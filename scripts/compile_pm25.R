# =============================================================================
# compile_pm25.R
#
# Purpose : Compile PM2.5 daily averages from CARB PM25HR site files
#           (PM25HR_SITE1YR_*.csv) and merge into Tahoe_PM25_Master.csv.
#           Fills missing Tahoe_City values from new files, rebuilds the
#           composite PM25 column, interpolates short gaps (<3 days), and
#           saves a processed master to data/processed/tahoe_pm25_compiled.csv.
#
# Inputs  : data/pm2.5/PM25HR_SITE1YR_*.csv   — CARB ARB hourly-site annual files
#           data/pm2.5/Tahoe_PM25_Master.csv   — existing master record
#
# Outputs : data/processed/tahoe_pm25_compiled.csv  — compiled + interpolated
#           figures/pm25/pm25_compiled_timeseries.png
#           figures/pm25/pm25_compiled_density.png
#           figures/pm25/pm25_compiled_density_log10.png
#           figures/pm25/pm25_compiled_density_fire_season.png
#           figures/pm25/pm25_compiled_density_fire_panels.png
#
# Site mapping (PM25HR `name` → master column):
#   "Tahoe City"              → Tahoe_City
#   "South Lake Tahoe"        → South_Lake
#   "Bliss State Park"        → Bliss
#   "Lake Tahoe College" etc. → LT_College
#   "Truckee"                 → Truckee
#
# Author  : Dante A. Capone
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(zoo)
  library(here)
  library(glue)
  library(patchwork)
})

# =============================================================================
# SECTION 0: Configuration
# =============================================================================

RAW_PM_DIR <- here("data", "pm2.5")
MASTER_FILE <- file.path(RAW_PM_DIR, "Tahoe_PM25_Master.csv")
PROC_DIR   <- here("data", "processed")
FIG_DIR    <- here("figures", "pm25")
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

# Master column names for each monitoring site
MASTER_SITE_COLS <- c("Bliss", "LT_College", "South_Lake", "Tahoe_City", "Truckee")

# Map PM25HR `name` field (lowercase) → master column
# Add variants here if new site names are encountered.
SITE_MAP <- function(nm) {
  nm <- tolower(trimws(nm))
  dplyr::case_when(
    str_detect(nm, "tahoe.{0,5}city")              ~ "Tahoe_City",
    str_detect(nm, "south.{0,5}lake")              ~ "South_Lake",
    str_detect(nm, "bliss")                        ~ "Bliss",
    str_detect(nm, "college|lt.?col|community")    ~ "LT_College",
    str_detect(nm, "truckee")                      ~ "Truckee",
    TRUE                                           ~ NA_character_
  )
}

# Fire seasons to annotate (year → label)
FIRE_EVENTS <- c(
  "2020" = "2020 Wildfires",
  "2021" = "Caldor Fire"
)

# Composite PM25 column priority: Tahoe_City > ARB_Daily_Average > other sites
PRIORITY_ORDER <- c("Tahoe_City", "ARB_Daily_Average",
                    "South_Lake", "Bliss", "LT_College", "Truckee",
                    "ARB_Daily_Max")

# =============================================================================
# SECTION 1: Read and standardise PM25HR site files
# =============================================================================

pm_files <- list.files(RAW_PM_DIR, pattern = "^PM25HR_SITE1YR.*\\.csv$",
                       full.names = TRUE)
cat(glue("Found {length(pm_files)} PM25HR file(s):\n"))
cat(paste0("  ", basename(pm_files), "\n"))

pm_new_long <- map_dfr(pm_files, function(f) {
  read_csv(f, show_col_types = FALSE) |>
    mutate(
      date      = as.Date(summary_date),
      pm25_davg = suppressWarnings(as.numeric(pm25_davg))
    ) |>
    filter(!is.na(date), !is.na(pm25_davg), pm25_davg >= 0, pm25_davg <= 500) |>
    select(date, site, name, pm25_davg, obs_count)
})

# Report unique sites found
site_names <- pm_new_long |> distinct(site, name) |> arrange(name)
cat("\nUnique sites in PM25HR files:\n")
print(site_names)

# Map to master column names
pm_new_long <- pm_new_long |>
  mutate(master_col = SITE_MAP(name))

unmapped <- pm_new_long |> filter(is.na(master_col)) |> distinct(name)
if (nrow(unmapped) > 0) {
  cat("\n\033[33mWARNING: Sites with no master column match (will be skipped):\033[0m\n")
  print(unmapped)
}

pm_new_long <- pm_new_long |> filter(!is.na(master_col))

# Collapse to one value per date × master_col (mean if multiple monitors)
pm_new_wide <- pm_new_long |>
  group_by(date, master_col) |>
  summarise(pm25 = mean(pm25_davg, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = master_col, values_from = pm25)

cat(glue("\nNew PM25HR data: {nrow(pm_new_wide)} site-days across ",
         "{format(min(pm_new_wide$date))} – {format(max(pm_new_wide$date))}\n"))

# =============================================================================
# SECTION 2: Read and parse the master CSV
# =============================================================================

cat("\nReading master file...\n")
master_raw <- read_csv(MASTER_FILE, show_col_types = FALSE)

master <- master_raw |>
  mutate(
    Date  = as.Date(Date),
    across(any_of(c("PM25", "ARB_Daily_Average", "ARB_Daily_Max",
                    MASTER_SITE_COLS)),
           ~ suppressWarnings(as.numeric(.x)))
  ) |>
  rename(date = Date) |>
  arrange(date)

cat(glue("Master records: {nrow(master)} days ",
         "({format(min(master$date, na.rm=TRUE))} – ",
         "{format(max(master$date, na.rm=TRUE))})\n"))

# =============================================================================
# SECTION 3: Merge new PM25HR data into master
# =============================================================================

# Ensure master has all site columns
for (col in MASTER_SITE_COLS) {
  if (!col %in% names(master)) master[[col]] <- NA_real_
}

# Separate master rows that overlap with new data vs. those that don't
new_dates  <- pm_new_wide$date
in_master  <- master$date %in% new_dates
new_only   <- !(new_dates %in% master$date)

# 3a. Update existing rows: fill NA site columns from pm_new_wide
update_cols <- intersect(MASTER_SITE_COLS, setdiff(names(pm_new_wide), "date"))

updated_master <- master |>
  left_join(pm_new_wide |>
              rename_with(~ paste0(".new_", .x), -date),
            by = "date") |>
  mutate(across(
    all_of(update_cols),
    ~ {
      new_col <- get(paste0(".new_", cur_column()))
      if_else(is.na(.x) & !is.na(new_col), new_col, .x)
    }
  )) |>
  select(-starts_with(".new_"))

# 3b. Append genuinely new dates
if (sum(new_only) > 0) {
  new_rows <- pm_new_wide |>
    filter(date %in% new_dates[new_only]) |>
    mutate(
      Julian          = yday(date),
      PM25            = NA_real_,
      ARB_Daily_Average = NA_real_,
      ARB_Daily_Max   = NA_real_,
      Source          = "PM25HR"
    )
  # add any missing site cols
  for (col in MASTER_SITE_COLS) {
    if (!col %in% names(new_rows)) new_rows[[col]] <- NA_real_
  }
  updated_master <- bind_rows(updated_master, new_rows) |> arrange(date)
  cat(glue("  Appended {sum(new_only)} new date(s) not previously in master.\n"))
}

cat(glue("  Updated {sum(in_master)} existing master rows with new PM25HR data.\n"))

# =============================================================================
# SECTION 4: Rebuild composite PM25 column (priority-based)
# =============================================================================

# Use first non-NA value from PRIORITY_ORDER
composite_cols <- intersect(PRIORITY_ORDER, names(updated_master))

updated_master <- updated_master |>
  mutate(
    PM25 = reduce(
      map(composite_cols, ~ updated_master[[.x]]),
      ~ if_else(is.na(.x), .y, .x)
    ),
    Source = reduce(
      seq_along(composite_cols),
      function(src_acc, i) {
        col <- composite_cols[[i]]
        if_else(
          is.na(src_acc) & !is.na(updated_master[[col]]),
          col,
          src_acc
        )
      },
      .init = NA_character_
    )
  )

n_pm25_before <- sum(!is.na(updated_master$PM25))
cat(glue("\nPM25 composite: {n_pm25_before} days with data ",
         "({round(100*n_pm25_before/nrow(updated_master),1)}% coverage)\n"))

# =============================================================================
# SECTION 5: Interpolate short gaps (<3 days) in PM25 composite
# =============================================================================

# Fill the full date sequence so zoo::na.approx sees true consecutive days
full_seq <- tibble(date = seq(min(updated_master$date, na.rm = TRUE),
                              max(updated_master$date, na.rm = TRUE),
                              by = "day"))

compiled <- full_seq |>
  left_join(updated_master, by = "date") |>
  mutate(
    PM25_raw        = PM25,
    PM25            = na.approx(PM25, na.rm = FALSE, maxgap = 2),
    interpolated    = is.na(PM25_raw) & !is.na(PM25),
    Source          = if_else(interpolated, "interpolated", Source)
  ) |>
  select(-PM25_raw)

n_interp <- sum(compiled$interpolated, na.rm = TRUE)
cat(glue("Interpolated {n_interp} gap day(s) (gaps of 1–2 days).\n"))

# =============================================================================
# SECTION 6: Save processed output
# =============================================================================

write_csv(compiled, file.path(PROC_DIR, "tahoe_pm25_compiled.csv"))
cat(glue("\nSaved: {file.path(PROC_DIR, 'tahoe_pm25_compiled.csv')}\n"))
cat(glue("  Total rows: {nrow(compiled)}  |  PM25 coverage: ",
         "{sum(!is.na(compiled$PM25))} days\n"))

# =============================================================================
# SECTION 7: Diagnostic multisite time series figure
# =============================================================================

# --- 7a. Compile long-format for individual site panels ---
site_long <- compiled |>
  select(date, all_of(MASTER_SITE_COLS)) |>
  pivot_longer(-date, names_to = "site", values_to = "pm25") |>
  mutate(
    year = year(date),
    site = factor(site, levels = MASTER_SITE_COLS)
  ) |>
  filter(!is.na(pm25))

# --- 7b. Colour palette per site ---
SITE_COLS <- c(
  Tahoe_City  = "#1B6CA8",
  South_Lake  = "#C0392B",
  Bliss       = "#27AE60",
  LT_College  = "#8E44AD",
  Truckee     = "#E67E22"
)

# --- 7c. Background shading bands: which site(s) contributed each day ---
source_band <- compiled |>
  select(date, Source) |>
  mutate(
    src_site = case_when(
      Source %in% MASTER_SITE_COLS ~ Source,
      Source == "ARB_Daily_Average" ~ "ARB_Avg",
      Source == "PM25HR"            ~ "Tahoe_City",
      Source == "interpolated"      ~ "interpolated",
      TRUE                          ~ "other"
    )
  )

# Contiguous bands per source (for geom_rect background)
source_rects <- source_band |>
  filter(!is.na(src_site)) |>
  mutate(run = cumsum(c(TRUE, diff(as.integer(factor(src_site))) != 0))) |>
  group_by(run, src_site) |>
  summarise(xmin = min(date) - 0.5, xmax = max(date) + 0.5,
            .groups = "drop")

BAND_COLS <- c(
  Tahoe_City    = "#D6EAF8",
  South_Lake    = "#FADBD8",
  Bliss         = "#D5F5E3",
  LT_College    = "#E8DAEF",
  Truckee       = "#FDEBD0",
  ARB_Avg       = "#FDFEFE",
  interpolated  = "#F9F9F9",
  other         = "#F2F3F4"
)

# Fire-event vertical lines
fire_lines <- tibble(
  year  = as.integer(names(FIRE_EVENTS)),
  label = unname(FIRE_EVENTS),
  x     = as.Date(paste0(year, "-08-01"))
)

# --- 7d. Panel A: compiled PM25 time series with source background ---
p_compiled <- ggplot() +
  geom_rect(
    data = source_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = src_site),
    alpha = 0.55, inherit.aes = FALSE
  ) +
  geom_vline(
    data = fire_lines,
    aes(xintercept = x),
    colour = "grey40", linetype = "dashed", linewidth = 0.4
  ) +
  annotate("text",
           x = fire_lines$x, y = Inf,
           label = fire_lines$label,
           vjust = 1.4, hjust = -0.05, size = 2.8, colour = "grey30") +
  geom_line(
    data = compiled |> filter(!is.na(PM25)),
    aes(x = date, y = PM25),
    colour = "grey20", linewidth = 0.25, alpha = 0.8
  ) +
  geom_point(
    data = compiled |> filter(interpolated),
    aes(x = date, y = PM25),
    shape = 21, fill = "white", colour = "grey50", size = 0.9
  ) +
  scale_fill_manual(
    values = BAND_COLS,
    name   = "Source site",
    na.value = "white"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Compiled Lake Tahoe PM\u2082.\u2085 (daily average)",
    subtitle = "Background shading = data source site; open circles = interpolated gaps",
    x = NULL, y = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(4, "mm"),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
    axis.text.x        = element_text(size = 8.5, angle = 0)
  )

# --- 7e. Panels B: individual site time series ---
p_sites <- ggplot(site_long, aes(x = date, y = pm25, colour = site)) +
  geom_vline(
    data = fire_lines |> mutate(dummy = TRUE),
    aes(xintercept = x),
    colour = "grey60", linetype = "dashed", linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  geom_hline(yintercept = c(16, 35), linetype = "dotted",
             colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = SITE_COLS, name = "Site", guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ site, ncol = 1, strip.position = "right") +
  labs(
    title    = "Individual monitoring sites",
    subtitle = "Dotted lines: Medium (16) and High (35) \u03bcg m\u207b\u00b3 thresholds",
    x = NULL, y = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.3),
    strip.text         = element_text(size = 8, face = "bold"),
    strip.background   = element_blank(),
    plot.title         = element_text(face = "bold", size = 10),
    plot.subtitle      = element_text(size = 8, colour = "grey40")
  )

# --- 7f. Assemble and save ---
p_full <- p_compiled / p_sites +
  plot_layout(heights = c(2, 3)) +
  plot_annotation(
    title  = "Lake Tahoe PM\u2082.\u2085 Compilation Diagnostics",
    theme  = theme(plot.title = element_text(face = "bold", size = 13))
  )

out_fig <- file.path(FIG_DIR, "pm25_compiled_timeseries.png")
ggsave(out_fig, plot = p_full, width = 26, height = 22, dpi = 250, units = "cm")
cat(glue("\nSaved figure: {out_fig}\n"))

# =============================================================================
# SECTION 8: Probability density of PM2.5 values by year
# =============================================================================

pm_density <- compiled |>
  filter(!is.na(PM25), PM25 > 0) |>
  mutate(year = factor(year(date), levels = sort(unique(year(date)))))

# --- 8a. Raw-scale density ---
p_density_raw <- ggplot(pm_density, aes(x = PM25, colour = year, group = year)) +
  geom_density(linewidth = 0.8, alpha = 0.7) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_viridis_d(option = "turbo", begin = 0, end = 0.95, name = "Year") +
  labs(
    title = "Probability density of daily PM\u2082.\u2085 by year",
    x = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3)",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 8.5),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

density_raw_fig <- file.path(FIG_DIR, "pm25_compiled_density.png")
ggsave(density_raw_fig, plot = p_density_raw, width = 26, height = 16, dpi = 300, units = "cm")
cat(glue("\nSaved figure: {density_raw_fig}\n"))

# --- 8b. Log10-scale density ---
p_density_log <- ggplot(pm_density, aes(x = PM25, colour = year, group = year)) +
  geom_density(linewidth = 0.8, alpha = 0.7) +
  scale_x_log10(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_viridis_d(option = "turbo", begin = 0, end = 0.95, name = "Year") +
  labs(
    title = "Probability density of daily PM\u2082.\u2085 by year (log\u2081\u2080 scale)",
    x = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3, log\u2081\u2080 scale)",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 8.5),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

density_log_fig <- file.path(FIG_DIR, "pm25_compiled_density_log10.png")
ggsave(density_log_fig, plot = p_density_log, width = 26, height = 16, dpi = 300, units = "cm")
cat(glue("\nSaved figure: {density_log_fig}\n"))

# --- 8c. Fire-season density: fire years vs. long-term mean (Jun-Oct) ---
FIRE_YEARS <- c(2007, 2008, 2018, 2020, 2021, 2022)

pm_fire_season <- compiled |>
  filter(!is.na(PM25), PM25 > 0, month(date) %in% 6:10) |>
  mutate(year = year(date))

pm_fire_mean <- pm_fire_season |>
  mutate(year_label = "Mean")

pm_fire_years <- pm_fire_season |>
  filter(year %in% FIRE_YEARS) |>
  mutate(year_label = factor(as.character(year), levels = as.character(FIRE_YEARS)))

FIRE_COLOURS <- c(
  "Mean" = "#7F7F7F",
  "2007" = "#F781BF",
  "2008" = "#E41A1C",
  "2018" = "#FF7F00",
  "2020" = "#4DAF4A",
  "2021" = "#984EA3",
  "2022" = "#377EB8"
)

p_fire_density <- ggplot() +
  geom_density(
    data = pm_fire_mean,
    aes(x = PM25, colour = year_label, fill = year_label),
    linewidth = 0.9, alpha = 0.25
  ) +
  geom_density(
    data = pm_fire_years,
    aes(x = PM25, colour = year_label, fill = year_label),
    linewidth = 0.8, alpha = 0.35
  ) +
  scale_x_log10(
    breaks = c(1, 5, 10, 20, 50, 100, 200),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = FIRE_COLOURS, name = "Year") +
  scale_fill_manual(values = FIRE_COLOURS, name = "Year") +
  labs(
    title = "Fire-season (Jun\u2013Oct) PM\u2082.\u2085 probability densities",
    subtitle = "Fire years highlighted against the 2006\u20132025 mean",
    x = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3, log\u2081\u2080 scale)",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 8.5),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

fire_density_fig <- file.path(FIG_DIR, "pm25_compiled_density_fire_season.png")
ggsave(fire_density_fig, plot = p_fire_density, width = 26, height = 16, dpi = 300, units = "cm")
cat(glue("\nSaved figure: {fire_density_fig}\n"))

# --- 8d. Multi-panel fire-year vs. mean with KS significance ---
ks_results <- map_dfr(FIRE_YEARS, function(fy) {
  fire_vals <- pm_fire_season |> filter(year == fy) |> pull(PM25)
  mean_vals <- pm_fire_season |> filter(year != fy) |> pull(PM25)
  test <- ks.test(fire_vals, mean_vals)
  tibble(
    year = fy,
    D = as.numeric(test$statistic),
    p_value = test$p.value,
    significant = test$p.value < 0.05
  )
})

cat("\nKS tests: fire-year vs. mean Jun-Oct PM2.5 distribution\n")
print(ks_results)

# Build long-format data for multi-panel plot
pm_fire_panels <- map_dfr(FIRE_YEARS, function(fy) {
  bind_rows(
    pm_fire_season |>
      filter(year == fy) |>
      mutate(curve = as.character(fy), panel_year = as.character(fy)),
    pm_fire_season |>
      filter(year != fy) |>
      mutate(curve = "Mean", panel_year = as.character(fy))
  )
})

# Panel titles with significance asterisks
panel_titles <- ks_results |>
  mutate(
    panel_title = paste0(year, ifelse(significant, " *", "")),
    panel_year = as.character(year)
  )

pm_fire_panels <- pm_fire_panels |>
  left_join(panel_titles |> select(panel_year, panel_title), by = "panel_year")

panel_colors <- c("Mean" = "#7F7F7F", FIRE_COLOURS)

p_fire_panels <- ggplot(pm_fire_panels,
                        aes(x = PM25, colour = curve, fill = curve)) +
  geom_density(linewidth = 0.8, alpha = 0.3) +
  scale_x_log10(
    breaks = c(1, 5, 10, 20, 50, 100, 200),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = panel_colors, name = NULL) +
  scale_fill_manual(values = panel_colors, name = NULL) +
  facet_wrap(~ panel_title, ncol = 3) +
  labs(
    title = "Fire-year PM\u2082.\u2085 distributions vs. long-term mean (Jun\u2013Oct)",
    subtitle = "* = significantly different from mean (KS test, p < 0.05)",
    x = "PM\u2082.\u2085 (\u03bcg m\u207b\u00b3, log\u2081\u2080 scale)",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "none",
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

fire_panels_fig <- file.path(FIG_DIR, "pm25_compiled_density_fire_panels.png")
ggsave(fire_panels_fig, plot = p_fire_panels, width = 26, height = 18, dpi = 300, units = "cm")
cat(glue("\nSaved figure: {fire_panels_fig}\n"))

cat("\nDone.\n")
