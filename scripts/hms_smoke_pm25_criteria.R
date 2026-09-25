# =============================================================================
# hms_smoke_pm25_criteria.R
#
# Purpose : Lake Tahoe smoke-day classification using a dual criterion:
#             (1) NOAA HMS smoke detected over the lake (Low / Medium / High)
#             (2) Surface EPA PM2.5 confirms smoke at the surface
#
# HMS density → estimated surface PM2.5 concentrations (NOAA HMS AOD codes):
#   Low    (AOD  5) ≈  5 μg m⁻³   |  formerly "Light"
#   Medium (AOD 16) ≈ 16 μg m⁻³   |  formerly "Medium"
#   High   (AOD 27) ≈ 27 μg m⁻³   |  formerly "Heavy"
#
# Smoke day criterion (any of three levels):
#   Low smoke day    : HMS ≥ Low    AND PM2.5 ≥  5 μg m⁻³
#   Medium smoke day : HMS ≥ Medium AND PM2.5 ≥ 16 μg m⁻³
#   High smoke day   : HMS = High   AND PM2.5 ≥ 27 μg m⁻³
#
# A "confirmed smoke day" (primary metric) = HMS ≥ Low AND PM2.5 ≥ 5 μg m⁻³.
# This runs in parallel to hms_smoke_swrad_criteria.R for comparison.
#
# Outputs : data/processed/tahoe_hms_pm25_daily.csv
#           data/processed/tahoe_hms_pm25_annual.csv
#           figures/hms_smoke/tahoe_smoke_pm25_sw_parallel.png  (parallel comparison)
#
# Author  : Dante A. Capone
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(here)
  library(glue)
  library(patchwork)
})

# =============================================================================
# SECTION 0: Configuration
# =============================================================================
PROC_DIR <- here("data", "processed")
FIG_DIR  <- here("figures", "hms_smoke")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

YEARS          <- 2006:2025
DENSITY_LEVELS <- c("Low", "Medium", "High")

# Estimated surface PM2.5 concentrations per HMS density category (μg m⁻³)
HMS_PM_CONC <- c(Low = 5, Medium = 16, High = 27)

# PM2.5 confirmation thresholds tied to HMS categories
THR_LOW    <-  5.0   # μg m⁻³ — Low density category
THR_MEDIUM <- 16.0   # μg m⁻³ — Medium density category
THR_HIGH   <- 27.0   # μg m⁻³ — High density category

# EPA AQS parameters
PARAM_FRM  <- "88101"
PARAM_NFRM <- "88502"
PM_MIN     <-   0.0
PM_MAX     <- 500.0
EXCL_COUNTY <- "Nevada"
EXCL_STATE  <- "California"

# =============================================================================
# SECTION 1: Load HMS daily smoke data
# =============================================================================
cat("Loading HMS daily data...\n")
daily_hms <- read_csv(
  file.path(PROC_DIR, "tahoe_hms_smoke_daily.csv"),
  show_col_types = FALSE
) |>
  mutate(
    date        = as.Date(date),
    max_density = factor(max_density, levels = DENSITY_LEVELS)
  ) |>
  filter(year %in% YEARS)

# =============================================================================
# SECTION 2: Load compiled PM2.5 daily data
# =============================================================================
# Source: data/processed/tahoe_pm25_compiled.csv  (run compile_pm25.R first)
# Falls back to tahoe_epa_pm_daily_all_years.parquet if compiled file missing.
cat("Loading PM2.5 data...\n")

compiled_pm_path <- file.path(PROC_DIR, "tahoe_pm25_compiled.csv")

if (file.exists(compiled_pm_path)) {
  cat("  Using compiled master: tahoe_pm25_compiled.csv\n")
  pm_compiled <- read_csv(compiled_pm_path, show_col_types = FALSE) |>
    mutate(date = as.Date(date))

  daily_pm <- pm_compiled |>
    filter(!is.na(PM25), PM25 >= PM_MIN, PM25 <= PM_MAX) |>
    mutate(
      year  = year(date),
      month = month(date)
    ) |>
    filter(year %in% YEARS) |>
    group_by(date, year, month) |>
    summarise(pm25_max = max(PM25, na.rm = TRUE), .groups = "drop")

} else {
  cat("\033[33m  WARNING: tahoe_pm25_compiled.csv not found — run compile_pm25.R first.\033[0m\n")
  cat("  Falling back to tahoe_epa_pm_daily_all_years.parquet\n")
  library(arrow)
  pm_raw <- read_parquet(file.path(PROC_DIR, "tahoe_epa_pm_daily_all_years.parquet"))

  daily_pm <- pm_raw |>
    filter(
      parameter_code %in% c(PARAM_FRM, PARAM_NFRM),
      !is.na(arithmetic_mean),
      arithmetic_mean >= PM_MIN,
      arithmetic_mean <= PM_MAX,
      !(county == EXCL_COUNTY & state == EXCL_STATE)
    ) |>
    mutate(
      date  = as.Date(date_local),
      year  = year(date),
      month = month(date)
    ) |>
    filter(year %in% YEARS) |>
    group_by(date, year, month) |>
    summarise(pm25_max = max(arithmetic_mean, na.rm = TRUE), .groups = "drop")
}

cat(sprintf("  PM2.5 records loaded: %d days\n", nrow(daily_pm)))

# =============================================================================
# SECTION 3: Merge HMS and PM2.5
# =============================================================================
merged <- daily_hms |>
  left_join(daily_pm |> select(date, pm25_max), by = "date") |>
  mutate(
    # PM2.5 flags at each concentration threshold
    pm_gte_low    = !is.na(pm25_max) & pm25_max >= THR_LOW,
    pm_gte_medium = !is.na(pm25_max) & pm25_max >= THR_MEDIUM,
    pm_gte_high   = !is.na(pm25_max) & pm25_max >= THR_HIGH,

    # Dual-criteria smoke days at each level
    smoke_low    = smoke_day & pm_gte_low,    # HMS any + PM ≥ 5 μg m⁻³
    smoke_medium = smoke_day & pm_gte_medium, # HMS any + PM ≥ 16 μg m⁻³
    smoke_high   = smoke_day & pm_gte_high,   # HMS any + PM ≥ 27 μg m⁻³

    # Strict: HMS ≥ Medium/High AND PM at corresponding threshold
    smoke_med_pm = !is.na(max_density) & max_density %in% c("Medium", "High") &
                   pm_gte_medium,
    smoke_high_pm = !is.na(max_density) & max_density == "High" & pm_gte_high,

    # Primary confirmed smoke day: HMS any density AND PM2.5 ≥ Low threshold
    smoke_confirmed = smoke_day & pm_gte_low
  )

n_no_pm <- sum(is.na(merged$pm25_max) & merged$smoke_day)
cat(sprintf("  HMS smoke days without PM2.5 data: %d (will be unconfirmed)\n", n_no_pm))

# =============================================================================
# SECTION 4: Annual summary Jul 1 – Sep 30
# =============================================================================
summer <- merged |> filter(month(date) %in% 7:9)

annual_pm25 <- summer |>
  group_by(year) |>
  summarise(
    hms_any_days      = sum(smoke_day,        na.rm = TRUE),
    pm25_days_avail   = sum(!is.na(pm25_max), na.rm = TRUE),
    confirmed_days    = sum(smoke_confirmed,  na.rm = TRUE),   # HMS + PM ≥ 5 μg m⁻³
    medium_pm_days    = sum(smoke_med_pm,     na.rm = TRUE),   # HMS Med/High + PM ≥ 16
    high_pm_days      = sum(smoke_high_pm,    na.rm = TRUE),   # HMS High + PM ≥ 27
    .groups = "drop"
  )

cat("\nHMS × PM2.5 dual-criteria smoke days (Jul 1 – Sep 30):\n")
print(annual_pm25, n = Inf)

# =============================================================================
# SECTION 5: Compare to SW criterion (load parallel output)
# =============================================================================
swrad_exists <- file.exists(file.path(PROC_DIR, "tahoe_hms_swrad_annual.csv"))
if (swrad_exists) {
  annual_sw <- read_csv(
    file.path(PROC_DIR, "tahoe_hms_swrad_annual.csv"),
    show_col_types = FALSE
  )
  cat("\nSW radiation criterion (Jul 1 – Sep 30):\n")
  print(annual_sw |> select(year, hms_med_high_days, smits_smoke_days), n = Inf)
} else {
  cat("Note: tahoe_hms_swrad_annual.csv not found — run hms_smoke_swrad_criteria.R first.\n")
  annual_sw <- NULL
}

# =============================================================================
# SECTION 6: Save processed outputs
# =============================================================================
write_csv(
  merged |> select(date, year, month, smoke_day, max_density, pm25_max,
                   pm_gte_low, pm_gte_medium, pm_gte_high,
                   smoke_confirmed, smoke_med_pm, smoke_high_pm),
  file.path(PROC_DIR, "tahoe_hms_pm25_daily.csv")
)
write_csv(annual_pm25, file.path(PROC_DIR, "tahoe_hms_pm25_annual.csv"))
cat("\nSaved: tahoe_hms_pm25_daily.csv and tahoe_hms_pm25_annual.csv\n")

# =============================================================================
# SECTION 7: Parallel comparison figure — PM2.5 criterion vs SW criterion
# =============================================================================

# --- Panel A: HMS × PM2.5 annual smoke days (Jul-Sep) ---
LBL_HMS  <- "HMS any density"
LBL_LOW  <- paste0("HMS + PM\u2082.\u2085 \u2265 5 \u03bcg m\u207b\u00b3")
LBL_MED  <- paste0("HMS + PM\u2082.\u2085 \u2265 16 \u03bcg m\u207b\u00b3")
LBL_SW   <- paste0("HMS + SW deficit >20 W m\u207b\u00b2")

pm25_bar_data <- annual_pm25 |>
  select(year, hms_any_days, confirmed_days, medium_pm_days) |>
  pivot_longer(-year, names_to = "criterion", values_to = "n_days") |>
  mutate(criterion = recode(criterion,
    hms_any_days   = LBL_HMS,
    confirmed_days = LBL_LOW,
    medium_pm_days = LBL_MED
  ),
  criterion = factor(criterion, levels = c(LBL_HMS, LBL_LOW, LBL_MED)))

p_pm25 <- ggplot(pm25_bar_data,
                 aes(x = factor(year), y = n_days, fill = criterion)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = setNames(c("#FED98E", "#FB8C00", "#B71C1C"),
                       c(LBL_HMS, LBL_LOW, LBL_MED)),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "HMS \u00d7 PM\u2082.\u2085 dual-criteria smoke days (Jul\u2013Sep, 2006\u20132025)",
    subtitle = "HMS any density vs confirmed by surface PM\u2082.\u2085 at 5 and 16 \u03bcg m\u207b\u00b3 thresholds",
    x = NULL, y = "Smoke days"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 8.5),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.4)
  )

# --- Panel B: Parallel comparison — PM2.5 vs SW criterion (if available) ---
if (!is.null(annual_sw)) {
  comp_data <- left_join(
    annual_pm25 |> select(year, confirmed_days, medium_pm_days),
    annual_sw   |> select(year, smits_smoke_days),
    by = "year"
  ) |>
    pivot_longer(-year, names_to = "criterion", values_to = "n_days") |>
    mutate(criterion = recode(criterion,
      confirmed_days   = LBL_LOW,
      medium_pm_days   = LBL_MED,
      smits_smoke_days = LBL_SW
    ))

  p_parallel <- ggplot(comp_data,
                       aes(x = factor(year), y = n_days, fill = criterion)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75,
             colour = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = setNames(c("#FB8C00", "#B71C1C", "#1565C0"),
                         c(LBL_LOW, LBL_MED, LBL_SW)),
      name = NULL
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = "Parallel smoke-day criteria comparison (Jul\u2013Sep)",
      subtitle = "PM\u2082.\u2085-based (orange/red) vs SW radiation-based (blue)",
      x = NULL, y = "Smoke days"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "top",
      legend.text        = element_text(size = 8.5),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9, colour = "grey40"),
      panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.4)
    )

  p_out <- p_pm25 / p_parallel +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.margin = margin()))
} else {
  p_out <- p_pm25
}

out_fig <- file.path(FIG_DIR, "tahoe_smoke_pm25_sw_parallel.png")
ggsave(out_fig, plot = p_out, width = 20, height = if (!is.null(annual_sw)) 13 else 7,
       dpi = 200, units = "cm")
cat(glue("Saved: {out_fig}\n"))

cat("\nDone.\n")
