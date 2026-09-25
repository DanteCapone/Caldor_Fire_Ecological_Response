# =============================================================================
# 04_usgs_tahoe_tributary_streamflow.R
#
# Purpose: Download USGS daily mean discharge for two Lake Tahoe tributaries
#          and compare July 2021-July 2022 with a 2015-2025 climatology.
#
# Gauges:
#   10336610  Upper Truckee River at South Lake Tahoe, CA
#   10336660  Blackwood Creek near Tahoe City, CA
#
# USGS parameter/statistic:
#   00060 = discharge; 00003 = daily mean
#
# The climatology uses all daily observations from 2015-2025 except dates in
# the focal interval (2021-07-01 through 2022-07-31), so the focal observations
# do not contribute to their own reference distribution.
#
# Outputs:
#   data/processed/usgs_tahoe_tributary_streamflow_daily_2015_2025.csv
#   data/processed/usgs_tahoe_tributary_streamflow_daily_climatology.csv
#   data/processed/usgs_tahoe_tributary_streamflow_monthly_comparison.csv
#   data/processed/usgs_tahoe_tributary_streamflow_summary.csv
#   figures/streamflow/tahoe_tributary_streamflow_daily_jul2021_jul2022.png
#   figures/streamflow/tahoe_tributary_streamflow_monthly_jul2021_jul2022.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(ragg)
})

# ---- Project paths -----------------------------------------------------------
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
  proj_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
} else {
  proj_root <- normalizePath(".", winslash = "/")
}
setwd(proj_root)
source("scripts/figure_aesthetics.R")

proc_dir <- file.path(proj_root, "data", "processed")
fig_dir  <- file.path(proj_root, "figures", "streamflow")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- Analysis configuration -------------------------------------------------
start_date  <- as.Date("2015-01-01")
end_date    <- as.Date("2025-12-31")
focus_start <- as.Date("2021-07-01")
focus_end   <- as.Date("2022-07-31")
fire_start  <- as.Date("2021-08-14")
fire_end    <- as.Date("2021-10-21")

sites <- tribble(
  ~site_no,    ~site,
  "10336610",  "Upper Truckee River",
  "10336660",  "Blackwood Creek"
) %>%
  mutate(
    site = factor(site, levels = c("Upper Truckee River", "Blackwood Creek"))
  )

site_cols <- c(
  "Upper Truckee River" = "#1565C0",
  "Blackwood Creek"     = "#2E7D32"
)

# ---- Download USGS daily values ---------------------------------------------
# Use the public NWIS daily-values service directly to keep this script
# reproducible without requiring the optional dataRetrieval package.
usgs_url <- paste0(
  "https://waterservices.usgs.gov/nwis/dv/?format=rdb",
  "&sites=", paste(sites$site_no, collapse = "%2C"),
  "&startDT=", start_date,
  "&endDT=", end_date,
  "&parameterCd=00060&statCd=00003&siteStatus=all"
)

message("Downloading USGS daily mean discharge...")
message(usgs_url)

cache_file <- file.path(
  proc_dir, "usgs_tahoe_tributary_streamflow_daily_2015_2025.csv"
)

raw <- tryCatch(
  readr::read_tsv(
    usgs_url,
    comment = "#",
    col_types = cols(.default = col_character()),
    progress = FALSE,
    name_repair = "minimal",
    show_col_types = FALSE
  ),
  error = function(e) {
    if (!file.exists(cache_file)) stop(e)
    warning(
      "USGS download failed; using cached 2015-2025 daily data: ",
      conditionMessage(e),
      call. = FALSE
    )
    read_csv(cache_file, show_col_types = FALSE) %>%
      transmute(
        agency_cd = "USGS",
        site_no = as.character(site_no),
        datetime = as.character(date),
        `00060_00003` = as.character(flow_cfs),
        `00060_00003_cd` = as.character(qualifier)
      )
  }
)

flow_col <- names(raw)[
  grepl("00060_00003$", names(raw)) &
    !grepl("_cd$", names(raw))
]
qual_col <- names(raw)[grepl("00060_00003_cd$", names(raw))]

if (length(flow_col) != 1L) {
  stop(
    "Could not uniquely identify the USGS daily-discharge column. Columns: ",
    paste(names(raw), collapse = ", ")
  )
}

flow <- raw %>%
  filter(agency_cd == "USGS", site_no %in% sites$site_no) %>%
  transmute(
    site_no,
    date = as.Date(datetime),
    flow_cfs = readr::parse_double(.data[[flow_col]]),
    qualifier = if (length(qual_col) == 1L) .data[[qual_col]] else NA_character_
  ) %>%
  left_join(sites, by = "site_no") %>%
  mutate(
    flow_cms = flow_cfs * 0.028316846592,
    year = year(date),
    month = month(date),
    month_day = format(date, "%m-%d"),
    water_year = if_else(month >= 10L, year + 1L, year)
  ) %>%
  arrange(site, date)

if (nrow(flow) == 0L || all(is.na(flow$flow_cfs))) {
  stop("USGS returned no usable daily discharge observations.")
}

write_csv(
  flow,
  file.path(proc_dir, "usgs_tahoe_tributary_streamflow_daily_2015_2025.csv"),
  na = ""
)

# ---- Daily focal series and 30-day moving-median climatology ----------------
moving_median <- function(x, window = 30L) {
  left <- floor((window - 1L) / 2L)
  right <- window - 1L - left
  vapply(seq_along(x), function(i) {
    idx <- seq.int(max(1L, i - left), min(length(x), i + right))
    median(x[idx], na.rm = TRUE)
  }, numeric(1))
}

focus_dates <- tibble(
  plot_date = seq(focus_start, focus_end, by = "day"),
  month_day = format(plot_date, "%m-%d")
)

daily_clim_stats <- flow %>%
  filter(
    !is.na(flow_cms),
    !(date >= focus_start & date <= focus_end)
  ) %>%
  group_by(site_no, site, month_day) %>%
  summarise(
    n_years = n_distinct(year),
    clim_daily_mean_cms = mean(flow_cms),
    clim_daily_sd_cms = sd(flow_cms),
    .groups = "drop"
  ) %>%
  mutate(
    clim_daily_se_cms = clim_daily_sd_cms / sqrt(n_years),
    clim_daily_ci_cms = if_else(
      n_years > 1L,
      qt(0.975, df = n_years - 1L) * clim_daily_se_cms,
      NA_real_
    ),
    clim_daily_lo_cms = pmax(clim_daily_mean_cms - clim_daily_ci_cms, 0),
    clim_daily_hi_cms = clim_daily_mean_cms + clim_daily_ci_cms
  )

daily_clim <- crossing(
  sites %>% select(site_no, site),
  focus_dates
) %>%
  left_join(
    daily_clim_stats,
    by = c("site_no", "site", "month_day")
  ) %>%
  group_by(site_no, site) %>%
  arrange(plot_date, .by_group = TRUE) %>%
  mutate(
    clim_center_cms = moving_median(clim_daily_mean_cms, 30L),
    clim_lo_cms = moving_median(clim_daily_lo_cms, 30L),
    clim_hi_cms = moving_median(clim_daily_hi_cms, 30L)
  ) %>%
  ungroup()

focus_daily <- flow %>%
  filter(date >= focus_start, date <= focus_end) %>%
  select(site_no, site, date, flow_cfs, flow_cms, qualifier) %>%
  left_join(
    daily_clim %>%
      select(site_no, plot_date, n_years, starts_with("clim_")),
    by = c("site_no", "date" = "plot_date")
  ) %>%
  mutate(anomaly_cms = flow_cms - clim_center_cms)

write_csv(
  focus_daily,
  file.path(proc_dir, "usgs_tahoe_tributary_streamflow_daily_climatology.csv"),
  na = ""
)
# ---- Monthly focal series and monthly climatology ---------------------------
# First average daily values within each site x year x month. The monthly
# climatology then treats each year-month mean as one independent replicate.
flow_monthly <- flow %>%
  filter(!is.na(flow_cms)) %>%
  mutate(month_date = floor_date(date, "month")) %>%
  group_by(site_no, site, year, month, month_date) %>%
  summarise(
    n_days = n(),
    flow_mean_cms = mean(flow_cms),
    .groups = "drop"
  )

monthly_clim <- flow_monthly %>%
  filter(!(month_date >= floor_date(focus_start, "month") &
             month_date <= floor_date(focus_end, "month"))) %>%
  group_by(site_no, site, month) %>%
  summarise(
    clim_n_years = n(),
    clim_mean_cms = mean(flow_mean_cms),
    clim_sd_cms = sd(flow_mean_cms),
    .groups = "drop"
  ) %>%
  mutate(
    clim_se_cms = clim_sd_cms / sqrt(clim_n_years),
    clim_ci_cms = if_else(
      clim_n_years > 1L,
      qt(0.975, df = clim_n_years - 1L) * clim_se_cms,
      NA_real_
    ),
    clim_lo_cms = pmax(clim_mean_cms - clim_ci_cms, 0),
    clim_hi_cms = clim_mean_cms + clim_ci_cms
  )

focus_monthly <- flow_monthly %>%
  filter(
    month_date >= floor_date(focus_start, "month"),
    month_date <= floor_date(focus_end, "month")
  ) %>%
  left_join(
    monthly_clim,
    by = c("site_no", "site", "month")
  ) %>%
  mutate(
    anomaly_cms = flow_mean_cms - clim_mean_cms,
    anomaly_pct = 100 * anomaly_cms / clim_mean_cms
  ) %>%
  arrange(site, month_date)

write_csv(
  focus_monthly,
  file.path(proc_dir, "usgs_tahoe_tributary_streamflow_monthly_comparison.csv"),
  na = ""
)

# ---- Compact analysis summary -----------------------------------------------
analysis_summary <- focus_daily %>%
  filter(!is.na(flow_cms)) %>%
  group_by(site_no, site) %>%
  summarise(
    n_observed_days = n(),
    mean_observed_cms = mean(flow_cms),
    mean_climatology_cms = mean(clim_center_cms, na.rm = TRUE),
    mean_anomaly_cms = mean(anomaly_cms, na.rm = TRUE),
    mean_anomaly_pct = 100 * mean_anomaly_cms / mean_climatology_cms,
    peak_observed_cms = max(flow_cms),
    peak_observed_date = date[which.max(flow_cms)],
    .groups = "drop"
  ) %>%
  left_join(
    focus_monthly %>%
      group_by(site_no, site) %>%
      summarise(
        months_above_climatology = sum(anomaly_cms > 0, na.rm = TRUE),
        months_below_climatology = sum(anomaly_cms < 0, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("site_no", "site")
  )

write_csv(
  analysis_summary,
  file.path(proc_dir, "usgs_tahoe_tributary_streamflow_summary.csv"),
  na = ""
)

# ---- Shared figure theme -----------------------------------------------------
stream_theme <- lo_theme(base_size = 9, family = LO_FONT) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.2, "mm"),
    legend.key.width = unit(6, "mm"),
    legend.key.height = unit(3, "mm"),
    axis.text = element_text(size = 8.5),
    axis.title.y = element_text(size = 8.5),
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

add_common_layers <- function(p) {
  p +
    annotate(
      "rect",
      xmin = fire_start, xmax = fire_end,
      ymin = -Inf, ymax = Inf,
      fill = "firebrick", alpha = 0.07
    ) +
    geom_vline(
      xintercept = fire_start,
      colour = "firebrick", linetype = "dashed", linewidth = LO_LW_MID
    ) +
    annotate(
      "text",
      x = fire_start + (fire_end - fire_start) / 2,
      y = Inf,
      label = "Caldor Fire",
      colour = "firebrick",
      family = LO_FONT,
      size = 3,
      vjust = 1.25
    )
}

# ---- Daily figure ------------------------------------------------------------
p_daily <- ggplot()
p_daily <- add_common_layers(p_daily) +
  geom_ribbon(
    data = daily_clim,
    aes(
      x = plot_date, ymin = clim_lo_cms, ymax = clim_hi_cms,
      fill = site, group = site
    ),
    alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = daily_clim,
    aes(
      x = plot_date, y = clim_center_cms, colour = site,
      linetype = "Climatology", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = focus_daily,
    aes(
      x = date, y = flow_cms, colour = site,
      linetype = "2021-2022", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.95, na.rm = TRUE
  ) +
  scale_colour_manual(values = site_cols, name = NULL) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_linetype_manual(
    values = c("2021-2022" = "solid", "Climatology" = "solid"),
    breaks = c("2021-2022", "Climatology"),
    name = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(alpha = 1, linewidth = 0.9)
    ),
    linetype = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(
        colour = "grey30", alpha = c(0.95, 0.34), linewidth = 0.9
      )
    )
  ) +
  scale_x_date(
    limits = c(focus_start, focus_end),
    date_breaks = "1 month",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    x = NULL,
    y = expression("Discharge (m"^3~"s"^{-1}*")")
  ) +
  stream_theme
save_lo_fig(
  p_daily,
  file.path(fig_dir, "tahoe_tributary_streamflow_daily_jul2021_jul2022.png"),
  width_type = "double",
  height_cm = 9.5
)

# ---- Monthly figure (matching Figure 2's monthly climatology convention) ----
p_monthly <- ggplot()
p_monthly <- add_common_layers(p_monthly) +
  geom_ribbon(
    data = focus_monthly,
    aes(
      x = month_date, ymin = clim_lo_cms, ymax = clim_hi_cms,
      fill = site, group = site
    ),
    alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = focus_monthly,
    aes(
      x = month_date, y = clim_mean_cms, colour = site,
      linetype = "Climatology", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = focus_monthly,
    aes(
      x = month_date, y = flow_mean_cms, colour = site,
      linetype = "2021-2022", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.95, na.rm = TRUE
  ) +
  geom_point(
    data = focus_monthly,
    aes(x = month_date, y = flow_mean_cms, colour = site),
    size = 1.6, alpha = 0.95, na.rm = TRUE
  ) +
  scale_colour_manual(values = site_cols, name = NULL) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_linetype_manual(
    values = c("2021-2022" = "solid", "Climatology" = "solid"),
    breaks = c("2021-2022", "Climatology"),
    name = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(alpha = 1, linewidth = 0.9, shape = NA)
    ),
    linetype = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(
        colour = "grey30", alpha = c(0.95, 0.34), linewidth = 0.9
      )
    )
  ) +
  scale_x_date(
    limits = c(floor_date(focus_start, "month"), floor_date(focus_end, "month")),
    date_breaks = "1 month",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    x = NULL,
    y = expression("Discharge (m"^3~"s"^{-1}*")")
  ) +
  stream_theme
save_lo_fig(
  p_monthly,
  file.path(fig_dir, "tahoe_tributary_streamflow_monthly_jul2021_jul2022.png"),
  width_type = "double",
  height_cm = 9.5
)

# ---- Console diagnostics -----------------------------------------------------
cat("\nUSGS records downloaded:\n")
print(
  flow %>%
    group_by(site_no, site) %>%
    summarise(
      first_date = min(date),
      last_date = max(date),
      n_days = n(),
      missing_flow = sum(is.na(flow_cms)),
      .groups = "drop"
    )
)

cat("\nFocal-period summary:\n")
print(analysis_summary)
cat("\nOutputs written to:\n  ", proc_dir, "\n  ", fig_dir, "\n", sep = "")
