# =============================================================================
# 04_usgs_streamflow_nutrient_timeseries.R
#
# Purpose: Combine the July 2021-July 2022 Upper Truckee/Blackwood discharge
#          comparison with discrete Upper Truckee nutrient time series and
#          pre-fire monthly climatologies.
#
# Inputs:
#   data/processed/usgs_tahoe_tributary_streamflow_daily_climatology.csv
#   data/processed/upper_truckee_water_quality.csv
#
# Outputs:
#   data/processed/upper_truckee_nutrient_timeseries_2021_2022.csv
#   data/processed/upper_truckee_nutrient_timeseries_availability.csv
#   figures/streamflow/tahoe_tributary_streamflow_nutrients_2021_2022.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(dataRetrieval)
  library(patchwork)
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
fig_dir <- file.path(proj_root, "figures", "streamflow")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

flow_file <- file.path(
  proc_dir, "usgs_tahoe_tributary_streamflow_daily_climatology.csv"
)
wq_file <- file.path(proc_dir, "upper_truckee_water_quality.csv")
wq_2022_cache <- file.path(
  proc_dir, "upper_truckee_water_quality_2022_plot_cache.csv"
)
plot_data_file <- file.path(
  proc_dir, "upper_truckee_nutrient_timeseries_2021_2022.csv"
)
availability_file <- file.path(
  proc_dir, "upper_truckee_nutrient_timeseries_availability.csv"
)
figure_file <- file.path(
  fig_dir, "tahoe_tributary_streamflow_nutrients_2021_2022.png"
)

if (!file.exists(flow_file)) {
  stop("Missing streamflow climatology: ", flow_file)
}
if (!file.exists(wq_file)) {
  stop("Missing Upper Truckee water-quality data: ", wq_file)
}

# ---- Analysis configuration -------------------------------------------------
focus_start <- as.Date("2021-07-01")
focus_end <- as.Date("2022-07-31")
fire_start <- as.Date("2021-08-14")
fire_end <- as.Date("2021-10-21")
baseline_start <- as.Date("2000-01-01")
baseline_end <- as.Date("2020-12-31")

site_cols <- c(
  "Upper Truckee River" = "#1565C0",
  "Blackwood Creek" = "#2E7D32"
)
upper_blue <- unname(site_cols[["Upper Truckee River"]])

# Use the same canonical USGS parameter/fraction combinations as the runoff
# comparison. All retained values are concentrations in mg/L.
nutrient_spec <- tribble(
  ~requested_variable, ~usgs_parameter_code, ~sample_fraction, ~result_units, ~plot_label,
  "Nitrate + nitrite", "00631", "Filtered field and/or lab", "mg/L", "NO3 + NO2 as N",
  "Ammonium/ammonia", "00608", "Filtered field and/or lab", "mg/L", "NH4/NH3 as N",
  "Total nitrogen", "00600", "Unfiltered", "mg/L", "Total nitrogen",
  "Kjeldahl/organic nitrogen", "00625", "Unfiltered", "mg/L", "TKN/organic N",
  "Orthophosphate", "00671", "Filtered field and/or lab", "mg/L", "Orthophosphate as P",
  "Total phosphorus", "00665", "Unfiltered", "mg/L", "Total phosphorus",
  "Dissolved phosphorus", "00666", "Filtered field and/or lab", "mg/L", "Dissolved phosphorus"
) %>%
  mutate(
    plot_label = factor(plot_label, levels = plot_label),
    facet_label = paste0(plot_label, "\n(mg L^-1)")
  )

# ---- Existing streamflow data -----------------------------------------------
flow_plot <- read_csv(
  flow_file,
  col_types = cols(
    site_no = col_character(),
    site = col_character(),
    date = col_date(),
    qualifier = col_character(),
    .default = col_double()
  ),
  show_col_types = FALSE
) %>%
  mutate(
    site = factor(site, levels = names(site_cols)),
    plot_date = date
  ) %>%
  filter(date >= focus_start, date <= focus_end)

if (nrow(flow_plot) == 0L) {
  stop("The processed streamflow climatology has no focal-period records.")
}

# ---- Extend the existing WQ table through July 2022 -------------------------
# The main two-site script intentionally ends at 2021. Only the small 2022
# extension needed for this figure is downloaded here and cached locally.
pull_or_default <- function(x, column, default = NA_character_) {
  if (column %in% names(x)) x[[column]] else rep(default, nrow(x))
}

download_2022_wq <- function() {
  raw <- as_tibble(
    read_waterdata_samples(
      monitoringLocationIdentifier = "USGS-10336610",
      activityStartDateLower = "2022-01-01",
      activityStartDateUpper = as.character(focus_end),
      dataProfile = "fullphyschem"
    )
  )

  tibble(
    sample_date = as.Date(pull_or_default(raw, "Activity_StartDate")),
    sample_datetime_utc = as.POSIXct(
      pull_or_default(raw, "Activity_StartDateTime"), tz = "UTC"
    ),
    requested_variable = NA_character_,
    usgs_parameter_code = str_pad(
      as.character(pull_or_default(raw, "USGSpcode")),
      width = 5, side = "left", pad = "0"
    ),
    sample_fraction = as.character(
      pull_or_default(raw, "Result_SampleFraction")
    ),
    result_value = suppressWarnings(as.numeric(
      pull_or_default(raw, "Result_Measure")
    )),
    result_units = as.character(pull_or_default(raw, "Result_MeasureUnit")),
    result_qualifier = as.character(
      pull_or_default(raw, "Result_MeasureQualifierCode")
    ),
    detection_condition = as.character(
      pull_or_default(raw, "Result_ResultDetectionCondition")
    )
  ) %>%
    mutate(
      result_units = if_else(
        str_to_lower(result_units) == "mg/l", "mg/L", result_units
      ),
      nondetect_flag = str_detect(
        str_to_lower(coalesce(detection_condition, "")),
        "not detected|non-detect|below|less than"
      ) | str_detect(coalesce(result_qualifier, ""), "<")
    ) %>%
    filter(usgs_parameter_code %in% nutrient_spec$usgs_parameter_code)
}

wq_2022 <- tryCatch(
  {
    message("Downloading Upper Truckee nutrient results for 2022...")
    downloaded <- download_2022_wq()
    write_csv(downloaded, wq_2022_cache, na = "")
    downloaded
  },
  error = function(e) {
    if (!file.exists(wq_2022_cache)) stop(e)
    warning(
      "The 2022 USGS sample download failed; using cached plot data: ",
      conditionMessage(e), call. = FALSE
    )
    read_csv(
      wq_2022_cache,
      col_types = cols(
        sample_date = col_date(),
        sample_datetime_utc = col_datetime(),
        usgs_parameter_code = col_character(),
        nondetect_flag = col_logical(),
        .default = col_character()
      ),
      show_col_types = FALSE
    ) %>%
      mutate(result_value = suppressWarnings(as.numeric(result_value)))
  }
)

wq_existing <- read_csv(
  wq_file,
  col_types = cols(
    sample_date = col_date(),
    sample_datetime_utc = col_datetime(),
    usgs_parameter_code = col_character(),
    result_value = col_double(),
    nondetect_flag = col_logical(),
    .default = col_character()
  ),
  show_col_types = FALSE
) %>%
  select(
    sample_date, sample_datetime_utc, requested_variable,
    usgs_parameter_code, sample_fraction, result_value, result_units,
    result_qualifier, detection_condition, nondetect_flag
  )

nutrient_wq <- bind_rows(wq_existing, wq_2022) %>%
  mutate(
    usgs_parameter_code = str_pad(
      usgs_parameter_code, width = 5, side = "left", pad = "0"
    ),
    result_units = if_else(
      str_to_lower(result_units) == "mg/l", "mg/L", result_units
    )
  ) %>%
  inner_join(
    nutrient_spec %>%
      select(-requested_variable) %>%
      mutate(plot_label = as.character(plot_label)),
    by = c(
      "usgs_parameter_code", "sample_fraction", "result_units"
    )
  ) %>%
  filter(
    sample_date >= baseline_start,
    sample_date <= focus_end,
    !coalesce(nondetect_flag, FALSE),
    is.finite(result_value)
  ) %>%
  distinct(
    sample_date, sample_datetime_utc, usgs_parameter_code,
    sample_fraction, result_value, result_units, .keep_all = TRUE
  )

# ---- Pre-fire monthly climatology and focal discrete observations -----------
historical_year_month <- nutrient_wq %>%
  filter(sample_date >= baseline_start, sample_date <= baseline_end) %>%
  mutate(year = year(sample_date), month = month(sample_date)) %>%
  group_by(plot_label, facet_label, year, month, sample_date) %>%
  summarise(date_mean = mean(result_value), .groups = "drop_last") %>%
  summarise(year_month_mean = mean(date_mean), .groups = "drop")

monthly_climatology <- historical_year_month %>%
  group_by(plot_label, facet_label, month) %>%
  summarise(
    n_years = n(),
    climatology_mean = mean(year_month_mean),
    climatology_sd = sd(year_month_mean),
    climatology_se = climatology_sd / sqrt(n_years),
    climatology_ci = if_else(
      n_years > 1L,
      qt(0.975, df = n_years - 1L) * climatology_se,
      NA_real_
    ),
    climatology_lo = pmax(climatology_mean - climatology_ci, 0),
    climatology_hi = climatology_mean + climatology_ci,
    .groups = "drop"
  ) %>%
  filter(n_years >= 3L)

focus_months <- tibble(
  month_start = seq(
    floor_date(focus_start, "month"),
    floor_date(focus_end, "month"),
    by = "month"
  )
) %>%
  mutate(
    plot_date = month_start + days(14),
    month = month(month_start)
  )

nutrient_climatology <- crossing(
  nutrient_spec %>%
    transmute(
      plot_label = as.character(plot_label),
      facet_label
    ),
  focus_months
) %>%
  left_join(
    monthly_climatology,
    by = c("plot_label", "facet_label", "month")
  )

focus_nutrients <- nutrient_wq %>%
  filter(sample_date >= focus_start, sample_date <= focus_end) %>%
  group_by(plot_label, facet_label, sample_date) %>%
  summarise(
    concentration_mg_l = mean(result_value),
    n_results = n(),
    .groups = "drop"
  )

availability <- nutrient_spec %>%
  transmute(
    plot_label = as.character(plot_label),
    facet_label,
    usgs_parameter_code,
    sample_fraction,
    result_units
  ) %>%
  left_join(
    focus_nutrients %>%
      group_by(plot_label, facet_label) %>%
      summarise(
        focal_dates = n_distinct(sample_date),
        focal_first = min(sample_date),
        focal_last = max(sample_date),
        .groups = "drop"
      ),
    by = c("plot_label", "facet_label")
  ) %>%
  left_join(
    monthly_climatology %>%
      group_by(plot_label, facet_label) %>%
      summarise(
        climatology_months = n_distinct(month),
        minimum_years_per_month = min(n_years),
        .groups = "drop"
      ),
    by = c("plot_label", "facet_label")
  ) %>%
  mutate(
    across(c(focal_dates, climatology_months), ~replace_na(.x, 0L)),
    usable = focal_dates >= 2L & climatology_months >= 6L,
    omission_reason = case_when(
      usable ~ NA_character_,
      focal_dates < 2L ~ "Fewer than two detected focal-period sample dates",
      climatology_months < 6L ~ "Pre-fire climatology available for fewer than six months",
      TRUE ~ "Unavailable"
    )
  )

write_csv(availability, availability_file, na = "")

usable_labels <- availability %>%
  filter(usable) %>%
  pull(plot_label)

if (length(usable_labels) == 0L) {
  stop("No nutrient series met the focal-data and climatology criteria.")
}

label_levels <- nutrient_spec %>%
  filter(as.character(plot_label) %in% usable_labels) %>%
  pull(facet_label)

focus_nutrients <- focus_nutrients %>%
  filter(plot_label %in% usable_labels) %>%
  mutate(facet_label = factor(facet_label, levels = label_levels))

nutrient_climatology <- nutrient_climatology %>%
  filter(plot_label %in% usable_labels) %>%
  mutate(facet_label = factor(facet_label, levels = label_levels))

plot_export <- bind_rows(
  focus_nutrients %>%
    transmute(
      series = "2021-2022 observations", plot_label, facet_label,
      date = sample_date, value = concentration_mg_l,
      lower = NA_real_, upper = NA_real_, n = n_results
    ),
  nutrient_climatology %>%
    transmute(
      series = "2000-2020 monthly climatology", plot_label, facet_label,
      date = plot_date, value = climatology_mean,
      lower = climatology_lo, upper = climatology_hi, n = n_years
    )
)
write_csv(plot_export, plot_data_file, na = "")

# ---- Shared aesthetics -------------------------------------------------------
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
    axis.title = element_text(size = 8.5),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

fire_layers <- list(
  annotate(
    "rect", xmin = fire_start, xmax = fire_end,
    ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.07
  ),
  geom_vline(
    xintercept = fire_start, colour = "firebrick",
    linetype = "dashed", linewidth = LO_LW_MID
  )
)

series_scale <- scale_linetype_manual(
  values = c("2021-2022" = "solid", "Climatology" = "solid"),
  breaks = c("2021-2022", "Climatology"),
  name = NULL
)

# ---- Panel a: discharge ------------------------------------------------------
p_flow <- ggplot() +
  fire_layers +
  annotate(
    "text",
    x = fire_start + (fire_end - fire_start) / 2,
    y = Inf,
    label = "Caldor Fire",
    colour = "firebrick",
    family = LO_FONT,
    size = 3,
    vjust = 1.25
  ) +
  geom_ribbon(
    data = flow_plot,
    aes(
      x = plot_date, ymin = clim_lo_cms, ymax = clim_hi_cms,
      fill = site, group = site
    ),
    alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = flow_plot,
    aes(
      x = plot_date, y = clim_center_cms, colour = site,
      linetype = "Climatology", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = flow_plot,
    aes(
      x = date, y = flow_cms, colour = site,
      linetype = "2021-2022", group = site
    ),
    linewidth = LO_LW_THICK, alpha = 0.95, na.rm = TRUE
  ) +
  scale_colour_manual(values = site_cols, name = NULL) +
  scale_fill_manual(values = site_cols, guide = "none") +
  series_scale +
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
    date_breaks = "1 month", date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(x = NULL, y = expression("Discharge (m"^3~"s"^{-1}*")")) +
  stream_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(1.5, 2, 0.5, 2, "mm")
  )

# ---- Panel b: discrete nutrient series --------------------------------------
p_nutrients <- ggplot() +
  fire_layers +
  geom_ribbon(
    data = nutrient_climatology,
    aes(
      x = plot_date, ymin = climatology_lo, ymax = climatology_hi,
      group = facet_label
    ),
    fill = upper_blue, alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = nutrient_climatology,
    aes(
      x = plot_date, y = climatology_mean,
      linetype = "Climatology", group = facet_label
    ),
    colour = upper_blue, linewidth = LO_LW_THICK,
    alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = focus_nutrients,
    aes(
      x = sample_date, y = concentration_mg_l,
      linetype = "2021-2022", group = facet_label
    ),
    colour = upper_blue, linewidth = LO_LW_THICK,
    alpha = 0.95, na.rm = TRUE
  ) +
  geom_point(
    data = focus_nutrients,
    aes(x = sample_date, y = concentration_mg_l),
    colour = upper_blue, fill = "white", shape = 21,
    stroke = 0.45, size = 1.55, alpha = 0.95, na.rm = TRUE
  ) +
  facet_wrap(vars(facet_label), ncol = 3, scales = "free_y") +
  series_scale +
  guides(
    linetype = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(
        colour = "grey30", alpha = c(0.95, 0.34), linewidth = 0.9
      )
    )
  ) +
  scale_x_date(
    limits = c(focus_start, focus_end),
    date_breaks = "3 months", date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(x = NULL, y = "Concentration") +
  stream_theme +
  theme(
    strip.text = element_text(size = 8, face = "plain"),
    plot.margin = margin(0.5, 2, 1.5, 2, "mm")
  )

combined <- (p_flow / p_nutrients) +
  plot_layout(heights = c(0.9, 1.55), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        family = LO_FONT, face = "bold", size = 9
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center"
  )

save_lo_fig(
  combined,
  figure_file,
  width_type = "double",
  height_cm = 15.2
)

# ---- Console summary ---------------------------------------------------------
cat("\n=== Combined streamflow and nutrient time-series figure ===\n")
cat("Focal interval:", as.character(focus_start), "to", as.character(focus_end), "\n")
cat("Pre-fire nutrient climatology:", as.character(baseline_start), "to", as.character(baseline_end), "\n")
cat("Nutrients plotted:\n")
print(
  availability %>%
    filter(usable) %>%
    select(plot_label, usgs_parameter_code, focal_dates, climatology_months),
  n = Inf
)
cat("Nutrients omitted:\n")
print(
  availability %>%
    filter(!usable) %>%
    select(plot_label, usgs_parameter_code, omission_reason),
  n = Inf
)
cat("Outputs:\n")
cat("  ", plot_data_file, "\n", sep = "")
cat("  ", availability_file, "\n", sep = "")
cat("  ", figure_file, "\n", sep = "")
