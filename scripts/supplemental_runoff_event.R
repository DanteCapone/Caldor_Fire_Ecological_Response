# =============================================================================
# supplemental_runoff_event.R
#
# Recreate the supplemental tributary discharge + nutrient time-series figure
# and quantify the 24-31 October 2021 post-Caldor runoff pulse.
#
# Event nutrient loads are estimates for the Upper Truckee River only. Daily
# concentrations are linearly interpolated between the 24 and 26 October USGS
# samples and carried forward after 26 October. "Excess" is relative to the
# median concentration on historical October-November sample dates (2000-2020)
# for the identical USGS parameter, fraction, and units.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
} else {
  project_root <- normalizePath(".", winslash = "/")
}
setwd(project_root)
source("scripts/figure_aesthetics.R")

processed_dir <- file.path("data", "processed")
supplement_dir <- Sys.getenv(
  "CALDOR_RUNOFF_SUPPLEMENT_DIR",
  unset = file.path("figures", "supplemental")
)
figure_only <- identical(Sys.getenv("CALDOR_RUNOFF_FIGURE_ONLY"), "1")
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

flow_file <- file.path(
  processed_dir, "usgs_tahoe_tributary_streamflow_daily_climatology.csv"
)
site_volume_file <- file.path(
  processed_dir, "usgs_tahoe_inflow_outflow_daily_site_volume.csv"
)
nutrient_plot_file <- file.path(
  processed_dir, "upper_truckee_nutrient_timeseries_2021_2022.csv"
)
availability_file <- file.path(
  processed_dir, "upper_truckee_nutrient_timeseries_availability.csv"
)
water_quality_file <- file.path(
  processed_dir, "upper_truckee_water_quality.csv"
)
hypsography_file <- file.path(
  processed_dir, "ctd", "lake_tahoe_hypsography_usgs.csv"
)

required_files <- c(
  flow_file, site_volume_file, nutrient_plot_file, availability_file,
  water_quality_file, hypsography_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input(s):\n", paste(missing_files, collapse = "\n"))
}

focus_start <- as.Date("2021-07-01")
focus_end <- as.Date("2022-07-31")
fire_start <- as.Date("2021-08-14")
fire_end <- as.Date("2021-10-21")
event_start <- as.Date("2021-10-24")
event_end <- as.Date("2021-10-31")
event_dates <- seq(event_start, event_end, by = "day")
historical_start <- as.Date("2000-01-01")
historical_end <- as.Date("2020-12-31")

runoff_month_labels <- function(x) {
  labels <- format(x, "%b")
  if (length(x) > 0L) {
    labels[1] <- paste0(labels[1], "\n", year(x[1]))
  }
  labels
}

site_colours <- c(
  "Upper Truckee River" = "#1565C0",
  "Blackwood Creek" = "#2E7D32"
)
upper_blue <- unname(site_colours[["Upper Truckee River"]])

# ---- Discharge and whole-lake event volume ---------------------------------
flow <- read_csv(
  flow_file,
  col_types = cols(
    site_no = col_character(), site = col_character(), date = col_date(),
    qualifier = col_character(), .default = col_double()
  ),
  show_col_types = FALSE
) %>%
  filter(date >= focus_start, date <= focus_end) %>%
  mutate(site = factor(site, levels = names(site_colours)))

site_volume <- read_csv(
  site_volume_file,
  col_types = cols(
    site_no = col_character(), date = col_date(), .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
  mutate(
    daily_volume_m3 = as.numeric(daily_volume_m3),
    flow_role = as.character(flow_role)
  ) %>%
  filter(date >= event_start, date <= event_end, is.finite(daily_volume_m3))

hypsography <- read_csv(hypsography_file, show_col_types = FALSE) %>%
  arrange(depth_m)
lake_volume_m3 <- sum(
  diff(hypsography$depth_m) *
    (head(hypsography$area_m2, -1L) + tail(hypsography$area_m2, -1L)) / 2
)

site_event_volume <- site_volume %>%
  group_by(site_no, site, flow_role) %>%
  summarise(
    n_days = n_distinct(date),
    event_volume_m3 = sum(daily_volume_m3),
    .groups = "drop"
  ) %>%
  transmute(
    component = site,
    role = flow_role,
    n_days,
    event_volume_m3
  )

upper_volume <- site_event_volume %>%
  filter(str_detect(component, regex("Upper Truckee", ignore_case = TRUE))) %>%
  summarise(value = sum(event_volume_m3)) %>%
  pull(value)
blackwood_volume <- site_event_volume %>%
  filter(str_detect(component, regex("Blackwood", ignore_case = TRUE))) %>%
  summarise(value = sum(event_volume_m3)) %>%
  pull(value)
outlet_volume <- site_event_volume %>%
  filter(role == "Outlet") %>%
  summarise(value = sum(event_volume_m3)) %>%
  pull(value)
monitored_inflow_volume <- upper_volume + blackwood_volume

event_volume_summary <- bind_rows(
  site_event_volume,
  tibble(
    component = c(
      "Upper Truckee River + Blackwood Creek",
      "Monitored inflow minus Tahoe outlet"
    ),
    role = c("Combined monitored inflow", "Partial monitored net"),
    n_days = length(event_dates),
    event_volume_m3 = c(
      monitored_inflow_volume,
      monitored_inflow_volume - outlet_volume
    )
  )
) %>%
  mutate(
    event_start = event_start,
    event_end = event_end,
    lake_volume_m3 = lake_volume_m3,
    lake_volume_km3 = lake_volume_m3 / 1e9,
    percent_of_lake_volume = 100 * event_volume_m3 / lake_volume_m3,
    volume_equation = "sum(daily mean discharge * 86400 s/day)",
    lake_volume_method = paste(
      "Trapezoidal integration of area by depth from",
      "USGS Lake Tahoe bathymetry-derived hypsography"
    )
  )

volume_summary_file <- file.path(
  processed_dir, "upper_truckee_blackwood_event_volume_2021_10_24_31.csv"
)
if (!figure_only) write_csv(event_volume_summary, volume_summary_file, na = "")

# ---- Upper Truckee event nutrient flux/load estimates ----------------------
availability <- read_csv(
  availability_file,
  col_types = cols(usgs_parameter_code = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  filter(usable) %>%
  mutate(
    usgs_parameter_code = str_pad(
      usgs_parameter_code, width = 5, side = "left", pad = "0"
    )
  )

wq <- read_csv(
  water_quality_file,
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
  mutate(
    usgs_parameter_code = str_pad(
      usgs_parameter_code, width = 5, side = "left", pad = "0"
    ),
    result_units = if_else(
      str_to_lower(result_units) == "mg/l", "mg/L", result_units
    )
  ) %>%
  inner_join(
    availability %>%
      select(
        plot_label, usgs_parameter_code, sample_fraction, result_units
      ),
    by = c("usgs_parameter_code", "sample_fraction", "result_units")
  ) %>%
  filter(!coalesce(nondetect_flag, FALSE), is.finite(result_value))

historical_baseline <- wq %>%
  filter(
    sample_date >= historical_start,
    sample_date <= historical_end,
    month(sample_date) %in% c(10L, 11L)
  ) %>%
  group_by(
    plot_label, usgs_parameter_code, sample_fraction, result_units,
    sample_date
  ) %>%
  summarise(date_concentration_mg_l = mean(result_value), .groups = "drop") %>%
  group_by(plot_label, usgs_parameter_code, sample_fraction, result_units) %>%
  summarise(
    baseline_median_mg_l = median(date_concentration_mg_l),
    baseline_mean_mg_l = mean(date_concentration_mg_l),
    baseline_n_dates = n(),
    baseline_n_years = n_distinct(year(sample_date)),
    .groups = "drop"
  )

# Samples beginning 21 October provide the nearest bracketing chemistry; for
# the 24-31 October load window, all usable analytes have 24 and/or 26 Oct data.
event_chemistry <- wq %>%
  filter(
    sample_date >= as.Date("2021-10-21"),
    sample_date <= event_end
  ) %>%
  group_by(
    plot_label, usgs_parameter_code, sample_fraction, result_units,
    sample_date
  ) %>%
  summarise(concentration_mg_l = mean(result_value), .groups = "drop")

interpolate_concentration <- function(dates, concentrations, targets) {
  if (length(unique(dates)) == 1L) {
    return(rep(mean(concentrations), length(targets)))
  }
  approx(
    x = as.numeric(dates), y = concentrations,
    xout = as.numeric(targets), method = "linear", rule = 2, ties = mean
  )$y
}

event_daily_concentration <- event_chemistry %>%
  group_by(plot_label, usgs_parameter_code, sample_fraction, result_units) %>%
  group_modify(~ tibble(
    date = event_dates,
    event_concentration_mg_l = interpolate_concentration(
      .x$sample_date, .x$concentration_mg_l, event_dates
    ),
    chemistry_first_date = min(.x$sample_date),
    chemistry_last_date = max(.x$sample_date),
    chemistry_n_dates = n_distinct(.x$sample_date)
  )) %>%
  ungroup()

upper_event_flow <- site_volume %>%
  filter(str_detect(site, regex("Upper Truckee", ignore_case = TRUE))) %>%
  transmute(date, discharge_cms = daily_volume_m3 / 86400)

event_flux_daily <- event_daily_concentration %>%
  left_join(historical_baseline, by = c(
    "plot_label", "usgs_parameter_code", "sample_fraction", "result_units"
  )) %>%
  left_join(upper_event_flow, by = "date") %>%
  mutate(
    daily_volume_m3 = discharge_cms * 86400,
    event_load_kg_day = event_concentration_mg_l * discharge_cms * 86.4,
    baseline_load_kg_day = baseline_median_mg_l * discharge_cms * 86.4,
    net_added_load_kg_day =
      (event_concentration_mg_l - baseline_median_mg_l) * discharge_cms * 86.4,
    positive_excess_load_kg_day = pmax(net_added_load_kg_day, 0)
  )

event_flux_daily_file <- file.path(
  processed_dir, "upper_truckee_event_nutrient_flux_daily_2021_10_24_31.csv"
)
if (!figure_only) write_csv(event_flux_daily, event_flux_daily_file, na = "")

event_flux_summary <- event_flux_daily %>%
  group_by(
    plot_label, usgs_parameter_code, sample_fraction, result_units,
    baseline_median_mg_l, baseline_mean_mg_l, baseline_n_dates,
    baseline_n_years, chemistry_first_date, chemistry_last_date,
    chemistry_n_dates
  ) %>%
  summarise(
    event_start = min(date),
    event_end = max(date),
    n_event_days = n(),
    upper_truckee_event_volume_m3 = sum(daily_volume_m3),
    event_mean_concentration_mg_l = mean(event_concentration_mg_l),
    event_total_load_kg = sum(event_load_kg_day),
    historical_concentration_load_kg = sum(baseline_load_kg_day),
    net_added_load_kg = sum(net_added_load_kg_day),
    positive_excess_load_kg = sum(positive_excess_load_kg_day),
    mean_event_flux_kg_day = mean(event_load_kg_day),
    mean_net_added_flux_kg_day = mean(net_added_load_kg_day),
    .groups = "drop"
  ) %>%
  mutate(
    concentration_method = paste(
      "Linear interpolation between event sample dates;",
      "nearest-sample carry before/after sampled dates"
    ),
    baseline_method = paste(
      "Median of historical October-November sample-date means,",
      "2000-2020; identical USGS parameter/fraction/units"
    ),
    scope_note = paste(
      "Upper Truckee River only; approximate event load from daily mean",
      "discharge and sparse discrete chemistry; not extrapolated to Blackwood Creek"
    )
  )

event_flux_summary_file <- file.path(
  processed_dir, "upper_truckee_event_nutrient_flux_summary_2021_10_24_31.csv"
)
if (!figure_only) write_csv(event_flux_summary, event_flux_summary_file, na = "")

# ---- Supplemental panel figure ---------------------------------------------
nutrient_plot <- read_csv(
  nutrient_plot_file,
  col_types = cols(date = col_date(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  filter(plot_label %in% availability$plot_label) %>%
  mutate(
    plot_label = factor(plot_label, levels = availability$plot_label),
    series = factor(
      series,
      levels = c("2021-2022 observations", "2000-2020 monthly climatology")
    ),
    # USGS chemistry is stored in mg/L; convert both observations and the
    # climatological summaries to ug/L for the manuscript-facing display.
    across(any_of(c("value", "lower", "upper")), ~ .x * 1000)
  )

observations <- nutrient_plot %>% filter(series == "2021-2022 observations")
climatology <- nutrient_plot %>%
  filter(series == "2000-2020 monthly climatology")

nutrient_facet_labels <- c(
  "NO3 + NO2 as N" = "NO[3]^'-' + NO[2]^'-'",
  "NH4/NH3 as N" = "NH[4]^'+'",
  "Total nitrogen" = "Total~N",
  "TKN/organic N" = "TKN",
  "Orthophosphate as P" = "Orthophosphate~as~P",
  "Total phosphorus" = "Total~P"
)

common_theme <- lo_theme(base_size = 8, family = LO_FONT) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.5),
    axis.text = element_text(size = 7.5),
    axis.title = element_text(size = 8),
    plot.title = element_blank()
  )

event_layers <- list(
  annotate(
    "rect", xmin = event_start, xmax = event_end + 1,
    ymin = -Inf, ymax = Inf, fill = "#4FC3F7", alpha = 0.13
  ),
  geom_vline(
    xintercept = event_start, colour = "#0277BD",
    linetype = "dotted", linewidth = LO_LW_MID
  )
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

p_discharge <- ggplot() +
  fire_layers + event_layers +
  geom_ribbon(
    data = flow,
    aes(date, ymin = clim_lo_cms, ymax = clim_hi_cms, fill = site, group = site),
    alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = flow,
    aes(date, clim_center_cms, colour = site, group = site),
    linewidth = LO_LW_THICK, alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = flow,
    aes(date, flow_cms, colour = site, group = site),
    linewidth = LO_LW_THICK, alpha = 0.95, na.rm = TRUE
  ) +
  scale_colour_manual(values = site_colours, name = NULL) +
  scale_fill_manual(values = site_colours, guide = "none") +
  guides(
    colour = guide_legend(order = 1, nrow = 1)
  ) +
  scale_x_date(
    limits = c(focus_start, focus_end),
    breaks = seq(focus_start, focus_end, by = "2 months"),
    labels = runoff_month_labels, expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(x = NULL, y = expression("Discharge (m"^3~"s"^{-1}*")")) +
  common_theme +
  theme(
    axis.text.x = element_text(size = 6.5),
    plot.margin = margin(1.5, 2, 0.5, 2, "mm")
  )

p_nutrients <- ggplot() +
  fire_layers + event_layers +
  geom_ribbon(
    data = climatology,
    aes(date, ymin = lower, ymax = upper, group = plot_label),
    fill = upper_blue, alpha = 0.09, colour = NA, na.rm = TRUE
  ) +
  geom_line(
    data = climatology,
    aes(date, value, group = plot_label),
    colour = upper_blue, linewidth = LO_LW_THICK, alpha = 0.34, na.rm = TRUE
  ) +
  geom_line(
    data = observations,
    aes(date, value, group = plot_label),
    colour = upper_blue, linewidth = LO_LW_THICK, alpha = 0.95, na.rm = TRUE
  ) +
  facet_wrap(
    vars(plot_label), ncol = 3, scales = "free_y",
    labeller = labeller(
      plot_label = as_labeller(nutrient_facet_labels, label_parsed)
    )
  ) +
  scale_x_date(
    limits = c(focus_start, focus_end),
    breaks = seq(focus_start, focus_end, by = "4 months"),
    labels = runoff_month_labels, expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(x = NULL, y = expression("Concentration ("*mu*"g L"^{-1}*")")) +
  common_theme +
  theme(
    strip.text = element_text(size = 7.5, face = "plain"),
    plot.margin = margin(0.5, 2, 1.5, 2, "mm")
  )

combined_plot <- (p_discharge / p_nutrients) +
  plot_layout(heights = c(0.85, 1.55), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = LO_FONT, face = "bold", size = 8.5)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = unit(0.2, "mm")
  )

figure_file <- file.path(
  supplement_dir, "supplemental_runoff_discharge_nutrients_2021_2022.png"
)
save_lo_fig(combined_plot, figure_file, width_type = "double", height_cm = 15.2)
ggsave(
  sub("\\.png$", ".pdf", figure_file), combined_plot,
  width = 12.7, height = 15.2, units = "cm", device = cairo_pdf,
  bg = "white"
)

cat("\n=== 24-31 October 2021 runoff summary ===\n")
print(
  event_volume_summary %>%
    select(component, role, event_volume_m3, percent_of_lake_volume),
  n = Inf
)
cat("\nUpper Truckee nutrient event estimates (kg):\n")
print(
  event_flux_summary %>%
    select(
      plot_label, event_total_load_kg, net_added_load_kg,
      positive_excess_load_kg, chemistry_n_dates
    ),
  n = Inf
)
cat("\nOutputs:\n")
cat(paste0("  ", c(
  figure_file, volume_summary_file, event_flux_daily_file,
  event_flux_summary_file
)), sep = "\n")
cat("\n")
