# =============================================================================
# Pickle-based CTD profiles, Schmidt stability, and Lake Number
#
# Inputs:
#   data/lake_environmental_data/ctd/01_ctd/{ltp,mltp}_ctd.pkl
#   data/lake_environmental_data/00_met/TB4.csv
#
# Focal profile window: June 2021-April 2022
# Stability window: January 2021-December 2022
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
  library(ggh4x)
  library(scales)
  library(patchwork)
  library(ragg)
  library(rLakeAnalyzer)
})

proj_root <- normalizePath(".", winslash = "/")
source("scripts/main_analysis/shared_aesthetics.R")

extractor <- file.path(
  proj_root, "scripts", "main_analysis", "ctd_pickle_extract.py"
)
ctd_pickle_dir <- file.path(
  proj_root, "data", "lake_environmental_data", "ctd", "01_ctd"
)
met_file <- file.path(
  proj_root, "data", "lake_environmental_data", "00_met", "TB4.csv"
)
proc_dir <- file.path(proj_root, "data", "processed", "ctd")
out_dir <- Sys.getenv(
  "CALDOR_CTD_OUT_DIR",
  unset = file.path(proj_root, "figures", "supplemental", "ctd_profiles")
)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

profile_summary_file <- file.path(
  proc_dir, "ctd_profile_summary_jun2021_apr2022_from_pkl.csv"
)
profile_availability_file <- file.path(
  proc_dir, "ctd_profile_availability_jun2021_apr2022_from_pkl.csv"
)
temperature_file <- file.path(
  proc_dir, "mltp_temperature_profiles_1m_from_pkl.csv"
)
hypsography_file <- file.path(proc_dir, "lake_tahoe_hypsography_usgs.csv")
metrics_file <- file.path(
  proc_dir, "mltp_stability_metrics_from_pkl_tb4.csv"
)
metric_plot_file <- file.path(
  proc_dir, "mltp_stability_metrics_2019_2023_plot_data.csv"
)
methods_file <- file.path(proc_dir, "mltp_stability_methods.csv")

profile_mltp_file <- file.path(
  out_dir, "ctd_profiles_mltp_jun2021_apr2022.png"
)
profile_ltp_file <- file.path(
  out_dir, "ctd_profiles_ltp_jun2021_apr2022.png"
)
schmidt_file <- file.path(
  out_dir, "ctd_schmidt_stability_mltp_2019_2023.png"
)
schmidt_pdf <- file.path(
  out_dir, "ctd_schmidt_stability_mltp_2019_2023.pdf"
)
lake_number_file <- file.path(
  out_dir, "ctd_lake_number_mltp_2019_2023.png"
)
combined_stability_file <- file.path(
  out_dir, "ctd_stability_metrics_mltp_2019_2023.png"
)

required_inputs <- c(
  extractor,
  file.path(ctd_pickle_dir, "ltp_ctd.pkl"),
  file.path(ctd_pickle_dir, "mltp_ctd.pkl"),
  met_file,
  hypsography_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required CTD inputs: ", paste(missing_inputs, collapse = ", "))
}

s7_figure_only <- identical(Sys.getenv("CALDOR_S7_FIGURE_ONLY"), "1")

# ---- Extract compact products directly from the pickle matrices ------------
python_candidates <- unique(c(
  Sys.getenv("CTD_PYTHON", unset = NA_character_),
  file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "geomap", "python.exe"),
  file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "geomap2", "python.exe"),
  file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "holoview", "python.exe"),
  Sys.which("python")
))
python_candidates <- python_candidates[
  !is.na(python_candidates) & nzchar(python_candidates) & file.exists(python_candidates)
]

python_ok <- vapply(python_candidates, function(candidate) {
  status <- suppressWarnings(system2(
    candidate,
    c("-c", shQuote("import pandas, numpy")),
    stdout = FALSE, stderr = FALSE
  ))
  identical(status, 0L)
}, logical(1))
if (!any(python_ok)) {
  stop(
    "No Python environment with pandas and numpy was found. ",
    "Set CTD_PYTHON to a compatible python.exe."
  )
}
python <- python_candidates[which(python_ok)[1]]
if (!s7_figure_only) {
  message("Extracting pickle CTD products with: ", python)
  extract_status <- system2(
    python,
    c(shQuote(extractor), "--project-root", shQuote(proj_root))
  )
  if (!identical(extract_status, 0L)) {
    stop("The CTD pickle extractor failed with status ", extract_status)
  }
} else {
  message("CALDOR_S7_FIGURE_ONLY=1; using existing extracted CTD products.")
}

# ---- June 2021-April 2022 profile figures ----------------------------------
profile_start <- as.Date("2021-06-01")
profile_end <- as.Date("2022-04-30")
stability_start <- as.Date("2019-01-01")
stability_end <- as.Date("2023-12-31")
fire_start <- as.Date("2021-08-14")
fire_end <- as.Date("2021-10-21")
month_levels <- c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
                  "Jan", "Feb", "Mar", "Apr")

variable_levels <- c(
  "Temperature", "Potential_Density", "Conductivity", "Dissolved_Oxygen",
  "Chl_Fluorescence", "Turbidity", "PAR"
)
variable_labels <- c(
  "Temperature" = "Temperature\n(deg C)",
  "Potential_Density" = "Potential density\n(kg m^-3)",
  "Conductivity" = "Conductivity\n(microS cm^-1)",
  "Dissolved_Oxygen" = "Dissolved oxygen\n(mg L^-1)",
  "Chl_Fluorescence" = "Chl fluorescence\n(mg L^-1)",
  "Turbidity" = "Turbidity\n(FTU)",
  "PAR" = "PAR\n(micromol m^-2 s^-1)"
)
obs_cols <- c(
  "Temperature" = "#D73027",
  "Potential_Density" = "#6A3D9A",
  "Conductivity" = "#2166AC",
  "Dissolved_Oxygen" = "#0072B2",
  "Chl_Fluorescence" = "#009E73",
  "Turbidity" = "#B8860B",
  "PAR" = "#E66101"
)
clim_cols <- c(
  "Temperature" = "#F4A6A1",
  "Potential_Density" = "#CAB2D6",
  "Conductivity" = "#9ECAE1",
  "Dissolved_Oxygen" = "#A6CEE3",
  "Chl_Fluorescence" = "#99D8C9",
  "Turbidity" = "#FDDC8A",
  "PAR" = "#FDBB84"
)

profile_summary <- read_csv(profile_summary_file, show_col_types = FALSE) %>%
  mutate(
    variable = factor(variable, levels = variable_levels),
    month_label = factor(month_label, levels = month_levels)
  )
profile_availability <- read_csv(
  profile_availability_file, show_col_types = FALSE
) %>%
  mutate(
    variable = factor(variable, levels = variable_levels),
    month_label = factor(month_label, levels = month_levels)
  )

station_alias <- c("Mid-lake" = "MLTP", "Index" = "LTP")
station_depth_cap <- c("Mid-lake" = 500, "Index" = 160)

make_profile_plot <- function(station_name) {
  station_data <- profile_summary %>% filter(Station_ID == station_name)
  missing_focal <- profile_availability %>%
    filter(Station_ID == station_name, focal_casts == 0L) %>%
    mutate(
      depth_m = unname(station_depth_cap[[station_name]]) * 0.53,
      note = "No focal\ndata"
    )
  depth_cap <- unname(station_depth_cap[[station_name]])
  depth_breaks <- if (station_name == "Mid-lake") {
    seq(0, 500, 100)
  } else {
    c(0, 30, 60, 90, 120, 150)
  }

  ggplot() +
    geom_ribbon(
      data = station_data %>% filter(is.finite(clim_mean)),
      aes(y = depth_m, xmin = clim_lo, xmax = clim_hi, fill = variable,
          group = interaction(variable, month)),
      orientation = "y", alpha = 0.25, colour = NA, na.rm = TRUE
    ) +
    geom_path(
      data = station_data %>% filter(is.finite(clim_mean)),
      aes(x = clim_mean, y = depth_m, colour = variable,
          group = interaction(variable, month),
          linetype = "Historical climatology"),
      linewidth = LO_LW_MID, alpha = 0.55, na.rm = TRUE
    ) +
    geom_errorbar(
      data = station_data %>% filter(is.finite(obs_mean)),
      aes(y = depth_m, xmin = obs_lo, xmax = obs_hi, colour = variable),
      orientation = "y", height = 3, linewidth = 0.3 * LO_FIG_SCALE,
      alpha = 0.65, na.rm = TRUE
    ) +
    geom_path(
      data = station_data %>% filter(is.finite(obs_mean)),
      aes(x = obs_mean, y = depth_m, colour = variable,
          group = interaction(variable, month),
          linetype = "June 2021-April 2022"),
      linewidth = LO_LW_THICK, na.rm = TRUE
    ) +
    geom_point(
      data = station_data %>% filter(is.finite(obs_mean)),
      aes(x = obs_mean, y = depth_m, colour = variable),
      size = 0.95 * LO_FIG_SCALE, stroke = 0.2 * LO_FIG_SCALE,
      na.rm = TRUE
    ) +
    geom_text(
      data = missing_focal,
      aes(x = Inf, y = depth_m, label = note),
      hjust = 1.08, vjust = 0.5, colour = "grey45",
      family = LO_FONT, size = 2.0, lineheight = 0.9
    ) +
    ggh4x::facet_grid2(
      rows = vars(variable), cols = vars(month_label),
      scales = "free_x", independent = "x", switch = "y", drop = FALSE,
      labeller = labeller(variable = as_labeller(variable_labels))
    ) +
    scale_y_reverse(
      limits = c(depth_cap, 0), breaks = depth_breaks,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_x_continuous(
      labels = label_number(accuracy = 0.01, big.mark = ","),
      expand = expansion(mult = c(0.05, 0.16))
    ) +
    scale_colour_manual(values = obs_cols, guide = "none") +
    scale_fill_manual(values = clim_cols, guide = "none") +
    scale_linetype_manual(
      values = c(
        "June 2021-April 2022" = "solid",
        "Historical climatology" = "dashed"
      ),
      breaks = c("June 2021-April 2022", "Historical climatology"),
      name = NULL
    ) +
    labs(
      title = paste0(station_name, " (", station_alias[[station_name]], ")"),
      x = "CTD value (units shown in row labels)", y = "Depth (m)"
    ) +
    theme_classic(base_size = 7, base_family = LO_FONT) +
    theme(
      text = element_text(family = LO_FONT),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 5.2),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = LO_LW_THIN),
      panel.grid.minor.y = element_line(colour = "grey95", linewidth = LO_LW_THIN),
      strip.background = element_blank(), strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 6.2),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing.x = unit(1.3, "mm"), panel.spacing.y = unit(1.4, "mm"),
      legend.position = "bottom", legend.key.width = unit(7, "mm"),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

p_mltp <- make_profile_plot("Mid-lake")
p_ltp <- make_profile_plot("Index")
if (!s7_figure_only) {
  ggsave(
    profile_mltp_file, p_mltp, width = 43, height = 30, units = "cm",
    dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE
  )
  ggsave(
    profile_ltp_file, p_ltp, width = 43, height = 30, units = "cm",
    dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE
  )
}

# ---- Schmidt stability and Lake Number -------------------------------------
hypsography <- read_csv(hypsography_file, show_col_types = FALSE) %>%
  filter(is.finite(depth_m), is.finite(area_m2)) %>%
  arrange(depth_m)
bth_depth <- hypsography$depth_m
bth_area <- hypsography$area_m2
minimum_stability_depth_m <- 400
wind_height_m <- 10
wind_window_hours <- 6
minimum_wind_speed_m_s <- 0.1

temperature_profiles <- fread(temperature_file, showProgress = FALSE)
temperature_profiles[, dateutc := as.POSIXct(dateutc, tz = "UTC")]
temperature_profiles[, date := as.IDate(date)]
setorder(temperature_profiles, cast_id, depth_m)

calculate_cast_metrics <- function(depth_m, temperature_c) {
  keep <- is.finite(depth_m) & is.finite(temperature_c)
  depth_m <- depth_m[keep]
  temperature_c <- temperature_c[keep]
  ord <- order(depth_m)
  depth_m <- depth_m[ord]
  temperature_c <- temperature_c[ord]
  max_depth <- max(depth_m)
  n_depths <- length(unique(depth_m))
  if (max_depth < minimum_stability_depth_m || n_depths < 50L) {
    return(list(
      maximum_depth_m = max_depth, n_depths = n_depths,
      schmidt_stability_j_m2 = NA_real_, meta_top_m = NA_real_,
      meta_bottom_m = NA_real_, average_epi_density_kg_m3 = NA_real_,
      average_hypo_density_kg_m3 = NA_real_
    ))
  }
  stability <- tryCatch(
    as.numeric(rLakeAnalyzer::schmidt.stability(
      temperature_c, depth_m, bth_area, bth_depth, sal = 0
    )),
    error = function(e) NA_real_
  )
  meta <- tryCatch(
    as.numeric(rLakeAnalyzer::meta.depths(
      temperature_c, depth_m, seasonal = TRUE
    )),
    error = function(e) c(NA_real_, NA_real_)
  )
  if (length(meta) != 2L || any(!is.finite(meta)) || meta[2] <= meta[1]) {
    meta <- c(NA_real_, NA_real_)
  }
  epi_density <- if (all(is.finite(meta))) {
    tryCatch(
      rLakeAnalyzer::layer.density(
        0, meta[1], temperature_c, depth_m, bth_area, bth_depth
      ),
      error = function(e) NA_real_
    )
  } else NA_real_
  hypo_density <- if (all(is.finite(meta)) && meta[2] < max(bth_depth)) {
    tryCatch(
      rLakeAnalyzer::layer.density(
        meta[2], max(bth_depth), temperature_c, depth_m, bth_area, bth_depth
      ),
      error = function(e) NA_real_
    )
  } else NA_real_

  list(
    maximum_depth_m = max_depth, n_depths = n_depths,
    schmidt_stability_j_m2 = stability,
    meta_top_m = meta[1], meta_bottom_m = meta[2],
    average_epi_density_kg_m3 = epi_density,
    average_hypo_density_kg_m3 = hypo_density
  )
}

message("Calculating Schmidt stability and thermal layers from MLTP casts...")
cast_metrics <- temperature_profiles[
  , calculate_cast_metrics(depth_m, temperature_c),
  by = .(cast_id, dateutc, date)
]

message("Reading TB4 wind and calculating centered 6-hour wind forcing...")
wind <- fread(
  met_file, select = c("dateutc", "WindSpeed", "WindDir"),
  showProgress = FALSE
)
wind[, dateutc := as.POSIXct(dateutc, tz = "UTC")]
wind <- wind[
  is.finite(WindSpeed) & WindSpeed >= 0 &
    is.finite(WindDir) & between(WindDir, 0, 360),
  .(dateutc, wind_speed_m_s = WindSpeed, wind_dir_deg = WindDir)
]
setkey(wind, dateutc)
cast_windows <- unique(cast_metrics[, .(cast_id, cast_time = dateutc)])
cast_windows[, `:=`(
  window_start = cast_time - hours(wind_window_hours / 2),
  window_end = cast_time + hours(wind_window_hours / 2)
)]
wind_matches <- wind[
  cast_windows,
  on = .(dateutc >= window_start, dateutc <= window_end),
  allow.cartesian = TRUE, nomatch = 0L,
  .(
    cast_id = i.cast_id, cast_time = i.cast_time,
    wind_speed_observation_m_s = x.wind_speed_m_s,
    wind_dir_observation_deg = x.wind_dir_deg
  )
]
wind_matches[, `:=`(
  wind_u_m_s = -wind_speed_observation_m_s * sin(wind_dir_observation_deg * pi / 180),
  wind_v_m_s = -wind_speed_observation_m_s * cos(wind_dir_observation_deg * pi / 180)
)]
wind_matches <- wind_matches[
  , .(
    n_wind_observations = .N,
    wind_speed_m_s = sqrt(mean(wind_speed_observation_m_s^2)),
    mean_wind_u_m_s = mean(wind_u_m_s),
    mean_wind_v_m_s = mean(wind_v_m_s)
  ),
  by = .(cast_id, cast_time)
]
wind_matches[, wind_dir_deg := (
  atan2(-mean_wind_u_m_s, -mean_wind_v_m_s) * 180 / pi
) %% 360]

cast_metrics <- merge(cast_metrics, wind_matches, by = "cast_id", all.x = TRUE)
cast_metrics[, water_friction_velocity_m_s := fifelse(
  is.finite(wind_speed_m_s) & wind_speed_m_s >= minimum_wind_speed_m_s &
    is.finite(average_epi_density_kg_m3),
  rLakeAnalyzer::uStar(wind_speed_m_s, wind_height_m,
                       average_epi_density_kg_m3),
  NA_real_
)]
cast_metrics[, lake_number := mapply(
  function(u_star, stability, meta_top, meta_bottom, hypo_density) {
    if (!all(is.finite(c(u_star, stability, meta_top, meta_bottom, hypo_density))) ||
        u_star <= 0 || stability <= 0) return(NA_real_)
    tryCatch(
      as.numeric(rLakeAnalyzer::lake.number(
        bth_area, bth_depth, u_star, stability,
        meta_top, meta_bottom, hypo_density
      )),
      error = function(e) NA_real_
    )
  },
  water_friction_velocity_m_s, schmidt_stability_j_m2,
  meta_top_m, meta_bottom_m, average_hypo_density_kg_m3
)]
cast_metrics[, `:=`(
  date = as.Date(date),
  year = lubridate::year(as.Date(date)),
  month = lubridate::month(as.Date(date))
)]
setorder(cast_metrics, date, cast_id)
if (!s7_figure_only) fwrite(cast_metrics, metrics_file, na = "")

if (!s7_figure_only) write_csv(
  tribble(
    ~item, ~value,
    "CTD source", "MLTP pickle temperature matrix, binned to 1 m",
    "Schmidt stability", "rLakeAnalyzer::schmidt.stability; freshwater salinity = 0",
    "Metalimnion", "rLakeAnalyzer::meta.depths; seasonal thermocline",
    "Wind source", "TB4 WindSpeed and WindDir; centered 6-hour forcing window",
    "Wind-speed aggregation", "Root-mean-square WindSpeed across the 6-hour window; WindDir is the speed-weighted circular vector direction",
    "Wind height", paste(wind_height_m, "m (10 m reference assumed)"),
    "Calm-wind cutoff", paste(minimum_wind_speed_m_s, "m s-1; Lake Number set NA below cutoff"),
    "Lake Number", "rLakeAnalyzer::lake.number; WindDir retained as a diagnostic but the standard formula uses wind-speed magnitude"
  ),
  methods_file
)

metric_daily <- cast_metrics %>%
  as_tibble() %>%
  group_by(date) %>%
  summarise(
    schmidt_stability_j_m2 = mean(schmidt_stability_j_m2, na.rm = TRUE),
    lake_number = mean(lake_number, na.rm = TRUE),
    n_casts = n(), .groups = "drop"
  ) %>%
  mutate(
    across(c(schmidt_stability_j_m2, lake_number), ~if_else(is.nan(.x), NA_real_, .x)),
    year = year(date), month = month(date)
  )

make_metric_data <- function(value_col, metric_name) {
  climatology_source <- metric_daily %>%
    filter(!year %in% c(2021L, 2022L), is.finite(.data[[value_col]]),
           if (value_col == "lake_number") .data[[value_col]] > 0 else TRUE) %>%
    group_by(year, month) %>%
    summarise(year_month_value = mean(.data[[value_col]]), .groups = "drop")
  if (value_col == "lake_number") {
    climatology <- climatology_source %>%
      group_by(month) %>%
      summarise(
        n = n(), log_mean = mean(log10(year_month_value)),
        log_sd = sd(log10(year_month_value)), .groups = "drop"
      ) %>%
      filter(n >= 3L) %>%
      mutate(
        clim_mean = 10^log_mean,
        clim_lo = 10^(log_mean - replace_na(log_sd, 0)),
        clim_hi = 10^(log_mean + replace_na(log_sd, 0))
      )
  } else {
    climatology <- climatology_source %>%
      group_by(month) %>%
      summarise(
        n = n(), clim_mean = mean(year_month_value),
        clim_sd = sd(year_month_value), .groups = "drop"
      ) %>%
      filter(n >= 3L) %>%
      mutate(
        clim_lo = pmax(clim_mean - replace_na(clim_sd, 0), 0),
        clim_hi = clim_mean + replace_na(clim_sd, 0)
      )
  }

  months <- tibble(
    month_start = seq(floor_date(stability_start, "month"),
                      floor_date(stability_end, "month"), by = "month")
  ) %>%
    mutate(date = month_start + days(14), month = month(month_start)) %>%
    left_join(climatology, by = "month")
  focal <- metric_daily %>%
    filter(date >= stability_start, date <= stability_end,
           is.finite(.data[[value_col]])) %>%
    transmute(metric = metric_name, series = "2019-2023 observations", date,
              value = .data[[value_col]], lower = NA_real_, upper = NA_real_,
              n = n_casts) %>%
    arrange(date) %>%
    mutate(
      segment = if (value_col == "lake_number") {
        cumsum(replace_na(as.numeric(date - lag(date)) > 90, TRUE))
      } else 1L
    )
  clim_export <- months %>%
    transmute(metric = metric_name, series = "Other-year monthly climatology",
              date, value = clim_mean, lower = clim_lo, upper = clim_hi, n)
  list(focal = focal, climatology = months,
       export = bind_rows(focal, clim_export))
}

schmidt_data <- make_metric_data(
  "schmidt_stability_j_m2", "Schmidt stability"
)
lake_data <- make_metric_data("lake_number", "Lake Number")
if (!s7_figure_only) {
  write_csv(bind_rows(schmidt_data$export, lake_data$export), metric_plot_file, na = "")
}

make_metric_plot <- function(metric_data, y_label, title,
                             obs_colour, clim_colour, log_y = FALSE) {
  p <- ggplot() +
    annotate("rect", xmin = fire_start, xmax = fire_end, ymin = -Inf, ymax = Inf,
             fill = "firebrick", alpha = 0.07) +
    geom_vline(xintercept = fire_start, colour = "firebrick", linetype = "dashed",
               linewidth = LO_LW_MID) +
    annotate("text", x = fire_start + (fire_end - fire_start) / 2,
             y = Inf, label = "Caldor Fire", colour = "firebrick",
             family = LO_FONT, size = 3, vjust = 1.25) +
    geom_ribbon(
      data = metric_data$climatology,
      aes(x = date, ymin = clim_lo, ymax = clim_hi),
      fill = clim_colour, alpha = 0.28, colour = NA, na.rm = TRUE
    ) +
    geom_line(
      data = metric_data$climatology,
      aes(x = date, y = clim_mean, linetype = "Climatology"),
      colour = clim_colour, linewidth = LO_LW_THICK, na.rm = TRUE
    ) +
    geom_line(
      data = metric_data$focal,
      aes(x = date, y = value, linetype = "2019-2023", group = segment),
      colour = obs_colour, linewidth = LO_LW_THICK, na.rm = TRUE
    ) +
    geom_point(
      data = metric_data$focal, aes(x = date, y = value),
      colour = obs_colour, fill = "white", shape = 21,
      size = 1.8, stroke = 0.5, na.rm = TRUE
    ) +
    scale_linetype_manual(
      values = c("2019-2023" = "solid", "Climatology" = "dashed"),
      breaks = c("2019-2023", "Climatology"), name = NULL
    ) +
    scale_x_date(
      limits = c(stability_start, stability_end), date_breaks = "1 year",
      date_labels = "%Y", expand = expansion(mult = c(0, 0))
    ) +
    labs(title = title, subtitle = NULL, x = NULL, y = y_label) +
    lo_theme(base_size = 9, family = LO_FONT) +
    theme(
      panel.grid.minor = element_blank(), legend.position = "bottom",
      legend.title = element_blank(), legend.key.width = unit(8, "mm"),
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_blank()
    )
  if (log_y) {
    p <- p + scale_y_log10(labels = label_scientific(digits = 2))
  } else {
    p <- p + scale_y_continuous(
      labels = label_number(accuracy = 1, big.mark = ","),
      expand = expansion(mult = c(0.02, 0.08))
    )
  }
  p
}

p_schmidt <- make_metric_plot(
  schmidt_data,
  expression("Schmidt stability (J m"^{-2}*")"),
  "Mid-lake (MLTP) thermal stratification",
  "#6A3D9A", "#CAB2D6"
)
p_lake_number <- make_metric_plot(
  lake_data, "Lake Number (dimensionless)",
  "Mid-lake (MLTP) resistance to wind mixing",
  "#0072B2", "#9ECAE1", log_y = TRUE
)

ggsave(
  schmidt_file, p_schmidt, width = LO_WIDTH_MAX, height = 8.0, units = "cm",
  dpi = LO_DPI_LINE, device = ragg::agg_png
)
ggsave(
  schmidt_pdf, p_schmidt, width = LO_WIDTH_MAX, height = 8.0, units = "cm",
  device = cairo_pdf
)
if (!s7_figure_only) {
  ggsave(
    lake_number_file, p_lake_number, width = LO_WIDTH_MAX, height = 8.0, units = "cm",
    dpi = LO_DPI_LINE, device = ragg::agg_png
  )
  ggsave(
    combined_stability_file, p_schmidt / p_lake_number,
    width = LO_WIDTH_MAX, height = LO_HEIGHT_MAX, units = "cm", dpi = LO_DPI_LINE,
    device = ragg::agg_png
  )
}

cat("\n=== Pickle CTD profiles and stability ===\n")
cat("Profile window:", as.character(profile_start), "to", as.character(profile_end), "\n")
cat("Focal profile casts by station:\n")
print(profile_availability %>%
        group_by(Station_ID) %>%
        summarise(focal_casts = sum(focal_casts), .groups = "drop"), n = Inf)
cat("Schmidt observations in 2019-2023:", nrow(schmidt_data$focal), "\n")
cat("Lake Number observations in 2019-2023:", nrow(lake_data$focal), "\n")
cat("Outputs:\n", paste(c(
  profile_mltp_file, profile_ltp_file, schmidt_file, lake_number_file,
  combined_stability_file, metrics_file, metric_plot_file, methods_file
), collapse = "\n"), "\n")
