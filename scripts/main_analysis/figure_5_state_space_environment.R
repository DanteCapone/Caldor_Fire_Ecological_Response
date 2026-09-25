# Stage 1 of the Figure 5C compound-disturbance state-space analysis.
# This script uses smoke, optical, deposition, radiation, and CTD observations.
# It freezes temporal windows and environmental PCA variable lists before any
# biological response data are read by stage 2.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(yaml)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")
source("scripts/main_analysis/figure_5_state_space_helpers.R")

PROC_DIR <- file.path("data", "processed", "figure_5_state_space")
FIG_DIR <- file.path("figures", "supplemental", "figure_5_state_space")
CONFIG_FILE <- file.path("project_state", "FIGURE_5_STATE_SPACE_CONFIG.yml")
ENVIRONMENT_FILE <- file.path(PROC_DIR, "environmental_state_definitions.rds")
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PRIMARY_SMOKE_DEFINITION <- "fixed_aug01_oct31"
PRIMARY_ANTECEDENT_RULE <- "fixed_jul01_aug31"
PRIMARY_ANTECEDENT_START_DOY <- standard_doy(as.Date("2001-07-01"))
PRIMARY_ANTECEDENT_END_DOY <- standard_doy(as.Date("2001-08-31"))
SMOKE_CLIMATOLOGY_EXCLUDED_YEARS <- c(2020L, 2021L)
OPTICAL_R2_MINIMUM <- 0.90
SUSTAIN_DAYS <- 14L
MINIMUM_COVERAGE <- 0.80
REDUNDANCY_CUTOFF <- 0.90

smoke_windows_spec <- tribble(
  ~smoke_definition, ~lower_probability, ~upper_probability,
  "p10_p90", 0.10, 0.90,
  "p05_p95", 0.05, 0.95,
  "p025_p975", 0.025, 0.975
)
antecedent_diagnostic_metadata <- tribble(
  ~variable, ~label, ~direction, ~priority,
  "schmidt_stability_j_m2", "Schmidt stability (J m^-2)", 1, 1,
  "delta_t_mixed_layer_c", "Mixed-layer thermal contrast (deg C)", 1, 2,
  "mixed_layer_depth_m", "Mixed-layer depth (m)", -1, 3,
  "surface_temperature_0_10m_c", "Surface temperature, 0-10 m (deg C)", 1, 4
)
# Whole-lake Schmidt stability requires casts reaching at least 400 m; the
# July-August 2023 casts reach only 101-106 m. Keep Schmidt stability as a
# separate antecedent diagnostic and fit the no-imputation PC to the three
# profile metrics observed in every analysis year.
state_metadata <- tribble(
  ~variable, ~label, ~direction, ~priority,
  "delta_t_mixed_layer_c", "Mixed-layer thermal contrast (deg C)", 1, 1,
  "mixed_layer_depth_m", "Mixed-layer depth (m)", -1, 2,
  "surface_temperature_0_10m_c", "Surface temperature, 0-10 m (deg C)", 1, 3
)
wildfire_metadata <- tribble(
  ~variable, ~label, ~direction, ~priority,
  "integrated_pm25_excess_ug_d_m3", "Integrated PM2.5 excess (ug d m^-3)", 1, 1
)

pm25_priority <- c(
  "Tahoe_City", "ARB_Daily_Average", "South_Lake", "Bliss",
  "LT_College", "Truckee", "ARB_Daily_Max"
)
pm25_compiled <- read_csv(
  file.path("data", "processed", "tahoe_pm25_compiled.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date))
pm25_priority <- intersect(pm25_priority, names(pm25_compiled))
pm25_corrected_daily <- pm25_compiled %>%
  transmute(
    date,
    pm25_max = pmap_dbl(select(., all_of(pm25_priority)), function(...) {
      values <- suppressWarnings(as.numeric(c(...)))
      valid <- is.finite(values) & values >= 0 & values <= 500
      if (any(valid)) values[which(valid)[1]] else NA_real_
    }),
    pm25_source = pmap_chr(select(., all_of(pm25_priority)), function(...) {
      values <- suppressWarnings(as.numeric(c(...)))
      valid <- is.finite(values) & values >= 0 & values <= 500
      if (any(valid)) pm25_priority[which(valid)[1]] else NA_character_
    })
  )
smoke_daily <- read_csv(
  file.path("data", "processed", "tahoe_hms_pm25_daily.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    date = as.Date(date), smoke_day, max_density,
    pm25_existing = as.numeric(pm25_max)
  ) %>%
  left_join(
    pm25_corrected_daily %>%
      rename(pm25_priority_fallback = pm25_max,
             pm25_priority_fallback_source = pm25_source),
    by = "date"
  ) %>%
  mutate(
    pm25_max = coalesce(pm25_existing, pm25_priority_fallback),
    pm25_source = if_else(
      is.finite(pm25_existing), "tahoe_hms_pm25_daily",
      pm25_priority_fallback_source
    ),
    year = year(date), month = month(date), doy = standard_doy(date),
    smoke_confirmed = coalesce(smoke_day, FALSE) & is.finite(pm25_max) & pm25_max >= 9.1
  ) %>%
  select(-pm25_existing, -pm25_priority_fallback, -pm25_priority_fallback_source)
write_csv(
  smoke_daily %>% select(date, pm25_max, pm25_source),
  file.path(PROC_DIR, "pm25_resolved_daily.csv")
)
shortwave_daily <- read_csv(
  file.path("data", "processed", "tahoe_hms_swrad_daily.csv"),
  show_col_types = FALSE
) %>%
  transmute(date = as.Date(date), sw_diff)
smoke_radiation_daily <- smoke_daily %>%
  left_join(shortwave_daily, by = "date")

optical_depth <- read_csv(
  file.path("data", "lake_environmental_data", "uv", "LTP_Kd_Thermocline_depth_results.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    date = as.Date(Date),
    year = year(date),
    doy = standard_doy(date),
    par_1pct_depth_m = if_else(
      coalesce(R_PAR, 0) >= OPTICAL_R2_MINIMUM & coalesce(Kd_PAR, 0) > 0,
      4.605 / Kd_PAR, NA_real_
    ),
    uv320_1pct_depth_m = if_else(
      coalesce(R_320, 0) >= OPTICAL_R2_MINIMUM & coalesce(Kd_320, 0) > 0,
      4.605 / Kd_320, NA_real_
    )
  )

deposition <- read_csv(
  file.path("data", "lake_environmental_data", "deposition",
            "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    start_date = as.Date(Start_Datetime),
    end_date = as.Date(End_Datetime),
    tn_daily_load_mg_m2_d = as.numeric(TN_Daily_Load),
    tp_daily_load_mg_m2_d = as.numeric(TP_Daily_Load)
  ) %>%
  filter(!is.na(start_date), !is.na(end_date), end_date > start_date)

source_overlap_years <- max(min(optical_depth$year, na.rm = TRUE),
                            year(min(deposition$start_date, na.rm = TRUE))):
  min(max(optical_depth$year, na.rm = TRUE),
      year(max(deposition$end_date, na.rm = TRUE)))
optical_complete_years <- optical_depth %>%
  filter(doy >= standard_doy(as.Date("2001-08-01")),
         doy <= standard_doy(as.Date("2001-10-31"))) %>%
  group_by(year) %>%
  summarise(
    par_available = any(is.finite(par_1pct_depth_m)),
    uv320_available = any(is.finite(uv320_1pct_depth_m)),
    .groups = "drop"
  ) %>%
  filter(par_available, uv320_available) %>%
  pull(year)
deposition_complete_years <- summarize_deposition_window(
  deposition, source_overlap_years,
  standard_doy(as.Date("2001-08-01")),
  standard_doy(as.Date("2001-10-31")),
  minimum_interval_coverage = 0
) %>%
  filter(deposition_interval_coverage >= MINIMUM_COVERAGE) %>%
  pull(year)
pm25_complete_years <- smoke_radiation_daily %>%
  filter(doy >= standard_doy(as.Date("2001-08-01")),
         doy <= standard_doy(as.Date("2001-10-31"))) %>%
  group_by(year) %>%
  summarise(coverage = mean(is.finite(pm25_max)), .groups = "drop") %>%
  filter(coverage >= MINIMUM_COVERAGE) %>%
  pull(year)
COMPLETE_SOURCE_YEARS <- Reduce(
  intersect,
  list(optical_complete_years, deposition_complete_years, pm25_complete_years)
)
ANALYSIS_YEARS <- 2006:2025
PCA_FOCAL_YEARS <- c(2011L, 2020L, 2021L)

historical_smoke_days <- smoke_radiation_daily %>%
  filter(!(year %in% SMOKE_CLIMATOLOGY_EXCLUDED_YEARS),
         smoke_confirmed %in% TRUE, is.finite(pm25_max))
if (nrow(historical_smoke_days) == 0L) stop("No historical operational smoke days found.")

smoke_windows <- smoke_windows_spec %>%
  mutate(window = pmap(
    list(lower_probability, upper_probability),
    ~quantile_window(historical_smoke_days$doy, ..1, ..2)
  )) %>%
  select(smoke_definition, window) %>%
  unnest(window) %>%
  mutate(
    start_month_day = format(date_from_standard_doy(2001L, start_doy), "%m-%d"),
    end_month_day = format(date_from_standard_doy(2001L, end_doy), "%m-%d"),
    n_historical_smoke_days = nrow(historical_smoke_days),
    historical_year_start = min(historical_smoke_days$year),
    historical_year_end = max(historical_smoke_days$year),
    excluded_years = paste(SMOKE_CLIMATOLOGY_EXCLUDED_YEARS, collapse = ";")
  )
smoke_windows <- bind_rows(
  smoke_windows,
  tibble(
    smoke_definition = PRIMARY_SMOKE_DEFINITION,
    lower_probability = NA_real_, upper_probability = NA_real_,
    raw_start_doy = as.numeric(standard_doy(as.Date("2001-08-01"))),
    raw_end_doy = as.numeric(standard_doy(as.Date("2001-10-31"))),
    start_doy = standard_doy(as.Date("2001-08-01")),
    end_doy = standard_doy(as.Date("2001-10-31")),
    start_month_day = "08-01", end_month_day = "10-31",
    n_historical_smoke_days = nrow(historical_smoke_days),
    historical_year_start = min(historical_smoke_days$year),
    historical_year_end = max(historical_smoke_days$year),
    excluded_years = paste(SMOKE_CLIMATOLOGY_EXCLUDED_YEARS, collapse = ";")
  )
)
write_csv(smoke_windows, file.path(PROC_DIR, "smoke_season_definitions.csv"))

smoke_climatology <- smoke_radiation_daily %>%
  filter(!(year %in% SMOKE_CLIMATOLOGY_EXCLUDED_YEARS), is.finite(pm25_max)) %>%
  group_by(doy) %>%
  summarise(
    years_observed = n_distinct(year),
    qualifying_smoke_days = sum(smoke_confirmed %in% TRUE),
    smoke_occurrence_frequency = qualifying_smoke_days / years_observed,
    .groups = "drop"
  )
write_csv(smoke_climatology, file.path(PROC_DIR, "historical_smoke_day_climatology.csv"))

stability <- read_csv(
  file.path("data", "processed", "ctd", "mltp_schmidt_stability_2005_2025.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date), year = year(date), doy = standard_doy(date)) %>%
  filter(is.finite(schmidt_stability_j_m2))
stability_historical <- stability %>%
  filter(year != 2021L) %>%
  add_count(year, name = "casts_in_year") %>%
  mutate(equal_year_weight = 1 / casts_in_year)

stability_gam <- gam(
  schmidt_stability_j_m2 ~ s(doy, bs = "cc", k = 20),
  data = stability_historical,
  weights = equal_year_weight,
  method = "REML",
  knots = list(doy = c(0.5, 365.5))
)
stability_curve <- tibble(doy = 1:365) %>%
  mutate(
    fitted_stability = pmax(
      as.numeric(predict(stability_gam, newdata = pick(everything()), type = "response")),
      0
    )
  )
stability_onsets <- map_dfr(
  c(0.40, 0.50, 0.60),
  ~find_sustained_onset(stability_curve, .x, SUSTAIN_DAYS)
)
write_csv(stability_curve, file.path(PROC_DIR, "schmidt_stability_climatology_gam.csv"))
write_csv(stability_onsets, file.path(PROC_DIR, "stratification_onset_definitions.csv"))

# Reuse the project's quality-control and 0.1 kg m^-3 density-threshold MLD
# algorithm, applied here to every MLTP profile without year-specific choices.
source("R/ctd_mixing_helpers.R")
CTD_STATION <- "Mid-lake"
CTD_STATION_ALIAS <- "MLTP"
mltp_profiles <- load_ltp_ctd_profiles(".")
cast_metrics <- mltp_profiles %>%
  group_by(CTD_ID, Event_ID, date) %>%
  group_modify(~calculate_mltp_cast_metrics(.x, threshold_depth)) %>%
  ungroup() %>%
  mutate(date = as.Date(date), year = year(date), doy = standard_doy(date)) %>%
  filter(max_profile_depth_m >= MIN_PROFILE_DEPTH_M)
write_csv(cast_metrics, file.path(PROC_DIR, "mltp_cast_state_metrics.csv"))

# This diagnostic is descriptive and is not used to select variables, orient a
# PCA, or define the July-August window.
antecedent_all_years <- summarize_antecedent_state(
  stability, cast_metrics,
  PRIMARY_ANTECEDENT_START_DOY,
  standard_doy(as.Date("2001-08-01")),
  PRIMARY_ANTECEDENT_END_DOY
)
antecedent_2021_diagnostic <- antecedent_all_years %>%
  select(year, all_of(antecedent_diagnostic_metadata$variable)) %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  filter(is.finite(value)) %>%
  group_by(variable) %>%
  group_modify(~{
    focal_value <- .x$value[.x$year == 2021L][1]
    historical <- .x$value[.x$year != 2021L]
    historical_mad <- mad(historical, constant = 1.4826)
    historical_interval <- quantile(historical, c(0.025, 0.975), type = 8)
    tibble(
      year_2021_value = focal_value,
      n_historical_years = length(historical),
      historical_median = median(historical),
      historical_mean = mean(historical),
      historical_sd = sd(historical),
      historical_p025 = historical_interval[[1]],
      historical_p975 = historical_interval[[2]],
      z_score = (focal_value - mean(historical)) / sd(historical),
      robust_z_score = (focal_value - median(historical)) / historical_mad,
      historical_percentile = mean(historical <= focal_value),
      outside_historical_95_interval = focal_value < historical_interval[[1]] |
        focal_value > historical_interval[[2]]
    )
  }) %>%
  ungroup() %>%
  left_join(
    antecedent_diagnostic_metadata %>% select(variable, label, direction),
    by = "variable"
  ) %>%
  relocate(label, .after = variable)
write_csv(
  antecedent_2021_diagnostic,
  file.path(PROC_DIR, "antecedent_stratification_2021_diagnostic.csv")
)

primary_smoke <- smoke_windows %>% filter(smoke_definition == PRIMARY_SMOKE_DEFINITION)
primary_onset <- tibble(
  stability_fraction = NA_real_,
  onset_doy = PRIMARY_ANTECEDENT_START_DOY,
  onset_month_day = "07-01"
)
stopifnot(nrow(primary_smoke) == 1L, nrow(primary_onset) == 1L)

# The requested definitions can logically yield no antecedent interval. Freeze
# and report that outcome before attempting annual summaries or reading biology.
if (PRIMARY_ANTECEDENT_RULE == "stability_onset" &&
    primary_onset$onset_doy > primary_smoke$start_doy - 1L) {
  blocked_reason <- paste0(
    "The objective 50% stratification onset (DOY ", primary_onset$onset_doy,
    ") occurs after the primary 5th-percentile smoke-season start (DOY ",
    primary_smoke$start_doy, "), so the prescribed antecedent interval is empty."
  )
  blocked_config <- list(
    analysis = "Figure 5C multivariate compound-disturbance state space",
    analysis_status = "blocked_empty_primary_antecedent_window",
    blocking_reason = blocked_reason,
    biological_data_used_to_define_windows_or_pcas = FALSE,
    standard_day_of_year = "365-day calendar; leap-day offset removed after February",
    smoke_day_definition = "NOAA HMS smoke over Lake Tahoe and composite daily PM2.5 >= 9.1 ug m^-3",
    smoke_climatology_years = sort(unique(historical_smoke_days$year)),
    excluded_window_definition_year = 2021L,
    quantile_type = 8L,
    integer_window_rule = "floor lower percentile and ceiling upper percentile",
    primary_smoke_definition = PRIMARY_SMOKE_DEFINITION,
    primary_smoke_start_doy = primary_smoke$start_doy,
    primary_smoke_end_doy = primary_smoke$end_doy,
    primary_smoke_start_month_day = primary_smoke$start_month_day,
    primary_smoke_end_month_day = primary_smoke$end_month_day,
    primary_antecedent_rule = PRIMARY_ANTECEDENT_RULE,
    primary_stratification_onset_doy = primary_onset$onset_doy,
    primary_stratification_onset_month_day = primary_onset$onset_month_day,
    sustained_days_required = SUSTAIN_DAYS,
    requested_primary_antecedent_end_doy = primary_smoke$start_doy - 1L,
    candidate_antecedent_variables = state_metadata$variable,
    candidate_wildfire_variables = wildfire_metadata$variable,
    antecedent_variables_retained = character(),
    wildfire_variables_retained = character(),
    sensitivity_smoke_definitions = smoke_windows$smoke_definition,
    sensitivity_stability_fractions = stability_onsets$stability_fraction
  )
  write_yaml(blocked_config, CONFIG_FILE)

  p_smoke_blocked <- ggplot(smoke_climatology, aes(doy, smoke_occurrence_frequency)) +
    annotate("rect", xmin = primary_smoke$start_doy, xmax = primary_smoke$end_doy,
             ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.12) +
    geom_col(fill = "grey35", width = 1) +
    geom_vline(xintercept = c(primary_smoke$start_doy, primary_smoke$end_doy),
               linetype = 2, colour = "#D55E00", linewidth = 0.45) +
    scale_x_continuous(
      breaks = standard_doy(as.Date(paste0("2001-", sprintf("%02d", 1:12), "-01"))),
      labels = month.abb
    ) +
    labs(x = NULL, y = "Historical smoke-day frequency") +
    theme_classic(base_size = 9, base_family = "Times New Roman")
  onset_lines_blocked <- stability_onsets %>%
    mutate(label = paste0(round(100 * stability_fraction), "%"))
  p_stability_blocked <- ggplot(stability_curve, aes(doy, fitted_stability)) +
    geom_line(colour = "#0072B2", linewidth = 0.7) +
    geom_hline(data = onset_lines_blocked,
               aes(yintercept = threshold_j_m2, colour = label),
               linetype = 3, linewidth = 0.35) +
    geom_vline(data = onset_lines_blocked,
               aes(xintercept = onset_doy, colour = label),
               linetype = 2, linewidth = 0.45) +
    scale_colour_manual(
      values = c("40%" = "#009E73", "50%" = "#D55E00", "60%" = "#CC79A7"),
      name = "Seasonal amplitude"
    ) +
    scale_x_continuous(
      breaks = standard_doy(as.Date(paste0("2001-", sprintf("%02d", 1:12), "-01"))),
      labels = month.abb
    ) +
    labs(x = NULL, y = expression("Fitted Schmidt stability (J m"^{-2}*")")) +
    theme_classic(base_size = 9, base_family = "Times New Roman") +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "smoke_day_climatology_and_window.png"), p_smoke_blocked,
         width = 12.7, height = 7.5, units = "cm", dpi = 600, device = ragg::agg_png)
  ggsave(file.path(FIG_DIR, "smoke_day_climatology_and_window.pdf"), p_smoke_blocked,
         width = 12.7, height = 7.5, units = "cm", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, "stability_climatology_and_onsets.png"), p_stability_blocked,
         width = 12.7, height = 9, units = "cm", dpi = 600, device = ragg::agg_png)
  ggsave(file.path(FIG_DIR, "stability_climatology_and_onsets.pdf"), p_stability_blocked,
         width = 12.7, height = 9, units = "cm", device = cairo_pdf)

  sensitivity_intervals <- crossing(
    smoke_windows %>% select(smoke_definition, smoke_start_doy = start_doy),
    stability_onsets %>% select(stability_fraction, onset_doy)
  ) %>%
    mutate(
      antecedent_end_doy = smoke_start_doy - 1L,
      antecedent_window_days = antecedent_end_doy - onset_doy + 1L,
      interval_exists = antecedent_window_days > 0,
      n_years_with_mltp_profiles = map2_int(onset_doy, smoke_start_doy, ~{
        cast_metrics %>% filter(doy >= .x, doy < .y) %>% summarise(n = n_distinct(year)) %>% pull(n)
      }),
      focal_years_with_mltp_profiles = map2_chr(onset_doy, smoke_start_doy, ~{
        years <- cast_metrics %>% filter(doy >= .x, doy < .y, year %in% c(2011L, 2020L, 2021L)) %>%
          distinct(year) %>% arrange(year) %>% pull(year)
        if (length(years) == 0L) "none" else paste(years, collapse = ";")
      }),
      n_years_with_stability_casts = map2_int(onset_doy, smoke_start_doy, ~{
        stability %>% filter(doy >= .x, doy < .y) %>% summarise(n = n_distinct(year)) %>% pull(n)
      }),
      focal_years_with_stability_casts = map2_chr(onset_doy, smoke_start_doy, ~{
        years <- stability %>% filter(doy >= .x, doy < .y, year %in% c(2011L, 2020L, 2021L)) %>%
          distinct(year) %>% arrange(year) %>% pull(year)
        if (length(years) == 0L) "none" else paste(years, collapse = ";")
      })
    )
  write_csv(sensitivity_intervals, file.path(PROC_DIR, "window_compatibility_diagnostic.csv"))
  analysis_record <- c(
    "# Figure 5C state-space analysis record",
    "",
    "This record is written in the required chronological order. No biological response data were read.",
    "",
    "## A. Historical smoke climatology",
    paste0("Operational smoke days pooled from ", min(historical_smoke_days$year), "-",
           max(historical_smoke_days$year), ", excluding 2021 (n = ", nrow(historical_smoke_days), ")."),
    "",
    "## B. Derived wildfire-season dates",
    paste0("Primary 5th-95th percentile window: standard DOY ", primary_smoke$start_doy,
           "-", primary_smoke$end_doy, " (", primary_smoke$start_month_day, " to ",
           primary_smoke$end_month_day, ")."),
    "",
    "## C. Historical Schmidt-stability climatology",
    paste0("Equal-year-weighted cyclic GAM fitted to MLTP casts excluding 2021; sustained threshold requirement = ",
           SUSTAIN_DAYS, " days."),
    "",
    "## D. Derived antecedent-state dates",
    paste0("Primary 50% onset: standard DOY ", primary_onset$onset_doy, " (",
           primary_onset$onset_month_day, "); requested antecedent end: DOY ",
           primary_smoke$start_doy - 1L, "."),
    paste0("BLOCKED: ", blocked_reason),
    "",
    "## E. Variables retained from coverage and redundancy screens",
    "Not estimable because the primary antecedent interval is empty.",
    "",
    "## F. PCA loadings",
    "Not fit because the primary antecedent interval is empty.",
    "",
    "## G. Leptolyngbya response and Figure 5C",
    "Not calculated; the analysis stopped before biological data were read."
  )
  writeLines(analysis_record, file.path(PROC_DIR, "analysis_record.md"), useBytes = TRUE)
  stop(blocked_reason, " Resolve the window definition before fitting environmental PCAs.")
}

primary <- build_environment_scenario(
  stability, cast_metrics, smoke_radiation_daily,
  primary_onset$onset_doy, primary_smoke$start_doy, primary_smoke$end_doy,
  state_metadata, wildfire_metadata,
  antecedent_end_doy = PRIMARY_ANTECEDENT_END_DOY,
  optical_depth = optical_depth,
  deposition = deposition,
  analysis_years = ANALYSIS_YEARS,
  focal_years = PCA_FOCAL_YEARS,
  required_wildfire_variables = wildfire_metadata$variable,
  wildfire_index_variable = "integrated_pm25_excess_ug_d_m3",
  minimum_coverage = MINIMUM_COVERAGE,
  redundancy_cutoff = REDUNDANCY_CUTOFF
)
if (!all(PCA_FOCAL_YEARS %in% primary$scores$year)) {
  stop("The primary environmental PCA cohort does not include 2011, 2020, and 2021.")
}

# Leave-one-year-out refits quantify whether a single annual observation drives
# the displayed PC1. Align each refit to the full-data loading before comparing
# loading direction and overlapping-year scores.
full_state_fit <- primary$state_pca$fit
full_state_values <- primary$state_pca$completed %>%
  select(all_of(primary$state_variables))
full_pc1_loading <- full_state_fit$rotation[, "PC1"]
full_pc1_scores <- full_state_fit$x[, "PC1"]
primary_pca_leave_one_year_out <- map_dfr(
  seq_len(nrow(full_state_values)),
  function(omitted_i) {
    refit <- prcomp(
      full_state_values[-omitted_i, , drop = FALSE],
      center = TRUE, scale. = TRUE
    )
    refit_loading <- refit$rotation[, "PC1"]
    refit_scores <- refit$x[, "PC1"]
    if (sum(refit_loading * full_pc1_loading) < 0) {
      refit_loading <- -refit_loading
      refit_scores <- -refit_scores
    }
    tibble(
      omitted_year = primary$state_pca$completed$year[omitted_i],
      pc1_loading_cosine_similarity = sum(refit_loading * full_pc1_loading) /
        sqrt(sum(refit_loading^2) * sum(full_pc1_loading^2)),
      overlapping_score_correlation = cor(
        refit_scores, full_pc1_scores[-omitted_i]
      ),
      pc1_variance_explained = refit$sdev[1]^2 / sum(refit$sdev^2)
    )
  }
)
write_csv(
  primary_pca_leave_one_year_out,
  file.path(PROC_DIR, "primary_antecedent_pca_leave_one_year_out.csv")
)

scenario_grid <- tibble(
  smoke_definition = PRIMARY_SMOKE_DEFINITION,
  smoke_start_doy = primary_smoke$start_doy,
  smoke_end_doy = primary_smoke$end_doy,
  stability_fraction = NA_real_,
  onset_doy = PRIMARY_ANTECEDENT_START_DOY,
  antecedent_end_doy = PRIMARY_ANTECEDENT_END_DOY,
  scenario_id = "fixed_aug01_oct31_with_jul01_aug31_antecedent"
)
sensitivity_scenarios <- pmap(
  scenario_grid,
  function(smoke_definition, smoke_start_doy, smoke_end_doy,
           stability_fraction, onset_doy, antecedent_end_doy, scenario_id) {
    scenario <- build_environment_scenario(
      stability, cast_metrics, smoke_radiation_daily,
      onset_doy, smoke_start_doy, smoke_end_doy,
      state_metadata, wildfire_metadata,
      antecedent_end_doy = antecedent_end_doy,
      optical_depth = optical_depth,
      deposition = deposition,
      analysis_years = ANALYSIS_YEARS,
      focal_years = PCA_FOCAL_YEARS,
      fixed_state_variables = primary$state_variables,
      fixed_wildfire_variables = primary$wildfire_variables,
      required_wildfire_variables = wildfire_metadata$variable,
      wildfire_index_variable = "integrated_pm25_excess_ug_d_m3",
      minimum_coverage = MINIMUM_COVERAGE,
      redundancy_cutoff = REDUNDANCY_CUTOFF
    )
    scenario$metadata <- tibble(
      smoke_definition = smoke_definition,
      smoke_start_doy = smoke_start_doy,
      smoke_end_doy = smoke_end_doy,
      stability_fraction = stability_fraction,
      onset_doy = onset_doy,
      antecedent_end_doy = antecedent_end_doy,
      scenario_id = scenario_id
    )
    scenario
  }
)
names(sensitivity_scenarios) <- scenario_grid$scenario_id

availability_long <- bind_rows(
  primary$state_screen$availability %>%
    pivot_longer(-year, names_to = "variable", values_to = "available") %>%
    mutate(block = "Antecedent lake state"),
  primary$wildfire_screen$availability %>%
    pivot_longer(-year, names_to = "variable", values_to = "available") %>%
    mutate(block = "Wildfire exposure")
)
write_csv(availability_long, file.path(PROC_DIR, "year_variable_availability.csv"))
write_csv(
  bind_rows(
    primary$state_screen$coverage %>% mutate(block = "Antecedent lake state"),
    primary$wildfire_screen$coverage %>% mutate(block = "Wildfire exposure")
  ),
  file.path(PROC_DIR, "environmental_variable_screen.csv")
)
write_csv(
  bind_rows(
    matrix_to_long(primary$state_screen$correlation, "Antecedent lake state"),
    matrix_to_long(primary$wildfire_screen$correlation, "Wildfire exposure")
  ),
  file.path(PROC_DIR, "environmental_variable_correlations.csv")
)
write_csv(primary$annual, file.path(PROC_DIR, "primary_environmental_annual_metrics.csv"))
write_csv(primary$scores, file.path(PROC_DIR, "primary_environmental_pca_scores.csv"))
write_csv(
  primary$state_pca$completed %>%
    left_join(
      primary$annual %>%
        select(year, integrated_pm25_excess_ug_d_m3),
      by = "year"
    ),
  file.path(PROC_DIR, "primary_environmental_pca_completed_inputs.csv")
)
write_csv(
  primary$state_pca$loadings,
  file.path(PROC_DIR, "primary_environmental_pca_loadings.csv")
)
write_csv(
  tibble(
    variable = names(primary$state_pca$fit$center),
    center_2006_2025 = unname(primary$state_pca$fit$center),
    sample_sd_2006_2025 = unname(primary$state_pca$fit$scale)
  ) %>%
    left_join(
      as_tibble(primary$state_pca$fit$rotation, rownames = "variable") %>%
        rename_with(~paste0(tolower(.x), "_loading"), starts_with("PC")),
      by = "variable"
    ) %>%
    left_join(state_metadata, by = "variable"),
  file.path(PROC_DIR, "primary_antecedent_pca_parameters.csv")
)
write_csv(
  primary$state_pca$variance,
  file.path(PROC_DIR, "primary_environmental_pca_variance.csv")
)
write_csv(
  map_dfr(sensitivity_scenarios, ~cross_join(.x$metadata, .x$scores)),
  file.path(PROC_DIR, "sensitivity_environmental_pca_scores.csv")
)

config <- list(
  analysis = "Figure 5C multivariate compound-disturbance state space",
  biological_data_used_to_define_windows_or_pcas = FALSE,
  standard_day_of_year = "365-day calendar; leap-day offset removed after February",
  smoke_day_definition = "NOAA HMS smoke over Lake Tahoe and composite daily PM2.5 >= 9.1 ug m^-3",
  smoke_climatology_years = sort(unique(historical_smoke_days$year)),
  excluded_smoke_climatology_years = SMOKE_CLIMATOLOGY_EXCLUDED_YEARS,
  quantile_type = 8L,
  integer_window_rule = "floor lower percentile and ceiling upper percentile",
  primary_smoke_definition = PRIMARY_SMOKE_DEFINITION,
  primary_smoke_start_doy = primary_smoke$start_doy,
  primary_smoke_end_doy = primary_smoke$end_doy,
  primary_smoke_start_month_day = primary_smoke$start_month_day,
  primary_smoke_end_month_day = primary_smoke$end_month_day,
  primary_antecedent_rule = PRIMARY_ANTECEDENT_RULE,
  primary_antecedent_start_doy = PRIMARY_ANTECEDENT_START_DOY,
  primary_antecedent_start_month_day = "07-01",
  primary_antecedent_end_doy = PRIMARY_ANTECEDENT_END_DOY,
  primary_antecedent_end_month_day = "08-31",
  antecedent_and_wildfire_windows_overlap = TRUE,
  sustained_days_required = SUSTAIN_DAYS,
  minimum_annual_daily_coverage = MINIMUM_COVERAGE,
  redundancy_correlation_cutoff = REDUNDANCY_CUTOFF,
  antecedent_variables_retained = primary$state_variables,
  wildfire_variables_retained = primary$wildfire_variables,
  pca_years = primary$scores$year,
  requested_analysis_years = ANALYSIS_YEARS,
  complete_case_all_auxiliary_forcing_years = COMPLETE_SOURCE_YEARS,
  focal_years_required = PCA_FOCAL_YEARS,
  wildfire_axis = "integrated_pm25_excess_ug_d_m3; native units, not a principal component",
  pm25_daily_source_rule = paste(
    "retain the operational tahoe_hms_pm25_daily value when finite; otherwise use the first valid",
    "0-500 ug m^-3 value in priority order Tahoe_City, ARB_Daily_Average, South_Lake,",
    "Bliss, LT_College, Truckee, ARB_Daily_Max"
  ),
  environmental_imputation_used = FALSE,
  antecedent_pca_calculation = paste(
    "For each variable x, z=(x-mean_2006_2025)/sample_sd_2006_2025;",
    "PC1=sum(z_j * loading_j). PCA is stats::prcomp with center=TRUE and scale.=TRUE,",
    "fit by singular-value decomposition. Sign is selected by positive alignment with",
    "the predeclared direction vector (1,-1,1); individual loadings are not constrained."
  ),
  optical_depth_definition = "1% depth = 4.605 / Kd; fitted profile R2 >= 0.90 and Kd > 0",
  nutrient_deposition_definition = paste(
    "sum over exact overlapping collector days of TN_Daily_Load + TP_Daily_Load",
    "in mg m^-2 d^-1; annual integrated unit mg m^-2"
  ),
  biological_response_window = "August 1-December 31 of the same calendar year",
  plotted_biological_response_metric = paste(
    "maximum monthly August-December depth-integrated Leptolyngbya abundance",
    "in cells m^-2; bubble area proportional to peak abundance"
  ),
  zero_and_missing_response_symbols = paste(
    "x = observed zero abundance; open circle = biological observations absent",
    "while environmental coordinates remain available"
  ),
  thermal_contrast_definition = paste(
    "median across casts of mean temperature at depths <= density-threshold MLD",
    "minus mean temperature >MLD through MLD+10 m"
  ),
  annual_summary = "median of all qualifying casts in each fixed antecedent window",
  climatological_smoke_definitions_retained_as_diagnostics = smoke_windows_spec$smoke_definition,
  climatological_stability_fractions_retained_as_diagnostics = stability_onsets$stability_fraction
)
write_yaml(config, CONFIG_FILE)

environmental_record <- list(
  config = config,
  smoke_windows = smoke_windows,
  smoke_climatology = smoke_climatology,
  stability_curve = stability_curve,
  stability_onsets = stability_onsets,
  state_metadata = state_metadata,
  antecedent_diagnostic_metadata = antecedent_diagnostic_metadata,
  wildfire_metadata = wildfire_metadata,
  primary = primary,
  scenario_grid = scenario_grid,
  sensitivity_scenarios = sensitivity_scenarios
)
saveRDS(environmental_record, ENVIRONMENT_FILE)

primary_rect <- primary_smoke
p_smoke <- ggplot(smoke_climatology, aes(doy, smoke_occurrence_frequency)) +
  annotate("rect", xmin = primary_rect$start_doy, xmax = primary_rect$end_doy,
           ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.12) +
  geom_col(fill = "grey35", width = 1) +
  geom_vline(xintercept = c(primary_rect$start_doy, primary_rect$end_doy),
             linetype = 2, colour = "#D55E00", linewidth = 0.45) +
  scale_x_continuous(
    breaks = standard_doy(as.Date(paste0("2001-", sprintf("%02d", 1:12), "-01"))),
    labels = month.abb
  ) +
  labs(x = NULL, y = "Historical smoke-day frequency") +
  theme_classic(base_size = 9, base_family = "Times New Roman")

onset_lines <- stability_onsets %>%
  mutate(label = paste0(round(100 * stability_fraction), "%"))
p_stability <- ggplot(stability_curve, aes(doy, fitted_stability)) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_hline(data = onset_lines, aes(yintercept = threshold_j_m2, colour = label),
             linetype = 3, linewidth = 0.35) +
  geom_vline(data = onset_lines, aes(xintercept = onset_doy, colour = label),
             linetype = 2, linewidth = 0.45) +
  scale_colour_manual(values = c("40%" = "#009E73", "50%" = "#D55E00", "60%" = "#CC79A7"),
                      name = "Seasonal amplitude") +
  scale_x_continuous(
    breaks = standard_doy(as.Date(paste0("2001-", sprintf("%02d", 1:12), "-01"))),
    labels = month.abb
  ) +
  labs(x = NULL, y = expression("Fitted Schmidt stability (J m"^{-2}*")")) +
  theme_classic(base_size = 9, base_family = "Times New Roman") +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "smoke_day_climatology_and_window.png"), p_smoke,
       width = 12.7, height = 7.5, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "smoke_day_climatology_and_window.pdf"), p_smoke,
       width = 12.7, height = 7.5, units = "cm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "stability_climatology_and_onsets.png"), p_stability,
       width = 12.7, height = 9, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "stability_climatology_and_onsets.pdf"), p_stability,
       width = 12.7, height = 9, units = "cm", device = cairo_pdf)

analysis_record <- c(
  "# Figure 5C state-space analysis record",
  "",
  "This record is written in the required chronological order. Sections A-F are generated without reading biological data.",
  "",
  "## A. Historical smoke climatology",
  paste0("Operational smoke days pooled from ", min(historical_smoke_days$year), "-",
         max(historical_smoke_days$year), ", excluding 2020 and 2021 (n = ",
         nrow(historical_smoke_days), ")."),
  "",
  "## B. Derived wildfire-season dates",
  paste0("User-specified primary window: standard DOY ", primary_smoke$start_doy,
         "-", primary_smoke$end_doy, " (", primary_smoke$start_month_day, " to ",
         primary_smoke$end_month_day, ")."),
  "",
  "## C. Historical Schmidt-stability climatology",
  paste0("Equal-year-weighted cyclic GAM fitted to MLTP casts excluding 2021; sustained threshold requirement = ",
         SUSTAIN_DAYS, " days."),
  "",
  "## D. Fixed antecedent-state dates",
  paste0("User-specified antecedent window: standard DOY ",
         PRIMARY_ANTECEDENT_START_DOY, "-", PRIMARY_ANTECEDENT_END_DOY,
         " (07-01 to 08-31). This intentionally overlaps the August portion of the wildfire window."),
  paste0(
    "The independent 2021 diagnostic classifies July-August Schmidt stability as elevated but not extreme: ",
    round(antecedent_2021_diagnostic$year_2021_value[
      antecedent_2021_diagnostic$variable == "schmidt_stability_j_m2"
    ]), " J m^-2; historical percentile ",
    round(100 * antecedent_2021_diagnostic$historical_percentile[
      antecedent_2021_diagnostic$variable == "schmidt_stability_j_m2"
    ]), "%; outside historical 95% interval = ",
    antecedent_2021_diagnostic$outside_historical_95_interval[
      antecedent_2021_diagnostic$variable == "schmidt_stability_j_m2"
    ], "."
  ),
  "",
  "## E. Variables retained from coverage and redundancy screens",
  paste0("Antecedent: ", paste(primary$state_variables, collapse = ", "), "."),
  "Wildfire exposure is not a PCA: it is the native-unit integrated smoke-day PM2.5 excess above 9.1 ug m^-3 during August-October.",
  "Invalid negative monitor values are removed before applying the fixed monitor-priority rule. This corrects the 2019 coverage artifact and retains directly observed PM2.5 exposure for 2006-2025 without environmental imputation.",
  "",
  "## F. Antecedent PCA calculation and loadings",
  paste0(
    "Each retained July-August annual variable is centered by its 2006-2025 mean and divided by its sample standard deviation. ",
    "PCA is fit by singular-value decomposition of that standardized year-by-variable matrix. The antecedent index is each year's PC1 score, ",
    "equal to the sum of its standardized variable values multiplied by the corresponding PC1 loadings. The PC1 sign is oriented toward ",
    "larger thermal contrast, shallower mixed layers, and warmer surface water. Whole-lake Schmidt stability remains a separate diagnostic because July-August 2023 lacks a >=400 m cast. See primary_environmental_pca_loadings.csv and primary_environmental_pca_variance.csv."
  ),
  paste0(
    "Leave-one-year-out robustness: loading cosine similarity range ",
    sprintf("%.4f-%.4f", min(primary_pca_leave_one_year_out$pc1_loading_cosine_similarity),
            max(primary_pca_leave_one_year_out$pc1_loading_cosine_similarity)),
    "; overlapping-score correlation range ",
    sprintf("%.4f-%.4f", min(primary_pca_leave_one_year_out$overlapping_score_correlation),
            max(primary_pca_leave_one_year_out$overlapping_score_correlation)),
    "; PC1 variance-explained range ",
    sprintf("%.2f-%.2f%%", 100 * min(primary_pca_leave_one_year_out$pc1_variance_explained),
            100 * max(primary_pca_leave_one_year_out$pc1_variance_explained)), "."
  ),
  "",
  "The frozen configuration is project_state/FIGURE_5_STATE_SPACE_CONFIG.yml. Biological response is calculated only by stage 2."
)
writeLines(analysis_record, file.path(PROC_DIR, "analysis_record.md"), useBytes = TRUE)

message(
  "Environmental definitions frozen: smoke DOY ", primary_smoke$start_doy, "-",
  primary_smoke$end_doy, "; antecedent onset DOY ", primary_onset$onset_doy,
  "; antecedent window DOY ", PRIMARY_ANTECEDENT_START_DOY, "-",
  PRIMARY_ANTECEDENT_END_DOY, "; smoke climatology excludes 2020 and 2021. ",
  "Configuration written before biological analysis."
)
