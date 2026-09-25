# Reproducible statistics for the manuscript's radiative-pathway Results text.
#
# Figure 1 optical and temperature values use the same inputs, quality filters,
# depth calculations, and monthly reference periods as figure_2.R, which is
# sourced by caldor_forcing_v2.R. Recovery dates are reported separately from
# the pathway synthesis because that analysis uses a stricter non-disturbance
# reference and requires two consecutive observations inside its envelope.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
R2_MIN <- 0.90

KD_FILE <- file.path(
  "data", "lake_environmental_data", "uv",
  "LTP_Kd_Thermocline_depth_results.csv"
)
SECCHI_FILE <- file.path(
  "data", "lake_environmental_data", "secchi", "Secchi_LTP.csv"
)
CTD_FILE <- file.path(
  "data", "lake_environmental_data", "ctd", "ctd_casts",
  "ctd_data_2005_2025.csv"
)
SMOKE_DAILY_FILE <- file.path(
  "data", "processed", "tahoe_hms_pm25_daily.csv"
)
PM25_FILE <- file.path("data", "processed", "tahoe_pm25_compiled.csv")
OPTICAL_TEST_FILE <- file.path(
  "figures", "figure_2_radiative", "figure_2_radiative_tests.csv"
)
RECOVERY_MARKER_FILE <- file.path(
  "figures", "supplemental", "resistance_and_resilience",
  "pathway_resistance_resilience_markers.csv"
)
OUT_DIR <- file.path("figures", "manuscript_statistics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

monthly_climatology <- function(data, value_col, start_year) {
  data %>%
    filter(
      year >= start_year, year <= 2025L, year != 2021L,
      is.finite(.data[[value_col]])
    ) %>%
    group_by(month) %>%
    summarise(
      n_hist = n(),
      expected = mean(.data[[value_col]], na.rm = TRUE),
      scale_sd = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ci_half_width = qt(0.975, df = pmax(n_hist - 1L, 1L)) *
        scale_sd / sqrt(pmax(n_hist, 1L)),
      reference_lo = pmax(expected - ci_half_width, 0),
      reference_hi = expected + ci_half_width
    )
}

join_reference <- function(observations, reference) {
  observations %>%
    left_join(reference, by = "month") %>%
    mutate(
      difference = value - expected,
      percent_difference = 100 * difference / expected,
      anomaly_sd = difference / scale_sd,
      within_reference_ci = between(value, reference_lo, reference_hi)
    )
}

kd_raw <- read_csv(KD_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date), year = year(date), month = month(date)) %>%
  filter(!is.na(date))

kd_qc <- kd_raw %>%
  mutate(
    par_1pct_m = if_else(
      coalesce(R_PAR, 0) >= R2_MIN & coalesce(Kd_PAR, 0) > 0,
      4.605 / Kd_PAR, NA_real_
    ),
    uv320_1pct_m = if_else(
      coalesce(R_320, 0) >= R2_MIN & coalesce(Kd_320, 0) > 0,
      4.605 / Kd_320, NA_real_
    )
  )

par_reference <- monthly_climatology(kd_qc, "par_1pct_m", 2015L)
uv_reference <- monthly_climatology(kd_qc, "uv320_1pct_m", 2015L)

par_2021 <- kd_qc %>%
  filter(year == 2021L, is.finite(par_1pct_m)) %>%
  transmute(date, month, value = par_1pct_m) %>%
  join_reference(par_reference)

uv_2021 <- kd_qc %>%
  filter(year == 2021L, is.finite(uv320_1pct_m)) %>%
  transmute(date, month, value = uv320_1pct_m) %>%
  join_reference(uv_reference)

secchi_raw <- read_csv(SECCHI_FILE, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(Date_Time_Local), year = year(date), month = month(date),
    Secchi = as.numeric(Secchi)
  ) %>%
  filter(!is.na(date), is.finite(Secchi))

secchi_reference <- monthly_climatology(secchi_raw, "Secchi", 2015L)
secchi_2021 <- secchi_raw %>%
  filter(year == 2021L) %>%
  transmute(date, month, value = Secchi) %>%
  join_reference(secchi_reference)

ctd_surface <- read_csv(
  CTD_FILE,
  col_select = c(
    CTD_ID, Station_ID, Cast_Date_Time_Local, Cast_Flag, Depth, Temperature
  ),
  na = c("", "NA", "NULL"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(date = as.Date(Cast_Date_Time_Local)) %>%
  filter(
    Station_ID == "Index", Cast_Flag == 1,
    is.finite(Depth), between(Depth, 0, 10), is.finite(Temperature),
    !is.na(date)
  ) %>%
  group_by(CTD_ID, date) %>%
  summarise(surface_temp_c = mean(Temperature, na.rm = TRUE), .groups = "drop") %>%
  mutate(year = year(date), month = month(date))

temp_reference <- monthly_climatology(ctd_surface, "surface_temp_c", 2005L)
temp_2021 <- ctd_surface %>%
  filter(year == 2021L) %>%
  transmute(date, month, value = surface_temp_c) %>%
  join_reference(temp_reference)

fire_minimum <- function(data) {
  data %>%
    filter(between(date, FIRE_START, FIRE_END)) %>%
    slice_min(value, n = 1L, with_ties = FALSE)
}

par_fire_min <- fire_minimum(par_2021)
par_2021_min <- par_2021 %>% slice_min(value, n = 1L, with_ties = FALSE)
uv_fire_min <- fire_minimum(uv_2021)
uv_2021_min <- uv_2021 %>% slice_min(value, n = 1L, with_ties = FALSE)
secchi_fire_min <- fire_minimum(secchi_2021)

# The August-December window is tested on month-standardized anomalies so that
# normal seasonal changes in light penetration do not drive the comparison.
par_aug_dec <- par_2021 %>%
  filter(between(date, as.Date("2021-08-01"), as.Date("2021-12-31")))
par_aug_dec_test_two_sided <- wilcox.test(
  par_aug_dec$anomaly_sd, mu = 0, alternative = "two.sided", exact = FALSE
)
par_aug_dec_test_less <- wilcox.test(
  par_aug_dec$anomaly_sd, mu = 0, alternative = "less", exact = FALSE
)
par_aug_dec_monthly <- par_aug_dec %>%
  group_by(month) %>%
  summarise(
    n_profiles = n(),
    mean_anomaly_sd = mean(anomaly_sd),
    median_anomaly_sd = median(anomaly_sd),
    .groups = "drop"
  )
par_aug_dec_monthly_test <- wilcox.test(
  par_aug_dec_monthly$mean_anomaly_sd,
  mu = 0, alternative = "two.sided", exact = TRUE
)

aug_dec_window <- function(data) {
  data %>%
    filter(between(date, as.Date("2021-08-01"), as.Date("2021-12-31")))
}

monthly_anomaly_summary <- function(data) {
  data %>%
    group_by(month) %>%
    summarise(
      n_profiles = n(),
      mean_anomaly_sd = mean(anomaly_sd),
      median_anomaly_sd = median(anomaly_sd),
      .groups = "drop"
    )
}

uv_aug_dec <- aug_dec_window(uv_2021)
secchi_aug_dec <- aug_dec_window(secchi_2021)
uv_aug_dec_test <- wilcox.test(
  uv_aug_dec$anomaly_sd, mu = 0, alternative = "two.sided", exact = FALSE
)
secchi_aug_dec_test <- wilcox.test(
  secchi_aug_dec$anomaly_sd, mu = 0, alternative = "two.sided", exact = FALSE
)
uv_aug_dec_monthly <- monthly_anomaly_summary(uv_aug_dec)
secchi_aug_dec_monthly <- monthly_anomaly_summary(secchi_aug_dec)
uv_aug_dec_monthly_test <- wilcox.test(
  uv_aug_dec_monthly$mean_anomaly_sd,
  mu = 0, alternative = "two.sided", exact = TRUE
)
secchi_aug_dec_monthly_test <- wilcox.test(
  secchi_aug_dec_monthly$mean_anomaly_sd,
  mu = 0, alternative = "two.sided", exact = TRUE
)

temp_fire <- temp_2021 %>% filter(between(date, FIRE_START, FIRE_END))
temp_fire_summary <- temp_fire %>%
  summarise(
    n_casts = n(),
    first_date = min(date),
    last_date = max(date),
    mean_difference_c = mean(difference),
    sd_difference_c = sd(difference),
    minimum_difference_c = min(difference),
    minimum_date = date[which.min(difference)],
    maximum_difference_c = max(difference),
    maximum_date = date[which.max(difference)]
  )

smoke_daily <- read_csv(SMOKE_DAILY_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))
pm25 <- read_csv(PM25_FILE, show_col_types = FALSE) %>%
  transmute(date = as.Date(date), pm25_ug_m3 = as.numeric(PM25))

fire_smoke <- smoke_daily %>%
  select(date, smoke_day) %>%
  left_join(pm25, by = "date") %>%
  filter(
    between(date, FIRE_START, FIRE_END), smoke_day,
    is.finite(pm25_ug_m3), pm25_ug_m3 >= 9.1
  )

smoke_summary <- fire_smoke %>%
  summarise(
    smoke_days = n(),
    smoke_days_ge_35_5 = sum(pm25_ug_m3 >= 35.5),
    maximum_pm25_ug_m3 = max(pm25_ug_m3),
    maximum_pm25_date = date[which.max(pm25_ug_m3)]
  )

optical_tests <- read_csv(OPTICAL_TEST_FILE, show_col_types = FALSE)
recovery_markers <- read_csv(RECOVERY_MARKER_FILE, show_col_types = FALSE) %>%
  mutate(
    peak_date = as.Date(peak_date),
    recovery_date = as.Date(recovery_date),
    recovery_days_after_ignition = as.integer(recovery_date - FIRE_START)
  )

par_test <- optical_tests %>% filter(metric == "PAR 1% depth")
uv_test <- optical_tests %>% filter(metric == "UV-320 1% depth")
secchi_test <- optical_tests %>% filter(metric == "Secchi depth")
par_recovery <- recovery_markers %>%
  filter(Variable == "PAR 1% depth (m)")
secchi_recovery <- recovery_markers %>%
  filter(Variable == "Secchi depth (m)")

statistic_rows <- function(metric, row, context) {
  tibble(
    metric,
    context,
    date = row$date,
    observed = row$value,
    expected = row$expected,
    difference = row$difference,
    percent_difference = row$percent_difference,
    anomaly_sd = row$anomaly_sd,
    historical_n = row$n_hist,
    reference_ci_lower = row$reference_lo,
    reference_ci_upper = row$reference_hi
  )
}

focal_observations <- bind_rows(
  statistic_rows("PAR 1% depth", par_fire_min,
                 "minimum during 2021-08-14 through 2021-10-21"),
  statistic_rows("PAR 1% depth", par_2021_min,
                 "minimum in calendar year 2021"),
  statistic_rows("UV-320 1% depth", uv_fire_min,
                 "minimum during 2021-08-14 through 2021-10-21"),
  statistic_rows("UV-320 1% depth", uv_2021_min,
                 "minimum in calendar year 2021"),
  statistic_rows("Secchi depth", secchi_fire_min,
                 "minimum during 2021-08-14 through 2021-10-21")
) %>%
  mutate(
    units = "m",
    reference_definition = paste(
      "monthly mean and 95% t confidence interval from individual observations;",
      "2015-2025 excluding 2021; optical profiles require R2 >= 0.90"
    )
  )

summary_statistics <- tribble(
  ~metric, ~statistic, ~estimate, ~units, ~method_detail,
  "PM2.5", "maximum during fire window", smoke_summary$maximum_pm25_ug_m3,
  "ug m^-3", format(smoke_summary$maximum_pm25_date, "%Y-%m-%d"),
  "PM2.5", "qualifying smoke days during fire window", smoke_summary$smoke_days,
  "days", "NOAA HMS smoke and composite PM2.5 >= 9.1 ug m^-3",
  "PM2.5", "qualifying smoke days >= 35.5 ug m^-3", smoke_summary$smoke_days_ge_35_5,
  "days", "subset of qualifying smoke days",
  "PAR 1% depth", "Wilcoxon V", par_test$statistic, "test statistic", par_test$test,
  "PAR 1% depth", "Wilcoxon p", par_test$p_value, "p value", par_test$test,
  "PAR 1% depth", "2021 median anomaly", par_test$median_anomaly_sd, "SD", par_test$test,
  "PAR 1% depth", "2021 observations in test", par_test$n_2021, "observations", par_test$test,
  "PAR 1% depth", "August-December Wilcoxon V",
  unname(par_aug_dec_test_two_sided$statistic), "test statistic",
  "one-sample Wilcoxon signed-rank test of month-standardized anomalies; 2021-08-01 through 2021-12-31",
  "PAR 1% depth", "August-December two-sided Wilcoxon p",
  par_aug_dec_test_two_sided$p.value, "p value",
  "primary window test; two-sided",
  "PAR 1% depth", "August-December one-sided Wilcoxon p",
  par_aug_dec_test_less$p.value, "p value",
  "directional alternative is shallower; use only if specified a priori",
  "PAR 1% depth", "August-December median anomaly",
  median(par_aug_dec$anomaly_sd), "SD", paste(nrow(par_aug_dec), "observations"),
  "PAR 1% depth", "August-December mean anomaly",
  mean(par_aug_dec$anomaly_sd), "SD", paste(nrow(par_aug_dec), "observations"),
  "PAR 1% depth", "August-December observations in test",
  nrow(par_aug_dec), "observations", "2021-08-01 through 2021-12-31",
  "PAR 1% depth", "August-December monthly-mean sensitivity Wilcoxon V",
  unname(par_aug_dec_monthly_test$statistic), "test statistic",
  "five calendar-month means treated as independent units",
  "PAR 1% depth", "August-December monthly-mean sensitivity two-sided p",
  par_aug_dec_monthly_test$p.value, "p value",
  "exact Wilcoxon test; five calendar-month means treated as independent units",
  "PAR 1% depth", "first of two consecutive observations inside recovery envelope",
  par_recovery$recovery_days_after_ignition, "days after ignition",
  paste("date", par_recovery$recovery_date,
        "; non-disturbance reference excludes 2011, 2012, 2020, and 2021"),
  "UV-320 1% depth", "Wilcoxon V", uv_test$statistic, "test statistic", uv_test$test,
  "UV-320 1% depth", "Wilcoxon p", uv_test$p_value, "p value", uv_test$test,
  "UV-320 1% depth", "2021 median anomaly", uv_test$median_anomaly_sd, "SD", uv_test$test,
  "UV-320 1% depth", "2021 observations in test", uv_test$n_2021, "observations", uv_test$test,
  "UV-320 1% depth", "August-December Wilcoxon V",
  unname(uv_aug_dec_test$statistic), "test statistic",
  "one-sample Wilcoxon signed-rank test of month-standardized anomalies; 2021-08-01 through 2021-12-31",
  "UV-320 1% depth", "August-December two-sided Wilcoxon p",
  uv_aug_dec_test$p.value, "p value", "primary window test; two-sided",
  "UV-320 1% depth", "August-December median anomaly",
  median(uv_aug_dec$anomaly_sd), "SD", paste(nrow(uv_aug_dec), "observations"),
  "UV-320 1% depth", "August-December observations in test",
  nrow(uv_aug_dec), "observations", "2021-08-01 through 2021-12-31",
  "UV-320 1% depth", "August-December monthly-mean sensitivity two-sided p",
  uv_aug_dec_monthly_test$p.value, "p value",
  paste("exact Wilcoxon test;", nrow(uv_aug_dec_monthly), "calendar-month means treated as independent units"),
  "UV-320 1% depth", "recovery status", NA_real_, "not applicable",
  "no significant effect in pathway recovery analysis",
  "Secchi depth", "Wilcoxon V", secchi_test$statistic, "test statistic", secchi_test$test,
  "Secchi depth", "Wilcoxon p", secchi_test$p_value, "p value", secchi_test$test,
  "Secchi depth", "2021 median anomaly", secchi_test$median_anomaly_sd, "SD", secchi_test$test,
  "Secchi depth", "2021 observations in test", secchi_test$n_2021, "observations", secchi_test$test,
  "Secchi depth", "August-December Wilcoxon V",
  unname(secchi_aug_dec_test$statistic), "test statistic",
  "one-sample Wilcoxon signed-rank test of month-standardized anomalies; 2021-08-01 through 2021-12-31",
  "Secchi depth", "August-December two-sided Wilcoxon p",
  secchi_aug_dec_test$p.value, "p value", "primary window test; two-sided",
  "Secchi depth", "August-December median anomaly",
  median(secchi_aug_dec$anomaly_sd), "SD", paste(nrow(secchi_aug_dec), "observations"),
  "Secchi depth", "August-December observations in test",
  nrow(secchi_aug_dec), "observations", "2021-08-01 through 2021-12-31",
  "Secchi depth", "August-December monthly-mean sensitivity two-sided p",
  secchi_aug_dec_monthly_test$p.value, "p value",
  paste("exact Wilcoxon test;", nrow(secchi_aug_dec_monthly), "calendar-month means treated as independent units"),
  "Secchi depth", "first of two consecutive observations inside recovery envelope",
  secchi_recovery$recovery_days_after_ignition, "days after ignition",
  paste("date", secchi_recovery$recovery_date,
        "; non-disturbance reference excludes 2011, 2012, 2020, and 2021"),
  "Surface temperature (0-10 m)", "mean fire-window departure from monthly climatology",
  temp_fire_summary$mean_difference_c, "deg C", paste(temp_fire_summary$n_casts, "casts"),
  "Surface temperature (0-10 m)", "SD of fire-window departures",
  temp_fire_summary$sd_difference_c, "deg C", paste(temp_fire_summary$n_casts, "casts"),
  "Surface temperature (0-10 m)", "minimum fire-window departure",
  temp_fire_summary$minimum_difference_c, "deg C", format(temp_fire_summary$minimum_date, "%Y-%m-%d"),
  "Surface temperature (0-10 m)", "maximum fire-window departure",
  temp_fire_summary$maximum_difference_c, "deg C", format(temp_fire_summary$maximum_date, "%Y-%m-%d")
)

write_csv(
  focal_observations,
  file.path(OUT_DIR, "radiative_pathway_focal_observations.csv")
)
write_csv(
  summary_statistics,
  file.path(OUT_DIR, "radiative_pathway_manuscript_statistics.csv")
)
write_csv(
  par_aug_dec,
  file.path(OUT_DIR, "radiative_pathway_par_aug_dec_observations.csv")
)
write_csv(
  par_aug_dec_monthly,
  file.path(OUT_DIR, "radiative_pathway_par_aug_dec_monthly_sensitivity.csv")
)
write_csv(
  bind_rows(
    par_aug_dec %>% mutate(metric = "PAR 1% depth", .before = 1),
    uv_aug_dec %>% mutate(metric = "UV-320 1% depth", .before = 1),
    secchi_aug_dec %>% mutate(metric = "Secchi depth", .before = 1)
  ),
  file.path(OUT_DIR, "radiative_pathway_aug_dec_optical_observations.csv")
)
write_csv(
  bind_rows(
    par_aug_dec_monthly %>% mutate(metric = "PAR 1% depth", .before = 1),
    uv_aug_dec_monthly %>% mutate(metric = "UV-320 1% depth", .before = 1),
    secchi_aug_dec_monthly %>% mutate(metric = "Secchi depth", .before = 1)
  ),
  file.path(OUT_DIR, "radiative_pathway_aug_dec_monthly_sensitivity.csv")
)

cat("Smoke summary\n")
print(smoke_summary)
cat("\nFocal optical observations\n")
print(focal_observations, n = Inf)
cat("\nFire-window temperature summary\n")
print(temp_fire_summary)
cat("\nOptical tests and recovery\n")
print(summary_statistics, n = Inf)
