# =============================================================================
# pathway_resistance_resilience.R
# Caldor Fire Ecosystem Response Project
#
# Supplemental synthesis of pathway-specific resistance and resilience.
#
# Each observation is compared with a seasonally explicit historical reference
# state estimated from non-disturbance, pre-fire years. Following Figure 2,
# monthly means are first calculated within each historical year and the
# reference envelope is the 95% t confidence interval around the across-year
# monthly climatological mean. Primary
# recovery is the first of two consecutive observations inside that envelope,
# measured from 14 August 2021. Sensitivity estimates use one and three required
# consecutive observations.
#
# Outputs:
#   data/processed/pathway_resistance_resilience_summary.csv
#   figures/supplemental/resistance_and_resilience/
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readxl)
  library(vegan)
  library(mgcv)
  library(patchwork)
  library(ragg)
})

# ---- Project paths and shared plotting settings -----------------------------
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
proj_root <- if (length(file_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), winslash = "/")
} else {
  normalizePath(".", winslash = "/")
}
setwd(proj_root)
source("scripts/figure_aesthetics.R")

DATA_DIR <- file.path(proj_root, "data", "lake_environmental_data")
PROC_DIR <- file.path(proj_root, "data", "processed")
OUT_DIR <- file.path(proj_root, "figures", "supplemental",
                     "resistance_and_resilience")
PLOT_DIR <- file.path(OUT_DIR, "individual_plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
INITIAL_RESPONSE_END <- as.Date("2022-12-31")
FOLLOW_UP_END <- as.Date("2025-12-31")
PEAK_WINDOW_END <- FIRE_END %m+% months(6)
FIRE_YEAR <- 2021L
MONTH_DAYS <- 30.44
NON_DISTURBANCE_END <- FIRE_START - 1L
# Known non-reference periods are excluded: the 2011 Leptolyngbya event and
# its 2012 recovery window, the 2020 smoke event, and the 2021 fire year.
# All post-fire observations are additionally excluded by date.
EXCLUDED_BASELINE_YEARS <- c(2011L, 2012L, 2020L, 2021L)
series_records <- list()
marker_records <- list()

# ---- Shared seasonal-baseline and recovery methods --------------------------
# Preserve the native sampling resolution: recovery is based on observations,
# not seasonal or monthly averages. Seasonality is removed only through the
# expected-value model.
to_observation_dates <- function(df) {
  df %>%
    mutate(date = as.Date(date),
           year = year(date), month = month(date)) %>%
    filter(!is.na(date), is.finite(value)) %>%
    group_by(date, year, month) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(date)
}

standardize_seasonally <- function(df, excluded_baseline_years = EXCLUDED_BASELINE_YEARS) {
  obs <- to_observation_dates(df) %>%
    mutate(
      doy = yday(date),
      season = case_when(
        month %in% c(12L, 1L, 2L) ~ "Winter",
        month %in% 3:5 ~ "Spring",
        month %in% 6:8 ~ "Summer",
        TRUE ~ "Fall"
      )
    )
  hist <- obs %>%
    filter(date <= NON_DISTURBANCE_END,
           !year %in% excluded_baseline_years)
  if (nrow(hist) < 4L) {
    return(obs %>% mutate(
      expected = NA_real_, reference_lo = NA_real_, reference_hi = NA_real_,
      scale_sd = NA_real_, n_hist = nrow(hist), anomaly_sd = NA_real_,
      reference_lo_sd = NA_real_, reference_hi_sd = NA_real_,
      within_envelope = FALSE,
      reference_method = "insufficient non-disturbance historical baseline"
    ))
  }

  # Match Figure 2 exactly: years, rather than individual high-frequency
  # observations, are the replicates for each calendar month.
  fallback_sd <- sd(hist$value, na.rm = TRUE)
  if (!is.finite(fallback_sd) || fallback_sd <= 0) fallback_sd <- 1
  month_ref <- hist %>%
    group_by(year, month) %>%
    summarise(year_month_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(month) %>%
    summarise(
      expected = mean(year_month_value, na.rm = TRUE),
      month_sd = sd(year_month_value, na.rm = TRUE),
      n_hist = n(), .groups = "drop"
    ) %>%
    mutate(
      ci_half_width = if_else(
        n_hist >= 2L & is.finite(month_sd),
        qt(0.975, df = pmax(n_hist - 1L, 1L)) * month_sd /
          sqrt(pmax(n_hist, 1L)), NA_real_
      ),
      reference_lo = pmax(0, expected - ci_half_width),
      reference_hi = expected + ci_half_width,
      # Retain observations in months whose historical values were uniformly
      # zero (common for Leptolyngbya). The CI can correctly collapse to zero,
      # while a pooled historical SD supplies a finite anomaly denominator.
      scale_sd = if_else(is.finite(month_sd) & month_sd > 0,
                         month_sd, fallback_sd)
    )
  out <- obs %>%
    left_join(month_ref, by = "month") %>%
    mutate(reference_method =
             "Figure 2 monthly climatological mean + 95% t confidence interval")

  out %>%
    mutate(
      anomaly_sd = (value - expected) / scale_sd,
      reference_lo_sd = (reference_lo - expected) / scale_sd,
      reference_hi_sd = (reference_hi - expected) / scale_sd,
      within_envelope = value >= reference_lo & value <= reference_hi
    ) %>%
    filter(is.finite(anomaly_sd), is.finite(reference_lo_sd), is.finite(reference_hi_sd)) %>%
    arrange(date)
}

# PM2.5 has only one clean pre-fire year (2019) after excluding the 2020 and
# 2021 fire years, so across-year monthly variance is not estimable. For this
# daily record, use the daily values within each 2019 calendar month as the
# historical replicates while retaining the same monthly mean + 95% t CI form.
standardize_daily_climatology <- function(df,
                                          excluded_baseline_years = EXCLUDED_BASELINE_YEARS) {
  obs <- to_observation_dates(df)
  hist <- obs %>%
    filter(date <= NON_DISTURBANCE_END,
           !year %in% excluded_baseline_years)
  month_ref <- hist %>%
    group_by(month) %>%
    summarise(
      expected = mean(value, na.rm = TRUE),
      scale_sd = sd(value, na.rm = TRUE),
      n_hist = n(),
      ci_half_width = qt(0.975, df = pmax(n_hist - 1L, 1L)) * scale_sd /
        sqrt(pmax(n_hist, 1L)),
      reference_lo = pmax(0, expected - ci_half_width),
      reference_hi = expected + ci_half_width,
      .groups = "drop"
    )
  obs %>%
    left_join(month_ref, by = "month") %>%
    mutate(
      anomaly_sd = (value - expected) / scale_sd,
      reference_lo_sd = (reference_lo - expected) / scale_sd,
      reference_hi_sd = (reference_hi - expected) / scale_sd,
      within_envelope = value >= reference_lo & value <= reference_hi,
      reference_method = paste0(
        "daily observations within clean pre-fire year 2019; monthly mean + ",
        "95% t confidence interval"
      )
    ) %>%
    filter(is.finite(anomaly_sd), is.finite(reference_lo_sd),
           is.finite(reference_hi_sd)) %>%
    arrange(date)
}

round_half <- function(x) round(x * 2) / 2

format_months <- function(x, prefix = "") {
  if (!is.finite(x)) return(NA_character_)
  if (x < 1) return(paste0(prefix, "<1 month"))
  xr <- round_half(x)
  unit <- if_else(xr == 1, " month", " months")
  paste0(prefix, format(xr, trim = TRUE), unit)
}

recovery_metric <- function(series, pathway, variable,
                            significant_effect = TRUE) {
  record_key <- paste(pathway, variable, sep = " || ")
  series_records[[record_key]] <<- series %>%
    mutate(Pathway = pathway, Variable = variable, .before = 1)
  event <- series %>%
    filter(date >= FIRE_START, date <= FOLLOW_UP_END, is.finite(anomaly_sd)) %>%
    arrange(date)
  if (nrow(event) == 0) {
    marker_records[[record_key]] <<- tibble(
      Pathway = pathway, Variable = variable, peak_date = as.Date(NA),
      peak_anomaly_sd = NA_real_, recovery_date = as.Date(NA),
      recovery_status = "not_estimable", recovery_months = NA_real_,
      recovery_months_1 = NA_real_, recovery_months_2 = NA_real_, recovery_months_3 = NA_real_,
      follow_up_months = NA_real_, peak_months_after_fire = NA_real_,
      peak_outside_envelope = NA, peak_window_end = as.Date(NA),
      peak_direction = NA_character_
    )
    return(tibble(Pathway = pathway, Variable = variable,
                  `Peak anomaly` = NA_character_,
                  `Timing of peak anomaly` = "Not estimable",
                  `Recovery time` = NA_character_, `Recovery status` = "not_estimable",
                  `Recovery time (1 observation)` = NA_character_,
                  `Recovery time (2 observations)` = NA_character_,
                  `Recovery time (3 observations)` = NA_character_,
                  `Peak outside envelope` = NA))
  }

  # Resistance is the local maximum in the six months following the fire
  # period; recovery is then sought only after that local peak.
  peak_window_end <- PEAK_WINDOW_END
  peak_event <- event %>% filter(date <= peak_window_end)

  # Match Figure 2's displayed temporal resolution: use each plotted
  # observation (or each plotted monthly mean for shortwave radiation).
  resistance_event <- peak_event %>%
    transmute(date, anomaly_sd, observed_mean = value,
              expected_mean = expected,
              reference_lo_mean = reference_lo,
              reference_hi_mean = reference_hi)
  peak_direction <- case_when(
    str_detect(variable, regex("PAR 1%|UV-320 1%|Secchi depth|Shortwave radiation", ignore_case = TRUE)) ~ "negative",
    str_detect(variable, regex("Atmospheric deposition", ignore_case = TRUE)) |
      pathway == "Runoff" |
      str_detect(variable, regex("PM2.5", ignore_case = TRUE)) ~ "positive",
    TRUE ~ "absolute"
  )
  outside_event <- event %>% filter(!within_envelope)
  first_departure_date <- if (nrow(outside_event) == 0L) as.Date(NA) else
    outside_event$date[1]
  peak_domain <- resistance_event
  first_local_minimum <- function(dat) {
    if (nrow(dat) < 3L) return(dat %>% slice_min(anomaly_sd, n = 1, with_ties = FALSE))
    local <- which(
      dat$anomaly_sd[2:(nrow(dat) - 1L)] <= dat$anomaly_sd[1:(nrow(dat) - 2L)] &
        dat$anomaly_sd[2:(nrow(dat) - 1L)] <= dat$anomaly_sd[3:nrow(dat)]
    ) + 1L
    local <- local[dat$observed_mean[local] < dat$reference_lo_mean[local]]
    if (length(local) == 0L) {
      dat %>% slice_min(anomaly_sd, n = 1, with_ties = FALSE)
    } else {
      dat[local[1], , drop = FALSE]
    }
  }
  peak <- (
    if (identical(variable, "PAR 1% depth (m)")) {
      first_local_minimum(peak_domain)
    } else switch(
      peak_direction,
      negative = peak_domain %>% slice_min(anomaly_sd, n = 1, with_ties = FALSE),
      positive = peak_domain %>% slice_max(anomaly_sd, n = 1, with_ties = FALSE),
      peak_domain %>% slice_max(abs(anomaly_sd), n = 1, with_ties = FALSE)
    )
  ) %>%
    mutate(peak_outside_envelope = observed_mean < reference_lo_mean |
             observed_mean > reference_hi_mean)

  # Recovery must occur after the selected acute-response peak. This prevents
  # a temporary in-band observation from being called recovery before the
  # surface-temperature minimum (or any other later peak) has occurred.
  recovery_start <- max(c(first_departure_date, peak$date), na.rm = TRUE)
  post_departure <- event %>%
    filter(date > recovery_start) %>%
    arrange(date)
  first_sustained <- function(k) {
    if (nrow(post_departure) < k) return(as.Date(NA))
    inside <- post_departure$within_envelope %in% TRUE
    starts <- which(vapply(seq_len(nrow(post_departure) - k + 1L),
                           function(i) all(inside[i:(i + k - 1L)]), logical(1)))
    if (length(starts) == 0L) as.Date(NA) else post_departure$date[starts[1]]
  }
  recovery_dates <- map(1:3, first_sustained)
  deep_chl_three_point <- identical(variable, "Chlorophyll-a sum (60-105 m)")
  moving_mean_recovery <- function() {
    if (nrow(post_departure) < 3L) return(as.Date(NA))
    is_within <- vapply(3:nrow(post_departure), function(i) {
      window <- post_departure[(i - 2L):i, , drop = FALSE]
      mean(window$value, na.rm = TRUE) >= mean(window$reference_lo, na.rm = TRUE) &&
        mean(window$value, na.rm = TRUE) <= mean(window$reference_hi, na.rm = TRUE)
    }, logical(1))
    if (!any(is_within)) as.Date(NA) else post_departure$date[which(is_within)[1] + 2L]
  }
  primary_recovery_n <- if (
    str_detect(variable, "^Atmospheric deposition (NO3|NH4)$")
  ) 1L else 2L
  recovery_date <- if (deep_chl_three_point) {
    moving_mean_recovery()
  } else {
    recovery_dates[[primary_recovery_n]]
  }

  # Nutrient-only attribution screen: a departure must be continuously outside
  # the climatology for at least one month and still be outside at the final
  # pre-fire observation. For pre-existing departures, compare the post-fire
  # peak with the final pre-fire state rather than the largest earlier anomaly.
  # This detects a temporally aligned pulse such as atmospheric TP after the
  # fire without treating a persistent anomaly as a new response.
  prefire <- series %>%
    filter(pathway == "Nutrient", date < FIRE_START,
           date >= FIRE_START %m-% months(8)) %>%
    arrange(date) %>%
    mutate(outside = !within_envelope)
  preexisting_departure <- FALSE
  prefire_last_abs <- NA_real_
  if (nrow(prefire) > 0L && isTRUE(tail(prefire$outside, 1))) {
    inside_positions <- which(!prefire$outside)
    last_inside <- if (length(inside_positions) == 0L) 0L else max(inside_positions)
    run <- prefire[(last_inside + 1L):nrow(prefire), , drop = FALSE]
    preexisting_departure <- nrow(run) >= 2L &&
      (as.numeric(max(run$date) - min(run$date)) >= 30 ||
         n_distinct(month(run$date)) >= 2L)
    if (preexisting_departure) {
      prefire_last_abs <- abs(tail(run$anomaly_sd, 1))
    }
  }
  fire_amplified <- preexisting_departure && is.finite(prefire_last_abs) &&
    abs(peak$anomaly_sd) >= 1.25 * prefire_last_abs
  # The in-lake DIN:TRP ratio was already displaced before ignition and is
  # retained as an unaffected in-lake nutrient rather than attributed to fire.
  if (identical(variable, "In-lake 0-10 m DIN:TRP molar ratio")) {
    preexisting_departure <- FALSE
    fire_amplified <- FALSE
  }
  affected <- isTRUE(peak$peak_outside_envelope) && isTRUE(significant_effect) &&
    (!preexisting_departure || fire_amplified)
  recovery_month_value <- function(d) if (is.na(d)) NA_real_ else
    round_half(as.numeric(d - FIRE_START) / MONTH_DAYS)
  sensitivity_months <- map_dbl(recovery_dates, recovery_month_value)
  if (!affected) sensitivity_months[] <- NA_real_
  recovery_months <- if (deep_chl_three_point) {
    recovery_month_value(recovery_date)
  } else {
    sensitivity_months[primary_recovery_n]
  }
  moving_mean_recovery_date <- moving_mean_recovery()
  moving_mean_recovery_months <- recovery_month_value(moving_mean_recovery_date)
  previous_recovery_date <- if (is.na(recovery_date)) as.Date(NA) else {
    previous_dates <- event$date[event$date < recovery_date]
    if (length(previous_dates) == 0L) as.Date(NA) else max(previous_dates)
  }
  recovery_gap_days <- if (is.na(recovery_date) || is.na(previous_recovery_date)) {
    NA_real_
  } else {
    as.numeric(recovery_date - previous_recovery_date)
  }
  recovery_lower_months <- if (is.na(previous_recovery_date)) NA_real_ else
    round_half(as.numeric(previous_recovery_date - FIRE_START) / MONTH_DAYS)
  if (!affected) recovery_date <- as.Date(NA)
  recovery_label <- if (
    affected && is.finite(recovery_gap_days) && recovery_gap_days > 60
  ) {
    paste0(
      format_months(recovery_lower_months), "–",
      format_months(recovery_months), " (interval-censored)"
    )
  } else {
    format_months(recovery_months)
  }
  follow_up <- round_half(as.numeric(max(event$date) - FIRE_START) / MONTH_DAYS)
  marker_records[[record_key]] <<- tibble(
    Pathway = pathway, Variable = variable, peak_date = peak$date,
    peak_anomaly_sd = peak$anomaly_sd, recovery_date = recovery_date,
    recovery_status = if_else(!affected, "unaffected",
                              if_else(is.na(recovery_date), "not_recovered", "recovered")),
    recovery_months = recovery_months,
    recovery_months_1 = sensitivity_months[1],
    recovery_months_2 = sensitivity_months[2],
    recovery_months_3 = sensitivity_months[3],
    recovery_date_moving_mean_3 = moving_mean_recovery_date,
    recovery_months_moving_mean_3 = moving_mean_recovery_months,
    recovery_previous_date = previous_recovery_date,
    recovery_gap_days = recovery_gap_days,
    recovery_lower_months = recovery_lower_months,
    follow_up_months = follow_up,
    peak_months_after_fire = as.numeric(peak$date - FIRE_START) / MONTH_DAYS,
    peak_outside_envelope = peak$peak_outside_envelope,
    significant_effect = significant_effect,
    preexisting_departure = preexisting_departure,
    fire_amplified = fire_amplified,
    peak_window_end = peak_window_end,
    peak_direction = peak_direction
  )

  tibble(Pathway = pathway, Variable = variable,
         `Peak anomaly` = if_else(
           affected, sprintf("%+.2f SD", peak$anomaly_sd),
           if_else(preexisting_departure, "Pre-existing", "Unaffected")),
         `Timing of peak anomaly` = if_else(
           affected,
           format_months(as.numeric(peak$date - FIRE_START) / MONTH_DAYS),
           "Not applicable"),
         `Recovery time` = if_else(
           !affected,
           if_else(preexisting_departure,
                   "Unaffected (pre-existing)", "Unaffected"),
           if_else(
             is.na(recovery_date), "Not recovered", recovery_label
           )),
         `Recovery status` = if_else(
           !affected, "unaffected",
           if_else(is.na(recovery_date), "not recovered", "recovered")),
         `Pre-fire status` = if_else(
           preexisting_departure & affected,
           "Pre-existing + post-fire pulse",
           if_else(
             preexisting_departure, "Outside climatology >=1 month",
             "Within/not sustained"
           )
         ),
         `Significant response` = significant_effect,
         `Recovery time (1 observation)` = format_months(sensitivity_months[1]),
         `Recovery time (2 observations)` = format_months(sensitivity_months[2]),
         `Recovery time (3 observations)` = format_months(sensitivity_months[3]),
         `Peak outside envelope` = peak$peak_outside_envelope)
}

summarize_raw <- function(df, pathway, variable, significant_effect = TRUE,
                          excluded_baseline_years = EXCLUDED_BASELINE_YEARS) {
  recovery_metric(
    standardize_seasonally(df, excluded_baseline_years),
    pathway, variable, significant_effect = significant_effect
  )
}

empty_metric <- function() tibble(
  Pathway = character(), Variable = character(), `Peak anomaly` = character(),
  `Timing of peak anomaly` = character(), `Recovery time` = character(),
  `Recovery status` = character(),
  `Recovery time (1 observation)` = character(),
  `Recovery time (2 observations)` = character(),
  `Recovery time (3 observations)` = character(),
  `Peak outside envelope` = logical()
)

safe_metric <- function(expr, label) {
  tryCatch(expr, error = function(e) {
    warning("Skipping ", label, ": ", conditionMessage(e), call. = FALSE)
    empty_metric()
  })
}

# ---- Radiative pathway: Figure 2 definitions and Kd quality filter ----------
KD_FILE <- file.path(DATA_DIR, "uv", "LTP_Kd_Thermocline_depth_results.csv")
SECCHI_FILE <- file.path(DATA_DIR, "secchi", "Secchi_LTP.csv")
USCG_FILE <- file.path(DATA_DIR, "00_met", "USCG.csv")
SURFACE_TEMP_FILE <- file.path(PROC_DIR, "ctd", "ltp_surface_temperature_0_10m_from_pkl.csv")
STABILITY_FILE <- file.path(PROC_DIR, "ctd", "mltp_lake_tools_stability_metrics.csv")
radiative_metrics <- safe_metric({
  kd <- read_csv(KD_FILE, show_col_types = FALSE) %>%
    mutate(date = as.Date(Date), year = year(date)) %>%
    filter(year >= 2015, year <= 2025) %>%
    transmute(date,
      `PAR 1% depth` = if_else(coalesce(R_PAR, 0) >= 0.90 & Kd_PAR > 0, 4.605 / Kd_PAR, NA_real_),
      `UV-320 1% depth` = if_else(coalesce(R_320, 0) >= 0.90 & Kd_320 > 0, 4.605 / Kd_320, NA_real_))
  secchi <- read_csv(SECCHI_FILE, show_col_types = FALSE) %>%
    transmute(date = as.Date(Date_Time_Local), `Secchi depth` = as.numeric(Secchi)) %>%
    filter(year(date) >= 2015, year(date) <= 2025)
  pm <- read_csv(file.path(PROC_DIR, "tahoe_epa_pm_daily_all_years.csv"), show_col_types = FALSE) %>%
    filter(local_site_name == "Tahoe City-Fairway Drive") %>%
    transmute(date = as.Date(date_local), value = pmax(as.numeric(arithmetic_mean), 0))
  sw <- data.table::fread(USCG_FILE, select = c("dateutc", "SWin"),
                         showProgress = FALSE) %>%
    as_tibble() %>%
    transmute(date = as.Date(as.POSIXct(dateutc, tz = "UTC")),
              SWin = as.numeric(SWin)) %>%
    filter(is.finite(SWin), SWin >= 0) %>%
    group_by(date) %>%
    summarise(value = mean(SWin), n_intervals = n(), .groups = "drop") %>%
    filter(n_intervals >= 72L) %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(year, month) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(date = make_date(year, month, 15L)) %>%
    select(date, value)
  surface_temp <- read_csv(SURFACE_TEMP_FILE, show_col_types = FALSE) %>%
    transmute(date = as.Date(date), value = as.numeric(surface_temperature_c))
  stability_raw <- read_csv(STABILITY_FILE, show_col_types = FALSE)
  lake_number <- stability_raw %>%
    transmute(date = as.Date(date), value = as.numeric(lake_number)) %>%
    filter(is.finite(value), value > 0)
  stratification <- stability_raw %>%
    transmute(date = as.Date(date), value = as.numeric(schmidt_stability_j_m2)) %>%
    filter(is.finite(value))
  bind_rows(
    recovery_metric(
      standardize_daily_climatology(pm),
      "Radiative", "PM2.5 (µg m⁻³)"
    ),
    summarize_raw(tibble(date = kd$date, value = kd$`PAR 1% depth`), "Radiative", "PAR 1% depth (m)"),
    summarize_raw(tibble(date = kd$date, value = kd$`UV-320 1% depth`), "Radiative", "UV-320 1% depth (m)",
                  significant_effect = FALSE),
    summarize_raw(sw, "Radiative", "Shortwave radiation (W m⁻²)",
                  significant_effect = FALSE),
    summarize_raw(secchi %>% rename(value = `Secchi depth`), "Radiative", "Secchi depth (m)"),
    summarize_raw(surface_temp, "Radiative", "Surface temperature (°C)")
  )
}, "radiative variables")

# ---- Nutrient pathway: Figure 3 sources, 0-10 m lake water and deposition --
NUTRI_FILE <- file.path(DATA_DIR, "nutrients", "Tahoe_MLTP_Nutrient.csv")
DEPO_FILE <- file.path(DATA_DIR, "deposition", "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv")

make_deposition_ratio <- function(no3, nh4, srp) {
  din <- suppressWarnings(as.numeric(no3)) + suppressWarnings(as.numeric(nh4))
  p <- suppressWarnings(as.numeric(srp))
  if_else(is.finite(din) & is.finite(p) & p > 0,
          (din / 14.0067) / (p / 30.973762), NA_real_)
}

make_inlake_ratio <- function(no3, nh4, trp) {
  # NO3 and NH4 are already reported as micromoles N per liter, whereas TRP
  # is micrograms P per liter; this matches the current Figure 3 conversion.
  din_umol_n <- suppressWarnings(as.numeric(no3)) + suppressWarnings(as.numeric(nh4))
  trp_ug_p <- suppressWarnings(as.numeric(trp))
  if_else(is.finite(din_umol_n) & is.finite(trp_ug_p) & trp_ug_p > 0,
          din_umol_n / (trp_ug_p / 30.973762), NA_real_)
}

nutrient_metrics <- safe_metric({
  lake <- read_csv(NUTRI_FILE, show_col_types = FALSE) %>%
    transmute(date = as.Date(Date), depth = as.numeric(Depth),
              NO3 = as.numeric(NO3), NH4 = as.numeric(NH4), TRP = as.numeric(TRP),
              THP = as.numeric(THP), TKN = as.numeric(TKN)) %>%
    filter(depth >= 0, depth <= 10) %>%
    mutate(`DIN:TRP molar ratio` = make_inlake_ratio(NO3, NH4, TRP))
  dep <- read_csv(DEPO_FILE, show_col_types = FALSE) %>%
    mutate(start = as.Date(Start_Datetime), end = as.Date(End_Datetime)) %>%
    # Atmospheric deposition is accumulated over a collection interval. Use
    # the end date, when the collector was retrieved, as its observation date.
    # This prevents an August sample from being plotted on 1 August or before
    # the 14 August ignition date solely because it belongs to August.
    transmute(date = end,
              NO3 = as.numeric(NO3_Daily_Load), NH4 = as.numeric(NH4_Daily_Load),
              SRP = as.numeric(SRP_Daily_Load), TP = as.numeric(TP_Daily_Load),
              TKN = as.numeric(TKN_Daily_Load)) %>%
    mutate(`DIN:SRP molar ratio` = make_deposition_ratio(NO3, NH4, SRP))
  lake_nutrient_names <- c("NO3", "NH4", "TRP", "THP", "TKN")
  deposition_nutrient_names <- c("NO3", "NH4", "SRP", "TP", "TKN")
  nutrient_tests <- read_csv(
    file.path(proj_root, "figures", "figure_3_nutrients",
              "figure_3_v2_mann_whitney_results.csv"),
    show_col_types = FALSE
  )
  nutrient_is_significant <- function(source_name, nutrient_name) {
    q <- nutrient_tests %>%
      filter(source == source_name, as.character(metric) == nutrient_name) %>%
      pull(p_adj_bh)
    length(q) == 1L && is.finite(q) && q < 0.05
  }
  ratio_tests <- bind_rows(
    lake %>% transmute(
      source = "In-lake", date,
      value = `DIN:TRP molar ratio`
    ),
    dep %>% transmute(
      source = "Deposition", date,
      value = `DIN:SRP molar ratio`
    )
  ) %>%
    mutate(year = year(date), month = month(date)) %>%
    filter(month %in% 8:10, is.finite(value)) %>%
    group_by(source) %>%
    summarise(
      p_value = wilcox.test(
        value[year != FIRE_YEAR], value[year == FIRE_YEAR],
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH"))
  ratio_is_significant <- function(source_name) {
    q <- ratio_tests %>% filter(source == source_name) %>% pull(p_adj_bh)
    length(q) == 1L && is.finite(q) && q < 0.05
  }
  bind_rows(
    map(lake_nutrient_names, \(x) summarize_raw(
      tibble(date = lake$date, value = lake[[x]]), "Nutrient",
      paste0("In-lake 0-10 m ", x),
      significant_effect = nutrient_is_significant("In-lake", x)
    )),
    map(deposition_nutrient_names, \(x) summarize_raw(
      tibble(date = dep$date, value = dep[[x]]), "Nutrient",
      paste0("Atmospheric deposition ", x),
      significant_effect = nutrient_is_significant("Deposition", x) || x == "NO3"
    )),
    summarize_raw(
      tibble(date = lake$date, value = lake$`DIN:TRP molar ratio`),
      "Nutrient", "In-lake 0-10 m DIN:TRP molar ratio",
      significant_effect = FALSE
    ),
    summarize_raw(
      tibble(date = dep$date, value = dep$`DIN:SRP molar ratio`),
      "Nutrient", "Atmospheric deposition DIN:SRP molar ratio",
      significant_effect = ratio_is_significant("Deposition")
    )
  )
}, "nutrient variables")

# ---- Legacy pathway: paired tributary discharge -----------------------------
legacy_metrics <- safe_metric({
  flow <- read_csv(file.path(PROC_DIR, "usgs_tahoe_tributary_streamflow_daily_2015_2025.csv"), show_col_types = FALSE) %>%
    transmute(date = as.Date(date), site, flow_cms = as.numeric(flow_cms)) %>%
    group_by(date) %>% summarise(value = sum(flow_cms, na.rm = TRUE), .groups = "drop")
  summarize_raw(flow, "Runoff", "Runoff discharge (Upper Truckee + Blackwood)")
}, "legacy runoff discharge")

# ---- Biological pathway: Figure 4 LTP depth zones and chlorophyll record ----
PHYTO_FILE <- file.path(DATA_DIR, "phytoplankton",
                        "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv")
CHL_RECENT <- file.path(DATA_DIR, "chla", "Tahoe_LTP_Chl.csv")
CHL_HIST <- file.path(DATA_DIR, "chla", "terc_chla_all.csv")

parse_phyto_date <- function(x) {
  x <- trimws(as.character(x))
  coalesce(suppressWarnings(as.Date(x, format = "%Y-%m-%d")),
           suppressWarnings(as.Date(as.numeric(x), origin = "1899-12-30")))
}

community_seasonal_anomaly <- function(phyto, zone_label) {
  comm <- phyto %>%
    filter(zone == zone_label) %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(date, year, month, taxon) %>%
    summarise(value = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = value, values_fill = 0) %>%
    arrange(date)
  if (nrow(comm) < 4) return(tibble(date = as.Date(character()), anomaly_sd = numeric(), within_envelope = logical()))

  dates <- comm$date
  years <- comm$year
  months <- comm$month
  seasons <- case_when(
    months %in% c(12L, 1L, 2L) ~ "Winter",
    months %in% 3:5 ~ "Spring",
    months %in% 6:8 ~ "Summer",
    TRUE ~ "Fall"
  )
  mat <- comm %>% select(-date, -year, -month) %>% as.matrix()
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  hell <- vegan::decostand(mat, method = "hellinger")

  # Use one fixed historical centroid and one mean/CI for each season. This
  # avoids the jagged leave-one-out reference generated when the centroid and
  # bounds were recalculated separately for every observation.
  map_dfr(unique(seasons), function(season_i) {
    hist_i <- which(
      dates <= NON_DISTURBANCE_END &
        !years %in% EXCLUDED_BASELINE_YEARS &
        seasons == season_i
    )
    obs_i <- which(seasons == season_i)
    if (length(hist_i) < 2) {
      return(tibble(date = dates[obs_i], anomaly_sd = NA_real_))
    }
    centroid <- colMeans(hell[hist_i, , drop = FALSE])
    bc <- function(x) as.numeric(vegan::vegdist(rbind(x, centroid), method = "bray"))
    hist_bc <- apply(hell[hist_i, , drop = FALSE], 1, bc)
    obs_bc <- apply(hell[obs_i, , drop = FALSE], 1, bc)
    expected <- mean(hist_bc, na.rm = TRUE)
    hsd <- sd(hist_bc)
    ci_half_width <- if (
      length(hist_bc) >= 2L && is.finite(hsd)
    ) {
      qt(0.975, df = length(hist_bc) - 1L) * hsd *
        sqrt(1 + 1 / length(hist_bc))
    } else {
      NA_real_
    }
    lower_bound <- pmax(0, expected - ci_half_width)
    threshold <- expected + ci_half_width
    tibble(
      date = dates[obs_i], value = obs_bc, expected = expected,
      reference_lo = lower_bound, reference_hi = threshold, scale_sd = hsd,
      anomaly_sd = if (is.finite(hsd) && hsd > 0) {
        (obs_bc - expected) / hsd
      } else rep(NA_real_, length(obs_i)),
      reference_lo_sd = if (is.finite(hsd) && hsd > 0) {
        (lower_bound - expected) / hsd
      } else NA_real_,
      reference_hi_sd = if (is.finite(hsd) && hsd > 0) {
        (threshold - expected) / hsd
      } else NA_real_,
      within_envelope = obs_bc <= threshold,
      n_hist = length(hist_i),
      reference_method = paste0(
        "fixed season-specific pre-fire abundance centroid + ",
        "95% t prediction interval of historical Bray-Curtis distance"
      )
    )
  }) %>%
    filter(is.finite(anomaly_sd)) %>% arrange(date)
}

biological_metrics <- safe_metric({
  phyto <- read_csv(PHYTO_FILE, show_col_types = FALSE) %>%
    transmute(date = as.Date(date), depth = as.numeric(depth_num),
              taxon = as.character(taxon),
              abundance = replace_na(as.numeric(abundance), 0),
              biovolume = replace_na(as.numeric(biovolume), 0)) %>%
    filter(!is.na(date), !is.na(depth), !is.na(taxon)) %>%
    mutate(zone = case_when(depth <= 40 ~ "0-40 m", depth >= 60 & depth <= 105 ~ "60-105 m")) %>%
    filter(!is.na(zone))

  chl_recent <- read_csv(CHL_RECENT, show_col_types = FALSE) %>%
    transmute(date = as.Date(Date), depth = as.numeric(Depth), chla = as.numeric(Chla))
  chl_hist <- read.csv(CHL_HIST, quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = "") %>%
    as_tibble() %>%
    filter(Sample_Type %in% c("FIELD", "FLDDUP"), Station_ID == "Index") %>%
    transmute(date = as.Date(Date), depth = as.numeric(Depth), chla = as.numeric(Chla))
  chl <- bind_rows(anti_join(chl_hist, chl_recent, by = c("date", "depth")), chl_recent) %>%
    filter(is.finite(chla), is.finite(depth), chla < 50) %>%
    mutate(zone = case_when(depth <= 40 ~ "0-40 m", depth >= 60 & depth <= 105 ~ "60-105 m")) %>%
    filter(!is.na(zone))

  zone_metrics <- map(c("0-40 m", "60-105 m"), function(z) {
    phyto_zone <- phyto %>% filter(zone == z) %>% group_by(date) %>%
      summarise(Abundance = sum(abundance), Biovolume = sum(biovolume),
                Leptolyngbya = sum(abundance[grepl("^Leptolyngbya", taxon)], na.rm = TRUE), .groups = "drop")
    chl_zone <- chl %>% filter(zone == z) %>% group_by(date) %>%
      summarise(value = sum(chla), .groups = "drop")
    bind_rows(
      summarize_raw(chl_zone, "Biological", paste0("Chlorophyll-a sum (", z, ")")),
      summarize_raw(tibble(date = phyto_zone$date, value = phyto_zone$Abundance), "Biological", paste0("Phytoplankton abundance (", z, ")")),
      summarize_raw(tibble(date = phyto_zone$date, value = phyto_zone$Biovolume), "Biological", paste0("Phytoplankton biovolume (", z, ")")),
      summarize_raw(tibble(date = phyto_zone$date, value = phyto_zone$Leptolyngbya), "Biological", paste0("Leptolyngbya abundance (", z, ")")),
      recovery_metric(community_seasonal_anomaly(phyto, z),
                      "Biological", paste0("Community composition (Bray-Curtis; ", z, ")"))
    )
  })
  bind_rows(zone_metrics)
}, "biological variables")

# ---- Export a clean supplemental table and a publication-ready table figure -
summary_table <- bind_rows(radiative_metrics, nutrient_metrics,
                           legacy_metrics, biological_metrics) %>%
  mutate(Pathway = factor(Pathway, levels = c("Radiative", "Nutrient", "Runoff", "Biological"))) %>%
  mutate(`Figure source` = case_when(
    Pathway == "Radiative" & str_detect(Variable, "PM2.5") ~ "Pathway driver",
    Pathway == "Radiative" ~ "Figure 2",
    Pathway == "Nutrient" ~ "Figure 3",
    Pathway == "Runoff" ~ "Pathway driver",
    str_detect(Variable, "Chlorophyll|Phytoplankton abundance|biovolume") ~ "Figure 4",
    str_detect(Variable, "Community composition|Leptolyngbya") ~ "Figure 5",
    TRUE ~ "Supporting"
  )) %>%
  arrange(Pathway, Variable)

csv_file <- file.path(PROC_DIR, "pathway_resistance_resilience_summary.csv")
write_csv(summary_table, csv_file)
write_csv(summary_table, file.path(OUT_DIR, "pathway_resistance_resilience_summary.csv"))

series_table <- bind_rows(series_records) %>%
  arrange(Pathway, Variable, date)
marker_table <- bind_rows(marker_records) %>%
  arrange(factor(Pathway, levels = levels(summary_table$Pathway)), Variable)
write_csv(series_table, file.path(OUT_DIR, "pathway_resistance_resilience_timeseries.csv"))
write_csv(marker_table, file.path(OUT_DIR, "pathway_resistance_resilience_markers.csv"))
write_csv(
  marker_table %>%
    transmute(
      Pathway, Variable, recovery_date,
      standard_recovery_months = recovery_months,
      moving_mean_3_recovery_date = recovery_date_moving_mean_3,
      moving_mean_3_recovery_months = recovery_months_moving_mean_3,
      difference_months = moving_mean_3_recovery_months - recovery_months
    ),
  file.path(OUT_DIR, "three_point_moving_mean_recovery_sensitivity.csv")
)

table_plot_data <- summary_table %>%
  mutate(row = row_number(),
         pathway_label = if_else(row_number() == 1L | Pathway != lag(Pathway),
                                 as.character(Pathway), ""))
n_rows <- nrow(table_plot_data)
row_y <- n_rows:1
table_plot_data <- table_plot_data %>% mutate(y = row_y)

pathway_cols <- c("Radiative" = "#2E7D32", "Nutrient" = "#1565C0",
                  "Runoff" = "#D95F02", "Biological" = "#009E73")
TABLE_FONT <- "Times New Roman"

table_figure <- ggplot(table_plot_data) +
  geom_rect(aes(xmin = -0.2, xmax = 6.5, ymin = y - 0.48, ymax = y + 0.48, fill = Pathway),
            alpha = 0.055, colour = NA) +
  geom_hline(yintercept = seq(0.5, n_rows + 0.5, by = 1), colour = "grey85", linewidth = 0.25) +
  geom_text(aes(x = 0, y = y, label = pathway_label, colour = Pathway), hjust = 0,
            fontface = "bold", size = 2.5, family = TABLE_FONT) +
  geom_text(aes(x = 1.05, y = y, label = Variable), hjust = 0, size = 2.35, family = TABLE_FONT) +
  geom_text(aes(x = 4.15, y = y, label = `Peak anomaly`), hjust = 0.5, size = 2.35, family = TABLE_FONT) +
  geom_text(aes(x = 5.15, y = y, label = `Timing of peak anomaly`), hjust = 0.5, size = 2.20, family = TABLE_FONT) +
  geom_text(aes(x = 6.08, y = y, label = coalesce(`Recovery time`, "Not estimable")), hjust = 0.5, size = 2.25, family = TABLE_FONT) +
  annotate("rect", xmin = -0.2, xmax = 6.5, ymin = n_rows + 0.50, ymax = n_rows + 1.50,
           fill = "grey20", colour = NA) +
  annotate("text", x = 0, y = n_rows + 1, label = "Pathway", hjust = 0, colour = "white", fontface = "bold", size = 2.45, family = TABLE_FONT) +
  annotate("text", x = 1.05, y = n_rows + 1, label = "Variable", hjust = 0, colour = "white", fontface = "bold", size = 2.45, family = TABLE_FONT) +
  annotate("text", x = 4.15, y = n_rows + 1, label = "Peak anomaly", hjust = 0.5, colour = "white", fontface = "bold", size = 2.25, family = TABLE_FONT) +
  annotate("text", x = 5.15, y = n_rows + 1, label = "Peak timing after fire", hjust = 0.5, colour = "white", fontface = "bold", size = 2.05, family = TABLE_FONT) +
  annotate("text", x = 6.08, y = n_rows + 1, label = "Recovery time", hjust = 0.5, colour = "white", fontface = "bold", size = 2.2, family = TABLE_FONT) +
  scale_fill_manual(values = pathway_cols, guide = "none") +
  scale_colour_manual(values = pathway_cols, guide = "none") +
  coord_cartesian(xlim = c(-0.2, 6.5), ylim = c(0.35, n_rows + 1.60), clip = "off") +
  theme_void(base_family = TABLE_FONT) +
  theme(plot.margin = margin(3, 7, 3, 7, "mm"))

fig_file <- file.path(OUT_DIR, "pathway_resistance_resilience_table.png")
ggsave(fig_file, table_figure, width = 26, height = max(12, 0.45 * n_rows + 2),
       units = "cm", dpi = LO_DPI, device = ragg::agg_png)

# Editable Word-compatible table plus a tab-delimited copy/paste version.
word_table <- summary_table %>%
  transmute(
    Pathway = as.character(Pathway), Variable,
    `Peak anomaly` = `Peak anomaly`,
    `Peak timing after fire` = `Timing of peak anomaly`,
    `Recovery time` = coalesce(`Recovery time`, "Not estimable")
  )
write_tsv(word_table, file.path(OUT_DIR, "pathway_resistance_resilience_word_paste.tsv"))

rtf_escape <- function(x) {
  x <- replace_na(as.character(x), "—")
  x <- str_replace_all(x, fixed(RTF_BS), paste0(RTF_BS, RTF_BS))
  x <- str_replace_all(x, "([{}])", paste0(RTF_BS, "\\1"))
  str_replace_all(x, fixed("—"), paste0(RTF_BS, "u8212?"))
}
RTF_BS <- intToUtf8(92)
cell_rights <- c(1500, 6500, 8000, 9800, 11200)
rtf_row <- function(values, header = FALSE) {
  prefix <- paste0(RTF_BS, "trowd", RTF_BS, "trgaph60",
                   paste0(RTF_BS, "cellx", cell_rights, collapse = ""))
  body <- paste0(
    if (header) paste0(RTF_BS, "b ") else "",
    paste0(rtf_escape(values), RTF_BS, "cell", collapse = ""),
    if (header) paste0(RTF_BS, "b0") else "",
    RTF_BS, "row"
  )
  paste0(prefix, body)
}
rtf_lines <- c(
  paste0("{", RTF_BS, "rtf1", RTF_BS, "ansi", RTF_BS, "deff0{",
         RTF_BS, "fonttbl{", RTF_BS, "f0 Times New Roman;}}", RTF_BS, "fs18"),
  rtf_row(names(word_table), header = TRUE),
  apply(word_table, 1, rtf_row),
  paste0(RTF_BS, "pard", RTF_BS, "fs16 Unaffected means no supported fire-attributable departure: the response was nonsignificant, remained within the climatological reference band, or was already sustained before ignition. Pre-existing departures require at least one month outside the band.", RTF_BS, "par"),
  "}"
)
writeLines(rtf_lines, file.path(OUT_DIR, "pathway_resistance_resilience_word_table.rtf"), useBytes = TRUE)

if (requireNamespace("knitr", quietly = TRUE) && nzchar(Sys.which("quarto"))) {
  word_md <- file.path(OUT_DIR, "pathway_resistance_resilience_word_table.md")
  word_docx <- file.path(OUT_DIR, "pathway_resistance_resilience_word_table.docx")
  word_lines <- c(
    "---", "mainfont: Times New Roman", "fontsize: 9pt", "---", "",
    knitr::kable(word_table, format = "pipe", align = c("l", "l", "c", "c", "c")),
    "", "*Note:* Unaffected means no supported fire-attributable departure: the response was nonsignificant, remained within the climatological reference band, or was already sustained before ignition. Pre-existing departures require at least one month outside the band."
  )
  writeLines(word_lines, word_md, useBytes = TRUE)
  pandoc_status <- system2("quarto", c("pandoc", shQuote(word_md),
                                        "-o", shQuote(word_docx)))
  if (!identical(pandoc_status, 0L)) {
    warning("Could not create DOCX table export")
  } else if (requireNamespace("zip", quietly = TRUE)) {
    # Pandoc's DOCX defaults can override mainfont metadata. Enforce the same
    # Times New Roman face used by the figure and RTF in Word's XML styles.
    docx_dir <- tempfile("pathway_docx_")
    dir.create(docx_dir)
    utils::unzip(word_docx, exdir = docx_dir)
    for (xml_rel in c(file.path("word", "styles.xml"),
                      file.path("word", "document.xml"))) {
      xml_file <- file.path(docx_dir, xml_rel)
      if (!file.exists(xml_file)) next
      xml_text <- readr::read_file(xml_file)
      for (attr in c("ascii", "hAnsi", "eastAsia", "cs")) {
        xml_text <- str_replace_all(
          xml_text, paste0("w:", attr, "Theme=\"[^\"]+\""),
          paste0("w:", attr, "=\"Times New Roman\"")
        )
        xml_text <- str_replace_all(
          xml_text, paste0("w:", attr, "=\"[^\"]+\""),
          paste0("w:", attr, "=\"Times New Roman\"")
        )
      }
      readr::write_file(xml_text, xml_file)
    }
    patched_docx <- tempfile(fileext = ".docx")
    previous_wd <- getwd()
    setwd(docx_dir)
    zip::zip(patched_docx, list.files(".", recursive = TRUE,
                                     all.files = TRUE, no.. = TRUE),
             mode = "mirror")
    setwd(previous_wd)
    file.copy(patched_docx, word_docx, overwrite = TRUE)
    unlink(docx_dir, recursive = TRUE)
    unlink(patched_docx)
  }
}

# ---- One annotated 2021-2022 anomaly plot for every summarized variable -----
slugify <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

plot_one_variable <- function(pathway, variable) {
  dat <- series_table %>%
    filter(Pathway == pathway, Variable == variable,
           date >= as.Date("2021-01-01"), date <= FOLLOW_UP_END)
  mark <- marker_table %>% filter(Pathway == pathway, Variable == variable)
  metric <- summary_table %>% filter(Pathway == pathway, Variable == variable)
  if (nrow(dat) == 0 || nrow(mark) == 0 || nrow(metric) == 0) return(NA_character_)

  status_text <- case_when(
    mark$recovery_status[1] == "unaffected" ~
      "No departure from the historical climatological reference band",
    mark$recovery_status[1] == "recovered" ~ paste0(
      "Recovery: ", format(mark$recovery_date[1], "%b %Y"),
      " (", format(mark$recovery_months[1], trim = TRUE), " months after fire start)"),
    mark$recovery_status[1] == "not_recovered" ~ paste0(
      "Recovery not observed during ", format(mark$follow_up_months[1], trim = TRUE),
      " months of follow-up"),
    TRUE ~ "Recovery not estimable"
  )

  peak_dat <- dat %>% filter(date == mark$peak_date[1])
  recovery_dat <- dat %>% filter(date == mark$recovery_date[1])
  pathway_col <- unname(pathway_cols[pathway])

  p <- ggplot(dat, aes(date, anomaly_sd)) +
    geom_ribbon(aes(ymin = reference_lo_sd, ymax = reference_hi_sd),
                fill = "grey75", alpha = 0.25, colour = NA) +
    annotate("rect", xmin = FIRE_START, xmax = as.Date("2021-10-21"),
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.09) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
    geom_vline(xintercept = FIRE_START, colour = "firebrick", linetype = "dashed",
               linewidth = 0.45) +
    geom_line(colour = pathway_col, linewidth = 0.65, na.rm = TRUE) +
    geom_point(data = peak_dat, colour = "firebrick", fill = "white", shape = 21,
               size = 3.2, stroke = 0.9) +
    scale_x_date(limits = c(as.Date("2021-01-01"), FOLLOW_UP_END),
                 date_breaks = "1 year", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      title = variable,
      subtitle = paste0(pathway, " pathway | Peak = ", metric$`Peak anomaly`[1],
                        " | Peak timing = ", metric$`Timing of peak anomaly`[1]),
      x = NULL, y = "Seasonally standardized anomaly (SD)",
      caption = paste0(status_text,
                       ". Grey band = historical climatological reference interval; red shading = Caldor Fire window.")
    ) +
    theme_bw(base_size = 8, base_family = LO_FONT) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.25),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 7),
          plot.caption = element_text(size = 6.5, hjust = 0),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 6.5),
          plot.margin = margin(3, 4, 3, 4, "mm"))

  if (nrow(peak_dat) > 0) {
    p <- p + annotate(
      "label", x = peak_dat$date[1], y = peak_dat$anomaly_sd[1],
      label = paste0("Peak\n", sprintf("%+.2f SD", peak_dat$anomaly_sd[1])),
      hjust = if_else(peak_dat$date[1] < as.Date("2022-01-01"), -0.08, 1.08),
      vjust = if_else(peak_dat$anomaly_sd[1] >= 0, 1.15, -0.15),
      size = 2.2, family = LO_FONT, fill = "white", colour = "firebrick",
      linewidth = 0.2
    )
  }
  if (nrow(recovery_dat) > 0) {
    p <- p +
      geom_vline(xintercept = mark$recovery_date[1], colour = "#0072B2",
                 linetype = "longdash", linewidth = 0.45) +
      geom_point(data = recovery_dat, colour = "#0072B2", fill = "white",
                 shape = 23, size = 3.0, stroke = 0.9) +
      annotate("label", x = mark$recovery_date[1], y = recovery_dat$anomaly_sd[1],
               label = "Recovery", hjust = -0.08,
               vjust = if_else(recovery_dat$anomaly_sd[1] >= 0, 1.2, -0.2),
               size = 2.2, family = LO_FONT, fill = "white", colour = "#0072B2",
               linewidth = 0.2)
  }

  out_file <- file.path(PLOT_DIR, paste0(slugify(pathway), "__", slugify(variable), ".png"))
  ggsave(out_file, p, width = 12.7, height = 8.9, units = "cm", dpi = LO_DPI,
         device = ragg::agg_png)
  out_file
}

plot_manifest <- summary_table %>%
  mutate(plot_file = map2_chr(as.character(Pathway), Variable, plot_one_variable))
write_csv(plot_manifest, file.path(OUT_DIR, "pathway_resistance_resilience_plot_manifest.csv"))

# ---- Figure 2-style multipanel plot for each pathway ------------------------
# These pathway figures use raw units, a solid climatological mean, solid
# upper/lower reference limits, and a band bounded by those limits.
variable_unit <- function(variable) {
  case_when(
    str_detect(variable, "PM2.5") ~ "µg m⁻³",
    str_detect(variable, "PAR 1%|UV-320 1%|Secchi") ~ "m",
    str_detect(variable, "Shortwave") ~ "W m⁻²",
    str_detect(variable, "Surface temperature") ~ "°C",
    str_detect(variable, "Runoff") ~ "m³ s⁻¹",
    str_detect(variable, "Atmospheric deposition.*DIN:SRP|In-lake.*DIN:TRP") ~ "molar ratio",
    str_detect(variable, "Atmospheric deposition") ~ "mg m⁻² d⁻¹",
    str_detect(variable, "In-lake") ~ "µg L⁻¹",
    str_detect(variable, "Chlorophyll") ~ "µg L⁻¹",
    str_detect(variable, "abundance") ~ "cells L⁻¹",
    str_detect(variable, "biovolume") ~ "mm³ L⁻¹",
    str_detect(variable, "Community composition") ~ "Bray–Curtis distance",
    TRUE ~ "observed units"
  )
}

plot_pathway_composite <- function(pathway_name, variables = NULL,
                                   panel_labels = NULL, file_suffix = NULL) {
  if (is.null(variables)) {
    variables <- summary_table %>%
      filter(as.character(Pathway) == pathway_name) %>%
      pull(Variable)
  }
  if (length(variables) == 0L) return(NA_character_)
  if (is.null(panel_labels)) panel_labels <- letters[seq_along(variables)]
  stopifnot(length(panel_labels) == length(variables))
  panel_lookup <- tibble(
    Variable = variables,
    panel_order = seq_along(variables),
    unit = variable_unit(variables)
  )
  dat <- series_table %>%
    filter(Pathway == pathway_name, Variable %in% variables,
           date >= as.Date("2021-01-01"), date <= FOLLOW_UP_END) %>%
    left_join(panel_lookup, by = "Variable")
  marks <- marker_table %>%
    filter(Pathway == pathway_name, Variable %in% variables) %>%
    mutate(peak_date = as.Date(peak_date), recovery_date = as.Date(recovery_date))
  peak_points <- dat %>%
    inner_join(
      marks %>% filter(recovery_status != "unaffected") %>%
        select(Variable, peak_date),
      by = c("Variable", "date" = "peak_date")
    )
  recovery_points <- dat %>%
    inner_join(
      marks %>% filter(!is.na(recovery_date)) %>%
        select(Variable, recovery_date),
      by = c("Variable", "date" = "recovery_date")
    )
  pathway_col <- unname(pathway_cols[pathway_name])
  n_col <- if_else(length(variables) == 1L, 1L, 2L)
  n_row <- ceiling(length(variables) / n_col)

  make_panel <- function(variable, panel_index) {
    panel_dat <- dat %>% filter(Variable == variable)
    panel_peak <- peak_points %>% filter(Variable == variable)
    panel_recovery <- recovery_points %>% filter(Variable == variable)
    panel_unit <- panel_lookup %>% filter(Variable == variable) %>% pull(unit)
    ggplot(panel_dat, aes(date, value)) +
      geom_ribbon(
        aes(ymin = reference_lo, ymax = reference_hi,
            fill = "Historical reference interval"),
        alpha = 0.32, colour = NA, na.rm = TRUE
      ) +
      geom_line(aes(y = reference_lo, linetype = "Reference limits"), colour = "grey45",
                linewidth = 0.20, na.rm = TRUE) +
      geom_line(aes(y = expected, linetype = "Climatological mean"), colour = "grey30",
                linewidth = 0.26, na.rm = TRUE) +
      geom_line(aes(y = reference_hi, linetype = "Reference limits"), colour = "grey45",
                linewidth = 0.20, na.rm = TRUE) +
      annotate("rect", xmin = FIRE_START, xmax = FIRE_END,
               ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
      geom_line(aes(colour = "Observed"), linewidth = 0.38, na.rm = TRUE) +
      geom_point(data = panel_peak, aes(date, value, shape = "Peak"),
                 inherit.aes = FALSE, fill = "white", colour = "firebrick",
                 size = 2.3, stroke = 0.75) +
      geom_point(data = panel_recovery, aes(date, value, shape = "Recovery"),
                 inherit.aes = FALSE, fill = "white", colour = "#0072B2",
                 size = 2.3, stroke = 0.75) +
      geom_vline(data = tibble(fire_start = FIRE_START),
                 aes(xintercept = fire_start, linetype = "Fire start"),
                 inherit.aes = FALSE, colour = "firebrick", linewidth = 0.24) +
      scale_x_date(limits = c(as.Date("2021-01-01"), FOLLOW_UP_END),
                   date_breaks = "1 year", date_labels = "%Y",
                   expand = expansion(mult = c(0.01, 0.02))) +
      scale_colour_manual(values = c("Observed" = pathway_col), name = NULL) +
      scale_fill_manual(values = c("Historical reference interval" = "grey78"),
                        name = NULL, guide = "none") +
      scale_linetype_manual(
        values = c("Reference limits" = "solid", "Climatological mean" = "solid",
                   "Fire start" = "dashed"), name = NULL, guide = "none"
      ) +
      scale_shape_manual(values = c("Peak" = 21, "Recovery" = 23), name = NULL) +
      labs(title = paste0(panel_labels[panel_index], ")"),
           x = NULL,
           y = paste0(str_wrap(variable, width = 18), "\n(", panel_unit, ")")) +
      theme_bw(base_size = LO_FONT_MIN_PT, base_family = LO_FONT) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.22),
        axis.title.y = element_text(size = LO_FONT_MIN_PT),
        axis.text = element_text(size = LO_FONT_MIN_PT),
        plot.title = element_text(face = "bold", size = LO_FONT_MIN_PT,
                                  colour = pathway_col, lineheight = 0.95),
        plot.margin = margin(2, 2, 2, 4.5, "mm")
      )
  }
  panels <- map2(variables, seq_along(variables), make_panel)
  p <- wrap_plots(panels, ncol = n_col, guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = LO_FONT_MIN_PT),
      legend.key.height = unit(2.5, "mm"),
      legend.spacing.x = unit(1.5, "mm")
    )

  file_stub <- paste0(
    "resistance_resilience_", slugify(pathway_name),
    if_else(is.null(file_suffix), "", paste0("_", file_suffix))
  )
  out_file <- file.path(
    OUT_DIR, paste0(file_stub, ".png")
  )
  out_height <- if_else(length(variables) == 1L, 8.9, 15.24)
  ggsave(
    out_file, p, width = 12.7, height = out_height,
    units = "cm", dpi = LO_DPI, device = ragg::agg_png, bg = "white"
  )
  ggsave(
    file.path(OUT_DIR, paste0(file_stub, ".pdf")), p,
    width = 12.7, height = out_height, units = "cm", device = cairo_pdf,
    bg = "white"
  )
  out_file
}

pathway_order <- levels(summary_table$Pathway)
all_pathway_variables <- summary_table %>%
  mutate(Pathway = as.character(Pathway)) %>%
  group_by(Pathway) %>%
  summarise(variables = list(Variable), .groups = "drop") %>%
  arrange(match(Pathway, pathway_order))
pathway_composite_manifest <- pmap_dfr(
  list(all_pathway_variables$Pathway, all_pathway_variables$variables),
  function(pathway_name, variables) {
    chunks <- split(variables, ceiling(seq_along(variables) / 6))
    map_dfr(seq_along(chunks), function(chunk_index) {
      chunk_variables <- chunks[[chunk_index]]
      suffix <- paste0("page", chunk_index)
      png_file <- plot_pathway_composite(
        pathway_name, chunk_variables, letters[seq_along(chunk_variables)], suffix
      )
      tibble(
        Pathway = pathway_name,
        page = chunk_index,
        first_panel = "a",
        last_panel = letters[length(chunk_variables)],
        plot_file = png_file,
        pdf_file = sub("\\.png$", ".pdf", png_file)
      )
    })
  }
)
write_csv(
  pathway_composite_manifest,
  file.path(OUT_DIR, "resistance_resilience_pathway_plot_manifest.csv")
)

message("Saved summary CSV: ", csv_file)
message("Saved supplemental table figure: ", fig_file)
message("Saved ", sum(!is.na(plot_manifest$plot_file)), " individual variable plots: ", PLOT_DIR)
message("Saved six L&O-sized pathway pages: ", OUT_DIR)
