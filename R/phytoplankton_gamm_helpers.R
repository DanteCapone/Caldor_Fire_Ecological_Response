# =============================================================================
# phytoplankton_gamm_helpers.R
# Caldor Fire Ecosystem Response Project
#
# Purpose: Shared helper functions for the leave-one-year-out GAMM anomaly
#          analysis of Leptolyngbya spp. at the LTP station (Analysis A).
#
# Data source: data/lake_environmental_data/phytoplankton/
#              caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv
#              (built by scripts/01b_full_phytoplankton_timeseries.R)
#
# Response-variable decision: Leptolyngbya spp. biovolume (mm3 m-3) is used
# in preference to abundance because biovolume is available for every
# enumerated sample and integrates cell-size variation that abundance alone
# does not capture. Raw units are preserved; log1p(biovolume) is used only
# as a diagnostic transform, never as the modeled response.
#
# Zero handling: the source file is a long-format enumeration table (one row
# per taxon identified in a sample). A sample date-depth combination with no
# Leptolyngbya spp. row is a true analytical zero, not a missing value,
# provided the sample was demonstrably enumerated for Cyanobacteria (i.e. at
# least one Cyanobacteria taxon was identified somewhere in the record for
# that station, and the sample itself yielded a minimum number of identified
# taxa). Samples below the minimum taxon-count threshold are flagged as
# low-effort counts and excluded from the primary model (see
# MIN_TAXA_PER_SAMPLE below).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(mgcv)
})

LEPTO_TAXON <- "Leptolyngbya spp."
MIN_DEPTH_TIME_POINTS <- 5L   # minimum distinct sample dates required to keep a depth
MIN_TAXA_PER_SAMPLE <- 5L     # minimum taxa identified per sample to trust a zero
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END <- as.Date("2021-10-21")
GAMM_SEED <- 20260804L

phyto_source_file <- function(proj_root) {
  file.path(
    proj_root, "data", "lake_environmental_data", "phytoplankton",
    "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"
  )
}

#' Build the zero-filled, analysis-ready Leptolyngbya spp. sample table.
#'
#' Every LTP date-depth sample that meets the minimum-taxon-count and
#' minimum-depth-coverage rules contributes exactly one row. Samples in
# which Leptolyngbya spp. was not identified receive abundance = 0 and
#' biovolume = 0 (true zero), not NA.
build_leptolyngbya_samples <- function(proj_root) {
  src <- phyto_source_file(proj_root)
  if (!file.exists(src)) {
    stop("Required phytoplankton file not found: ", src,
         "\nRun scripts/01b_full_phytoplankton_timeseries.R first.")
  }

  phyto <- read_csv(src, show_col_types = FALSE) %>%
    filter(station_id == "LTP", !is.na(date), !is.na(depth_num))

  required_cols <- c(
    "station_id", "date", "depth_num", "taxon", "abundance", "biovolume",
    "year", "month", "season", "fire_period", "sample_id"
  )
  missing_cols <- setdiff(required_cols, names(phyto))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in phytoplankton source: ",
         paste(missing_cols, collapse = ", "))
  }

  # Sample-level enumeration effort (used to distinguish true zeros from
  # samples too sparsely counted to trust an absence).
  sample_effort <- phyto %>%
    group_by(sample_id, date, depth_num) %>%
    summarise(n_taxa_identified = n_distinct(taxon), .groups = "drop") %>%
    mutate(low_effort_sample = n_taxa_identified < MIN_TAXA_PER_SAMPLE)

  cyanobacteria_ever_identified <- phyto %>%
    filter(phylum == "Cyanobacteria") %>%
    nrow() > 0
  if (!cyanobacteria_ever_identified) {
    stop("No Cyanobacteria taxa found in the LTP record; cannot treat ",
         "missing Leptolyngbya rows as true zeros.")
  }

  depth_coverage <- phyto %>%
    distinct(depth_num, date) %>%
    count(depth_num, name = "n_time_points") %>%
    mutate(depth_included = n_time_points >= MIN_DEPTH_TIME_POINTS)
  eligible_depths <- depth_coverage %>% filter(depth_included) %>% pull(depth_num)

  lepto_present <- phyto %>%
    filter(taxon == LEPTO_TAXON) %>%
    transmute(
      sample_id, date, depth_num,
      abundance_obs = abundance,
      biovolume_obs = biovolume
    )

  samples <- sample_effort %>%
    filter(depth_num %in% eligible_depths) %>%
    left_join(lepto_present, by = c("sample_id", "date", "depth_num")) %>%
    mutate(
      taxon_evaluated = TRUE,
      abundance = coalesce(abundance_obs, 0),
      biovolume = coalesce(biovolume_obs, 0),
      analytical_zero = is.na(abundance_obs)
    ) %>%
    select(-abundance_obs, -biovolume_obs)

  station_meta <- phyto %>%
    distinct(sample_id, date, depth_num, station_id, year, month, season, fire_period)

  samples <- samples %>%
    left_join(station_meta, by = c("sample_id", "date", "depth_num")) %>%
    mutate(
      doy = yday(date),
      # Cyclic day-of-year on a 366-day circle so that Dec 31 / Jan 1 are adjacent.
      doy_cyclic = doy,
      year_factor = factor(year),
      depth_factor = factor(depth_num),
      log1p_biovolume = log1p(biovolume),
      log1p_abundance = log1p(abundance)
    ) %>%
    arrange(date, depth_num)

  attr(samples, "eligible_depths") <- eligible_depths
  attr(samples, "excluded_depths") <- depth_coverage %>%
    filter(!depth_included) %>% pull(depth_num)
  samples
}

#' Summarize sampling effort by year / season / depth (diagnostic table).
summarize_sampling_effort <- function(samples) {
  samples %>%
    group_by(year, season, depth_num) %>%
    summarise(
      n_samples = n(),
      n_low_effort = sum(low_effort_sample),
      n_dates = n_distinct(date),
      .groups = "drop"
    ) %>%
    arrange(year, season, depth_num)
}

#' Identify years with insufficient temporal coverage for the LOYO analysis.
flag_insufficient_years <- function(samples, min_dates_per_year = 4L) {
  samples %>%
    group_by(year) %>%
    summarise(n_dates = n_distinct(date), .groups = "drop") %>%
    mutate(sufficient_coverage = n_dates >= min_dates_per_year)
}

#' Historical (non-2021) bloom threshold: upper seasonal percentile of
#' log1p(biovolume) computed while excluding the focal year, evaluated at
#' shallow depths (<=20 m) where the bloom concentrates.
compute_bloom_threshold <- function(samples, exclude_year = 2021L,
                                     percentile = 0.90, max_depth_m = 20) {
  ref <- samples %>%
    filter(year != exclude_year, depth_num <= max_depth_m, !low_effort_sample)
  quantile(ref$biovolume, probs = percentile, na.rm = TRUE, type = 7)
}

#' Fit the primary Tweedie bam() GAMM to a training set (all years but the
#' held-out focal year when used inside the LOYO loop).
fit_lepto_gamm <- function(train_data, k_doy = 10, k_year = 8, k_depth = 5) {
  bam(
    biovolume ~
      s(doy_cyclic, bs = "cc", k = k_doy) +
      s(year, k = k_year) +
      s(depth_num, k = k_depth) +
      s(year_factor, bs = "re"),
    data = train_data,
    family = tw(link = "log"),
    knots = list(doy_cyclic = c(0.5, 366.5)),
    method = "fREML",
    discrete = TRUE
  )
}

#' Predict the expected seasonal trajectory (response scale) for a focal
#' year on a regular grid of days, holding depth at the observed sampling
#' depths for that year (so predictions match the true sampling structure).
#'
#' Uncertainty: `ci_mean` from predict.gam (parameter uncertainty only,
#' excludes the year random effect which cannot be predicted for an unseen
#' year). `pi_obs` is a simulation-based predictive interval for a new
#' observation that adds Tweedie observation-level variance and one
#' simulated draw from the year random-effect distribution, using
#' `simulate_new_year_re()`.
predict_focal_year <- function(model, focal_data, n_sim = 1000,
                                seed = GAMM_SEED) {
  set.seed(seed)
  grid <- focal_data %>% distinct(doy_cyclic, year, depth_num)

  # The focal year's random-effect level was never observed by this model
  # (it was excluded from training). A dummy in-sample level is supplied so
  # model.frame() can be built, and s(year_factor) is excluded from the
  # linear predictor so its contribution is exactly zero for the focal year.
  train_year_levels <- levels(model$model$year_factor)
  grid$year_factor <- factor(train_year_levels[1], levels = train_year_levels)

  lp <- predict(
    model, newdata = grid, type = "link", se.fit = TRUE,
    exclude = "s(year_factor)"
  )
  grid$fit_link <- lp$fit
  grid$se_link <- lp$se.fit
  grid$expected_mean <- exp(grid$fit_link)
  grid$ci_lower <- exp(grid$fit_link - 1.96 * grid$se_link)
  grid$ci_upper <- exp(grid$fit_link + 1.96 * grid$se_link)

  re_sd <- gam.vcomp(model)
  year_re_sd <- suppressWarnings(re_sd["s(year_factor)", "std.dev"])
  if (is.na(year_re_sd) || is.null(year_re_sd)) year_re_sd <- 0

  p_param <- rmvn(n_sim, coef(model), model$Vp)
  Xp <- predict(
    model, newdata = grid, type = "lpmatrix", exclude = "s(year_factor)"
  )
  sim_link <- Xp %*% t(p_param)
  sim_link <- sim_link + matrix(
    rnorm(n_sim * nrow(grid), mean = 0, sd = year_re_sd),
    nrow = nrow(grid), ncol = n_sim
  )
  sim_mu <- exp(sim_link)

  phi <- model$sig2
  p_tw <- model$family$getTheta(TRUE)
  sim_obs <- matrix(
    tweedie::rtweedie(
      length(sim_mu), mu = as.vector(sim_mu),
      phi = phi, power = p_tw
    ),
    nrow = nrow(grid), ncol = n_sim
  )

  grid$pi_lower <- apply(sim_obs, 1, quantile, probs = 0.025, na.rm = TRUE)
  grid$pi_upper <- apply(sim_obs, 1, quantile, probs = 0.975, na.rm = TRUE)

  grid %>% left_join(
    focal_data %>% select(date, doy_cyclic, year, depth_num, biovolume, abundance),
    by = c("doy_cyclic", "year", "depth_num")
  )
}

#' Annual anomaly metrics comparing observed values with the LOYO expected
#' trajectory. All metrics are computed identically for every year.
compute_annual_anomaly_metrics <- function(pred_obs, bloom_threshold) {
  pred_obs %>%
    group_by(year) %>%
    summarise(
      n_obs = n(),
      max_std_anomaly = max((biovolume - expected_mean) / pmax(se_link, 1e-6), na.rm = TRUE),
      integrated_positive_anomaly = sum(pmax(biovolume - expected_mean, 0), na.rm = TRUE),
      integrated_anomaly_above_pi = sum(pmax(biovolume - pi_upper, 0), na.rm = TRUE),
      n_above_pi = sum(biovolume > pi_upper, na.rm = TRUE),
      prop_above_pi = mean(biovolume > pi_upper, na.rm = TRUE),
      annual_peak_biovolume = max(biovolume, na.rm = TRUE),
      peak_date = date[which.max(biovolume)],
      cumulative_biovolume = sum(biovolume, na.rm = TRUE),
      bloom_days = sum(biovolume >= bloom_threshold, na.rm = TRUE),
      bloom_initiation = if (any(biovolume >= bloom_threshold)) {
        min(date[biovolume >= bloom_threshold])
      } else {
        as.Date(NA)
      },
      bloom_termination = if (any(biovolume >= bloom_threshold)) {
        max(date[biovolume >= bloom_threshold])
      } else {
        as.Date(NA)
      },
      .groups = "drop"
    ) %>%
    mutate(
      bloom_duration_days = if_else(
        bloom_days > 0, as.numeric(bloom_termination - bloom_initiation), NA_real_
      )
    )
}

#' Rank years for a metric and report the empirical percentile of 2021.
rank_2021 <- function(metrics, metric_col, focal_year = 2021L) {
  m <- metrics %>%
    mutate(rank = rank(-.data[[metric_col]], ties.method = "min")) %>%
    arrange(rank)
  focal_rank <- m$rank[m$year == focal_year]
  percentile <- 100 * (1 - (focal_rank - 1) / nrow(m))
  list(table = m, focal_rank = focal_rank, focal_percentile = percentile)
}
