# Shared helpers for the Figure 5C compound-disturbance state-space analysis.
# These functions contain no Leptolyngbya-specific logic so the environmental
# state definitions can be generated and frozen before biological data are read.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

standard_doy <- function(date) {
  date <- as.Date(date)
  yday(date) - if_else(leap_year(date) & month(date) > 2L, 1L, 0L)
}

date_from_standard_doy <- function(year, doy) {
  as.Date(sprintf("%d-01-01", as.integer(year))) + as.integer(doy) - 1L +
    if_else(leap_year(as.integer(year)) & as.integer(doy) >= 60L, 1L, 0L)
}

quantile_window <- function(doy, lower_prob, upper_prob) {
  raw_start <- as.numeric(quantile(doy, lower_prob, type = 8, na.rm = TRUE))
  raw_end <- as.numeric(quantile(doy, upper_prob, type = 8, na.rm = TRUE))
  tibble(
    lower_probability = lower_prob,
    upper_probability = upper_prob,
    raw_start_doy = raw_start,
    raw_end_doy = raw_end,
    start_doy = floor(raw_start),
    end_doy = ceiling(raw_end)
  )
}

find_sustained_onset <- function(fitted_curve, fraction, sustain_days = 14L) {
  stopifnot(all(c("doy", "fitted_stability") %in% names(fitted_curve)))
  x <- fitted_curve %>% arrange(doy)
  minimum_i <- which.min(x$fitted_stability)
  circular_order <- c(minimum_i:nrow(x), seq_len(max(minimum_i - 1L, 0L)))
  maximum_position <- which.max(x$fitted_stability[circular_order])
  rising_i <- circular_order[seq_len(maximum_position)]
  minimum_value <- min(x$fitted_stability)
  maximum_value <- max(x$fitted_stability)
  threshold <- minimum_value + fraction * (maximum_value - minimum_value)

  qualifies <- vapply(rising_i, function(i) {
    next_i <- ((i - 1L + seq_len(sustain_days) - 1L) %% nrow(x)) + 1L
    all(x$fitted_stability[next_i] >= threshold)
  }, logical(1))
  if (!any(qualifies)) {
    stop("No sustained stratification onset found for threshold fraction ", fraction)
  }
  onset_i <- rising_i[which(qualifies)[1]]
  tibble(
    stability_fraction = fraction,
    onset_doy = x$doy[onset_i],
    onset_month_day = format(date_from_standard_doy(2001L, x$doy[onset_i]), "%m-%d"),
    sustain_days = sustain_days,
    baseline_minimum_j_m2 = minimum_value,
    annual_maximum_j_m2 = maximum_value,
    seasonal_amplitude_j_m2 = maximum_value - minimum_value,
    threshold_j_m2 = threshold
  )
}

trapz <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

calculate_mltp_cast_metrics <- function(profile, mld_function) {
  temp <- profile %>% filter(is.finite(depth_m), is.finite(Temperature)) %>% arrange(depth_m)
  density <- profile %>% filter(is.finite(depth_m), is.finite(Density)) %>% arrange(depth_m)
  max_depth <- suppressWarnings(max(profile$depth_m, na.rm = TRUE))
  mld <- mld_function(density$depth_m, density$Density, 0.1, "increase")
  if (!is.finite(max_depth)) max_depth <- NA_real_

  epilimnetic_temperature <- if (is.finite(mld)) {
    mean(temp$Temperature[temp$depth_m <= mld], na.rm = TRUE)
  } else NA_real_
  below_mixed_layer_temperature <- if (is.finite(mld)) {
    mean(temp$Temperature[temp$depth_m > mld & temp$depth_m <= mld + 10], na.rm = TRUE)
  } else NA_real_
  surface_temperature <- mean(temp$Temperature[temp$depth_m <= 10], na.rm = TRUE)

  tibble(
    max_profile_depth_m = max_depth,
    mixed_layer_depth_m = mld,
    epilimnetic_temperature_c = if_else(is.nan(epilimnetic_temperature), NA_real_, epilimnetic_temperature),
    below_mixed_layer_temperature_c = if_else(is.nan(below_mixed_layer_temperature), NA_real_, below_mixed_layer_temperature),
    delta_t_mixed_layer_c = epilimnetic_temperature_c - below_mixed_layer_temperature_c,
    surface_temperature_0_10m_c = if_else(is.nan(surface_temperature), NA_real_, surface_temperature)
  )
}

summarize_antecedent_state <- function(stability, cast_metrics, onset_doy, smoke_start_doy,
                                       antecedent_end_doy = NULL) {
  end_doy <- antecedent_end_doy %||% (smoke_start_doy - 1L)
  if (onset_doy > end_doy) stop("Antecedent window starts after it ends.")

  stability_annual <- stability %>%
    filter(doy >= onset_doy, doy <= end_doy, is.finite(schmidt_stability_j_m2)) %>%
    group_by(year) %>%
    summarise(
      schmidt_stability_j_m2 = median(schmidt_stability_j_m2),
      n_stability_casts = n(),
      stability_first_date = min(date),
      stability_last_date = max(date),
      stability_median_date = as.Date(median(as.numeric(date)), origin = "1970-01-01"),
      .groups = "drop"
    )

  profile_annual <- cast_metrics %>%
    filter(doy >= onset_doy, doy <= end_doy) %>%
    group_by(year) %>%
    summarise(
      mixed_layer_depth_m = median(mixed_layer_depth_m, na.rm = TRUE),
      delta_t_mixed_layer_c = median(delta_t_mixed_layer_c, na.rm = TRUE),
      surface_temperature_0_10m_c = median(surface_temperature_0_10m_c, na.rm = TRUE),
      epilimnetic_temperature_c = median(epilimnetic_temperature_c, na.rm = TRUE),
      n_profile_casts = n(),
      n_mixed_layer_depth = sum(is.finite(mixed_layer_depth_m)),
      n_delta_t = sum(is.finite(delta_t_mixed_layer_c)),
      n_surface_temperature = sum(is.finite(surface_temperature_0_10m_c)),
      profile_first_date = min(date),
      profile_last_date = max(date),
      profile_median_date = as.Date(median(as.numeric(date)), origin = "1970-01-01"),
      .groups = "drop"
    ) %>%
    mutate(across(c(mixed_layer_depth_m, delta_t_mixed_layer_c,
                    surface_temperature_0_10m_c, epilimnetic_temperature_c),
                  ~if_else(is.nan(.x), NA_real_, .x)))

  full_join(stability_annual, profile_annual, by = "year") %>%
    mutate(antecedent_start_doy = onset_doy, antecedent_end_doy = end_doy)
}

summarize_deposition_window <- function(deposition, years, start_doy, end_doy,
                                        minimum_interval_coverage = 0.8) {
  map_dfr(years, function(target_year) {
    window_start <- date_from_standard_doy(target_year, start_doy)
    window_end <- date_from_standard_doy(target_year, end_doy)
    window_end_exclusive <- window_end + 1L
    window_days <- as.numeric(window_end_exclusive - window_start)

    overlapping <- deposition %>%
      mutate(
        overlap_start = as.Date(
          pmax(as.numeric(start_date), as.numeric(window_start)), origin = "1970-01-01"
        ),
        overlap_end_exclusive = as.Date(
          pmin(as.numeric(end_date), as.numeric(window_end_exclusive)), origin = "1970-01-01"
        ),
        overlap_days = pmax(0, as.numeric(overlap_end_exclusive - overlap_start)),
        combined_tn_tp_daily_load_mg_m2_d = tn_daily_load_mg_m2_d + tp_daily_load_mg_m2_d
      ) %>%
      filter(overlap_days > 0)

    covered_days <- sum(
      overlapping$overlap_days[is.finite(overlapping$combined_tn_tp_daily_load_mg_m2_d)],
      na.rm = TRUE
    )
    integrated_load <- sum(
      overlapping$overlap_days * overlapping$combined_tn_tp_daily_load_mg_m2_d,
      na.rm = TRUE
    )
    coverage <- min(covered_days / window_days, 1)
    tibble(
      year = target_year,
      deposition_window_days = window_days,
      deposition_days_covered = covered_days,
      deposition_interval_coverage = coverage,
      n_deposition_intervals = sum(
        is.finite(overlapping$combined_tn_tp_daily_load_mg_m2_d)
      ),
      integrated_total_nutrient_deposition_mg_m2 = if_else(
        coverage >= minimum_interval_coverage, integrated_load, NA_real_
      )
    )
  })
}

summarize_wildfire_exposure <- function(smoke_radiation_daily, start_doy, end_doy,
                                        optical_depth = NULL,
                                        deposition = NULL,
                                        pm_threshold = 9.1,
                                        high_pm_threshold = 35.5,
                                        minimum_daily_coverage = 0.8) {
  smoke_summary <- smoke_radiation_daily %>%
    filter(doy >= start_doy, doy <= end_doy) %>%
    group_by(year) %>%
    summarise(
      wildfire_window_days = n(),
      pm25_days_available = sum(is.finite(pm25_max)),
      shortwave_days_available = sum(is.finite(sw_diff)),
      pm25_daily_coverage = pm25_days_available / wildfire_window_days,
      shortwave_daily_coverage = shortwave_days_available / wildfire_window_days,
      smoke_day_count_raw = sum(smoke_confirmed %in% TRUE, na.rm = TRUE),
      high_pm25_smoke_day_count_raw = sum(
        smoke_confirmed %in% TRUE & is.finite(pm25_max) & pm25_max >= high_pm_threshold,
        na.rm = TRUE
      ),
      integrated_pm25_excess_ug_d_m3_raw = sum(
        if_else(smoke_confirmed %in% TRUE & is.finite(pm25_max),
                pmax(pm25_max - pm_threshold, 0), 0),
        na.rm = TRUE
      ),
      smoke_days_with_shortwave = sum(smoke_confirmed %in% TRUE & is.finite(sw_diff)),
      integrated_smoke_shortwave_deficit_w_m2_d_raw = sum(
        if_else(smoke_confirmed %in% TRUE & is.finite(sw_diff), pmax(sw_diff, 0), 0),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      smoke_shortwave_coverage = if_else(
        smoke_day_count_raw > 0,
        smoke_days_with_shortwave / smoke_day_count_raw,
        shortwave_daily_coverage
      ),
      smoke_day_count = if_else(pm25_daily_coverage >= minimum_daily_coverage,
                                as.numeric(smoke_day_count_raw), NA_real_),
      high_pm25_smoke_day_count = if_else(pm25_daily_coverage >= minimum_daily_coverage,
                                          as.numeric(high_pm25_smoke_day_count_raw), NA_real_),
      integrated_pm25_excess_ug_d_m3 = if_else(
        pm25_daily_coverage >= minimum_daily_coverage,
        integrated_pm25_excess_ug_d_m3_raw, NA_real_
      ),
      integrated_smoke_shortwave_deficit_w_m2_d = if_else(
        pm25_daily_coverage >= minimum_daily_coverage &
          smoke_shortwave_coverage >= minimum_daily_coverage,
        integrated_smoke_shortwave_deficit_w_m2_d_raw, NA_real_
      ),
      wildfire_start_doy = start_doy,
      wildfire_end_doy = end_doy
    )

  if (!is.null(optical_depth)) {
    optical_summary <- optical_depth %>%
      filter(doy >= start_doy, doy <= end_doy) %>%
      group_by(year) %>%
      summarise(
        n_par_depth_observations = sum(is.finite(par_1pct_depth_m)),
        n_uv320_depth_observations = sum(is.finite(uv320_1pct_depth_m)),
        par_1pct_depth_m = median(par_1pct_depth_m, na.rm = TRUE),
        uv320_1pct_depth_m = median(uv320_1pct_depth_m, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        par_1pct_depth_m = if_else(is.nan(par_1pct_depth_m), NA_real_, par_1pct_depth_m),
        uv320_1pct_depth_m = if_else(is.nan(uv320_1pct_depth_m), NA_real_, uv320_1pct_depth_m)
      )
    smoke_summary <- smoke_summary %>% full_join(optical_summary, by = "year")
  }

  if (!is.null(deposition)) {
    deposition_summary <- summarize_deposition_window(
      deposition, smoke_summary$year, start_doy, end_doy,
      minimum_interval_coverage = minimum_daily_coverage
    )
    smoke_summary <- smoke_summary %>% left_join(deposition_summary, by = "year")
  }
  smoke_summary
}

variable_screen <- function(data, metadata, focal_years = c(2011L, 2020L, 2021L),
                            minimum_coverage = 0.8, redundancy_cutoff = 0.9,
                            required_variables = character(),
                            minimum_variables = 2L) {
  variables <- metadata$variable
  availability <- data %>%
    select(year, all_of(variables)) %>%
    mutate(across(all_of(variables), is.finite))

  coverage <- tibble(
    variable = variables,
    coverage_fraction = map_dbl(variables, ~mean(is.finite(data[[.x]]))),
    all_focal_years_available = map_lgl(
      variables,
      ~all(is.finite(data[[.x]][match(focal_years, data$year)]))
    )
  ) %>%
    left_join(metadata, by = "variable") %>%
    arrange(priority)

  correlation <- data %>%
    select(all_of(variables)) %>%
    cor(use = "pairwise.complete.obs")

  decisions <- coverage %>%
    mutate(
      adequate_coverage = coverage_fraction >= minimum_coverage & all_focal_years_available,
      retained = FALSE,
      decision_reason = if_else(
        adequate_coverage, "eligible after coverage screen",
        "removed: inadequate coverage or missing focal year"
      )
    )

  retained <- character()
  eligible_variables <- decisions$variable[
    decisions$adequate_coverage | decisions$variable %in% required_variables
  ]
  for (v in eligible_variables) {
    if (v %in% required_variables) {
      retained <- c(retained, v)
      decisions$retained[decisions$variable == v] <- TRUE
      decisions$decision_reason[decisions$variable == v] <- if (
        decisions$adequate_coverage[decisions$variable == v]
      ) {
        "retained: required by frozen specification"
      } else {
        "retained with explicit missing-data estimation: required by frozen specification"
      }
      next
    }
    if (length(retained) == 0L) {
      retained <- c(retained, v)
      decisions$retained[decisions$variable == v] <- TRUE
      decisions$decision_reason[decisions$variable == v] <- "retained"
      next
    }
    r <- abs(correlation[v, retained, drop = TRUE])
    if (all(!is.finite(r) | r < redundancy_cutoff)) {
      retained <- c(retained, v)
      decisions$retained[decisions$variable == v] <- TRUE
      decisions$decision_reason[decisions$variable == v] <- "retained"
    } else {
      redundant_with <- retained[which.max(replace(r, !is.finite(r), -Inf))]
      decisions$decision_reason[decisions$variable == v] <- paste0(
        "removed: |r| >= ", redundancy_cutoff, " with ", redundant_with
      )
    }
  }

  if (length(retained) < minimum_variables) {
    stop("Fewer than ", minimum_variables,
         " variables survived the environmental screen: ",
         paste(retained, collapse = ", "))
  }
  list(
    retained = retained,
    availability = availability,
    coverage = decisions,
    correlation = correlation
  )
}

fit_oriented_pca <- function(data, variables, directions, index_name, block_name,
                             impute_missing = FALSE, imputation_ncp = 1L) {
  raw_values <- data %>% select(all_of(variables)) %>%
    mutate(across(everything(), ~if_else(is.finite(.x), .x, NA_real_)))
  if (impute_missing) {
    if (!requireNamespace("missMDA", quietly = TRUE)) {
      stop("Package 'missMDA' is required for explicit PCA missing-data estimation.")
    }
    if (any(map_lgl(raw_values, ~all(is.na(.x))))) {
      stop("At least one ", block_name, " variable is missing for every analysis year.")
    }
    maximum_ncp <- min(nrow(raw_values) - 2L, ncol(raw_values) - 1L)
    ncp_used <- min(as.integer(imputation_ncp), maximum_ncp)
    if (anyNA(raw_values)) {
      set.seed(2021)
      completed_values <- missMDA::imputePCA(
        as.data.frame(raw_values), ncp = ncp_used, scale = TRUE,
        method = "Regularized"
      )$completeObs %>%
        as_tibble()
    } else {
      completed_values <- raw_values
    }
    model_data <- bind_cols(data %>% select(year), completed_values)
  } else {
    model_data <- data %>% filter(if_all(all_of(variables), is.finite))
    completed_values <- model_data %>% select(all_of(variables))
    raw_values <- completed_values
    ncp_used <- NA_integer_
  }
  if (nrow(model_data) < max(5L, length(variables) + 1L)) {
    stop("Insufficient complete years for ", block_name, " PCA.")
  }
  fit <- prcomp(completed_values, center = TRUE, scale. = TRUE)
  orientation <- sum(fit$rotation[, 1] * directions[variables])
  if (orientation < 0) {
    fit$rotation[, 1] <- -fit$rotation[, 1]
    fit$x[, 1] <- -fit$x[, 1]
  }
  if (ncol(fit$x) >= 2L) {
    anchor <- which.max(abs(fit$rotation[, 2]))
    if (fit$rotation[anchor, 2] < 0) {
      fit$rotation[, 2] <- -fit$rotation[, 2]
      fit$x[, 2] <- -fit$x[, 2]
    }
  }
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  scores <- tibble(
    year = model_data$year,
    !!index_name := fit$x[, 1],
    pc2 = if (ncol(fit$x) >= 2L) fit$x[, 2] else NA_real_,
    n_variables_observed = rowSums(is.finite(as.matrix(raw_values))),
    n_values_imputed = length(variables) - n_variables_observed,
    imputed_variables = pmap_chr(raw_values, function(...) {
      missing <- variables[!is.finite(c(...))]
      if (length(missing) == 0L) "none" else paste(missing, collapse = ";")
    })
  )
  loadings <- as_tibble(fit$rotation, rownames = "variable") %>%
    mutate(block = block_name, .before = 1)
  variance_table <- tibble(
    block = block_name,
    component = paste0("PC", seq_along(variance)),
    variance_explained = variance
  )
  completed <- bind_cols(model_data %>% select(year), completed_values)
  list(
    fit = fit, scores = scores, loadings = loadings,
    variance = variance_table, years = model_data$year,
    completed = completed, imputation_ncp = ncp_used
  )
}

build_environment_scenario <- function(stability, cast_metrics, smoke_radiation_daily,
                                       onset_doy, smoke_start_doy, smoke_end_doy,
                                       state_metadata, wildfire_metadata,
                                       antecedent_end_doy = NULL,
                                       optical_depth = NULL,
                                       deposition = NULL,
                                       analysis_years = NULL,
                                       focal_years = c(2011L, 2020L, 2021L),
                                       fixed_state_variables = NULL,
                                       fixed_wildfire_variables = NULL,
                                       required_state_variables = character(),
                                       required_wildfire_variables = character(),
                                       wildfire_index_variable = NULL,
                                       impute_missing_pca = FALSE,
                                       imputation_ncp = 1L,
                                       minimum_coverage = 0.8,
                                       redundancy_cutoff = 0.9) {
  state <- summarize_antecedent_state(
    stability, cast_metrics, onset_doy, smoke_start_doy, antecedent_end_doy
  )
  wildfire <- summarize_wildfire_exposure(
    smoke_radiation_daily, smoke_start_doy, smoke_end_doy,
    optical_depth = optical_depth,
    deposition = deposition,
    minimum_daily_coverage = minimum_coverage
  )
  if (is.null(analysis_years)) {
    year_min <- max(min(stability$year), min(smoke_radiation_daily$year))
    year_max <- min(max(stability$year), max(smoke_radiation_daily$year))
    analysis_years <- year_min:year_max
  }
  annual <- tibble(year = as.integer(analysis_years)) %>%
    left_join(state, by = "year") %>%
    left_join(wildfire, by = "year")

  state_screen <- variable_screen(
    annual, state_metadata, focal_years = focal_years,
    minimum_coverage = minimum_coverage,
    redundancy_cutoff = redundancy_cutoff,
    required_variables = required_state_variables
  )
  wildfire_screen <- variable_screen(
    annual, wildfire_metadata, focal_years = focal_years,
    minimum_coverage = minimum_coverage,
    redundancy_cutoff = redundancy_cutoff,
    required_variables = required_wildfire_variables,
    minimum_variables = if (is.null(wildfire_index_variable)) 2L else 1L
  )
  state_variables <- fixed_state_variables %||% state_screen$retained
  wildfire_variables <- fixed_wildfire_variables %||% wildfire_screen$retained

  complete_case <- annual %>%
    filter(if_all(all_of(c(state_variables, wildfire_variables)), is.finite))
  model_data <- if (impute_missing_pca) annual else complete_case
  state_pca <- fit_oriented_pca(
    model_data, state_variables,
    setNames(state_metadata$direction, state_metadata$variable),
    "antecedent_lake_state_index", "Antecedent lake state",
    impute_missing = impute_missing_pca, imputation_ncp = imputation_ncp
  )
  wildfire_pca <- if (is.null(wildfire_index_variable)) {
    fit_oriented_pca(
      model_data, wildfire_variables,
      setNames(wildfire_metadata$direction, wildfire_metadata$variable),
      "wildfire_exposure_index", "Wildfire exposure",
      impute_missing = impute_missing_pca, imputation_ncp = imputation_ncp
    )
  } else NULL
  scores <- model_data %>%
    select(year) %>%
    left_join(
      state_pca$scores %>%
        select(-pc2) %>%
        rename_with(~paste0("antecedent_", .x),
                    c(n_variables_observed, n_values_imputed, imputed_variables)),
      by = "year"
    ) %>%
    {
      if (is.null(wildfire_index_variable)) {
        left_join(
          .,
          wildfire_pca$scores %>%
            select(-pc2) %>%
            rename_with(~paste0("wildfire_", .x),
                        c(n_variables_observed, n_values_imputed, imputed_variables)),
          by = "year"
        )
      } else {
        left_join(
          .,
          annual %>%
            transmute(
              year,
              wildfire_exposure_index = .data[[wildfire_index_variable]],
              wildfire_n_variables_observed = as.integer(
                is.finite(.data[[wildfire_index_variable]])
              ),
              wildfire_n_values_imputed = 0L,
              wildfire_imputed_variables = "none"
            ),
          by = "year"
        )
      }
    }
  list(
    annual = annual,
    common = model_data,
    complete_case = complete_case,
    scores = scores,
    state_screen = state_screen,
    wildfire_screen = wildfire_screen,
    state_variables = state_variables,
    wildfire_variables = wildfire_variables,
    state_pca = state_pca,
    wildfire_pca = wildfire_pca
  )
}

matrix_to_long <- function(x, block) {
  as.data.frame(x) %>%
    rownames_to_column("variable_1") %>%
    pivot_longer(-variable_1, names_to = "variable_2", values_to = "correlation") %>%
    mutate(block = block, .before = 1)
}
