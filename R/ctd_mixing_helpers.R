# =============================================================================
# ctd_mixing_helpers.R
# Caldor Fire Ecosystem Response Project
#
# Purpose: Shared helper functions for the mixed-layer dilution and
# depth-integrated biomass analysis (Analysis B). All physical metrics are
# calculated at the Index (LTP) CTD station, the same station where the
# Leptolyngbya spp. phytoplankton samples were collected, so that mixed-layer
# depth and stratification metrics describe the water column actually
# sampled by the phytoplankton net/bottle casts.
#
# Density is taken directly from the CTD Density channel (kg m-3) rather
# than recomputed from temperature alone, because Lake Tahoe has measurable
# (if small) specific conductance that the instrument's own density
# calculation already accounts for. No oceanic salinity correction is
# applied; Practical_Salinity in this record is a freshwater conductivity-
# based value, consistent with USGS/UC Davis TERC processing of this
# station.
#
# Two mixed-layer-depth (MLD) definitions are calculated for every cast:
#   - Density threshold: shallowest depth where density exceeds the
#     near-surface reference density by DELTA_RHO_MLD (kg m-3). Primary
#     definition, used because Density is a directly measured CTD channel.
#   - Temperature threshold: shallowest depth where temperature drops below
#     the near-surface reference temperature by DELTA_T_MLD (deg C).
#     Reported as a sensitivity check.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
  library(rLakeAnalyzer)
})

CTD_STATION <- "Index"           # Station_ID value for LTP in the CTD record
CTD_STATION_ALIAS <- "LTP"
DELTA_RHO_MLD <- 0.1              # kg m-3, primary MLD definition
DELTA_T_MLD <- 1.0                # deg C, sensitivity MLD definition
MIN_PROFILE_DEPTH_M <- 30         # casts shallower than this cannot resolve
                                   # epilimnion/metalimnion structure reliably
MIN_STRATIFICATION_RANGE_KGM3 <- 0.15  # below this, the near-surface profile
                                         # is treated as effectively unstratified
SURFACE_LAYER_SEARCH_M <- 50       # thermocline/N2 search window. Lake Tahoe
                                     # is deep enough (>400 m) that a whole-
                                     # profile gradient search can lock onto
                                     # unrelated deep structure; restricting
                                     # to the upper water column keeps these
                                     # metrics relevant to the euphotic,
                                     # bloom-relevant part of the profile.
PROFILE_BIN_M <- 1
G_ACCEL <- 9.81                   # m s-2

ctd_paths <- function(proj_root) {
  ctd_dir <- file.path(proj_root, "data", "lake_environmental_data", "ctd")
  list(
    data_file = file.path(ctd_dir, "ctd_casts", "ctd_data_2005_2025.csv"),
    info_file = file.path(ctd_dir, "ctd_casts", "ctd_info_2005_2025.csv"),
    kd_file = file.path(ctd_dir, "LTP_Kd_Thermocline_depth_results.csv"),
    schmidt_file = file.path(proj_root, "data", "processed", "ctd", "mltp_schmidt_stability_2005_2025.csv")
  )
}

#' Load quality-controlled Index (LTP) CTD profiles at 1 m depth bins.
#'
#' Cast-quality rules mirror scripts/supplemental_ctd_profiles_stability.R:
#' duplicate casts and casts flagged 4 (do not use) or 5 (missing) are
#' dropped, per-variable flags of 4 or 5 null only that variable.
load_ltp_ctd_profiles <- function(proj_root) {
  paths <- ctd_paths(proj_root)
  if (!file.exists(paths$data_file)) stop("Missing CTD data: ", paths$data_file)
  if (!file.exists(paths$info_file)) stop("Missing CTD cast info: ", paths$info_file)

  variables <- c("Temperature", "Density", "Chl_Fluorescence", "Phycocyanin", "PAR")
  ctd <- fread(
    paths$data_file,
    select = c(
      "CTD_ID", "Event_ID", "Station_ID", "Cast_Date_Time_Local",
      "Cast_Duplicate", "Cast_Flag", "Depth", variables
    ),
    na.strings = c("", "NA", "NULL"),
    showProgress = FALSE
  )
  ctd <- ctd[Station_ID == CTD_STATION]
  if (nrow(ctd) == 0) {
    stop("No CTD rows found for Station_ID == '", CTD_STATION, "'.")
  }

  flag_cols <- paste0(variables, "_Flag")
  ctd_info <- read_csv(paths$info_file, show_col_types = FALSE, progress = FALSE, lazy = FALSE) %>%
    transmute(
      CTD_ID,
      info_duplicate = as.integer(Duplicate),
      info_flag = as.integer(Flag),
      depth_flag = as.integer(Depth_Flag),
      across(all_of(flag_cols), ~suppressWarnings(as.integer(.x)))
    ) %>%
    as.data.table()

  ctd[, date := as.IDate(substr(Cast_Date_Time_Local, 1L, 10L))]
  ctd[, Cast_Duplicate := as.integer(Cast_Duplicate)]
  ctd[, Cast_Flag := as.integer(Cast_Flag)]
  ctd <- merge(ctd, ctd_info, by = "CTD_ID", all.x = TRUE, sort = FALSE)
  ctd[, duplicate_use := fcoalesce(Cast_Duplicate, info_duplicate, 0L)]
  ctd[, cast_flag_use := fcoalesce(info_flag, Cast_Flag)]
  ctd[, valid_cast :=
        duplicate_use == 0L &
        (is.na(cast_flag_use) | cast_flag_use <= 3L) &
        (is.na(depth_flag) | depth_flag <= 3L)]

  ctd <- ctd[valid_cast == TRUE & is.finite(Depth) & Depth >= 0]

  long_list <- lapply(variables, function(v) {
    fc <- paste0(v, "_Flag")
    x <- ctd[
      (is.na(get(fc)) | get(fc) <= 3L) & is.finite(get(v)),
      .(CTD_ID, Event_ID, date, Depth, value = get(v))
    ]
    if (nrow(x) == 0L) return(NULL)
    x[, depth_m := round(Depth / PROFILE_BIN_M) * PROFILE_BIN_M]
    out <- x[, .(value = mean(value)), by = .(CTD_ID, Event_ID, date, depth_m)]
    out[, variable := v]
    out
  })
  long <- rbindlist(long_list, use.names = TRUE, fill = TRUE)

  profiles <- long %>%
    as_tibble() %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    arrange(CTD_ID, depth_m) %>%
    mutate(date = as.Date(date))

  profiles
}

#' Per-cast maximum profile depth and a not-deep-enough flag.
flag_cast_depth <- function(profiles) {
  profiles %>%
    group_by(CTD_ID, Event_ID, date) %>%
    summarise(max_depth_m = max(depth_m), n_depths = n(), .groups = "drop") %>%
    mutate(insufficient_depth = max_depth_m < MIN_PROFILE_DEPTH_M)
}

#' Shallowest depth at which a threshold is exceeded relative to a
#' near-surface reference value (first two valid bins averaged).
threshold_depth <- function(depth_m, value, delta, direction = c("increase", "decrease")) {
  direction <- match.arg(direction)
  keep <- is.finite(depth_m) & is.finite(value)
  depth_m <- depth_m[keep]; value <- value[keep]
  if (length(depth_m) < 3) return(NA_real_)
  ord <- order(depth_m)
  depth_m <- depth_m[ord]; value <- value[ord]
  ref <- mean(value[1:min(2, length(value))])
  exceed <- if (direction == "increase") value - ref >= delta else ref - value >= delta
  if (!any(exceed)) return(NA_real_)
  depth_m[which(exceed)[1]]
}

#' Depth of the maximum absolute vertical gradient of a profile variable
#' (used for thermocline depth from density and, as a cross-check, from
#' rLakeAnalyzer's temperature-based thermo.depth()).
max_gradient_depth <- function(depth_m, value) {
  keep <- is.finite(depth_m) & is.finite(value)
  depth_m <- depth_m[keep]; value <- value[keep]
  if (length(depth_m) < 4) return(NA_real_)
  ord <- order(depth_m)
  depth_m <- depth_m[ord]; value <- value[ord]
  dz <- diff(depth_m)
  dv <- diff(value)
  grad <- dv / dz
  mid_depth <- depth_m[-1] - dz / 2
  mid_depth[which.max(abs(grad))]
}

#' Fluorescence-weighted depth centroid: sum(z * C) / sum(C), restricted to
#' non-negative signal (raw fluorescence occasionally reads small negative
#' values near the noise floor, which are set to zero for weighting only).
fluorescence_centroid <- function(depth_m, value) {
  keep <- is.finite(depth_m) & is.finite(value)
  depth_m <- depth_m[keep]; value <- pmax(value[keep], 0)
  if (sum(value) <= 0) return(NA_real_)
  sum(depth_m * value) / sum(value)
}

fluorescence_max_depth <- function(depth_m, value) {
  keep <- is.finite(depth_m) & is.finite(value)
  depth_m <- depth_m[keep]; value <- value[keep]
  if (length(value) == 0) return(NA_real_)
  depth_m[which.max(value)]
}

#' Maximum buoyancy frequency (N2, s-2) from the measured density profile.
max_buoyancy_frequency <- function(depth_m, density_kgm3) {
  keep <- is.finite(depth_m) & is.finite(density_kgm3)
  depth_m <- depth_m[keep]; density_kgm3 <- density_kgm3[keep]
  if (length(depth_m) < 4) return(NA_real_)
  ord <- order(depth_m)
  depth_m <- depth_m[ord]; density_kgm3 <- density_kgm3[ord]
  dz <- diff(depth_m)
  drho <- diff(density_kgm3)
  rho0 <- mean(density_kgm3)
  n2 <- (G_ACCEL / rho0) * (drho / dz)
  max(n2, na.rm = TRUE)
}

#' Compute the full per-cast metric table for Analysis B.
compute_ctd_cast_metrics <- function(profiles, kd_table) {
  depth_flags <- flag_cast_depth(profiles)

  metrics <- profiles %>%
    group_by(CTD_ID, Event_ID, date) %>%
    summarise(
      mld_density_m = threshold_depth(depth_m, Density, DELTA_RHO_MLD, "increase"),
      mld_temperature_m = threshold_depth(depth_m, Temperature, DELTA_T_MLD, "decrease"),
      thermocline_depth_density_m = max_gradient_depth(
        depth_m[depth_m <= SURFACE_LAYER_SEARCH_M], Density[depth_m <= SURFACE_LAYER_SEARCH_M]
      ),
      thermocline_depth_rla_m = tryCatch(
        suppressWarnings(rLakeAnalyzer::thermo.depth(
          Temperature[order(depth_m)][depth_m[order(depth_m)] <= SURFACE_LAYER_SEARCH_M],
          sort(depth_m)[sort(depth_m) <= SURFACE_LAYER_SEARCH_M]
        )),
        error = function(e) NA_real_
      ),
      n2_max_s2 = max_buoyancy_frequency(
        depth_m[depth_m <= SURFACE_LAYER_SEARCH_M], Density[depth_m <= SURFACE_LAYER_SEARCH_M]
      ),
      upper_density_diff_kgm3 = {
        surf <- mean(Density[depth_m <= 2], na.rm = TRUE)
        d20 <- Density[which.min(abs(depth_m - 20))]
        d20 - surf
      },
      integrated_density_gradient = {
        ord <- order(depth_m)
        dz <- diff(depth_m[ord]); drho <- diff(Density[ord])
        sum(abs(drho), na.rm = TRUE)
      },
      chl_max_depth_m = fluorescence_max_depth(depth_m, Chl_Fluorescence),
      chl_centroid_depth_m = fluorescence_centroid(depth_m, Chl_Fluorescence),
      phyco_max_depth_m = fluorescence_max_depth(depth_m, Phycocyanin),
      phyco_centroid_depth_m = fluorescence_centroid(depth_m, Phycocyanin),
      mixed_layer_mean_par = {
        mld <- threshold_depth(depth_m, Density, DELTA_RHO_MLD, "increase")
        if (is.na(mld)) NA_real_ else mean(PAR[depth_m <= mld], na.rm = TRUE)
      },
      density_range_kgm3 = diff(range(
        Density[depth_m <= SURFACE_LAYER_SEARCH_M], na.rm = TRUE
      )),
      .groups = "drop"
    ) %>%
    left_join(depth_flags, by = c("CTD_ID", "Event_ID", "date")) %>%
    left_join(kd_table, by = c("Event_ID" = "Event", "date" = "Date")) %>%
    mutate(
      euphotic_depth_m = if_else(is.finite(Kd_PAR) & Kd_PAR > 0, log(100) / Kd_PAR, NA_real_),
      zmix_over_zeu = mld_density_m / euphotic_depth_m,
      # A profile with a very small top-to-bottom density range is
      # effectively unstratified; a "depth of maximum gradient" under those
      # conditions reflects measurement noise rather than a true
      # thermocline, so thermocline/N2 metrics are suppressed.
      weakly_stratified = density_range_kgm3 < MIN_STRATIFICATION_RANGE_KGM3,
      across(
        c(thermocline_depth_density_m, thermocline_depth_rla_m, n2_max_s2),
        ~if_else(weakly_stratified, NA_real_, .x)
      ),
      across(
        c(mld_density_m, mld_temperature_m, thermocline_depth_density_m,
          thermocline_depth_rla_m, n2_max_s2, chl_max_depth_m, chl_centroid_depth_m,
          phyco_max_depth_m, phyco_centroid_depth_m, mixed_layer_mean_par,
          euphotic_depth_m, zmix_over_zeu),
        ~if_else(insufficient_depth, NA_real_, .x)
      )
    ) %>%
    arrange(date)

  metrics
}

#' Load the pre-computed LTP Kd (light-attenuation) regression results and
#' return the PAR attenuation coefficient per cast, keyed by Event/Date to
#' match compute_ctd_cast_metrics().
load_kd_table <- function(proj_root) {
  paths <- ctd_paths(proj_root)
  if (!file.exists(paths$kd_file)) {
    warning("Kd/euphotic-depth source not found: ", paths$kd_file,
            "; euphotic depth will be NA.")
    return(tibble(Event = character(), Date = as.Date(character()), Kd_PAR = numeric()))
  }
  read_csv(paths$kd_file, show_col_types = FALSE) %>%
    transmute(Event = Event, Date = as.Date(Date), Kd_PAR = as.numeric(Kd_PAR))
}

#' Change in mixed-layer depth and Schmidt-proxy stability between
#' consecutive casts, ordered by date.
compute_cast_to_cast_change <- function(metrics) {
  metrics %>%
    arrange(date) %>%
    mutate(
      delta_mld_density_m = mld_density_m - lag(mld_density_m),
      delta_n2_max_s2 = n2_max_s2 - lag(n2_max_s2),
      days_since_prior_cast = as.numeric(date - lag(date))
    )
}
