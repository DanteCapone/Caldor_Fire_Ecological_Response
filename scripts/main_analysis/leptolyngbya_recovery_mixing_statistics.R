# =============================================================================
# Quantitative support for the association between winter mixing and the
# post-peak decline of Leptolyngbya in 2021-2022.
#
# This analysis estimates the timing and rate of the biological decline and
# quantifies its association with Schmidt stability. It does not estimate a
# causal effect of mixing: the two series are observational, seasonally
# structured, sampled approximately monthly, and measured at different lake
# stations (LTP phytoplankton; MLTP whole-lake stability).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(ragg)
})

source("scripts/figure_aesthetics.R")

LEPTO_FILE <- file.path("data", "processed", "phytoplankton_lepto_samples.csv")
FULL_PHYTO_FILE <- file.path(
  "data", "lake_environmental_data", "phytoplankton",
  "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"
)
STABILITY_FILE <- file.path(
  "data", "processed", "ctd", "mltp_lake_tools_stability_metrics.csv"
)
COMMUNITY_FILE <- file.path(
  "figures", "figure_5_nmds_braycurtis", "figure_5_panel_a_data.csv"
)
OUT_DIR <- file.path("figures", "manuscript_statistics", "leptolyngbya_recovery_mixing")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(LEPTO_FILE, FULL_PHYTO_FILE, STABILITY_FILE, COMMUNITY_FILE)
if (any(!file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
MODEL_START <- as.Date("2021-04-01")
MODEL_END <- as.Date("2022-06-30")
DISPLAY_START <- as.Date("2021-08-01")
DISPLAY_END <- as.Date("2022-06-30")
POST_SPRING_START <- as.Date("2022-05-31")
PRIMARY_MATCH_TOLERANCE_DAYS <- 14
MATCH_SENSITIVITY_DAYS <- c(7, 14, 21, 30)
BOOTSTRAP_REPLICATES <- 2000L
BOOTSTRAP_SEED <- 20260917L

integrate_profile <- function(depth_m, concentration) {
  keep <- is.finite(depth_m) & is.finite(concentration)
  depth_m <- depth_m[keep]
  concentration <- concentration[keep]
  if (length(depth_m) < 2L) {
    return(tibble(
      integrated_value = NA_real_, minimum_depth_m = NA_real_,
      maximum_depth_m = NA_real_, sampled_depths_n = length(depth_m)
    ))
  }
  ord <- order(depth_m)
  depth_m <- depth_m[ord]
  concentration <- concentration[ord]
  tibble(
    integrated_value = sum(
      diff(depth_m) * (head(concentration, -1L) + tail(concentration, -1L)) / 2
    ),
    minimum_depth_m = min(depth_m), maximum_depth_m = max(depth_m),
    sampled_depths_n = length(depth_m)
  )
}

lepto_samples <- read_csv(LEPTO_FILE, show_col_types = FALSE) %>%
  transmute(
    date = as.Date(date), depth_m = as.numeric(depth_num),
    abundance_cells_l = replace_na(as.numeric(abundance), 0),
    biovolume = replace_na(as.numeric(biovolume), 0)
  ) %>%
  filter(depth_m >= 5, depth_m <= 90)

lepto_profiles <- lepto_samples %>%
  group_by(date) %>%
  group_modify(~{
    abundance <- integrate_profile(.x$depth_m, .x$abundance_cells_l)
    biovolume <- integrate_profile(.x$depth_m, .x$biovolume)
    tibble(
      abundance_cells_m2 = abundance$integrated_value * 1000,
      biovolume_areal = biovolume$integrated_value,
      minimum_depth_m = abundance$minimum_depth_m,
      maximum_depth_m = abundance$maximum_depth_m,
      sampled_depths_n = abundance$sampled_depths_n
    )
  }) %>%
  ungroup() %>%
  arrange(date)

model_data <- lepto_profiles %>%
  filter(date >= MODEL_START, date <= MODEL_END, abundance_cells_m2 >= 0) %>%
  mutate(
    time_days = as.numeric(date - min(date)),
    log10_abundance = log10(abundance_cells_m2 + 1)
  )

aicc <- function(model) {
  n <- nobs(model)
  k <- length(coef(model)) + 1L
  AIC(model) + 2 * k * (k + 1) / (n - k - 1)
}

candidate_breaks <- model_data$date[4:(nrow(model_data) - 3L)]
breakpoint_fits <- map_dfr(candidate_breaks, function(candidate_date) {
  candidate_day <- as.numeric(candidate_date - min(model_data$date))
  fit_data <- model_data %>%
    mutate(post_break_days = pmax(0, time_days - candidate_day))
  fit <- lm(log10_abundance ~ time_days + post_break_days, data = fit_data)
  beta <- coef(fit)
  vc <- vcov(fit)
  post_slope <- beta[["time_days"]] + beta[["post_break_days"]]
  post_slope_se <- sqrt(
    vc["time_days", "time_days"] + vc["post_break_days", "post_break_days"] +
      2 * vc["time_days", "post_break_days"]
  )
  tibble(
    breakpoint_date = candidate_date,
    breakpoint_day = candidate_day,
    aicc = aicc(fit),
    pre_slope_log10_per_day = beta[["time_days"]],
    post_slope_log10_per_day = post_slope,
    post_slope_se = post_slope_se,
    r_squared = summary(fit)$r.squared
  )
}) %>%
  mutate(delta_aicc = aicc - min(aicc)) %>%
  arrange(aicc)

best_break <- breakpoint_fits %>% slice_head(n = 1)
best_break_day <- best_break$breakpoint_day
best_data <- model_data %>%
  mutate(post_break_days = pmax(0, time_days - best_break_day))
best_model <- lm(log10_abundance ~ time_days + post_break_days, data = best_data)
linear_model <- lm(log10_abundance ~ time_days, data = model_data)

best_beta <- coef(best_model)
best_vcov <- vcov(best_model)
pre_slope <- best_beta[["time_days"]]
post_slope <- best_beta[["time_days"]] + best_beta[["post_break_days"]]
post_slope_se <- sqrt(
  best_vcov["time_days", "time_days"] +
    best_vcov["post_break_days", "post_break_days"] +
    2 * best_vcov["time_days", "post_break_days"]
)
post_slope_df <- df.residual(best_model)
post_slope_ci <- post_slope + c(-1, 1) * qt(0.975, post_slope_df) * post_slope_se
post_slope_p <- 2 * pt(-abs(post_slope / post_slope_se), df = post_slope_df)
half_life_days <- if (post_slope < 0) log10(0.5) / post_slope else NA_real_

set.seed(BOOTSTRAP_SEED)
bootstrap_breaks <- map_dfr(seq_len(BOOTSTRAP_REPLICATES), function(iteration) {
  simulated <- model_data %>%
    mutate(
      log10_abundance = fitted(best_model) +
        sample(residuals(best_model), replace = TRUE)
    )
  fits <- map_dfr(candidate_breaks, function(candidate_date) {
    candidate_day <- as.numeric(candidate_date - min(simulated$date))
    d <- simulated %>% mutate(post_break_days = pmax(0, time_days - candidate_day))
    fit <- lm(log10_abundance ~ time_days + post_break_days, data = d)
    beta <- coef(fit)
    tibble(
      breakpoint_date = candidate_date,
      aicc = aicc(fit),
      post_slope = beta[["time_days"]] + beta[["post_break_days"]]
    )
  })
  fits %>% slice_min(aicc, n = 1, with_ties = FALSE) %>% mutate(iteration = iteration)
})

bootstrap_break_ci <- quantile(
  as.numeric(bootstrap_breaks$breakpoint_date), c(0.025, 0.975), na.rm = TRUE
) %>%
  as.Date(origin = "1970-01-01")
bootstrap_half_life_ci <- bootstrap_breaks %>%
  filter(post_slope < 0) %>%
  transmute(half_life_days = log10(0.5) / post_slope) %>%
  summarise(
    lower = quantile(half_life_days, 0.025, na.rm = TRUE),
    upper = quantile(half_life_days, 0.975, na.rm = TRUE)
  )

segmented_summary <- tibble(
  model_start = min(model_data$date), model_end = max(model_data$date),
  observations_n = nrow(model_data), response = "5-90 m depth-integrated Leptolyngbya abundance",
  transformation = "log10(cells m^-2 + 1)",
  breakpoint_date = best_break$breakpoint_date,
  breakpoint_bootstrap_95_lower = bootstrap_break_ci[1],
  breakpoint_bootstrap_95_upper = bootstrap_break_ci[2],
  pre_slope_log10_per_30_days = pre_slope * 30,
  post_slope_log10_per_30_days = post_slope * 30,
  post_slope_95_lower_per_30_days = post_slope_ci[1] * 30,
  post_slope_95_upper_per_30_days = post_slope_ci[2] * 30,
  post_slope_p = post_slope_p,
  half_life_days = half_life_days,
  half_life_bootstrap_95_lower = bootstrap_half_life_ci$lower,
  half_life_bootstrap_95_upper = bootstrap_half_life_ci$upper,
  segmented_aicc = aicc(best_model), linear_aicc = aicc(linear_model),
  delta_aicc_linear_minus_segmented = aicc(linear_model) - aicc(best_model),
  segmented_r_squared = summary(best_model)$r.squared,
  inference = paste(
    "Breakpoint selected by minimum AICc across observed interior dates;",
    "conditional slope p-value and residual-bootstrap intervals are descriptive",
    "for this single observational event"
  )
)

stability <- read_csv(STABILITY_FILE, show_col_types = FALSE) %>%
  transmute(
    stability_date = as.Date(date),
    schmidt_stability_j_m2 = as.numeric(schmidt_stability_j_m2)
  ) %>%
  filter(is.finite(schmidt_stability_j_m2)) %>%
  arrange(stability_date)

match_to_stability <- function(profile_data, tolerance_days) {
  map_dfr(seq_len(nrow(profile_data)), function(i) {
    differences <- abs(as.numeric(stability$stability_date - profile_data$date[i]))
    nearest <- which.min(differences)
    if (length(nearest) == 0L || differences[nearest] > tolerance_days) return(tibble())
    bind_cols(
      profile_data[i, ],
      stability[nearest, ],
      tibble(match_difference_days = differences[nearest])
    )
  })
}

peak_profile <- lepto_profiles %>%
  filter(date >= FIRE_START, date <= MODEL_END) %>%
  slice_max(abundance_cells_m2, n = 1, with_ties = FALSE)

association_by_tolerance <- map_dfr(MATCH_SENSITIVITY_DAYS, function(tolerance_days) {
  matched <- lepto_profiles %>%
    filter(date >= peak_profile$date, date <= as.Date("2022-12-31")) %>%
    match_to_stability(tolerance_days) %>%
    filter(is.finite(abundance_cells_m2), is.finite(schmidt_stability_j_m2))
  if (nrow(matched) < 4L) {
    return(tibble(
      tolerance_days = tolerance_days, matched_n = nrow(matched),
      spearman_rho = NA_real_, spearman_p = NA_real_,
      first_difference_rho = NA_real_, first_difference_p = NA_real_
    ))
  }
  level_test <- suppressWarnings(cor.test(
    matched$abundance_cells_m2, matched$schmidt_stability_j_m2,
    method = "spearman", exact = FALSE
  ))
  differenced <- matched %>%
    arrange(date) %>%
    mutate(
      abundance_rate = (log10(abundance_cells_m2 + 1) -
                          lag(log10(abundance_cells_m2 + 1))) /
        as.numeric(date - lag(date)),
      stability_rate = (log10(schmidt_stability_j_m2) -
                          lag(log10(schmidt_stability_j_m2))) /
        as.numeric(date - lag(date))
    ) %>%
    filter(is.finite(abundance_rate), is.finite(stability_rate))
  difference_test <- if (nrow(differenced) >= 4L) {
    suppressWarnings(cor.test(
      differenced$abundance_rate, differenced$stability_rate,
      method = "spearman", exact = FALSE
    ))
  } else NULL
  tibble(
    tolerance_days = tolerance_days,
    matched_n = nrow(matched),
    spearman_rho = unname(level_test$estimate), spearman_p = level_test$p.value,
    first_difference_n = nrow(differenced),
    first_difference_rho = if (is.null(difference_test)) NA_real_ else unname(difference_test$estimate),
    first_difference_p = if (is.null(difference_test)) NA_real_ else difference_test$p.value
  )
}) %>%
  mutate(
    inference = paste(
      "Exploratory paired-date association; no multiplicity correction;",
      "first differences reduce shared trend but do not remove seasonal confounding or autocorrelation"
    )
  )

primary_matches <- lepto_profiles %>%
  filter(date >= peak_profile$date, date <= as.Date("2022-12-31")) %>%
  match_to_stability(PRIMARY_MATCH_TOLERANCE_DAYS) %>%
  filter(is.finite(abundance_cells_m2), is.finite(schmidt_stability_j_m2)) %>%
  arrange(date)

oct_stability <- stability %>%
  filter(year(stability_date) == 2021L, month(stability_date) == 10L) %>%
  slice_min(abs(day(stability_date) - 1L), n = 1)
dec_stability <- stability %>%
  filter(year(stability_date) == 2021L, month(stability_date) == 12L) %>%
  slice_min(abs(day(stability_date) - 15L), n = 1)
stability_decline_percent <- 100 * (
  oct_stability$schmidt_stability_j_m2 - dec_stability$schmidt_stability_j_m2
) / oct_stability$schmidt_stability_j_m2

full_phyto <- read_csv(FULL_PHYTO_FILE, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date), depth_m = as.numeric(depth_num),
    abundance = replace_na(as.numeric(abundance), 0),
    biovolume = replace_na(as.numeric(biovolume), 0),
    is_leptolyngbya = str_detect(taxon, regex("^Leptolyngbya", ignore_case = TRUE))
  )

upper_community <- full_phyto %>%
  filter(depth_bin == "0-40 m") %>%
  group_by(date) %>%
  summarise(
    total_abundance = sum(abundance), total_biovolume = sum(biovolume),
    lepto_abundance = sum(abundance[is_leptolyngbya]),
    lepto_biovolume = sum(biovolume[is_leptolyngbya]),
    .groups = "drop"
  ) %>%
  mutate(
    lepto_abundance_percent = 100 * lepto_abundance / total_abundance,
    lepto_biovolume_percent = 100 * lepto_biovolume / total_biovolume
  )

post_spring_upper_max <- upper_community %>%
  filter(date > POST_SPRING_START, lepto_abundance > 0) %>%
  summarise(
    abundance_max_percent = max(lepto_abundance_percent, na.rm = TRUE),
    abundance_max_date = date[which.max(lepto_abundance_percent)],
    biovolume_max_percent = max(lepto_biovolume_percent, na.rm = TRUE),
    biovolume_max_date = date[which.max(lepto_biovolume_percent)]
  ) %>%
  mutate(
    scope = "sum across observed 0-40 m LTP depths on each sampling date",
    inference = "descriptive maximum; no significance test"
  )

community_recovery <- read_csv(COMMUNITY_FILE, show_col_types = FALSE) %>%
  filter(metric == "Community departure", stratum == "0-40 m") %>%
  transmute(
    community_peak_date = as.Date(peak_date),
    community_recovery_date = as.Date(recovery_date),
    community_recovery_months = as.numeric(recovery_months),
    recovery_rule
  ) %>%
  distinct() %>%
  slice_head(n = 1)

event_summary <- tibble(
  statistic = c(
    "October 2021 Schmidt stability", "December 2021 Schmidt stability",
    "October-December Schmidt stability decline", "Peak depth-integrated Leptolyngbya date",
    "Segmented-regression decline breakpoint", "Post-break abundance half-life",
    "Upper-community departure peak date", "Upper-community recovery date",
    "Maximum post-spring upper-zone abundance share",
    "Maximum post-spring upper-zone biovolume share"
  ),
  value = c(
    oct_stability$schmidt_stability_j_m2, dec_stability$schmidt_stability_j_m2,
    stability_decline_percent, as.numeric(peak_profile$date),
    as.numeric(best_break$breakpoint_date), half_life_days,
    as.numeric(community_recovery$community_peak_date),
    as.numeric(community_recovery$community_recovery_date),
    post_spring_upper_max$abundance_max_percent,
    post_spring_upper_max$biovolume_max_percent
  ),
  value_date = c(
    as.character(oct_stability$stability_date), as.character(dec_stability$stability_date),
    NA_character_, as.character(peak_profile$date), as.character(best_break$breakpoint_date),
    NA_character_, as.character(community_recovery$community_peak_date),
    as.character(community_recovery$community_recovery_date),
    as.character(post_spring_upper_max$abundance_max_date),
    as.character(post_spring_upper_max$biovolume_max_date)
  ),
  units_or_scope = c(
    "J m^-2", "J m^-2", "%", "date; 5-90 m trapezoidal integration",
    "date; minimum-AICc continuous segmented regression", "days",
    "date; season-specific Bray-Curtis distance from historical centroid",
    paste0("date; first two consecutive observations in reference envelope; ",
           community_recovery$community_recovery_months, " months after ignition"),
    "% of summed 0-40 m community abundance", "% of summed 0-40 m community biovolume"
  )
)

write_csv(lepto_profiles, file.path(OUT_DIR, "leptolyngbya_depth_integrated_timeseries.csv"))
write_csv(breakpoint_fits, file.path(OUT_DIR, "segmented_breakpoint_candidate_fits.csv"))
write_csv(segmented_summary, file.path(OUT_DIR, "segmented_decline_summary.csv"))
write_csv(bootstrap_breaks, file.path(OUT_DIR, "segmented_decline_bootstrap.csv"))
write_csv(primary_matches, file.path(OUT_DIR, "stability_abundance_primary_matches.csv"))
write_csv(association_by_tolerance, file.path(OUT_DIR, "stability_association_sensitivity.csv"))
write_csv(post_spring_upper_max, file.path(OUT_DIR, "post_spring_upper_community_share.csv"))
write_csv(event_summary, file.path(OUT_DIR, "manuscript_event_summary.csv"))

display_abundance_min <- lepto_profiles %>%
  filter(date >= DISPLAY_START, date <= DISPLAY_END, abundance_cells_m2 > 0) %>%
  summarise(value = min(abundance_cells_m2) * 0.8) %>%
  pull(value)

plot_stability <- stability %>%
  filter(stability_date >= DISPLAY_START, stability_date <= DISPLAY_END) %>%
  ggplot(aes(stability_date, schmidt_stability_j_m2 / 1000)) +
  annotate(
    "rect", xmin = FIRE_START, xmax = FIRE_END, ymin = -Inf, ymax = Inf,
    fill = "firebrick", alpha = 0.08
  ) +
  geom_line(colour = "#6A3D9A", linewidth = LO_LW_THICK) +
  geom_point(shape = 21, fill = "white", colour = "#6A3D9A", size = 1.8, stroke = 0.4) +
  geom_vline(xintercept = best_break$breakpoint_date, linetype = "dashed", linewidth = LO_LW_MID) +
  scale_x_date(limits = c(DISPLAY_START, DISPLAY_END), date_breaks = "2 months", date_labels = "%b\n%Y") +
  labs(x = NULL, y = expression("Schmidt stability (kJ m"^{-2}*")")) +
  lo_theme(base_size = 8.5) +
  labs(tag = "a") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.tag = element_text(size = LO_FS_TAG, face = "bold", family = LO_FONT),
    plot.tag.position = "topleft"
  )

plot_abundance <- lepto_profiles %>%
  filter(date >= DISPLAY_START, date <= DISPLAY_END) %>%
  ggplot(aes(date, abundance_cells_m2)) +
  annotate(
    "rect", xmin = FIRE_START, xmax = FIRE_END,
    ymin = display_abundance_min, ymax = Inf,
    fill = "firebrick", alpha = 0.08
  ) +
  annotate(
    "rect", xmin = bootstrap_break_ci[1], xmax = bootstrap_break_ci[2],
    ymin = display_abundance_min, ymax = Inf, fill = "grey50", alpha = 0.10
  ) +
  geom_line(colour = "#009E73", linewidth = LO_LW_THICK) +
  geom_point(shape = 21, fill = "white", colour = "#009E73", size = 1.8, stroke = 0.4) +
  geom_vline(xintercept = best_break$breakpoint_date, linetype = "dashed", linewidth = LO_LW_MID) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_x_date(limits = c(DISPLAY_START, DISPLAY_END), date_breaks = "2 months", date_labels = "%b\n%Y") +
  labs(
    x = NULL, y = "Leptolyngbya abundance\n(cells m⁻²; 5–90 m)", tag = "b"
  ) +
  lo_theme(base_size = 8.5) +
  theme(
    plot.tag = element_text(size = LO_FS_TAG, face = "bold", family = LO_FONT),
    plot.tag.position = "topleft"
  )

combined_plot <- plot_stability / plot_abundance

png_file <- file.path(OUT_DIR, "leptolyngbya_recovery_mixing_timing.png")
pdf_file <- file.path(OUT_DIR, "leptolyngbya_recovery_mixing_timing.pdf")
save_lo_fig(combined_plot, png_file, width_type = "1half", height_cm = 12.5)
ggsave(
  pdf_file, combined_plot, width = LO_WIDTH_1HALF, height = 12.5, units = "cm",
  device = cairo_pdf, family = LO_FONT
)

caption <- paste0(
  "Temporal association between water-column stability and the 2021-2022 decline of ",
  "Leptolyngbya. (a) Schmidt stability calculated from deep MLTP CTD profiles. ",
  "(b) Trapezoidally integrated Leptolyngbya abundance across the sampled 5-90 m LTP ",
  "profile. Firebrick shading marks the Caldor Fire window. The dashed line is the ",
  "minimum-AICc breakpoint in a continuous segmented regression of log10 abundance; ",
  "grey shading is its residual-bootstrap 95% interval. Sampling dates differ between ",
  "panels, and the temporal association does not by itself establish a causal mixing effect."
)
writeLines(caption, file.path(OUT_DIR, "leptolyngbya_recovery_mixing_caption.txt"), useBytes = TRUE)

cat("\n=== Leptolyngbya recovery and mixing statistics ===\n")
print(segmented_summary)
cat("\nStability association sensitivity:\n")
print(association_by_tolerance)
cat("\nPost-spring upper-community maxima:\n")
print(post_spring_upper_max)
cat("\nOutputs written to:", OUT_DIR, "\n")
