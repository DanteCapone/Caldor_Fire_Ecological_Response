# =============================================================================
# phytoplankton_gamm_anomaly.R
# Caldor Fire Ecosystem Response Project
#
# Analysis A: leave-one-year-out (LOYO) GAMM anomaly analysis quantifying how
# anomalous the 2021 Leptolyngbya spp. bloom was relative to the LTP
# 2005-2025 record.
#
# Response variable: Leptolyngbya spp. biovolume (mm3 m-3), Tweedie bam().
# See R/phytoplankton_gamm_helpers.R for the zero-handling and modeling
# rationale.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(tweedie)
  library(ragg)
})

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
source("R/phytoplankton_gamm_helpers.R")

data_dir <- file.path(proj_root, "data", "processed")
model_dir <- file.path(proj_root, "models", "phytoplankton_gamm")
fig_dir <- file.path(proj_root, "figures", "phytoplankton_anomaly")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples_file <- file.path(data_dir, "phytoplankton_lepto_samples.csv")
predictions_file <- file.path(data_dir, "phytoplankton_gamm_predictions.rds")
metrics_file <- file.path(data_dir, "annual_bloom_anomaly_metrics.csv")
model_spec_file <- file.path(model_dir, "model_specification.txt")
full_model_file <- file.path(model_dir, "full_model_2005_2025.rds")
loyo_models_dir <- file.path(model_dir, "loyo_fits")
dir.create(loyo_models_dir, recursive = TRUE, showWarnings = FALSE)

FOCAL_YEARS <- c(2021L, 2011L)
FIT_YEARS <- 2005:2025

# ---- Data preparation --------------------------------------------------------
samples <- build_leptolyngbya_samples(proj_root)

stopifnot(
  "Expected years 2005-2025 in the record" =
    all(FIT_YEARS %in% unique(samples$year))
)

sampling_effort <- summarize_sampling_effort(samples)
write_csv(sampling_effort, file.path(data_dir, "phytoplankton_sampling_effort_summary.csv"))

year_coverage <- flag_insufficient_years(samples)
write_csv(year_coverage, file.path(data_dir, "phytoplankton_year_coverage.csv"))
insufficient_years <- year_coverage %>% filter(!sufficient_coverage) %>% pull(year)
if (length(insufficient_years) > 0) {
  message(
    "Years with fewer than 4 sample dates (retained in model, flagged in table): ",
    paste(insufficient_years, collapse = ", ")
  )
}

write_csv(samples, samples_file)

bloom_threshold <- compute_bloom_threshold(samples)
message(
  "Historical (non-2021, <=20 m, 90th pct) bloom threshold: ",
  signif(bloom_threshold, 4), " mm3 m-3"
)

cat("\n=== Analysis-ready sample summary ===\n")
cat("Rows:", nrow(samples), "| distinct dates:", n_distinct(samples$date),
    "| eligible depths:", paste(attr(samples, "eligible_depths"), collapse = ","), "m\n")
cat("Excluded depths (< 5 sample dates):",
    paste(attr(samples, "excluded_depths"), collapse = ","), "\n")
cat("Zero fraction (biovolume == 0):", round(mean(samples$biovolume == 0), 3), "\n")
cat("Low-effort samples flagged:", sum(samples$low_effort_sample), "of", nrow(samples), "\n")

# ---- Full model (all years) for diagnostics and model comparison -----------
full_model <- fit_lepto_gamm(samples)
saveRDS(full_model, full_model_file)

alt_model_depth_re <- bam(
  biovolume ~
    s(doy_cyclic, bs = "cc", k = 10) +
    s(year, k = 8) +
    depth_factor +
    s(year_factor, bs = "re"),
  data = samples, family = tw(link = "log"),
  knots = list(doy_cyclic = c(0.5, 366.5)),
  method = "fREML", discrete = TRUE
)

model_comparison <- tibble(
  model = c("primary: s(depth_num) smooth", "alternative: depth_factor"),
  AIC = c(AIC(full_model), AIC(alt_model_depth_re)),
  deviance_explained = c(
    summary(full_model)$dev.expl, summary(alt_model_depth_re)$dev.expl
  )
)
write_csv(model_comparison, file.path(model_dir, "model_comparison.csv"))

writeLines(
  c(
    "Primary model formula:",
    "biovolume ~ s(doy_cyclic, bs='cc', k=10) + s(year, k=8) + s(depth_num, k=5) + s(year_factor, bs='re')",
    "Family: Tweedie (log link), fit via mgcv::bam(method='fREML', discrete=TRUE)",
    paste("Bloom threshold (90th pct, non-2021, <=20 m):", signif(bloom_threshold, 4), "mm3 m-3"),
    paste("Random seed:", GAMM_SEED),
    paste("mgcv version:", as.character(packageVersion("mgcv"))),
    paste("R version:", R.version.string)
  ),
  model_spec_file
)

# ---- Diagnostics --------------------------------------------------------------
diagnostics_file <- file.path(fig_dir, "gamm_diagnostics.png")
agg_png(diagnostics_file, width = 20, height = 16, units = "cm", res = LO_DPI)
par(mfrow = c(2, 2))
gam.check(full_model, pch = 19, cex = 0.4)
dev.off()

concurvity_table <- concurvity(full_model, full = FALSE)$estimate
write.csv(concurvity_table, file.path(model_dir, "concurvity_estimates.csv"))

cat("\nOverdispersion / basis-dimension check (see gam.check output above).\n")

# ---- Leave-one-year-out procedure -------------------------------------------
run_loyo_year <- function(focal_year) {
  train <- samples %>% filter(year != focal_year)
  focal <- samples %>% filter(year == focal_year)
  if (nrow(focal) == 0) return(NULL)

  model <- fit_lepto_gamm(train)
  saveRDS(model, file.path(loyo_models_dir, paste0("loyo_", focal_year, ".rds")))

  pred <- predict_focal_year(model, focal)
  pred$focal_year <- focal_year
  pred
}

message("Running leave-one-year-out fits for ", length(FIT_YEARS), " years...")
loyo_predictions <- map_dfr(FIT_YEARS, run_loyo_year)
saveRDS(loyo_predictions, predictions_file)

# ---- Annual anomaly metrics ---------------------------------------------------
annual_metrics <- compute_annual_anomaly_metrics(loyo_predictions, bloom_threshold)

rank_integrated <- rank_2021(annual_metrics, "integrated_positive_anomaly")
rank_duration <- rank_2021(annual_metrics, "bloom_duration_days")
rank_peak <- rank_2021(annual_metrics, "annual_peak_biovolume")
rank_above_pi <- rank_2021(annual_metrics, "integrated_anomaly_above_pi")

annual_metrics_ranked <- annual_metrics %>%
  left_join(
    rank_integrated$table %>% select(year, rank_integrated_anomaly = rank),
    by = "year"
  ) %>%
  left_join(
    rank_duration$table %>% select(year, rank_bloom_duration = rank),
    by = "year"
  ) %>%
  left_join(
    rank_peak$table %>% select(year, rank_peak_biovolume = rank),
    by = "year"
  ) %>%
  left_join(
    rank_above_pi$table %>% select(year, rank_anomaly_above_pi = rank),
    by = "year"
  )
write_csv(annual_metrics_ranked, metrics_file)

cat("\n=== 2021 percentile ranks across LOYO anomaly metrics ===\n")
cat("Integrated positive anomaly percentile:", round(rank_integrated$focal_percentile, 1), "\n")
cat("Integrated anomaly above 95% PI percentile:", round(rank_above_pi$focal_percentile, 1), "\n")
cat("Bloom duration percentile:", round(rank_duration$focal_percentile, 1), "\n")
cat("Annual peak biovolume percentile:", round(rank_peak$focal_percentile, 1), "\n")

# ---- Figures ------------------------------------------------------------------
plot_focal_year <- function(focal_year) {
  d <- loyo_predictions %>% filter(focal_year == !!focal_year)
  agg_by_date <- d %>%
    group_by(date) %>%
    summarise(
      biovolume = sum(biovolume, na.rm = TRUE),
      expected_mean = sum(expected_mean, na.rm = TRUE),
      ci_lower = sum(ci_lower, na.rm = TRUE),
      ci_upper = sum(ci_upper, na.rm = TRUE),
      pi_lower = sum(pi_lower, na.rm = TRUE),
      pi_upper = sum(pi_upper, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(agg_by_date, aes(x = date)) +
    geom_ribbon(aes(ymin = pi_lower, ymax = pi_upper), fill = "grey75", alpha = 0.4) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#3182BD", alpha = 0.35) +
    geom_line(aes(y = expected_mean), colour = "#3182BD", linewidth = LO_LW_THICK) +
    geom_point(aes(y = biovolume), colour = "#009E73", size = 1.6) +
    geom_line(aes(y = biovolume), colour = "#009E73", linewidth = LO_LW_MID) +
    labs(
      title = paste0(focal_year, " observed vs. leave-one-year-out expected trajectory"),
      x = NULL,
      y = expression("Leptolyngbya spp. biovolume (mm"^3*" m"^-3*", summed over depths)")
    ) +
    lo_theme(base_size = 8)
  save_lo_fig(
    p, file.path(fig_dir, paste0("loyo_trajectory_", focal_year, ".png")),
    width_type = "1half", height_cm = 9
  )
}
walk(FOCAL_YEARS, plot_focal_year)

overview_data <- loyo_predictions %>%
  group_by(focal_year, date) %>%
  summarise(
    biovolume = sum(biovolume, na.rm = TRUE),
    expected_mean = sum(expected_mean, na.rm = TRUE),
    pi_lower = sum(pi_lower, na.rm = TRUE),
    pi_upper = sum(pi_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(doy = yday(date))

p_overview <- ggplot(overview_data, aes(x = doy)) +
  geom_ribbon(aes(ymin = pi_lower, ymax = pi_upper), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = expected_mean), colour = "#3182BD", linewidth = LO_LW_MID) +
  geom_point(aes(y = biovolume), colour = "#009E73", size = 0.8) +
  facet_wrap(~focal_year, ncol = 5) +
  labs(
    title = "Observed vs. leave-one-year-out expected Leptolyngbya spp. biovolume, 2005-2025",
    x = "Day of year",
    y = expression("Biovolume (mm"^3*" m"^-3*")")
  ) +
  lo_theme(base_size = 7)
save_lo_fig(
  p_overview, file.path(fig_dir, "loyo_overview_all_years.png"),
  width_type = "double", height_cm = 15
)

p_rank_anomaly <- ggplot(
  rank_integrated$table,
  aes(x = reorder(factor(year), integrated_positive_anomaly), y = integrated_positive_anomaly)
) +
  geom_col(aes(fill = year == 2021)) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey60"), guide = "none") +
  labs(
    title = "Ranked annual integrated positive anomaly",
    x = NULL,
    y = expression("Integrated positive anomaly (mm"^3*" m"^-3*")")
  ) +
  lo_theme(base_size = 8)
save_lo_fig(
  p_rank_anomaly, file.path(fig_dir, "ranked_integrated_anomaly.png"),
  width_type = "1half", height_cm = 10
)

p_rank_duration_peak <- annual_metrics %>%
  mutate(is_2021 = year == 2021) %>%
  pivot_longer(
    cols = c(bloom_duration_days, annual_peak_biovolume),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      bloom_duration_days = "Bloom duration (days)",
      annual_peak_biovolume = "Peak biovolume (mm3 m-3)"
    )
  ) %>%
  ggplot(aes(x = reorder(factor(year), value), y = value)) +
  geom_col(aes(fill = is_2021)) +
  coord_flip() +
  facet_wrap(~metric, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey60"), guide = "none") +
  labs(title = "Ranked bloom duration and peak magnitude", x = NULL, y = NULL) +
  lo_theme(base_size = 8)
save_lo_fig(
  p_rank_duration_peak, file.path(fig_dir, "ranked_duration_peak.png"),
  width_type = "double", height_cm = 10
)

cat("\n=== Analysis A outputs written ===\n")
cat(" ", samples_file, "\n")
cat(" ", predictions_file, "\n")
cat(" ", metrics_file, "\n")
cat(" ", model_spec_file, "\n")
cat(" ", full_model_file, "\n")
cat(" ", fig_dir, "/*.png\n")
