# =============================================================================
# Final manuscript analysis of the Leptolyngbya rise-to-decline transition.
#
# A REML GAM estimates a smooth abundance trajectory and its turning date. The
# aligned Schmidt-stability panel is descriptive: it shows temporal coincidence
# with winter destabilization but is not a causal test of mixing.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(MASS)
  library(patchwork)
  library(ragg)
})

source("scripts/figure_aesthetics.R")

PROFILE_FILE <- file.path(
  "figures", "manuscript_statistics", "leptolyngbya_recovery_mixing",
  "leptolyngbya_depth_integrated_timeseries.csv"
)
STABILITY_FILE <- file.path(
  "data", "processed", "ctd", "mltp_lake_tools_stability_metrics.csv"
)
MIXED_LAYER_FILE <- file.path(
  "data", "processed", "mltp_mixed_layer_depth_density_2005_2025.csv"
)
OUT_DIR <- file.path(
  "figures", "supplemental", "leptolyngbya_gam_winter_destabilization"
)
FINAL_DIR <- file.path("figures", "supplemental_final")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PROFILE_FILE) || !file.exists(STABILITY_FILE) ||
    !file.exists(MIXED_LAYER_FILE)) {
  stop("Missing required input. Run the Leptolyngbya integration workflow first.")
}

MODEL_START <- as.Date("2021-04-01")
MODEL_END <- as.Date("2022-06-30")
DISPLAY_START <- as.Date("2021-08-01")
DISPLAY_END <- as.Date("2022-06-30")
FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
SIMULATIONS <- 4000L
SIMULATION_SEED <- 20260919L

model_data <- read_csv(PROFILE_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= MODEL_START, date <= MODEL_END, abundance_cells_m2 >= 0) %>%
  arrange(date) %>%
  mutate(
    time_days = as.numeric(date - min(date)),
    log10_abundance = log10(abundance_cells_m2 + 1)
  )

# k = 5 limits flexibility for the 15-date series. REML estimates smoothing.
gam_model <- gam(
  log10_abundance ~ s(time_days, k = 5, bs = "tp"),
  data = model_data,
  method = "REML"
)

prediction_grid <- tibble(
  date = seq(min(model_data$date), max(model_data$date), by = "1 day")
) %>%
  mutate(time_days = as.numeric(date - min(model_data$date)))

prediction <- predict(gam_model, newdata = prediction_grid, se.fit = TRUE)
prediction_grid <- prediction_grid %>%
  mutate(
    estimate_log10 = as.numeric(prediction$fit),
    se_log10 = as.numeric(prediction$se.fit),
    lower_log10 = estimate_log10 - 1.96 * se_log10,
    upper_log10 = estimate_log10 + 1.96 * se_log10,
    estimate_cells_m2 = pmax(10^estimate_log10 - 1, 0),
    lower_cells_m2 = pmax(10^lower_log10 - 1, 0),
    upper_cells_m2 = pmax(10^upper_log10 - 1, 0)
  )

turning_row <- prediction_grid %>%
  slice_max(estimate_log10, n = 1, with_ties = FALSE)
turning_date <- turning_row$date

# Simulate coefficient uncertainty from the fitted GAM, preserving its smooth
# structure. The resulting interval describes uncertainty in the fitted peak
# date conditional on the chosen model basis and smoothing specification.
set.seed(SIMULATION_SEED)
coefficient_draws <- mvrnorm(
  n = SIMULATIONS,
  mu = coef(gam_model),
  Sigma = vcov(gam_model, unconditional = TRUE)
)
prediction_matrix <- predict(gam_model, newdata = prediction_grid, type = "lpmatrix")
simulated_curves <- prediction_matrix %*% t(coefficient_draws)
simulated_peak_indices <- apply(simulated_curves, 2, which.max)
simulated_turning_dates <- prediction_grid$date[simulated_peak_indices]
turning_interval <- quantile(
  as.numeric(simulated_turning_dates), c(0.025, 0.975), na.rm = TRUE
) %>%
  as.Date(origin = "1970-01-01")

stability <- read_csv(STABILITY_FILE, show_col_types = FALSE) %>%
  transmute(
    date = as.Date(date),
    schmidt_stability_j_m2 = as.numeric(schmidt_stability_j_m2)
  ) %>%
  filter(is.finite(schmidt_stability_j_m2)) %>%
  arrange(date)

# Mixed-layer depth uses the project's density-threshold definition: the first
# depth where density is 0.1 kg m^-3 greater than its near-surface reference.
# It is plotted downward from the top of the secondary axis so deeper mixing is
# visually lower, consistent with the other manuscript-facing depth panels.
mixed_layer <- read_csv(MIXED_LAYER_FILE, show_col_types = FALSE) %>%
  transmute(date = as.Date(date), mixed_layer_depth_m = as.numeric(mld_m)) %>%
  filter(is.finite(mixed_layer_depth_m)) %>%
  arrange(date)

stability_display <- stability %>%
  filter(date >= DISPLAY_START, date <= DISPLAY_END) %>%
  mutate(schmidt_stability_kj_m2 = schmidt_stability_j_m2 / 1000)

mixed_layer_display <- mixed_layer %>%
  filter(date >= DISPLAY_START, date <= DISPLAY_END)

if (!identical(stability_display$date, mixed_layer_display$date)) {
  stop("Expected mixed-layer depth on the same dates as Schmidt stability.")
}

STABILITY_AXIS_MAX_KJ_M2 <- 72
MIXED_LAYER_AXIS_MAX_M <- 30
if (any(mixed_layer_display$mixed_layer_depth_m < 0 |
        mixed_layer_display$mixed_layer_depth_m > MIXED_LAYER_AXIS_MAX_M)) {
  stop("Mixed-layer depth lies outside the configured 0-30 m secondary axis.")
}
mixed_layer_display <- mixed_layer_display %>%
  mutate(
    plot_y = STABILITY_AXIS_MAX_KJ_M2 *
      (1 - mixed_layer_depth_m / MIXED_LAYER_AXIS_MAX_M)
  )

october_stability <- stability %>% filter(date == as.Date("2021-10-01"))
december_stability <- stability %>% filter(date == as.Date("2021-12-21"))
if (nrow(october_stability) != 1L || nrow(december_stability) != 1L) {
  stop("Expected Schmidt-stability observations on 2021-10-01 and 2021-12-21.")
}
stability_decline_percent <- 100 * (
  october_stability$schmidt_stability_j_m2 -
    december_stability$schmidt_stability_j_m2
) / october_stability$schmidt_stability_j_m2

end_row <- prediction_grid %>% slice_tail(n = 1)
fitted_decline_percent <- 100 * (
  turning_row$estimate_cells_m2 - end_row$estimate_cells_m2
) / turning_row$estimate_cells_m2

smooth_table <- summary(gam_model)$s.table
gam_summary <- tibble(
  model = "Gaussian GAM of log10(cells m^-2 + 1)",
  smoothing_method = "REML",
  smooth_basis = "thin-plate regression spline",
  basis_dimension_k = 5L,
  observations_n = nrow(model_data),
  model_start = min(model_data$date),
  model_end = max(model_data$date),
  smooth_edf = unname(smooth_table[1, "edf"]),
  smooth_f = unname(smooth_table[1, "F"]),
  smooth_p = unname(smooth_table[1, "p-value"]),
  gam_identified_inflection_date = turning_date,
  inflection_date_simulation_95_lower = turning_interval[1],
  inflection_date_simulation_95_upper = turning_interval[2],
  coefficient_simulations = SIMULATIONS,
  simulation_seed = SIMULATION_SEED,
  fitted_peak_abundance_cells_m2 = turning_row$estimate_cells_m2,
  fitted_end_abundance_cells_m2 = end_row$estimate_cells_m2,
  fitted_peak_to_end_decline_percent = fitted_decline_percent,
  october_to_december_stability_decline_percent = stability_decline_percent,
  mixed_layer_depth_definition = paste(
    "First depth where density is 0.1 kg m^-3 above the near-surface reference;",
    "displayed on a 0-30 m secondary axis with depth increasing downward."
  ),
  inference = paste(
    "The turning date is the maximum of the daily fitted GAM curve;",
    "its simulation interval is conditional on the GAM specification.",
    "The aligned stability trajectory shows temporal coincidence, not causation."
  )
)

write_csv(gam_summary, file.path(OUT_DIR, "leptolyngbya_gam_transition_summary.csv"))
write_csv(prediction_grid, file.path(OUT_DIR, "leptolyngbya_gam_predictions.csv"))
write_csv(
  stability_display %>%
    dplyr::select(date, schmidt_stability_kj_m2) %>%
    left_join(
      mixed_layer_display %>% dplyr::select(date, mixed_layer_depth_m),
      by = "date"
    ),
  file.path(OUT_DIR, "schmidt_stability_mixed_layer_depth_plot_data.csv")
)
write_csv(
  tibble(simulation = seq_len(SIMULATIONS), turning_date = simulated_turning_dates),
  file.path(OUT_DIR, "leptolyngbya_gam_turning_date_simulations.csv")
)

transition_band <- tibble(
  xmin = turning_interval[1], xmax = turning_interval[2]
)

fire_midpoint <- FIRE_START + floor(as.numeric(FIRE_END - FIRE_START) / 2)

plot_stability <- ggplot() +
  annotate(
    "rect", xmin = FIRE_START, xmax = FIRE_END, ymin = -Inf, ymax = Inf,
    fill = "firebrick", alpha = 0.07
  ) +
  geom_rect(
    data = transition_band,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey55", alpha = 0.10
  ) +
  geom_line(
    data = stability_display,
    aes(date, schmidt_stability_kj_m2),
    colour = "#6A3D9A", linewidth = LO_LW_THICK
  ) +
  geom_point(
    data = stability_display,
    aes(date, schmidt_stability_kj_m2),
    shape = 21, fill = "white", colour = "#6A3D9A",
    size = 1.8, stroke = 0.4
  ) +
  geom_line(
    data = mixed_layer_display,
    aes(date, plot_y),
    colour = "#D55E00", linetype = "22", linewidth = LO_LW_THICK
  ) +
  geom_point(
    data = mixed_layer_display,
    aes(date, plot_y),
    shape = 24, fill = "white", colour = "#D55E00",
    size = 2, stroke = 0.45
  ) +
  geom_vline(
    xintercept = turning_date, linetype = "dashed", linewidth = LO_LW_MID
  ) +
  annotate(
    "text", x = fire_midpoint, y = 69,
    label = "Caldor Fire\nwindow", size = 2.25, family = LO_FONT,
    lineheight = 0.9
  ) +
  scale_y_continuous(
    limits = c(0, STABILITY_AXIS_MAX_KJ_M2),
    breaks = c(0, 20, 40, 60),
    name = expression("Schmidt stability (kJ m"^{-2}*")"),
    sec.axis = sec_axis(
      ~ (STABILITY_AXIS_MAX_KJ_M2 - .) *
        MIXED_LAYER_AXIS_MAX_M / STABILITY_AXIS_MAX_KJ_M2,
      name = "Mixed-layer depth (m)", breaks = c(0, 10, 20, 30)
    )
  ) +
  scale_x_date(
    limits = c(DISPLAY_START, DISPLAY_END),
    date_breaks = "2 months", date_labels = "%b\n%Y"
  ) +
  labs(x = NULL, tag = "a") +
  lo_theme(base_size = 8.5) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y.left = element_text(colour = "#6A3D9A"),
    axis.text.y.left = element_text(colour = "#6A3D9A"),
    axis.ticks.y.left = element_line(colour = "#6A3D9A"),
    axis.title.y.right = element_text(colour = "#D55E00"),
    axis.text.y.right = element_text(colour = "#D55E00"),
    axis.ticks.y.right = element_line(colour = "#D55E00"),
    plot.tag.position = c(0.985, 0.97),
    plot.tag = element_text(hjust = 1, vjust = 1, face = "bold", family = LO_FONT)
  )

plot_abundance <- ggplot() +
  annotate(
    "rect", xmin = FIRE_START, xmax = FIRE_END, ymin = 1e7, ymax = Inf,
    fill = "firebrick", alpha = 0.07
  ) +
  geom_rect(
    data = transition_band,
    aes(xmin = xmin, xmax = xmax, ymin = 1e7, ymax = Inf),
    inherit.aes = FALSE, fill = "grey55", alpha = 0.10
  ) +
  geom_ribbon(
    data = filter(prediction_grid, date >= DISPLAY_START, date <= DISPLAY_END),
    aes(date, ymin = lower_cells_m2, ymax = upper_cells_m2),
    fill = "#009E73", alpha = 0.16
  ) +
  geom_line(
    data = filter(prediction_grid, date >= DISPLAY_START, date <= DISPLAY_END),
    aes(date, estimate_cells_m2), colour = "#009E73", linewidth = 0.9
  ) +
  geom_point(
    data = filter(model_data, date >= DISPLAY_START, date <= DISPLAY_END),
    aes(date, abundance_cells_m2), shape = 21, fill = "white",
    colour = "black", size = 2, stroke = 0.4
  ) +
  geom_vline(
    xintercept = turning_date, linetype = "dashed", linewidth = LO_LW_MID
  ) +
  scale_y_log10(
    breaks = 10^(7:11),
    labels = function(x) parse(text = paste0("10^", round(log10(x))))
  ) +
  scale_x_date(
    limits = c(DISPLAY_START, DISPLAY_END),
    date_breaks = "2 months", date_labels = "%b\n%Y"
  ) +
  labs(
    x = NULL,
    y = expression(italic("Leptolyngbya") * " abundance" ~ ("cells " * m^{-2}) * "; 5-90 m"),
    tag = "b"
  ) +
  lo_theme(base_size = 8.5) +
  theme(
    plot.tag.position = c(0.985, 0.97),
    plot.tag = element_text(hjust = 1, vjust = 1, face = "bold", family = LO_FONT)
  )

supplemental_figure <- plot_stability / plot_abundance +
  plot_layout(heights = c(0.85, 1.15))

png_file <- file.path(OUT_DIR, "leptolyngbya_gam_winter_destabilization.png")
pdf_file <- file.path(OUT_DIR, "leptolyngbya_gam_winter_destabilization.pdf")
save_lo_fig(supplemental_figure, png_file, width_type = "1half", height_cm = 12.5)
ggsave(
  pdf_file, supplemental_figure, width = LO_WIDTH_1HALF, height = 12.5,
  units = "cm", device = cairo_pdf, family = LO_FONT
)

caption <- paste0(
  "Winter destabilization and the Leptolyngbya abundance transition. ",
  "(a) Schmidt stability (purple solid line and circles; left axis) and ",
  "density-threshold mixed-layer depth (orange dashed line and triangles; ",
  "right axis, with depth increasing downward) calculated from MLTP CTD ",
  "profiles. (b) Observed 5-90 m ",
  "depth-integrated LTP Leptolyngbya abundance and the REML GAM fit (green line; ",
  "pointwise 95% confidence band). The dashed line marks the GAM-identified ",
  "inflection point, operationally defined as the fitted maximum, on ",
  format(turning_date, "%d %B %Y"), " (+20/-16 days). Grey shading is its ",
  "simulation-based 95% interval; red shading marks the Caldor Fire window. ",
  "Phytoplankton was measured at LTP, whereas stability and mixed-layer depth ",
  "were measured at MLTP; their alignment indicates temporal coincidence rather ",
  "than causation."
)
writeLines(caption, file.path(OUT_DIR, "leptolyngbya_gam_winter_destabilization_caption.txt"))

final_png <- file.path(FINAL_DIR, "Figure_S12_leptolyngbya_gam_winter_destabilization.png")
final_pdf <- file.path(FINAL_DIR, "Figure_S12_leptolyngbya_gam_winter_destabilization.pdf")
final_caption <- file.path(
  FINAL_DIR, "Figure_S12_leptolyngbya_gam_winter_destabilization_caption.txt"
)
if (!file.copy(png_file, final_png, overwrite = TRUE, copy.date = TRUE) ||
    !file.copy(pdf_file, final_pdf, overwrite = TRUE, copy.date = TRUE) ||
    !file.copy(
      file.path(OUT_DIR, "leptolyngbya_gam_winter_destabilization_caption.txt"),
      final_caption, overwrite = TRUE, copy.date = TRUE
    )) {
  stop("Could not copy the GAM supplemental figure into supplemental_final.")
}

cat("\n=== Final GAM transition summary ===\n")
print(gam_summary)
cat("\nOutputs written to: ", OUT_DIR, "\n", sep = "")
