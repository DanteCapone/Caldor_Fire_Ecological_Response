# Stage 2 of the Figure 5C compound-disturbance state-space analysis.
# The frozen environmental configuration is read before Leptolyngbya data are
# loaded. Temporal windows and PCA variables are never selected using biology.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(yaml)
  library(ggrepel)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")
source("scripts/main_analysis/figure_5_state_space_helpers.R")

PROC_DIR <- file.path("data", "processed", "figure_5_state_space")
FIG_DIR <- file.path("figures", "supplemental", "figure_5_state_space")
MAIN_FIG_DIR <- "figures_v2"
CONFIG_FILE <- file.path("project_state", "FIGURE_5_STATE_SPACE_CONFIG.yml")
ENVIRONMENT_FILE <- file.path(PROC_DIR, "environmental_state_definitions.rds")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAIN_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(CONFIG_FILE) || !file.exists(ENVIRONMENT_FILE)) {
  stop("Run figure_5_state_space_environment.R before the biological stage.")
}

frozen_config <- read_yaml(CONFIG_FILE)
environmental <- readRDS(ENVIRONMENT_FILE)
if (!identical(frozen_config$biological_data_used_to_define_windows_or_pcas, FALSE)) {
  stop("Frozen configuration does not certify biology-independent definitions.")
}

CORE_DEPTHS_M <- c(5, 20, 40, 60, 75, 90)
MAX_MONTHLY_INTEGRATION_GAP_DAYS <- 62L
MINIMUM_SAMPLED_MONTHS <- 3L
RESPONSE_MONTHS <- 8:12
RESPONSE_START_DOY <- standard_doy(as.Date("2001-08-01"))
RESPONSE_END_DOY <- standard_doy(as.Date("2001-12-31"))

integrate_monthly_positive_anomaly <- function(date, positive_anomaly,
                                                max_gap_days = 62L) {
  keep <- is.finite(as.numeric(date)) & is.finite(positive_anomaly)
  date <- as.Date(date[keep])
  positive_anomaly <- positive_anomaly[keep]
  if (length(date) < 2L) {
    return(tibble(integrated_positive_anomaly_cells_d_m2 = NA_real_,
                  n_integrated_segments = 0L, n_long_gaps = 0L,
                  maximum_gap_days = NA_real_))
  }
  ord <- order(date)
  date <- date[ord]
  positive_anomaly <- positive_anomaly[ord]
  gaps <- as.numeric(diff(date))
  valid <- gaps <= max_gap_days
  areas <- gaps * (head(positive_anomaly, -1L) + tail(positive_anomaly, -1L)) / 2
  tibble(
    integrated_positive_anomaly_cells_d_m2 = if (any(valid)) sum(areas[valid]) else NA_real_,
    n_integrated_segments = sum(valid),
    n_long_gaps = sum(!valid),
    maximum_gap_days = max(gaps)
  )
}

lepto_samples <- read_csv(
  file.path("data", "processed", "phytoplankton_lepto_samples.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  filter(!low_effort_sample, depth_num %in% CORE_DEPTHS_M)

depth_availability <- lepto_samples %>%
  distinct(date, depth_num) %>%
  count(depth_num, name = "n_dates") %>%
  mutate(core_depth = depth_num %in% CORE_DEPTHS_M)
write_csv(depth_availability, file.path(PROC_DIR, "leptolyngbya_depth_availability.csv"))

date_profiles <- lepto_samples %>%
  group_by(date) %>%
  filter(all(CORE_DEPTHS_M %in% depth_num)) %>%
  summarise(
    water_column_abundance_cells_m2 = 1000 * trapz(depth_num, abundance),
    water_column_biovolume_mm3_m2 = trapz(depth_num, biovolume),
    n_depths = n_distinct(depth_num),
    minimum_depth_m = min(depth_num),
    maximum_depth_m = max(depth_num),
    .groups = "drop"
  ) %>%
  mutate(year = year(date), month = month(date), month_id = year * 12L + month)

make_response <- function(environmental_years) {
  monthly <- date_profiles %>%
    filter(month %in% RESPONSE_MONTHS) %>%
    mutate(response_year = year) %>%
    group_by(response_year, year, month, month_id) %>%
    summarise(
      monthly_date = as.Date(median(as.numeric(date)), origin = "1970-01-01"),
      monthly_abundance_cells_m2 = median(water_column_abundance_cells_m2),
      monthly_biovolume_mm3_m2 = median(water_column_biovolume_mm3_m2),
      n_sample_dates = n(),
      .groups = "drop"
    )

  climatology <- map_dfr(sort(unique(environmental_years)), function(target_year) {
    monthly %>%
      filter(response_year != target_year, response_year != 2021L) %>%
      group_by(month) %>%
      summarise(
        historical_month_median_cells_m2 = median(monthly_abundance_cells_m2),
        historical_month_p95_cells_m2 = as.numeric(
          quantile(monthly_abundance_cells_m2, 0.95, type = 8)
        ),
        n_historical_response_years = n_distinct(response_year),
        .groups = "drop"
      ) %>%
      mutate(response_year = target_year)
  })

  monthly_anomaly <- monthly %>%
    filter(response_year %in% environmental_years) %>%
    left_join(climatology, by = c("response_year", "month")) %>%
    mutate(
      seasonal_anomaly_cells_m2 = monthly_abundance_cells_m2 - historical_month_median_cells_m2,
      positive_seasonal_anomaly_cells_m2 = pmax(seasonal_anomaly_cells_m2, 0),
      exceeds_historical_p95 = monthly_abundance_cells_m2 > historical_month_p95_cells_m2
    )

  observed_min <- floor_date(min(date_profiles$date), unit = "month")
  observed_max <- ceiling_date(max(date_profiles$date), unit = "month") - days(1)
  annual <- monthly_anomaly %>%
    group_by(response_year) %>%
    group_modify(~{
      integration <- integrate_monthly_positive_anomaly(
        .x$monthly_date, .x$positive_seasonal_anomaly_cells_m2,
        MAX_MONTHLY_INTEGRATION_GAP_DAYS
      )
      tibble(
        response_start_date = date_from_standard_doy(.y$response_year, RESPONSE_START_DOY),
        response_end_date = date_from_standard_doy(.y$response_year, RESPONSE_END_DOY),
        n_sampled_months = n_distinct(.x$month_id),
        n_sample_dates = sum(.x$n_sample_dates),
        first_sample_date = min(.x$monthly_date),
        last_sample_date = max(.x$monthly_date),
        median_sample_date = as.Date(median(as.numeric(.x$monthly_date)), origin = "1970-01-01"),
        maximum_monthly_anomaly_cells_m2 = max(.x$seasonal_anomaly_cells_m2),
        maximum_monthly_abundance_cells_m2 = max(.x$monthly_abundance_cells_m2),
        fraction_sampled_months_above_historical_p95 = mean(.x$exceeds_historical_p95),
        integration
      )
    }) %>%
    ungroup() %>%
    mutate(
      complete_response_interval_observable = response_start_date >= observed_min &
        response_end_date <= observed_max,
      adequate_monthly_sampling = n_sampled_months >= MINIMUM_SAMPLED_MONTHS,
      included_in_state_space = complete_response_interval_observable & adequate_monthly_sampling,
      response_start_doy = RESPONSE_START_DOY,
      response_end_doy = RESPONSE_END_DOY
    )
  list(monthly = monthly_anomaly, annual = annual)
}

smoke_definitions <- environmental$smoke_windows %>%
  select(smoke_definition, smoke_start_doy = start_doy, smoke_end_doy = end_doy)
response_by_smoke_definition <- setNames(
  map(smoke_definitions$smoke_start_doy,
      ~make_response(environmental$primary$scores$year)),
  smoke_definitions$smoke_definition
)
primary_response <- response_by_smoke_definition[[frozen_config$primary_smoke_definition]]

primary_joint_all <- environmental$primary$annual %>%
  select(year, everything()) %>%
  left_join(environmental$primary$scores, by = "year") %>%
  left_join(primary_response$annual, by = c("year" = "response_year"))
primary_joint <- primary_joint_all %>%
  filter(is.finite(antecedent_lake_state_index), is.finite(wildfire_exposure_index)) %>%
  mutate(response_available = coalesce(included_in_state_space, FALSE))
if (!all(2011:2025 %in% primary_joint$year)) {
  stop("The primary plotted cohort does not include every year from 2011-2025.")
}

write_csv(primary_response$monthly, file.path(PROC_DIR, "primary_leptolyngbya_monthly_anomalies.csv"))
write_csv(primary_response$annual, file.path(PROC_DIR, "primary_leptolyngbya_response_metrics.csv"))
write_csv(primary_joint_all, file.path(PROC_DIR, "primary_state_space_annual_summary.csv"))
# Compatibility paths used by the existing manifest and downstream review notes.
write_csv(primary_joint, file.path("data", "processed", "historical_joint_state_annual_summary.csv"))
write_csv(
  primary_joint_all %>%
    filter(!(included_in_state_space %in% TRUE)) %>%
    transmute(
      year,
      missing_reason = case_when(
        !is.finite(antecedent_lake_state_index) | !is.finite(wildfire_exposure_index) ~
          "incomplete environmental PCA inputs",
        is.na(n_sampled_months) ~
          "no August-December biological samples; environmental point retained without response size",
        !(complete_response_interval_observable %in% TRUE) ~ "incomplete response interval",
        !(adequate_monthly_sampling %in% TRUE) ~ "fewer than three sampled months",
        TRUE ~ "not included"
      )
    ),
  file.path("data", "processed", "historical_joint_state_dropped_years.csv")
)

state_variance <- environmental$primary$state_pca$variance %>%
  filter(component == "PC1") %>% pull(variance_explained)
highlight_years <- c(2011L, 2020L, 2021L, 2022L, 2023L, 2024L, 2025L)
year_colours <- c(
  "2011" = "#009E73",
  "2020" = "#0072B2",
  "2021" = "#D55E00",
  "2022" = "#CC79A7",
  "2023" = "#E69F00",
  "2024" = "#56B4E9",
  "2025" = "#882255",
  "Other years" = "grey72"
)
plot_data <- primary_joint %>%
  mutate(
    year_group = if_else(year %in% highlight_years, as.character(year), "Other years"),
    focal_label = if_else(year %in% highlight_years, as.character(year), NA_character_),
    leptolyngbya_present = is.finite(maximum_monthly_abundance_cells_m2) &
      maximum_monthly_abundance_cells_m2 > 0
  )

p_joint_abundance <- ggplot(
  plot_data,
  aes(antecedent_lake_state_index, wildfire_exposure_index)
) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  # A cross means observed zero abundance. The one year lacking biological
  # observations remains an open environmental point and is not called zero.
  geom_point(
    data = plot_data %>% filter(!is.finite(maximum_monthly_abundance_cells_m2)),
    shape = 21, size = 2.0, fill = NA, colour = "grey60", stroke = 0.65,
    show.legend = FALSE
  ) +
  geom_point(
    data = plot_data %>%
      filter(is.finite(maximum_monthly_abundance_cells_m2), !leptolyngbya_present),
    aes(colour = year_group),
    shape = 4, size = 2.3, stroke = 0.7, show.legend = FALSE
  ) +
  geom_point(
    data = plot_data %>% filter(leptolyngbya_present),
    aes(fill = year_group, colour = year_group),
    shape = 21, size = 2.0, stroke = 0.45, show.legend = FALSE
  ) +
  geom_point(
    data = plot_data %>% filter(leptolyngbya_present),
    aes(size = maximum_monthly_abundance_cells_m2,
        fill = year_group, colour = year_group),
    shape = 21, alpha = 0.86, stroke = 0.45
  ) +
  geom_text_repel(
    data = plot_data %>% filter(!is.na(focal_label)),
    aes(label = focal_label), family = "Times New Roman", size = lo_geom_text_size(8.5),
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.35,
    seed = 2021, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = year_colours, guide = "none") +
  scale_colour_manual(values = year_colours, guide = "none") +
  scale_size_area(
    max_size = 7.5,
    limits = c(0, 5e10),
    breaks = c(1e9, 5e9, 2e10, 5e10),
    labels = c("1", "5", "20", "50"),
    name = expression("Peak Aug-Dec "*italic("Leptolyngbya")*
                        " (10"^9*" cells m"^{-2}*")"),
    guide = guide_legend(
      direction = "horizontal", nrow = 1, byrow = TRUE,
      title.position = "top", label.position = "bottom",
      override.aes = list(fill = "grey60", colour = "grey20")
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.09))) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.20)),
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    x = "Lake Stratification State (PC1)",
    y = expression("Integrated PM"[2.5]*" excess ("*mu*"g d m"^{-3}*")")
  ) +
  theme_classic(base_size = 9, base_family = "Times New Roman") +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.1),
    legend.text = element_text(size = 8.1),
    legend.key.width = grid::unit(5.5, "mm"),
    legend.key.height = grid::unit(4.5, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0, "mm"),
    plot.margin = margin(2.5, 3, 1, 2, "mm")
  )
p_joint <- p_joint_abundance

ggsave(file.path(MAIN_FIG_DIR, "historical_joint_state_abundance.png"), p_joint_abundance,
       width = 12.7, height = 10.5, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(MAIN_FIG_DIR, "historical_joint_state_abundance.pdf"), p_joint_abundance,
       width = 12.7, height = 10.5, units = "cm", device = cairo_pdf)
ggsave(file.path(MAIN_FIG_DIR, "historical_joint_state.png"), p_joint,
       width = 12.7, height = 10.5, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(MAIN_FIG_DIR, "historical_joint_state.pdf"), p_joint,
       width = 12.7, height = 10.5, units = "cm", device = cairo_pdf)

loading_data <- bind_rows(
  environmental$primary$state_pca$loadings
) %>%
  left_join(bind_rows(
    environmental$state_metadata %>% mutate(block = "Antecedent lake state"),
    environmental$wildfire_metadata %>% mutate(block = "Wildfire exposure")
  ), by = c("block", "variable")) %>%
  mutate(label = coalesce(label, variable))
p_loadings <- ggplot(loading_data, aes(reorder(label, PC1), PC1, fill = PC1 > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  coord_flip() +
  facet_wrap(~block, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00")) +
  labs(x = NULL, y = "PC1 loading") +
  theme_classic(base_size = 9, base_family = "Times New Roman") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(FIG_DIR, "primary_environmental_pca_loadings.png"), p_loadings,
       width = 12.7, height = 11, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "primary_environmental_pca_loadings.pdf"), p_loadings,
       width = 12.7, height = 11, units = "cm", device = cairo_pdf)

sensitivity_joint <- map_dfr(environmental$sensitivity_scenarios, function(scenario) {
  smoke_name <- scenario$metadata$smoke_definition
  response <- response_by_smoke_definition[[smoke_name]]$annual
  scenario$scores %>%
    left_join(response, by = c("year" = "response_year")) %>%
    filter(is.finite(antecedent_lake_state_index), is.finite(wildfire_exposure_index)) %>%
    mutate(
      smoke_definition = smoke_name,
      smoke_start_doy = scenario$metadata$smoke_start_doy,
      smoke_end_doy = scenario$metadata$smoke_end_doy,
      stability_fraction = scenario$metadata$stability_fraction,
      onset_doy = scenario$metadata$onset_doy,
      scenario_id = scenario$metadata$scenario_id,
      focal = year %in% c(2011L, 2020L, 2021L)
    )
})
write_csv(sensitivity_joint, file.path(PROC_DIR, "sensitivity_state_space_scores_and_response.csv"))
write_csv(
  sensitivity_joint %>% filter(focal) %>%
    mutate(
      x_quadrant = sign(antecedent_lake_state_index),
      y_quadrant = sign(wildfire_exposure_index)
    ),
  file.path(PROC_DIR, "sensitivity_focal_year_diagnostics.csv")
)

smoke_display <- c(
  "fixed_aug01_oct31" = "Fixed August 1-October 31",
  "p10_p90" = "10th-90th percentile",
  "p05_p95" = "5th-95th percentile",
  "p025_p975" = "2.5th-97.5th percentile"
)
for (smoke_name in intersect(names(smoke_display), unique(sensitivity_joint$smoke_definition))) {
  plot_sensitivity <- sensitivity_joint %>%
    filter(smoke_definition == smoke_name) %>%
    mutate(
      threshold_label = factor(
        if_else(is.na(stability_fraction), "Fixed July-August antecedent",
                paste0(round(100 * stability_fraction), "% onset")),
        levels = c("Fixed July-August antecedent", "40% onset", "50% onset", "60% onset")
      ),
      year_group = if_else(focal, as.character(year), "Other years")
    )
  p_sensitivity <- ggplot(
    plot_sensitivity,
    aes(antecedent_lake_state_index, wildfire_exposure_index)
  ) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
    geom_point(shape = 21, size = 1.8, fill = "grey82", colour = "grey25",
               stroke = 0.3, show.legend = FALSE) +
    geom_point(
               data = plot_sensitivity %>%
                 filter(is.finite(integrated_positive_anomaly_cells_d_m2)),
               aes(size = integrated_positive_anomaly_cells_d_m2, fill = year_group),
               shape = 21, colour = "grey20", alpha = 0.85, stroke = 0.3,
               show.legend = FALSE) +
    geom_text_repel(
      data = plot_sensitivity %>% filter(focal),
      aes(label = year), family = "Times New Roman", size = lo_geom_text_size(8.5),
      min.segment.length = 0, seed = 2021, max.overlaps = Inf,
      show.legend = FALSE
    ) +
    facet_wrap(~threshold_label, ncol = 3) +
    scale_fill_manual(values = year_colours) +
    scale_size_area(max_size = 6) +
    labs(
      x = "Antecedent Lake-State Index",
      y = "Integrated PM2.5 excess (ug d m^-3)",
      title = smoke_display[[smoke_name]]
    ) +
    theme_classic(base_size = 9, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  output_stub <- paste0("state_space_sensitivity_", smoke_name)
  ggsave(file.path(FIG_DIR, paste0(output_stub, ".png")), p_sensitivity,
         width = 12.7, height = 8.5, units = "cm", dpi = 600, device = ragg::agg_png)
  ggsave(file.path(FIG_DIR, paste0(output_stub, ".pdf")), p_sensitivity,
         width = 12.7, height = 8.5, units = "cm", device = cairo_pdf)
}

record_file <- file.path(PROC_DIR, "analysis_record.md")
cat(
  paste0(
    "\n## G. Leptolyngbya response and Figure 5C\n",
    "After reading the frozen configuration, samples were aggregated to monthly 5-90 m ",
    "water-column abundance using fixed depths (5, 20, 40, 60, 75, and 90 m) during August-December of each calendar year. ",
    "Month-specific leave-one-response-year-out medians excluded 2021 from every historical baseline. ",
    "Positive anomalies were trapezoidally integrated only across gaps <= ",
    MAX_MONTHLY_INTEGRATION_GAP_DAYS, " days. Years required at least ",
    MINIMUM_SAMPLED_MONTHS, " sampled months and an observable complete response interval. Figure 5C bubble area maps the maximum monthly August-December depth-integrated abundance. Crosses indicate observed zero abundance; a year without biological observations remains an open environmental point and is not classified as zero. Years 2011 and 2020-2025 are colored and labeled for comparison.\n",
    "The primary figure and fixed-window diagnostic were generated without reselecting dates, variables, the antecedent PCA orientation, or the PM2.5 metric from biological results.\n"
  ),
  file = record_file, append = TRUE
)

message(
  "Biological response added after configuration freeze; Figure 5C cohort n = ",
  nrow(primary_joint), "."
)
