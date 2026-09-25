# =============================================================================
# Figure 5 quantitative synthesis
#
# Manuscript-facing synthesis of (a) community and Leptolyngbya trajectories,
# (b) first observed statistical return across measured pathways, and
# (c) descriptive historical context for antecedent lake state, smoke exposure,
# and the fixed-window Leptolyngbya response. This workflow does not fit a
# biological interaction model and does not interpret observational contrasts
# as causal, necessary, or sufficient conditions.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(ragg)
  library(yaml)
})

source("scripts/main_analysis/shared_aesthetics.R")
source("scripts/main_analysis/figure_5_state_space_helpers.R")

# ---- Fixed definitions ------------------------------------------------------
FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
MONTH_DAYS <- 30.4375
ANALYSIS_YEARS <- 2006:2025
FOCAL_YEARS <- c(2011L, 2020L, 2021L)
REFERENCE_YEARS <- setdiff(ANALYSIS_YEARS, FOCAL_YEARS)
ANTECEDENT_START_MD <- "07-01"
ANTECEDENT_END_MD <- "08-13"
EXPOSURE_START_MD <- "08-14"
EXPOSURE_END_MD <- "10-21"
RESPONSE_START_MD <- "08-01"
RESPONSE_END_MD <- "12-31"
PM_THRESHOLD <- 9.1
MINIMUM_PM_COVERAGE <- 0.80
CORE_DEPTHS_M <- c(5, 20, 40, 60, 75, 90)
MINIMUM_RESPONSE_DATES <- 3L
MINIMUM_RESPONSE_SPAN_DAYS <- 60
BASE_FAMILY <- "Times New Roman"

RENDER_ONLY <- identical(Sys.getenv("CALDOR_FIGURE_5_RENDER_ONLY", unset = "0"), "1")
OUT_DIR <- Sys.getenv(
  "CALDOR_FIGURE_5_OUT_DIR",
  unset = file.path("figures", "figure_5_nmds_braycurtis")
)
MIRROR_DIR <- Sys.getenv("CALDOR_FIGURE_5_MIRROR_DIR", unset = "figures_v2")
PROC_DIR <- file.path("data", "processed", "figure_5_state_space")
STATE_FIG_DIR <- Sys.getenv(
  "CALDOR_FIGURE_5_STATE_FIG_DIR",
  unset = file.path("figures", "supplemental", "figure_5_state_space")
)
CONFIG_FILE <- file.path("project_state", "FIGURE_5_STATE_SPACE_CONFIG.yml")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MIRROR_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATE_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path("figures", "supplemental", "resistance_and_resilience",
            "pathway_resistance_resilience_markers.csv"),
  file.path("figures", "supplemental", "resistance_and_resilience",
            "pathway_resistance_resilience_timeseries.csv"),
  file.path(PROC_DIR, "mltp_cast_state_metrics.csv"),
  file.path("data", "processed", "ctd", "mltp_schmidt_stability_2005_2025.csv"),
  file.path(PROC_DIR, "pm25_resolved_daily.csv"),
  file.path("data", "processed", "tahoe_hms_pm25_daily.csv"),
  file.path("data", "processed", "phytoplankton_lepto_samples.csv"),
  file.path("data", "lake_environmental_data", "phytoplankton",
            "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing Figure 5 inputs: ", paste(missing_files, collapse = ", "))
}

# Retain the pre-revision figure products exactly once. Existing SIMPER and
# state-space products are never removed or renamed by this script.
archive_dir <- file.path(OUT_DIR, "archived", "pre_quantitative_synthesis_2026-09-16")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
archive_targets <- c(
  file.path(OUT_DIR, "figure_5_final.png"),
  file.path(OUT_DIR, "figure_5_final.pdf"),
  file.path(MIRROR_DIR, "figure_5_final.png"),
  file.path(MIRROR_DIR, "figure_5_final.pdf")
)
for (source_file in archive_targets[file.exists(archive_targets)]) {
  archive_name <- paste0(
    if_else(dirname(source_file) == MIRROR_DIR, "figures_v2__", "figures__"),
    basename(source_file)
  )
  archive_file <- file.path(archive_dir, archive_name)
  if (!file.exists(archive_file)) file.copy(source_file, archive_file)
}

window_start <- function(year, month_day) as.Date(paste0(year, "-", month_day))
window_end <- window_start
safe_min <- function(x) if (all(!is.finite(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)
scientific_axis_labels <- function(x) {
  vapply(x, function(value) {
    if (!is.finite(value)) return("")
    if (value == 0) return("0")
    exponent <- floor(log10(abs(value)))
    coefficient <- signif(value / (10^exponent), 2)
    coefficient_label <- format(coefficient, scientific = FALSE, trim = TRUE)
    exponent_label <- chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(exponent))
    paste0(coefficient_label, "×10", exponent_label)
  }, character(1))
}

figure_theme <- theme_bw(base_size = 8.5, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.22),
    axis.text = element_text(size = 8.5, colour = "grey20"),
    axis.title = element_text(size = 8.5, face = "bold"),
    strip.text = element_text(size = 8.5, face = "bold"),
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.tag = element_text(size = 10, face = "bold", hjust = 1, vjust = 1),
    plot.tag.position = c(0.99, 0.99),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(2.4, "mm"),
    legend.key.width = unit(4.0, "mm"),
    plot.margin = margin(1.2, 1.8, 1.2, 1.8, "mm")
  )

# ---- Panels A and B: shared recovery products ------------------------------
marker_file <- required_files[1]
series_file <- required_files[2]
markers <- read_csv(marker_file, show_col_types = FALSE) %>%
  mutate(
    peak_date = as.Date(peak_date),
    recovery_date = as.Date(recovery_date),
    recovery_date_moving_mean_3 = as.Date(recovery_date_moving_mean_3),
    recovery_previous_date = as.Date(recovery_previous_date)
  )
series <- read_csv(series_file, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

recovery_rule_for <- function(variable) {
  case_when(
    str_detect(variable, "^Atmospheric deposition (NO3|NH4)$") ~
      "first post-peak observation within the historical reference envelope",
    variable == "Chlorophyll-a sum (60-105 m)" ~
      "first trailing three-observation mean within the corresponding mean reference envelope",
    TRUE ~
      "first of two consecutive post-peak observations within the historical reference envelope"
  )
}

confirmation_for <- function(pathway, variable, recovery_date, rule) {
  if (is.na(recovery_date)) return(as.Date(NA))
  # A one-observation exception has no separate confirmation observation.
  if (str_detect(rule, "first post-peak")) return(as.Date(NA))
  # For the trailing-mean exception, the recovery date is the third and final
  # observation that makes the three-observation summary meet the criterion.
  if (str_detect(rule, "trailing three-observation")) return(recovery_date)
  candidate <- series %>%
    filter(Pathway == pathway, Variable == variable, date >= recovery_date,
           within_envelope %in% TRUE) %>%
    arrange(date) %>%
    slice_head(n = 2)
  if (nrow(candidate) < 2L || candidate$date[1] != recovery_date) as.Date(NA) else
    candidate$date[2]
}

recovery_export <- markers %>%
  mutate(
    classification = case_when(
      recovery_status %in% c("recovered", "not_recovered") ~ "affected",
      preexisting_departure %in% TRUE & !(fire_amplified %in% TRUE) ~ "pre-existing",
      recovery_status == "unaffected" ~ "unaffected",
      TRUE ~ "not estimable"
    ),
    recovery_rule = map_chr(Variable, recovery_rule_for),
    confirmation_date = as.Date(
      pmap_dbl(
        list(Pathway, Variable, recovery_date, recovery_rule),
        ~as.numeric(confirmation_for(..1, ..2, ..3, ..4))
      ),
      origin = "1970-01-01"
    ),
    sensitivity_1obs = recovery_months_1,
    sensitivity_2obs = recovery_months_2,
    sensitivity_3obs = recovery_months_3,
    sensitivity_min = pmap_dbl(
      list(sensitivity_1obs, sensitivity_2obs, sensitivity_3obs),
      ~safe_min(c(...))
    ),
    sensitivity_max = pmap_dbl(
      list(sensitivity_1obs, sensitivity_2obs, sensitivity_3obs),
      ~safe_max(c(...))
    ),
    selected_recovery_months = recovery_months,
    variable_specific_exception = str_detect(
      recovery_rule, "first post-peak|trailing three-observation"
    ),
    classification_note = case_when(
      preexisting_departure %in% TRUE & fire_amplified %in% TRUE ~
        "pre-existing departure retained because the post-fire peak met the current amplification rule",
      preexisting_departure %in% TRUE ~ "pre-existing departure; excluded from numerical ranking",
      recovery_status == "unaffected" ~ "unaffected under the current workflow",
      TRUE ~ NA_character_
    )
  ) %>%
  transmute(
    pathway = as.character(Pathway), variable = Variable, classification,
    selected_recovery_months, recovery_rule,
    sensitivity_1obs, sensitivity_2obs, sensitivity_3obs,
    sensitivity_min, sensitivity_max, peak_date, recovery_date,
    confirmation_date, recovery_previous_date, recovery_gap_days,
    recovery_lower_months, peak_anomaly_sd, peak_direction,
    return_direction = case_when(
      peak_anomaly_sd < 0 ~ "Negative anomaly returning upward",
      peak_anomaly_sd > 0 ~ "Positive anomaly returning downward",
      TRUE ~ "No directional anomaly"
    ),
    variable_specific_exception, classification_note, recovery_status,
    follow_up_months
  )

write_csv(recovery_export, file.path(OUT_DIR, "figure_5_recovery_synthesis.csv"))

trajectory_lookup <- tribble(
  ~metric, ~variable_pattern, ~stratum,
  "Community departure", "Community composition (Bray-Curtis; 0-40 m)", "0-40 m",
  "Community departure", "Community composition (Bray-Curtis; 60-105 m)", "60-105 m"
)

panel_a_start <- as.Date("2019-01-01")
panel_a_end <- as.Date("2025-12-31")
community_sampling_coverage <- read_csv(required_files[8], show_col_types = FALSE) %>%
  transmute(
    date = as.Date(date), stratum = depth_bin, depth_m = as.numeric(depth_num),
    taxon, abundance = as.numeric(abundance)
  ) %>%
  filter(stratum %in% c("0-40 m", "60-105 m")) %>%
  group_by(date, stratum) %>%
  summarise(
    sampled_depth_count = n_distinct(depth_m),
    sampled_depths_m = paste(sort(unique(depth_m)), collapse = ";"),
    positive_taxa_count = n_distinct(taxon[abundance > 0]),
    total_abundance_cells_l = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

panel_a_data_all <- series %>%
  inner_join(trajectory_lookup, by = c("Variable" = "variable_pattern")) %>%
  filter(date >= panel_a_start, date <= panel_a_end) %>%
  left_join(
    markers %>%
      inner_join(trajectory_lookup, by = c("Variable" = "variable_pattern")) %>%
      select(metric, stratum, peak_date, recovery_date, recovery_months),
    by = c("metric", "stratum")
  ) %>%
  left_join(community_sampling_coverage, by = c("date", "stratum")) %>%
  mutate(
    is_peak = date == peak_date,
    is_recovery = date == recovery_date,
    marker = case_when(is_peak ~ "Peak", is_recovery ~ "First observed return", TRUE ~ NA_character_),
    recovery_rule = map_chr(Variable, recovery_rule_for),
    included_minimum_3_depths = sampled_depth_count >= 3,
    exclusion_reason = if_else(
      included_minimum_3_depths,
      NA_character_,
      "fewer than three sampled depths in the depth zone"
    )
  ) %>%
  transmute(
    metric, stratum, variable = Variable, date, value, expected,
    reference_lo, reference_hi, n_historical = n_hist,
    reference_method, sampled_depth_count, sampled_depths_m,
    positive_taxa_count, total_abundance_cells_l,
    peak_date, recovery_date, recovery_months, recovery_rule, is_peak, is_recovery,
    included_minimum_3_depths, exclusion_reason
  )
panel_a_data <- panel_a_data_all %>% filter(included_minimum_3_depths)
if (any(!panel_a_data_all$included_minimum_3_depths &
        (panel_a_data_all$is_peak | panel_a_data_all$is_recovery))) {
  stop("A selected Panel A peak or recovery visit has fewer than three sampled depths.")
}
if (any(panel_a_data$sampled_depth_count < 3, na.rm = TRUE)) {
  stop("Panel A contains a visit with fewer than three sampled depths in its zone.")
}
write_csv(panel_a_data, file.path(OUT_DIR, "figure_5_panel_a_data.csv"))
write_csv(
  panel_a_data_all %>%
    select(date, stratum, value, expected, reference_lo, reference_hi,
           sampled_depth_count, sampled_depths_m, positive_taxa_count,
           total_abundance_cells_l, included_minimum_3_depths, exclusion_reason),
  file.path(OUT_DIR, "figure_5_panel_a_sampling_diagnostic.csv")
)

stratum_colours <- c("0-40 m" = "#7BC8A4", "60-105 m" = "#7FB3D5")
recovery_lines <- panel_a_data %>%
  filter(is_recovery) %>%
  distinct(recovery_date)

p_community <- ggplot(
  filter(panel_a_data, metric == "Community departure"),
  aes(date, value, colour = stratum, fill = stratum)
) +
  annotate("rect", xmin = FIRE_START, xmax = FIRE_END,
           ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
  geom_ribbon(aes(ymin = reference_lo, ymax = reference_hi, group = stratum),
              alpha = 0.11, colour = NA) +
  geom_line(aes(y = reference_lo, group = stratum), linewidth = 0.30, alpha = 0.52) +
  geom_line(aes(y = reference_hi, group = stratum), linewidth = 0.30, alpha = 0.52) +
  geom_line(aes(y = expected, group = stratum), linewidth = 0.28,
            linetype = "dashed", alpha = 0.8) +
  geom_line(aes(group = stratum), linewidth = 0.55) +
  geom_vline(
    xintercept = FIRE_START, colour = "firebrick", linetype = "dashed",
    linewidth = 0.45, alpha = 0.85
  ) +
  geom_vline(
    data = recovery_lines, aes(xintercept = recovery_date),
    colour = "black", linetype = "dashed", linewidth = 0.45
  ) +
  scale_colour_manual(values = stratum_colours, name = NULL) +
  scale_fill_manual(values = stratum_colours, guide = "none") +
  scale_x_date(
    limits = c(panel_a_start, panel_a_end),
    breaks = as.Date(paste0(2019:2025, "-07-01")),
    date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) +
  labs(
    tag = "a", title = NULL,
    x = NULL,
    y = "Bray-Curtis\ncommunity departure"
  ) +
  figure_theme +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.box = "horizontal", legend.margin = margin(-1, 0, 0, 0, "mm"),
    axis.text.x = element_text(size = 7.2, angle = 45, hjust = 1),
    plot.margin = margin(1.2, 3.5, 1.2, 1.2, "mm")
  )
panel_a <- p_community

display_group_levels <- c(
    "Radiative", "Nutrient", "Runoff", "Biological"
)
short_variable <- function(x) {
  x %>%
    str_replace("Community composition \\(Bray-Curtis; ", "Community departure (") %>%
    str_replace("Leptolyngbya abundance \\(", "Leptolyngbya (") %>%
    str_replace("Phytoplankton abundance \\(", "Phyto. abundance (") %>%
    str_replace("Phytoplankton biovolume \\(", "Phyto. biovolume (") %>%
    str_replace("Chlorophyll-a sum \\(", "Chl-a sum (") %>%
    str_replace("Atmospheric deposition DIN:SRP molar ratio", "DIN:SRP molar ratio") %>%
    str_replace("Atmospheric deposition ", "") %>%
    str_replace("Runoff discharge \\(Upper Truckee \\+ Blackwood\\)", "Runoff discharge") %>%
    str_replace("In-lake 0-10 m TRP", "TRP (in-lake 0-10 m)") %>%
    str_replace("Secchi depth \\(m\\)", "Secchi depth") %>%
    str_replace("PAR 1% depth \\(m\\)", "PAR 1% depth") %>%
    str_replace("PM2.5 \\(µg m⁻³\\)", "PM2.5") %>%
    str_replace("Surface temperature \\(°C\\)", "Surface temperature")
}
italicize_leptolyngbya_axis <- function(labels) {
  lapply(labels, function(label) {
    if (str_starts(label, "Leptolyngbya")) {
      suffix <- str_remove(label, "^Leptolyngbya")
      bquote(italic(Leptolyngbya) * .(suffix))
    } else {
      label
    }
  })
}

recovery_plot_data <- recovery_export %>%
  filter(classification == "affected", is.finite(selected_recovery_months)) %>%
  mutate(
    display_group = case_when(
      pathway == "Radiative" ~ "Radiative",
      pathway == "Nutrient" ~ "Nutrient",
      pathway == "Runoff" ~ "Runoff",
      pathway == "Biological" ~ "Biological"
    ),
    display_group = factor(display_group, levels = display_group_levels),
    row_label = short_variable(variable),
    sensitivity_n = rowSums(is.finite(as.matrix(
      pick(sensitivity_1obs, sensitivity_2obs, sensitivity_3obs)
    ))),
    show_sensitivity_interval = sensitivity_n >= 2L,
    return_direction = factor(
      return_direction,
      levels = c(
        "Negative anomaly returning upward",
        "Positive anomaly returning downward",
        "No directional anomaly"
      )
    )
  ) %>%
  arrange(display_group, selected_recovery_months, row_label) %>%
  mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

fire_end_months <- as.numeric(FIRE_END - FIRE_START) / MONTH_DAYS
recovery_colours <- c(
  "Radiative" = "#0072B2",
  "Nutrient" = "#CC79A7",
  "Runoff" = "#E69F00",
  "Biological" = "#009E73"
)
p_recovery <- ggplot(
  recovery_plot_data,
  aes(selected_recovery_months, row_label, colour = display_group)
) +
  geom_vline(xintercept = fire_end_months, colour = "firebrick",
             linetype = "dashed", linewidth = 0.42) +
  geom_segment(
    data = filter(recovery_plot_data, show_sensitivity_interval),
    aes(x = sensitivity_min, xend = sensitivity_max,
        y = row_label, yend = row_label),
    linewidth = 0.42, alpha = 0.75
  ) +
  geom_point(aes(shape = return_direction), size = 2.25, fill = "white", stroke = 0.75) +
  ggh4x::facet_grid2(
    display_group ~ ., scales = "free_y", space = "free_y", switch = "y",
    strip = ggh4x::strip_themed(
      text_y = ggh4x::elem_list_text(
        colour = unname(recovery_colours[display_group_levels]),
        face = "bold", size = 8.5, angle = 0, hjust = 1
      )
    )
  ) +
  scale_colour_manual(values = recovery_colours, guide = "none") +
  scale_shape_manual(
    values = c(
      "Negative anomaly returning upward" = 24,
      "Positive anomaly returning downward" = 25,
      "No directional anomaly" = 21
    ),
    labels = c(
      "Negative anomaly returning upward" = "Negative ↑",
      "Positive anomaly returning downward" = "Positive ↓",
      "No directional anomaly" = "No directional anomaly"
    ),
    name = NULL,
    guide = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(size = 1.65, stroke = 0.60)
    )
  ) +
  scale_x_continuous(
    limits = c(0, max(recovery_plot_data$sensitivity_max, na.rm = TRUE) + 1),
    breaks = seq(0, 24, by = 3), expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_discrete(labels = italicize_leptolyngbya_axis) +
  labs(
    tag = "c", title = "Ecosystem Recovery",
    x = "Months since Caldor Fire",
    y = NULL
  ) +
  figure_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.3, "mm"),
    panel.grid.major.y = element_blank(),
    legend.position = "inside", legend.position.inside = c(0.985, 0.985),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.86), colour = NA),
    legend.text = element_text(size = 6.8),
    legend.key.width = unit(3.0, "mm"),
    legend.spacing.x = unit(0.6, "mm"),
    legend.margin = margin(0.5, 0.8, 0.5, 0.8, "mm"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(1, 2, 1, 9, "mm")
  )

# ---- Panel B: fixed-window environmental state and biological response ------
cast_metrics <- read_csv(file.path(PROC_DIR, "mltp_cast_state_metrics.csv"),
                         show_col_types = FALSE) %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  filter(
    year %in% ANALYSIS_YEARS,
    date >= window_start(year, ANTECEDENT_START_MD),
    date <= window_end(year, ANTECEDENT_END_MD)
  )
if (any(cast_metrics$date > window_end(cast_metrics$year, ANTECEDENT_END_MD))) {
  stop("Antecedent-state input contains an observation after August 13.")
}

antecedent_annual <- cast_metrics %>%
  group_by(year) %>%
  summarise(
    mixed_layer_thermal_contrast_c = median(delta_t_mixed_layer_c, na.rm = TRUE),
    mixed_layer_depth_m = median(mixed_layer_depth_m, na.rm = TRUE),
    mixed_layer_shallowness_m = -1 * mixed_layer_depth_m,
    surface_temperature_0_10m_c = median(surface_temperature_0_10m_c, na.rm = TRUE),
    antecedent_n_casts = n(),
    antecedent_first_date = min(date), antecedent_last_date = max(date),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~if_else(is.nan(.x), NA_real_, .x)))

# Schmidt stability is summarized over the identical July 1-August 13 window
# for external validation of the thermal-state PC1. It is not a PCA input, so
# incomplete full-depth coverage cannot exclude otherwise usable years.
stability_window_annual <- read_csv(
  file.path("data", "processed", "ctd", "mltp_schmidt_stability_2005_2025.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  filter(
    deep_enough %in% TRUE, is.finite(schmidt_stability_j_m2),
    year %in% ANALYSIS_YEARS,
    date >= window_start(year, ANTECEDENT_START_MD),
    date <= window_end(year, ANTECEDENT_END_MD)
  ) %>%
  group_by(year) %>%
  summarise(
    schmidt_stability_j_m2 = median(schmidt_stability_j_m2),
    schmidt_n_casts = n(),
    schmidt_first_date = min(date), schmidt_last_date = max(date),
    .groups = "drop"
  )

antecedent_annual <- antecedent_annual %>%
  left_join(stability_window_annual, by = "year")
focal_antecedent_dates <- antecedent_annual %>% filter(year %in% FOCAL_YEARS)
if (
  nrow(focal_antecedent_dates) != length(FOCAL_YEARS) ||
    any(
      focal_antecedent_dates$antecedent_first_date <
        window_start(focal_antecedent_dates$year, ANTECEDENT_START_MD) |
        focal_antecedent_dates$antecedent_last_date >
        window_end(focal_antecedent_dates$year, ANTECEDENT_END_MD)
    )
) {
  stop("A focal-year antecedent summary includes a date outside July 1-August 13.")
}

state_variables <- c(
  "mixed_layer_thermal_contrast_c",
  "mixed_layer_shallowness_m",
  "surface_temperature_0_10m_c"
)
pca_input <- tibble(year = ANALYSIS_YEARS) %>%
  left_join(antecedent_annual, by = "year")
pca_reference <- pca_input %>%
  filter(year %in% REFERENCE_YEARS, if_all(all_of(state_variables), is.finite))
if (any(pca_reference$year %in% FOCAL_YEARS)) {
  stop("A focal year accidentally entered the reference-only PCA fit.")
}
complete_antecedent_years <- pca_input %>%
  filter(if_all(all_of(state_variables), is.finite)) %>%
  pull(year)
if (!all(FOCAL_YEARS %in% complete_antecedent_years)) {
  stop("At least one focal year lacks complete antecedent-state data.")
}
if (nrow(pca_reference) < 8L) stop("Insufficient background years for reference PCA.")

reference_fit <- prcomp(
  pca_reference %>% select(all_of(state_variables)), center = TRUE, scale. = TRUE
)
if (sum(reference_fit$rotation[, "PC1"]) < 0) {
  reference_fit$rotation[, "PC1"] <- -reference_fit$rotation[, "PC1"]
  reference_fit$x[, "PC1"] <- -reference_fit$x[, "PC1"]
}

project_pca <- function(data, fit) {
  values <- as.matrix(data %>% select(all_of(state_variables)))
  standardized <- sweep(values, 2, fit$center, "-")
  standardized <- sweep(standardized, 2, fit$scale, "/")
  scores <- standardized %*% fit$rotation
  bind_cols(data %>% select(year), as_tibble(scores))
}

pca_complete <- pca_input %>% filter(if_all(all_of(state_variables), is.finite))
reference_scores <- project_pca(pca_complete, reference_fit) %>%
  rename(
    antecedent_pc1 = PC1, antecedent_pc2 = PC2,
    antecedent_pc3 = PC3
  )
pca_variance <- reference_fit$sdev^2 / sum(reference_fit$sdev^2)

loo_diagnostics <- map_dfr(pca_reference$year, function(omitted_year) {
  refit_data <- pca_reference %>% filter(year != omitted_year)
  refit <- prcomp(refit_data %>% select(all_of(state_variables)),
                  center = TRUE, scale. = TRUE)
  if (sum(refit$rotation[, "PC1"] * reference_fit$rotation[, "PC1"]) < 0) {
    refit$rotation[, "PC1"] <- -refit$rotation[, "PC1"]
    refit$x[, "PC1"] <- -refit$x[, "PC1"]
  }
  refit_scores <- project_pca(refit_data, refit)
  full_scores <- reference_scores %>% filter(year %in% refit_data$year)
  tibble(
    omitted_year,
    pc1_loading_cosine_similarity = sum(
      refit$rotation[, "PC1"] * reference_fit$rotation[, "PC1"]
    ) / sqrt(sum(refit$rotation[, "PC1"]^2) *
               sum(reference_fit$rotation[, "PC1"]^2)),
    overlapping_score_correlation = cor(
      refit_scores$PC1[match(full_scores$year, refit_scores$year)],
      full_scores$antecedent_pc1
    ),
    pc1_variance_explained = refit$sdev[1]^2 / sum(refit$sdev^2)
  )
})

all_years_fit <- prcomp(
  pca_complete %>% select(all_of(state_variables)), center = TRUE, scale. = TRUE
)
if (sum(all_years_fit$rotation[, "PC1"] * reference_fit$rotation[, "PC1"]) < 0) {
  all_years_fit$rotation[, "PC1"] <- -all_years_fit$rotation[, "PC1"]
  all_years_fit$x[, "PC1"] <- -all_years_fit$x[, "PC1"]
}
all_years_scores <- project_pca(pca_complete, all_years_fit)
reference_vs_all_score_correlation <- cor(
  reference_scores$antecedent_pc1[
    match(all_years_scores$year, reference_scores$year)
  ],
  all_years_scores$PC1
)
reference_vs_all_loading_cosine <- sum(
  reference_fit$rotation[, "PC1"] * all_years_fit$rotation[, "PC1"]
) / sqrt(sum(reference_fit$rotation[, "PC1"]^2) *
           sum(all_years_fit$rotation[, "PC1"]^2))

pca_parameters <- tibble(
  variable = state_variables,
  center_reference_years = unname(reference_fit$center[state_variables]),
  scale_reference_years = unname(reference_fit$scale[state_variables]),
  pc1_loading = unname(reference_fit$rotation[state_variables, "PC1"]),
  pc2_loading = unname(reference_fit$rotation[state_variables, "PC2"]),
  pc3_loading = unname(reference_fit$rotation[state_variables, "PC3"]),
  pc1_variance_explained = pca_variance[1],
  pc2_variance_explained = pca_variance[2],
  pc3_variance_explained = pca_variance[3],
  fit_years = paste(pca_reference$year, collapse = ";"),
  excluded_focal_years = paste(FOCAL_YEARS, collapse = ";")
)
write_csv(pca_parameters, file.path(OUT_DIR, "figure_5_state_space_pca_parameters.csv"))

# Schmidt stability is an external validation variable and is evaluated only
# where a qualifying full-depth cast exists in the fixed antecedent window.
stability_annual <- pca_complete %>%
  select(year, schmidt_stability_j_m2, schmidt_n_casts) %>%
  filter(is.finite(schmidt_stability_j_m2)) %>%
  inner_join(reference_scores %>% select(year, antecedent_pc1), by = "year")
if (nrow(stability_annual) < 3L) {
  stop("Insufficient qualifying years for external Schmidt-stability validation.")
}
schmidt_pearson <- cor.test(
  stability_annual$antecedent_pc1, stability_annual$schmidt_stability_j_m2,
  method = "pearson"
)
schmidt_spearman <- suppressWarnings(cor.test(
  stability_annual$antecedent_pc1, stability_annual$schmidt_stability_j_m2,
  method = "spearman", exact = FALSE
))
schmidt_summary <- tibble(
  n = nrow(stability_annual),
  pearson_correlation = unname(schmidt_pearson$estimate),
  pearson_p_value = schmidt_pearson$p.value,
  spearman_correlation = unname(schmidt_spearman$estimate),
  spearman_p_value = schmidt_spearman$p.value,
  independent_validation = TRUE,
  interpretation = "external validation; Schmidt stability is not included in the PCA"
)
write_csv(stability_annual, file.path(OUT_DIR, "figure_5_state_space_schmidt_validation_data.csv"))
write_csv(schmidt_summary, file.path(OUT_DIR, "figure_5_state_space_schmidt_validation_summary.csv"))
write_csv(stability_annual, file.path(OUT_DIR, "figure_5_state_space_schmidt_pc1_association_data.csv"))
write_csv(schmidt_summary, file.path(OUT_DIR, "figure_5_state_space_schmidt_pc1_association_summary.csv"))

p_schmidt <- ggplot(stability_annual, aes(antecedent_pc1, schmidt_stability_j_m2)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey45", fill = "grey80",
              linewidth = 0.45) +
  geom_point(size = 2, fill = "white", shape = 21) +
  labs(
    x = "Antecedent thermal-state PC1",
    y = expression("Schmidt stability (J m"^{-2}*")")
  ) +
  figure_theme
ggsave(file.path(STATE_FIG_DIR, "antecedent_pc1_schmidt_validation.png"), p_schmidt,
       width = 8.5, height = 7.0, units = "cm", dpi = 1200,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(STATE_FIG_DIR, "antecedent_pc1_schmidt_validation.pdf"), p_schmidt,
       width = 8.5, height = 7.0, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(STATE_FIG_DIR, "antecedent_pc1_schmidt_input_association.png"), p_schmidt,
       width = 8.5, height = 7.0, units = "cm", dpi = 1200,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(STATE_FIG_DIR, "antecedent_pc1_schmidt_input_association.pdf"), p_schmidt,
       width = 8.5, height = 7.0, units = "cm", device = cairo_pdf, bg = "white")

# Event-matched smoke exposure and two operational-threshold sensitivities.
resolved_pm <- read_csv(file.path(PROC_DIR, "pm25_resolved_daily.csv"),
                        show_col_types = FALSE) %>%
  transmute(date = as.Date(date), pm25 = as.numeric(pm25_max), pm25_source)
hms_daily_all <- read_csv(file.path("data", "processed", "tahoe_hms_pm25_daily.csv"),
                          show_col_types = FALSE) %>%
  transmute(date = as.Date(date), hms_smoke_day = coalesce(as.logical(smoke_day), FALSE)) %>%
  left_join(resolved_pm, by = "date") %>%
  mutate(year = year(date), month_day = format(date, "%m-%d")) %>%
  filter(year %in% ANALYSIS_YEARS)
hms_daily <- hms_daily_all %>%
  filter(
    date >= window_start(year, EXPOSURE_START_MD),
    date <= window_end(year, EXPOSURE_END_MD)
  )
if (any(
  hms_daily$date < window_start(hms_daily$year, EXPOSURE_START_MD) |
    hms_daily$date > window_end(hms_daily$year, EXPOSURE_END_MD)
)) {
  stop("Wildfire exposure includes dates outside August 14-October 21.")
}

calendar_pm_baseline <- hms_daily %>%
  filter(year %in% REFERENCE_YEARS, is.finite(pm25)) %>%
  group_by(month_day) %>%
  summarise(historical_calendar_pm25_median = median(pm25),
            n_reference_years = n_distinct(year), .groups = "drop")
hms_daily <- hms_daily %>% left_join(calendar_pm_baseline, by = "month_day")
expected_exposure_days <- as.integer(as.Date("2001-10-21") - as.Date("2001-08-14")) + 1L
exposure_annual <- hms_daily %>%
  group_by(year) %>%
  summarise(
    exposure_window_first_date = min(date),
    exposure_window_last_date = max(date),
    exposure_window_days = n_distinct(date),
    pm25_days_available = sum(is.finite(pm25)),
    pm25_daily_coverage = pm25_days_available / expected_exposure_days,
    hms_smoke_days = sum(hms_smoke_day & is.finite(pm25)),
    qualifying_smoke_days = sum(
      hms_smoke_day & is.finite(pm25) & pm25 >= PM_THRESHOLD
    ),
    integrated_pm25_excess_ug_d_m3_raw = sum(
      if_else(hms_smoke_day & is.finite(pm25) & pm25 >= PM_THRESHOLD,
              pm25 - PM_THRESHOLD, 0), na.rm = TRUE
    ),
    cumulative_pm25_hms_smoke_days_ug_d_m3_raw = sum(
      if_else(hms_smoke_day & is.finite(pm25), pm25, 0), na.rm = TRUE
    ),
    cumulative_pm25_excess_calendar_baseline_ug_d_m3_raw = sum(
      if_else(
        hms_smoke_day & is.finite(pm25) & is.finite(historical_calendar_pm25_median),
        pmax(pm25 - historical_calendar_pm25_median, 0), 0
      ), na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    integrated_pm25_excess_ug_d_m3 = if_else(
      pm25_daily_coverage >= MINIMUM_PM_COVERAGE,
      integrated_pm25_excess_ug_d_m3_raw, NA_real_
    ),
    cumulative_pm25_hms_smoke_days_ug_d_m3 = if_else(
      pm25_daily_coverage >= MINIMUM_PM_COVERAGE,
      cumulative_pm25_hms_smoke_days_ug_d_m3_raw, NA_real_
    ),
    cumulative_pm25_excess_calendar_baseline_ug_d_m3 = if_else(
      pm25_daily_coverage >= MINIMUM_PM_COVERAGE,
      cumulative_pm25_excess_calendar_baseline_ug_d_m3_raw, NA_real_
    )
  )
if (!all(exposure_annual$exposure_window_days == expected_exposure_days)) {
  stop("At least one wildfire-exposure year does not contain the exact event-matched date window.")
}
if (any(
  exposure_annual$exposure_window_first_date !=
    window_start(exposure_annual$year, EXPOSURE_START_MD) |
    exposure_annual$exposure_window_last_date !=
      window_end(exposure_annual$year, EXPOSURE_END_MD)
)) {
  stop("An annual wildfire-exposure summary includes a date outside August 14-October 21.")
}

# Panel B integrates PM2.5 over the identical August 1-December 31 window used
# for the biological AUC. Values are summed on threshold-qualified HMS smoke
# days without subtracting the operational 9.1 ug m^-3 threshold.
expected_panel_b_pm_days <- as.integer(
  as.Date("2001-12-31") - as.Date("2001-08-01")
) + 1L
panel_b_pm_daily <- hms_daily_all %>%
  filter(
    date >= window_start(year, RESPONSE_START_MD),
    date <= window_end(year, RESPONSE_END_MD)
  )
if (any(
  panel_b_pm_daily$date < window_start(panel_b_pm_daily$year, RESPONSE_START_MD) |
    panel_b_pm_daily$date > window_end(panel_b_pm_daily$year, RESPONSE_END_MD)
)) {
  stop("Panel B PM2.5 includes dates outside August 1-December 31.")
}
panel_b_pm_annual <- panel_b_pm_daily %>%
  group_by(year) %>%
  summarise(
    panel_b_pm_window_first_date = min(date),
    panel_b_pm_window_last_date = max(date),
    panel_b_pm_window_days = n_distinct(date),
    panel_b_pm_days_available = sum(is.finite(pm25)),
    panel_b_pm_daily_coverage = panel_b_pm_days_available / expected_panel_b_pm_days,
    panel_b_qualifying_smoke_days = sum(
      hms_smoke_day & is.finite(pm25) & pm25 >= PM_THRESHOLD
    ),
    integrated_pm25_aug_dec_ug_d_m3_raw = sum(
      if_else(
        hms_smoke_day & is.finite(pm25) & pm25 >= PM_THRESHOLD,
        pm25, 0
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    integrated_pm25_aug_dec_ug_d_m3 = if_else(
      panel_b_pm_daily_coverage >= MINIMUM_PM_COVERAGE,
      integrated_pm25_aug_dec_ug_d_m3_raw,
      NA_real_
    )
  )
if (!all(panel_b_pm_annual$panel_b_pm_window_days == expected_panel_b_pm_days)) {
  stop("At least one Panel B PM2.5 year lacks the exact August 1-December 31 window.")
}

# Date-level depth integration followed by temporal trapezoidal AUC.
lepto_profiles <- read_csv(
  file.path("data", "processed", "phytoplankton_lepto_samples.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date), year = year(date), depth_num = as.numeric(depth_num)) %>%
  filter(
    !low_effort_sample, depth_num %in% CORE_DEPTHS_M,
    date >= window_start(year, RESPONSE_START_MD),
    date <= window_end(year, RESPONSE_END_MD)
  ) %>%
  group_by(date, year, depth_num) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(date, year) %>%
  filter(all(CORE_DEPTHS_M %in% depth_num)) %>%
  summarise(
    depth_integrated_leptolyngbya_cells_m2 = 1000 * trapz(depth_num, abundance),
    n_depths = n_distinct(depth_num), .groups = "drop"
  )

response_annual_observed <- lepto_profiles %>%
  arrange(year, date) %>%
  group_by(year) %>%
  summarise(
    biological_n_dates = n_distinct(date),
    biological_first_date = min(date), biological_last_date = max(date),
    biological_span_days = as.numeric(biological_last_date - biological_first_date),
    biological_temporal_auc_cells_d_m2_raw = trapz(
      as.numeric(date), depth_integrated_leptolyngbya_cells_m2
    ),
    biological_maximum_cells_m2_raw = max(depth_integrated_leptolyngbya_cells_m2),
    biological_median_positive_cells_m2_raw = if (
      any(depth_integrated_leptolyngbya_cells_m2 > 0)
    ) median(depth_integrated_leptolyngbya_cells_m2[
      depth_integrated_leptolyngbya_cells_m2 > 0
    ]) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    biological_coverage_adequate = biological_n_dates >= MINIMUM_RESPONSE_DATES &
      biological_span_days >= MINIMUM_RESPONSE_SPAN_DAYS,
    biological_observed_zero = biological_coverage_adequate &
      biological_maximum_cells_m2_raw == 0,
    biological_temporal_auc_cells_d_m2 = if_else(
      biological_coverage_adequate, biological_temporal_auc_cells_d_m2_raw, NA_real_
    ),
    biological_maximum_cells_m2 = if_else(
      biological_coverage_adequate, biological_maximum_cells_m2_raw, NA_real_
    ),
    biological_median_positive_cells_m2 = if_else(
      biological_coverage_adequate,
      biological_median_positive_cells_m2_raw, NA_real_
    )
  )

state_space_annual <- tibble(year = ANALYSIS_YEARS) %>%
  left_join(pca_input, by = "year") %>%
  left_join(reference_scores, by = "year") %>%
  left_join(exposure_annual, by = "year") %>%
  left_join(panel_b_pm_annual, by = "year") %>%
  left_join(response_annual_observed, by = "year") %>%
  mutate(
    environmental_point_available = is.finite(antecedent_pc1) &
      is.finite(integrated_pm25_aug_dec_ug_d_m3),
    biological_response_available = biological_coverage_adequate %in% TRUE,
    biological_observed_zero = coalesce(biological_observed_zero, FALSE)
  )
if (any(
  !state_space_annual$biological_response_available &
    state_space_annual$biological_temporal_auc_cells_d_m2 == 0,
  na.rm = TRUE
)) {
  stop("A missing biological year was incorrectly classified as zero.")
}

write_csv(state_space_annual,
          file.path(OUT_DIR, "figure_5_state_space_annual_summary.csv"))
write_csv(
  state_space_annual %>%
    select(
      year, antecedent_n_casts, antecedent_first_date, antecedent_last_date,
      schmidt_n_casts, schmidt_first_date, schmidt_last_date,
      pm25_days_available, pm25_daily_coverage, hms_smoke_days,
      qualifying_smoke_days,
      panel_b_pm_window_first_date, panel_b_pm_window_last_date,
      panel_b_pm_window_days, panel_b_pm_days_available,
      panel_b_pm_daily_coverage, panel_b_qualifying_smoke_days,
      biological_n_dates, biological_first_date, biological_last_date,
      biological_span_days, biological_coverage_adequate,
      environmental_point_available, biological_response_available,
      biological_observed_zero
    ),
  file.path(OUT_DIR, "figure_5_state_space_coverage.csv")
)

smoke_sensitivity <- state_space_annual %>%
  select(
    year,
    panel_b_aug_dec_integrated_pm25 = integrated_pm25_aug_dec_ug_d_m3,
    event_smoke_day_count = qualifying_smoke_days,
    event_primary_excess = integrated_pm25_excess_ug_d_m3,
    event_cumulative_hms = cumulative_pm25_hms_smoke_days_ug_d_m3,
    event_calendar_baseline_excess = cumulative_pm25_excess_calendar_baseline_ug_d_m3
  ) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(rank_descending = min_rank(desc(value))) %>%
  ungroup()
biology_sensitivity <- state_space_annual %>%
  select(
    year,
    temporal_auc = biological_temporal_auc_cells_d_m2,
    maximum = biological_maximum_cells_m2,
    median_positive = biological_median_positive_cells_m2
  ) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(rank_descending = min_rank(desc(value))) %>%
  ungroup()

sensitivity_export <- bind_rows(
  loo_diagnostics %>%
    transmute(
      diagnostic = "reference_pca_leave_one_year_out",
      omitted_year, metric = "PC1", value = pc1_variance_explained,
      loading_cosine_similarity = pc1_loading_cosine_similarity,
      score_correlation = overlapping_score_correlation
    ),
  tibble(
    diagnostic = "reference_pca_vs_all_years_pca", omitted_year = NA_integer_,
    metric = "PC1", value = pca_variance[1],
    loading_cosine_similarity = reference_vs_all_loading_cosine,
    score_correlation = reference_vs_all_score_correlation
  ),
  smoke_sensitivity %>%
    transmute(
      diagnostic = "wildfire_exposure_metric", year, metric, value,
      rank_descending
    ),
  biology_sensitivity %>%
    transmute(
      diagnostic = "biological_response_metric", year, metric, value,
      rank_descending
    )
)
write_csv(sensitivity_export,
          file.path(OUT_DIR, "figure_5_state_space_sensitivity.csv"))

smoke_focal <- smoke_sensitivity %>% filter(year %in% c(2020L, 2021L)) %>%
  select(year, metric, value) %>% pivot_wider(names_from = year, values_from = value)
event_intensity_order_robust <- nrow(filter(
  smoke_focal,
  metric %in% c(
    "event_primary_excess", "event_cumulative_hms",
    "event_calendar_baseline_excess"
  )
)) == 3L &&
  all(is.finite(smoke_focal$`2021`), is.finite(smoke_focal$`2020`)) &&
  all(filter(
    smoke_focal,
    metric %in% c(
      "event_primary_excess", "event_cumulative_hms",
      "event_calendar_baseline_excess"
    )
  )$`2021` > filter(
    smoke_focal,
    metric %in% c(
      "event_primary_excess", "event_cumulative_hms",
      "event_calendar_baseline_excess"
    )
  )$`2020`)
smoke_day_count_2020_exceeds_2021 <- smoke_focal %>%
  filter(metric == "event_smoke_day_count") %>%
  summarise(result = n() == 1L && `2020` > `2021`) %>%
  pull(result)
panel_b_integrated_2021_exceeds_2020 <- smoke_focal %>%
  filter(metric == "panel_b_aug_dec_integrated_pm25") %>%
  summarise(result = n() == 1L && `2021` > `2020`) %>%
  pull(result)
biological_focal <- biology_sensitivity %>%
  filter(year %in% FOCAL_YEARS) %>%
  select(year, metric, value) %>%
  pivot_wider(names_from = year, values_from = value)
auc_and_maximum <- biological_focal %>% filter(metric %in% c("temporal_auc", "maximum"))
median_positive_focal <- biological_focal %>% filter(metric == "median_positive")
year_2020_observed_zero <- state_space_annual %>%
  filter(year == 2020L) %>%
  summarise(confirmed = n() == 1L && biological_coverage_adequate %in% TRUE &&
              biological_observed_zero %in% TRUE) %>%
  pull(confirmed)
biological_interpretation_robust <-
  nrow(auc_and_maximum) == 2L &&
  all(is.finite(as.matrix(auc_and_maximum[, c("2011", "2020", "2021")]))) &&
  all(auc_and_maximum$`2021` > auc_and_maximum$`2011`) &&
  all(auc_and_maximum$`2021` > auc_and_maximum$`2020`) &&
  nrow(median_positive_focal) == 1L &&
  is.finite(median_positive_focal$`2011`) &&
  is.finite(median_positive_focal$`2021`) &&
  median_positive_focal$`2021` > median_positive_focal$`2011` &&
  year_2020_observed_zero

plot_state <- state_space_annual %>% filter(environmental_point_available) %>%
  mutate(
    biological_positive = biological_response_available &
      !biological_observed_zero & biological_temporal_auc_cells_d_m2 > 0,
    year_group = case_when(
      year %in% c(2020L, 2021L) ~ "2020/2021",
      biological_positive ~ "Leptolyngbya present",
      TRUE ~ "Background"
    ),
    year_label = if_else(biological_positive | year == 2020L,
                         as.character(year), NA_character_),
    integrated_leptolyngbya_abundance_trillion =
      biological_temporal_auc_cells_d_m2 / 1e12
  )
state_colours <- c(
  "Background" = "grey68", "Leptolyngbya present" = "#009E73",
  "2020/2021" = "#E69F00"
)
positive_auc <- plot_state$biological_temporal_auc_cells_d_m2[
  is.finite(plot_state$biological_temporal_auc_cells_d_m2) &
    plot_state$biological_temporal_auc_cells_d_m2 > 0
]
auc_breaks <- if (length(positive_auc) >= 2L) {
  unique(c(
    as.numeric(quantile(positive_auc, c(0.25, 0.60), type = 8)),
    max(positive_auc)
  ) / 1e12)
} else pretty(positive_auc / 1e12, n = 3)

p_state <- ggplot(plot_state, aes(antecedent_pc1, integrated_pm25_aug_dec_ug_d_m3)) +
  geom_point(
    data = filter(plot_state, !biological_response_available),
    aes(colour = year_group), shape = 21, fill = NA, size = 2.1, stroke = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(plot_state, biological_response_available, biological_observed_zero),
    aes(colour = year_group), shape = 4, size = 2.2, stroke = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(
      plot_state, biological_response_available, !biological_observed_zero,
      is.finite(integrated_leptolyngbya_abundance_trillion)
    ),
    aes(size = integrated_leptolyngbya_abundance_trillion,
        fill = year_group, colour = year_group),
    shape = 21, alpha = 0.82, stroke = 0.55
  ) +
  geom_text_repel(
    data = filter(plot_state, !is.na(year_label), year != 2021L),
    aes(label = year_label, colour = year_group),
    size = lo_geom_text_size(7.5), family = BASE_FAMILY,
    min.segment.length = 0, box.padding = 0.22, point.padding = 0.25,
    seed = 2021, max.overlaps = Inf, show.legend = FALSE
  ) +
  geom_text(
    data = filter(plot_state, year == 2021L, !is.na(year_label)),
    aes(label = year_label, colour = year_group),
    nudge_y = -430, hjust = 0.5, vjust = 0.5,
    size = lo_geom_text_size(7.5), family = BASE_FAMILY,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = state_colours, guide = "none") +
  scale_colour_manual(values = state_colours, guide = "none") +
  scale_size_area(
    max_size = 5.5, breaks = auc_breaks,
    labels = function(x) {
      numeric_labels <- sub("\\.?0+$", "", formatC(x, format = "f", digits = 2))
      as.expression(lapply(seq_along(x), function(i) {
        if (i == length(x)) {
          bquote(.(numeric_labels[i])~"(10"^12~"cells d m"^{-2}*")")
        } else {
          numeric_labels[i]
        }
      }))
    },
    name = expression("Integrated"~italic(Leptolyngbya)~"abundance"),
    guide = guide_legend(
      direction = "horizontal", title.position = "top", title.hjust = 0.5,
      label.position = "right", label.hjust = 0,
      nrow = 1, byrow = TRUE,
      override.aes = list(
        size = c(1.3, 2.1, 3.8), fill = "grey75", colour = "grey30"
      )
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = scientific_axis_labels,
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(
    tag = "b", title = NULL,
    x = "Pre-existing thermal state (PC1)",
    y = expression("Integrated PM"[2.5]~"("*mu*"g d m"^{-3}*")")
  ) +
  figure_theme +
  theme(
    legend.position = "bottom",
    legend.justification = c(0.70, 0.5),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.title = ggplot2::element_text(size = 8.1, face = "plain", lineheight = 0.9,
                                         margin = margin(0, 0, 0.2, 0, "mm")),
    legend.text = ggplot2::element_text(size = 8.1),
    legend.key.width = unit(1.7, "mm"),
    legend.key.height = unit(2.3, "mm"),
    legend.spacing.x = unit(0.05, "mm"),
    legend.margin = margin(-1.5, 0, -1.0, 0, "mm"),
    axis.text.y = ggplot2::element_text(
      size = 7.8, margin = margin(0, 0.4, 0, 0, "mm")
    ),
    axis.title.y = ggplot2::element_text(
      size = 8.0, face = "bold", margin = margin(0, 0.8, 0, 0, "mm")
    ),
    plot.tag.position = c(0.985, 0.975),
    plot.tag.location = "panel",
    plot.tag = element_text(size = 10, face = "bold", hjust = 1, vjust = 1),
    plot.margin = margin(0.8, 1.2, 0.4, 0.2, "mm")
  )

# Explicitly fail if a secondary y axis is ever added to either trajectory.
has_secondary_y <- function(plot) {
  any(vapply(plot$scales$scales, function(scale) {
    "y" %in% scale$aesthetics && !is.null(scale$secondary.axis) &&
      !inherits(scale$secondary.axis, "waiver")
  }, logical(1)))
}
if (has_secondary_y(p_community)) {
  stop("A secondary y-axis was introduced into Panel A.")
}

# Confirm that every plotted Panel C recovery estimate is exactly the value written to
# the machine-readable export, rather than comparing two aliases in memory.
recovery_export_written <- read_csv(
  file.path(OUT_DIR, "figure_5_recovery_synthesis.csv"),
  show_col_types = FALSE
)
plot_export_check <- recovery_plot_data %>%
  select(pathway, variable, plotted = selected_recovery_months) %>%
  left_join(
    recovery_export_written %>%
      select(pathway, variable, exported = selected_recovery_months),
    by = c("pathway", "variable")
  )
if (
  any(!is.finite(plot_export_check$exported)) ||
    any(!near(plot_export_check$plotted, plot_export_check$exported), na.rm = TRUE)
) {
  stop("A plotted recovery value disagrees with figure_5_recovery_synthesis.csv.")
}

# ---- Figure, current configuration, caption, methods, and diagnostics -------
draw_figure_5 <- function() {
  grid::grid.newpage()
  top_row_weight <- 4.6
  separator_row_weight <- 0.30
  bottom_row_weight <- 7.4
  layout <- grid::grid.layout(
    nrow = 3, ncol = 2,
    heights = grid::unit(
      c(top_row_weight, separator_row_weight, bottom_row_weight), "null"
    ),
    widths = grid::unit(c(1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(p_community, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_state, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_recovery, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  separator_y <- (bottom_row_weight + separator_row_weight / 2) /
    (top_row_weight + separator_row_weight + bottom_row_weight)
  grid::grid.lines(
    x = grid::unit(c(0.015, 0.985), "npc"),
    y = grid::unit(c(separator_y, separator_y), "npc"),
    gp = grid::gpar(col = "grey55", lwd = 0.7)
  )
  grid::upViewport()
  invisible()
}

main_png <- file.path(OUT_DIR, "figure_5_final.png")
main_pdf <- file.path(OUT_DIR, "figure_5_final.pdf")
ragg::agg_png(main_png, width = 12.2, height = 14.6, units = "cm",
              res = LO_DPI_LINE, background = "white")
draw_figure_5()
invisible(dev.off())
grDevices::cairo_pdf(main_pdf, width = 12.2 / 2.54, height = 14.6 / 2.54,
                     bg = "white")
draw_figure_5()
invisible(dev.off())
message("Saved (12.2 x 14.6 cm, 1200 dpi): figure_5_final.png")

# Standalone direct Bray-Curtis version requested for metric verification. The
# main Panel A uses this same authoritative distance; this export makes that
# identity explicit without substituting the archived PCoA-coordinate metric.
ggsave(
  file.path(OUT_DIR, "figure_5_panel_a_direct_bray_curtis.png"), p_community,
  width = 12.7, height = 8.0, units = "cm", dpi = LO_DPI_LINE,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  file.path(OUT_DIR, "figure_5_panel_a_direct_bray_curtis.pdf"), p_community,
  width = 12.7, height = 8.0, units = "cm", device = cairo_pdf, bg = "white"
)
if (!file.copy(main_png, file.path(MIRROR_DIR, "figure_5_final.png"), overwrite = TRUE)) {
  stop("Failed to update the canonical figures_v2 PNG mirror.")
}
if (!file.copy(main_pdf, file.path(MIRROR_DIR, "figure_5_final.pdf"), overwrite = TRUE)) {
  stop("Failed to update the canonical figures_v2 PDF mirror.")
}
if (RENDER_ONLY) {
  message("Render-only Figure 5 build completed; analytical and project-state exports were not rewritten.")
  quit(save = "no", status = 0, runLast = FALSE)
}

config <- list(
  analysis = "Figure 5 quantitative synthesis",
  biological_data_used_to_define_windows_or_pcas = FALSE,
  primary_antecedent_rule = "fixed_jul01_aug13",
  primary_antecedent_start_month_day = ANTECEDENT_START_MD,
  primary_antecedent_end_month_day = ANTECEDENT_END_MD,
  panel_b_pm_window_rule = "fixed_aug01_dec31_matched_to_biological_response",
  panel_b_pm_start_month_day = RESPONSE_START_MD,
  panel_b_pm_end_month_day = RESPONSE_END_MD,
  event_sensitivity_window_rule = "fixed_aug14_oct21",
  event_sensitivity_start_month_day = EXPOSURE_START_MD,
  event_sensitivity_end_month_day = EXPOSURE_END_MD,
  smoke_day_definition = paste0(
    "NOAA HMS smoke over Lake Tahoe with composite daily PM2.5 >= ",
    PM_THRESHOLD, " ug m^-3"
  ),
  minimum_annual_daily_coverage = MINIMUM_PM_COVERAGE,
  primary_exposure_metric = paste(
    "sum of daily PM2.5 on NOAA HMS smoke days with PM2.5 >= 9.1 ug m^-3",
    "during August 1-December 31; threshold is not subtracted"
  ),
  panel_b_exposure_axis_transform = "linear native-unit scale",
  exposure_sensitivity_metrics = c(
    "qualifying smoke-day count during August 14-October 21",
    "sum(max(PM2.5 - 9.1, 0)) on qualifying HMS smoke days",
    "cumulative PM2.5 on HMS smoke days",
    "cumulative PM2.5 excess above the reference-year calendar-day median on HMS smoke days"
  ),
  antecedent_variables = state_variables,
  mixed_layer_depth_transform = "mixed_layer_shallowness_m = -1 * mixed_layer_depth_m",
  pca_fit_years = pca_reference$year,
  pca_excluded_focal_years = FOCAL_YEARS,
  pca_projection_years = pca_complete$year,
  pca_orientation = paste(
    "positive net alignment with warmer surface water, greater mixed-layer thermal contrast,",
    "and shallower mixed layer"
  ),
  schmidt_stability_role = paste(
    "external validation using years with a qualifying full-depth cast;",
    "not a PCA input; no interpolation or deep-water extrapolation"
  ),
  biological_response_window = "August 1-December 31",
  biological_depth_integration_m = CORE_DEPTHS_M,
  plotted_biological_response_metric = "trapezoidal temporal AUC of date-level 5-90 m depth-integrated Leptolyngbya abundance",
  biological_minimum_dates = MINIMUM_RESPONSE_DATES,
  biological_minimum_span_days = MINIMUM_RESPONSE_SPAN_DAYS,
  recovery_time_origin = as.character(FIRE_START),
  active_fire_end = as.character(FIRE_END),
  final_figure = main_png
)
write_yaml(config, CONFIG_FILE)

pc1_loadings_text <- paste(
  paste0(pca_parameters$variable, " = ", sprintf("%.6f", pca_parameters$pc1_loading)),
  collapse = "; "
)
insufficient_biology_years <- state_space_annual %>%
  filter(!biological_response_available) %>% pull(year)
insufficient_environment_years <- state_space_annual %>%
  filter(!environmental_point_available) %>% pull(year)
schmidt_agrees <- schmidt_summary$pearson_correlation > 0 &
  schmidt_summary$spearman_correlation > 0
excluded_panel_a_rows <- panel_a_data_all %>%
  filter(!included_minimum_3_depths)
deep_2022_excluded <- excluded_panel_a_rows %>%
  filter(stratum == "60-105 m", year(date) == 2022L) %>%
  slice_max(value, n = 1, with_ties = FALSE)

caption_lines <- c(
  "Fig 5. Persistence and historical context of biological responses following the Caldor Fire.",
  paste0(
    "a) Phytoplankton community departure, expressed as Bray-Curtis dissimilarity between Hellinger-transformed abundance and a season-matched pre-fire centroid, in the upper (0-40 m) and lower (60-105 m) water column. Only dates with at least three sampled depths per zone are shown. Shaded bands and their boundaries show the 95% historical prediction intervals. Red shading and the red dashed line denote the active-fire period (14 August-21 October 2021) and fire ignition, respectively; black dashed lines mark the first observed statistical return."
  ),
  paste0(
    "b) Historical comparison of antecedent lake thermal state and integrated PM2.5 exposure. Antecedent state (PC1) summarizes mixed-layer thermal contrast, mixed-layer shallowness, and 0-10 m temperature from 1 July-13 August; Schmidt stability is evaluated separately where qualifying full-depth casts are available. Integrated PM2.5 is the sum of daily PM2.5 on NOAA HMS smoke days with PM2.5 ≥ 9.1 µg m⁻³ from 1 August-31 December. Bubble area is proportional to the temporal integral of depth-integrated *Leptolyngbya* abundance over the same August-December window. Green symbols identify years with positive *Leptolyngbya* abundance other than 2021, orange identifies the focal smoke years 2020 and 2021, and gray denotes other background years. Open points lack sufficient biological coverage, whereas crosses are adequately sampled observed zeros."
  ),
  paste0(
    "c) Recovery time from fire ignition for affected radiative, nutrient, runoff, and biological variables. Points show the first observed statistical return to the historical reference range, and horizontal lines span estimates under the one-, two-, and three-observation return criteria; these ranges are not confidence intervals. Upward triangles indicate return from a negative anomaly, and downward triangles indicate decline from a positive anomaly. The red dashed line marks the end of the active-fire period. Recovery estimates are interval-censored by sampling frequency and do not represent exact process termination times."
  )
)
writeLines(caption_lines, file.path(OUT_DIR, "figure_5_caption.txt"), useBytes = TRUE)

methods_lines <- c(
  "# Figure 5 quantitative-synthesis methods",
  "",
  "## Panel A: community departure",
  "Panel A reads the authoritative pathway resistance/recovery time series without recalculating its community endpoint and displays 2019-2025. The metric is direct Bray-Curtis dissimilarity between Hellinger-transformed taxon abundance and a fixed season-specific pre-fire centroid; it is not the archived additive-corrected PCoA-coordinate distance. A visit is displayed only when at least three distinct depths were sampled within its depth zone. The coverage diagnostic retains excluded rows and records the exclusion reason. The 25 July 2022 deep value is excluded because only 75 m and two positive taxa were represented. The band and its faint solid bounds are the 95% historical prediction interval. The first-observed-return date comes from pathway_resistance_resilience_markers.csv.",
  "",
  "## Panel B: historical state space",
  "Antecedent state uses 1 July-13 August in every year. Annual medians are calculated for mixed-layer thermal contrast, mixed-layer shallowness (-1 times mixed-layer depth), and mean 0-10 m temperature. The three variables are standardized with centers and scales from background years only. PCA excludes the predeclared focal contrast years 2011, 2020, and 2021 from fitting and projects those years into the fixed space. The PC1 sign is oriented toward positive net alignment with warmer surface water, greater thermal contrast, and a shallower mixed layer. Biological observations are not used to select variables, windows, fit years, or orientation.",
  "Mixed-layer shallowness is defined as -1 times mixed-layer depth, so a shallower observed mixed layer has a less-negative input value. Schmidt stability is not included in the PCA. It is compared with PC1 as an external validation only for years with a qualifying full-depth cast in the fixed antecedent window; profiles are not extrapolated or interpolated to fill missing stability values.",
  paste0("Panel B PM2.5 and the biological response use the identical 1 August-31 December window in every year. Integrated PM2.5 is the sum of daily PM2.5 on NOAA HMS smoke days with daily PM2.5 >= ", PM_THRESHOLD, " ug m^-3; the operational threshold qualifies days but is not subtracted from the sum. At least ", 100 * MINIMUM_PM_COVERAGE, "% daily PM2.5 coverage is required. The plotted axis is linear and retains native units. Event-window sensitivities retain the 14 August-21 October qualifying smoke-day count, integrated excess above 9.1 ug m^-3, cumulative PM2.5 on all HMS smoke days, and positive PM2.5 excess above the reference-year calendar-day median."),
  paste0("At each biological date, *Leptolyngbya* abundance is trapezoidally integrated across fixed depths ", paste(CORE_DEPTHS_M, collapse = ", "), " m. The displayed response is temporal trapezoidal AUC from 1 August-31 December between observed dates. A year requires at least ", MINIMUM_RESPONSE_DATES, " dates spanning at least ", MINIMUM_RESPONSE_SPAN_DAYS, " days. Inadequately sampled years retain their environmental point and are not assigned zero. Sensitivities are maximum depth-integrated abundance and median positive depth-integrated abundance."),
  "",
  "## Panel C: cross-pathway recovery",
  paste0("Time is measured from ", FIRE_START, ". The standard estimate is the first of two consecutive post-peak observations inside the relevant reference envelope. Atmospheric NO3 and NH4 deposition use the first post-peak observation; deep chlorophyll uses the first trailing three-observation mean inside the corresponding mean envelope. Point shape records whether the selected peak anomaly was negative and returned upward or positive and declined toward the reference envelope; it no longer encodes the recovery-rule exception. Thin horizontal error bars show the minimum and maximum of available one-, two-, and three-observation estimates and are omitted when fewer than two sensitivity values exist. They are sensitivity ranges, not confidence intervals. These are observation-based, sampling-interval-censored estimates, not exact process termination times. Unaffected and pre-existing variables are exported but not ranked."),
  "",
  "No interaction, causal, necessary-condition, or sufficient-condition model is fit."
)
writeLines(methods_lines, file.path(OUT_DIR, "figure_5_methods.md"), useBytes = TRUE)

diagnostic_lines <- c(
  "# Figure 5 diagnostics",
  "",
  "## Validation checks",
  "- PASS: all antecedent observations are within July 1-August 13.",
  "- PASS: all exposure observations are within August 14-October 21.",
  "- PASS: all Panel B integrated-PM2.5 observations are within August 1-December 31.",
  "- PASS: inadequately sampled biological years are not classified as observed zero.",
  "- PASS: 2011, 2020, and 2021 are excluded from the reference PCA fit.",
  "- PASS: plotted Panel C recovery estimates equal the recovery export.",
  "- PASS: Panel A has no secondary y-axis.",
  "- PASS: every displayed Panel A visit has at least three sampled depths in its zone.",
  "",
  "## PCA",
  paste0("- Reference fit years: ", paste(pca_reference$year, collapse = ", "), "."),
  paste0("- PC1 loadings: ", pc1_loadings_text, "."),
  paste0("- PC1 variance explained: ", sprintf("%.2f%%", 100 * pca_variance[1]), "."),
  paste0("- Leave-one-background-year-out loading cosine range: ", sprintf("%.4f-%.4f", min(loo_diagnostics$pc1_loading_cosine_similarity), max(loo_diagnostics$pc1_loading_cosine_similarity)), "."),
  paste0("- Leave-one-background-year-out score-correlation range: ", sprintf("%.4f-%.4f", min(loo_diagnostics$overlapping_score_correlation), max(loo_diagnostics$overlapping_score_correlation)), "."),
  paste0("- Reference-only versus all-years score correlation: ", sprintf("%.4f", reference_vs_all_score_correlation), "; loading cosine similarity: ", sprintf("%.4f", reference_vs_all_loading_cosine), "."),
  "",
  "## Schmidt stability validation",
  paste0("- n = ", schmidt_summary$n, "; Pearson r = ", sprintf("%.3f", schmidt_summary$pearson_correlation), " (p = ", format.pval(schmidt_summary$pearson_p_value, digits = 3), "); Spearman rho = ", sprintf("%.3f", schmidt_summary$spearman_correlation), " (p = ", format.pval(schmidt_summary$spearman_p_value, digits = 3), ")."),
  paste0("- PC1 is positively associated with independently calculated Schmidt stability: ", if_else(schmidt_agrees, "yes", "no"), ". Schmidt stability is not a PCA input; the validation uses only qualifying full-depth casts in the fixed antecedent window. The main axis retains the neutral label 'Antecedent thermal state (PC1)'."),
  "",
  "## Coverage and sensitivity",
  paste0("- Years excluded from Panel B for incomplete environmental coordinates: ", if_else(length(insufficient_environment_years) == 0L, "none", paste(insufficient_environment_years, collapse = ", ")), "."),
  paste0("- Years lacking the required biological coverage: ", if_else(length(insufficient_biology_years) == 0L, "none", paste(insufficient_biology_years, collapse = ", ")), "."),
  paste0("- The coverage-limited deep-community observation on ", deep_2022_excluded$date, " (Bray-Curtis = ", sprintf("%.3f", deep_2022_excluded$value), "; ", deep_2022_excluded$sampled_depth_count, " sampled depth; ", deep_2022_excluded$positive_taxa_count, " positive taxa) is excluded by the three-depth Panel A criterion."),
  paste0("- Panel B August-December integrated PM2.5 is higher in 2021 than 2020: ", panel_b_integrated_2021_exceeds_2020, "."),
  paste0("- In the shorter event window, 2020 has more qualifying smoke days than 2021: ", smoke_day_count_2020_exceeds_2021, "."),
  paste0("- 2021 exceeds 2020 under all three event-window PM2.5 intensity/excess metrics: ", event_intensity_order_robust, "."),
  paste0("- The focal-year biological interpretation is unchanged across AUC, maximum, and median-positive summaries where estimable: ", biological_interpretation_robust, ". The median-positive metric is undefined for the adequately sampled observed-zero year 2020; its zero classification is confirmed by AUC and maximum abundance."),
  "- SIMPER is not displayed in the main figure; its prior outputs remain unchanged for Supporting Information.",
  "",
  "## Publication-readiness checks",
  "- Figure dimensions are 12.2 x 14.6 cm; the PNG is 1200 dpi and the PDF is vector.",
  "- Manual final-size inspection is required after each regeneration; the script does not auto-certify visual legibility."
)
writeLines(diagnostic_lines, file.path(OUT_DIR, "figure_5_diagnostics.md"), useBytes = TRUE)

# Compatibility copies keep the existing processed state-space location useful
# without deleting any earlier exports.
write_csv(state_space_annual,
          file.path(PROC_DIR, "primary_state_space_annual_summary.csv"))
write_csv(pca_parameters,
          file.path(PROC_DIR, "primary_antecedent_pca_parameters.csv"))
write_csv(loo_diagnostics,
          file.path(PROC_DIR, "primary_antecedent_pca_leave_one_year_out.csv"))

message("Saved revised quantitative-synthesis Figure 5 and supporting exports to ", OUT_DIR)
