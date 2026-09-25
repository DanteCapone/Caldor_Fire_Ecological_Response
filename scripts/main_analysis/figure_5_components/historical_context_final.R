# =============================================================================
# figure_5_ltp_abundance_historical_context_final.R
#
# Final Figure 5 historical-context variants.
#
# Produces:
#   1. Upper + Lower version with paired background violins.
#   2. Upper-only version with no displayed depth-zone label.
#
# Both versions use a seven-month August-February window, a full-width
# historical 95th-percentile line, and a faint grey region below that line.
# Panel C has no title. Panel tags a, b, and c share the same native size.
# =============================================================================

# Rebuild the native core figure objects and calculate the depth-resolved
# historical data. The sourced scripts also retain their standalone outputs.
source("scripts/main_analysis/figure_5_components/figure_5_revised.R")
source("scripts/main_analysis/figure_5_components/historical_context_depth.R")

# ---- Final output paths ------------------------------------------------------
UPPER_PANEL_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_upper_only_panel.png"
)
UPPER_COMPOSITE_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_upper_only.png"
)

# ---- Shared final styling ----------------------------------------------------
FINAL_TAG_SIZE <- 10
REFERENCE_FILL <- "grey75"
REFERENCE_ALPHA <- 0.13
REFERENCE_LINE <- "grey35"

full_width_q95 <- as.numeric(
  quantile(
    background_maxima$maximum_departure,
    probs = 0.95,
    na.rm = TRUE,
    type = 8
  )
)

# Reuse the existing grouped positioning so Upper and Lower remain paired.
final_depth_dodge <- position_dodge(width = 0.50)
final_raw_jitter_dodge <- position_jitterdodge(
  jitter.width = 0.055,
  jitter.height = 0,
  dodge.width = 0.50,
  seed = 42
)

# =============================================================================
# VERSION 1: Upper + Lower
# =============================================================================

final_depth_panel <- ggplot() +
  # The light-grey region denotes values below the pooled historical 95th
  # percentile of background August-February maxima.
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = full_width_q95,
    fill = REFERENCE_FILL,
    alpha = REFERENCE_ALPHA
  ) +
  geom_hline(
    yintercept = full_width_q95,
    colour = REFERENCE_LINE,
    linetype = "dashed",
    linewidth = 0.55 * LO_FIG_SCALE
  ) +
  annotate(
    "text",
    x = 1.08,
    y = full_width_q95,
    label = "Historical 95th percentile",
    hjust = 0,
    vjust = -0.45,
    size = 1.75 * LO_FIG_SCALE,
    family = BASE_FAMILY,
    colour = REFERENCE_LINE
  ) +
  geom_violin(
    data = background_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = interaction(comparison_group, depth_zone),
      fill = event_depth
    ),
    position = final_depth_dodge,
    width = 0.82,
    scale = "width",
    trim = FALSE,
    linewidth = 0.45 * LO_FIG_SCALE,
    colour = "grey25",
    alpha = 0.78
  ) +
  geom_point(
    data = background_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = depth_zone,
      shape = depth_zone,
      fill = event_depth
    ),
    position = position_jitterdodge(
      jitter.width = 0.04,
      jitter.height = 0,
      dodge.width = 0.50,
      seed = 42
    ),
    size = 1.2 * LO_FIG_SCALE,
    stroke = 0.2 * LO_FIG_SCALE,
    colour = "grey35",
    alpha = 0.38
  ) +
  geom_point(
    data = focal_values,
    aes(
      x = comparison_group,
      y = departure,
      group = depth_zone,
      shape = depth_zone,
      fill = event_depth
    ),
    position = final_raw_jitter_dodge,
    size = 1.55 * LO_FIG_SCALE,
    stroke = 0.25 * LO_FIG_SCALE,
    colour = "grey20",
    alpha = 0.18
  ) +
  geom_point(
    data = focal_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = depth_zone,
      shape = depth_zone,
      fill = event_depth
    ),
    position = final_depth_dodge,
    size = 3.35 * LO_FIG_SCALE,
    stroke = 0.75 * LO_FIG_SCALE,
    colour = "grey10"
  ) +
  geom_text(
    data = focal_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = depth_zone,
      label = sprintf("%.2f", maximum_departure)
    ),
    position = final_depth_dodge,
    vjust = -1.10,
    size = 1.85 * LO_FIG_SCALE,
    family = BASE_FAMILY,
    fontface = "bold",
    colour = "grey10"
  ) +
  scale_shape_manual(
    values = DEPTH_SHAPES,
    name = "Depth zone",
    labels = c("Upper", "Lower")
  ) +
  scale_fill_manual(values = EVENT_DEPTH_COLORS, guide = "none") +
  scale_x_discrete(
    limits = GROUP_LEVELS,
    labels = c(
      "Background non-fire years" = "Background non-fire\nyears",
      "2011 Leptolyngbya" = "2011\nLeptolyngbya",
      "2020 smoke" = "2020\nsmoke exposure",
      "2021 Caldor" = "2021\nCaldor fire"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    tag = "c",
    title = NULL,
    x = NULL,
    y = "Maximum Bray-Curtis distance\nto seasonal background centroid"
  ) +
  guides(
    shape = guide_legend(
      direction = "horizontal",
      override.aes = list(
        size = 2.4 * LO_FIG_SCALE,
        fill = c("grey78", "grey45"),
        alpha = 1
      )
    )
  ) +
  theme_bw(base_size = 8, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      colour = "#EBEBEB",
      linewidth = LO_LW_THIN
    ),
    axis.text = element_text(size = LO_FS_AXIS_TEXT, colour = "grey20"),
    axis.title.y = element_text(size = LO_FS_AXIS_TITLE, face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.tag = element_text(size = FINAL_TAG_SIZE, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(size = 6.1 * LO_FIG_SCALE, face = "bold"),
    legend.text = element_text(size = 5.8 * LO_FIG_SCALE),
    legend.key.width = unit(4.5 * LO_FIG_SCALE, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 3, 1, 2, "mm")
  )

# Keep the existing native core plots but enlarge their panel tags. Because the
# core is combined natively here, a, b, and c use the same physical font size.
final_full_core <- results$abundance$fig &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = FINAL_TAG_SIZE,
      family = BASE_FAMILY
    )
  )

final_full_composite <- wrap_elements(full = final_full_core) /
  final_depth_panel +
  plot_layout(heights = c(13.5, 7.0))

agg_png(PANEL_FILE, width = 17.8, height = 7.0, units = "cm", res = LO_DPI)
print(final_depth_panel)
invisible(dev.off())

agg_png(COMPOSITE_FILE, width = 17.8, height = 20.5,
        units = "cm", res = LO_DPI)
print(final_full_composite)
invisible(dev.off())

# =============================================================================
# VERSION 2: Upper-only, with no displayed depth-zone label
# =============================================================================

# ---- Upper-only panel a: resistance trajectory -------------------------------
upper_distance <- results$abundance$distance %>%
  filter(depth_bin == "Surface", n_ref >= MIN_PREFIRE_REFS, !is.na(dist)) %>%
  arrange(date)

upper_lepto <- ltp_site %>%
  add_depth_bin() %>%
  filter(depth_bin == "Surface", taxon == LEPTO_TAXON) %>%
  group_by(date) %>%
  summarise(lepto_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Figure 5a and the pathway table use the same seasonal-mean recovery result.
# Reading the audited table value here prevents the former mismatch caused by
# Figure 5a reporting time since ignition while the table reported time since
# the seasonal peak.
recovery_summary_file <- file.path(
  "data", "processed", "pathway_resistance_resilience_summary.csv"
)
if (!file.exists(recovery_summary_file)) {
  stop("Run scripts/pathway_resistance_resilience.R before Figure 5 so recovery labels are synchronized.")
}
recovery_summary_shared <- read_csv(recovery_summary_file, show_col_types = FALSE)
shared_recovery_label <- function(depth_zone) {
  value <- recovery_summary_shared %>%
    filter(Variable == paste0("Community composition (Bray-Curtis; ", depth_zone, ")")) %>%
    pull(`Recovery time`)
  if (length(value) != 1L) stop("Missing shared recovery result for ", depth_zone)
  if (is.na(value) || value == "NA") "Recovery: not observed" else
    paste0("Recovery: ", value)
}
upper_recovery_label <- shared_recovery_label("0-40 m")
lower_recovery_label <- shared_recovery_label("60-105 m")

upper_distance_plot <- make_distance_subplot(
  upper_distance,
  show_y = TRUE,
  lepto_dat = upper_lepto,
  recovery_label_override = upper_recovery_label
) +
  labs(title = NULL) +
  theme(plot.title = element_blank())

upper_y_title <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = "Bray-Curtis community dissimilarity\nto pre-fire centroid",
    angle = 90, size = 2.1, fontface = "bold", family = BASE_FAMILY
  ) +
  theme_void()

upper_lepto_title <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = expression(ln(italic("Leptolyngbya")~" abundance (cells L"^-1*")")),
    angle = -90, size = 2.1, colour = LEPTO_COL,
    family = BASE_FAMILY
  ) +
  theme_void()

upper_panel_a <- upper_y_title + upper_distance_plot + upper_lepto_title +
  plot_layout(widths = c(0.105, 1, 0.055))

# Alternate panel a: compare every sample with one pooled pre-fire abundance
# centroid and show a constant historical mean and 95% CI. Fire-affected years
# are excluded from both the centroid and its reference distribution.
make_constant_prefire_centroid <- function(zone_name = "Surface") {
  comm <- ltp_site %>%
    add_depth_bin() %>%
    filter(depth_bin == zone_name) %>%
    group_by(date, taxon) %>%
    summarise(value = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = value, values_fill = 0) %>%
    arrange(date)
  dates <- comm$date
  mat <- comm %>% select(-date) %>% as.matrix()
  keep <- rowSums(mat) > 0
  dates <- dates[keep]
  mat <- mat[keep, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  hell <- vegan::decostand(mat, method = "hellinger")
  hist <- dates < CALDOR_START &
    !lubridate::year(dates) %in% c(2011L, 2012L, 2020L, 2021L)
  if (sum(hist) < 2L) stop("Insufficient undisturbed pre-fire data for constant Figure 5 reference.")
  centroid <- colMeans(hell[hist, , drop = FALSE])
  bray_to_centroid <- function(x) {
    as.numeric(vegan::vegdist(rbind(x, centroid), method = "bray"))
  }
  distances <- apply(hell, 1, bray_to_centroid)
  hist_distances <- distances[hist]
  expected <- mean(hist_distances)
  ci_half <- stats::qt(0.975, df = sum(hist) - 1L) *
    stats::sd(hist_distances) * sqrt(1 + 1 / sum(hist))
  recovery_source <- if (zone_name == "Surface") upper_distance else lower_distance
  shared_recovery <- recovery_source %>%
    summarise(
      recovery_date = first(na.omit(recovery_date), default = as.Date(NA)),
      recovery_months = first(na.omit(recovery_months), default = NA_real_)
    )
  tibble(
    date = dates,
    dist = distances,
    pre_05 = pmax(0, expected - ci_half),
    pre_median = expected,
    pre_95 = expected + ci_half,
    n_ref = sum(hist),
    depth_bin = zone_name,
    recovery_date = shared_recovery$recovery_date,
    recovery_months = shared_recovery$recovery_months
  )
}

upper_distance_constant <- make_constant_prefire_centroid("Surface")
upper_distance_plot_constant <- make_distance_subplot(
  upper_distance_constant,
  show_y = TRUE,
  lepto_dat = upper_lepto,
  recovery_label_override = upper_recovery_label
) +
  labs(title = NULL) +
  theme(plot.title = element_blank())
upper_panel_a_constant <- upper_y_title + upper_distance_plot_constant +
  upper_lepto_title + plot_layout(widths = c(0.105, 1, 0.055))


# Standalone panel-a versions with explicit calendar-year limits. The plotting
# function uses a dynamic y range, so observations above the former fixed 0.9
# ceiling remain visible and no longer create artificial gaps in the line.
make_upper_panel_a_window <- function(start_date, end_date) {
  distance_plot <- make_distance_subplot(
    upper_distance, show_y = TRUE, lepto_dat = upper_lepto,
    date_start = as.Date(start_date), date_end = as.Date(end_date),
    recovery_label_override = upper_recovery_label
  ) + labs(title = NULL) + theme(plot.title = element_blank())
  upper_y_title + distance_plot + upper_lepto_title +
    plot_layout(widths = c(0.105, 1, 0.055))
}

upper_panel_a_2005_2025 <- make_upper_panel_a_window("2005-01-01", "2025-12-31")
upper_panel_a_2015_2025 <- make_upper_panel_a_window("2015-01-01", "2025-12-31")

for (spec in list(
  list(plot = upper_panel_a_2005_2025,
       file = file.path("figures", "figure_5_nmds_braycurtis", "figure_5a_2005_2025.png")),
  list(plot = upper_panel_a_2015_2025,
       file = file.path("figures", "figure_5_nmds_braycurtis", "figure_5a_2015_2025.png"))
)) {
  ragg::agg_png(spec$file, width = 17.8, height = 7.0, units = "cm", res = LO_DPI)
  print(spec$plot)
  invisible(dev.off())
}

# Supplemental panel a containing both LTP depth zones over the complete
# curated record.
lower_distance <- results$abundance$distance %>%
  filter(depth_bin == "Deep", n_ref >= MIN_PREFIRE_REFS, !is.na(dist)) %>%
  arrange(date)
lower_lepto <- ltp_site %>%
  add_depth_bin() %>%
  filter(depth_bin == "Deep", taxon == LEPTO_TAXON) %>%
  group_by(date) %>%
  summarise(lepto_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

lower_distance_constant <- make_constant_prefire_centroid("Deep")

supp_upper_plot <- make_distance_subplot(
  upper_distance_constant, show_y = TRUE, lepto_dat = upper_lepto,
  date_start = as.Date("2005-01-01"), date_end = as.Date("2025-12-31"),
  recovery_label_override = upper_recovery_label,
  date_breaks_override = "2 years", log_lepto = TRUE, show_bc_points = FALSE,
  y_lower = 0, recovery_position = "upper_left"
) + labs(tag = "a")
supp_lower_plot <- make_distance_subplot(
  lower_distance_constant, show_y = TRUE, lepto_dat = lower_lepto,
  date_start = as.Date("2005-01-01"), date_end = as.Date("2025-12-31"),
  recovery_label_override = lower_recovery_label,
  date_breaks_override = "2 years", log_lepto = TRUE, show_bc_points = FALSE,
  y_lower = 0, recovery_position = "upper_left"
) + labs(tag = "b")
supp_upper_lower <- upper_y_title + (supp_upper_plot / supp_lower_plot) +
  upper_lepto_title + plot_layout(widths = c(0.085, 1, 0.050))
supp_fig5_dir <- file.path("figures", "supplemental", "figure_5")
dir.create(supp_fig5_dir, recursive = TRUE, showWarnings = FALSE)
supp_upper_lower_file <- file.path(
  supp_fig5_dir, "supplemental_figure_5a_upper_lower_2005_2025.png"
)
ragg::agg_png(supp_upper_lower_file, width = 12.7, height = 13.5,
              units = "cm", res = LO_DPI)
print(supp_upper_lower)
invisible(dev.off())
ggsave(
  sub("\\.png$", ".pdf", supp_upper_lower_file), supp_upper_lower,
  width = 12.7, height = 13.5, units = "cm", device = cairo_pdf,
  bg = "white"
)

# ---- Upper-only panel b: SIMPER contributors --------------------------------
upper_simper <- results$abundance$simper %>%
  filter(depth_bin == "Surface") %>%
  slice_min(rank, n = 8, with_ties = FALSE) %>%
  arrange(contribution_pct) %>%
  mutate(
    y = row_number(),
    signed_contribution = if_else(
      higher_in == "Higher pre-fire",
      -contribution_pct,
      contribution_pct
    )
  )

upper_simper_max <- max(abs(upper_simper$signed_contribution), na.rm = TRUE) * 1.12

upper_panel_b <- ggplot(
  upper_simper,
  aes(x = signed_contribution, y = y, fill = higher_in)
) +
  geom_vline(
    xintercept = 0,
    colour = "grey25",
    linewidth = 0.30 * LO_FIG_SCALE
  ) +
  geom_col(
    width = 0.72,
    colour = "grey25",
    linewidth = 0.18 * LO_FIG_SCALE,
    orientation = "y"
  ) +
  scale_y_continuous(
    breaks = upper_simper$y,
    labels = upper_simper$taxon_label,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(
    labels = function(x) label_number(accuracy = 1)(abs(x)),
    limits = c(-upper_simper_max, upper_simper_max),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c(
      "Higher post-fire" = "firebrick",
      "Higher pre-fire" = "#1B4F72"
    ),
    guide = "none"
  ) +
  labs(
    x = "Contribution to dissimilarity (%)",
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  base_theme +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 4.7, face = "italic"),
    axis.text.x = element_text(size = 5.5),
    axis.title.x = element_text(size = 6, face = "bold"),
    plot.margin = margin(1, 2, 1, 1, "mm")
  )

upper_core <- wrap_elements(full = upper_panel_a) +
  wrap_elements(full = upper_panel_b) +
  plot_layout(widths = c(0.59, 0.41)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = FINAL_TAG_SIZE,
      family = BASE_FAMILY
    ),
    plot.tag.position = c(0.03, 0.99)
  )

upper_core_constant <- wrap_elements(full = upper_panel_a_constant) +
  wrap_elements(full = upper_panel_b) +
  plot_layout(widths = c(0.59, 0.41)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = FINAL_TAG_SIZE,
      family = BASE_FAMILY
    ),
    plot.tag.position = c(0.03, 0.99)
  )

# ---- Upper-only panel c ------------------------------------------------------
upper_background_maxima <- background_maxima %>%
  filter(depth_zone == "Upper")
upper_focal_values <- focal_values %>%
  filter(depth_zone == "Upper")
upper_focal_maxima <- focal_maxima %>%
  filter(depth_zone == "Upper")

upper_q95 <- as.numeric(
  quantile(
    upper_background_maxima$maximum_departure,
    probs = 0.95,
    na.rm = TRUE,
    type = 8
  )
)

upper_event_colors <- c(
  "Background non-fire years" = "#D8D8D8",
  "2011 Leptolyngbya" = "#63C7A6",
  "2020 smoke" = "#F3BD58",
  "2021 Caldor" = "#E67878"
)

upper_only_panel <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = upper_q95,
    fill = REFERENCE_FILL,
    alpha = REFERENCE_ALPHA
  ) +
  geom_hline(
    yintercept = upper_q95,
    colour = REFERENCE_LINE,
    linetype = "dashed",
    linewidth = 0.55 * LO_FIG_SCALE
  ) +
  annotate(
    "text",
    x = 1.08,
    y = upper_q95,
    label = "Historical 95th percentile",
    hjust = 0,
    vjust = -0.45,
    size = 1.75 * LO_FIG_SCALE,
    family = BASE_FAMILY,
    colour = REFERENCE_LINE
  ) +
  geom_violin(
    data = upper_background_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      fill = as.character(comparison_group)
    ),
    width = 0.55,
    scale = "width",
    trim = FALSE,
    linewidth = 0.45 * LO_FIG_SCALE,
    colour = "grey25",
    alpha = 0.78
  ) +
  geom_jitter(
    data = upper_background_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      fill = as.character(comparison_group)
    ),
    width = 0.055,
    height = 0,
    shape = 21,
    size = 1.2 * LO_FIG_SCALE,
    stroke = 0.2 * LO_FIG_SCALE,
    colour = "grey35",
    alpha = 0.38
  ) +
  geom_jitter(
    data = upper_focal_values,
    aes(
      x = comparison_group,
      y = departure,
      fill = as.character(comparison_group)
    ),
    width = 0.055,
    height = 0,
    shape = 21,
    size = 1.55 * LO_FIG_SCALE,
    stroke = 0.25 * LO_FIG_SCALE,
    colour = "grey20",
    alpha = 0.18
  ) +
  geom_point(
    data = upper_focal_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      fill = as.character(comparison_group)
    ),
    shape = 21,
    size = 3.35 * LO_FIG_SCALE,
    stroke = 0.75 * LO_FIG_SCALE,
    colour = "grey10"
  ) +
  geom_text(
    data = upper_focal_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      label = sprintf("%.2f", maximum_departure)
    ),
    vjust = -1.10,
    size = 1.85 * LO_FIG_SCALE,
    family = BASE_FAMILY,
    fontface = "bold",
    colour = "grey10"
  ) +
  scale_fill_manual(values = upper_event_colors, guide = "none") +
  scale_x_discrete(
    limits = GROUP_LEVELS,
    labels = c(
      "Background non-fire years" = "Background non-fire\nyears",
      "2011 Leptolyngbya" = "2011\nLeptolyngbya",
      "2020 smoke" = "2020\nsmoke exposure",
      "2021 Caldor" = "2021\nCaldor fire"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    tag = "c",
    title = NULL,
    x = NULL,
    y = "Maximum Bray-Curtis distance\nto seasonal background centroid"
  ) +
  theme_bw(base_size = 8, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      colour = "#EBEBEB",
      linewidth = LO_LW_THIN
    ),
    axis.text = element_text(size = LO_FS_AXIS_TEXT, colour = "grey20"),
    axis.title.y = element_text(size = LO_FS_AXIS_TITLE, face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.tag = element_text(size = FINAL_TAG_SIZE, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "none",
    plot.margin = margin(2, 3, 1, 2, "mm")
  )

upper_only_composite <- wrap_elements(full = upper_core) /
  upper_only_panel +
  plot_layout(heights = c(13.5, 6.5))

agg_png(UPPER_PANEL_FILE, width = 17.8, height = 6.5,
        units = "cm", res = LO_DPI)
print(upper_only_panel)
invisible(dev.off())

agg_png(UPPER_COMPOSITE_FILE, width = 17.8, height = 20.0,
        units = "cm", res = LO_DPI)
print(upper_only_composite)
invisible(dev.off())

# ---- Diagnostics -------------------------------------------------------------
cat("\nPooled full-width historical 95th percentile:",
    round(full_width_q95, 3), "\n")
cat("Upper-only historical 95th percentile:", round(upper_q95, 3), "\n")
cat("Saved final Upper + Lower composite:", COMPOSITE_FILE, "\n")
cat("Saved final Upper-only composite:   ", UPPER_COMPOSITE_FILE, "\n")
