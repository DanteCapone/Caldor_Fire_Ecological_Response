# =============================================================================
# figure_5_ltp_abundance_historical_context_depth.R
# Caldor Fire Ecosystem Response Project
#
# Purpose
#   Create an Upper/Lower-depth historical-context panel for the existing core
#   LTP abundance Figure 5. The comparison window is seven months beginning in
#   August (1 August through the end of February in the following calendar year).
#
# Design
#   * Background pre-fire non-fire years: paired Upper and Lower violin plots of
#     each window-year's maximum community departure.
#   * 2011 Leptolyngbya, 2020 smoke exposure, and 2021 Caldor: all observed
#     August-February distances appear with high transparency; each depth-zone
#     maximum is emphasized with a large outlined point and bold value label.
#   * Upper uses a circle and a lighter event shade; Lower uses a triangle and a
#     darker event shade.
#
# Metric
#   Matches scripts/figure_5_revised.R: Hellinger transformation, Bray-Curtis
#   dissimilarity represented in corrected PCoA space, and Euclidean distance
#   from each sample to its season-matched background centroid.
#
# Outputs
#   figures/figure_5_nmds_braycurtis/revised/
#     figure_5_ltp_abundance_historical_context_depth_panel.png
#     figure_5_ltp_abundance_historical_context_depth.png
#     figure_5_ltp_abundance_historical_context_depth_maxima.csv
#     figure_5_ltp_abundance_historical_context_depth_samples.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
  library(patchwork)
  library(ragg)
  library(png)
})

source("scripts/main_analysis/shared_aesthetics.R")

# ---- Paths ------------------------------------------------------------------
DATA_FILE <- file.path(
  "data", "lake_environmental_data", "phytoplankton",
  "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"
)
OUT_DIR <- file.path("figures", "figure_5_nmds_braycurtis", "archived",
                     "revised_analysis")
CORE_FIGURE <- file.path(OUT_DIR, "figure_5_ltp_abundance.png")
PANEL_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_depth_panel.png"
)
COMPOSITE_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_depth.png"
)
MAXIMA_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_depth_maxima.csv"
)
SAMPLE_FILE <- file.path(
  OUT_DIR, "figure_5_ltp_abundance_historical_context_depth_samples.csv"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(DATA_FILE)) stop("Historical dataset not found: ", DATA_FILE)
if (!file.exists(CORE_FIGURE)) stop("Core Figure 5 not found: ", CORE_FIGURE)

# ---- Analysis definitions ----------------------------------------------------
BASE_FAMILY <- "Times New Roman"
FOCAL_WINDOW_YEARS <- c(2011L, 2020L, 2021L)
EXCLUDED_BACKGROUND_WINDOW_YEARS <- c(2011L, 2014L, 2020L, 2021L)
MIN_REFERENCE_SAMPLES <- 2L

DEPTH_LOOKUP <- c(
  "0-40 m" = "Upper",
  "60-105 m" = "Lower"
)
DEPTH_LEVELS <- c("Upper", "Lower")
DEPTH_SHAPES <- c("Upper" = 21, "Lower" = 24)

GROUP_LEVELS <- c(
  "Background non-fire years",
  "2011 Leptolyngbya",
  "2020 smoke",
  "2021 Caldor"
)

# Lighter shades represent Upper; darker shades represent Lower.
EVENT_DEPTH_COLORS <- c(
  "Background non-fire years | Upper" = "#D8D8D8",
  "Background non-fire years | Lower" = "#777777",
  "2011 Leptolyngbya | Upper" = "#63C7A6",
  "2011 Leptolyngbya | Lower" = "#007A59",
  "2020 smoke | Upper" = "#F3BD58",
  "2020 smoke | Lower" = "#B96B00",
  "2021 Caldor | Upper" = "#E67878",
  "2021 Caldor | Lower" = "#981F1F"
)

# Seven-month comparison window: August through February, inclusive.
window_year_from_date <- function(dates) {
  case_when(
    month(dates) >= 8L ~ year(dates),
    month(dates) <= 2L ~ year(dates) - 1L,
    TRUE ~ NA_integer_
  )
}

get_season <- function(dates) {
  month_number <- month(dates)
  case_when(
    month_number %in% c(12, 1, 2)  ~ "Winter",
    month_number %in% c(3, 4, 5)   ~ "Spring",
    month_number %in% c(6, 7, 8)   ~ "Summer",
    month_number %in% c(9, 10, 11) ~ "Fall",
    TRUE ~ NA_character_
  )
}

# ---- Load the long-term LTP community record --------------------------------
# Data through February 2022 are required to complete the 2021 comparison
# window. Later recovery samples are outside the requested seven-month window.
community_long <- read_csv(DATA_FILE, show_col_types = FALSE) %>%
  transmute(
    station_id = as.character(station_id),
    date = as.Date(date),
    calendar_year = year(date),
    window_year = window_year_from_date(date),
    season = get_season(date),
    depth_original = as.character(depth_bin),
    depth_zone = unname(DEPTH_LOOKUP[depth_original]),
    taxon = str_squish(as.character(taxon)),
    abundance = as.numeric(abundance)
  ) %>%
  filter(
    station_id == "LTP",
    date <= as.Date("2022-02-28"),
    !is.na(date),
    !is.na(season),
    !is.na(depth_zone),
    !is.na(taxon),
    taxon != "",
    is.finite(abundance),
    abundance >= 0
  ) %>%
  mutate(depth_zone = factor(depth_zone, levels = DEPTH_LEVELS)) %>%
  group_by(
    station_id, depth_zone, date, calendar_year, window_year, season, taxon
  ) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Background is defined by August-February window year, not calendar year. This
# prevents Jan-Feb 2012 (part of the 2011 Leptolyngbya window), Jan-Feb 2015
# (part of the excluded 2014 fire window), and Jan-Feb 2021 (part of the 2020
# smoke window) from leaking into the reference state.
window_coverage <- community_long %>%
  filter(!is.na(window_year)) %>%
  distinct(window_year, depth_zone, date) %>%
  mutate(calendar_month = month(date)) %>%
  group_by(window_year, depth_zone) %>%
  summarise(
    has_august_to_december = any(calendar_month >= 8L),
    has_january_to_february = any(calendar_month <= 2L),
    .groups = "drop"
  )

# A background window is retained only when both depth zones contain samples in
# both calendar-year portions of the seven-month interval. This removes partial
# edge/gap windows that would otherwise have fewer opportunities to produce a
# maximum than the focal years.
background_window_years <- window_coverage %>%
  group_by(window_year) %>%
  summarise(
    n_depth_zones = n_distinct(depth_zone),
    complete_window =
      all(has_august_to_december) & all(has_january_to_february),
    .groups = "drop"
  ) %>%
  filter(
    window_year < 2021L,
    !window_year %in% EXCLUDED_BACKGROUND_WINDOW_YEARS,
    n_depth_zones == length(DEPTH_LEVELS),
    complete_window
  ) %>%
  arrange(window_year) %>%
  pull(window_year)

if (length(background_window_years) < 3L) {
  stop("Fewer than three usable background August-February windows were found.")
}

# ---- Distance-to-seasonal-centroid calculation by depth zone ----------------
calculate_depth_departure <- function(depth_data) {
  depth_name <- as.character(first(depth_data$depth_zone))

  community_wide <- depth_data %>%
    select(date, calendar_year, window_year, season, taxon, abundance) %>%
    pivot_wider(
      names_from = taxon,
      values_from = abundance,
      values_fill = 0,
      values_fn = sum
    ) %>%
    arrange(date)

  metadata <- community_wide %>%
    select(date, calendar_year, window_year, season)
  abundance_matrix <- community_wide %>%
    select(-date, -calendar_year, -window_year, -season) %>%
    as.matrix()
  storage.mode(abundance_matrix) <- "double"

  keep_samples <- rowSums(abundance_matrix, na.rm = TRUE) > 0
  metadata <- metadata[keep_samples, , drop = FALSE]
  abundance_matrix <- abundance_matrix[keep_samples, , drop = FALSE]

  keep_taxa <- colSums(abundance_matrix, na.rm = TRUE) > 0
  abundance_matrix <- abundance_matrix[, keep_taxa, drop = FALSE]

  if (nrow(abundance_matrix) < 6L || ncol(abundance_matrix) < 2L) {
    stop("Insufficient community data for ", depth_name)
  }

  hellinger_matrix <- decostand(abundance_matrix, method = "hellinger")
  bray_matrix <- vegdist(hellinger_matrix, method = "bray")
  pcoa <- cmdscale(
    bray_matrix,
    k = nrow(hellinger_matrix) - 1L,
    eig = TRUE,
    add = TRUE
  )

  eigenvalues <- pcoa$eig[seq_len(ncol(pcoa$points))]
  coordinates <- pcoa$points[, eigenvalues > 1e-12, drop = FALSE]
  if (ncol(coordinates) == 0L) stop("No positive PCoA axes for ", depth_name)

  map_dfr(seq_len(nrow(metadata)), function(sample_index) {
    sample_season <- metadata$season[sample_index]

    reference_indices <- which(
      metadata$window_year %in% background_window_years &
        metadata$season == sample_season
    )
    # Leave a background sample out of its own seasonal centroid.
    reference_indices <- setdiff(reference_indices, sample_index)

    if (length(reference_indices) < MIN_REFERENCE_SAMPLES) {
      return(tibble(
        date = metadata$date[sample_index],
        calendar_year = metadata$calendar_year[sample_index],
        window_year = metadata$window_year[sample_index],
        season = sample_season,
        depth_zone = depth_name,
        n_reference = length(reference_indices),
        departure = NA_real_
      ))
    }

    seasonal_centroid <- colMeans(
      coordinates[reference_indices, , drop = FALSE]
    )

    tibble(
      date = metadata$date[sample_index],
      calendar_year = metadata$calendar_year[sample_index],
      window_year = metadata$window_year[sample_index],
      season = sample_season,
      depth_zone = depth_name,
      n_reference = length(reference_indices),
      departure = sqrt(sum((coordinates[sample_index, ] - seasonal_centroid)^2))
    )
  })
}

sample_departure <- community_long %>%
  group_by(depth_zone) %>%
  group_split(.keep = TRUE) %>%
  map_dfr(calculate_depth_departure) %>%
  filter(
    is.finite(departure),
    !is.na(window_year),
    window_year %in% c(background_window_years, FOCAL_WINDOW_YEARS)
  ) %>%
  mutate(
    depth_zone = factor(depth_zone, levels = DEPTH_LEVELS),
    comparison_group = case_when(
      window_year == 2011L ~ "2011 Leptolyngbya",
      window_year == 2020L ~ "2020 smoke",
      window_year == 2021L ~ "2021 Caldor",
      window_year %in% background_window_years ~ "Background non-fire years",
      TRUE ~ NA_character_
    ),
    comparison_group = factor(comparison_group, levels = GROUP_LEVELS),
    event_depth = paste(comparison_group, depth_zone, sep = " | ")
  ) %>%
  filter(!is.na(comparison_group))

# ---- Maximum departure for every window year and depth zone -----------------
depth_maxima <- sample_departure %>%
  group_by(window_year, comparison_group, depth_zone, event_depth) %>%
  slice_max(departure, n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    window_year,
    comparison_group,
    depth_zone,
    event_depth,
    maximum_departure = departure,
    date_of_maximum = date,
    season_of_maximum = season,
    n_reference
  ) %>%
  arrange(window_year, depth_zone)

background_maxima <- depth_maxima %>%
  filter(comparison_group == "Background non-fire years")
focal_maxima <- depth_maxima %>%
  filter(window_year %in% FOCAL_WINDOW_YEARS)
focal_values <- sample_departure %>%
  filter(window_year %in% FOCAL_WINDOW_YEARS)

# Retain a separate historical 95th-percentile reference for the Upper and
# Lower background distributions. These are shown as short dashed lines across
# the paired background violins rather than a single pooled threshold.
background_q95_by_depth <- background_maxima %>%
  group_by(depth_zone) %>%
  summarise(
    q95 = quantile(maximum_departure, 0.95, na.rm = TRUE, type = 8),
    .groups = "drop"
  ) %>%
  mutate(
    x_center = if_else(depth_zone == "Upper", 0.875, 1.125),
    x_start = x_center - 0.16,
    x_end = x_center + 0.16
  )

# Confirm that all six emphasized points are the true depth-specific maxima.
maximum_audit <- focal_values %>%
  group_by(window_year, depth_zone) %>%
  summarise(recalculated_maximum = max(departure), .groups = "drop") %>%
  left_join(
    focal_maxima %>%
      select(window_year, depth_zone, plotted_maximum = maximum_departure),
    by = c("window_year", "depth_zone")
  ) %>%
  mutate(matches = near(recalculated_maximum, plotted_maximum))

if (nrow(maximum_audit) != 6L || !all(maximum_audit$matches)) {
  stop("One or more emphasized Upper/Lower points are not true window maxima.")
}

write_csv(depth_maxima, MAXIMA_FILE)
write_csv(sample_departure, SAMPLE_FILE)

# ---- Grouped violin and focal-observation panel ------------------------------
# Positioning is shared across layers so Upper and Lower remain paired within
# each x-axis category.
depth_dodge <- position_dodge(width = 0.50)
raw_jitter_dodge <- position_jitterdodge(
  jitter.width = 0.055,
  jitter.height = 0,
  dodge.width = 0.50,
  seed = 42
)

historical_depth_panel <- ggplot() +
  # Separate background violins show the distribution of annual maxima for the
  # Upper and Lower zones while sharing one background x-axis category.
  geom_violin(
    data = background_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = interaction(comparison_group, depth_zone),
      fill = event_depth
    ),
    position = depth_dodge,
    width = 0.82,
    scale = "width",
    trim = FALSE,
    linewidth = 0.45 * LO_FIG_SCALE,
    colour = "grey25",
    alpha = 0.78
  ) +
  # Depth-specific historical 95th-percentile reference lines.
  geom_segment(
    data = background_q95_by_depth,
    aes(x = x_start, xend = x_end, y = q95, yend = q95),
    inherit.aes = FALSE,
    colour = "grey20",
    linetype = "dashed",
    linewidth = 0.5 * LO_FIG_SCALE
  ) +
  # Faint background maxima reveal the sample support behind each violin.
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
  # All focal observations remain visible but highly transparent.
  geom_point(
    data = focal_values,
    aes(
      x = comparison_group,
      y = departure,
      group = depth_zone,
      shape = depth_zone,
      fill = event_depth
    ),
    position = raw_jitter_dodge,
    size = 1.55 * LO_FIG_SCALE,
    stroke = 0.25 * LO_FIG_SCALE,
    colour = "grey20",
    alpha = 0.18
  ) +
  # The true Upper and Lower maxima are large, opaque, and outlined.
  geom_point(
    data = focal_maxima,
    aes(
      x = comparison_group,
      y = maximum_departure,
      group = depth_zone,
      shape = depth_zone,
      fill = event_depth
    ),
    position = depth_dodge,
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
    position = depth_dodge,
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
    title = "Historical context: maximum August-February community departure",
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
    plot.title = element_text(size = LO_FS_TITLE, face = "bold"),
    plot.subtitle = element_blank(),
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_text(size = 6.1 * LO_FIG_SCALE, face = "bold"),
    legend.text = element_text(size = 5.8 * LO_FIG_SCALE),
    legend.key.width = unit(4.5 * LO_FIG_SCALE, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 3, 2, 2, "mm")
  )

# ---- Export the standalone panel and unchanged-core composite ----------------
agg_png(PANEL_FILE, width = 17.8, height = 6.5, units = "cm", res = LO_DPI)
print(historical_depth_panel)
invisible(dev.off())

core_raster <- png::readPNG(CORE_FIGURE)
core_panel <- wrap_elements(
  full = grid::rasterGrob(core_raster, interpolate = TRUE)
)

figure_5_ltp_abundance_historical_depth <- core_panel /
  historical_depth_panel +
  plot_layout(heights = c(13.5, 6.5))

agg_png(COMPOSITE_FILE, width = 17.8, height = 20.0,
        units = "cm", res = LO_DPI)
print(figure_5_ltp_abundance_historical_depth)
invisible(dev.off())

# ---- Execution diagnostics --------------------------------------------------
cat("\nBackground August-February window years:\n")
cat(paste(background_window_years, collapse = ", "), "\n")
cat("\nVerified focal maxima by window year and depth zone:\n")
print(
  focal_maxima %>%
    select(
      window_year, depth_zone, maximum_departure,
      date_of_maximum, season_of_maximum
    ),
  n = Inf
)
cat("\nMaximum audit:\n")
print(maximum_audit, n = Inf)
cat("\nSaved panel:    ", PANEL_FILE, "\n", sep = "")
cat("Saved composite:", COMPOSITE_FILE, "\n", sep = "")
