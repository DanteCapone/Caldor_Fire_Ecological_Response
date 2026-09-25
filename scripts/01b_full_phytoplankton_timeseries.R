# =============================================================================
# 01b_full_phytoplankton_timeseries.R
# Caldor Fire Ecosystem Response Project
#
# Purpose: Prepare and plot the complete LTP phytoplankton time series from the
#          curated 2005-2025 workbook. Depths with fewer than five distinct
#          sampling dates are excluded from every output.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(scales)
  library(ragg)
})
source("scripts/figure_aesthetics.R")

PHYTO_FILE <- file.path(
  "data/lake_environmental_data/phytoplankton",
  "PhytoData_LTP_2005-2025_counts_biovolume_size_v1.xlsx"
)
PHYTO_SHEET <- "LTP_2005-2025"
OUT_DIR <- "figures/phytoplankton/full_time_series"
DATA_EXPORT <- file.path(
  dirname(PHYTO_FILE),
  "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CALDOR_START <- as.Date("2021-08-14")
CALDOR_END <- as.Date("2021-10-21")
SERIES_START <- as.Date("2005-01-01")
SERIES_END <- as.Date("2025-12-31")
MIN_TIME_POINTS <- 5L
TOP_N <- 10L
SIZE_COL <- "Life form size classes by longest linear dimension (LLD)"
SIZE_LEVELS <- c("2-20 µm", ">20 µm")
SIZE_COLS <- c("2-20 µm" = "#56B4E9", ">20 µm" = "#D55E00")

parse_excel_date <- function(x) {
  x <- trimws(as.character(x))
  d_iso <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  n <- suppressWarnings(as.numeric(x))
  d_serial <- as.Date(n, origin = as.Date("1899-12-30"))
  coalesce(d_iso, d_serial)
}

harmonize_phyto_taxon <- function(x) {
  cleaned <- str_squish(x)
  cleaned <- str_remove(cleaned, regex("^cf\\.\\s*", ignore_case = TRUE))
  original_genus <- str_extract(cleaned, "^[A-Za-z-]+")
  corrected_genus <- recode(
    original_genus,
    "Kephrion" = "Kephyrion",
    "Planktonema" = "Planctonema",
    "Scendesmus" = "Scenedesmus",
    "Peduastrum" = "Pediastrum",
    "Staurasturm" = "Staurastrum",
    .default = original_genus
  )
  remainder <- str_remove(cleaned, "^[A-Za-z-]+\\s*")
  normalized <- if_else(
    !is.na(corrected_genus) & corrected_genus != original_genus,
    str_c(corrected_genus, if_else(remainder == "", "", " "), remainder),
    cleaned
  )
  generic_prefixes <- c(
    "Unidentified", "Unknown", "Flagellate", "Flagellates", "Coccoid",
    "Colonial", "Non-motile", "Motile", "Small", "Large", "Centric",
    "Pennate", "Diatom", "Diatoms", "Chlorococcales", "Chrysophyte",
    "Cryptophyte", "Dinoflagellate", "Green", "Blue-green", "Other", "Cyst"
  )
  has_taxonomic_epithet <- str_detect(normalized, "^[A-Z][A-Za-z-]+\\s+\\S+")
  genus_level <- has_taxonomic_epithet &
    !is.na(corrected_genus) &
    !corrected_genus %in% generic_prefixes
  if_else(genus_level, str_c(corrected_genus, " spp."), normalized)
}

if (!file.exists(PHYTO_FILE)) stop("Missing phytoplankton workbook: ", PHYTO_FILE)

# Date is deliberately read as text because 2005-2024 are Excel serial dates,
# while 2025 is stored as ISO text in the same column.
raw <- read_excel(
  PHYTO_FILE,
  sheet = PHYTO_SHEET,
  col_types = c(
    "text", "text", "text", "text", "text", "numeric", "text", "text",
    "numeric", "text", "numeric", "text", "text"
  )
)

phyto_all <- raw %>%
  transmute(
    station_id = .data[["Station_ID"]],
    event_id = .data[["Event_ID"]],
    date = parse_excel_date(.data[["Date"]]),
    depth_num = as.numeric(.data[["Depth_(m)"]]),
    phylum = .data[["Phylum"]],
    taxon_original = str_squish(.data[["Taxon"]]),
    taxon = harmonize_phyto_taxon(taxon_original),
    abundance = as.numeric(.data[["Abundance_(units/L)"]]),
    unit_biovolume_um3 = as.numeric(.data[["Unit_biovolume_(µm3)"]]),
    # Equivalent to mm3 m-3: (um3 L-1) / 1e6.
    biovolume = abundance * unit_biovolume_um3 / 1e6,
    size_class = str_replace_all(.data[[SIZE_COL]], "μ", "µ"),
    counting_unit = .data[["Counting_unit"]],
    notes = .data[["Notes"]]
  ) %>%
  mutate(
    size_class = str_squish(size_class),
    size_class = if_else(size_class %in% SIZE_LEVELS, size_class, NA_character_),
    year = year(date),
    month = month(date),
    season = factor(
      case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% 3:5 ~ "Spring",
        month %in% 6:8 ~ "Summer",
        month %in% 9:11 ~ "Fall"
      ),
      levels = c("Winter", "Spring", "Summer", "Fall")
    ),
    fire_period = if_else(date < CALDOR_START, "Pre-fire", "Post-fire")
  ) %>%
  filter(
    station_id == "LTP",
    !is.na(date),
    date >= SERIES_START,
    date <= SERIES_END,
    !is.na(depth_num),
    !is.na(taxon),
    taxon != "",
    is.finite(abundance),
    abundance >= 0
  )

depth_coverage <- phyto_all %>%
  distinct(depth_num, date) %>%
  count(depth_num, name = "n_time_points") %>%
  arrange(depth_num) %>%
  mutate(
    minimum_required = MIN_TIME_POINTS,
    included = n_time_points >= MIN_TIME_POINTS
  )
write_csv(depth_coverage, file.path(OUT_DIR, "ltp_depth_time_point_coverage.csv"))

eligible_depths <- depth_coverage %>% filter(included) %>% pull(depth_num)
excluded_depths <- depth_coverage %>% filter(!included) %>% pull(depth_num)
depth_levels <- paste0(eligible_depths, " m")
bin_levels <- c("0-40 m", "60-105 m")

phyto <- phyto_all %>%
  filter(depth_num %in% eligible_depths) %>%
  mutate(
    depth_category = factor(paste0(depth_num, " m"), levels = depth_levels),
    depth_bin = case_when(
      depth_num >= 0 & depth_num <= 40 ~ "0-40 m",
      depth_num >= 60 & depth_num <= 105 ~ "60-105 m",
      TRUE ~ NA_character_
    ),
    depth_bin = factor(depth_bin, levels = bin_levels),
    depth_m = as.character(depth_num),
    depth_sort = depth_num,
    sample_id = paste0(station_id, "_", format(date, "%Y-%m-%d"), "_", depth_m)
  )

taxon_crosswalk <- phyto %>%
  group_by(taxon_original, taxon) %>%
  summarise(
    Phylum = paste(sort(unique(na.omit(phylum))), collapse = " | "),
    first_date = min(date),
    last_date = max(date),
    n_dates = n_distinct(date),
    total_abundance = sum(abundance, na.rm = TRUE),
    present_2005_2018 = any(year <= 2018),
    present_2019_2025 = any(year >= 2019),
    .groups = "drop"
  ) %>%
  arrange(taxon, taxon_original)
write_csv(taxon_crosswalk, file.path(OUT_DIR, "ltp_taxon_harmonization_crosswalk.csv"))

phyto_export <- phyto %>%
  select(
    station_id, date, depth_m, depth_num, depth_sort, depth_category, depth_bin,
    taxon_original, taxon, abundance, biovolume, size_class, year, month,
    season, fire_period, sample_id, event_id, phylum, counting_unit, notes
  )
write_csv(phyto_export, DATA_EXPORT)

full_ts_top_taxa <- phyto %>%
  group_by(taxon) %>%
  summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total, n = TOP_N, with_ties = FALSE) %>%
  arrange(desc(total)) %>%
  pull(taxon)

all_taxa <- sort(unique(phyto$taxon))
all_taxon_cols <- setNames(hcl.colors(length(all_taxa), palette = "Dark 3"), all_taxa)
overrides <- c(
  "Leptolyngbya spp." = "#009E73",
  "Cryptomonas spp." = "#9467BD",
  "Cyclotella spp." = "#0072B2",
  "Synedra spp." = "#D55E00"
)
all_taxon_cols[intersect(names(all_taxon_cols), names(overrides))] <-
  overrides[intersect(names(all_taxon_cols), names(overrides))]
taxon_cols <- c(all_taxon_cols[full_ts_top_taxa], Other = "#BBBBBB")

full_ts_theme <- theme_classic(base_size = 7, base_family = LO_FONT) +
  theme(
    text = element_text(family = LO_FONT),
    axis.title = element_text(size = 7, family = LO_FONT),
    axis.text = element_text(size = 6, family = LO_FONT),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 6, family = LO_FONT),
    legend.text = element_text(size = 6, family = LO_FONT),
    legend.key.size = unit(1.8, "mm"),
    panel.grid.major.y = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold")
  )

make_taxon_plot <- function(metric_col, metric_label, group_col, relative = FALSE) {
  group_levels <- if (group_col == "depth_category") depth_levels else bin_levels
  dat <- phyto %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[metric_col]])) %>%
    mutate(
      plot_group = factor(as.character(.data[[group_col]]), levels = group_levels),
      taxon_group = if_else(taxon %in% full_ts_top_taxa, taxon, "Other")
    ) %>%
    group_by(plot_group, date, taxon_group) %>%
    summarise(value = sum(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(plot_group, date) %>%
    mutate(
      total_value = sum(value, na.rm = TRUE),
      plot_value = if (relative) if_else(total_value > 0, 100 * value / total_value, 0) else value
    ) %>%
    ungroup() %>%
    mutate(taxon_group = factor(taxon_group, levels = rev(c(full_ts_top_taxa, "Other"))))

  y_scale <- if (relative) {
    scale_y_continuous(limits = c(0, 100), labels = label_number(accuracy = 1), expand = c(0, 0))
  } else {
    scale_y_continuous(labels = label_scientific(), expand = expansion(mult = c(0, 0.12)))
  }
  y_label <- if (relative) {
    if (metric_col == "abundance") "Relative abundance (%)" else "Relative biovolume (%)"
  } else metric_label

  p <- ggplot(dat, aes(date, plot_value, fill = taxon_group)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_col(width = 29, colour = NA) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = LO_LW_MID) +
    facet_grid(rows = vars(plot_group), scales = if (relative) "fixed" else "free_y", switch = "y") +
    scale_fill_manual(
      values = taxon_cols,
      breaks = c(full_ts_top_taxa, "Other"),
      labels = lo_abbreviate_taxon,
      name = "Taxon",
      guide = guide_legend(ncol = 2, reverse = TRUE, byrow = TRUE)
    ) +
    scale_x_date(
      limits = c(SERIES_START, SERIES_END),
      date_breaks = "2 years", date_labels = "%Y",
      expand = expansion(mult = c(0.002, 0.002))
    ) +
    y_scale +
    labs(x = NULL, y = y_label) +
    full_ts_theme

  metric_stub <- if_else(metric_col == "abundance", "abundance", "biovolume")
  group_stub <- if_else(group_col == "depth_bin", "by_figure4_depth_bins", "by_depth")
  value_stub <- if_else(relative, "relative", "absolute")
  out_path <- file.path(
    OUT_DIR,
    paste0("full_time_series_", metric_stub, "_", value_stub, "_", group_stub, ".png")
  )
  fig_height <- if (group_col == "depth_category") max(18, 3.25 * length(depth_levels)) else 14
  ggsave(out_path, p, width = 24, height = fig_height, units = "cm",
         dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
  message("Saved: ", out_path)
}

for (metric_col in c("abundance", "biovolume")) {
  metric_label <- if (metric_col == "abundance") {
    expression("Abundance (cells L"^-1*")")
  } else {
    expression("Biovolume (mm"^3*" m"^-3*")")
  }
  for (group_col in c("depth_category", "depth_bin")) {
    make_taxon_plot(metric_col, metric_label, group_col, relative = FALSE)
    make_taxon_plot(metric_col, metric_label, group_col, relative = TRUE)
  }
}

make_size_plot <- function(group_col) {
  group_levels <- if (group_col == "depth_category") depth_levels else bin_levels
  classified <- phyto %>%
    filter(!is.na(.data[[group_col]]), !is.na(size_class)) %>%
    mutate(plot_group = factor(as.character(.data[[group_col]]), levels = group_levels)) %>%
    group_by(plot_group, date, size_class) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(plot_group, date) %>%
    mutate(
      classified_abundance = sum(abundance, na.rm = TRUE),
      percent = if_else(classified_abundance > 0, 100 * abundance / classified_abundance, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = SIZE_LEVELS))

  ratio <- classified %>%
    select(plot_group, date, size_class, abundance, percent) %>%
    pivot_wider(names_from = size_class, values_from = c(abundance, percent), values_fill = 0) %>%
    mutate(
      large_to_small_abundance_ratio = if_else(
        .data[["abundance_2-20 µm"]] > 0,
        .data[["abundance_>20 µm"]] / .data[["abundance_2-20 µm"]],
        NA_real_
      ),
      group_type = group_col
    )

  p <- ggplot(classified, aes(date, percent, fill = size_class)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_col(width = 29, colour = NA) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = LO_LW_MID) +
    facet_grid(rows = vars(plot_group), switch = "y") +
    scale_fill_manual(values = SIZE_COLS, breaks = SIZE_LEVELS,
                      name = "LLD size class") +
    scale_x_date(
      limits = c(SERIES_START, SERIES_END), date_breaks = "2 years",
      date_labels = "%Y", expand = expansion(mult = c(0.002, 0.002))
    ) +
    scale_y_continuous(limits = c(0, 100), labels = label_percent(scale = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Classified abundance (%)") +
    full_ts_theme

  stub <- if_else(group_col == "depth_category", "by_depth", "by_figure4_depth_bins")
  out_path <- file.path(OUT_DIR, paste0("full_time_series_size_class_relative_", stub, ".png"))
  fig_height <- if (group_col == "depth_category") max(18, 3.25 * length(depth_levels)) else 14
  ggsave(out_path, p, width = 24, height = fig_height, units = "cm",
         dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
  message("Saved: ", out_path)
  ratio
}

size_ratio_by_depth <- make_size_plot("depth_category")
size_ratio_by_bin <- make_size_plot("depth_bin")
write_csv(
  bind_rows(size_ratio_by_depth, size_ratio_by_bin),
  file.path(OUT_DIR, "full_time_series_size_class_ratios.csv")
)

size_class_coverage <- phyto %>%
  summarise(
    total_rows = n(),
    classified_rows = sum(!is.na(size_class)),
    total_abundance = sum(abundance, na.rm = TRUE),
    classified_abundance = sum(abundance[!is.na(size_class)], na.rm = TRUE),
    classified_row_percent = 100 * classified_rows / total_rows,
    classified_abundance_percent = 100 * classified_abundance / total_abundance
  )
write_csv(size_class_coverage, file.path(OUT_DIR, "ltp_size_class_coverage.csv"))

make_leptolyngbya_plot <- function(group_col) {
  group_levels <- if (group_col == "depth_category") depth_levels else bin_levels
  sample_dates <- phyto %>%
    filter(!is.na(.data[[group_col]])) %>%
    transmute(plot_group = as.character(.data[[group_col]]), date) %>%
    distinct()
  dat <- phyto %>%
    filter(taxon == "Leptolyngbya spp.", !is.na(.data[[group_col]])) %>%
    transmute(plot_group = as.character(.data[[group_col]]), date, abundance) %>%
    group_by(plot_group, date) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    right_join(sample_dates, by = c("plot_group", "date")) %>%
    mutate(
      abundance = replace_na(abundance, 0),
      plot_group = factor(plot_group, levels = group_levels)
    )
  p <- ggplot(dat, aes(date, abundance)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_line(colour = "#009E73", linewidth = LO_LW_MID, alpha = 0.8) +
    geom_point(colour = "#009E73", size = 1, alpha = 0.85) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = LO_LW_MID) +
    facet_grid(rows = vars(plot_group), scales = "free_y", switch = "y") +
    scale_x_date(limits = c(SERIES_START, SERIES_END), date_breaks = "2 years",
                 date_labels = "%Y", expand = expansion(mult = c(0.002, 0.002))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10), labels = label_scientific(),
                       expand = expansion(mult = c(0.03, 0.12))) +
    labs(x = NULL, y = expression("Leptolyngbya spp. abundance (cells L"^-1*")")) +
    full_ts_theme + theme(legend.position = "none")
  stub <- if_else(group_col == "depth_category", "by_depth", "by_depth_bin")
  out_path <- file.path(OUT_DIR, paste0("full_time_series_leptolyngbya_abundance_", stub, ".png"))
  fig_height <- if (group_col == "depth_category") max(18, 3.25 * length(depth_levels)) else 13
  ggsave(out_path, p, width = 24, height = fig_height, units = "cm",
         dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
  message("Saved: ", out_path)
}

make_leptolyngbya_plot("depth_category")
make_leptolyngbya_plot("depth_bin")

cat("\nInput date range:", as.character(min(phyto$date)), "to", as.character(max(phyto$date)), "\n")
cat("Included depths (>=", MIN_TIME_POINTS, " dates):", paste(eligible_depths, collapse = ", "), "m\n")
cat("Excluded depths (<", MIN_TIME_POINTS, " dates):", paste(excluded_depths, collapse = ", "), "m\n")
cat("Analysis export:", DATA_EXPORT, "\n")
cat("All full time-series outputs written to:", OUT_DIR, "\n")
