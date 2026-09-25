# =============================================================================
# nutrients_2020_historical_comparison.R
#
# Compare August-October 2020 atmospheric deposition and surface-water
# nutrients with same-calendar-month observations from all other available
# years. Analytical definitions intentionally match Figure 3 v2.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")

FOCAL_YEAR <- 2020L
SEASON_MONTHS <- 8:10
MAX_DEPTH_M <- 10
N_ATOMIC_MASS <- 14.0067
P_ATOMIC_MASS <- 30.973762

NUTRI_MLTP <- file.path(
  "data", "lake_environmental_data", "nutrients", "Tahoe_MLTP_Nutrient.csv"
)
DEPO_FILE <- file.path(
  "data", "lake_environmental_data", "deposition",
  "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv"
)
OUT_DIR <- Sys.getenv(
  "CALDOR_2020_OUT_DIR",
  unset = file.path("figures", "supplemental", "nutrients_2020_historical")
)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NUT_LEVELS_DEP <- c("NO3", "NH4", "TKN", "SRP", "TP")
NUT_LEVELS_LAKE <- c("NO3", "NH4", "TKN", "TRP", "THP")
METRIC_LEVELS <- c("NO3", "NH4", "TKN", "SRP", "TP", "TRP", "THP",
                   "DIN:SRP", "DIN:TRP")
METRIC_DISPLAY <- c(
  NO3 = "NO3", NH4 = "NH4", TKN = "TKN", SRP = "SRP", TP = "TP",
  TRP = "TRP", THP = "THP", "DIN:SRP" = "DIN:SRP", "DIN:TRP" = "DIN:TRP"
)
COL_HIST <- "#BDBDBD"
COL_DEP_2020 <- "#D55E00"
COL_LAKE_2020 <- "#0072B2"

safe_log2 <- function(x) log2(pmax(x, .Machine$double.eps))

# Figure 3 v2 treats individual shallow observations as the analysis units.
# Preserve that convention here so the 2020 comparison is directly comparable.
lake_long <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Month %in% SEASON_MONTHS) %>%
  select(Date, Year, Month, NO3, NH4, TRP, THP, TKN) %>%
  pivot_longer(c(NO3, NH4, TRP, THP, TKN),
               names_to = "metric", values_to = "value") %>%
  mutate(
    source = "In-lake",
    period = if_else(Year == FOCAL_YEAR, "2020", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

dep_raw <- read_csv(DEPO_FILE, show_col_types = FALSE) %>%
  filter(!is.na(Start_Datetime), !is.na(End_Datetime)) %>%
  mutate(
    Start_Date = as.Date(Start_Datetime),
    End_Date = as.Date(End_Datetime),
    Date = Start_Date + as.integer(End_Date - Start_Date) %/% 2L,
    Year = year(Date), Month = month(Date)
  )

dep_long <- dep_raw %>%
  filter(Month %in% SEASON_MONTHS) %>%
  select(Date, Year, Month,
         NO3 = NO3_Daily_Load, NH4 = NH4_Daily_Load,
         TKN = TKN_Daily_Load, SRP = SRP_Daily_Load, TP = TP_Daily_Load) %>%
  pivot_longer(c(NO3, NH4, TKN, SRP, TP),
               names_to = "metric", values_to = "value") %>%
  mutate(
    source = "Deposition",
    period = if_else(Year == FOCAL_YEAR, "2020", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

# Ratio calculations reproduce Figure 3 v2 exactly. In particular, the
# in-lake DIN:TRP expression is retained as implemented there rather than
# silently changing an established manuscript-facing definition.
lake_ratio <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Month %in% SEASON_MONTHS) %>%
  transmute(
    Date, Year, Month,
    value = (as.numeric(NO3) + as.numeric(NH4)) /
      (as.numeric(TRP) / P_ATOMIC_MASS),
    metric = "DIN:TRP", source = "In-lake",
    period = if_else(Year == FOCAL_YEAR, "2020", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

dep_ratio <- dep_raw %>%
  filter(Month %in% SEASON_MONTHS) %>%
  transmute(
    Date, Year, Month,
    value = ((as.numeric(NO3_Daily_Load) + as.numeric(NH4_Daily_Load)) /
               N_ATOMIC_MASS) /
      (as.numeric(SRP_Daily_Load) / P_ATOMIC_MASS),
    metric = "DIN:SRP", source = "Deposition",
    period = if_else(Year == FOCAL_YEAR, "2020", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

analysis_input <- bind_rows(lake_long, dep_long, lake_ratio, dep_ratio) %>%
  mutate(
    metric = factor(metric, levels = METRIC_LEVELS),
    input_file = if_else(source == "Deposition", DEPO_FILE, NUTRI_MLTP),
    value_units = case_when(
      source == "Deposition" & metric %in% c("NO3", "NH4", "TKN") ~
        "mg N m^-2 d^-1",
      source == "Deposition" & metric %in% c("SRP", "TP") ~
        "mg P m^-2 d^-1",
      source == "In-lake" & metric %in% c("NO3", "NH4", "TKN") ~
        "ug N L^-1",
      source == "In-lake" & metric %in% c("TRP", "THP") ~
        "ug P L^-1",
      TRUE ~ "dimensionless ratio"
    ),
    analytical_definition = case_when(
      source == "Deposition" & metric == "DIN:SRP" ~
        "((NO3_Daily_Load + NH4_Daily_Load) / 14.0067) / (SRP_Daily_Load / 30.973762)",
      source == "In-lake" & metric == "DIN:TRP" ~
        "(NO3 + NH4) / (TRP / 30.973762), retained from Figure 3 v2",
      TRUE ~ "Positive finite source-column value; individual observations retained"
    )
  )

expected_metrics <- tribble(
  ~source, ~metric,
  "Deposition", "NO3", "Deposition", "NH4", "Deposition", "TKN",
  "Deposition", "SRP", "Deposition", "TP", "Deposition", "DIN:SRP",
  "In-lake", "NO3", "In-lake", "NH4", "In-lake", "TKN",
  "In-lake", "TRP", "In-lake", "THP", "In-lake", "DIN:TRP"
)
observed_metrics <- analysis_input %>%
  distinct(source, metric) %>%
  mutate(metric = as.character(metric))
if (nrow(anti_join(expected_metrics, observed_metrics, by = c("source", "metric")))) {
  stop("One or more expected nutrient/source combinations are absent.")
}

month_reference <- analysis_input %>%
  filter(period == "Historical") %>%
  group_by(source, metric, Month) %>%
  summarise(
    historical_month_median = median(value, na.rm = TRUE),
    n_historical_month = n(),
    .groups = "drop"
  )
if (any(!is.finite(month_reference$historical_month_median) |
        month_reference$historical_month_median <= 0)) {
  stop("A historical calendar-month median is missing, non-finite, or non-positive.")
}

plot_data <- analysis_input %>%
  left_join(month_reference, by = c("source", "metric", "Month")) %>%
  mutate(log_fold_difference = safe_log2(value / historical_month_median)) %>%
  filter(is.finite(log_fold_difference))

count_check <- plot_data %>%
  count(source, metric, period, name = "n") %>%
  filter(n < 2)
if (nrow(count_check)) {
  stop("At least one nutrient/source/period group has fewer than two observations.")
}

tests <- plot_data %>%
  group_by(source, metric) %>%
  summarise(
    n_historical = sum(period == "Historical"),
    n_2020 = sum(period == "2020"),
    input_file = first(input_file),
    value_units = first(value_units),
    analytical_definition = first(analytical_definition),
    wilcox_w_2020 = {
      historical <- log_fold_difference[period == "Historical"]
      focal <- log_fold_difference[period == "2020"]
      as.numeric(wilcox.test(focal, historical, exact = FALSE)$statistic)
    },
    p_value = {
      historical <- log_fold_difference[period == "Historical"]
      focal <- log_fold_difference[period == "2020"]
      wilcox.test(historical, focal, exact = FALSE)$p.value
    },
    median_2020_log_fold = median(log_fold_difference[period == "2020"]),
    .groups = "drop"
  ) %>%
  mutate(
    mann_whitney_u_2020 = wilcox_w_2020 - n_2020 * (n_2020 + 1) / 2,
    p_adjusted_bh = p.adjust(p_value, method = "BH"),
    effect_symbol = case_when(
      p_adjusted_bh < .05 & median_2020_log_fold > 0 ~ "(+)",
      p_adjusted_bh < .05 & median_2020_log_fold < 0 ~ "(-)",
      TRUE ~ ""
    ),
    effect_direction = case_when(
      effect_symbol == "(+)" ~ "higher",
      effect_symbol == "(-)" ~ "lower",
      TRUE ~ "not significant"
    ),
    effect_colour = case_when(
      effect_symbol == "(+)" ~ "firebrick3",
      effect_symbol == "(-)" ~ "#0072B2",
      TRUE ~ "transparent"
    ),
    focal_window = "August-October 2020",
    historical_reference = paste0(
      "All other available years; normalized to source-, metric-, and ",
      "calendar-month-specific medians"
    ),
    maximum_in_lake_depth_m = if_else(source == "In-lake", MAX_DEPTH_M, NA_real_),
    test = "two-sided Wilcoxon rank-sum (Mann-Whitney)",
    multiplicity_correction = "Benjamini-Hochberg across 12 source-metric comparisons"
  )

reference <- plot_data %>%
  filter(period == "Historical") %>%
  group_by(source, metric) %>%
  summarise(
    historical_median = median(log_fold_difference),
    historical_2.5_percentile = quantile(log_fold_difference, .025),
    historical_97.5_percentile = quantile(log_fold_difference, .975),
    reference_year_min = min(Year),
    reference_year_max = max(Year),
    reference_years = paste(sort(unique(Year)), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(tests, by = c("source", "metric"))

shared_x <- range(
  c(plot_data$log_fold_difference,
    reference$historical_2.5_percentile,
    reference$historical_97.5_percentile),
  na.rm = TRUE
)

make_panel <- function(source_name, focal_colour, title, metric_levels, tag) {
  panel_data <- plot_data %>% filter(source == source_name)
  historical <- panel_data %>% filter(period == "Historical")
  focal <- panel_data %>% filter(period == "2020")
  panel_reference <- reference %>% filter(source == source_name)
  label_x <- shared_x[2] + 0.25

  ggplot() +
    geom_vline(xintercept = 0, colour = "grey35",
               linewidth = .42 * LO_FIG_SCALE) +
    geom_violin(
      data = historical,
      aes(x = log_fold_difference, y = metric),
      orientation = "y", fill = COL_HIST, colour = "grey55", alpha = .70,
      linewidth = .32 * LO_FIG_SCALE, trim = FALSE
    ) +
    geom_segment(
      data = panel_reference,
      aes(x = historical_2.5_percentile, xend = historical_97.5_percentile,
          y = metric, yend = metric),
      colour = "grey20", linewidth = 1.0 * LO_FIG_SCALE
    ) +
    geom_point(
      data = panel_reference,
      aes(x = historical_median, y = metric),
      shape = 21, fill = "white", colour = "grey15",
      size = 2.2 * LO_FIG_SCALE, stroke = .35
    ) +
    geom_point(
      data = focal,
      aes(x = log_fold_difference, y = metric),
      position = position_jitter(height = .10, width = 0, seed = FOCAL_YEAR),
      shape = 21, fill = focal_colour, colour = focal_colour,
      size = 2.35 * LO_FIG_SCALE, stroke = .25, alpha = .95
    ) +
    geom_text(
      data = panel_reference %>% filter(effect_symbol != ""),
      aes(x = label_x, y = metric, label = effect_symbol,
          colour = effect_colour),
      hjust = 0, vjust = .35, size = 9 / LO_PT_PER_MM,
      fontface = "bold", family = LO_FONT, show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_y_discrete(
      limits = rev(metric_levels),
      labels = unname(METRIC_DISPLAY[rev(metric_levels)])
    ) +
    scale_x_continuous(
      breaks = pretty(shared_x, n = 6),
      expand = expansion(mult = c(.03, .16))
    ) +
    coord_cartesian(xlim = c(shared_x[1] - .18, label_x + .20), clip = "off") +
    labs(
      tag = tag, title = title,
      x = "Log2-fold difference (2020 vs. historical)", y = NULL
    ) +
    theme_classic(base_size = 9, base_family = LO_FONT) +
    theme(
      axis.title.x = element_text(size = 9),
      axis.text = element_text(size = 8.5, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 9, face = "bold"),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(.01, .99),
      plot.margin = margin(2, 8, 2, 2, "mm")
    )
}

p_dep <- make_panel(
  "Deposition", COL_DEP_2020, "Atmospheric nutrient deposition",
  c(NUT_LEVELS_DEP, "DIN:SRP"), "a"
)
p_lake <- make_panel(
  "In-lake", COL_LAKE_2020, "Surface-water nutrients (0-10 m)",
  c(NUT_LEVELS_LAKE, "DIN:TRP"), "b"
)
supplemental_2020_nutrients <- p_dep / p_lake + plot_layout(heights = c(1, 1))

coverage <- plot_data %>%
  count(source, metric, period, Year, Month, name = "n_observations") %>%
  arrange(source, metric, period, Year, Month)

figure_only <- identical(Sys.getenv("CALDOR_2020_FIGURE_ONLY"), "1")
if (!figure_only) {
  write_csv(plot_data %>% arrange(source, metric, Date),
            file.path(OUT_DIR, "supplemental_2020_nutrients_logfold_data.csv"))
  write_csv(reference,
            file.path(OUT_DIR, "supplemental_2020_nutrients_reference.csv"))
  write_csv(tests,
            file.path(OUT_DIR, "supplemental_2020_nutrients_mann_whitney_results.csv"))
  write_csv(coverage,
            file.path(OUT_DIR, "supplemental_2020_nutrients_coverage.csv"))
}

caption <- paste0(
  "August-October 2020 nutrient conditions relative to same-calendar-month ",
  "observations from all other available years. (a) Atmospheric deposition and ",
  "(b) surface-water nutrients at 0-10 m. Gray violins show historical ",
  "log2-fold differences after each observation was divided by its source-, ",
  "metric-, and calendar-month-specific historical median; black bars and open ",
  "circles show the historical 95% interval and median. Colored points are 2020 ",
  "observations. Red (+) and blue (-) denote significantly higher or lower 2020 ",
  "values, respectively, using two-sided Wilcoxon rank-sum tests with ",
  "Benjamini-Hochberg correction across 12 source-by-metric comparisons ",
  "(adjusted p < 0.05); unlabeled rows were not significant. Deposition panels ",
  "use SRP, TP, and DIN:SRP; surface-water panels use TRP, THP, and DIN:TRP."
)
if (!figure_only) {
  writeLines(caption, file.path(OUT_DIR, "supplemental_2020_nutrients_caption.txt"))
}

png_file <- file.path(OUT_DIR, "supplemental_2020_nutrients_historical_comparison.png")
pdf_file <- file.path(OUT_DIR, "supplemental_2020_nutrients_historical_comparison.pdf")
S2020_WIDTH_CM <- 7.62
S2020_HEIGHT_CM <- S2020_WIDTH_CM * 14.5 / LO_WIDTH_DOUBLE
ggsave(
  png_file, supplemental_2020_nutrients,
  width = S2020_WIDTH_CM, height = S2020_HEIGHT_CM, units = "cm", dpi = LO_DPI_LINE,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  pdf_file, supplemental_2020_nutrients,
  width = S2020_WIDTH_CM, height = S2020_HEIGHT_CM, units = "cm",
  device = grDevices::cairo_pdf, bg = "white"
)
message("Saved 2020 nutrient comparison: ", png_file, " and ", pdf_file)
