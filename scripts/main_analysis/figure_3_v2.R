# =============================================================================
# figure_3_v2.R -- Caldor Fire nutrient pathways
# Keeps the Figure 3 Aug-Oct filters and Wilcoxon tests, but displays each
# response relative to the seasonal historical baseline.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(ragg)
})
source('scripts/mixed_layer_nutrient_budget.R', local = FALSE)
source("scripts/main_analysis/shared_aesthetics.R")

NUTRI_MLTP <- file.path("data", "lake_environmental_data", "nutrients", "Tahoe_MLTP_Nutrient.csv")
DEPO_FILE <- file.path("data", "lake_environmental_data", "deposition", "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv")
ASH_FILE <- file.path("data", "experiments", "Table 1- Individual jar regression stats DO rates of change.xlsx")
OUT_DIR <- file.path("figures", "figure_3_nutrients")
MIRROR_DIR <- file.path("figures_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MIRROR_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY <- "Times New Roman"
FIRE_YEAR <- 2021
FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
SEASON_MONTHS <- 8:10
MAX_DEPTH_M <- 10
N_ATOMIC_MASS <- 14.0067
P_ATOMIC_MASS <- 30.973762
COL_DEP_HIST <- "#BDBDBD"
COL_DEP_2021 <- "#D55E00"
COL_LAK_2021 <- "#0072B2"
NUT_LEVELS_DEP <- c("NO3", "NH4", "TKN", "SRP", "TP")
NUT_LEVELS_LAKE <- c("NO3", "NH4", "TKN", "TRP", "THP")
NUT_LEVELS <- c("NO3", "NH4", "TKN", "SRP", "TP", "TRP", "THP")
NUT_DISPLAY <- c(
  NO3 = "NO3", NH4 = "NH4", TKN = "TKN", SRP = "SRP", TP = "TP",
  TRP = "TRP", THP = "THP"
)

# ---- Established Figure 3 data processing -----------------------------------
lake_long <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Month %in% SEASON_MONTHS) %>%
  select(Date, Year, Month, NO3, NH4, TRP, THP, TKN) %>%
  pivot_longer(c(NO3, NH4, TRP, THP, TKN), names_to = "nut_raw", values_to = "value") %>%
  mutate(nutrient = recode(nut_raw, TRP = "TRP", THP = "THP"),
         nutrient = factor(nutrient, levels = NUT_LEVELS_LAKE), source = "In-lake",
         period = if_else(Year == FIRE_YEAR, "2021", "Historical")) %>%
  filter(!is.na(value))

dep_raw <- read_csv(DEPO_FILE, show_col_types = FALSE) %>%
  filter(!is.na(Start_Datetime), !is.na(End_Datetime)) %>%
  mutate(Start_Date = as.Date(Start_Datetime), End_Date = as.Date(End_Datetime),
         Mid_Date = Start_Date + as.integer(End_Date - Start_Date) %/% 2L,
         Year = year(Mid_Date), Month = month(Mid_Date))

dep_long <- dep_raw %>%
  filter(Month %in% SEASON_MONTHS) %>%
  select(Date = Mid_Date, Year, Month, NO3 = NO3_Daily_Load, NH4 = NH4_Daily_Load,
         TKN = TKN_Daily_Load, SRP = SRP_Daily_Load, TP = TP_Daily_Load) %>%
  pivot_longer(c(NO3, NH4, TKN, SRP, TP), names_to = "nutrient", values_to = "value") %>%
  mutate(nutrient = factor(nutrient, levels = NUT_LEVELS_DEP), source = "Deposition",
         period = if_else(Year == FIRE_YEAR, "2021", "Historical")) %>%
  filter(!is.na(value))

box_data <- bind_rows(lake_long, dep_long)

# The original two-sided Aug-Oct Wilcoxon comparisons are retained verbatim in
# scope, then BH-adjusted across the 10 source-by-nutrient comparisons.
wilcox_results <- box_data %>%
  group_by(source, nutrient) %>%
  summarise(
    n_hist = sum(period == "Historical"), n_2021 = sum(period == "2021"),
    p_val = {
      h <- value[period == "Historical"]
      f <- value[period == "2021"]
      if (length(h) >= 2 && length(f) >= 2) wilcox.test(h, f, exact = FALSE)$p.value else NA_real_
    }, .groups = "drop"
  ) %>%
  mutate(p_adj_bh = p.adjust(p_val, method = "BH"),
         sig_label = case_when(is.na(p_adj_bh) ~ "", p_adj_bh < .001 ~ "***",
                               p_adj_bh < .01 ~ "**", p_adj_bh < .05 ~ "*", TRUE ~ "ns"))

make_fold_summary <- function(data, source_name) {
  data %>%
    filter(source == source_name) %>%
    group_by(nutrient) %>%
    summarise(historical_median = median(value[period == "Historical"], na.rm = TRUE),
              historical_lo = quantile(value[period == "Historical"], .025, na.rm = TRUE),
              historical_hi = quantile(value[period == "Historical"], .975, na.rm = TRUE),
              estimate_2021 = median(value[period == "2021"], na.rm = TRUE), .groups = "drop") %>%
    left_join(wilcox_results %>% filter(source == source_name) %>%
                select(nutrient, p_val, p_adj_bh, sig_label), by = "nutrient") %>%
    mutate(fold = estimate_2021 / historical_median,
           interval_lo = historical_lo / historical_median,
           interval_hi = historical_hi / historical_median,
           nutrient = factor(nutrient, levels = if (source_name == "Deposition") NUT_LEVELS_DEP else NUT_LEVELS_LAKE))
}
dep_fold <- make_fold_summary(box_data, "Deposition")
lake_fold <- make_fold_summary(box_data, "In-lake")
safe_log2 <- function(x) log2(pmax(x, .Machine$double.eps))

make_enrichment_panel <- function(d, colour, title, ylab) {
  vals <- safe_log2(c(d$fold, d$interval_lo, d$interval_hi))
  label_x <- max(vals, na.rm = TRUE) + .20
  ggplot(d, aes(y = nutrient)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = .42 * LO_FIG_SCALE) +
    geom_segment(aes(x = safe_log2(interval_lo), xend = safe_log2(interval_hi), yend = nutrient),
                 colour = "grey45", linewidth = .72 * LO_FIG_SCALE) +
    geom_point(aes(x = safe_log2(fold)), shape = 21, fill = colour, colour = colour,
               size = 2.7 * LO_FIG_SCALE, stroke = .3 * LO_FIG_SCALE) +
    geom_text(aes(x = label_x, label = sig_label), hjust = 0, vjust = .35,
              size = 2.4 * LO_FIG_SCALE, fontface = "bold", family = BASE_FAMILY) +
    scale_y_discrete(limits = rev(NUT_LEVELS), labels = NUT_DISPLAY) +
    scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4),
                       labels = function(x) paste0(formatC(2^x, format = "fg", digits = 3), "Ã—"),
                       expand = expansion(mult = c(.03, .18))) +
    coord_cartesian(xlim = c(min(vals, na.rm = TRUE) - .18, label_x + .20), clip = "off") +
    labs(title = title, x = "2021 / historical Aug-Oct median (fold enrichment)", y = ylab) +
    theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
    theme(axis.title = element_text(size = 7.2 * LO_FIG_SCALE),
          axis.text = element_text(size = 6.4 * LO_FIG_SCALE),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 8 * LO_FIG_SCALE, face = "bold"),
          plot.margin = margin(1.5, 5.5, 1.5, 1.5, "mm"))
}

# ---- Panel B: seasonal atmospheric DIN:SRP anomaly --------------------------
dep_ratio_all <- dep_raw %>%
  mutate(DIN_v = suppressWarnings(as.numeric(NO3_Daily_Load)) + suppressWarnings(as.numeric(NH4_Daily_Load)),
         SRP_v = suppressWarnings(as.numeric(SRP_Daily_Load)), Date = Mid_Date,
         din_srp_ratio = (DIN_v / N_ATOMIC_MASS) / (SRP_v / P_ATOMIC_MASS)) %>%
  filter(!is.na(din_srp_ratio), is.finite(din_srp_ratio), din_srp_ratio > 0, Month %in% SEASON_MONTHS)
ratio_hist <- dep_ratio_all %>% filter(Year != FIRE_YEAR)
ratio_month_ref <- ratio_hist %>% group_by(Month) %>%
  summarise(monthly_median = median(din_srp_ratio, na.rm = TRUE), .groups = "drop")
ratio_anomaly <- dep_ratio_all %>% left_join(ratio_month_ref, by = "Month") %>%
  mutate(log2_anomaly = safe_log2(din_srp_ratio / monthly_median))
ratio_hist_anomaly <- ratio_anomaly %>% filter(Year != FIRE_YEAR)
ratio_caldor <- ratio_anomaly %>% filter(Year == FIRE_YEAR, Date >= FIRE_START, Date <= FIRE_END)
stoich_reference <- ratio_hist_anomaly %>%
  summarise(historical_median = median(log2_anomaly, na.rm = TRUE),
            historical_lo = quantile(log2_anomaly, .025, na.rm = TRUE),
            historical_hi = quantile(log2_anomaly, .975, na.rm = TRUE))

p_stoich <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = .42 * LO_FIG_SCALE) +
  geom_violin(data = ratio_hist_anomaly, aes(x = log2_anomaly, y = "Historical\nAug-Oct"),
              orientation = "y", fill = COL_DEP_HIST, colour = "grey55", alpha = .70,
              linewidth = .32 * LO_FIG_SCALE, trim = FALSE) +
  geom_segment(data = stoich_reference,
               aes(x = historical_lo, xend = historical_hi, y = "Historical\nAug-Oct", yend = "Historical\nAug-Oct"),
               colour = "grey20", linewidth = 1.05 * LO_FIG_SCALE) +
  geom_point(data = stoich_reference, aes(x = historical_median, y = "Historical\nAug-Oct"),
             shape = 21, fill = "white", colour = "grey15", size = 2.25 * LO_FIG_SCALE) +
  geom_point(data = ratio_caldor, aes(x = log2_anomaly, y = "Caldor\n2021"),
             position = position_jitter(height = .09, width = 0), shape = 21,
             fill = COL_DEP_2021, colour = COL_DEP_2021, size = 2.35 * LO_FIG_SCALE, alpha = .95) +
  scale_y_discrete(limits = c("Caldor\n2021", "Historical\nAug-Oct")) +
  scale_x_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), formatC(x, format = "f", digits = 1)),
                     breaks = pretty(c(ratio_hist_anomaly$log2_anomaly, ratio_caldor$log2_anomaly), n = 5)) +
  labs(title = "Atmospheric DIN:SRP was not unusually high",
       x = "logâ‚‚(DIN:SRP / historical monthly median)", y = NULL) +
  theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
  theme(axis.title = element_text(size = 7.2 * LO_FIG_SCALE),
        axis.text = element_text(size = 6.4 * LO_FIG_SCALE),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8 * LO_FIG_SCALE, face = "bold"),
        plot.subtitle = element_text(size = 6.2 * LO_FIG_SCALE, colour = "grey30"),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

p_dep <- make_enrichment_panel(dep_fold, COL_DEP_2021,
  "Atmospheric nutrient deposition increased during Caldor", "Nutrient")
p_lake <- make_enrichment_panel(lake_fold, COL_LAK_2021,
  "Surface-water nutrient response was comparatively weak", "Nutrient")

# ---- Figure 3 v3-style log-fold distributions ------------------------------
# Use the panel-B convention for every nutrient: the historical distribution is
# shown as a violin and all 2021 observations are overlaid on that same row.
# Values are log2-fold differences from the same-calendar-month historical
# median, so monthly seasonality is removed before comparing 2021 to history.
METRIC_LEVELS <- c(NUT_LEVELS, "DIN:SRP", "DIN:TRP")
METRIC_DISPLAY <- c(NUT_DISPLAY, "DIN:SRP" = "DIN:SRP", "DIN:TRP" = "DIN:TRP")

lake_ratio_fold <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Month %in% SEASON_MONTHS) %>%
  transmute(
    Date, Year, Month,
    value = (as.numeric(NO3) + as.numeric(NH4)) / (as.numeric(TRP) / P_ATOMIC_MASS),
    metric = "DIN:TRP", source = "In-lake",
    period = if_else(Year == FIRE_YEAR, "2021", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

dep_ratio_fold <- dep_raw %>%
  filter(Month %in% SEASON_MONTHS) %>%
  transmute(
    Date = Mid_Date, Year, Month,
    value = ((as.numeric(NO3_Daily_Load) + as.numeric(NH4_Daily_Load)) / N_ATOMIC_MASS) /
      (as.numeric(SRP_Daily_Load) / P_ATOMIC_MASS),
    metric = "DIN:SRP", source = "Deposition",
    period = if_else(Year == FIRE_YEAR, "2021", "Historical")
  ) %>%
  filter(is.finite(value), value > 0)

fold_plot_input <- bind_rows(
  box_data %>%
    transmute(Date, Year, Month, value, metric = as.character(nutrient), source, period),
  lake_ratio_fold,
  dep_ratio_fold
) %>%
  filter(is.finite(value), value > 0) %>%
  mutate(metric = factor(metric, levels = METRIC_LEVELS))

fold_month_reference <- fold_plot_input %>%
  filter(period == "Historical") %>%
  group_by(source, metric, Month) %>%
  summarise(historical_month_median = median(value, na.rm = TRUE), .groups = "drop")

fold_plot_data <- fold_plot_input %>%
  left_join(fold_month_reference, by = c("source", "metric", "Month")) %>%
  mutate(log_fold_difference = safe_log2(value / historical_month_median)) %>%
  filter(is.finite(log_fold_difference))

fold_tests <- fold_plot_data %>%
  group_by(source, metric) %>%
  summarise(
    n_hist = sum(period == "Historical"),
    n_2021 = sum(period == "2021"),
    wilcox_w_2021 = {
      historical <- log_fold_difference[period == "Historical"]
      focal <- log_fold_difference[period == "2021"]
      if (length(historical) >= 2L && length(focal) >= 2L) {
        as.numeric(wilcox.test(focal, historical, exact = FALSE)$statistic)
      } else NA_real_
    },
    p_val = {
      historical <- log_fold_difference[period == "Historical"]
      focal <- log_fold_difference[period == "2021"]
      if (length(historical) >= 2L && length(focal) >= 2L) {
        wilcox.test(historical, focal, exact = FALSE)$p.value
      } else NA_real_
    },
    median_2021_log_fold = median(log_fold_difference[period == "2021"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mann_whitney_u_2021 = wilcox_w_2021 - n_2021 * (n_2021 + 1) / 2,
    p_adj_bh = p.adjust(p_val, method = "BH"),
    sig_label = case_when(
      is.na(p_adj_bh) ~ "",
      p_adj_bh < .001 ~ "***",
      p_adj_bh < .01 ~ "**",
      p_adj_bh < .05 ~ "*",
      TRUE ~ ""
    ),
    effect_symbol = case_when(
      p_adj_bh < .05 & median_2021_log_fold > 0 ~ "(+)",
      p_adj_bh < .05 & median_2021_log_fold < 0 ~ "(-)",
      TRUE ~ ""
    ),
    effect_colour = case_when(
      effect_symbol == "(+)" ~ "firebrick3",
      effect_symbol == "(-)" ~ "#0072B2",
      TRUE ~ "transparent"
    )
  )

fold_reference <- fold_plot_data %>%
  filter(period == "Historical") %>%
  group_by(source, metric) %>%
  summarise(
    historical_median = median(log_fold_difference, na.rm = TRUE),
    historical_lo = quantile(log_fold_difference, .025, na.rm = TRUE),
    historical_hi = quantile(log_fold_difference, .975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(fold_tests, by = c("source", "metric"))

shared_fold_x <- range(c(
  fold_plot_data$log_fold_difference,
  fold_reference$historical_lo,
  fold_reference$historical_hi
), na.rm = TRUE)

make_logfold_distribution_panel <- function(source_name, focal_colour, title,
                                            metric_levels) {
  panel_data <- fold_plot_data %>% filter(source == source_name)
  hist_data <- panel_data %>% filter(period == "Historical")
  focal_data <- panel_data %>% filter(period == "2021")
  panel_ref <- fold_reference %>% filter(source == source_name)
  label_x <- shared_fold_x[2] + 0.25

  ggplot() +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = .42 * LO_FIG_SCALE) +
    geom_violin(
      data = hist_data,
      aes(x = log_fold_difference, y = metric),
      orientation = "y", fill = COL_DEP_HIST, colour = "grey55", alpha = .70,
      linewidth = .32 * LO_FIG_SCALE, trim = FALSE
    ) +
    geom_segment(
      data = panel_ref,
      aes(x = historical_lo, xend = historical_hi, y = metric, yend = metric),
      colour = "grey20", linewidth = 1.0 * LO_FIG_SCALE
    ) +
    geom_point(
      data = panel_ref,
      aes(x = historical_median, y = metric),
      shape = 21, fill = "white", colour = "grey15", size = 2.2 * LO_FIG_SCALE
    ) +
    geom_point(
      data = focal_data,
      aes(x = log_fold_difference, y = metric),
      position = position_jitter(height = .10, width = 0),
      shape = 21, fill = focal_colour, colour = focal_colour,
      size = 2.35 * LO_FIG_SCALE, alpha = .95
    ) +
    geom_text(
      data = panel_ref,
      aes(x = label_x, y = metric, label = effect_symbol, colour = effect_colour),
      hjust = 0, vjust = .35, size = 2.4 * LO_FIG_SCALE,
      fontface = "bold", family = BASE_FAMILY, show.legend = FALSE
    ) +
    scale_y_discrete(limits = rev(metric_levels), labels = unname(METRIC_DISPLAY[rev(metric_levels)])) +
    scale_x_continuous(
      breaks = pretty(shared_fold_x, n = 5),
      expand = expansion(mult = c(.03, .18))
    ) +
    coord_cartesian(xlim = c(shared_fold_x[1] - .18, label_x + .20), clip = "off") +
    scale_colour_identity() +
    labs(title = title, x = "Log2-fold difference (2021 vs. historical)", y = NULL) +
    theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
    theme(
      axis.title = element_text(size = 7.2 * LO_FIG_SCALE),
      axis.text = element_text(size = 6.4 * LO_FIG_SCALE),
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      plot.title = element_text(size = 8 * LO_FIG_SCALE, face = "bold"),
      plot.margin = margin(1.5, 5.5, 1.5, 1.5, "mm")
    )
}

p_dep <- make_logfold_distribution_panel(
  "Deposition", COL_DEP_2021,
  "Atmospheric nutrient deposition", c(NUT_LEVELS_DEP, "DIN:SRP")
)
p_lake <- make_logfold_distribution_panel(
  "In-lake", COL_LAK_2021,
  "Surface-water nutrients (0-10 m)", c(NUT_LEVELS_LAKE, "DIN:TRP")
)


COMPOSITION_FILE <- file.path(
  "data", "Figure 3", "Figure 3", "Close_Far_Composition_Analysis",
  "group_mean_nutrient_percentages.csv"
)
ASH_CHLA_FILE <- file.path("data", "Figure 3", "Figure 3", "ashexpt_chla.csv")
ASH_PPR_FILE <- file.path("data", "Figure 3", "Figure 3", "ashexpt_PPR.csv")

composition_colours <- c(
  Ca = "#4E79A7", K = "#F28E2B", Fe = "#E15759", Mg = "#76B7B2",
  TP = "#59A14F", Na = "#EDC948", Mn = "#B07AA1", Zn = "#FF9DA7",
  Cu = "#9C755F", Co = "#BAB0AC", Pb = "#86BCB6", Cr = "#FFBE7D",
  Ni = "#8CD17D", Cd = "#B699C0", As = "#F1CE63"
)
composition_order <- list(
  "Major nutrients" = c("Ca", "K", "Fe", "Mg", "TP", "Na", "Mn"),
  "Minor nutrients" = c("Zn", "Cu", "Co"),
  "Trace metals" = c("Pb", "Cr", "Ni", "Cd", "As")
)
composition_backgrounds <- c(
  "Major nutrients" = "#EDF3F8", "Minor nutrients" = "#FDF1F2",
  "Trace metals" = "#EDF7EE"
)

composition_data <- read_csv(COMPOSITION_FILE, show_col_types = FALSE) %>%
  pivot_longer(c(Close_mean_percent, Far_mean_percent),
               names_to = "distance_group", values_to = "percent") %>%
  mutate(
    distance_group = recode(
      distance_group,
      Close_mean_percent = "Close\nAsh",
      Far_mean_percent = "Far\nAsh"
    ),
    Category = factor(Category, levels = names(composition_order)),
    Element = factor(Element, levels = unlist(composition_order, use.names = FALSE))
  )

make_composition_subplot <- function(category, y_breaks, y_max, y_label = NULL) {
  category_data <- composition_data %>%
    filter(Category == category) %>%
    mutate(Element = factor(Element, levels = composition_order[[category]]))

  ggplot(category_data, aes(x = distance_group, y = percent, fill = Element)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = composition_backgrounds[[category]], colour = NA) +
    geom_col(width = .64, position = position_stack(reverse = TRUE),
             colour = "white", linewidth = .22 * LO_FIG_SCALE) +
    scale_fill_manual(values = composition_colours[composition_order[[category]]], name = NULL) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(title = category, x = NULL, y = y_label) +
    theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
    theme(
      axis.title = element_text(size = 7.2 * LO_FIG_SCALE),
      axis.text = element_text(size = 6.4 * LO_FIG_SCALE),
      axis.line = element_line(linewidth = LO_LW_THIN),
      axis.ticks = element_line(linewidth = LO_LW_THIN),
      plot.title = element_text(size = 7.4 * LO_FIG_SCALE, face = "bold", hjust = .5),
      legend.position = "bottom",
      legend.text = element_text(size = 5.1 * LO_FIG_SCALE),
      legend.key.height = unit(1.8, "mm"),
      legend.key.width = unit(2.2, "mm"),
      legend.spacing.x = unit(.25, "mm"),
      legend.margin = margin(t = -3, b = -2, unit = "pt"),
      plot.margin = margin(1.5, 0, 0, 0, "mm")
    )
}

p_composition_major <- make_composition_subplot(
  "Major nutrients", seq(0, 100, 25), 100,
  "Mean contribution\nto total deposition (%)"
)
p_composition_minor <- make_composition_subplot(
  "Minor nutrients", seq(0, .56, .14), .56
)
p_composition_trace <- make_composition_subplot(
  "Trace metals", seq(0, .12, .03), .12
)
p_composition <- wrap_plots(
  p_composition_major + labs(tag = "c"), p_composition_minor, p_composition_trace,
  ncol = 3, widths = c(1.25, .85, .95), guides = "collect"
) +
  plot_annotation(theme = theme(
    legend.position = "bottom",
    legend.box.just = "center",
    legend.direction = "horizontal"
  ))
make_ppr_normalized_panel <- function() {
  treatment_labels <- c(
    "a.control" = "Control",
    "b.NP" = "Inorganic\nnutrient",
    "h.BA1" = "Far Ash",
    "i.BA2" = "Close Ash"
  )
  chla <- read_csv(ASH_CHLA_FILE, show_col_types = FALSE) %>%
    filter(lake == "Tahoe", smoke_trt == "no smoke", treatment == "a.control")
  control_chla_mean <- mean(chla$ChlorophyllA, na.rm = TRUE)
  if (!is.finite(control_chla_mean) || control_chla_mean <= 0) {
    stop("Cannot normalize PPR: Tahoe no-smoke control chlorophyll-a mean is invalid.")
  }
  ppr_data <- read_csv(ASH_PPR_FILE, show_col_types = FALSE) %>%
    filter(lake == "Tahoe", smoke_trt == "no smoke", treatment %in% names(treatment_labels)) %>%
    transmute(
      sample, treatment = factor(treatment, levels = names(treatment_labels)),
      net_ppr = as.numeric(NET_PPR),
      normalized_ppr = net_ppr / control_chla_mean
    ) %>%
    filter(is.finite(normalized_ppr))
  ppr_summary <- ppr_data %>%
    group_by(treatment) %>%
    summarise(
      n = n(), mean_normalized_ppr = mean(normalized_ppr),
      sd_normalized_ppr = sd(normalized_ppr), .groups = "drop"
    )
  fit <- aov(normalized_ppr ~ treatment, data = ppr_data)
  tukey <- TukeyHSD(fit)$treatment %>% as.data.frame() %>% rownames_to_column("comparison")
  letters <- multcompView::multcompLetters4(fit, TukeyHSD(fit))$treatment
  ppr_letters <- tibble(
    treatment = factor(names(letters$Letters), levels = levels(ppr_data$treatment)),
    significance_group = unname(letters$Letters)
  )
  ppr_summary <- ppr_summary %>%
    left_join(ppr_letters, by = "treatment") %>%
    mutate(label_y = mean_normalized_ppr + sd_normalized_ppr + 0.30)

  write_csv(ppr_data, file.path(OUT_DIR, "figure_3_v2_ppr_normalized_data.csv"))
  write_csv(ppr_summary, file.path(OUT_DIR, "figure_3_v2_ppr_normalized_summary.csv"))
  write_csv(tukey, file.path(OUT_DIR, "figure_3_v2_ppr_normalized_tukey.csv"))

  ggplot(ppr_data, aes(x = treatment, y = normalized_ppr)) +
    geom_point(position = position_jitter(width = 0.07, height = 0, seed = 42),
               shape = 21, fill = "#56B4E9", colour = "grey20", size = 2.4,
               stroke = 0.35) +
    geom_errorbar(data = ppr_summary,
                  aes(x = treatment, ymin = mean_normalized_ppr - sd_normalized_ppr,
                      ymax = mean_normalized_ppr + sd_normalized_ppr),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.42) +
    geom_point(data = ppr_summary, aes(x = treatment, y = mean_normalized_ppr),
               inherit.aes = FALSE, shape = 23, fill = "black", colour = "black", size = 3) +
    geom_text(data = ppr_summary, aes(x = treatment, y = label_y, label = significance_group),
              inherit.aes = FALSE, fontface = "bold", size = 3.6) +
    scale_x_discrete(labels = treatment_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(
      x = NULL,
      y = expression(atop("Net PPR / Chl "*italic(a),
                           "(mg C mg Chl "*italic(a)^-1*" h"^-1*")"))
    ) +
    theme_classic(base_size = LO_FONT_MIN_PT, base_family = BASE_FAMILY) +
    theme(
      axis.title = element_text(size = LO_FONT_MIN_PT),
      axis.text.x = element_text(size = 6.6 * LO_FIG_SCALE, lineheight = .85,
                                 angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = LO_FONT_MIN_PT),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}
p_ppr <- make_ppr_normalized_panel() + labs(tag = "d")

# ---- Poster-only ash composition and experiment extract ---------------------
# Major and minor composition bars share a common x axis. Minor-element values
# are scaled solely for plotting; the secondary axis reports their native
# percent contributions.
POSTER_OUT <- file.path("figures", "for_poster")
dir.create(POSTER_OUT, showWarnings = FALSE, recursive = TRUE)

poster_composition_data <- composition_data %>%
  filter(Category %in% c("Major nutrients", "Minor nutrients")) %>%
  mutate(
    Element = factor(Element,
                     levels = c(composition_order[["Major nutrients"]],
                                composition_order[["Minor nutrients"]])),
    distance_group = factor(
      gsub("\\n", " ", distance_group),
      levels = c("Close Ash", "Far Ash")
    )
  )
poster_major <- poster_composition_data %>% filter(Category == "Major nutrients")
poster_minor <- poster_composition_data %>% filter(Category == "Minor nutrients")
poster_major_axis_max <- ceiling(max(poster_major$percent, na.rm = TRUE) / 10) * 10
poster_minor_axis_max <- ceiling(max(poster_minor$percent, na.rm = TRUE) / 0.1) * 0.1
poster_minor_scale <- poster_major_axis_max / poster_minor_axis_max
poster_distance_colours <- c("Close Ash" = "#D55E00", "Far Ash" = "#0072B2")

poster_composition <- ggplot() +
  geom_col(
    data = poster_major,
    aes(x = Element, y = percent, fill = distance_group),
    position = position_dodge(width = 0.76), width = 0.68,
    colour = "white", linewidth = 0.35
  ) +
  geom_col(
    data = poster_minor,
    aes(x = Element, y = percent * poster_minor_scale, fill = distance_group),
    position = position_dodge(width = 0.76), width = 0.68,
    colour = "white", linewidth = 0.35
  ) +
  geom_vline(xintercept = length(composition_order[["Major nutrients"]]) + 0.5,
             colour = "grey50", linewidth = 0.45, linetype = "dashed") +
  annotate("text", x = 4, y = poster_major_axis_max * 0.97,
           label = "Major nutrients", fontface = "bold", size = 5.0,
           family = BASE_FAMILY) +
  annotate("text", x = 9, y = poster_major_axis_max * 0.97,
           label = "Minor nutrients", fontface = "bold", size = 5.0,
           family = BASE_FAMILY) +
  scale_fill_manual(values = poster_distance_colours, name = NULL) +
  scale_y_continuous(
    name = "Major-nutrient contribution (%)",
    limits = c(0, poster_major_axis_max),
    breaks = seq(0, poster_major_axis_max, by = 10),
    sec.axis = sec_axis(~ . / poster_minor_scale,
                        name = "Minor-nutrient contribution (%)",
                        breaks = seq(0, poster_minor_axis_max, by = 0.1))
  ) +
  labs(title = "Ash composition", x = NULL) +
  theme_classic(base_size = 12, base_family = BASE_FAMILY) +
  theme(
    axis.title.y.left = element_text(size = 16, face = "bold", margin = margin(r = 6)),
    axis.title.y.right = element_text(size = 16, face = "bold", margin = margin(l = 6)),
    axis.text = element_text(size = 13, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 0.55),
    axis.ticks = element_line(linewidth = 0.55),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 13, face = "bold"),
    legend.key.size = unit(5, "mm"),
    plot.margin = margin(5, 10, 4, 8, "mm")
  )

poster_ppr_data <- p_ppr$data
poster_ppr_summary <- p_ppr$layers[[2]]$data
poster_treatment_colours <- c(
  "a.control" = "black",
  "b.NP" = "#009E73",
  "h.BA1" = poster_distance_colours[["Far Ash"]],
  "i.BA2" = poster_distance_colours[["Close Ash"]]
)
poster_experiment <- ggplot(poster_ppr_data, aes(x = treatment, y = normalized_ppr)) +
  geom_point(aes(fill = treatment),
             position = position_jitter(width = 0.07, height = 0, seed = 42),
             shape = 21, colour = "grey20", size = 3.1, stroke = 0.45,
             show.legend = FALSE) +
  geom_errorbar(
    data = poster_ppr_summary,
    aes(x = treatment, ymin = mean_normalized_ppr - sd_normalized_ppr,
        ymax = mean_normalized_ppr + sd_normalized_ppr),
    inherit.aes = FALSE, width = 0.12, linewidth = 0.55
  ) +
  geom_point(
    data = poster_ppr_summary,
    aes(x = treatment, y = mean_normalized_ppr, fill = treatment),
    inherit.aes = FALSE, shape = 23, colour = "black", size = 4.0,
    show.legend = FALSE
  ) +
  geom_text(
    data = poster_ppr_summary,
    aes(x = treatment, y = label_y, label = significance_group),
    inherit.aes = FALSE, fontface = "bold", size = 5.0, family = BASE_FAMILY
  ) +
  scale_fill_manual(values = poster_treatment_colours) +
  scale_x_discrete(labels = c(
    "a.control" = "Control",
    "b.NP" = "Inorganic\nnutrient",
    "h.BA1" = "Far\nash",
    "i.BA2" = "Close\nash"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(
    title = "Phytoplankton growth",
    x = NULL,
    y = expression(atop("Net PPR / Chl "*italic(a),
                         "(mg C mg Chl "*italic(a)^-1*" h"^-1*")"))
  ) +
  theme_classic(base_size = 12, base_family = BASE_FAMILY) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 5)),
    axis.text.x = element_text(size = 10.5, face = "bold", lineheight = 0.9,
                               angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.line = element_line(linewidth = 0.55),
    axis.ticks = element_line(linewidth = 0.55),
    plot.margin = margin(5, 6, 5, 6, "mm")
  )

poster_figure_3 <- wrap_plots(
  poster_composition, plot_spacer(), poster_experiment,
  ncol = 3, widths = c(2.25, 0.16, 1)
) & theme(plot.tag = element_blank())
ggsave(file.path(POSTER_OUT, "figure_3_ash_composition_experiment.png"),
       poster_figure_3, width = 30.0, height = 12.5, units = "cm", dpi = 300,
       device = ragg::agg_png, bg = "white")

# The combined view uses one y-axis row per nutrient. Atmospheric deposition
# and surface-water distributions are vertically offset within that row and
# identified with a source legend.
grouped_fold_data <- fold_plot_data %>%
  mutate(
    metric_index = as.numeric(factor(metric, levels = METRIC_LEVELS)),
    source_offset = if_else(source == "Deposition", 0.16, -0.16),
    y_position = metric_index + source_offset
  )
grouped_reference <- fold_reference %>%
  mutate(
    metric_index = as.numeric(factor(metric, levels = METRIC_LEVELS)),
    source_offset = if_else(source == "Deposition", 0.16, -0.16),
    y_position = metric_index + source_offset
  )
source_fills <- c("Deposition" = COL_DEP_HIST, "In-lake" = "#7CB9E8")
source_colours <- c("Deposition" = COL_DEP_2021, "In-lake" = COL_LAK_2021)
source_labels <- c("Deposition" = "Atmospheric deposition", "In-lake" = "Surface water (0-10 m)")

make_grouped_logfold_panel <- function() {
  historical <- grouped_fold_data %>% filter(period == "Historical")
  focal <- grouped_fold_data %>% filter(period == "2021")
  x_values <- c(grouped_fold_data$log_fold_difference,
                grouped_reference$historical_lo, grouped_reference$historical_hi)
  label_x <- max(x_values, na.rm = TRUE) + .25

  ggplot() +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = .42 * LO_FIG_SCALE) +
    geom_violin(
      data = historical,
      aes(x = log_fold_difference, y = y_position, fill = source, colour = source,
          group = interaction(metric, source)),
      orientation = "y", alpha = .35, linewidth = .32 * LO_FIG_SCALE,
      trim = FALSE, show.legend = TRUE
    ) +
    geom_segment(
      data = grouped_reference,
      aes(x = historical_lo, xend = historical_hi, y = y_position, yend = y_position,
          colour = source),
      linewidth = 0.85 * LO_FIG_SCALE, show.legend = FALSE
    ) +
    geom_point(
      data = grouped_reference,
      aes(x = historical_median, y = y_position, colour = source),
      shape = 21, fill = "white", size = 2.0 * LO_FIG_SCALE, show.legend = FALSE
    ) +
    geom_point(
      data = focal,
      aes(x = log_fold_difference, y = y_position, fill = source, colour = source),
      position = position_jitter(height = .07, width = 0), shape = 21,
      size = 2.3 * LO_FIG_SCALE, alpha = .95, show.legend = TRUE
    ) +
    geom_text(
      data = grouped_reference %>% filter(effect_symbol == "(+)"),
      aes(x = label_x, y = y_position, label = effect_symbol),
      hjust = 0, vjust = .35, size = 5.2 * LO_FIG_SCALE,
      fontface = "bold", family = BASE_FAMILY, colour = "firebrick3", show.legend = FALSE
    ) +
    geom_text(
      data = grouped_reference %>% filter(effect_symbol == "(-)"),
      aes(x = label_x, y = y_position, label = effect_symbol),
      hjust = 0, vjust = .35, size = 5.2 * LO_FIG_SCALE,
      fontface = "bold", family = BASE_FAMILY, colour = "#0072B2", show.legend = FALSE
    ) +
    scale_fill_manual(values = source_fills, labels = source_labels, name = NULL) +
    scale_colour_manual(values = source_colours, labels = source_labels, name = NULL) +
    scale_y_reverse(
      breaks = seq_along(METRIC_LEVELS), labels = unname(METRIC_DISPLAY[METRIC_LEVELS]),
      limits = c(length(METRIC_LEVELS) + .5, .5)
    ) +
    scale_x_continuous(
      breaks = pretty(range(x_values, na.rm = TRUE), n = 5),
      expand = expansion(mult = c(.03, .18))
    ) +
    coord_cartesian(xlim = c(min(x_values, na.rm = TRUE) - .18, label_x + .20), clip = "off") +
    labs(x = "Log2-fold difference (2021 vs. historical)", y = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = .65)),
           colour = guide_legend(override.aes = list(shape = 21, size = 2.4))) +
    theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
    theme(
      axis.title = element_text(size = 7.2 * LO_FIG_SCALE),
      axis.text.y = element_text(size = 11.5 * LO_FIG_SCALE),
      axis.text.x = element_text(size = 6.4 * LO_FIG_SCALE),
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      legend.position = "bottom", legend.text = element_text(size = 6.7 * LO_FIG_SCALE),
      legend.key.height = unit(2.5, "mm"), legend.key.width = unit(4, "mm"),
      plot.margin = margin(1.5, 8, 1.5, 1.5, "mm")
    )
}

p_grouped <- make_grouped_logfold_panel()

# ---- Existing ash-addition experiment, retained without analytical changes ---
make_ash_addition_panel <- function() {
  if (!file.exists(ASH_FILE)) return(NULL)
  ash_summary <- readxl::read_excel(ASH_FILE, sheet = "Individual jar regression stats") %>%
    transmute(light_treatment = .data[["Light Treatment"]],
              nutrient_treatment = .data[["Nutrient treatment"]], lake = Lake,
              period = Period, estimate = as.numeric(Estimate)) %>%
    filter(!is.na(estimate), lake == "Tahoe") %>%
    mutate(light_treatment = factor(light_treatment, levels = c("Dark", "Smoke", "No Smoke")),
           nutrient_treatment = factor(nutrient_treatment,
             levels = c("Control", "BA-1", "BA-2", "GA1", "GA2", "GA3", "N", "P", "NP"))) %>%
    group_by(period, light_treatment, nutrient_treatment) %>%
    summarise(mean_estimate = mean(estimate, na.rm = TRUE),
              se_estimate = sd(estimate, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  ggplot(ash_summary, aes(y = nutrient_treatment, x = mean_estimate, fill = light_treatment)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = .32 * LO_FIG_SCALE) +
    geom_col(position = position_dodge(width = .72), width = .64, colour = "grey25", linewidth = .15 * LO_FIG_SCALE) +
    geom_errorbar(aes(xmin = mean_estimate - se_estimate, xmax = mean_estimate + se_estimate),
                  position = position_dodge(width = .72), width = .16, linewidth = .25 * LO_FIG_SCALE) +
    facet_wrap(~period, nrow = 1) +
    scale_fill_manual(values = c("Dark" = "#3B3B3B", "Smoke" = "#B07AA1", "No Smoke" = "#4E79A7"), name = "Light") +
    labs(title = "Ash-addition experiment", x = "DO slope estimate", y = NULL) +
    theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
    theme(axis.title = element_text(size = 7.2 * LO_FIG_SCALE), axis.text = element_text(size = 5.7 * LO_FIG_SCALE),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 6.2 * LO_FIG_SCALE, face = "bold"), legend.position = "bottom",
          legend.title = element_text(size = 6 * LO_FIG_SCALE), legend.text = element_text(size = 5.6 * LO_FIG_SCALE),
          legend.key.size = unit(2.5, "mm"), plot.title = element_text(size = 8 * LO_FIG_SCALE, face = "bold"),
          plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
}
# This analysis is not displayed in the current four-panel layout; avoid
# reading and constructing it during Figure 3 rendering.
p_ash_addition <- NULL

make_placeholder <- function(tag) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey97", colour = "grey62", linewidth = 0.45,
             linetype = "dashed") +
    annotate("text", x = 0.5, y = 0.5, label = "Plot in\npreparation",
             colour = "grey42", size = 2.5, lineheight = 0.95,
             family = BASE_FAMILY) +
    annotate("text", x = 0.035, y = 0.975, label = tag,
             hjust = 0, vjust = 1, colour = "black", size = 3.2,
             fontface = "bold", family = BASE_FAMILY) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = BASE_FAMILY) +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
}

top_row <- wrap_plots(p_dep + labs(tag = "a"), p_lake + labs(tag = "b"), ncol = 2)
bottom_row <- wrap_plots(p_composition, p_ppr, ncol = 2, widths = c(2, 1.2))
figure_3_v2 <- wrap_plots(top_row, plot_spacer(), bottom_row, ncol = 1,
                           heights = c(1, .055, .92))
figure_3_v2 <- figure_3_v2 &
  theme(
    plot.tag = element_text(face = "bold", size = 9 * LO_FIG_SCALE,
                            family = BASE_FAMILY),
    plot.tag.position = c(.015, .985),
    axis.text.x = element_text(lineheight = .85),
    plot.subtitle = element_blank()
  )

write_csv(bind_rows(mutate(dep_fold, panel = "Atmospheric deposition"),
                    mutate(lake_fold, panel = "In-lake concentration")),
          file.path(OUT_DIR, "figure_3_v2_enrichment_summary.csv"))
write_csv(composition_data, file.path(OUT_DIR, "figure_3_v2_composition_data.csv"))
write_csv(fold_plot_data, file.path(OUT_DIR, "figure_3_v2_logfold_data.csv"))
write_csv(fold_reference, file.path(OUT_DIR, "figure_3_v2_logfold_reference.csv"))
write_csv(fold_tests, file.path(OUT_DIR, "figure_3_v2_mann_whitney_results.csv"))
caption_text <- paste0(
  "Nutrient pathways. (a) Atmospheric deposition and (b) surface-water (0-10 m) nutrient responses during the August-October 2021 Caldor Fire period. Labels use source-column names: deposition panels show SRP and TP, whereas surface-water panels show TRP and THP; molar ratios are DIN:SRP and DIN:TRP, respectively. Historical distributions (violins) and 2021 observations (points) are expressed as log2-fold differences from the same-calendar-month historical median. Black bars and open circles denote the historical 95% interval and median, respectively. Red (+) and blue (-) indicate a significantly higher or lower 2021 response, respectively, based on two-sided Mann-Whitney U tests (Wilcoxon rank-sum tests) with Benjamini-Hochberg correction across the 12 nutrient-source comparisons (pBH < 0.05); unlabeled rows were not significant. Panel (c) compares the mean elemental composition of deposition collected close to and farther from the fire. Major nutrients, minor nutrients, and trace metals use independent y-axis scales. Panel (d) shows Tahoe no-smoke NET PPR normalized by the mean control chlorophyll a concentration; points are individual observations, diamonds are means, error bars are SD, and letters are Tukey HSD groups."
)
writeLines(caption_text, file.path(OUT_DIR, "figure_3_v2_caption.txt"))
ggsave(file.path(OUT_DIR, "figure_3_v2.png"), figure_3_v2,
       width = 17.8, height = 15.2, units = "cm", dpi = LO_DPI,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT_DIR, "figure_3_v2.pdf"), figure_3_v2,
       width = 17.8, height = 15.2, units = "cm",
       device = grDevices::cairo_pdf, bg = "white")

# Keep the manuscript-facing figures_v2 mirror byte-identical to the verified
# v2 export so an older Figure 3 file cannot be mistaken for the current build.
mirror_ok <- file.copy(
  file.path(OUT_DIR, c("figure_3_v2.png", "figure_3_v2.pdf")),
  file.path(MIRROR_DIR, c("fig03_nutrient_pathways.png", "fig03_nutrient_pathways.pdf")),
  overwrite = TRUE
)
if (!all(mirror_ok)) stop("Could not update the canonical Figure 3 v2 mirror.")

message(
  "Saved Figure 3 v2 with Times New Roman: ",
  file.path(OUT_DIR, "figure_3_v2.png"),
  " and mirrored it to ", file.path(MIRROR_DIR, "fig03_nutrient_pathways.png")
)
