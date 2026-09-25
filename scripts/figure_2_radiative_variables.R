# =============================================================================
# figure_2_radiative_variables.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Figure 2 - Light transparency at LTP: 2021 vs. multi-year climatology
#   3-panel figure (stacked):
#   Panel 1: PAR 1% depth       [4.605 / Kd_PAR]  (m)
#   Panel 2: UV-320 1% depth    [4.605 / Kd_320]  (m)
#   Panel 3: Secchi depth       (m)
#
#   Each panel: climatology = monthly means +/- 95% CI (all years except 2021,
#   NA removed); 2021 = all individual observations as points.
#   Quality filter for Kd: R-squared >= 0.90.
#   Caldor Fire window (Aug 14 - Oct 21, 2021) shaded.
#
# Data:
#   data/lake_environmental_data/uv/LTP_Kd_Thermocline_depth_results.csv
#   data/lake_environmental_data/secchi/Secchi_LTP.csv
# Output: figures/figure_2_radiative/figure_2_radiative.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})
source("scripts/figure_aesthetics.R")

# ---- Paths ------------------------------------------------------------------
KD_FILE <- file.path(
  "data", "lake_environmental_data", "uv",
  "LTP_Kd_Thermocline_depth_results.csv"
)
SECCHI_FILE <- file.path(
  "data", "lake_environmental_data", "secchi", "Secchi_LTP.csv"
)
OUT_DIR <- file.path("figures", "figure_2_radiative")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY  <- "Times New Roman"
R2_MIN       <- 0.90
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")

# ---- Colours ----------------------------------------------------------------
COL_PAR    <- "#2E7D32"   # PAR 1% depth    — dark green
COL_UV320  <- "#E65100"   # UV-320 1% depth — dark orange
COL_SECCHI <- "#1565C0"   # Secchi depth    — dark blue
COL_CLIM   <- "grey55"    # climatology line/ribbon

# =============================================================================
# SECTION 1: Load and derive variables
# =============================================================================

kd_raw <- read_csv(KD_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date), year = year(date), month = month(date)) %>%
  filter(!is.na(date))

kd <- kd_raw %>%
  mutate(
    par_1pct   = if_else(coalesce(R_PAR, 0) >= R2_MIN & coalesce(Kd_PAR, 0) > 0,
                         4.605 / Kd_PAR, NA_real_),
    uv320_1pct = if_else(coalesce(R_320, 0) >= R2_MIN & coalesce(Kd_320, 0) > 0,
                         4.605 / Kd_320, NA_real_)
  )

secchi_raw <- read_csv(SECCHI_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date_Time_Local),
         year  = year(date),
         month = month(date)) %>%
  filter(!is.na(date), !is.na(Secchi))

cat("Kd rows:", nrow(kd),
    "| valid PAR 1%:", sum(!is.na(kd$par_1pct)),
    "| valid UV320 1%:", sum(!is.na(kd$uv320_1pct)), "\n")
cat("Secchi rows:", nrow(secchi_raw), "\n")

# =============================================================================
# SECTION 2: Monthly climatology (all years except 2021)
# =============================================================================

make_clim <- function(df, col, excl_yr = 2021) {
  df %>%
    filter(year >= 2015, year <= 2025, year != excl_yr, !is.na(.data[[col]])) %>%
    group_by(month) %>%
    summarise(
      n         = n(),
      clim_mean = mean(.data[[col]], na.rm = TRUE),
      clim_sd   = sd(.data[[col]],   na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      clim_se  = clim_sd / sqrt(pmax(n, 1)),
      clim_ci  = qt(0.975, df = pmax(n - 1, 1)) * clim_se,
      clim_lo  = pmax(clim_mean - clim_ci, 0),
      clim_hi  = clim_mean + clim_ci,
      plot_date = as.Date(paste0("2021-", sprintf("%02d", month), "-15"))
    ) %>%
    filter(!is.na(clim_mean))
}

clim_par   <- make_clim(kd,         "par_1pct")
clim_uv320 <- make_clim(kd,         "uv320_1pct")
clim_sec   <- make_clim(secchi_raw, "Secchi")

make_anomaly_test <- function(df, value_col, metric_label) {
  dat <- df %>%
    filter(year >= 2015, year <= 2025, !is.na(.data[[value_col]])) %>%
    transmute(year, month, value = .data[[value_col]])
  hist <- dat %>%
    filter(year != 2021) %>%
    group_by(month) %>%
    summarise(hist_mean = mean(value), hist_sd = sd(value), .groups = "drop")
  z_2021 <- dat %>%
    filter(year == 2021) %>%
    left_join(hist, by = "month") %>%
    filter(is.finite(hist_sd), hist_sd > 0) %>%
    mutate(anomaly_sd = (value - hist_mean) / hist_sd)
  wt <- if (nrow(z_2021) >= 2) {
    suppressWarnings(wilcox.test(z_2021$anomaly_sd, mu = 0, exact = FALSE))
  } else {
    NULL
  }
  tibble(
    metric = metric_label,
    test = "one-sample Wilcoxon signed-rank test of month-standardized 2021 anomalies",
    statistic = if (is.null(wt)) NA_real_ else unname(wt$statistic),
    p_value = if (is.null(wt)) NA_real_ else wt$p.value,
    n_2021 = nrow(z_2021),
    median_anomaly_sd = median(z_2021$anomaly_sd, na.rm = TRUE),
    mean_anomaly_sd = mean(z_2021$anomaly_sd, na.rm = TRUE),
    sd_anomaly_sd = sd(z_2021$anomaly_sd, na.rm = TRUE)
  )
}

figure_2_tests <- bind_rows(
  make_anomaly_test(kd, "par_1pct", "PAR 1% depth"),
  make_anomaly_test(kd, "uv320_1pct", "UV-320 1% depth"),
  make_anomaly_test(secchi_raw, "Secchi", "Secchi depth")
)
write_csv(figure_2_tests, file.path(OUT_DIR, "figure_2_radiative_tests.csv"))

# 2021 individual observations — ALL rows (quality filter shown as shape)
dat21_par <- kd_raw %>%
  filter(year == 2021, !is.na(Kd_PAR), Kd_PAR > 0) %>%
  mutate(par_1pct = 4.605 / Kd_PAR,
         qc_pass  = coalesce(R_PAR, 0) >= R2_MIN)

dat21_uv320 <- kd_raw %>%
  filter(year == 2021, !is.na(Kd_320), Kd_320 > 0) %>%
  mutate(uv320_1pct = 4.605 / Kd_320,
         qc_pass   = coalesce(R_320, 0) >= R2_MIN)

dat21_sec <- secchi_raw %>%
  filter(year == 2021) %>%
  mutate(qc_pass = TRUE)

cat("Climatology months — PAR:", paste(clim_par$month,   collapse=","),
    "| UV320:", paste(clim_uv320$month, collapse=","),
    "| Secchi:", paste(clim_sec$month,  collapse=","), "\n")
cat("2021 obs — PAR:", nrow(dat21_par),
    "| UV320:", nrow(dat21_uv320),
    "| Secchi:", nrow(dat21_sec), "\n")

# ---- Diagnostic: why July 2021 is missing -----------------------------------
cat("\n--- 2021 Kd data: all rows with quality flags ---\n")
kd_2021_diag <- kd_raw %>%
  filter(year == 2021) %>%
  mutate(par_pass  = coalesce(R_PAR, 0) >= R2_MIN & coalesce(Kd_PAR, 0) > 0,
         uv320_pass = coalesce(R_320, 0) >= R2_MIN & coalesce(Kd_320, 0) > 0) %>%
  select(date, month, Kd_PAR, R_PAR, par_pass, Kd_320, R_320, uv320_pass) %>%
  arrange(date)
print(as.data.frame(kd_2021_diag), digits = 3)
cat("\n2021 Secchi months available:", sort(unique(secchi_raw$month[secchi_raw$year == 2021])), "\n\n")

# =============================================================================
# SECTION 3: Shared theme
# =============================================================================

base_theme <- theme_bw(base_size = 8 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
    axis.text        = ggplot2::element_text(size = 5),
    axis.title.y     = ggplot2::element_text(size = 6),
    axis.title.x     = element_blank(),
    plot.title       = element_text(size = LO_FS_TITLE, face = "bold")
  )

# =============================================================================
# SECTION 4: Panel-building function
# =============================================================================

make_panel <- function(clim, dat21, y_col, y_lab, col21,
                       add_fire_label = FALSE, y_lo = NULL) {

  dat21_pass <- dat21 %>% filter(qc_pass)
  dat21_fail <- dat21 %>% filter(!qc_pass)

  p <- ggplot() +
    # Caldor Fire window
    annotate("rect",
             xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf,
             fill = "firebrick", alpha = 0.08) +
    geom_vline(xintercept = CALDOR_START,
               linetype = "dashed", colour = "firebrick", linewidth = 0.45 * LO_FIG_SCALE) +
    # Climatology CI ribbon + line
    geom_ribbon(data = clim,
                aes(x = plot_date, ymin = clim_lo, ymax = clim_hi),
                fill = COL_CLIM, alpha = 0.28, na.rm = TRUE) +
    geom_line(data = clim,
              aes(x = plot_date, y = clim_mean),
              colour = COL_CLIM, linewidth = 0.90 * LO_FIG_SCALE, na.rm = TRUE) +
    geom_point(data = clim,
               aes(x = plot_date, y = clim_mean),
               colour = COL_CLIM, size = 1.6 * LO_FIG_SCALE, shape = 16, na.rm = TRUE) +
    # 2021: quality-passed — filled circles + connecting line
    geom_line(data = dat21_pass,
              aes(x = date, y = .data[[y_col]]),
              colour = col21, linewidth = 0.70 * LO_FIG_SCALE, alpha = 0.65, na.rm = TRUE) +
    geom_point(data = dat21_pass,
               aes(x = date, y = .data[[y_col]]),
               colour = col21, size = 2.2 * LO_FIG_SCALE, shape = 16, alpha = 0.90, na.rm = TRUE) +
    # 2021: quality-failed — hollow open circles (R² < 0.90)
    geom_point(data = dat21_fail,
               aes(x = date, y = .data[[y_col]]),
               colour = col21, size = 2.2 * LO_FIG_SCALE, shape = 1, stroke = 0.7 * LO_FIG_SCALE,
               alpha = 0.70, na.rm = TRUE) +
    scale_x_date(limits      = c(as.Date("2021-01-01"), as.Date("2021-12-31")),
                 date_breaks = "1 month", date_labels = "%b",
                 expand      = expansion(mult = 0.01)) +
    scale_y_reverse() +
    labs(y = y_lab) +
    base_theme

  if (add_fire_label) {
    p <- p + annotate("text",
                      x = CALDOR_START + 19, y = Inf,
                      label = "Caldor Fire", hjust = 0.5, vjust = 1.5,
                      colour = "firebrick", size = 1.3 * LO_FIG_SCALE,
                      family = BASE_FAMILY)
  }
  p
}

# =============================================================================
# SECTION 5: Build panels
# =============================================================================

p1 <- make_panel(clim_par,   dat21_par,   "par_1pct",   "PAR 1% depth (m)",    COL_PAR,
                 add_fire_label = TRUE)
p2 <- make_panel(clim_uv320, dat21_uv320, "uv320_1pct", "UV-320 1% depth (m)", COL_UV320)
p3 <- make_panel(clim_sec,   dat21_sec,   "Secchi",     "Secchi depth (m)",    COL_SECCHI)

p1 <- p1 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2 <- p2 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# =============================================================================
# SECTION 6: Combine and save
# =============================================================================

fig2 <- p1 / p2 / p3 &
  theme(plot.margin = margin(0.8, 1.2, 0.8, 1.2, "mm"))

ragg::agg_png(
  file.path(OUT_DIR, "figure_2_radiative.png"),
  width = 5.08, height = 6.35, units = "cm", res = LO_DPI
)
print(fig2)
invisible(dev.off())

message("Saved: figure_2_radiative.png")
message("Done. Output in: ", OUT_DIR)
