# =============================================================================
# Figure 2 final: daily/cast observations for calendar year 2021
#
# Two otherwise identical versions are produced:
#   1. Lake Number (Lake-Tools)
#   2. Schmidt stability (Lake-Tools)
#
# Observations are never reduced to monthly means when finer-resolution data
# exist. Historical temperature/monthly metric ribbons and the daily shortwave
# ribbon are 95% t confidence intervals among years.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")

KD_FILE <- file.path("data", "lake_environmental_data", "uv",
                     "LTP_Kd_Thermocline_depth_results.csv")
SECCHI_FILE <- file.path("data", "lake_environmental_data", "secchi",
                         "Secchi_LTP.csv")
SURFACE_TEMP_FILE <- file.path("data", "processed", "ctd",
                               "ltp_surface_temperature_0_10m_from_pkl.csv")
STABILITY_FILE <- file.path("data", "processed", "ctd",
                            "mltp_lake_tools_stability_metrics.csv")
USCG_FILE <- file.path("data", "lake_environmental_data", "00_met", "USCG.csv")
OUT_DIR <- file.path("figures", "figure_2_radiative")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LAKE_FILE <- file.path(OUT_DIR, "figure_2_2021_lake_number_final.png")
SCHMIDT_FILE <- file.path(OUT_DIR, "figure_2_2021_schmidt_stability_final.png")
CANONICAL_FILE <- file.path(OUT_DIR, "figure_2_final.png")
PLOT_DATA_FILE <- file.path(OUT_DIR, "figure_2_2021_final_plot_data.csv")
SUMMARY_FILE <- file.path(OUT_DIR, "figure_2_2021_final_summary.csv")

required <- c(KD_FILE, SECCHI_FILE, SURFACE_TEMP_FILE, STABILITY_FILE, USCG_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing Figure 2 inputs: ", paste(missing, collapse = ", "))

BASE_FAMILY <- "Times New Roman"
R2_MIN <- 0.90
FOCAL_START <- as.Date("2021-01-01")
FOCAL_END <- as.Date("2021-12-31")
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END <- as.Date("2021-10-21")
COL_PAR <- "#2E7D32"
COL_UV320 <- "#E65100"
COL_SECCHI <- "#1565C0"
COL_TEMP <- "#7B3294"
COL_SW <- "#D89C00"
COL_LAKE_NUMBER <- "#0072B2"
COL_SCHMIDT <- "#008B8B"
COL_CLIM <- "grey52"

t_ci <- function(mean_value, sd_value, n_value) {
  if_else(
    n_value >= 2L & is.finite(sd_value),
    qt(0.975, df = n_value - 1L) * sd_value / sqrt(n_value),
    NA_real_
  )
}

monthly_climatology <- function(data, value_col, start_year = NULL,
                                log_transform = FALSE) {
  source_data <- data %>%
    filter(year != 2021L, is.finite(.data[[value_col]]))
  if (!is.null(start_year)) source_data <- source_data %>% filter(year >= start_year)
  annual_month <- source_data %>%
    group_by(year, month) %>%
    summarise(year_month_value = mean(.data[[value_col]], na.rm = TRUE),
              .groups = "drop")
  if (log_transform) {
    annual_month %>%
      filter(year_month_value > 0) %>%
      mutate(transformed = log10(year_month_value)) %>%
      group_by(month) %>%
      summarise(n = n(), center = mean(transformed), spread = sd(transformed),
                .groups = "drop") %>%
      mutate(ci = t_ci(center, spread, n),
             clim_mean = 10^center,
             clim_lo = 10^(center - ci), clim_hi = 10^(center + ci)) %>%
      filter(n >= 2L, is.finite(clim_mean))
  } else {
    annual_month %>%
      group_by(month) %>%
      summarise(n = n(), clim_mean = mean(year_month_value),
                clim_sd = sd(year_month_value), .groups = "drop") %>%
      mutate(ci = t_ci(clim_mean, clim_sd, n),
             clim_lo = clim_mean - ci, clim_hi = clim_mean + ci) %>%
      filter(n >= 2L, is.finite(clim_mean))
  }
}

repeat_monthly_climatology <- function(clim) {
  tibble(month = 1:12,
         plot_date = as.Date(sprintf("2021-%02d-15", month))) %>%
    left_join(clim, by = "month")
}

# ---- Optical observations and historical monthly confidence intervals ------
kd_raw <- read_csv(KD_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date), year = year(date), month = month(date)) %>%
  filter(!is.na(date))
kd <- kd_raw %>%
  mutate(
    par_1pct = if_else(coalesce(R_PAR, 0) >= R2_MIN & coalesce(Kd_PAR, 0) > 0,
                       4.605 / Kd_PAR, NA_real_),
    uv320_1pct = if_else(coalesce(R_320, 0) >= R2_MIN & coalesce(Kd_320, 0) > 0,
                         4.605 / Kd_320, NA_real_)
  )
secchi <- read_csv(SECCHI_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date_Time_Local), year = year(date), month = month(date)) %>%
  filter(!is.na(date), is.finite(Secchi))

clim_par <- repeat_monthly_climatology(
  monthly_climatology(kd, "par_1pct", 2015L)
) %>% mutate(clim_lo = pmax(clim_lo, 0))
clim_uv320 <- repeat_monthly_climatology(
  monthly_climatology(kd, "uv320_1pct", 2015L)
) %>% mutate(clim_lo = pmax(clim_lo, 0))
clim_secchi <- repeat_monthly_climatology(
  monthly_climatology(secchi, "Secchi", 2015L)
) %>% mutate(clim_lo = pmax(clim_lo, 0))

dat_par <- kd_raw %>%
  filter(date >= FOCAL_START, date <= FOCAL_END,
         is.finite(Kd_PAR), Kd_PAR > 0) %>%
  transmute(date, value = 4.605 / Kd_PAR,
            qc_pass = coalesce(R_PAR, 0) >= R2_MIN)
dat_uv320 <- kd_raw %>%
  filter(date >= FOCAL_START, date <= FOCAL_END,
         is.finite(Kd_320), Kd_320 > 0) %>%
  transmute(date, value = 4.605 / Kd_320,
            qc_pass = coalesce(R_320, 0) >= R2_MIN)
dat_secchi <- secchi %>%
  filter(date >= FOCAL_START, date <= FOCAL_END) %>%
  transmute(date, value = Secchi, qc_pass = TRUE)

# ---- Every LTP CTD cast; monthly climatology CI uses years as replicates ----
ctd_surface <- read_csv(SURFACE_TEMP_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), year = year(date), month = month(date)) %>%
  filter(is.finite(surface_temperature_c))
clim_temp <- repeat_monthly_climatology(
  monthly_climatology(ctd_surface, "surface_temperature_c")
)
dat_temp <- ctd_surface %>%
  filter(date >= FOCAL_START, date <= FOCAL_END) %>%
  transmute(date, value = surface_temperature_c, qc_pass = TRUE)

# ---- USCG shortwave monthly means; historical and focal 95% CIs -------------
message("Reading USCG SWin and calculating monthly means and 95% CIs...")
sw <- fread(USCG_FILE, select = c("dateutc", "SWin"), showProgress = FALSE)
sw[, dateutc := as.POSIXct(dateutc, tz = "UTC")]
sw <- sw[is.finite(SWin) & SWin >= 0]
sw[, date := as.IDate(dateutc)]
sw_daily <- sw[, .(sw_w_m2 = mean(SWin), n_intervals = .N), by = date][
  n_intervals >= 72L
] %>%
  as_tibble() %>%
  mutate(date = as.Date(date), year = year(date), month = month(date))

sw_climatology <- repeat_monthly_climatology(
  monthly_climatology(sw_daily, "sw_w_m2", 2005L)
) %>% mutate(clim_lo = pmax(clim_lo, 0))
sw_focal <- sw_daily %>%
  filter(date >= FOCAL_START, date <= FOCAL_END) %>%
  group_by(year, month) %>%
  summarise(
    value = mean(sw_w_m2), daily_sd = sd(sw_w_m2), n_days = n(),
    .groups = "drop"
  ) %>%
  mutate(
    date = make_date(year, month, 15L),
    obs_ci = qt(0.975, pmax(n_days - 1L, 1L)) * daily_sd / sqrt(n_days),
    obs_lo = pmax(value - obs_ci, 0), obs_hi = value + obs_ci,
    qc_pass = n_days >= 2L
  )

# ---- Every MLTP cast; Lake-Tools metrics and historical monthly CIs ---------
stability <- read_csv(STABILITY_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), year = year(date), month = month(date))
stability_2021 <- stability %>% filter(year == 2021L)
if (!nrow(stability_2021)) stop("No 2021 Lake-Tools metrics were found.")
if (any(!is.finite(stability_2021$schmidt_stability_j_m2))) {
  stop("A 2021 MLTP cast has missing Schmidt stability; see the audit CSV.")
}

lake_obs <- stability_2021 %>%
  filter(is.finite(lake_number), lake_number > 0) %>%
  transmute(date, value = lake_number, qc_pass = TRUE)
schmidt_obs <- stability_2021 %>%
  transmute(date, value = schmidt_stability_j_m2, qc_pass = TRUE)
lake_clim <- repeat_monthly_climatology(
  monthly_climatology(stability, "lake_number", 2005L, log_transform = TRUE)
)
schmidt_clim <- repeat_monthly_climatology(
  monthly_climatology(stability, "schmidt_stability_j_m2", 2005L)
) %>% mutate(clim_lo = pmax(clim_lo, 0))

# ---- Shared plotting --------------------------------------------------------
base_theme <- theme_bw(base_size = 7.2, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.22),
    axis.text = element_text(size = 5.7, colour = "grey20"),
    axis.title = element_text(size = 6.5, face = "bold"),
    axis.title.y = element_text(size = 6.5, face = "bold"),
    axis.title.x = element_blank(),
    plot.tag = element_text(size = 9, face = "bold"),
    plot.tag.position = c(0.015, 0.985),
    plot.margin = margin(1.2, 1.8, 1.2, 1.8, "mm")
  )

fire_layers <- list(
  annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
           ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.07),
  geom_vline(xintercept = CALDOR_START, linetype = "dashed",
             colour = "firebrick", linewidth = 0.35)
)
focal_axis <- scale_x_date(
  limits = c(FOCAL_START, FOCAL_END),
  breaks = as.Date(c("2021-01-01", "2021-03-01", "2021-05-01",
                     "2021-07-01", "2021-09-01", "2021-11-01",
                     "2021-12-31")),
  labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov", "Dec"),
  expand = expansion(mult = 0)
)

make_panel <- function(clim, obs, y_lab, colour, tag,
                       reverse_y = FALSE, fire_label = FALSE,
                       log_y = FALSE, daily_climatology = FALSE,
                       y_limits = NULL) {
  if (!"segment" %in% names(obs)) obs <- obs %>% mutate(segment = 1L)
  p <- ggplot() +
    fire_layers +
    geom_ribbon(
      data = clim, aes(plot_date, ymin = clim_lo, ymax = clim_hi),
      fill = COL_CLIM, alpha = 0.24, na.rm = TRUE
    ) +
    geom_line(data = clim, aes(plot_date, clim_mean),
              colour = COL_CLIM,
              linewidth = if (daily_climatology) 0.35 else 0.65,
              na.rm = TRUE) +
    {if (!daily_climatology) geom_point(
      data = clim, aes(plot_date, clim_mean),
      colour = COL_CLIM, size = 1.15, na.rm = TRUE
    )} +
    {if (all(c("obs_lo", "obs_hi") %in% names(obs))) geom_errorbar(
      data = obs %>% filter(qc_pass),
      aes(date, ymin = obs_lo, ymax = obs_hi),
      colour = colour, width = 7, linewidth = 0.35, alpha = 0.8,
      na.rm = TRUE
    )} +
    geom_line(data = obs %>% filter(qc_pass), aes(date, value, group = segment),
              colour = colour,
              linewidth = if (daily_climatology) 0.35 else 0.55,
              alpha = 0.75, na.rm = TRUE) +
    geom_point(data = obs %>% filter(qc_pass), aes(date, value),
               colour = colour,
               size = if (daily_climatology) 0.55 else 1.45,
               alpha = if (daily_climatology) 0.65 else 0.9, na.rm = TRUE) +
    geom_point(data = obs %>% filter(!qc_pass), aes(date, value),
               colour = colour, shape = 1, stroke = 0.45, size = 1.45,
               alpha = 0.75, na.rm = TRUE) +
    focal_axis + labs(y = y_lab, tag = tag) + base_theme
  if (reverse_y) p <- p + scale_y_reverse()
  if (log_y) p <- p + scale_y_log10(labels = label_scientific(digits = 2))
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  if (fire_label) {
    p <- p + annotate("text", x = as.Date("2021-09-17"), y = Inf,
                      label = "Caldor Fire", colour = "firebrick",
                      vjust = 1.4, size = 2.0, family = BASE_FAMILY)
  }
  p
}

p_a <- make_panel(clim_par, dat_par, "PAR 1% depth (m)",
                  COL_PAR, "a", TRUE, TRUE)
p_b <- make_panel(clim_uv320, dat_uv320, "UV-320 1% depth (m)",
                  COL_UV320, "b", TRUE)
p_c <- make_panel(clim_secchi, dat_secchi, "Secchi depth (m)",
                  COL_SECCHI, "c", TRUE)
p_d <- make_panel(
  clim_temp, dat_temp, expression("Surface temperature ("*degree*"C)"),
  COL_TEMP, "d"
)
p_e <- make_panel(
  sw_climatology, sw_focal, expression("Shortwave radiation (W m"^-2*")"),
  COL_SW, "e"
)

compose_figure <- function(metric_panel) {
  (p_a | p_d) / (p_b | p_e) / (p_c | metric_panel) +
    plot_layout(heights = c(1, 1, 1))
}

p_lake <- make_panel(
  lake_clim, lake_obs, "Lake Number (dimensionless)",
  COL_LAKE_NUMBER, "f", log_y = TRUE, y_limits = c(0.5, 2000)
)
p_schmidt <- make_panel(
  schmidt_clim, schmidt_obs,
  expression("Schmidt stability (J m"^-2*")"), COL_SCHMIDT, "f"
)
fig_lake <- compose_figure(p_lake)
fig_schmidt <- compose_figure(p_schmidt)

for (file in c(LAKE_FILE, CANONICAL_FILE)) {
  ggsave(file, fig_lake, width = 5, height = 6, units = "in", dpi = LO_DPI,
         device = ragg::agg_png, bg = "white")
}
ggsave(SCHMIDT_FILE, fig_schmidt, width = 5, height = 6, units = "in",
       dpi = LO_DPI, device = ragg::agg_png, bg = "white")

export_series <- function(dat, panel, series, date_col = "date",
                          value_col = "value") {
  dat %>% transmute(panel = panel, series = series,
                    date = .data[[date_col]], value = .data[[value_col]],
                    lower = if ("clim_lo" %in% names(dat)) clim_lo else if ("obs_lo" %in% names(dat)) obs_lo else NA_real_,
                    upper = if ("clim_hi" %in% names(dat)) clim_hi else if ("obs_hi" %in% names(dat)) obs_hi else NA_real_,
                    n = if ("n" %in% names(dat)) n else NA_integer_)
}
plot_export <- bind_rows(
  export_series(dat_par, "PAR 1% depth", "2021 observations"),
  export_series(dat_uv320, "UV-320 1% depth", "2021 observations"),
  export_series(dat_secchi, "Secchi depth", "2021 observations"),
  export_series(dat_temp, "Surface temperature", "2021 cast observations"),
  export_series(sw_focal, "USCG SWin", "2021 monthly means with 95% CI"),
  export_series(lake_obs, "Lake Number", "2021 cast observations"),
  export_series(schmidt_obs, "Schmidt stability", "2021 cast observations"),
  export_series(clim_par, "PAR 1% depth", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(clim_uv320, "UV-320 1% depth", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(clim_secchi, "Secchi depth", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(clim_temp, "Surface temperature", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(sw_climatology, "USCG SWin", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(lake_clim, "Lake Number", "historical monthly climatology", "plot_date", "clim_mean"),
  export_series(schmidt_clim, "Schmidt stability", "historical monthly climatology", "plot_date", "clim_mean")
)
write_csv(plot_export, PLOT_DATA_FILE, na = "")
write_csv(
  tribble(
    ~metric, ~n, ~note,
    "LTP surface-temperature casts", nrow(dat_temp),
    "All 2021 casts; 95% CI is among historical annual monthly means",
    "USCG monthly SWin values", nrow(sw_focal),
    "Monthly means of valid daily means; days require >=72 ten-minute records; error bars are 95% CIs among days and historical 95% CI is among annual monthly means",
    "Lake Number MLTP casts", nrow(lake_obs),
    "Finite 2021 cast values; current Lake-Tools with nearest raw TB4 wind and no averaging; exact-zero wind is undefined upstream",
    "Schmidt-stability MLTP casts", nrow(schmidt_obs),
    "All 2021 casts; no missing metric values; Lake-Tools",
    "Inference note", NA_integer_,
    "Grey ribbons are 95% confidence intervals for the historical monthly climatological mean. The resistance/recovery table now uses these same Figure 2 confidence bands to classify departures and recovery."
  ),
  SUMMARY_FILE
)

message("Saved: ", LAKE_FILE)
message("Saved: ", SCHMIDT_FILE)
message("Updated canonical Figure 2: ", CANONICAL_FILE)
