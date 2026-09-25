# =============================================================================
# Lake Number diagnostic, 2021
# Shows raw TB4 wind and raw-wind-resolution Lake Number alongside the
# cast-resolved hydrographic terms used to derive it.
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

TEMP_FILE <- file.path("data", "processed", "ctd", "mltp_temperature_profiles_1m_from_pkl.csv")
METRIC_FILE <- file.path("data", "processed", "ctd", "mltp_lake_tools_stability_metrics.csv")
RAW_LAKE_FILE <- file.path("data", "processed", "ctd", "mltp_lake_number_raw_wind_2021.csv")
OUT_DIR <- file.path("figures", "supplemental", "lake_number_diagnostic")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

YEAR_START <- as.POSIXct("2021-01-01 00:00:00", tz = "UTC")
YEAR_END <- as.POSIXct("2022-01-01 00:00:00", tz = "UTC")
DATE_START <- as.Date("2021-01-01")
DATE_END <- as.Date("2021-12-31")
FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
BASE_FAMILY <- "Times New Roman"

wind_lake <- read_csv(RAW_LAKE_FILE, show_col_types = FALSE) %>%
  transmute(
    datetime = as.POSIXct(dateutc, tz = "UTC"),
    wind_speed_m_s = as.numeric(wind_speed_m_s),
    wind_speed_forcing_m_s = as.numeric(wind_speed_forcing_m_s),
    lake_number_diagnostic = as.numeric(lake_number_diagnostic),
    interpolated_schmidt_stability_j_m2,
    interpolated_thermocline_depth_m,
    interpolated_epilimnion_density_kg_m3
  ) %>%
  arrange(datetime) %>%
  mutate(wind_24h_moving_mean_m_s = data.table::frollmean(
    wind_speed_m_s, n = 144L, align = "center", na.rm = TRUE
  ), lake_number_7d_moving_median = data.table::frollapply(
    lake_number_diagnostic, n = 1008L, FUN = median,
    align = "center", na.rm = TRUE
  ))

surface_temp <- read_csv(TEMP_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(year(date) == 2021L, depth_m >= 0, depth_m <= 10,
         is.finite(temperature_c)) %>%
  group_by(date) %>%
  summarise(surface_temperature_c = mean(temperature_c), .groups = "drop")

metrics <- read_csv(METRIC_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(year(date) == 2021L) %>%
  left_join(surface_temp, by = "date") %>%
  arrange(date)

if (nrow(metrics) != 12L || any(!is.finite(metrics$schmidt_stability_j_m2))) {
  stop("Expected 12 complete MLTP Schmidt-stability casts in 2021.")
}

fire_layer <- list(
  annotate("rect", xmin = FIRE_START, xmax = FIRE_END,
           ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.07),
  geom_vline(xintercept = FIRE_START, colour = "firebrick",
             linetype = "dashed", linewidth = 0.35)
)
lake_fire_layer <- list(
  annotate("rect", xmin = FIRE_START, xmax = FIRE_END,
           ymin = min(wind_lake$lake_number_diagnostic, na.rm = TRUE) * 0.8,
           ymax = max(wind_lake$lake_number_diagnostic, na.rm = TRUE) * 1.2,
           fill = "firebrick", alpha = 0.07),
  geom_vline(xintercept = FIRE_START, colour = "firebrick",
             linetype = "dashed", linewidth = 0.35)
)
x_scale <- scale_x_date(
  limits = c(DATE_START, DATE_END), date_breaks = "2 months",
  date_labels = "%b", expand = expansion(mult = 0)
)
panel_theme <- theme_bw(base_size = 8, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.22),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6),
    plot.tag = element_text(size = 9, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    plot.margin = margin(1, 2, 1, 2, "mm")
  )

p_wind <- ggplot(wind_lake, aes(as.Date(datetime))) +
  fire_layer +
  geom_line(aes(y = wind_speed_m_s, colour = "Raw TB4"),
            linewidth = 0.18, alpha = 0.28, na.rm = TRUE) +
  geom_line(aes(y = wind_24h_moving_mean_m_s, colour = "24-hour moving mean"),
            linewidth = 0.65, na.rm = TRUE) +
  geom_point(data = metrics,
             aes(x = date, y = wind_speed_forcing_m_s,
                 colour = "Nearest raw cast forcing"),
             size = 1.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c(
    "Raw TB4" = "grey65", "24-hour moving mean" = "#0072B2",
    "Nearest raw cast forcing" = "#D55E00"
  )) +
  x_scale + labs(x = NULL, y = expression("Wind speed (m s"^-1*")"), tag = "a") +
  panel_theme

p_temp <- ggplot(metrics, aes(date, surface_temperature_c)) +
  fire_layer + geom_line(colour = "#7B3294", linewidth = 0.6) +
  geom_point(colour = "#7B3294", size = 1.8) + x_scale +
  labs(x = NULL, y = expression("MLTP 0-10 m temperature ("*degree*"C)"), tag = "b") +
  panel_theme + theme(legend.position = "none")

p_schmidt <- ggplot(metrics, aes(date, schmidt_stability_j_m2)) +
  fire_layer + geom_line(colour = "#008B8B", linewidth = 0.6) +
  geom_point(colour = "#008B8B", size = 1.8) + x_scale +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = NULL, y = expression("Schmidt stability (J m"^-2*")"), tag = "c") +
  panel_theme + theme(legend.position = "none")

p_lake <- ggplot(wind_lake, aes(as.Date(datetime))) +
  lake_fire_layer +
  geom_line(aes(y = lake_number_diagnostic, colour = "Raw-wind resolution"),
            linewidth = 0.18, alpha = 0.28, na.rm = TRUE) +
  geom_line(aes(y = lake_number_7d_moving_median, colour = "7-day moving median"),
            linewidth = 0.65, na.rm = TRUE) +
  geom_point(data = metrics %>% filter(is.finite(lake_number)),
             aes(x = date, y = lake_number, colour = "Cast value"),
             inherit.aes = FALSE, size = 1.7) +
  scale_colour_manual(values = c(
    "Raw-wind resolution" = "#7FB3D5", "7-day moving median" = "#1B4F72",
    "Cast value" = "#D55E00"
  )) + x_scale +
  scale_y_log10(labels = label_number()) +
  labs(x = NULL, y = "Lake Number (log scale)", tag = "d") +
  panel_theme

diagnostic_plot <- p_wind / p_temp / p_schmidt / p_lake +
  plot_layout(heights = c(1.15, 1, 1, 1))

OUT_FILE <- file.path(OUT_DIR, "supplemental_lake_number_diagnostic_2021.png")
ggsave(OUT_FILE, diagnostic_plot, width = 12.7, height = 18.5, units = "cm",
       dpi = LO_DPI, device = ragg::agg_png, bg = "white")

write_csv(metrics, file.path(OUT_DIR, "lake_number_cast_diagnostics_2021.csv"))
write_csv(wind_lake, file.path(OUT_DIR, "tb4_raw_wind_and_lake_number_2021.csv"))
message("Saved: ", OUT_FILE)
