# Rebuild the 2021-2022 standalone stability figures from Lake-Tools metrics.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")

METRIC_FILE <- file.path("data", "processed", "ctd",
                         "mltp_lake_tools_stability_metrics.csv")
PLOT_DATA_FILE <- file.path("data", "processed", "ctd",
                            "mltp_stability_metrics_2021_2022_plot_data.csv")
OUT_DIR <- file.path("figures", "supplemental", "ctd_profiles")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FOCAL_START <- as.Date("2021-01-01")
FOCAL_END <- as.Date("2022-12-31")
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END <- as.Date("2021-10-21")

metrics <- read_csv(METRIC_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), year = year(date), month = month(date))
focal <- metrics %>%
  filter(date >= FOCAL_START, date <= FOCAL_END)
if (any(!is.finite(focal$schmidt_stability_j_m2)) ||
    any(!is.finite(focal$lake_number))) {
  stop("The 2021-2022 Lake-Tools stability series contains a missing value.")
}

monthly_climatology <- function(value_col, log_scale = FALSE) {
  annual_month <- metrics %>%
    filter(year >= 2005L, !year %in% c(2021L, 2022L),
           is.finite(.data[[value_col]])) %>%
    group_by(year, month) %>%
    summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  if (log_scale) {
    summary <- annual_month %>%
      filter(value > 0) %>% mutate(value = log10(value)) %>%
      group_by(month) %>%
      summarise(n = n(), center = mean(value), spread = sd(value), .groups = "drop") %>%
      mutate(ci = qt(0.975, n - 1) * spread / sqrt(n),
             clim_mean = 10^center, clim_lo = 10^(center - ci),
             clim_hi = 10^(center + ci))
  } else {
    summary <- annual_month %>%
      group_by(month) %>%
      summarise(n = n(), clim_mean = mean(value), spread = sd(value),
                .groups = "drop") %>%
      mutate(ci = qt(0.975, n - 1) * spread / sqrt(n),
             clim_lo = pmax(clim_mean - ci, 0), clim_hi = clim_mean + ci)
  }
  tibble(month_start = seq(as.Date("2021-01-01"), as.Date("2022-12-01"), by = "month")) %>%
    mutate(date = month_start + days(14), month = month(month_start)) %>%
    left_join(summary, by = "month")
}

schmidt_clim <- monthly_climatology("schmidt_stability_j_m2")
lake_clim <- monthly_climatology("lake_number", log_scale = TRUE)

make_metric_plot <- function(value_col, clim, y_label, colour, log_scale = FALSE) {
  p <- ggplot() +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = .Machine$double.xmin, ymax = Inf,
             fill = "firebrick", alpha = 0.07) +
    geom_ribbon(data = clim, aes(date, ymin = clim_lo, ymax = clim_hi),
                fill = "grey65", alpha = 0.25) +
    geom_line(data = clim, aes(date, clim_mean), colour = "grey45",
              linewidth = LO_LW_MID) +
    geom_point(data = clim, aes(date, clim_mean), colour = "grey45",
               size = 0.9 * LO_FIG_SCALE) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = LO_LW_MID) +
    geom_line(data = focal, aes(date, .data[[value_col]]), colour = colour,
              linewidth = LO_LW_MID, alpha = 0.75) +
    geom_point(data = focal, aes(date, .data[[value_col]]), colour = colour,
               size = 1.25 * LO_FIG_SCALE) +
    scale_x_date(limits = c(FOCAL_START, FOCAL_END), date_breaks = "3 months",
                 date_labels = "%b\n%Y", expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 7, base_family = LO_FONT) +
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"),
          plot.margin = margin(2, 2, 2, 2, "mm"))
  if (log_scale) p <- p + scale_y_log10(labels = label_scientific(digits = 2))
  p
}

p_schmidt <- make_metric_plot(
  "schmidt_stability_j_m2", schmidt_clim,
  expression("Schmidt stability (J m"^-2*")"), "#008B8B"
)
p_lake <- make_metric_plot(
  "lake_number", lake_clim, "Lake Number (dimensionless)", "#0072B2", TRUE
)

ggsave(file.path(OUT_DIR, "ctd_schmidt_stability_mltp_2021_2022.png"),
       p_schmidt, width = 12.7, height = 7.6, units = "cm", dpi = LO_DPI,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT_DIR, "ctd_schmidt_stability_mltp_2021_2022.pdf"),
       p_schmidt, width = 12.7, height = 7.6, units = "cm",
       device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, "ctd_lake_number_mltp_2021_2022.png"),
       p_lake, width = 12.7, height = 7.6, units = "cm", dpi = LO_DPI,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT_DIR, "ctd_stability_metrics_mltp_2021_2022.png"),
       p_schmidt / p_lake, width = 12.7, height = 14.5, units = "cm",
       dpi = LO_DPI, device = ragg::agg_png, bg = "white")

plot_data <- bind_rows(
  focal %>% transmute(metric = "Schmidt stability", series = "observations",
                      date, value = schmidt_stability_j_m2,
                      lower = NA_real_, upper = NA_real_, n = NA_integer_),
  focal %>% transmute(metric = "Lake Number", series = "observations",
                      date, value = lake_number,
                      lower = NA_real_, upper = NA_real_, n = NA_integer_),
  schmidt_clim %>% transmute(metric = "Schmidt stability", series = "climatology",
                            date, value = clim_mean, lower = clim_lo,
                            upper = clim_hi, n),
  lake_clim %>% transmute(metric = "Lake Number", series = "climatology",
                         date, value = clim_mean, lower = clim_lo,
                         upper = clim_hi, n)
)
write_csv(plot_data, PLOT_DATA_FILE, na = "")
message("Rebuilt standalone 2021-2022 stability figures with Lake-Tools.")
