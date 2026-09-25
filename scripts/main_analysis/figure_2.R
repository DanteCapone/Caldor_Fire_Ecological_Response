# =============================================================================
# Final Figure 2: optical and physical context during the Caldor Fire
#
# Panels a-c show LTP optical-depth observations against monthly climatology.
# Panels d-e show LTP 0-10 m CTD temperature and NASA POWER shortwave radiation.
# Output: figures/figure_2_radiative/figure_2_final.png (5 x 6 inches)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
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
CTD_FILE <- file.path("data", "lake_environmental_data", "ctd", "ctd_casts",
                      "ctd_data_2005_2025.csv")
SW_FILE <- file.path("data", "raw", "swrad",
                     "tahoe_nasa_power_sw_2006_2025.csv")
HMS_FILE <- file.path("data", "processed", "tahoe_hms_smoke_daily.csv")
OUT_DIR <- Sys.getenv(
  "CALDOR_FIG2_OUT_DIR",
  unset = file.path("figures", "figure_2_radiative")
)
OUT_FILE <- file.path(OUT_DIR, "figure_2_final.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY <- "Times New Roman"
R2_MIN <- 0.90
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END <- as.Date("2021-10-21")
TAHOE_LAT <- 39.0868
COL_PAR <- "#2E7D32"
COL_UV320 <- "#E65100"
COL_SECCHI <- "#1565C0"
COL_TEMP <- "#7B3294"
COL_SW <- "#D89C00"
COL_DEFICIT <- "#B2182B"
COL_CLIM <- "grey52"

monthly_climatology <- function(data, value_col, start_year = 2015L) {
  data |>
    filter(year >= start_year, year <= 2025L, year != 2021L,
           is.finite(.data[[value_col]])) |>
    group_by(month) |>
    summarise(
      n = n(),
      clim_mean = mean(.data[[value_col]], na.rm = TRUE),
      clim_sd = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      clim_se = clim_sd / sqrt(pmax(n, 1)),
      clim_ci = qt(0.975, df = pmax(n - 1, 1)) * clim_se,
      clim_lo = clim_mean - clim_ci,
      clim_hi = clim_mean + clim_ci,
      plot_date = as.Date(sprintf("2021-%02d-15", month))
    ) |>
    filter(is.finite(clim_mean))
}

# ---- Optical observations ---------------------------------------------------
kd_raw <- read_csv(KD_FILE, show_col_types = FALSE) |>
  mutate(date = as.Date(Date), year = year(date), month = month(date)) |>
  filter(!is.na(date))

kd <- kd_raw |>
  mutate(
    par_1pct = if_else(coalesce(R_PAR, 0) >= R2_MIN & coalesce(Kd_PAR, 0) > 0,
                       4.605 / Kd_PAR, NA_real_),
    uv320_1pct = if_else(coalesce(R_320, 0) >= R2_MIN & coalesce(Kd_320, 0) > 0,
                         4.605 / Kd_320, NA_real_)
  )

secchi <- read_csv(SECCHI_FILE, show_col_types = FALSE) |>
  mutate(date = as.Date(Date_Time_Local), year = year(date), month = month(date)) |>
  filter(!is.na(date), is.finite(Secchi))

clim_par <- monthly_climatology(kd, "par_1pct") |>
  mutate(clim_lo = pmax(clim_lo, 0))
clim_uv320 <- monthly_climatology(kd, "uv320_1pct") |>
  mutate(clim_lo = pmax(clim_lo, 0))
clim_secchi <- monthly_climatology(secchi, "Secchi") |>
  mutate(clim_lo = pmax(clim_lo, 0))

dat21_par <- kd_raw |>
  filter(year == 2021L, is.finite(Kd_PAR), Kd_PAR > 0) |>
  mutate(value = 4.605 / Kd_PAR, qc_pass = coalesce(R_PAR, 0) >= R2_MIN)
dat21_uv320 <- kd_raw |>
  filter(year == 2021L, is.finite(Kd_320), Kd_320 > 0) |>
  mutate(value = 4.605 / Kd_320, qc_pass = coalesce(R_320, 0) >= R2_MIN)
dat21_secchi <- secchi |>
  filter(year == 2021L) |>
  transmute(date, value = Secchi, qc_pass = TRUE)

# ---- LTP 0-10 m CTD temperature --------------------------------------------
ctd_surface <- read_csv(
  CTD_FILE,
  col_select = c(CTD_ID, Station_ID, Cast_Date_Time_Local, Cast_Flag,
                 Depth, Temperature),
  show_col_types = FALSE,
  progress = FALSE
) |>
  mutate(date = as.Date(Cast_Date_Time_Local)) |>
  filter(
    Station_ID == "Index", Cast_Flag == 1,
    is.finite(Depth), between(Depth, 0, 10), is.finite(Temperature),
    !is.na(date)
  ) |>
  group_by(CTD_ID, date) |>
  summarise(surface_temp = mean(Temperature, na.rm = TRUE), .groups = "drop") |>
  mutate(year = year(date), month = month(date))

clim_temp <- monthly_climatology(ctd_surface, "surface_temp", 2005L)
dat21_temp <- ctd_surface |>
  filter(year == 2021L) |>
  transmute(date, value = surface_temp, qc_pass = TRUE)

# ---- NASA POWER shortwave and clear-sky deficit -----------------------------
sw_lines <- readLines(SW_FILE, warn = FALSE)
sw_start <- which(grepl("^YEAR", sw_lines))[1]
power <- read.csv(
  text = paste(sw_lines[sw_start:length(sw_lines)], collapse = "\n"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("-999", "-999.0")
) |>
  as_tibble() |>
  rename(year = YEAR, month = MO, day = DY,
         allsky_kwh = ALLSKY_SFC_SW_DWN) |>
  mutate(
    date = as.Date(sprintf("%04d-%02d-%02d", year, month, day)),
    doy = yday(date),
    sw_meas = allsky_kwh * 1000 / 24
  ) |>
  filter(year %in% 2006:2025, is.finite(sw_meas))

calc_ra <- function(doy, lat_deg) {
  lat <- lat_deg * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
  delta <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmin(1, pmax(-1, -tan(lat) * tan(delta))))
  (1367 / pi) * dr * (ws * sin(lat) * sin(delta) +
                        cos(lat) * cos(delta) * sin(ws))
}

power <- power |>
  mutate(ra = calc_ra(doy, TAHOE_LAT))

hms <- read_csv(HMS_FILE, show_col_types = FALSE) |>
  mutate(date = as.Date(date))
power_hms <- hms |>
  left_join(power |> select(date, sw_meas, ra), by = "date")

tau_clear <- quantile(
  power_hms$sw_meas[!power_hms$smoke_day & month(power_hms$date) %in% 7:9] /
    power_hms$ra[!power_hms$smoke_day & month(power_hms$date) %in% 7:9],
  probs = 0.90, na.rm = TRUE
)

sw_daily <- power_hms |>
  mutate(
    year = year(date), month = month(date), doy = yday(date),
    sw_clear = ra * tau_clear,
    sw_deficit = sw_clear - sw_meas,
    smoke_criterion = !is.na(sw_deficit) & sw_deficit > 20
  )

sw_monthly_climatology <- sw_daily |>
  filter(year != 2021L, is.finite(sw_meas)) |>
  group_by(month) |>
  summarise(
    n = n(), clim_mean = mean(sw_meas), clim_sd = sd(sw_meas), .groups = "drop"
  ) |>
  mutate(
    clim_se = clim_sd / sqrt(n),
    clim_ci = qt(0.975, df = pmax(n - 1, 1)) * clim_se,
    clim_lo = pmax(clim_mean - clim_ci, 0),
    clim_hi = clim_mean + clim_ci,
    plot_date = as.Date(sprintf("2021-%02d-15", month))
  )

sw_2021_monthly <- sw_daily |>
  filter(year == 2021L) |>
  group_by(month) |>
  summarise(
    date = as.Date(sprintf("2021-%02d-15", first(month))),
    sw_meas = mean(sw_meas, na.rm = TRUE),
    sw_deficit = mean(sw_deficit, na.rm = TRUE),
    smoke_days = sum(smoke_criterion, na.rm = TRUE),
    .groups = "drop"
  )

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
           ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.07),
  geom_vline(xintercept = CALDOR_START, linetype = "dashed",
             colour = "firebrick", linewidth = 0.35)
)

year_axis <- scale_x_date(
  limits = c(as.Date("2021-01-01"), as.Date("2021-12-31")),
  date_breaks = "2 months", date_labels = "%b",
  expand = expansion(mult = 0.01)
)

make_monthly_panel <- function(clim, obs, y_lab, colour, tag,
                               reverse_y = FALSE, fire_label = FALSE) {
  p <- ggplot() +
    fire_layers +
    geom_ribbon(data = clim,
                aes(plot_date, ymin = clim_lo, ymax = clim_hi),
                fill = COL_CLIM, alpha = 0.24, na.rm = TRUE) +
    geom_line(data = clim, aes(plot_date, clim_mean),
              colour = COL_CLIM, linewidth = 0.65, na.rm = TRUE) +
    geom_point(data = clim, aes(plot_date, clim_mean),
               colour = COL_CLIM, size = 1.15, na.rm = TRUE) +
    geom_line(data = obs |> filter(qc_pass), aes(date, value),
              colour = colour, linewidth = 0.55, alpha = 0.75, na.rm = TRUE) +
    geom_point(data = obs |> filter(qc_pass), aes(date, value),
               colour = colour, size = 1.45, alpha = 0.9, na.rm = TRUE) +
    geom_point(data = obs |> filter(!qc_pass), aes(date, value),
               colour = colour, shape = 1, stroke = 0.45, size = 1.45,
               alpha = 0.75, na.rm = TRUE) +
    year_axis + labs(y = y_lab, tag = tag) + base_theme
  if (reverse_y) p <- p + scale_y_reverse()
  if (fire_label) {
    p <- p + annotate("text", x = as.Date("2021-09-17"), y = Inf,
                      label = "Caldor Fire", colour = "firebrick",
                      vjust = 1.4, size = 2.0, family = BASE_FAMILY)
  }
  p
}

p_a <- make_monthly_panel(clim_par, dat21_par, "PAR 1% depth (m)",
                          COL_PAR, "a", TRUE, TRUE)
p_b <- make_monthly_panel(clim_uv320, dat21_uv320, "UV-320 1% depth (m)",
                          COL_UV320, "b", TRUE)
p_c <- make_monthly_panel(clim_secchi, dat21_secchi, "Secchi depth (m)",
                          COL_SECCHI, "c", TRUE)
p_d <- make_monthly_panel(clim_temp, dat21_temp,
                          expression("Mean 0-10 m temperature ("*degree*"C)"),
                          COL_TEMP, "d")

p_e <- ggplot() +
  fire_layers +
  geom_ribbon(data = sw_monthly_climatology,
              aes(plot_date, ymin = clim_lo, ymax = clim_hi),
              fill = COL_CLIM, alpha = 0.24, na.rm = TRUE) +
  geom_line(data = sw_monthly_climatology, aes(plot_date, clim_mean),
            colour = COL_CLIM, linewidth = 0.55, na.rm = TRUE) +
  geom_point(data = sw_monthly_climatology, aes(plot_date, clim_mean),
             colour = COL_CLIM, size = 1.15, na.rm = TRUE) +
  geom_line(data = sw_2021_monthly, aes(date, sw_meas),
            colour = COL_SW, linewidth = 0.60, alpha = 0.85, na.rm = TRUE) +
  geom_point(data = sw_2021_monthly, aes(date, sw_meas),
             colour = COL_SW, size = 1.45, na.rm = TRUE) +
  year_axis +
  labs(y = expression("Shortwave radiation (W m"^-2*")"), tag = "e") +
  base_theme

fig2 <- (p_a | p_d) / (p_b | p_e) / (p_c | plot_spacer()) +
  plot_layout(heights = c(1, 1, 1))

ggsave(
  OUT_FILE, fig2, width = 5, height = 6, units = "in", dpi = LO_DPI,
  device = ragg::agg_png, bg = "white"
)

write_csv(
  tibble(
    metric = c("LTP 0-10 m CTD temperature casts", "NASA POWER monthly mean shortwave",
               "2021 days exceeding the 20 W m-2 deficit criterion"),
    n = c(nrow(dat21_temp), nrow(sw_2021_monthly), sum(sw_daily$smoke_criterion[sw_daily$year == 2021L], na.rm = TRUE)),
    note = c(
      "Mean of accepted Index-station CTD observations from 0 through 10 m",
      "ALLSKY_SFC_SW_DWN converted from kWh m-2 day-1 to W m-2",
      sprintf("Clear-sky transmissivity tau_90 = %.4f", tau_clear)
    )
  ),
  file.path(OUT_DIR, "figure_2_final_summary.csv")
)

message("Saved final Figure 2: ", OUT_FILE)
