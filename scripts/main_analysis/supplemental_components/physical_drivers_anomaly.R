# =============================================================================
# supplemental_physical_drivers_anomaly.R
# Caldor Fire Ecosystem Response Project
#
# Purpose: Compile a suite of physical/chemical drivers at LTP/mid-lake and
#          plot shaded monthly anomalies (opposing colors for positive vs.
#          negative anomalies) and absolute values around the two
#          Leptolyngbya spp. bloom windows identified in the LTP time series,
#          plus one anomaly figure spanning the full 2005-2025 record:
#            bloom 1: 2011-04-16 to 2012-05-11 (identified from consecutive
#                     LTP Leptolyngbya detections in
#                     data/lake_environmental_data/phytoplankton/
#                     caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv)
#            bloom 2: 2021-08 onward (Caldor Fire, Aug 14 - Oct 21 2021)
#
# Variables (16 physical/chemical drivers, split into two columns, plus a
# Leptolyngbya abundance panel for the two bloom-window figures):
#   Column 1 (physical/whole-driver):
#     1. Surface temperature   - LTP CTD, 0-10 m mean Temperature per cast
#     2. PM2.5                 - Tahoe-wide compiled daily PM2.5
#     3. Upper Truckee River discharge (USGS)
#     4. Deposition DIN:SRP molar ratio (mid-lake atmospheric deposition)
#     5. In-lake DIN:TRP molar ratio (MLTP, 0-10 m)
#     6. 1% PAR depth          - LTP Kd_PAR regression results
#     7. Integrated Chl-a      - LTP Chl-a, trapezoidal integration, 0-100 m
#     8. Stratification        - MLTP Schmidt stability index
#   Column 2 (individual nutrients, deposition + in-lake pairs):
#     9-10.  NO3  (deposition mg m-2 d-1 / in-lake ug L-1)
#    11-12.  NH4  (deposition mg m-2 d-1 / in-lake ug L-1)
#    13-14.  TKN  (deposition mg m-2 d-1 / in-lake ug L-1)
#    15-16.  SRP/TRP  (deposition mg m-2 d-1 / in-lake ug L-1)
#
# Output:
#   data/processed/physical_drivers_compiled.csv   (long format, native
#                                                    sampling resolution)
#   data/processed/physical_drivers_monthly.csv     (monthly means used for
#                                                    the anomaly/value plots)
#   figures/supplemental/physical_drivers/supplemental_S6_physical_drivers_anomaly_2011_2012.png
#   figures/supplemental/physical_drivers/supplemental_S6b_physical_drivers_values_2011_2012.png
#   figures/supplemental/physical_drivers/supplemental_S7_physical_drivers_anomaly_2021_2022.png
#   figures/supplemental/physical_drivers/supplemental_S7b_physical_drivers_values_2021_2022.png
#   figures/supplemental/physical_drivers/supplemental_S8_physical_drivers_anomaly_2005_2025.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
  library(data.table)
})
source("scripts/main_analysis/shared_aesthetics.R")
source("scripts/main_analysis/supplemental_components/ctd_mixing_helpers.R")

proj_root <- normalizePath(".", winslash = "/")
DATA_DIR <- file.path(proj_root, "data")
PROC_DIR <- file.path(DATA_DIR, "processed")
OUT_DIR  <- file.path(proj_root, "figures", "supplemental", "physical_drivers")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BASE_FAMILY  <- "Times New Roman"
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")
LEPTO_BLOOM1_START <- as.Date("2011-04-16")
LEPTO_BLOOM1_END   <- as.Date("2012-03-21")  # last 0-10 m Leptolyngbya spp. abundance detection in 2011-2012
LEPTO_BLOOM2_START <- as.Date("2021-04-01")
LEPTO_BLOOM2_END   <- as.Date("2022-07-31")
LEPTO_COL    <- "#1B9E77"
N_ATOMIC_MASS <- 14.0067
P_ATOMIC_MASS <- 30.973762
MAX_DEPTH_M  <- 10   # epilimnion filter for in-lake nutrients / surface temp
R2_MIN       <- 0.90 # Kd quality filter

# =============================================================================
# SECTION 1: Load and derive each physical driver (native resolution)
# =============================================================================

# ---- 1a. Surface temperature (LTP CTD, 0-10 m mean) -------------------------
message("Loading LTP CTD profiles for surface temperature...")
ctd_profiles <- load_ltp_ctd_profiles(proj_root)
surface_temp <- ctd_profiles %>%
  filter(depth_m <= MAX_DEPTH_M, !is.na(Temperature)) %>%
  group_by(date) %>%
  summarise(value = mean(Temperature, na.rm = TRUE), .groups = "drop") %>%
  transmute(date, variable = "surface_temp_c", value)

# ---- 1b. PM2.5 (Tahoe-wide compiled daily record) ---------------------------
pm25 <- read_csv(file.path(PROC_DIR, "tahoe_pm25_compiled.csv"), show_col_types = FALSE) %>%
  transmute(date = as.Date(date), value = suppressWarnings(as.numeric(PM25))) %>%
  filter(!is.na(value)) %>%
  transmute(date, variable = "pm25_ugm3", value)

# ---- 1c. Upper Truckee River discharge (USGS) -------------------------------
truckee_discharge <- read_csv(file.path(PROC_DIR, "upper_truckee_discharge_daily.csv"),
                              show_col_types = FALSE) %>%
  transmute(date = as.Date(date), value = suppressWarnings(as.numeric(discharge_cms))) %>%
  filter(!is.na(value)) %>%
  group_by(date) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  transmute(date, variable = "upper_truckee_discharge_cms", value)

# ---- 1d. Deposition nutrients (all years, all months) ----------------------
dep_raw <- read_csv(
  file.path(DATA_DIR, "lake_environmental_data", "deposition",
            "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(Start_Datetime), !is.na(End_Datetime), QA_Code == "OK") %>%
  mutate(
    Start_Date = as.Date(Start_Datetime),
    End_Date   = as.Date(End_Datetime),
    date       = Start_Date + as.integer(End_Date - Start_Date) %/% 2L
  )

din_srp_deposition <- dep_raw %>%
  mutate(
    DIN_v = suppressWarnings(as.numeric(NO3_Daily_Load)) +
            suppressWarnings(as.numeric(NH4_Daily_Load)),
    SRP_v = suppressWarnings(as.numeric(SRP_Daily_Load))
  ) %>%
  filter(!is.na(DIN_v), !is.na(SRP_v), SRP_v > 0) %>%
  transmute(date, value = (DIN_v / N_ATOMIC_MASS) / (SRP_v / P_ATOMIC_MASS)) %>%
  transmute(date, variable = "din_srp_deposition_molar", value)

dep_nut_map <- c(NO3 = "NO3_Daily_Load", NH4 = "NH4_Daily_Load",
                 TKN = "TKN_Daily_Load", SRP = "SRP_Daily_Load")
dep_individual <- imap(dep_nut_map, function(col, nut) {
  dep_raw %>%
    transmute(date, value = suppressWarnings(as.numeric(.data[[col]]))) %>%
    filter(!is.na(value)) %>%
    transmute(date, variable = paste0(tolower(nut), "_deposition_mgm2d"), value)
}) %>% bind_rows()

# ---- 1e. In-lake nutrients (MLTP, 0-10 m, all years) -----------------------
lake_nut <- read_csv(
  file.path(DATA_DIR, "lake_environmental_data", "nutrients", "Tahoe_MLTP_Nutrient.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(Date)) %>%
  filter(Depth <= MAX_DEPTH_M)

din_srp_inlake <- lake_nut %>%
  mutate(
    DIN_v = suppressWarnings(as.numeric(NO3)) + suppressWarnings(as.numeric(NH4)),
    SRP_v = suppressWarnings(as.numeric(TRP))
  ) %>%
  filter(!is.na(DIN_v), !is.na(SRP_v), SRP_v > 0) %>%
  group_by(date) %>%
  summarise(value = mean((DIN_v / N_ATOMIC_MASS) / (SRP_v / P_ATOMIC_MASS), na.rm = TRUE),
            .groups = "drop") %>%
  transmute(date, variable = "din_srp_inlake_molar", value)

lake_nut_map <- c(NO3 = "NO3", NH4 = "NH4", TKN = "TKN", SRP = "TRP")
lake_individual <- imap(lake_nut_map, function(col, nut) {
  lake_nut %>%
    transmute(date, value = suppressWarnings(as.numeric(.data[[col]]))) %>%
    filter(!is.na(value)) %>%
    group_by(date) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    transmute(date, variable = paste0(tolower(nut), "_inlake_ugl"), value)
}) %>% bind_rows()

# ---- 1f. 1% PAR depth (LTP Kd_PAR regression results) -----------------------
kd_raw <- read_csv(
  file.path(DATA_DIR, "lake_environmental_data", "uv", "LTP_Kd_Thermocline_depth_results.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(Date))

par_1pct_depth <- kd_raw %>%
  filter(coalesce(R_PAR, 0) >= R2_MIN, coalesce(Kd_PAR, 0) > 0) %>%
  transmute(date, value = 4.605 / Kd_PAR) %>%
  filter(is.finite(value)) %>%
  transmute(date, variable = "par_1pct_depth_m", value)

# ---- 1g. Integrated Chl-a (LTP, trapezoidal, 0-100 m) ----------------------
# Follows the same historical + recent merge as scripts/figure_4_phytoplankton_by_depth.R
# so that "the dataset used for figure 4" also covers the 2011-2012 window.
chl_ltp_recent <- read_csv(
  file.path(DATA_DIR, "lake_environmental_data", "chla", "Tahoe_LTP_Chl.csv"),
  show_col_types = FALSE
) %>%
  transmute(station = "LTP", date = as.Date(Date), depth = as.numeric(Depth),
            chla = as.numeric(Chla)) %>%
  filter(!is.na(chla), !is.na(depth))

chl_all_raw <- tryCatch(
  read.csv(file.path(DATA_DIR, "lake_environmental_data", "chla", "terc_chla_all.csv"),
           quote = "", stringsAsFactors = FALSE),
  error = function(e) read.csv(
    file.path(DATA_DIR, "lake_environmental_data", "chla", "terc_chla_all.csv"),
    quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = ""
  )
)

chl_hist_ltp <- chl_all_raw %>%
  as_tibble() %>%
  filter(Sample_Type %in% c("FIELD", "FLDDUP"), Station_ID == "Index") %>%
  transmute(station = "LTP", date = as.Date(Date),
            depth = suppressWarnings(as.numeric(Depth)),
            chla = suppressWarnings(as.numeric(Chla))) %>%
  filter(!is.na(chla), !is.na(depth), chla > 0)

chl_ltp <- bind_rows(
  anti_join(chl_hist_ltp, chl_ltp_recent, by = c("station", "date", "depth")),
  chl_ltp_recent
) %>%
  filter(depth <= 100, chla < 50)

trap_int <- function(depth, value) {
  ord <- order(depth); d <- depth[ord]; v <- value[ord]
  ok  <- !is.na(v);    d <- d[ok];     v <- v[ok]
  if (length(d) < 2) return(NA_real_)
  sum(diff(d) * (head(v, -1) + tail(v, -1)) / 2)
}

integrated_chla <- chl_ltp %>%
  group_by(date) %>%
  summarise(value = trap_int(depth, chla), .groups = "drop") %>%
  filter(!is.na(value)) %>%
  transmute(date, variable = "integrated_chla_mg_m2", value)

# ---- 1h. Stratification (MLTP Schmidt stability index) ---------------------
schmidt <- read_csv(file.path(PROC_DIR, "ctd", "mltp_schmidt_stability_2005_2025.csv"),
                    show_col_types = FALSE) %>%
  filter(deep_enough) %>%
  transmute(date = as.Date(date), value = schmidt_stability_j_m2) %>%
  filter(!is.na(value)) %>%
  transmute(date, variable = "schmidt_stability_j_m2", value)

# =============================================================================
# SECTION 2: Compile and save
# =============================================================================

physical_drivers_compiled <- bind_rows(
  surface_temp, pm25, truckee_discharge, din_srp_deposition,
  din_srp_inlake, par_1pct_depth, integrated_chla, schmidt,
  dep_individual, lake_individual
) %>%
  arrange(variable, date)

write_csv(physical_drivers_compiled, file.path(PROC_DIR, "physical_drivers_compiled.csv"))
message("Saved: ", file.path(PROC_DIR, "physical_drivers_compiled.csv"),
        " (", nrow(physical_drivers_compiled), " rows)")

physical_drivers_monthly <- physical_drivers_compiled %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(variable, year, month) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(month_date = as.Date(sprintf("%04d-%02d-15", year, month)))

write_csv(physical_drivers_monthly, file.path(PROC_DIR, "physical_drivers_monthly.csv"))
message("Saved: ", file.path(PROC_DIR, "physical_drivers_monthly.csv"),
        " (", nrow(physical_drivers_monthly), " rows)")

# ---- Leptolyngbya abundance by depth zone -----------------------------------
# Monthly means are used so the abundance panel has the same time resolution
# as the physical-driver anomaly panels.
lepto_abundance_monthly <- read_csv(
  file.path(DATA_DIR, "lake_environmental_data", "phytoplankton",
            "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    date = floor_date(as.Date(date), "month") + days(14),
    depth_zone = depth_bin,
    taxon = taxon,
    value = suppressWarnings(as.numeric(abundance))
  ) %>%
  filter(taxon == "Leptolyngbya spp.",
         depth_zone %in% c("0-40 m", "60-105 m"), !is.na(value)) %>%
  group_by(date, depth_zone) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

LEPTO_DEPTH_COLS <- c("0-40 m" = "#1B9E77", "60-105 m" = "#66C2A5")

# =============================================================================
# SECTION 3: Variable specification (title, units, colors, column assignment)
# =============================================================================

DEP_UNIT  <- "mg m^-2 d^-1"
LAKE_UNIT <- "\u00b5g L^-1"

var_spec <- tribble(
  ~variable,                      ~title,                                ~ylab_val,                          ~ylab_anom,                       ~main_col,  ~high_col,  ~low_col,   ~col,
  "surface_temp_c",               "Surface temperature (0-10 m)",        "Temperature (\u00b0C)",           "Temp. anomaly (\u00b0C)",         "#DE4E6A",  "#DE4E6A",  "#5698E8",  1,
  "pm25_ugm3",                    "PM2.5",                               "PM2.5 (\u00b5g m^-3)",            "PM2.5 anomaly (\u00b5g m^-3)",    "#5831B8",  "#5831B8",  "#A89FF3",  1,
  "upper_truckee_discharge_cms",  "Upper Truckee River discharge",       "Discharge (m^3 s^-1)",            "Discharge anomaly (m^3 s^-1)",    "#BC5090",  "#BC5090",  "#FFA600",  1,
  "din_srp_deposition_molar",     "Deposition DIN:SRP (molar)",          "DIN:SRP (molar)",                 "DIN:SRP anomaly",                 "#E65100",  "#E65100",  "#FFAB76",  1,
  "din_srp_inlake_molar",         "In-lake DIN:TRP (molar, 0-10 m)",     "DIN:TRP (molar)",                 "DIN:TRP anomaly",                 "#1565C0",  "#1565C0",  "#7CB9E8",  1,
  "par_1pct_depth_m",             "1% PAR depth",                        "1% PAR depth (m)",                "1% PAR depth anomaly (m)",        "#2E7D32",  "#2E7D32",  "#A5D6A7",  1,
  "integrated_chla_mg_m2",        "Integrated Chl-a (0-100 m)",          "Chl-a (mg m^-2)",                 "Chl-a anomaly (mg m^-2)",          "#00695C",  "#00695C",  "#80CBC4",  1,
  "schmidt_stability_j_m2",       "Stratification (Schmidt stability)",  "Stability (J m^-2)",              "Stability anomaly (J m^-2)",       "#5E35B1",  "#5E35B1",  "#D1C4E9",  1,
  "no3_deposition_mgm2d",         "Deposition NO3",                      paste0("NO3 (", DEP_UNIT, ")"),    "NO3 anomaly",                      "#E65100",  "#E65100",  "#FFCC80",  2,
  "no3_inlake_ugl",               "In-lake NO3 (0-10 m)",                paste0("NO3 (", LAKE_UNIT, ")"),   "NO3 anomaly",                      "#1565C0",  "#1565C0",  "#90CAF9",  2,
  "nh4_deposition_mgm2d",         "Deposition NH4",                      paste0("NH4 (", DEP_UNIT, ")"),    "NH4 anomaly",                      "#BF360C",  "#BF360C",  "#FFAB91",  2,
  "nh4_inlake_ugl",               "In-lake NH4 (0-10 m)",                paste0("NH4 (", LAKE_UNIT, ")"),   "NH4 anomaly",                      "#0D47A1",  "#0D47A1",  "#82B1FF",  2,
  "tkn_deposition_mgm2d",         "Deposition TKN",                      paste0("TKN (", DEP_UNIT, ")"),    "TKN anomaly",                      "#F57F17",  "#F57F17",  "#FFE082",  2,
  "tkn_inlake_ugl",               "In-lake TKN (0-10 m)",                paste0("TKN (", LAKE_UNIT, ")"),   "TKN anomaly",                      "#01579B",  "#01579B",  "#81D4FA",  2,
  "srp_deposition_mgm2d",         "Deposition SRP",                      paste0("SRP (", DEP_UNIT, ")"),    "SRP anomaly",                      "#AD1457",  "#AD1457",  "#F48FB1",  2,
  "srp_inlake_ugl",               "In-lake TRP (0-10 m)",                paste0("TRP (", LAKE_UNIT, ")"),   "TRP anomaly",                      "#4527A0",  "#4527A0",  "#B39DDB",  2
)

# Anomaly panels show the same unit as the absolute value panels (no "anomaly" label)
var_spec <- var_spec %>% mutate(ylab_anom = ylab_val)

# =============================================================================
# SECTION 4: Panel builders
# =============================================================================

add_event_shading <- function(p, show_caldor, show_lepto_bloom1) {
  if (show_lepto_bloom1) {
    p <- p +
      annotate("rect", xmin = LEPTO_BLOOM1_START, xmax = LEPTO_BLOOM1_END,
               ymin = -Inf, ymax = Inf, fill = LEPTO_COL, alpha = 0.10) +
      geom_vline(xintercept = LEPTO_BLOOM1_START, linetype = "dashed",
                 colour = LEPTO_COL, linewidth = LO_LW_MID)
  }
  if (show_caldor) {
    p <- p +
      annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
               ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
      geom_vline(xintercept = CALDOR_START, linetype = "dashed",
                 colour = "firebrick", linewidth = LO_LW_MID)
  }
  p
}

#' Add black dashed vertical lines for the end of each Leptolyngbya bloom
#' when the date falls inside the current x-axis window.
add_bloom_end_lines <- function(p, window_start, window_end) {
  term_dates <- c(LEPTO_BLOOM1_END, LEPTO_BLOOM2_END)
  in_window <- term_dates >= window_start & term_dates <= window_end
  for (d in term_dates[in_window]) {
    p <- p + geom_vline(xintercept = d, linetype = "dashed",
                        colour = "black", linewidth = LO_LW_MID)
  }
  p
}

base_panel_theme <- function(show_x) {
  th <- theme_bw(base_size = 8 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
    theme(panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 8 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY),
          axis.title.y = element_text(size = 7 * LO_FIG_SCALE, family = BASE_FAMILY),
          axis.text   = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY))
  if (!show_x) th <- th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  th
}

#' Monthly anomaly relative to a calendar-month climatology. When
#' `exclude_years` is supplied, those years are dropped from the climatology
#' baseline (used for the two bloom-window figures so the baseline doesn't
#' include the window itself). For the full-record figure, `exclude_years`
#' is NULL and the climatology uses the entire available record.
build_monthly_anomaly <- function(var_id, exclude_years = NULL) {
  dat <- physical_drivers_monthly %>% filter(variable == var_id)
  clim <- dat %>%
    filter(is.null(exclude_years) | !year %in% exclude_years) %>%
    group_by(month) %>%
    summarise(clim_mean = mean(value, na.rm = TRUE), .groups = "drop")
  dat %>%
    left_join(clim, by = "month") %>%
    mutate(anomaly = value - clim_mean)
}

make_anomaly_panel <- function(var_id, title, ylab, high_col, low_col,
                               window_start, window_end, exclude_years,
                               show_x, show_caldor, show_lepto_bloom1) {
  anom <- build_monthly_anomaly(var_id, exclude_years) %>%
    filter(month_date >= floor_date(window_start, "month"),
           month_date <= ceiling_date(window_end, "month") - days(1)) %>%
    filter(!is.na(anomaly))

  p <- ggplot(anom, aes(month_date, anomaly))
  p <- add_event_shading(p, show_caldor, show_lepto_bloom1)
  p <- add_bloom_end_lines(p, window_start, window_end)

  if (nrow(anom) == 0) {
    p <- p +
      annotate("text", x = mean(c(window_start, window_end)), y = 0,
               label = "No data", colour = "grey45", size = lo_geom_text_size())
  } else {
    p <- p +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = LO_LW_THIN) +
      geom_ribbon(aes(ymin = 0, ymax = pmax(anomaly, 0)), fill = high_col,
                  alpha = 0.65, colour = NA, na.rm = TRUE) +
      geom_ribbon(aes(ymin = pmin(anomaly, 0), ymax = 0), fill = low_col,
                  alpha = 0.65, colour = NA, na.rm = TRUE) +
      geom_line(colour = "grey25", linewidth = LO_LW_THIN, na.rm = TRUE)
  }

  date_breaks <- if (as.numeric(difftime(window_end, window_start, units = "days")) > 3000) {
    "2 years"
  } else {
    "2 months"
  }
  date_labels <- if (date_breaks == "2 years") "%Y" else "%b %Y"

  p <- p +
    scale_x_date(limits = c(window_start, window_end), date_breaks = date_breaks,
                 date_labels = date_labels, expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = ylab, title = title) +
    base_panel_theme(show_x)
  p
}

make_value_panel <- function(var_id, title, ylab, main_col,
                             window_start, window_end,
                             show_x, show_caldor, show_lepto_bloom1) {
  dat <- physical_drivers_compiled %>%
    filter(variable == var_id, date >= window_start, date <= window_end) %>%
    filter(!is.na(value))

  p <- ggplot(dat, aes(date, value))
  p <- add_event_shading(p, show_caldor, show_lepto_bloom1)
  p <- add_bloom_end_lines(p, window_start, window_end)

  if (nrow(dat) == 0) {
    p <- p +
      annotate("text", x = mean(c(window_start, window_end)), y = 0,
               label = "No data", colour = "grey45", size = lo_geom_text_size())
  } else {
    p <- p +
      geom_line(colour = main_col, linewidth = LO_LW_MID, alpha = 0.85, na.rm = TRUE)
    if (var_id %in% c("din_srp_deposition_molar", "din_srp_inlake_molar")) {
      p <- p +
        geom_hline(yintercept = 16, linetype = "dashed", colour = "black",
                   linewidth = LO_LW_THIN) +
        annotate("text", x = window_start + 10, y = 16,
                 label = "Redfield ratio", hjust = 0, vjust = -0.7,
                 colour = "black", size = lo_geom_text_size(), family = BASE_FAMILY)
    }
  }

  date_breaks <- if (as.numeric(difftime(window_end, window_start, units = "days")) > 3000) {
    "2 years"
  } else {
    "2 months"
  }
  date_labels <- if (date_breaks == "2 years") "%Y" else "%b %Y"

  p <- p +
    scale_x_date(limits = c(window_start, window_end), date_breaks = date_breaks,
                 date_labels = date_labels, expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = ylab, title = title) +
    base_panel_theme(show_x)
  p
}

make_lepto_panel <- function(window_start, window_end, show_x) {
  dat <- lepto_abundance_monthly %>%
    filter(date >= floor_date(window_start, "month"),
           date <= ceiling_date(window_end, "month") - days(1))

  p <- ggplot(dat, aes(date, value, colour = depth_zone, group = depth_zone)) +
    geom_line(linewidth = LO_LW_MID, na.rm = TRUE) +
    geom_point(size = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = LEPTO_DEPTH_COLS, name = "Depth zone") +
    scale_x_date(limits = c(window_start, window_end),
                 date_breaks = "2 months", date_labels = "%b %Y",
                 expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = expression("Abundance (cells L"^-1*")"),
         title = "Leptolyngbya abundance") +
    base_panel_theme(show_x) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY),
          legend.text = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY),
          legend.key.size = unit(2.5 * LO_FIG_SCALE, "mm"))
  p
}

# =============================================================================
# SECTION 5: Figure builders (2-column, 8-row layouts)
# =============================================================================

n_col1 <- sum(var_spec$col == 1)
n_col2 <- sum(var_spec$col == 2)

build_anomaly_figure <- function(window_start, window_end, exclude_years,
                                 fname, show_caldor, show_lepto_bloom1, title,
                                 add_lepto = FALSE) {
  ordered_vars <- bind_rows(
    var_spec %>% filter(col == 1),
    var_spec %>% filter(col == 2)
  )
  panels <- pmap(ordered_vars, function(variable, title, ylab_val, ylab_anom,
                                        main_col, high_col, low_col, col) {
    is_last_in_col <- variable == (ordered_vars %>% filter(col == !!col) %>% pull(variable) %>% tail(1))
    make_anomaly_panel(
      variable, title, ylab_anom, high_col, low_col,
      window_start = window_start, window_end = window_end,
      exclude_years = exclude_years,
      show_x = is_last_in_col, show_caldor = show_caldor,
      show_lepto_bloom1 = show_lepto_bloom1
    )
  })
  if (add_lepto) panels <- append(panels, list(make_lepto_panel(window_start, window_end, TRUE)))

  fig <- wrap_plots(panels, ncol = 2, byrow = FALSE) +
    plot_annotation(
      title = title,
      tag_levels = "a",
      theme = theme(plot.title = element_text(size = 9 * LO_FIG_SCALE, face = "bold",
                                              family = BASE_FAMILY),
                    plot.tag = element_text(size = LO_FS_TAG, face = "bold", family = BASE_FAMILY))
    )

  out_path <- file.path(OUT_DIR, fname)
  ragg::agg_png(out_path, width = lo_word_width(7.62) * 1.7, height = 24, units = "cm", res = LO_DPI)
  print(fig)
  invisible(dev.off())
  message("Saved: ", out_path)
}

build_value_figure <- function(window_start, window_end, fname,
                               show_caldor, show_lepto_bloom1, title,
                               add_lepto = FALSE) {
  ordered_vars <- bind_rows(
    var_spec %>% filter(col == 1),
    var_spec %>% filter(col == 2)
  )
  panels <- pmap(ordered_vars, function(variable, title, ylab_val, ylab_anom,
                                        main_col, high_col, low_col, col) {
    is_last_in_col <- variable == (ordered_vars %>% filter(col == !!col) %>% pull(variable) %>% tail(1))
    make_value_panel(
      variable, title, ylab_val, main_col,
      window_start = window_start, window_end = window_end,
      show_x = is_last_in_col, show_caldor = show_caldor,
      show_lepto_bloom1 = show_lepto_bloom1
    )
  })
  if (add_lepto) panels <- append(panels, list(make_lepto_panel(window_start, window_end, TRUE)))

  fig <- wrap_plots(panels, ncol = 2, byrow = FALSE) +
    plot_annotation(
      title = title,
      tag_levels = "a",
      theme = theme(plot.title = element_text(size = 9 * LO_FIG_SCALE, face = "bold",
                                              family = BASE_FAMILY),
                    plot.tag = element_text(size = LO_FS_TAG, face = "bold", family = BASE_FAMILY))
    )

  out_path <- file.path(OUT_DIR, fname)
  ragg::agg_png(out_path, width = lo_word_width(7.62) * 1.7, height = 24, units = "cm", res = LO_DPI)
  print(fig)
  invisible(dev.off())
  message("Saved: ", out_path)
}

# ---- Window 1: first Leptolyngbya bloom (no Caldor Fire overlap) -----------
build_anomaly_figure(
  window_start = as.Date("2011-01-01"), window_end = as.Date("2012-04-30"),
  exclude_years = c(2011, 2012),
  fname = "supplemental_S6_physical_drivers_anomaly_2011_2012.png",
  show_caldor = FALSE, show_lepto_bloom1 = TRUE,
  title = "Physical/chemical driver anomalies: Jan 2011 - Apr 2012",
  add_lepto = TRUE
)
build_value_figure(
  window_start = as.Date("2011-01-01"), window_end = as.Date("2012-04-30"),
  fname = "supplemental_S6b_physical_drivers_values_2011_2012.png",
  show_caldor = FALSE, show_lepto_bloom1 = TRUE,
  title = "Physical/chemical driver values: Jan 2011 - Apr 2012",
  add_lepto = TRUE
)

# ---- Window 2: second Leptolyngbya bloom (Caldor Fire window shaded) ------
build_anomaly_figure(
  window_start = as.Date("2021-04-01"), window_end = as.Date("2022-07-31"),
  exclude_years = c(2021, 2022),
  fname = "supplemental_S7_physical_drivers_anomaly_2021_2022.png",
  show_caldor = TRUE, show_lepto_bloom1 = FALSE,
  title = "Physical/chemical driver anomalies: Apr 2021 - Jul 2022",
  add_lepto = TRUE
)
build_value_figure(
  window_start = as.Date("2021-04-01"), window_end = as.Date("2022-07-31"),
  fname = "supplemental_S7b_physical_drivers_values_2021_2022.png",
  show_caldor = TRUE, show_lepto_bloom1 = FALSE,
  title = "Physical/chemical driver values: Apr 2021 - Jul 2022",
  add_lepto = TRUE
)

# ---- Full record: 2005-2025 anomaly, deseasonalized against full record ---
# All 16 variable panels are always drawn, even when a variable's native
# record only partially overlaps 2005-2025 (e.g., deposition nutrients begin
# in 2013); the panel simply shows data over its available span.
build_anomaly_figure(
  window_start = as.Date("2005-01-01"), window_end = as.Date("2025-12-31"),
  exclude_years = NULL,
  fname = "supplemental_S8_physical_drivers_anomaly_2005_2025.png",
  show_caldor = TRUE, show_lepto_bloom1 = TRUE,
  title = "Physical/chemical driver anomalies: 2005-2025"
)
build_value_figure(
  window_start = as.Date("2005-01-01"), window_end = as.Date("2025-12-31"),
  fname = "supplemental_S8b_physical_drivers_values_2005_2025.png",
  show_caldor = TRUE, show_lepto_bloom1 = TRUE,
  title = "Physical/chemical driver values: 2005-2025"
)

message("\nAll physical driver outputs saved to: ", OUT_DIR)
