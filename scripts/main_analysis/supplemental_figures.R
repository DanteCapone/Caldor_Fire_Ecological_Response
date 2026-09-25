# =============================================================================
# supplemental_figures.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Supplemental figures for manuscript
#
#   S1: Leptolyngbya and Cryptomonas depth profiles Aug 2021-Apr 2022
#       Separate compact 3 x 3 figures (nine monthly panels per taxon)
#
#   S2: Nutrient time series - deposition and in-lake, all years
#       NO3, NH4, SRP, TP, TKN for 2021 with historical reference
#
#   S3: Cross-correlation heatmap: phytoplankton taxa abundance vs. PM2.5
#       Taxa on y-axis, lag (months) on x-axis, ccf coefficient as color
#
# Output: figures/supplemental/
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")

# ---- Paths -------------------------------------------------------------------
PHYTO_FILE <- file.path("data", "lake_environmental_data", "phytoplankton",
                        "PhytoData_LTP_MLTP_2019-2025_counts_biovolume.xlsx")
NUTRI_MLTP <- file.path("data", "lake_environmental_data", "nutrients",
                        "Tahoe_MLTP_Nutrient.csv")
NUTRI_LTP  <- file.path("data", "lake_environmental_data", "nutrients",
                        "Tahoe_LTP_Nutrient.csv")
DEPO_FILE  <- file.path("data", "lake_environmental_data", "deposition",
                        "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv")
PM25_FILE  <- file.path("data", "processed", "tahoe_hms_pm25_daily.csv")
OUT_DIR    <- file.path("figures", "supplemental")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY  <- "Times New Roman"
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")
FIRE_YEAR    <- 2021
SEASON_MONTHS <- 8:10
MAX_DEPTH_M  <- 10

COL_LEPTO <- "#009E73"    # Leptolyngbya - teal/green
COL_CRYPT <- "#9467BD"    # Cryptomonas  - purple
COL_DEP_HIST <- "#FFAB76"
COL_DEP_2021 <- "#E65100"
COL_LAK_HIST <- "#7CB9E8"
COL_LAK_2021 <- "#1565C0"

# =============================================================================
# -- Load phytoplankton data --------------------------------------------------
# =============================================================================
parse_excel_date <- function(x) {
  x <- trimws(as.character(x))
  d_iso    <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  n        <- suppressWarnings(as.numeric(x))
  d_serial <- as.Date(n, origin = as.Date("1899-12-30"))
  coalesce(d_iso, d_serial)
}

ltp_raw <- read_excel(PHYTO_FILE, sheet = "LTP_2019-2025",
                      col_types = c("text","text","text","text","text",
                                    "numeric","text","text","numeric","text","numeric"))

clean_ltp <- ltp_raw %>%
  rename(station_id  = Station_ID,
         event_id    = Event_ID,
         depth_m     = `Depth_(m)`,
         phylum      = Phylum,
         taxon       = Taxon,
         abundance   = `Abundance_(units/L)`,
         unit_biovol = `Unit_biovolume_(µm3)`) %>%
  mutate(date        = parse_excel_date(Date),
         year        = year(date),
         month       = month(date),
         month_date  = floor_date(date, "month"),
         depth_m     = as.numeric(depth_m),
         abundance   = as.numeric(abundance),
         unit_biovol = as.numeric(unit_biovol),
         biovolume   = abundance * unit_biovol) %>%
  filter(!is.na(abundance), !is.na(depth_m), !is.na(date))

mltp_raw <- read_excel(PHYTO_FILE, sheet = "MLTP_2019-2025",
                       col_types = c("text", "text", "text", "text", "text",
                                     "text", "text", "text", "numeric",
                                     "text", "numeric"))

clean_mltp <- mltp_raw %>%
  rename(station_id  = Station_ID,
         event_id    = Event_ID,
         depth_m     = all_of("Collection_depth_range_(m)"),
         phylum      = Phylum,
         taxon       = Taxon,
         abundance   = all_of("Abundance_(units/L)")) %>%
  rename_with(~ "unit_biovol", matches("Unit_biovolume")) %>%
  mutate(date        = parse_excel_date(Date),
         year        = year(date),
         month       = month(date),
         month_date  = floor_date(date, "month"),
         depth_m     = as.character(depth_m),
         abundance   = as.numeric(abundance),
         unit_biovol = as.numeric(unit_biovol),
         biovolume   = abundance * unit_biovol) %>%
  filter(!is.na(abundance), !is.na(depth_m), !is.na(date))

all_taxa_levels <- sort(unique(c(clean_ltp$taxon, clean_mltp$taxon)))
taxon_all_cols <- setNames(hcl.colors(length(all_taxa_levels), palette = "Dark 3"), all_taxa_levels)
taxon_all_cols[intersect(names(taxon_all_cols), c("Leptolyngbya sp."))] <- COL_LEPTO
taxon_all_cols[intersect(names(taxon_all_cols), c("Cryptomonas sp."))] <- COL_CRYPT
taxon_all_cols[intersect(names(taxon_all_cols), c("Cyclotella sp."))] <- "#0072B2"
taxon_all_cols[intersect(names(taxon_all_cols), c("Synedra acus var. radians"))] <- "#D55E00"

cat("Phytoplankton rows:", nrow(clean_ltp),
    "| dates:", n_distinct(clean_ltp$date), "\n")

# =============================================================================
# -- S1: Leptolyngbya + Cryptomonas profiles, Aug 2021-Apr 2022 --------------
# =============================================================================
S1_MONTHS <- tibble(
  year = c(rep(2021L, 5), rep(2022L, 4)),
  month = c(8L, 9L, 10L, 11L, 12L, 1L, 2L, 3L, 4L)
)
S1_WIDTH_CM <- 8.5
S1_HEIGHT_CM <- 12.0

make_taxon_depth_panel <- function(taxon_name, taxon_col,
                                   target_year, mon_num, show_y = TRUE,
                                   show_x_title = TRUE,
                                   panel_title = NULL) {
  if (is.null(panel_title)) panel_title <- paste(month.abb[mon_num], target_year)
  dat_obs <- clean_ltp %>%
    filter(taxon == taxon_name, !is.na(depth_m), depth_m <= 105,
           year == target_year, month == mon_num) %>%
    mutate(depth_zone = factor(if_else(depth_m <= 40, "0-40 m", "60-105 m"),
                               levels = c("0-40 m", "60-105 m"))) %>%
    group_by(depth_m, depth_zone) %>%
    summarise(obs_mean = mean(abundance, na.rm = TRUE), .groups = "drop")

  if (nrow(dat_obs) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No sample", colour = "grey45") +
        labs(title = panel_title) +
        theme_void(base_size = 7, base_family = BASE_FAMILY) +
        theme(plot.title = element_text(face = "bold", size = 7))
    )
  }

  p <- ggplot(dat_obs, aes(x = obs_mean, y = depth_m)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25 * LO_FIG_SCALE) +
    geom_path(colour = taxon_col, linewidth = LO_LW_THICK, na.rm = TRUE) +
    geom_point(colour = taxon_col, fill = taxon_col, shape = 21,
               size = 1.25 * LO_FIG_SCALE, stroke = 0.35 * LO_FIG_SCALE, na.rm = TRUE) +
    scale_y_reverse(breaks = c(0, 5, 10, 20, 30, 40, 50, 60, 75, 90, 105)) +
    scale_x_continuous(labels = label_scientific(digits = 2),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = panel_title,
         y = "Depth (m)",
         x = expression("Abundance (cells L"^-1*")")) +
    theme_classic(base_size = 7, base_family = BASE_FAMILY) +
    theme(axis.title = element_text(size = 6.5, family = BASE_FAMILY),
          axis.text.x = element_text(size = 6, angle = 35, hjust = 1,
                                     family = BASE_FAMILY),
          axis.text.y = element_text(size = 6, family = BASE_FAMILY),
          plot.title = element_text(face = "bold", size = 7,
                                    family = BASE_FAMILY),
          panel.grid.major.y = element_line(colour = "#EBEBEB",
                                            linewidth = LO_LW_THIN),
          plot.margin = margin(0.8, 0.8, 0.8, 0.8, "mm"))

  if (!show_y)      p <- p + theme(axis.title.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
  if (!show_x_title) p <- p + theme(axis.title.x = element_blank())
  p
}
# NOTE: taxon strings in the data carry the " sp." suffix.
TAXON_LEPTO <- "Leptolyngbya sp."
TAXON_CRYPT <- "Cryptomonas sp."

# Leptolyngbya monthly panels
lepto_panels <- pmap(
  list(S1_MONTHS$year, S1_MONTHS$month, seq_len(nrow(S1_MONTHS))),
  function(yr, m, i) {
    make_taxon_depth_panel(
      TAXON_LEPTO, COL_LEPTO, yr, m,
      show_y       = (i %% 3 == 1),
      show_x_title = (i == 8),
      panel_title  = paste(month.abb[m], yr)
    )
  }
)

# Cryptomonas monthly panels
crypt_panels <- pmap(
  list(S1_MONTHS$year, S1_MONTHS$month, seq_len(nrow(S1_MONTHS))),
  function(yr, m, i) {
    make_taxon_depth_panel(
      TAXON_CRYPT, COL_CRYPT, yr, m,
      show_y       = (i %% 3 == 1),
      show_x_title = (i == 8),
      panel_title  = paste(month.abb[m], yr)
    )
  }
)

# Separate compact 3 x 3 canvases keep the full nine-month sequence readable.
unlink(file.path(OUT_DIR, "supplemental_S1_depth_profiles.png"))
fig_s1a <- wrap_plots(lepto_panels, ncol = 3, nrow = 3, byrow = TRUE) +
  plot_annotation(
    title      = "Leptolyngbya sp.",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 8, face = "bold",
                                             family = BASE_FAMILY),
                  plot.tag   = element_text(size = 7, face = "bold",
                                            family = BASE_FAMILY),
                  plot.margin = margin(2.5, 2.5, 2.0, 2.5, "mm"))
  )
fig_s1b <- wrap_plots(crypt_panels, ncol = 3, nrow = 3, byrow = TRUE) +
  plot_annotation(
    title      = "Cryptomonas sp.",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 8, face = "bold",
                                             family = BASE_FAMILY),
                  plot.tag   = element_text(size = 7, face = "bold",
                                            family = BASE_FAMILY),
                  plot.margin = margin(2.5, 2.5, 2.0, 2.5, "mm"))
  )

ragg::agg_png(file.path(OUT_DIR, "supplemental_S1a_depth_profiles_leptolyngbya.png"),
              width = S1_WIDTH_CM, height = S1_HEIGHT_CM, units = "cm", res = LO_DPI)
print(fig_s1a)
invisible(dev.off())
message("Saved: supplemental_S1a_depth_profiles_leptolyngbya.png")

ragg::agg_png(file.path(OUT_DIR, "supplemental_S1b_depth_profiles_cryptomonas.png"),
              width = S1_WIDTH_CM, height = S1_HEIGHT_CM, units = "cm", res = LO_DPI)
print(fig_s1b)
invisible(dev.off())
message("Saved: supplemental_S1b_depth_profiles_cryptomonas.png")

# =============================================================================
# -- S2: Nutrient time series 2021 -------------------------------------------
#    Two rows: Atmospheric deposition (top) and In-lake 0-10 m (bottom).
#    Columns = nutrients. Each panel: 2021 monthly observed on the LEFT axis;
#    monthly climatology (mean +/- 95% CI) on the RIGHT axis (similar to Fig 4).
# =============================================================================
lake_2021 <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Year == FIRE_YEAR)

dep_2021 <- read_csv(DEPO_FILE, show_col_types = FALSE) %>%
  filter(!is.na(Start_Datetime), !is.na(End_Datetime)) %>%
  mutate(Start_Date = as.Date(Start_Datetime),
         End_Date   = as.Date(End_Datetime),
         Mid_Date   = Start_Date + as.integer(End_Date - Start_Date) %/% 2L,
         Year = year(Mid_Date), Month = month(Mid_Date)) %>%
  filter(Year == FIRE_YEAR)

# Full historical record (all calendar months, all years != 2021) for climatology
lake_climsrc <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date), Year = year(Date), Month = month(Date)) %>%
  filter(Depth <= MAX_DEPTH_M, Year != FIRE_YEAR)

dep_climsrc <- read_csv(DEPO_FILE, show_col_types = FALSE) %>%
  filter(!is.na(Start_Datetime), !is.na(End_Datetime)) %>%
  mutate(Start_Date = as.Date(Start_Datetime),
         End_Date   = as.Date(End_Datetime),
         Mid_Date   = Start_Date + as.integer(End_Date - Start_Date) %/% 2L,
         Year = year(Mid_Date), Month = month(Mid_Date)) %>%
  filter(Year != FIRE_YEAR)

nut_map <- list(
  NO3 = list(lake_col = "NO3", dep_col = "NO3_Daily_Load", disp = "NO\u2083\u207B"),
  NH4 = list(lake_col = "NH4", dep_col = "NH4_Daily_Load", disp = "NH\u2084\u207A"),
  SRP = list(lake_col = "TRP", dep_col = "SRP_Daily_Load",
             disp_dep = "SRP", disp_lake = "TRP"),
  TP  = list(lake_col = "THP", dep_col = "TP_Daily_Load",
             disp_dep = "TP", disp_lake = "THP"),
  TKN = list(lake_col = "TKN", dep_col = "TKN_Daily_Load", disp = "TKN")
)

LAKE_UNIT_S2 <- "\u00b5g L\u207B\u00B9"
DEP_UNIT_S2  <- "mg m\u207B\u00B2 d\u207B\u00B9"

# Monthly climatology (mean +/- 95% CI) from a raw source column
monthly_clim <- function(src_df, col) {
  v <- suppressWarnings(as.numeric(src_df[[col]]))
  tibble(Month = src_df$Month, value = v) %>%
    filter(!is.na(value)) %>%
    group_by(Month) %>%
    summarise(n = n(), clim_mean = mean(value), clim_sd = sd(value),
              .groups = "drop") %>%
    mutate(clim_ci   = qt(0.975, df = pmax(n - 1, 1)) * clim_sd / sqrt(pmax(n, 1)),
           clim_lo   = pmax(clim_mean - clim_ci, 0),
           clim_hi   = clim_mean + clim_ci,
           plot_date = as.Date(paste0(FIRE_YEAR, "-", sprintf("%02d", Month), "-15")))
}

make_nut_panel_s2 <- function(nut_id, source,
                              show_left_title = TRUE, show_right_title = TRUE,
                              show_title = TRUE, show_x = TRUE,
                              compact_axis_titles = FALSE) {
  nm <- nut_map[[nut_id]]

  if (source == "In-lake") {
    col_obs <- COL_LAK_2021; unit <- LAKE_UNIT_S2
    obs <- tibble(Date = lake_2021$Date,
                  value = suppressWarnings(as.numeric(lake_2021[[nm$lake_col]]))) %>%
      filter(!is.na(value)) %>% group_by(Date) %>%
      summarise(val = mean(value), .groups = "drop")
    clim <- monthly_clim(lake_climsrc, nm$lake_col)
  } else {
    col_obs <- COL_DEP_2021; unit <- DEP_UNIT_S2
    obs <- tibble(Date = dep_2021$Mid_Date,
                  val = suppressWarnings(as.numeric(dep_2021[[nm$dep_col]]))) %>%
      filter(!is.na(val))
    clim <- monthly_clim(dep_climsrc, nm$dep_col)
  }

  # Scale climatology (right axis) onto the observed (left axis) coordinate space
  o_rng <- range(obs$val, na.rm = TRUE)
  c_rng <- range(c(clim$clim_lo, clim$clim_hi), na.rm = TRUE)
  if (!all(is.finite(o_rng))) o_rng <- c(0, 1)
  if (!all(is.finite(c_rng))) c_rng <- c(0, 1)
  if (diff(o_rng) < .Machine$double.eps) o_rng <- o_rng + c(-0.5, 0.5)
  if (diff(c_rng) < .Machine$double.eps) c_rng <- c_rng + c(-0.5, 0.5)
  if (source == "In-lake") {
    common_rng <- range(c(o_rng, c_rng), na.rm = TRUE)
    o_rng <- common_rng
    c_rng <- common_rng
  } else {
    c_mid <- mean(c_rng)
    c_half <- diff(c_rng) * 0.65
    c_rng <- c_mid + c(-c_half, c_half)
  }
  slope  <- diff(o_rng) / diff(c_rng)
  intcpt <- o_rng[1] - c_rng[1] * slope

  x_lim <- c(as.Date("2021-01-01"), as.Date("2021-12-31"))

  p <- ggplot() +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.07) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = 0.4 * LO_FIG_SCALE) +
    # Climatology (right axis): CI ribbon + mean line
    geom_ribbon(data = clim,
                aes(x = plot_date, ymin = clim_lo * slope + intcpt,
                    ymax = clim_hi * slope + intcpt),
                fill = "grey70", alpha = 0.40, na.rm = TRUE) +
    geom_line(data = clim, aes(x = plot_date, y = clim_mean * slope + intcpt),
              colour = "grey55", linetype = "dashed", linewidth = 0.6 * LO_FIG_SCALE,
              na.rm = TRUE) +
    # 2021 observed (left axis): line only
    geom_line(data = obs, aes(x = Date, y = val),
              colour = col_obs, linewidth = 0.80 * LO_FIG_SCALE, na.rm = TRUE) +
    scale_x_date(limits = x_lim, date_breaks = "2 months",
                 date_labels = "%b", expand = expansion(mult = 0.02)) +
    scale_y_continuous(
      name     = if (show_left_title) {
        if (compact_axis_titles) "2021" else paste0("2021 (", unit, ")")
      } else NULL,
      labels   = label_number(accuracy = NULL),
      expand   = expansion(mult = c(0.05, 0.15)),
      sec.axis = sec_axis(~ (. - intcpt) / slope,
                          name   = if (show_right_title) {
                            if (compact_axis_titles) "Historical" else
                              paste0("Climatology (", unit, ")")
                          } else NULL,
                          labels = label_number(accuracy = NULL))
    ) +
    labs(x = NULL, title = if (show_title) {
      if (source == "In-lake") {
        if (!is.null(nm$disp_lake)) nm$disp_lake else nm$disp
      } else {
        if (!is.null(nm$disp_dep)) nm$disp_dep else nm$disp
      }
    } else NULL) +
    theme_bw(base_size = 9 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
    theme(panel.grid.minor   = element_blank(),
          plot.title         = element_text(size = 9 * LO_FIG_SCALE, face = "bold",
                                             family = BASE_FAMILY),
          axis.title.y       = element_text(colour = col_obs,  size = 7.5 * LO_FIG_SCALE),
          axis.title.y.right = element_text(colour = "grey45", size = 7.5 * LO_FIG_SCALE),
          axis.text.y        = element_text(colour = col_obs,  size = 7 * LO_FIG_SCALE),
          axis.text.y.right  = element_text(colour = "grey45", size = 7 * LO_FIG_SCALE),
          axis.ticks.y       = element_line(colour = col_obs),
          axis.ticks.y.right = element_line(colour = "grey45"))
  if (!show_x)
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

nut_ids <- names(nut_map)
n_nut   <- length(nut_ids)

dep_row_panels <- imap(nut_ids, function(id, j)
  make_nut_panel_s2(id, "Deposition",
                    show_left_title  = (j == 1),
                    show_right_title = (j == n_nut),
                    show_title = TRUE,  show_x = FALSE))

lake_row_panels <- imap(nut_ids, function(id, j)
  make_nut_panel_s2(id, "In-lake",
                    show_left_title  = (j == 1),
                    show_right_title = (j == n_nut),
                    show_title = TRUE, show_x = TRUE))

dep_page_panels <- imap(nut_ids, function(id, j)
  make_nut_panel_s2(id, "Deposition",
                    show_left_title  = (j %% 2 == 1),
                    show_right_title = (j %% 2 == 0 || j == n_nut),
                    show_title = TRUE, show_x = TRUE))

lake_page_panels <- imap(nut_ids, function(id, j)
  make_nut_panel_s2(id, "In-lake",
                    show_left_title  = (j %% 2 == 1),
                    show_right_title = (j %% 2 == 0 || j == n_nut),
                    show_title = TRUE, show_x = TRUE))

fig_s2_dep <- wrap_plots(dep_page_panels, ncol = 2) +
  plot_annotation(
    title      = "Atmospheric deposition",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 9 * LO_FIG_SCALE, face = "bold",
                                             family = BASE_FAMILY),
                  plot.tag   = element_text(size = LO_FS_TAG, face = "bold",
                                            family = BASE_FAMILY)))

fig_s2_lake <- wrap_plots(lake_page_panels, ncol = 2) +
  plot_annotation(
    title      = "In-lake 0-10 m",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 9 * LO_FIG_SCALE, face = "bold",
                                             family = BASE_FAMILY),
                  plot.tag   = element_text(size = LO_FS_TAG, face = "bold",
                                            family = BASE_FAMILY)))

unlink(file.path(OUT_DIR, "supplemental_S2_nutrient_timeseries.png"))
ragg::agg_png(file.path(OUT_DIR, "supplemental_S2a_nutrient_timeseries_deposition.png"),
              width = lo_word_width(7.62), height = lo_word_height(7.62, 9.14), units = "cm", res = LO_DPI)
print(fig_s2_dep)
invisible(dev.off())
message("Saved: supplemental_S2a_nutrient_timeseries_deposition.png")

ragg::agg_png(file.path(OUT_DIR, "supplemental_S2b_nutrient_timeseries_inlake.png"),
              width = lo_word_width(7.62), height = lo_word_height(7.62, 9.14), units = "cm", res = LO_DPI)
print(fig_s2_lake)
invisible(dev.off())
message("Saved: supplemental_S2b_nutrient_timeseries_inlake.png")

# Final manuscript-facing S2: atmospheric and in-lake panels are paired
# horizontally by nutrient. A five-row layout preserves >8 pt text within the
# L&O 5 x 6 in limit; tags run consecutively across the complete figure.
dep_final_panels <- imap(nut_ids, function(id, j) {
  p <- make_nut_panel_s2(
    id, "Deposition", show_left_title = TRUE, show_right_title = TRUE,
    show_title = TRUE, show_x = (j == n_nut), compact_axis_titles = TRUE
  )
  if (j == 1) p <- p + labs(title = paste0(
    "Atmospheric deposition (", DEP_UNIT_S2, ")\n", nut_map[[id]]$disp
  ))
  p
})

lake_final_panels <- imap(nut_ids, function(id, j) {
  p <- make_nut_panel_s2(
    id, "In-lake", show_left_title = TRUE, show_right_title = TRUE,
    show_title = TRUE, show_x = (j == n_nut), compact_axis_titles = TRUE
  )
  if (j == 1) p <- p + labs(title = paste0(
    "In-lake 0-10 m (", LAKE_UNIT_S2, ")\n", nut_map[[id]]$disp
  ))
  p
})

s2_horizontal_panels <- unlist(
  map2(dep_final_panels, lake_final_panels, ~list(.x, .y)),
  recursive = FALSE
)

fig_s2_horizontal <- wrap_plots(s2_horizontal_panels, ncol = 2, byrow = TRUE) +
  plot_annotation(tag_levels = LO_TAG_LEVEL) &
  theme(
    plot.tag = element_text(size = LO_FS_TAG, face = "bold", family = BASE_FAMILY,
                            hjust = 1, vjust = 1),
    plot.tag.position = c(.98, .92)
  )

s2_horizontal_png <- file.path(OUT_DIR, "supplemental_S2_nutrient_timeseries_horizontal.png")
s2_horizontal_pdf <- file.path(OUT_DIR, "supplemental_S2_nutrient_timeseries_horizontal.pdf")
ragg::agg_png(s2_horizontal_png, width = LO_WIDTH_MAX, height = LO_HEIGHT_MAX,
              units = "cm", res = LO_DPI_LINE)
print(fig_s2_horizontal)
invisible(dev.off())
ggsave(s2_horizontal_pdf, fig_s2_horizontal, width = LO_WIDTH_MAX,
       height = LO_HEIGHT_MAX, units = "cm", device = cairo_pdf)
message("Saved L&O final: supplemental_S2_nutrient_timeseries_horizontal.png/.pdf")

if (identical(Sys.getenv("CALDOR_S2_ONLY"), "1")) {
  message("CALDOR_S2_ONLY=1; stopping after the final S2 export.")
  quit(save = "no", status = 0)
}

# =============================================================================
# -- S3: Cross-correlation heatmap: taxa vs. in-lake N:P by depth bin --------
np_monthly <- tryCatch(
  read_csv(NUTRI_LTP, show_col_types = FALSE) %>%
    mutate(date = as.Date(Date),
           depth = suppressWarnings(as.numeric(Depth)),
           depth_zone = case_when(depth >= 0 & depth <= 40 ~ "Surface community (0-40 m)",
                                  depth >= 60 & depth <= 105 ~ "Deep community (60-105 m)",
                                  TRUE ~ NA_character_),
           np_ratio = ((suppressWarnings(as.numeric(NO3)) +
                          suppressWarnings(as.numeric(NH4))) / 14.0067) /
             (suppressWarnings(as.numeric(TRP)) / 30.973762),
           month_date = floor_date(date, "month")) %>%
    filter(!is.na(depth_zone), is.finite(np_ratio), np_ratio > 0) %>%
    group_by(depth_zone, month_date) %>%
    summarise(NP_ratio = mean(np_ratio, na.rm = TRUE), .groups = "drop"),
  error = function(e) {
    cat("\033[31m  DIN:SRP ratio data not found - skipping S3\033[0m\n"); NULL }
)

if (!is.null(np_monthly)) {
  TOP_N_CCF <- 15
  MIN_LAG <- -6
  MAX_LAG <- 1
  LAG_MAX_ABS <- max(abs(c(MIN_LAG, MAX_LAG)))
  CCF_ZONES <- tibble(
    depth_zone = c("Surface community (0-40 m)", "Deep community (60-105 m)"),
    zmin = c(0, 60),
    zmax = c(40, 105)
  )

  compute_ccf_zone <- function(depth_zone, zmin, zmax) {
    zone_dat <- clean_ltp %>% filter(depth_m >= zmin, depth_m <= zmax)
    top_taxa_ccf <- zone_dat %>%
      group_by(taxon) %>%
      summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      slice_max(total, n = TOP_N_CCF, with_ties = FALSE) %>%
      pull(taxon)
    top_taxa_ccf <- unique(c(top_taxa_ccf, TAXON_LEPTO))

    phyto_monthly <- zone_dat %>%
      filter(taxon %in% top_taxa_ccf) %>%
      group_by(month_date, taxon) %>%
      summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop")

    np_zone <- np_monthly %>% filter(depth_zone == !!depth_zone)
    if (nrow(phyto_monthly) == 0 || nrow(np_zone) == 0) {
      return(tibble(depth_zone = depth_zone, taxon = character(),
                    lag = integer(), ccf_val = numeric(), n_months = integer()))
    }

    win_start <- max(min(phyto_monthly$month_date), min(np_zone$month_date))
    win_end <- min(max(phyto_monthly$month_date), max(np_zone$month_date))
    phyto_monthly <- phyto_monthly %>% filter(month_date >= win_start, month_date <= win_end)
    np_zone <- np_zone %>% filter(month_date >= win_start, month_date <= win_end)

    all_months <- seq(win_start, win_end, by = "month")
    full_grid <- expand.grid(month_date = all_months, taxon = top_taxa_ccf,
                             stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      mutate(month_date = as.Date(month_date, origin = "1970-01-01"))

    ccf_df <- full_grid %>%
      left_join(phyto_monthly, by = c("month_date", "taxon")) %>%
      left_join(np_zone, by = "month_date") %>%
      arrange(taxon, month_date)

    map_dfr(top_taxa_ccf, function(tx) {
      sub <- ccf_df %>% filter(taxon == tx) %>% arrange(month_date)
      cc_ok <- !is.na(sub$NP_ratio) & !is.na(sub$mean_abund)
      if (sum(cc_ok) < (2 * LAG_MAX_ABS + 3)) {
        return(tibble(depth_zone = depth_zone, taxon = tx,
                      lag = integer(0), ccf_val = numeric(0), n_months = integer(0)))
      }
      cc <- ccf(sub$NP_ratio[cc_ok], sub$mean_abund[cc_ok],
                lag.max = LAG_MAX_ABS, plot = FALSE)
      tibble(depth_zone = depth_zone, taxon = tx,
             lag = as.integer(cc$lag), ccf_val = as.numeric(cc$acf),
             n_months = sum(cc_ok)) %>%
        filter(lag >= MIN_LAG, lag <= MAX_LAG)
    })
  }

  ccf_results <- pmap_dfr(CCF_ZONES, compute_ccf_zone)

  if (nrow(ccf_results) > 0) {
    ccf_results <- ccf_results %>%
      group_by(depth_zone) %>%
      mutate(max_abs_taxon = ave(abs(ccf_val), taxon, FUN = max),
             taxa_order = dense_rank(max_abs_taxon)) %>%
      arrange(depth_zone, desc(max_abs_taxon), taxon, lag) %>%
      mutate(taxon_ord = factor(taxon, levels = rev(unique(taxon))),
             y_idx = as.integer(taxon_ord)) %>%
      ungroup()

    lim_val <- max(abs(ccf_results$ccf_val), na.rm = TRUE)
    sig_df <- ccf_results %>%
      mutate(sig_threshold = qnorm(0.975) / sqrt(pmax(n_months, 1))) %>%
      filter(abs(ccf_val) >= sig_threshold)

    ccf_diag <- ccf_results %>%
      group_by(depth_zone, taxon) %>%
      slice_max(abs(ccf_val), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(depth_zone, taxon, peak_lag_months = lag,
                peak_ccf = ccf_val, peak_abs_ccf = abs(ccf_val), n_months) %>%
      arrange(depth_zone, desc(peak_abs_ccf)) %>%
      group_by(depth_zone) %>%
      mutate(rank = row_number()) %>%
      ungroup()
    write_csv(ccf_diag, file.path(OUT_DIR, "supplemental_S3_ccf_diagnostics.csv"))

    lepto_diag <- ccf_diag %>% filter(taxon == TAXON_LEPTO)
    if (nrow(lepto_diag) > 0) {
      lepto_diag %>% pmap(function(depth_zone, taxon, peak_lag_months,
                                   peak_ccf, peak_abs_ccf, n_months, rank) {
        message("Leptolyngbya N:P CCF peak (", depth_zone, "): lag ",
                peak_lag_months, " months, r = ", round(peak_ccf, 3),
                "; rank ", rank, ".")
      })
    }

    y_breaks <- ccf_results %>%
      distinct(depth_zone, y_idx, taxon_ord) %>%
      arrange(depth_zone, y_idx)

    p_ccf <- ggplot(ccf_results, aes(x = lag, y = y_idx, fill = ccf_val)) +
      geom_tile(width = 0.98, height = 0.98, colour = "grey20", linewidth = 0.12 * LO_FIG_SCALE) +
      geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.45 * LO_FIG_SCALE,
                 linetype = "dashed") +
      geom_point(data = sig_df, aes(x = lag, y = y_idx), inherit.aes = FALSE,
                 shape = 21, fill = NA, colour = "grey20", size = 1.0 * LO_FIG_SCALE,
                 stroke = 0.35 * LO_FIG_SCALE) +
      facet_wrap(~ depth_zone, ncol = 1, scales = "free_y") +
      scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#e6550d",
                           midpoint = 0, limits = c(-lim_val, lim_val),
                           name = "CCF", labels = label_number(accuracy = 0.01)) +
      scale_x_continuous(breaks = seq(MIN_LAG, MAX_LAG, by = 1),
                         labels = as.character, expand = expansion(mult = 0.01)) +
      scale_y_continuous(breaks = y_breaks$y_idx,
                         labels = lo_abbreviate_taxon(as.character(y_breaks$taxon_ord)),
                         expand = expansion(add = 0.5)) +
      labs(x = "Lag (months, DIN:SRP molar ratio vs. phytoplankton)", y = NULL) +
      theme_bw(base_size = 9 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(face = "italic", size = 6.8 * LO_FIG_SCALE,
                                       family = BASE_FAMILY),
            axis.text.x = element_text(size = 7 * LO_FIG_SCALE, family = BASE_FAMILY),
            axis.title.x = element_text(size = 8 * LO_FIG_SCALE, family = BASE_FAMILY),
            strip.text = element_text(size = 7.4 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY),
            plot.title = element_text(size = 10 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY),
            legend.title = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY),
            legend.text = element_text(size = 5.8 * LO_FIG_SCALE, family = BASE_FAMILY),
            legend.key.size = unit(2 * LO_FIG_SCALE, "mm"))

    ragg::agg_png(file.path(OUT_DIR, "supplemental_S3_ccf_heatmap.png"),
                  width = lo_word_width(7.62), height = lo_word_height(7.62, 9.14), units = "cm", res = LO_DPI)
    print(p_ccf)
    invisible(dev.off())
    message("Saved: supplemental_S3_ccf_heatmap.png")
  } else {
    cat("\033[33m  No CCF results computed - check data overlap\033[0m\n")
  }
}

# =============================================================================# =============================================================================
# Okabe-Ito colorblind-friendly palette (as in exploratory analysis)
make_anom_depth_plot <- function(dat, zone_filter, zone_label) {
  z <- dat %>% filter(zone_filter(depth_m))

  top10 <- z %>%
    group_by(taxon) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    slice_max(total, n = 10, with_ties = FALSE) %>%
    arrange(desc(total)) %>%
    pull(taxon)

  monthly <- z %>%
    filter(taxon %in% top10) %>%
    group_by(taxon, month_date, month) %>%
    summarise(val = mean(abundance, na.rm = TRUE), .groups = "drop")

  clim <- monthly %>%
    filter(year(month_date) != FIRE_YEAR) %>%
    group_by(taxon, month) %>%
    summarise(clim_mean = mean(val, na.rm = TRUE),
              clim_sd   = sd(val,   na.rm = TRUE), .groups = "drop")

  anom <- monthly %>%
    left_join(clim, by = c("taxon", "month")) %>%
    mutate(anom_units = val - clim_mean,
           taxon_fac = factor(taxon, levels = top10))

  taxon_cols <- taxon_all_cols[top10]

  ggplot(anom, aes(x = month_date, colour = taxon_fac, fill = taxon_fac)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_hline(yintercept = 0, linewidth = 0.25 * LO_FIG_SCALE, colour = "grey40") +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = 0.45 * LO_FIG_SCALE) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(anom_units, 0)), alpha = 0.55,
                colour = NA, na.rm = TRUE) +
    geom_ribbon(aes(ymin = pmin(anom_units, 0), ymax = 0), alpha = 0.20,
                colour = NA, na.rm = TRUE) +
    geom_line(aes(y = anom_units), linewidth = 0.45 * LO_FIG_SCALE, na.rm = TRUE) +
    scale_fill_manual(values = taxon_cols, guide = "none") +
    scale_colour_manual(values = taxon_cols, guide = "none") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_scientific(digits = 2)) +
    facet_wrap(~ taxon_fac, ncol = 2, scales = "free_y",
               labeller = labeller(taxon_fac = lo_abbreviate_taxon)) +
    labs(title = zone_label,
         x = NULL, y = expression("Abundance anomaly (cells L"^-1*")")) +
    theme_bw(base_size = 7 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
    theme(panel.grid.minor = element_blank(),
          panel.spacing = unit(0.75 * LO_FIG_SCALE, "mm"),
          strip.text = element_text(size = 5.7 * LO_FIG_SCALE, face = "italic", family = BASE_FAMILY),
          axis.text.x = element_text(size = 5.5 * LO_FIG_SCALE, family = BASE_FAMILY),
          axis.text.y = element_text(size = 5.3 * LO_FIG_SCALE, family = BASE_FAMILY),
          axis.title.y = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY),
          plot.title = element_text(size = 8 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY))
}

anom_zones <- list(
  list(dat = clean_ltp,  zone_filter = function(d) d >= 0 & d <= 40,
       label = "LTP 0-40 m", key = "shallow"),
  list(dat = clean_ltp,  zone_filter = function(d) d >= 60 & d <= 105,
       label = "LTP 60-105 m", key = "deep"),
  list(dat = clean_mltp, zone_filter = function(d) d == "0-100",
       label = "MLTP 0-100 m", key = "mltp_shallow"),
  list(dat = clean_mltp, zone_filter = function(d) d == "150-450",
       label = "MLTP 150-450 m", key = "mltp_deep")
)

unlink(file.path(OUT_DIR, c("supplemental_S5_phyto_anomaly_shallow_part1.png",
                           "supplemental_S5_phyto_anomaly_shallow_part2.png",
                           "supplemental_S5_phyto_anomaly_deep_part1.png",
                           "supplemental_S5_phyto_anomaly_deep_part2.png")))
for (az in anom_zones) {
  p_anom <- make_anom_depth_plot(az$dat, az$zone_filter, az$label)
  fn <- file.path(OUT_DIR, paste0("supplemental_S5_phyto_anomaly_", az$key, ".png"))
  ragg::agg_png(fn, width = lo_word_width(7.62), height = lo_word_height(7.62, 9.14), units = "cm", res = LO_DPI)
  print(p_anom)
  invisible(dev.off())
  message("Saved: ", basename(fn))
}
make_combined_abundance_anomaly <- function() {
  z <- clean_ltp %>%
    filter(depth_m <= 105) %>%
    mutate(depth_zone = case_when(
      depth_m <= 40 ~ "0-40 m",
      depth_m >= 60 & depth_m <= 105 ~ "60-105 m",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(depth_zone))

  top_taxa_combined <- z %>%
    group_by(taxon) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    slice_max(total, n = 12, with_ties = FALSE) %>%
    pull(taxon)

  monthly <- z %>%
    filter(taxon %in% top_taxa_combined) %>%
    group_by(depth_zone, taxon, month_date, month) %>%
    summarise(val = mean(abundance, na.rm = TRUE), .groups = "drop")

  clim <- monthly %>%
    filter(year(month_date) != FIRE_YEAR) %>%
    group_by(depth_zone, taxon, month) %>%
    summarise(clim_mean = mean(val, na.rm = TRUE), .groups = "drop")

  anom <- monthly %>%
    left_join(clim, by = c("depth_zone", "taxon", "month")) %>%
    mutate(anom_units = val - clim_mean,
           taxon = factor(taxon, levels = rev(top_taxa_combined)),
           y_idx = as.integer(taxon)) %>%
    filter(!is.na(anom_units))

  lim <- quantile(abs(anom$anom_units), 0.98, na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- max(abs(anom$anom_units), na.rm = TRUE)

  ggplot(anom, aes(month_date, y_idx, fill = anom_units)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_tile(width = 31, height = 0.98, colour = "grey20", linewidth = 0.12 * LO_FIG_SCALE) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = 0.45 * LO_FIG_SCALE) +
    facet_wrap(~ depth_zone, ncol = 1, scales = "free_y") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = sort(unique(anom$y_idx)),
                       labels = lo_abbreviate_taxon(levels(anom$taxon)),
                       expand = expansion(add = 0.5)) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-lim, lim), oob = squish,
                         labels = label_scientific(),
                         name = expression("cells L"^-1)) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_bw(base_size = 8 * LO_FIG_SCALE, base_family = BASE_FAMILY) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.2 * LO_FIG_SCALE, face = "italic", family = BASE_FAMILY),
          axis.text.x = element_text(size = 6.5 * LO_FIG_SCALE, family = BASE_FAMILY),
          strip.text = element_text(size = 7.5 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY),
          plot.title = element_text(size = 9.5 * LO_FIG_SCALE, face = "bold", family = BASE_FAMILY),
          legend.position = "bottom",
          legend.title = element_text(size = 6.2 * LO_FIG_SCALE, family = BASE_FAMILY),
          legend.text = element_text(size = 5.8 * LO_FIG_SCALE, family = BASE_FAMILY),
          legend.key.width = unit(10 * LO_FIG_SCALE, "mm"),
          legend.key.height = unit(2 * LO_FIG_SCALE, "mm"))
}

p_s5_combined <- make_combined_abundance_anomaly()
ragg::agg_png(file.path(OUT_DIR, "supplemental_S5_phyto_abundance_anomaly_combined.png"),
              width = lo_word_width(7.62), height = lo_word_height(7.62, 9.14), units = "cm", res = LO_DPI)
print(p_s5_combined)
invisible(dev.off())
message("Saved: supplemental_S5_phyto_abundance_anomaly_combined.png")
message("\nAll supplemental outputs saved to: ", OUT_DIR)

# Additional supplemental suites are isolated in their own environments so
# their configuration variables cannot overwrite the core S1-S5 workflow.
supplemental_components <- c(
  "scripts/main_analysis/supplemental_components/ctd_profiles_stability_pkl.R",
  "scripts/main_analysis/supplemental_components/lake_number_diagnostic_2021.R",
  "scripts/main_analysis/supplemental_components/nutrient_profiles.R",
  "scripts/main_analysis/supplemental_components/nutrients_2020_historical_comparison.R",
  "scripts/main_analysis/supplemental_components/physical_drivers_anomaly.R",
  "scripts/main_analysis/supplemental_components/taxa_depth_profiles.R"
)
for (component in supplemental_components) {
  message("Running supplemental component: ", component)
  sys.source(component, envir = new.env(parent = globalenv()))
}
message("All supplemental figure suites completed.")
