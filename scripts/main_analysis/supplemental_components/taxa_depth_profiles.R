# =============================================================================
# supplemental_taxa_depth_profiles.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Extended taxon depth-profile supplemental figures (LTP station).
#   For each taxon-year combination below, two figure types are produced:
#     1) Monthly depth-profile panel grids (abundance and biovolume), with a
#        fixed x-axis scale across all months in that taxon-year-variable
#        combination for direct visual comparison.
#     2) A depth-time "temporal evolution" figure: sampled depths are
#        linearly interpolated within each month and rendered as a
#        continuous linear color gradient (depth x month heatmap).
#
#   Combinations:
#     - Leptolyngbya sp., Apr 2021 - Apr 2022 (abundance, biovolume)
#     - Cryptomonas sp.,  Apr 2021 - Apr 2022 (abundance, biovolume)
#     - Leptolyngbya sp., Apr 2011 - Apr 2012 (abundance, biovolume)
#
# Data source: PhytoData_LTP_2005-2025_counts_biovolume_size_v1.xlsx
#              (station LTP only; this workbook spans 2005-2025 and is the
#              only source covering the 2011-2012 window).
#
# Output: figures/supplemental/taxa_profiles/
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
                        "PhytoData_LTP_2005-2025_counts_biovolume_size_v1.xlsx")
PHYTO_SHEET <- "LTP_2005-2025"
OUT_DIR <- file.path("figures", "supplemental", "taxa_profiles")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY <- "Times New Roman"

# Consistent taxon color palette (matches scripts/supplemental_figures.R)
COL_LEPTO <- "#009E73"    # Leptolyngbya - teal/green
COL_CRYPT <- "#9467BD"    # Cryptomonas  - purple

DEPTH_BREAKS <- c(0, 5, 10, 20, 30, 40, 50, 60, 75, 90, 105)

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

raw <- read_excel(
  PHYTO_FILE, sheet = PHYTO_SHEET,
  col_types = c("text", "text", "text", "text", "text", "numeric", "text",
                "text", "numeric", "text", "numeric", "text", "text")
)

clean_ltp <- raw %>%
  transmute(
    station_id  = Station_ID,
    date        = parse_excel_date(Date),
    depth_m     = as.numeric(`Depth_(m)`),
    taxon       = str_squish(Taxon),
    abundance   = as.numeric(`Abundance_(units/L)`),
    unit_biovol = as.numeric(`Unit_biovolume_(µm3)`)
  ) %>%
  mutate(
    year      = year(date),
    month     = month(date),
    biovolume = abundance * unit_biovol / 1e9   # mm3 L^-1; 1 mm3 = 1e9 um3
  ) %>%
  filter(station_id == "LTP", !is.na(abundance), !is.na(depth_m),
         !is.na(date), depth_m <= 105)

cat("LTP phytoplankton rows loaded:", nrow(clean_ltp), "\n")

VAR_INFO <- list(
  abundance = list(
    col   = "abundance",
    label = expression("Abundance (cells L"^-1*")"),
    short = "Abundance"
  ),
  biovolume = list(
    col   = "biovolume",
    label = expression("Biovolume (mm"^3*" L"^-1*")"),
    short = "Biovolume"
  )
)

# April(year1) - April(year2), 13 months, chronological order.
make_month_seq <- function(year1) {
  tibble(year  = c(rep(year1, 9L), rep(year1 + 1L, 4L)),
         month = c(4:12, 1:4))
}

# =============================================================================
# -- Monthly depth-profile panel grid (fixed x-scale within combo) -----------
# =============================================================================
make_month_panel <- function(taxon_name, taxon_col, target_year, mon_num,
                             value_col, value_lab, x_limits, y_limits,
                             show_y = TRUE, show_x_title = TRUE) {
  panel_title <- paste(month.abb[mon_num], target_year)

  dat_obs <- clean_ltp %>%
    filter(taxon == taxon_name, year == target_year, month == mon_num) %>%
    group_by(depth_m) %>%
    summarise(val = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  if (nrow(dat_obs) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No sample", colour = "grey45") +
        labs(title = panel_title) +
        theme_void(base_size = 7, base_family = BASE_FAMILY) +
        theme(plot.title = element_text(face = "bold", size = 7))
    )
  }

  p <- ggplot(dat_obs, aes(x = val, y = depth_m)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25 * LO_FIG_SCALE) +
    geom_path(colour = taxon_col, linewidth = LO_LW_THICK, na.rm = TRUE) +
    geom_point(colour = taxon_col, fill = taxon_col, shape = 21,
               size = 1.25 * LO_FIG_SCALE, stroke = 0.35 * LO_FIG_SCALE, na.rm = TRUE) +
    scale_y_reverse(breaks = DEPTH_BREAKS, limits = y_limits) +
    scale_x_continuous(labels = label_scientific(digits = 2),
                       limits = x_limits,
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = panel_title, y = "Depth (m)", x = value_lab) +
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

  if (!show_y)       p <- p + theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank())
  if (!show_x_title) p <- p + theme(axis.title.x = element_blank())
  p
}

make_panel_grid_figure <- function(taxon_name, taxon_col, taxon_label,
                                   year1, variable, fname) {
  months_tbl <- make_month_seq(year1)
  vinfo <- VAR_INFO[[variable]]

  dat_all <- clean_ltp %>%
    filter(taxon == taxon_name) %>%
    filter((year == year1 & month >= 4) | (year == year1 + 1L & month <= 4)) %>%
    group_by(year, month, depth_m) %>%
    summarise(val = mean(.data[[vinfo$col]], na.rm = TRUE), .groups = "drop")

  x_max <- suppressWarnings(max(dat_all$val, na.rm = TRUE))
  if (!is.finite(x_max) || x_max <= 0) x_max <- 1
  x_limits <- c(0, x_max)
  y_limits <- c(0, max(dat_all$depth_m, na.rm = TRUE))

  n_panels <- nrow(months_tbl)
  ncol_grid <- 4L

  panels <- pmap(
    list(months_tbl$year, months_tbl$month, seq_len(n_panels)),
    function(yr, m, i) {
      make_month_panel(
        taxon_name, taxon_col, yr, m, vinfo$col, vinfo$label, x_limits, y_limits,
        show_y       = (i %% ncol_grid == 1),
        show_x_title = (i == n_panels)
      )
    }
  )

  fig <- wrap_plots(panels, ncol = ncol_grid, byrow = TRUE) +
    plot_annotation(
      title      = paste0(taxon_label, ": ", vinfo$short),
      tag_levels = "a",
      theme = theme(plot.title = element_text(size = 8, face = "bold",
                                               family = BASE_FAMILY),
                    plot.tag   = element_text(size = 7, face = "bold",
                                              family = BASE_FAMILY),
                    plot.margin = margin(2.5, 2.5, 2.0, 2.5, "mm"))
    )

  out_path <- file.path(OUT_DIR, fname)
  unlink(out_path)
  ragg::agg_png(out_path, width = 11.4, height = 14.0, units = "cm", res = LO_DPI)
  print(fig)
  invisible(dev.off())
  message("Saved: ", fname)
}

# =============================================================================
# -- Temporal evolution (depth x month) - linear color gradient --------------
# =============================================================================
make_heatmap_figure <- function(taxon_name, taxon_col, taxon_label,
                                year1, variable, fname) {
  months_tbl <- make_month_seq(year1) %>%
    mutate(month_label = factor(paste(month.abb[month], year),
                                levels = paste(month.abb[month], year)))
  vinfo <- VAR_INFO[[variable]]

  monthly_obs <- clean_ltp %>%
    filter(taxon == taxon_name) %>%
    filter((year == year1 & month >= 4) | (year == year1 + 1L & month <= 4)) %>%
    group_by(year, month, depth_m) %>%
    summarise(val = mean(.data[[vinfo$col]], na.rm = TRUE), .groups = "drop") %>%
    left_join(months_tbl, by = c("year", "month"))

  depth_grid <- seq(0, 105, by = 1)

  interp_month <- function(yr, m, label) {
    d <- monthly_obs %>% filter(year == yr, month == m) %>% arrange(depth_m)
    if (nrow(d) < 2) return(tibble(month_label = label, depth_m = depth_grid, val = NA_real_))
    out <- approx(x = d$depth_m, y = d$val, xout = depth_grid, rule = 1)
    tibble(month_label = label, depth_m = depth_grid, val = out$y)
  }

  grid_df <- pmap_dfr(
    list(months_tbl$year, months_tbl$month, months_tbl$month_label),
    interp_month
  ) %>%
    mutate(month_label = factor(month_label, levels = levels(months_tbl$month_label)))

  obs_points <- monthly_obs %>% select(month_label, depth_m)

  p <- ggplot(grid_df, aes(x = month_label, y = depth_m, fill = val)) +
    geom_raster(na.rm = FALSE) +
    geom_point(data = obs_points, aes(x = month_label, y = depth_m),
               inherit.aes = FALSE, shape = 3, size = 0.6, stroke = 0.3,
               colour = "grey20") +
    scale_y_reverse(breaks = DEPTH_BREAKS, expand = expansion(mult = c(0.02, 0.02))) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_gradient(low = "white", high = taxon_col, na.value = "grey92",
                       name = vinfo$label) +
    labs(title = paste0(taxon_label, ": ", vinfo$short, " depth-time distribution"),
         x = NULL, y = "Depth (m)") +
    theme_classic(base_size = 8, base_family = BASE_FAMILY) +
    theme(axis.title = element_text(size = 7.5, family = BASE_FAMILY),
          axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1,
                                     family = BASE_FAMILY),
          axis.text.y = element_text(size = 7, family = BASE_FAMILY),
          plot.title = element_text(face = "bold", size = 8, family = BASE_FAMILY),
          legend.title = element_text(size = 7, family = BASE_FAMILY),
          legend.text = element_text(size = 6.5, family = BASE_FAMILY),
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(3, "mm"),
          plot.margin = margin(2, 2, 2, 2, "mm"))

  out_path <- file.path(OUT_DIR, fname)
  unlink(out_path)
  ragg::agg_png(out_path, width = 11.4, height = 8.5, units = "cm", res = LO_DPI)
  print(p)
  invisible(dev.off())
  message("Saved: ", fname)
}

# =============================================================================
# -- Build all taxon-year-variable combinations ------------------------------
# =============================================================================
TAXON_LEPTO <- "Leptolyngbya sp."
TAXON_CRYPT <- "Cryptomonas sp."

combos <- list(
  list(taxon = TAXON_LEPTO, col = COL_LEPTO, label = "Leptolyngbya sp.",
       year1 = 2021, tag = "leptolyngbya_2021_2022"),
  list(taxon = TAXON_CRYPT, col = COL_CRYPT, label = "Cryptomonas sp.",
       year1 = 2021, tag = "cryptomonas_2021_2022"),
  list(taxon = TAXON_LEPTO, col = COL_LEPTO, label = "Leptolyngbya sp.",
       year1 = 2011, tag = "leptolyngbya_2011_2012")
)

for (cmb in combos) {
  for (variable in c("abundance", "biovolume")) {
    make_panel_grid_figure(
      cmb$taxon, cmb$col, cmb$label, cmb$year1, variable,
      fname = paste0("taxa_profiles_", cmb$tag, "_", variable, "_panels.png")
    )
    make_heatmap_figure(
      cmb$taxon, cmb$col, cmb$label, cmb$year1, variable,
      fname = paste0("taxa_profiles_", cmb$tag, "_", variable, "_heatmap.png")
    )
  }
}

message("All taxa depth-profile figures saved to: ", OUT_DIR)
