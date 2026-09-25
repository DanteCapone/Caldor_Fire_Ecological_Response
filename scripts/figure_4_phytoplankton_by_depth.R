# =============================================================================
# figure_4_phytoplankton_by_depth.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Composite 5-panel figure showing depth-resolved dynamics:
#   (a-c) Aug-Oct 2021 Chl-a depth profiles vs. historical climatology
#         Different point shapes: 0-40 m (circle) vs. 60-105 m (triangle)
#         Dashed grey lines at 40 m and 60 m depth zone boundaries
#   (d)   Total biovolume time series, depths 0-40 m, with Chl-a right axis
#   (e)   Total biovolume time series, depths 60-105 m, with Chl-a right axis
#         Stars mark Aug, Sep, Oct 2021 (months shown in panels a-c)
#
# Caldor Fire ignition: 14 August 2021 (red dashed vline)
# Output: figures/figure_4_phytoplankton/figure_4_phyto_by_depth.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(ragg)
})
source("scripts/figure_aesthetics.R")

# ---- Paths -------------------------------------------------------------------
PHYTO_FILE <- file.path(
  "data", "lake_environmental_data", "phytoplankton",
  "PhytoData_LTP_MLTP_2019-2025_counts_biovolume.xlsx"
)
CHL_LTP <- file.path(
  "data", "lake_environmental_data", "chla", "Tahoe_LTP_Chl.csv"
)
CHL_MLTP <- file.path(
  "data", "lake_environmental_data", "chla", "Tahoe_MLTP_Chl.csv"
)
CHL_ALL <- file.path(
  "data", "lake_environmental_data", "chla", "terc_chla_all.csv"
)
OUT_DIR <- file.path("figures", "figure_4_phytoplankton")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY  <- "Times New Roman"
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")
TOP_N        <- 10

# ---- Load phytoplankton data -------------------------------------------------
# The Excel Date column contains a mix of serial numbers (2019-2024) and
# ISO-style text dates (2025). Reading it as text lets us parse both forms.
ltp_raw <- read_excel(PHYTO_FILE, sheet = "LTP_2019-2025",
                      col_types = c("text","text","text","text","text",
                                    "numeric","text","text","numeric","text","numeric"))

parse_excel_date <- function(x) {
  x <- trimws(as.character(x))
  # Try ISO text date first (e.g. "2025-01-16")
  d_iso <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  # Fall back to Excel serial number (Windows origin; e.g. "43493")
  n <- suppressWarnings(as.numeric(x))
  d_serial <- as.Date(n, origin = as.Date("1899-12-30"))
  coalesce(d_iso, d_serial)
}

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
                                     "text", "numeric")) %>%
  rename_with(~ "unit_biovol", matches("Unit_biovolume")) %>%
  rename(station_id = Station_ID,
         event_id   = Event_ID,
         depth_m    = all_of("Collection_depth_range_(m)"),
         phylum     = Phylum,
         taxon      = Taxon,
         abundance  = all_of("Abundance_(units/L)"))

clean_mltp <- mltp_raw %>%
  mutate(date        = parse_excel_date(Date),
         year        = year(date),
         month       = month(date),
         month_date  = floor_date(date, "month"),
         depth_m     = as.character(depth_m),
         abundance   = as.numeric(abundance),
         unit_biovol = as.numeric(unit_biovol),
         biovolume   = abundance * unit_biovol) %>%
  filter(!is.na(abundance), !is.na(depth_m), !is.na(date))

# ---- Top taxa & color palette ------------------------------------------------
all_taxa_levels <- sort(unique(c(clean_ltp$taxon, clean_mltp$taxon)))
qual_cols <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#9467BD",
               "#56B4E9", "#F0E442", "#8C564B", "#17BECF", "#CC79A7",
               "#999933", "#44AA99")
taxon_all_cols <- setNames(rep(qual_cols, length.out = length(all_taxa_levels)),
                           all_taxa_levels)
taxon_all_cols[intersect(names(taxon_all_cols), c("Leptolyngbya sp."))] <- "#009E73"
taxon_all_cols[intersect(names(taxon_all_cols), c("Cryptomonas sp."))] <- "#9467BD"
taxon_all_cols[intersect(names(taxon_all_cols), c("Cyclotella sp."))] <- "#0072B2"
taxon_all_cols[intersect(names(taxon_all_cols), c("Synedra acus var. radians"))] <- "#D55E00"

top_taxa <- bind_rows(
  clean_ltp %>% select(taxon, abundance),
  clean_mltp %>% select(taxon, abundance)
) %>%
  group_by(taxon) %>%
  summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total, n = TOP_N, with_ties = FALSE) %>%
  arrange(desc(total)) %>%
  pull(taxon)

taxon_cols <- c(taxon_all_cols[top_taxa], "Other" = "#BBBBBB")
mltp_qual_cols <- c("#0072B2", "#009E73", "#9467BD", "#C44E52", "#56B4E9",
                    "#CC79A7", "#8C564B", "#17BECF", "#999933", "#44AA99")
taxon_cols_mltp <- c(
  setNames(rep(mltp_qual_cols, length.out = length(top_taxa)), top_taxa),
  "Other" = "#BBBBBB"
)
taxon_cols_mltp[intersect(names(taxon_cols_mltp), "Leptolyngbya sp.")] <- "#009E73"
FIG4_WIDTH_CM  <- 11.43
FIG4_HEIGHT_CM <- 13.716
FIG4_TS_START  <- as.Date("2019-01-01")
FIG4_TS_END    <- as.Date("2025-12-31")

recode_taxa <- function(x) {
  factor(if_else(x %in% top_taxa, x, "Other"),
         levels = rev(c(top_taxa, "Other")))
}

# ---- Shared bar theme (L&O publication spec) --------------------------------
bar_theme <- theme_classic(base_size = 7, base_family = BASE_FAMILY) +
  theme(
    text              = element_text(family = BASE_FAMILY),
    axis.title        = element_text(size = 7, family = BASE_FAMILY),
    axis.text         = element_text(size = 6, family = BASE_FAMILY),
    plot.title        = element_text(face = "bold", size = 7, family = BASE_FAMILY),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 6, family = BASE_FAMILY),
    legend.text       = element_text(size = 6, family = BASE_FAMILY),
    legend.key.size   = unit(1.2 * LO_FIG_SCALE, "mm"),
    legend.key.spacing.x = unit(0.27 * LO_FIG_SCALE, "mm"),
    legend.key.spacing.y = unit(0.18 * LO_FIG_SCALE, "mm"),
    legend.margin     = margin(0, 0, 0, 0, "mm"),
    legend.box.spacing = unit(0.3 * LO_FIG_SCALE, "mm"),
    axis.title.y      = element_text(vjust = 0.5, margin = margin(r = 0.1 * LO_FIG_SCALE, unit = "mm")),
    panel.grid.major.y = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
    strip.text        = element_text(face = "bold", size = LO_FS_STRIP, family = BASE_FAMILY),
    strip.background  = element_rect(fill = "#F0F0F0", colour = "grey70")
  )

# =============================================================================
# -- Load chlorophyll data  ----------------------------------------------------
# =============================================================================
read_chl_recent <- function(path, station_label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(station = station_label, date = as.Date(Date),
           depth = as.numeric(Depth), chla = as.numeric(Chla)) %>%
    filter(!is.na(chla), !is.na(depth)) %>%
    select(station, date, depth, chla)
}

chl_ltp_recent <- read_chl_recent(CHL_LTP, "LTP")
chl_mltp_recent <- read_chl_recent(CHL_MLTP, "MLTP")

chl_all_raw <- tryCatch(
  read.csv(CHL_ALL, quote = "", stringsAsFactors = FALSE),
  error = function(e) read.csv(CHL_ALL, quote = "", stringsAsFactors = FALSE,
                                fill = TRUE, comment.char = "")
)

chl_all <- chl_all_raw %>%
  as_tibble() %>%
  filter(Sample_Type %in% c("FIELD","FLDDUP"),
         Station_ID  %in% c("Index", "Mid-lake")) %>%
  mutate(station = if_else(Station_ID == "Index", "LTP", "MLTP"),
         date    = as.Date(Date),
         depth   = suppressWarnings(as.numeric(Depth)),
         chla    = suppressWarnings(as.numeric(Chla))) %>%
  filter(!is.na(chla), !is.na(depth), chla > 0) %>%
  select(station, date, depth, chla)

chl_hist_only <- anti_join(chl_all, bind_rows(chl_ltp_recent, chl_mltp_recent),
                           by = c("station","date","depth"))

chl <- bind_rows(chl_hist_only, chl_ltp_recent, chl_mltp_recent) %>%
  mutate(year = year(date), month = month(date),
         month_lab  = month(date, label = TRUE, abbr = TRUE),
         month_date = floor_date(date, "month")) %>%
  filter(depth <= 150, chla < 50)

chl_ltp <- chl %>% filter(station == "LTP")
chl_mltp <- chl %>% filter(station == "MLTP")

cat("LTP Chl rows:", nrow(chl_ltp),
    "| date range:", format(range(chl_ltp$date)), "\n")
cat("MLTP Chl rows:", nrow(chl_mltp),
    "| date range:", format(range(chl_mltp$date)), "\n")

# =============================================================================
# -- Monthly Chl-a sums by depth zone (for secondary axis) --------------------
# =============================================================================
# Shallow: 0-40 m   |   Deep: 60-105 m
CHL_COL <- "#4CAF50"  # right-axis Chl point / axis colour

chl_shallow_ts <- chl_ltp %>%
  filter(year >= 2019, year <= 2025, depth <= 40) %>%
  group_by(month_date) %>%
  summarise(chl_mean = mean(chla, na.rm = TRUE),
            chl_sd   = sd(chla,  na.rm = TRUE),
            .groups  = "drop")

chl_deep_ts <- chl_ltp %>%
  filter(year >= 2019, year <= 2025, depth >= 60, depth <= 105) %>%
  group_by(month_date) %>%
  summarise(chl_mean = mean(chla, na.rm = TRUE),
            chl_sd   = sd(chla,  na.rm = TRUE),
            .groups  = "drop")

# =============================================================================
# -- Panels (a-c): Aug-Oct 2021 Chl depth profiles vs. historical -------------
# =============================================================================
FIRE_MONTHS <- c(8L, 9L, 10L, 11L, 12L, 1L)

chl_clim_p1 <- chl_ltp %>%
  filter(month %in% FIRE_MONTHS) %>%
  mutate(target_year = if_else(month == 1L, 2022L, 2021L)) %>%
  filter(year != target_year) %>%
  group_by(month, depth) %>%
  summarise(clim_mean = mean(chla, na.rm = TRUE),
            clim_sd   = sd(chla,  na.rm = TRUE),
            clim_n    = n(), .groups = "drop") %>%
  filter(clim_n >= 3, month %in% FIRE_MONTHS)

chl_obs_p1 <- chl_ltp %>%
  filter((year == 2021 & month %in% 8:12) |
           (year == 2022 & month == 1L)) %>%
  group_by(month, depth) %>%
  summarise(obs_mean = mean(chla, na.rm = TRUE),
            obs_sd   = sd(chla,   na.rm = TRUE),
            obs_n    = n(), .groups = "drop")

stats_p1 <- chl_obs_p1 %>%
  left_join(chl_clim_p1, by = c("month","depth")) %>%
  filter(!is.na(clim_mean), !is.na(clim_sd), clim_sd > 0) %>%
  mutate(
    z_score    = (obs_mean - clim_mean) / (clim_sd / sqrt(clim_n)),
    p_approx   = 2 * pnorm(-abs(z_score)),
    sig_label  = case_when(p_approx < 0.01 ~ "**", p_approx < 0.05 ~ "*", TRUE ~ ""),
    depth_zone = factor(if_else(depth <= 40, "0-40 m", "60-105 m"),
                        levels = c("0-40 m", "60-105 m"))
  )

# ---- Trapezoidal integrated Chl-a for Aug/Sep/Oct 2021 vs. climatology ------
trap_int <- function(depth, value) {
  ord <- order(depth); d <- depth[ord]; v <- value[ord]
  ok  <- !is.na(v);    d <- d[ok];     v <- v[ok]
  if (length(d) < 2) return(NA_real_)
  sum(diff(d) * (head(v, -1) + tail(v, -1)) / 2)
}

integ_obs  <- stats_p1 %>% group_by(month) %>%
  summarise(int_obs = trap_int(depth, obs_mean), .groups = "drop")

integ_clim_yearly <- chl_ltp %>%
  filter(month %in% FIRE_MONTHS) %>%
  filter(!((month == 1L & year == 2022L) |
           (month != 1L & year == 2021L))) %>%
  group_by(year, month) %>%
  summarise(int_year = trap_int(depth, chla), .groups = "drop") %>%
  filter(!is.na(int_year))

integ_clim <- integ_clim_yearly %>%
  group_by(month) %>%
  summarise(int_clim    = mean(int_year, na.rm = TRUE),
            int_clim_sd = sd(int_year,  na.rm = TRUE),
            int_clim_n  = n(),
            .groups = "drop")

integ_p1 <- left_join(integ_obs, integ_clim, by = "month") %>%
  mutate(
    target_year = if_else(month == 1L, 2022L, 2021L),
    pred_se = int_clim_sd * sqrt(1 + 1 / pmax(int_clim_n, 1)),
    t_pred  = if_else(!is.na(pred_se) & pred_se > 0,
                      (int_obs - int_clim) / pred_se, NA_real_),
    p_pred  = 2 * pt(-abs(t_pred), df = pmax(int_clim_n - 1, 1)),
    sd_anom = if_else(!is.na(int_clim_sd) & int_clim_sd > 0,
                      (int_obs - int_clim) / int_clim_sd, NA_real_),
    p_text  = case_when(
      is.na(p_pred) ~ "p = NA",
      p_pred < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", formatC(p_pred, digits = 2, format = "g"))
    ),
    narrative = sprintf(
      "2021 integrated Chl-a: %.1f mg m^-2\nhistorical same-month: %.1f +/- %.1f\nanomaly: %+.1f SD; %s",
      int_obs, int_clim, int_clim_sd, sd_anom, p_text
    )
  )

write_csv(integ_p1, file.path(OUT_DIR, "figure_4_integrated_chla_tests.csv"))

chl_theme <- theme_classic(base_size = 7, base_family = BASE_FAMILY) +
  theme(
    text             = element_text(family = BASE_FAMILY),
    axis.title       = element_text(size = 7, family = BASE_FAMILY),
    axis.text        = element_text(size = 6, family = BASE_FAMILY),
    plot.title       = element_text(face = "bold", size = 7, hjust = 0.62, family = BASE_FAMILY),
    panel.grid.major = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = LO_FS_LEGEND_TTL, family = BASE_FAMILY),
    legend.text      = element_text(size = 6, family = BASE_FAMILY),
    legend.key.size  = unit(1.2 * LO_FIG_SCALE, "mm"),
    legend.spacing.x = unit(0.3 * LO_FIG_SCALE, "mm"),
    legend.margin    = margin(0, 0, 0, 0, "mm"),
    axis.title.y    = element_text(margin = margin(r = -0.35, unit = "mm")),
    axis.title.x    = element_text(margin = margin(t = 0.5 * LO_FIG_SCALE, unit = "mm"))
  )

make_chl_panel <- function(mon_num, show_legend = FALSE,
                           show_y = TRUE, show_x_title = TRUE,
                           clim_df = chl_clim_p1,
                           stats_df = stats_p1) {
  dat_clim <- clim_df %>% filter(month == mon_num)
  dat_obs  <- stats_df %>% filter(month == mon_num)
  p <- ggplot() +
    geom_vline(xintercept = c(40, 60),
               linetype = "dashed", colour = "grey55", linewidth = LO_LW_THIN) +
    geom_ribbon(data = dat_clim,
                aes(x = depth, ymin = pmax(clim_mean - clim_sd, 0),
                    ymax = clim_mean + clim_sd),
                alpha = 0.20, fill = "#5AAE61") +
    geom_path(data = dat_clim, aes(x = depth, y = clim_mean),
              colour = "#5AAE61", linewidth = LO_LW_MID, linetype = "dashed") +
    geom_errorbar(data = dat_obs,
                  aes(x = depth,
                      ymin = pmax(obs_mean - obs_sd, 0),
                      ymax = obs_mean + obs_sd),
                  colour = "#1B7837", width = 1.0, linewidth = 0.3 * LO_FIG_SCALE,
                  alpha = 0.6) +
    geom_path(data = dat_obs, aes(x = depth, y = obs_mean),
              colour = "#1B7837", linewidth = LO_LW_THICK) +
    geom_point(data = dat_obs,
               aes(x = depth, y = obs_mean),
               colour = "#1B7837", fill = "#1B7837", shape = 21,
               size = 1.2 * LO_FIG_SCALE, stroke = 0.35 * LO_FIG_SCALE)
  p <- p +
    scale_x_reverse(breaks = c(0, 5, 10, 20, 30, 40, 50, 60, 75, 90, 105)) +
    scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_pretty(n = 2),
                       guide = guide_axis(check.overlap = TRUE),
                       expand = expansion(mult = c(0, 0.22))) +
    coord_flip(clip = "off") +
    labs(title = month.abb[mon_num],
         y     = expression("Chlorophyll-a (µg L"^-1*")"),
         x     = "Depth (m)") +
    chl_theme +
    theme(axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5,
                                     family = BASE_FAMILY),
          axis.text.y = element_text(size = 6, margin = margin(r = 0, unit = "pt")))

  p <- p + theme(legend.position = "none")
  if (!show_y)        p <- p + theme(axis.title.y = element_blank(),
                                     axis.text.y  = element_blank(),
                                     axis.ticks.y = element_blank())
  if (!show_x_title)  p <- p + theme(axis.title.x = element_blank())
  p
}

p1_aug <- make_chl_panel(8,  show_y = TRUE,  show_x_title = FALSE)
p1_sep <- make_chl_panel(9,  show_y = FALSE, show_x_title = FALSE)
p1_oct <- make_chl_panel(10, show_y = FALSE, show_x_title = TRUE)
p1_nov <- make_chl_panel(11, show_y = FALSE, show_x_title = FALSE)
p1_dec <- make_chl_panel(12, show_y = FALSE, show_x_title = FALSE)
p1_jan <- make_chl_panel(1,  show_legend = TRUE, show_y = FALSE, show_x_title = FALSE)

chl_mltp_clim_p1 <- chl_mltp %>%
  filter(month %in% FIRE_MONTHS) %>%
  mutate(target_year = if_else(month == 1L, 2022L, 2021L)) %>%
  filter(year != target_year) %>%
  group_by(month, depth) %>%
  summarise(clim_mean = mean(chla, na.rm = TRUE),
            clim_sd   = sd(chla,  na.rm = TRUE),
            clim_n    = n(), .groups = "drop") %>%
  filter(clim_n >= 3, month %in% FIRE_MONTHS)

chl_mltp_obs_p1 <- chl_mltp %>%
  filter((year == 2021 & month %in% 8:12) |
           (year == 2022 & month == 1L)) %>%
  group_by(month, depth) %>%
  summarise(obs_mean = mean(chla, na.rm = TRUE),
            obs_sd   = sd(chla,   na.rm = TRUE),
            obs_n    = n(), .groups = "drop")

stats_mltp_p1 <- chl_mltp_obs_p1 %>%
  left_join(chl_mltp_clim_p1, by = c("month","depth")) %>%
  filter(!is.na(clim_mean), !is.na(clim_sd), clim_sd > 0) %>%
  mutate(
    z_score    = (obs_mean - clim_mean) / (clim_sd / sqrt(clim_n)),
    p_approx   = 2 * pnorm(-abs(z_score)),
    sig_label  = case_when(p_approx < 0.01 ~ "**", p_approx < 0.05 ~ "*", TRUE ~ ""),
    depth_zone = factor(if_else(depth <= 100, "0-100 m", ">100 m"),
                        levels = c("0-100 m", ">100 m"))
  )

p_mltp_aug <- make_chl_panel(8,  show_y = TRUE,  show_x_title = FALSE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)
p_mltp_sep <- make_chl_panel(9,  show_y = FALSE, show_x_title = FALSE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)
p_mltp_oct <- make_chl_panel(10, show_y = FALSE, show_x_title = TRUE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)
p_mltp_nov <- make_chl_panel(11, show_y = FALSE, show_x_title = FALSE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)
p_mltp_dec <- make_chl_panel(12, show_y = FALSE, show_x_title = FALSE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)
p_mltp_jan <- make_chl_panel(1,  show_legend = TRUE, show_y = FALSE,
                             show_x_title = FALSE,
                             clim_df = chl_mltp_clim_p1, stats_df = stats_mltp_p1)

# =============================================================================
# -- Panels (d-e): Total biovolume timeseries by depth zone -------------------
#    Biovolume in µm³/L (scientific notation); Chl-a on right y-axis
# =============================================================================
make_metric_panel <- function(depth_fn, zone_title, metric = c("biovolume", "abundance"),
                              show_y_title = TRUE) {
  metric <- match.arg(metric)
  value_col <- if (metric == "biovolume") "biovolume" else "abundance"
  y_lab <- if (metric == "biovolume") {
    expression("Biovolume ("*mu*"m"^3*" L"^-1*")")
  } else {
    expression("Abundance (cells L"^-1*")")
  }

  bv <- clean_ltp %>%
    filter(depth_fn(depth_m),
           month_date >= FIG4_TS_START, month_date <= FIG4_TS_END) %>%
    mutate(taxon_grp = recode_taxa(taxon)) %>%
    group_by(month_date, taxon_grp) %>%
    summarise(total_value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  bv_tot <- bv %>% group_by(month_date) %>%
    summarise(tot = sum(total_value), .groups = "drop")
  y_upper <- max(bv_tot$tot, na.rm = TRUE) * 1.03
  if (!is.finite(y_upper) || y_upper <= 0) y_upper <- 1

  ggplot(bv, aes(x = month_date, y = total_value, fill = taxon_grp)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_col(position = "stack", width = 27, colour = NA) +
    geom_vline(xintercept = CALDOR_START,
               linetype = "dashed", colour = "firebrick", linewidth = LO_LW_MID) +
    scale_fill_manual(values = taxon_cols, labels = lo_abbreviate_taxon, name = "Taxon",
                      guide = guide_legend(nrow = 3, reverse = TRUE,
                                           byrow = TRUE,
                                           override.aes = list(linewidth = 0))) +
    scale_x_date(breaks = seq(FIG4_TS_START, as.Date("2025-01-01"), by = "1 year"),
                 date_labels = "%Y",
                 limits = c(FIG4_TS_START - 15, FIG4_TS_END + 15),
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(
      labels = label_scientific(),
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = zone_title, x = NULL,
         y = if (show_y_title) y_lab else NULL) +
    bar_theme +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0, margin = margin(b = 0.2 * LO_FIG_SCALE, unit = "mm"), family = BASE_FAMILY),
      legend.text = element_text(face = "italic", size = 6,
                                 family = BASE_FAMILY),
      legend.title = element_text(face = "bold", size = 6,
                                  family = BASE_FAMILY)
    )
}

p2 <- make_metric_panel(
  function(d) d <= 40,
  "0-40m",
  metric = "biovolume"
)

p3 <- make_metric_panel(
  function(d) d >= 60 & d <= 105,
  "60-105m",
  metric = "biovolume",
  show_y_title = TRUE
)

p2_abund <- make_metric_panel(
  function(d) d <= 40,
  "0-40m",
  metric = "abundance"
)

p3_abund <- make_metric_panel(
  function(d) d >= 60 & d <= 105,
  "60-105m",
  metric = "abundance",
  show_y_title = TRUE
)
# =============================================================================
# axis_titles = "collect" merges the shared left (Biovolume) and right (Chl-a)
# titles into a single, vertically centred label spanning both panels.
# Tight top/bottom plot margins compact the gap between panels g and h.
make_fig4_composite <- function(lower_top, lower_bottom) {
  bottom_panels <- (lower_top / lower_bottom) +
    plot_layout(guides = "collect", axes = "collect_y",
                axis_titles = "collect") &
    theme(legend.position = "bottom",
          plot.margin = margin(0.5 * LO_FIG_SCALE, 1 * LO_FIG_SCALE,
                               0.5 * LO_FIG_SCALE, 1 * LO_FIG_SCALE, "mm"))

  wrap_elements(full = wrap_plots(p1_aug, p1_sep, p1_oct, p1_nov, p1_dec, p1_jan,
             nrow = 1, widths = c(1.12, rep(1, 5)))) / bottom_panels +
    plot_layout(heights = c(4, 7)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.margin = margin(1.5, 1 * LO_FIG_SCALE, 1 * LO_FIG_SCALE, 1 * LO_FIG_SCALE, "mm"),
          plot.tag = element_text(face = "bold", size = 7,
                                  family = BASE_FAMILY, margin = margin(0, 0, 0, 0)),
          plot.tag.position = c(0.045, 0.985))
}

fig4 <- make_fig4_composite(p2, p3)
fig4_path <- file.path(OUT_DIR, "figure_4_phyto_by_depth.png")
ggsave(fig4_path, plot = fig4,
       width = FIG4_WIDTH_CM, height = FIG4_HEIGHT_CM, dpi = LO_DPI, units = "cm",
       device = ragg::agg_png)
cat("Saved:", fig4_path, "\n")

fig4_abund <- make_fig4_composite(p2_abund, p3_abund)
fig4_abund_path <- file.path(OUT_DIR, "figure_4_phyto_by_depth_v2_abundance.png")
ggsave(fig4_abund_path, plot = fig4_abund,
       width = FIG4_WIDTH_CM, height = FIG4_HEIGHT_CM, dpi = LO_DPI, units = "cm",
       device = ragg::agg_png)
cat("Saved:", fig4_abund_path, "\n")

# =============================================================================
# -- MLTP Figure 4 variation: Chl profiles + phytoplankton depth-bin panels ---
# =============================================================================
make_metric_panel_mltp <- function(depth_filter, zone_title,
                                   metric = c("biovolume", "abundance"),
                                   show_y_title = TRUE) {
  metric <- match.arg(metric)
  value_col <- if (metric == "biovolume") "biovolume" else "abundance"
  y_lab <- if (metric == "biovolume") {
    expression("Biovolume ("*mu*"m"^3*" L"^-1*")")
  } else {
    expression("Abundance (cells L"^-1*")")
  }

  dat <- clean_mltp %>%
    filter(depth_filter(depth_m),
           month_date >= FIG4_TS_START, month_date <= FIG4_TS_END) %>%
    mutate(taxon_grp = recode_taxa(taxon)) %>%
    group_by(month_date, taxon_grp) %>%
    summarise(total_value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  y_upper <- dat %>%
    group_by(month_date) %>%
    summarise(total = sum(total_value), .groups = "drop") %>%
    summarise(y_upper = max(total, na.rm = TRUE) * 1.03) %>%
    pull(y_upper)
  if (!is.finite(y_upper) || y_upper <= 0) y_upper <- 1

  ggplot(dat, aes(x = month_date, y = total_value, fill = taxon_grp)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_col(position = "stack", width = 27, colour = NA) +
    geom_vline(xintercept = CALDOR_START,
               linetype = "dashed", colour = "firebrick", linewidth = LO_LW_MID) +
    scale_fill_manual(values = taxon_cols_mltp, labels = lo_abbreviate_taxon, name = "Taxon",
                      guide = guide_legend(nrow = 3, reverse = TRUE,
                                           byrow = TRUE,
                                           override.aes = list(linewidth = 0))) +
    scale_x_date(breaks = seq(FIG4_TS_START, as.Date("2025-01-01"), by = "1 year"),
                 date_labels = "%Y",
                 limits = c(FIG4_TS_START - 15, FIG4_TS_END + 15),
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(
      labels = label_scientific(),
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = zone_title, x = NULL,
         y = if (show_y_title) y_lab else NULL) +
    bar_theme +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0,
                                    margin = margin(b = 0.2 * LO_FIG_SCALE, unit = "mm"),
                                    family = BASE_FAMILY),
          legend.text = element_text(face = "italic", size = 6,
                                     family = BASE_FAMILY),
          legend.title = element_text(face = "bold", size = 6,
                                      family = BASE_FAMILY),
          axis.text.y = element_text(size = 6, margin = margin(r = 0, unit = "pt")))
}

make_fig4_mltp <- function(metric = c("biovolume", "abundance")) {
  metric <- match.arg(metric)
  p_top <- make_metric_panel_mltp(function(d) d == "0-100", "0-100m",
                                  metric = metric)
  p_bottom <- make_metric_panel_mltp(function(d) d == "150-450", "150-450m",
                                     metric = metric)
  chl_row <- wrap_elements(full = wrap_plots(p_mltp_aug, p_mltp_sep, p_mltp_oct,
                        p_mltp_nov, p_mltp_dec, p_mltp_jan,
                        nrow = 1, widths = c(1.12, rep(1, 5))))
  bottom_panels <- (p_top / p_bottom) +
    plot_layout(guides = "collect", axes = "collect_y",
                axis_titles = "collect")
  chl_row / bottom_panels +
    plot_layout(heights = c(4, 7)) +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom",
          plot.margin = margin(1.5, 1 * LO_FIG_SCALE,
                               0.5 * LO_FIG_SCALE, 1 * LO_FIG_SCALE, "mm"),
          plot.tag = element_text(face = "bold", size = 7,
                                  family = BASE_FAMILY, margin = margin(0, 0, 0, 0)),
          plot.tag.position = c(0.045, 0.985))
}

fig4_mltp <- make_fig4_mltp("biovolume")
fig4_mltp_path <- file.path(OUT_DIR, "figure_4_phyto_by_depth_mltp.png")
ggsave(fig4_mltp_path, plot = fig4_mltp,
       width = FIG4_WIDTH_CM, height = FIG4_HEIGHT_CM, dpi = LO_DPI, units = "cm",
       device = ragg::agg_png)
cat("Saved:", fig4_mltp_path, "\n")

fig4_mltp_abund <- make_fig4_mltp("abundance")
fig4_mltp_abund_path <- file.path(OUT_DIR, "figure_4_phyto_by_depth_mltp_v2_abundance.png")
ggsave(fig4_mltp_abund_path, plot = fig4_mltp_abund,
       width = FIG4_WIDTH_CM, height = FIG4_HEIGHT_CM, dpi = LO_DPI, units = "cm",
       device = ragg::agg_png)
cat("Saved:", fig4_mltp_abund_path, "\n")

# =============================================================================
make_depth_time_anom <- function(obs_df, value_col) {
  monthly_depth <- obs_df %>%
    mutate(month_date = floor_date(date, "month"),
           month = month(date), year = year(date)) %>%
    group_by(month_date, year, month, depth) %>%
    summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  full_grid <- expand_grid(
    month_date = seq(as.Date("2019-01-01"), as.Date("2025-12-01"), by = "month"),
    depth = sort(unique(monthly_depth$depth))
  ) %>%
    mutate(year = year(month_date), month = month(month_date))

  monthly_depth <- full_grid %>%
    left_join(monthly_depth, by = c("month_date", "year", "month", "depth"))

  clim <- monthly_depth %>%
    filter(year != 2021) %>%
    group_by(month, depth) %>%
    summarise(clim_mean = mean(value, na.rm = TRUE),
              clim_sd   = sd(value, na.rm = TRUE),
              .groups = "drop")

  monthly_depth %>%
    left_join(clim, by = c("month", "depth")) %>%
    mutate(anom_sd = if_else(!is.na(clim_sd) & clim_sd > 0,
                             (value - clim_mean) / clim_sd, NA_real_),
           anom_sd = pmax(pmin(anom_sd, 4), -4))
}

chl_heat <- chl_ltp %>%
  filter(year >= 2019, year <= 2025, depth <= 105) %>%
  select(date, depth, chla) %>%
  make_depth_time_anom("chla")

lepto_heat <- clean_ltp %>%
  filter(taxon == "Leptolyngbya sp.", year >= 2019, year <= 2025,
         depth_m <= 105, !depth_m %in% c(0, 10)) %>%
  transmute(date, depth = depth_m, abundance = log10(abundance + 1)) %>%
  make_depth_time_anom("abundance")

make_heat_panel <- function(df, title_txt, fill_title) {
  ggplot(df, aes(x = month_date, y = depth, fill = anom_sd)) +
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.10) +
    geom_tile(width = 31, height = 5, colour = "grey20", linewidth = 0.12 * LO_FIG_SCALE) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = LO_LW_MID) +
    scale_y_reverse(breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish,
                         name = fill_title, na.value = "grey85") +
    labs(x = NULL, y = "Depth (m)", title = title_txt) +
    lo_theme(base_size = 4.2) +
    theme(legend.position = "right",
          plot.title = element_text(size = LO_FS_TITLE, face = "bold",
                                    family = BASE_FAMILY))
}

p_heat_chl <- make_heat_panel(chl_heat, "Chlorophyll-a", "SD")
p_heat_lep <- make_heat_panel(lepto_heat, "Leptolyngbya sp.", "SD")
fig4_heat <- p_heat_chl / p_heat_lep +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 7,
                                  family = BASE_FAMILY, margin = margin(0, 0, 0, 0)),
          plot.tag.position = c(0.045, 0.985))

heat_path <- file.path(OUT_DIR, "figure_4_supporting_chla_lepto_heatmaps.png")
ggsave(heat_path, plot = fig4_heat,
       width = FIG4_WIDTH_CM, height = FIG4_HEIGHT_CM, dpi = LO_DPI, units = "cm",
       device = ragg::agg_png)
cat("Saved:", heat_path, "\n")
cat("Done.\n")
