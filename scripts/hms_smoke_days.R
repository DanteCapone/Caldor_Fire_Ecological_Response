# =============================================================================
# hms_smoke_days.R
#
# Purpose : Quantify Lake Tahoe "smoke days" directly from the NOAA Hazard
#           Mapping System (HMS) Smoke Product, replicating the lake
#           smoke-day method of Farruggia et al. (2024) for a single lake.
#
# Method  : Following Farruggia et al. (2024) and Paul et al. (2023):
#           A lake smoke-day is any day on which any portion of the lake
#           boundary intersects an area mapped as smoke by NOAA HMS. HMS
#           classifies daily smoke density as Low (~5 μg m⁻³), Medium (~16 μg m⁻³),
#           or High (~27 μg m⁻³) from satellite-derived aerosol optical depth (AOD).
#           Here each smoke-day is assigned the MAXIMUM density category that
#           intersected the lake on that day, so categories are mutually
#           exclusive and sum to the total number of smoke-days. Smoke-days
#           are then summed annually (and monthly).
#
# Data    : NOAA HMS Smoke Polygons (daily shapefiles), 2019-2025
#           https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/
#           Ruminski et al. (2006); product page:
#           https://www.ospo.noaa.gov/products/land/hms.html
#
# Lake    : Lake Tahoe shoreline polygon (same approximate outline used by the
#           project's satellite masks, scripts/satellite/tahoe_masks.py).
#
# Inputs  : downloaded on demand to data/raw/hms_smoke/
# Outputs : data/processed/tahoe_hms_smoke_daily.csv
#           data/processed/tahoe_hms_smoke_annual.csv
#           data/processed/tahoe_hms_smoke_monthly.csv
#           figures/hms_smoke/tahoe_hms_smoke_days_annual.png
#           figures/hms_smoke/tahoe_hms_smoke_days_monthly_heatmap.png
#
# Author  : Dante A. Capone
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(sf)
  library(here)
  library(glue)
  library(patchwork)
})

sf::sf_use_s2(FALSE)   # planar ops are fine at this scale and avoid s2 errors

# =============================================================================
# SECTION 0: Configuration
# =============================================================================

YEARS      <- 2006:2025
HMS_BASE   <- "https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile"

RAW_DIR    <- here("data", "raw", "hms_smoke")
PROC_DIR   <- here("data", "processed")
FIG_DIR    <- here("figures", "hms_smoke")
dir.create(RAW_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

# Caldor Fire dates (for annotation)
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")

# HMS smoke density categories with estimated surface PM2.5 concentrations
# (NOAA HMS AOD codes: 5 = Low ~5 μg m⁻³, 16 = Medium ~16 μg m⁻³, 27 = High ~27 μg m⁻³).
# Older files store density as text ("Light"/"Medium"/"Heavy"); both are handled.
DENSITY_LEVELS <- c("Low", "Medium", "High")
DENSITY_COLS   <- c(
  Low    = "#FED98E",   # pale yellow  (~5 μg m⁻³)
  Medium = "#FB8C00",   # orange       (~16 μg m⁻³)
  High   = "#B71C1C"    # deep red     (~27 μg m⁻³)
)

# Key fire/event years to annotate above bars
EVENT_YEARS <- c(
  "2020" = "2020 Wildfires",
  "2021" = "Caldor Fire"
)

# =============================================================================
# SECTION 1: Lake Tahoe boundary
# (matches scripts/satellite/tahoe_masks.py TAHOE_VERTICES; EPSG:4326)
# =============================================================================

tahoe_vertices <- matrix(
  c(-120.1320, 39.2520,
    -120.1640, 39.2300,
    -120.1900, 39.1900,
    -120.2050, 39.1500,
    -120.2050, 39.1100,
    -120.1950, 39.0700,
    -120.1750, 39.0400,
    -120.1500, 39.0100,
    -120.1100, 38.9700,
    -120.0700, 38.9500,
    -120.0400, 38.9350,
    -120.0050, 38.9300,
    -119.9700, 38.9350,
    -119.9450, 38.9500,
    -119.9300, 38.9700,
    -119.9200, 39.0000,
    -119.9150, 39.0400,
    -119.9150, 39.0800,
    -119.9200, 39.1200,
    -119.9300, 39.1600,
    -119.9450, 39.1900,
    -119.9700, 39.2100,
    -120.0000, 39.2250,
    -120.0350, 39.2400,
    -120.0700, 39.2500,
    -120.1000, 39.2540,
    -120.1320, 39.2520),  # close ring
  ncol = 2, byrow = TRUE
)

tahoe <- sf::st_sfc(sf::st_polygon(list(tahoe_vertices)), crs = 4326) |>
  sf::st_make_valid()

# =============================================================================
# SECTION 2: Download daily HMS smoke shapefiles
# =============================================================================

#' Download one day's HMS smoke shapefile (zip) into RAW_DIR/<year>/.
#' Returns the path to the extracted .shp, or NA if unavailable / no smoke file.
download_hms_day <- function(date) {
  ymd_str <- format(date, "%Y%m%d")
  yr      <- format(date, "%Y")
  mo      <- format(date, "%m")

  year_dir <- file.path(RAW_DIR, yr)
  dir.create(year_dir, showWarnings = FALSE, recursive = TRUE)

  zip_path <- file.path(year_dir, glue("hms_smoke{ymd_str}.zip"))
  shp_path <- file.path(year_dir, glue("hms_smoke{ymd_str}.shp"))

  # Already extracted
  if (file.exists(shp_path)) return(shp_path)

  url <- glue("{HMS_BASE}/{yr}/{mo}/hms_smoke{ymd_str}.zip")

  ok <- tryCatch({
    utils::download.file(url, destfile = zip_path, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!ok || !file.exists(zip_path) || file.info(zip_path)$size < 200) {
    if (file.exists(zip_path)) file.remove(zip_path)
    return(NA_character_)
  }

  extracted <- tryCatch(
    utils::unzip(zip_path, exdir = year_dir),
    error = function(e) character(0),
    warning = function(w) character(0)
  )

  if (!file.exists(shp_path)) return(NA_character_)
  shp_path
}

# =============================================================================
# SECTION 3: Extract lake smoke-day records from one daily shapefile
# =============================================================================

#' Read a daily HMS shapefile, intersect with the lake, and return the set of
#' density categories that touched the lake that day (character vector, may be
#' empty). Returns NULL when the file is unreadable.
extract_lake_densities <- function(shp_path) {
  if (is.na(shp_path) || !file.exists(shp_path)) return(NULL)

  sm <- tryCatch(
    suppressWarnings(sf::st_read(shp_path, quiet = TRUE)),
    error = function(e) NULL
  )
  if (is.null(sm) || nrow(sm) == 0) return(character(0))

  # HMS files are WGS84; force CRS if the .prj was missing.
  if (is.na(sf::st_crs(sm))) sf::st_crs(sm) <- 4326
  sm <- sf::st_transform(sm, 4326)
  sm <- tryCatch(sf::st_make_valid(sm), error = function(e) sm)
  sm <- sm[!sf::st_is_empty(sm), ]
  if (nrow(sm) == 0) return(character(0))

  # Keep only polygons that intersect the lake
  hits <- tryCatch(
    suppressMessages(sf::st_intersects(sm, tahoe, sparse = FALSE)[, 1]),
    error = function(e) rep(FALSE, nrow(sm))
  )
  sm   <- sm[hits, ]
  if (nrow(sm) == 0) return(character(0))

  # Resolve the density attribute (column name varies: Density / DENSITY)
  dcol <- names(sm)[tolower(names(sm)) == "density"]
  if (length(dcol) == 0) {
    # No density attribute -> treat as unspecified smoke (count as Light)
    return("Low")
  }
  raw <- sm[[dcol[1]]]

  classify_density(raw)
}

#' Map raw HMS density values (numeric AOD code or text) to category labels.
classify_density <- function(raw) {
  txt <- trimws(as.character(raw))
  num <- suppressWarnings(as.numeric(txt))

  out <- character(length(txt))
  # Numeric AOD codes: 5 = Low (~5 μg m⁻³), 16 = Medium (~16 μg m⁻³), 27 = High (~27 μg m⁻³)
  is_num <- !is.na(num)
  out[is_num] <- cut(num[is_num],
                     breaks = c(-Inf, 10, 21, Inf),
                     labels = DENSITY_LEVELS) |> as.character()
  # Text labels (older files; "Light"→"Low", "Heavy"→"High")
  is_txt <- !is_num
  lt <- tolower(txt[is_txt])
  out[is_txt] <- dplyr::case_when(
    str_detect(lt, "heav|high") ~ "High",
    str_detect(lt, "med")       ~ "Medium",
    str_detect(lt, "light|low") ~ "Low",
    TRUE                        ~ "Low"
  )
  unique(out[out %in% DENSITY_LEVELS])
}

# =============================================================================
# SECTION 4: Build the daily smoke-day table
# =============================================================================

all_dates <- seq(as.Date(glue("{min(YEARS)}-01-01")),
                 as.Date(glue("{max(YEARS)}-12-31")), by = "day")

cat(glue("Processing {length(all_dates)} days ({min(YEARS)}-{max(YEARS)}) ...\n\n"))

records <- vector("list", length(all_dates))
pb_step <- max(1, floor(length(all_dates) / 50))

for (i in seq_along(all_dates)) {
  d   <- all_dates[i]
  shp <- download_hms_day(d)
  dens <- extract_lake_densities(shp)

  file_ok <- !is.na(shp)
  max_cat <- NA_character_
  if (length(dens) > 0) {
    max_cat <- DENSITY_LEVELS[max(match(dens, DENSITY_LEVELS))]
  }

  records[[i]] <- tibble(
    date        = d,
    file_avail  = file_ok,
    smoke_day   = file_ok & length(dens) > 0,
    has_low     = "Low"    %in% dens,
    has_medium  = "Medium" %in% dens,
    has_high    = "High"   %in% dens,
    max_density = max_cat
  )

  if (i %% pb_step == 0 || i == length(all_dates)) {
    cat(glue("  {i}/{length(all_dates)}  ({format(d)})\r"))
  }
}
cat("\n\n")

daily <- bind_rows(records) |>
  mutate(
    year  = year(date),
    month = month(date),
    doy   = yday(date),
    max_density = factor(max_density, levels = DENSITY_LEVELS)
  )

write_csv(daily, file.path(PROC_DIR, "tahoe_hms_smoke_daily.csv"))
cat(glue("Saved: {file.path(PROC_DIR, 'tahoe_hms_smoke_daily.csv')}\n"))

# =============================================================================
# SECTION 5: Annual and monthly summaries
# =============================================================================

annual_avail <- daily |>
  group_by(year) |>
  summarise(n_files = sum(file_avail), .groups = "drop") |>
  mutate(days_in_year = as.integer(format(as.Date(glue("{year}-12-31")), "%j")),
         pct_complete = round(n_files / days_in_year * 100, 0))

annual <- daily |>
  filter(smoke_day) |>
  count(year, max_density, name = "n_smoke_days") |>
  complete(year = YEARS,
           max_density = factor(DENSITY_LEVELS, levels = DENSITY_LEVELS),
           fill = list(n_smoke_days = 0L)) |>
  left_join(annual_avail, by = "year")

annual_totals <- annual |>
  group_by(year) |>
  summarise(total_smoke_days = sum(n_smoke_days), .groups = "drop") |>
  left_join(annual_avail, by = "year")

write_csv(annual,        file.path(PROC_DIR, "tahoe_hms_smoke_annual.csv"))
cat(glue("Saved: {file.path(PROC_DIR, 'tahoe_hms_smoke_annual.csv')}\n"))

monthly <- daily |>
  filter(smoke_day) |>
  count(year, month, name = "n_smoke_days") |>
  complete(year = YEARS, month = 1:12, fill = list(n_smoke_days = 0L))

write_csv(monthly, file.path(PROC_DIR, "tahoe_hms_smoke_monthly.csv"))
cat(glue("Saved: {file.path(PROC_DIR, 'tahoe_hms_smoke_monthly.csv')}\n\n"))

cat("--- Annual lake smoke-day summary ---\n")
annual |>
  pivot_wider(id_cols = year, names_from = max_density,
              values_from = n_smoke_days, values_fill = 0L) |>
  left_join(annual_totals |> select(year, total_smoke_days, pct_complete),
            by = "year") |>
  arrange(year) |>
  print(n = Inf)

# =============================================================================
# SECTION 6: Annual stacked-bar figure
# =============================================================================

annual_plot <- annual |>
  mutate(max_density = factor(max_density, levels = rev(DENSITY_LEVELS)))  # heavy on top

bar_tops <- annual_totals |> rename(bar_top = total_smoke_days)

event_df <- tibble(
  year  = as.integer(names(EVENT_YEARS)),
  label = unname(EVENT_YEARS)
) |>
  filter(year %in% YEARS) |>
  left_join(bar_tops, by = "year")

p_bar <- ggplot(annual_plot,
                aes(x = factor(year), y = n_smoke_days, fill = max_density)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_segment(
    data = event_df,
    aes(x = factor(year), xend = factor(year),
        y = bar_top + 7, yend = bar_top + 1),
    colour = "#333333", linewidth = 0.5,
    arrow = arrow(length = unit(3, "pt"), type = "closed"),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = event_df,
    aes(x = factor(year), y = bar_top + 8, label = label),
    angle = 45, hjust = 0, vjust = 0, size = 3.2,
    colour = "#333333", fontface = "italic", inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = DENSITY_COLS,
    breaks = DENSITY_LEVELS,
    name   = "HMS smoke density"
  ) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Lake Tahoe smoke days from NOAA HMS",
    subtitle = "Days any portion of the lake intersected mapped smoke; binned by max daily density",
    x = NULL, y = "Smoke days per year"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "#555555"),
    legend.position    = "top",
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.4),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1))

p_comp <- ggplot(annual_avail, aes(x = factor(year), y = pct_complete)) +
  geom_col(fill = "#999999", width = 0.72) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "#666666") +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "HMS files available (%)") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.4),
    panel.grid.major.x = element_blank()
  )

p_combined <- p_bar / p_comp + plot_layout(heights = c(3, 1))

ggsave(file.path(FIG_DIR, "tahoe_hms_smoke_days_annual.png"),
       plot = p_combined, width = 8, height = 7, dpi = 200)
cat(glue("\nSaved: {file.path(FIG_DIR, 'tahoe_hms_smoke_days_annual.png')}\n"))

# =============================================================================
# SECTION 7: Monthly heatmap
# =============================================================================

p_heat <- ggplot(monthly,
                 aes(x = factor(month, labels = month.abb),
                     y = factor(year), fill = n_smoke_days)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_smoke_days > 0, n_smoke_days, "")),
            size = 3, colour = "grey20") +
  scale_fill_gradient(low = "#FFFFFF", high = "#B71C1C",
                      name = "Smoke days\nper month", limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  labs(title = "Lake Tahoe monthly smoke-day frequency (NOAA HMS)",
       x = "Month", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggsave(file.path(FIG_DIR, "tahoe_hms_smoke_days_monthly_heatmap.png"),
       plot = p_heat, width = 9, height = 6, dpi = 200)
cat(glue("Saved: {file.path(FIG_DIR, 'tahoe_hms_smoke_days_monthly_heatmap.png')}\n"))

cat("\nDone.\n")
