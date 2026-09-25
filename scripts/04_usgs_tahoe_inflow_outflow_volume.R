# =============================================================================
# 04_usgs_tahoe_inflow_outflow_volume.R
#
# Convert USGS daily mean discharge to daily water volume and compare two
# monitored Lake Tahoe tributary inputs with the Truckee River outlet.
#
# Daily volume (m3) = daily mean discharge (m3 s-1) * 86,400 s day-1.
# Only matched days with all three gauges are compared. Upper Truckee River +
# Blackwood Creek are a partial monitored inflow, not total Lake Tahoe inflow.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(ragg)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
  proj_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
} else {
  proj_root <- normalizePath(".", winslash = "/")
}
setwd(proj_root)
source("scripts/figure_aesthetics.R")

proc_dir <- file.path(proj_root, "data", "processed")
fig_dir <- file.path(proj_root, "figures", "streamflow")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

focus_start <- as.Date("2021-07-01")
focus_end <- as.Date("2022-07-31")
fire_start <- as.Date("2021-08-14")
fire_end <- as.Date("2021-10-21")

sites <- tribble(
  ~site_no,   ~site,                     ~flow_role,
  "10336610", "Upper Truckee River",     "Monitored inflow",
  "10336660", "Blackwood Creek",         "Monitored inflow",
  "10337500", "Truckee River at outlet", "Outlet"
)

cache_file <- file.path(
  proc_dir, "usgs_tahoe_inflow_outflow_daily_site_volume.csv"
)
usgs_url <- paste0(
  "https://waterservices.usgs.gov/nwis/dv/?format=rdb",
  "&sites=", paste(sites$site_no, collapse = "%2C"),
  "&startDT=", focus_start, "&endDT=", focus_end,
  "&parameterCd=00060&statCd=00003&siteStatus=all"
)

message("Downloading USGS daily mean discharge for volume comparison...")
raw <- tryCatch(
  readr::read_tsv(
    usgs_url, comment = "#", col_types = cols(.default = col_character()),
    progress = FALSE, name_repair = "minimal", show_col_types = FALSE
  ),
  error = function(e) {
    if (!file.exists(cache_file)) stop(e)
    warning("USGS download failed; using cached data: ", conditionMessage(e))
    NULL
  }
)

if (is.null(raw)) {
  daily_site <- read_csv(cache_file, show_col_types = FALSE) %>%
    mutate(date = as.Date(date))
} else {
  flow_col <- names(raw)[
    grepl("00060_00003$", names(raw)) & !grepl("_cd$", names(raw))
  ]
  qual_col <- names(raw)[grepl("00060_00003_cd$", names(raw))]
  if (length(flow_col) != 1L) {
    stop("Could not uniquely identify daily discharge: ",
         paste(names(raw), collapse = ", "))
  }

  daily_site <- raw %>%
    filter(agency_cd == "USGS", site_no %in% sites$site_no) %>%
    transmute(
      site_no,
      date = as.Date(datetime),
      flow_cfs = parse_double(.data[[flow_col]]),
      qualifier = if (length(qual_col) == 1L) .data[[qual_col]] else NA_character_
    ) %>%
    left_join(sites, by = "site_no") %>%
    mutate(
      flow_cms = flow_cfs * 0.028316846592,
      daily_volume_m3 = flow_cms * 86400
    ) %>%
    filter(date >= focus_start, date <= focus_end) %>%
    arrange(date, site)

  if (nrow(daily_site) == 0L || all(is.na(daily_site$daily_volume_m3))) {
    stop("USGS returned no usable daily discharge observations.")
  }
  write_csv(daily_site, cache_file, na = "")
}

# Complete the date x gauge grid so missing gauge-days cannot silently bias
# either side of the comparison.
volume_grid <- crossing(
  date = seq(focus_start, focus_end, by = "day"),
  sites
) %>%
  left_join(
    daily_site %>%
      select(site_no, date, flow_cfs, flow_cms, daily_volume_m3, qualifier),
    by = c("site_no", "date")
  )

daily_comparison <- volume_grid %>%
  group_by(date) %>%
  summarise(
    n_inflow_gauges = sum(
      flow_role == "Monitored inflow" & is.finite(daily_volume_m3)
    ),
    n_outlet_gauges = sum(
      flow_role == "Outlet" & is.finite(daily_volume_m3)
    ),
    monitored_inflow_volume_m3 = if (n_inflow_gauges == 2L) {
      sum(daily_volume_m3[flow_role == "Monitored inflow"], na.rm = TRUE)
    } else NA_real_,
    outlet_volume_m3 = if (n_outlet_gauges == 1L) {
      sum(daily_volume_m3[flow_role == "Outlet"], na.rm = TRUE)
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    complete_matched_day = !is.na(monitored_inflow_volume_m3) &
      !is.na(outlet_volume_m3),
    monitored_net_volume_m3 = monitored_inflow_volume_m3 - outlet_volume_m3,
    outlet_to_monitored_inflow_ratio = outlet_volume_m3 /
      monitored_inflow_volume_m3
  )

write_csv(
  daily_comparison,
  file.path(proc_dir, "usgs_tahoe_inflow_outflow_daily_comparison.csv"),
  na = ""
)

monthly_comparison <- daily_comparison %>%
  mutate(month_date = floor_date(date, "month")) %>%
  group_by(month_date) %>%
  summarise(
    n_calendar_days = n(),
    n_matched_days = sum(complete_matched_day),
    matched_day_coverage_pct = 100 * n_matched_days / n_calendar_days,
    monitored_inflow_volume_m3 = sum(
      monitored_inflow_volume_m3[complete_matched_day], na.rm = TRUE
    ),
    outlet_volume_m3 = sum(
      outlet_volume_m3[complete_matched_day], na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    full_month_coverage = n_matched_days == n_calendar_days,
    monitored_net_volume_m3 = monitored_inflow_volume_m3 - outlet_volume_m3,
    outlet_to_monitored_inflow_ratio = outlet_volume_m3 /
      monitored_inflow_volume_m3
  )

write_csv(
  monthly_comparison,
  file.path(proc_dir, "usgs_tahoe_inflow_outflow_monthly_comparison.csv"),
  na = ""
)

matched <- daily_comparison %>% filter(complete_matched_day)
volume_summary <- tibble(
  focus_start = focus_start,
  focus_end = focus_end,
  n_calendar_days = nrow(daily_comparison),
  n_matched_days = nrow(matched),
  matched_day_coverage_pct = 100 * n_matched_days / n_calendar_days,
  monitored_inflow_volume_m3 = sum(matched$monitored_inflow_volume_m3),
  outlet_volume_m3 = sum(matched$outlet_volume_m3),
  monitored_net_volume_m3 = monitored_inflow_volume_m3 - outlet_volume_m3,
  outlet_to_monitored_inflow_ratio = outlet_volume_m3 /
    monitored_inflow_volume_m3,
  comparison_scope = paste(
    "Upper Truckee River + Blackwood Creek versus Tahoe City outlet;",
    "partial monitored inflow, not a whole-lake water balance"
  )
)
write_csv(
  volume_summary,
  file.path(proc_dir, "usgs_tahoe_inflow_outflow_summary.csv"),
  na = ""
)

plot_data <- monthly_comparison %>%
  transmute(
    month_date,
    `Monitored tributary inflow` = monitored_inflow_volume_m3,
    `Tahoe City outlet` = outlet_volume_m3
  ) %>%
  pivot_longer(-month_date, names_to = "volume_type", values_to = "volume_m3") %>%
  mutate(volume_type = factor(
    volume_type,
    levels = c("Monitored tributary inflow", "Tahoe City outlet")
  ))

p_volume <- ggplot(plot_data, aes(month_date, volume_m3, fill = volume_type)) +
  annotate("rect", xmin = fire_start, xmax = fire_end,
           ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.07) +
  geom_col(position = position_dodge2(width = 25, preserve = "single"),
           width = 23, colour = NA) +
  geom_vline(xintercept = fire_start, colour = "firebrick",
             linetype = "dashed", linewidth = LO_LW_MID) +
  scale_fill_manual(
    values = c("Monitored tributary inflow" = "#2E7D32",
               "Tahoe City outlet" = "#1565C0"),
    name = NULL
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y",
               expand = expansion(mult = c(0.015, 0.02))) +
  scale_y_continuous(labels = label_number(scale = 1e-6, accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(
    x = NULL,
    y = expression("Matched-day water volume (million m"^3*")"),
    caption = paste(
      "Inflow includes Upper Truckee River and Blackwood Creek only.",
      "Daily mean Q was multiplied by 86,400 s; comparisons use matched days."
    )
  ) +
  lo_theme(base_size = 8, family = LO_FONT) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7),
        plot.caption = element_text(size = 6.5, hjust = 0))

save_lo_fig(
  p_volume,
  file.path(fig_dir, "tahoe_monitored_inflow_vs_outflow_volume.png"),
  width_type = "double", height_cm = 9.5
)

cat("\nMatched-day volume comparison:\n")
print(volume_summary)
cat("\nGuardrail: the two gauges represent partial monitored inflow only.\n")
