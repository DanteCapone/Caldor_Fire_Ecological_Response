# Mixed-layer nutrient inventories and Caldor atmospheric-load scaling.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")
source("R/ctd_mixing_helpers.R")

OUT_DIR <- "figures_v2"
PROC_DIR <- file.path("data", "processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
EXPERIMENT_ASH_MG_L <- 25
MAX_MLD_BRACKET_DAYS <- 60
MAX_NEAREST_CTD_DAYS <- 35

# Reuse the project's measured-density, 0.1 kg m-3 threshold definition at MLTP.
CTD_STATION <- "Mid-lake"
CTD_STATION_ALIAS <- "MLTP"
mltp_profiles <- load_ltp_ctd_profiles(".")

mld_casts <- mltp_profiles %>%
  group_by(CTD_ID, Event_ID, date) %>%
  summarise(
    mld_m = threshold_depth(depth_m, Density, DELTA_RHO_MLD, "increase"),
    max_profile_depth_m = max(depth_m),
    n_depth_bins = n_distinct(depth_m),
    .groups = "drop"
  ) %>%
  filter(max_profile_depth_m >= MIN_PROFILE_DEPTH_M, is.finite(mld_m)) %>%
  group_by(date) %>%
  summarise(
    mld_m = mean(mld_m), n_casts = n(),
    max_profile_depth_m = max(max_profile_depth_m), .groups = "drop"
  ) %>%
  arrange(date)

write_csv(mld_casts, file.path(PROC_DIR, "mltp_mixed_layer_depth_density_2005_2025.csv"))

match_mld <- function(target_dates, mld_data = mld_casts) {
  target_dates <- as.Date(target_dates)
  x <- as.numeric(mld_data$date)
  tx <- as.numeric(target_dates)
  interpolated <- approx(x, mld_data$mld_m, xout = tx, rule = 1)$y
  nearest_i <- vapply(tx, function(z) which.min(abs(x - z)), integer(1))
  left_i <- findInterval(tx, x)
  right_i <- pmin(left_i + 1L, length(x))
  left_i[left_i < 1L] <- 1L
  bracket_days <- x[right_i] - x[left_i]
  nearest_days <- abs(tx - x[nearest_i])
  tibble(
    date = target_dates,
    mld_m = interpolated,
    nearest_ctd_date = mld_data$date[nearest_i],
    nearest_ctd_days = nearest_days,
    bracket_days = bracket_days,
    reliable_mld_match = is.finite(interpolated) &
      nearest_days <= MAX_NEAREST_CTD_DAYS & bracket_days <= MAX_MLD_BRACKET_DAYS
  )
}

mld_2021 <- mld_casts %>% filter(year(date) == 2021)
p_mld <- ggplot(mld_2021, aes(date, mld_m)) +
  annotate("rect", xmin = FIRE_START, xmax = FIRE_END, ymin = -Inf, ymax = Inf,
           fill = "#E69F00", alpha = 0.18) +
  geom_line(colour = "#0072B2", linewidth = 0.65) +
  geom_point(shape = 21, fill = "white", colour = "#0072B2", size = 2) +
  scale_y_reverse() +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(x = NULL, y = "Mixed-layer depth (m)",
       title = "Mid-lake mixed-layer depth, 2021") +
  theme_classic(base_size = 8, base_family = "Times New Roman")

ggsave(file.path(OUT_DIR, "mixed_layer_depth_2021.png"), p_mld,
       width = 12.7, height = 8.5, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(OUT_DIR, "mixed_layer_depth_2021.pdf"), p_mld,
       width = 12.7, height = 8.5, units = "cm", device = cairo_pdf)

# Interval loads are existing integrated mg m-2 values; daily-load products are
# retained only for an algebraic duration check.
dep <- read_csv(
  file.path("data", "lake_environmental_data", "deposition",
            "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    start_date = as.Date(Start_Datetime), end_date = as.Date(End_Datetime),
    midpoint = start_date + as.integer(end_date - start_date) %/% 2L,
    interval_days = as.numeric(difftime(End_Datetime, Start_Datetime, units = "days")),
    overlap_start = pmax(start_date, FIRE_START),
    overlap_end = pmin(end_date, FIRE_END),
    overlap_days = pmax(as.numeric(overlap_end - overlap_start), 0),
    event_fraction = pmin(overlap_days / interval_days, 1)
  )

dep_map <- tribble(
  ~nutrient, ~load_col, ~daily_col, ~pool,
  "NO3", "NO3_Load", "NO3_Daily_Load", "N",
  "NH4", "NH4_Load", "NH4_Daily_Load", "N",
  "DIN", "DIN_Load", "DIN_Daily_Load", "N",
  "TKN", "TKN_Load", "TKN_Daily_Load", "N",
  "SRP", "SRP_Load", "SRP_Daily_Load", "P",
  "TP", "TP_Load", "TP_Daily_Load", "P"
)

dep_long <- dep_map %>%
  pmap_dfr(function(nutrient, load_col, daily_col, pool) {
    dep %>% transmute(
      Sample_ID, start_date, end_date, midpoint, interval_days, overlap_days,
      event_fraction, QA_Code, nutrient, pool,
      deposition_mg_m2 = as.numeric(.data[[load_col]]),
      daily_flux_mg_m2_d = as.numeric(.data[[daily_col]]),
      duration_check_ratio = deposition_mg_m2 /
        (daily_flux_mg_m2_d * interval_days)
    )
  }) %>%
  left_join(match_mld(unique(dep$midpoint)), by = c("midpoint" = "date")) %>%
  mutate(
    equivalent_concentration_mg_L = deposition_mg_m2 / (mld_m * 1000),
    event_load_mg_m2 = deposition_mg_m2 * event_fraction,
    event_equivalent_mg_L = event_load_mg_m2 / (mld_m * 1000)
  )

write_csv(dep_long, file.path(PROC_DIR, "caldor_deposition_mixed_layer_equivalent.csv"))

deposition_event <- dep_long %>%
  filter(overlap_days > 0, reliable_mld_match, is.finite(event_load_mg_m2)) %>%
  group_by(nutrient, pool) %>%
  summarise(
    atmospheric_load_mg_m2 = sum(event_load_mg_m2),
    mixed_layer_equivalent_mg_L = sum(event_equivalent_mg_L),
    n_intervals = n(), .groups = "drop"
  )

# Integrate only where the sampled profile brackets the surface-to-MLD layer.
integrate_profile_to_mld <- function(depth_m, concentration_ug_L, mld_m) {
  ok <- is.finite(depth_m) & is.finite(concentration_ug_L)
  d <- depth_m[ok]
  v <- concentration_ug_L[ok]
  ord <- order(d)
  d <- d[ord]; v <- v[ord]
  if (length(d) < 2 || min(d) > 2 || max(d) < mld_m || mld_m <= min(d)) {
    return(tibble(inventory_mg_m2 = NA_real_,
                  max_sample_depth_m = if (length(d)) max(d) else NA_real_,
                  n_profile_depths = length(d), profile_brackets_mld = FALSE))
  }
  boundary <- approx(d, v, xout = mld_m, rule = 1)$y
  keep <- d < mld_m
  d_use <- c(d[keep], mld_m)
  v_use <- c(v[keep], boundary)
  inventory <- sum(diff(d_use) * (head(v_use, -1) + tail(v_use, -1)) / 2)
  tibble(inventory_mg_m2 = inventory, max_sample_depth_m = max(d),
         n_profile_depths = length(d_use), profile_brackets_mld = TRUE)
}

nut_raw <- read_csv(
  file.path("data", "lake_environmental_data", "nutrients", "Tahoe_MLTP_Nutrient.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(Date), depth_m = as.numeric(Depth), DIN = as.numeric(NO3) + as.numeric(NH4))

nut_long <- nut_raw %>%
  transmute(date, depth_m, NO3, NH4, DIN, TKN, SRP = TRP, TP = THP) %>%
  pivot_longer(-c(date, depth_m), names_to = "nutrient", values_to = "concentration_ug_L") %>%
  mutate(pool = if_else(nutrient %in% c("SRP", "TP"), "P", "N")) %>%
  left_join(match_mld(unique(nut_raw$date)), by = "date")

inventory <- nut_long %>%
  filter(reliable_mld_match) %>%
  group_by(date, nutrient, pool, mld_m, nearest_ctd_date, nearest_ctd_days, bracket_days) %>%
  group_modify(~integrate_profile_to_mld(.x$depth_m, .x$concentration_ug_L,
                                         unique(.y$mld_m))) %>%
  ungroup() %>%
  mutate(year = year(date), month = month(date))

inventory_climatology <- inventory %>%
  filter(year != 2021, profile_brackets_mld, is.finite(inventory_mg_m2)) %>%
  group_by(nutrient, month) %>%
  summarise(
    expected_inventory_mg_m2 = median(inventory_mg_m2),
    historical_n = n(), .groups = "drop"
  )

inventory_2021 <- inventory %>%
  filter(year == 2021) %>%
  left_join(inventory_climatology, by = c("nutrient", "month")) %>%
  mutate(
    inventory_anomaly_mg_m2 = inventory_mg_m2 - expected_inventory_mg_m2,
    proportional_anomaly = inventory_anomaly_mg_m2 / expected_inventory_mg_m2
  )

write_csv(inventory, file.path(PROC_DIR, "mltp_mixed_layer_nutrient_inventories.csv"))
write_csv(inventory_2021, file.path(PROC_DIR, "mltp_mixed_layer_nutrient_anomalies_2021.csv"))

budget_summary <- inventory_2021 %>%
  filter(date >= FIRE_START, date <= FIRE_END, profile_brackets_mld) %>%
  group_by(nutrient, pool) %>%
  summarise(
    historical_inventory_mg_m2 = median(expected_inventory_mg_m2, na.rm = TRUE),
    observed_2021_inventory_mg_m2 = median(inventory_mg_m2, na.rm = TRUE),
    inventory_anomaly_mg_m2 = observed_2021_inventory_mg_m2 - historical_inventory_mg_m2,
    proportional_anomaly = inventory_anomaly_mg_m2 / historical_inventory_mg_m2,
    n_2021_profiles = n(), .groups = "drop"
  ) %>%
  left_join(deposition_event, by = c("nutrient", "pool")) %>%
  mutate(
    atmospheric_pct_historical_inventory = 100 * atmospheric_load_mg_m2 /
      historical_inventory_mg_m2,
    anomaly_direction = if_else(inventory_anomaly_mg_m2 >= 0,
                                "apparent excess", "apparent deficit"),
    atmospheric_to_inventory_anomaly = atmospheric_load_mg_m2 / inventory_anomaly_mg_m2,
    atmospheric_to_abs_inventory_anomaly = atmospheric_load_mg_m2 /
      abs(inventory_anomaly_mg_m2)
  )

write_csv(budget_summary, file.path(PROC_DIR, "caldor_mixed_layer_nutrient_budget_summary.csv"))

# No field ash-deposition measurements or experimental ash chemistry are present
# in the project. The supported inverse calculation is retained for every 2021 MLD.
ash_scaling <- mld_2021 %>%
  filter(date >= FIRE_START, date <= FIRE_END) %>%
  transmute(
    date, mld_m,
    required_ash_mg_m2 = EXPERIMENT_ASH_MG_L * mld_m * 1000,
    required_ash_g_m2 = required_ash_mg_m2 / 1000,
    field_ash_equivalent_mg_L = NA_real_,
    experiment_to_field_ratio = NA_real_,
    limitation = "Field ash-mass deposition dataset not available in project"
  )
write_csv(ash_scaling, file.path(PROC_DIR, "ash_experiment_required_areal_deposition.csv"))

budget_plot_data <- budget_summary %>%
  select(nutrient, pool, atmospheric_load_mg_m2, historical_inventory_mg_m2,
         observed_2021_inventory_mg_m2) %>%
  pivot_longer(ends_with("mg_m2"), names_to = "component", values_to = "value") %>%
  mutate(
    component = recode(component,
      atmospheric_load_mg_m2 = "Caldor atmospheric input",
      historical_inventory_mg_m2 = "Historical mixed-layer inventory",
      observed_2021_inventory_mg_m2 = "Observed 2021 inventory"),
    component = factor(component, levels = c("Historical mixed-layer inventory",
      "Observed 2021 inventory", "Caldor atmospheric input"))
  )

p_budget <- ggplot(budget_plot_data, aes(component, value, fill = component)) +
  geom_col(width = 0.72) +
  facet_grid(pool ~ nutrient, scales = "free_y", space = "free_x") +
  scale_fill_manual(values = c(
    "Historical mixed-layer inventory" = "grey65",
    "Observed 2021 inventory" = "#0072B2",
    "Caldor atmospheric input" = "#D55E00"), guide = "none") +
  labs(x = NULL, y = expression("Areal nutrient mass (mg m"^{-2}*")"),
       title = "Caldor atmospheric input relative to mixed-layer inventories") +
  theme_classic(base_size = 7, base_family = "Times New Roman") +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, size = 5.8),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold"))

ash_q <- quantile(ash_scaling$required_ash_g_m2, c(0, 0.5, 1), na.rm = TRUE)
p_ash <- ggplot(tibble(x = 1, y = ash_q[2], ymin = ash_q[1], ymax = ash_q[3]), aes(x, y)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), linewidth = 0.8, colour = "#009E73") +
  geom_point(shape = 21, fill = "white", colour = "#009E73", size = 2.5) +
  annotate("text", x = 1, y = ash_q[3], label = "MLD range", vjust = -0.6, size = 2.1) +
  scale_x_continuous(NULL, breaks = 1, labels = "25 mg ash L^-1") +
  labs(y = expression("Required ash deposition (g m"^{-2}*")"),
       title = "Experimental ash dose: inverse field scaling") +
  theme_classic(base_size = 7, base_family = "Times New Roman")

mixed_layer_budget_figure <- p_budget / p_ash + plot_layout(heights = c(3.4, 1.4))
ggsave(file.path(OUT_DIR, "mixed_layer_nutrient_budget.png"), mixed_layer_budget_figure,
       width = 17.8, height = 15.2, units = "cm", dpi = 600, device = ragg::agg_png)
ggsave(file.path(OUT_DIR, "mixed_layer_nutrient_budget.pdf"), mixed_layer_budget_figure,
       width = 17.8, height = 15.2, units = "cm", device = cairo_pdf)

duration_error_mg_m2 <- dep_long %>%
  filter(is.finite(duration_check_ratio)) %>%
  summarise(max_abs_error = max(abs(deposition_mg_m2 -
    daily_flux_mg_m2_d * interval_days))) %>% pull()
# Daily fluxes are rounded in the source table; compare in absolute load units.
if (duration_error_mg_m2 > 0.2) {
  warning("Integrated and rounded daily deposition loads differ by >0.2 mg m-2.")
}
if (any(abs(dep$DIN_Load - dep$NO3_Load - dep$NH4_Load) > 0.05, na.rm = TRUE)) {
  warning("DIN_Load is not always NO3_Load + NH4_Load within rounding tolerance.")
}

message("Mixed-layer budget outputs written to figures_v2/ and data/processed/.")
invisible(TRUE)
