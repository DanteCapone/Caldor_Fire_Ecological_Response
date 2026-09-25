# Reproducible values requested for manuscript Results text.
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

PHYTO_FILE <- file.path("data", "lake_environmental_data", "phytoplankton",
                        "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv")
CHL_RECENT <- file.path("data", "lake_environmental_data", "chla", "Tahoe_LTP_Chl.csv")
CHL_HIST <- file.path("data", "lake_environmental_data", "chla", "terc_chla_all.csv")
STABILITY_FILE <- file.path("data", "processed", "ctd", "mltp_lake_tools_stability_metrics.csv")
FIG4_DEEP_TEST <- file.path("figures", "figure_4_phytoplankton",
                            "figure_4_ltp_chla_depth_bin_sum_tests.csv")
FIG5_STATE_SPACE <- file.path("figures", "figure_5_nmds_braycurtis",
                              "figure_5_state_space_annual_summary.csv")
SMOKE_ANNUAL <- file.path("data", "processed",
                          "supplemental_smoke_days_annual_summary.csv")
NUTRIENT_2020_TESTS <- file.path(
  "figures", "supplemental", "nutrients_2020_historical",
  "supplemental_2020_nutrients_mann_whitney_results.csv"
)
OUT_DIR <- Sys.getenv(
  "CALDOR_MANUSCRIPT_STATS_OUT_DIR",
  unset = file.path("figures", "manuscript_statistics")
)
WRITE_EXISTING_OUTPUTS <- tolower(Sys.getenv(
  "CALDOR_MANUSCRIPT_STATS_WRITE_EXISTING",
  unset = "true"
)) %in% c("true", "1", "yes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

trap_int <- function(depth, value) {
  ord <- order(depth)
  d <- depth[ord]
  v <- value[ord]
  ok <- is.finite(d) & is.finite(v)
  d <- d[ok]
  v <- v[ok]
  if (length(d) < 2L) return(NA_real_)
  sum(diff(d) * (head(v, -1) + tail(v, -1)) / 2)
}

# ---- Chlorophyll: deep October anomaly and Aug 2021-Feb 2022 integration ----
deep_tests <- read_csv(FIG4_DEEP_TEST, show_col_types = FALSE) %>%
  mutate(
    sd_anomaly = (obs_sum - clim_sum) / clim_sd,
    t_simple = t_pred,
    p_simple = p_value,
    q_fdr = p.adjust(p_simple, method = "BH")
  )
deep_oct <- deep_tests %>% filter(month == 10L, depth_zone == "60-105 m")

chl_recent <- read_csv(CHL_RECENT, show_col_types = FALSE) %>%
  transmute(station = "LTP", date = as.Date(Date), depth = as.numeric(Depth),
            chla = as.numeric(Chla)) %>%
  filter(is.finite(depth), is.finite(chla), chla > 0)
chl_hist <- read.csv(CHL_HIST, quote = "", stringsAsFactors = FALSE,
                     fill = TRUE, comment.char = "") %>%
  as_tibble() %>%
  filter(Sample_Type %in% c("FIELD", "FLDDUP"), Station_ID == "Index") %>%
  transmute(station = "LTP", date = as.Date(Date),
            depth = suppressWarnings(as.numeric(Depth)),
            chla = suppressWarnings(as.numeric(Chla))) %>%
  filter(is.finite(depth), is.finite(chla), chla > 0)
chl <- bind_rows(anti_join(chl_hist, chl_recent, by = c("station", "date", "depth")),
                 chl_recent) %>%
  mutate(year = year(date), month = month(date)) %>%
  filter(depth <= 150, chla < 50)

target_months <- tibble(
  year = c(rep(2021L, 5), rep(2022L, 2)),
  month = c(8:12, 1:2)
)
obs_integrated <- chl %>%
  inner_join(target_months, by = c("year", "month")) %>%
  group_by(year, month, depth) %>%
  summarise(chla = mean(chla, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, month) %>%
  summarise(integrated_chla_mg_m2 = trap_int(depth, chla), .groups = "drop")

hist_integrated <- chl %>%
  anti_join(target_months, by = c("year", "month")) %>%
  filter(month %in% 1:2 | month %in% 8:12) %>%
  group_by(year, month, depth) %>%
  summarise(chla = mean(chla, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, month) %>%
  summarise(int_year = trap_int(depth, chla), .groups = "drop") %>%
  filter(is.finite(int_year))

integrated_tests <- hist_integrated %>%
  group_by(month) %>%
  summarise(hist_mean = mean(int_year), hist_sd = sd(int_year), hist_n = n(),
            .groups = "drop") %>%
  right_join(obs_integrated, by = "month") %>%
  mutate(
    comparison_se = hist_sd * sqrt(1 + 1 / hist_n),
    t_simple = (integrated_chla_mg_m2 - hist_mean) / comparison_se,
    p_simple = 2 * pt(-abs(t_simple), df = hist_n - 1),
    anomaly_sd = (integrated_chla_mg_m2 - hist_mean) / hist_sd
  ) %>%
  mutate(q_fdr = p.adjust(p_simple, method = "BH")) %>%
  arrange(year, month)

paired_window <- t.test(integrated_tests$integrated_chla_mg_m2,
                        integrated_tests$hist_mean, paired = TRUE,
                        alternative = "greater")
paired_window_two_sided <- t.test(integrated_tests$integrated_chla_mg_m2,
                                  integrated_tests$hist_mean, paired = TRUE,
                                  alternative = "two.sided")

# ---- Leptolyngbya maxima and proportional contribution ---------------------
phyto <- read_csv(PHYTO_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), depth = as.numeric(depth_num),
         abundance = replace_na(as.numeric(abundance), 0),
         biovolume = replace_na(as.numeric(biovolume), 0),
         is_lepto = str_detect(taxon, regex("^Leptolyngbya", ignore_case = TRUE)))

sample_totals <- phyto %>%
  group_by(date, depth_bin, depth) %>%
  summarise(
    total_abundance = sum(abundance), total_biovolume = sum(biovolume),
    lepto_abundance = sum(abundance[is_lepto]),
    lepto_biovolume = sum(biovolume[is_lepto]), .groups = "drop"
  ) %>%
  mutate(abundance_pct = 100 * lepto_abundance / total_abundance,
         biovolume_pct = 100 * lepto_biovolume / total_biovolume)

peak_upper <- sample_totals %>%
  filter(depth_bin == "0-40 m", date >= as.Date("2021-08-14")) %>%
  slice_max(lepto_abundance, n = 1, with_ties = FALSE)

upper_date_totals <- phyto %>%
  filter(depth_bin == "0-40 m") %>%
  group_by(date) %>%
  summarise(total_abundance = sum(abundance), total_biovolume = sum(biovolume),
            lepto_abundance = sum(abundance[is_lepto]),
            lepto_biovolume = sum(biovolume[is_lepto]), .groups = "drop") %>%
  mutate(abundance_pct = 100 * lepto_abundance / total_abundance,
         biovolume_pct = 100 * lepto_biovolume / total_biovolume)
peak_upper_date <- upper_date_totals %>%
  filter(date == peak_upper$date) %>% slice_head(n = 1)

# The manuscript's reported 2011 values are sums across the observed depths in
# each LTP depth zone on a sampling date. Retain that definition when adding the
# corresponding dates and community-biovolume contributions.
historical_2011_maxima <- phyto %>%
  filter(year(date) == 2011L, depth_bin %in% c("0-40 m", "60-105 m")) %>%
  group_by(date, depth_bin) %>%
  summarise(
    sampled_depths_n = n_distinct(depth),
    total_abundance = sum(abundance),
    total_biovolume = sum(biovolume),
    lepto_abundance = sum(abundance[is_lepto]),
    lepto_biovolume = sum(biovolume[is_lepto]),
    .groups = "drop"
  ) %>%
  mutate(
    abundance_pct = 100 * lepto_abundance / total_abundance,
    biovolume_pct = 100 * lepto_biovolume / total_biovolume
  ) %>%
  group_by(depth_bin) %>%
  slice_max(lepto_abundance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(factor(depth_bin, levels = c("0-40 m", "60-105 m"))) %>%
  transmute(
    year = year(date), date, depth_zone = depth_bin, sampled_depths_n,
    lepto_abundance_depth_sum_cells_l = lepto_abundance,
    community_abundance_depth_sum_cells_l = total_abundance,
    lepto_abundance_percent = abundance_pct,
    lepto_biovolume_depth_sum_mm3_l = lepto_biovolume,
    community_biovolume_depth_sum_mm3_l = total_biovolume,
    lepto_biovolume_percent = biovolume_pct,
    aggregation_definition = paste(
      "sum across observed LTP depths in the stated zone on each date;",
      "row is the 2011 date with maximum Leptolyngbya abundance"
    ),
    inference = "descriptive; no significance test"
  )

historical_2011_shallow <- historical_2011_maxima %>%
  filter(depth_zone == "0-40 m")
historical_2011_deep <- historical_2011_maxima %>%
  filter(depth_zone == "60-105 m")

# ---- Figure 5b focal-year environmental and biological contrasts -----------
figure5b_annual <- read_csv(FIG5_STATE_SPACE, show_col_types = FALSE)
smoke_annual <- read_csv(SMOKE_ANNUAL, show_col_types = FALSE)
nutrient_2020_tests <- read_csv(NUTRIENT_2020_TESTS, show_col_types = FALSE)

figure5b_focal_years <- figure5b_annual %>%
  filter(year %in% c(2011L, 2020L, 2021L)) %>%
  select(
    year, mixed_layer_thermal_contrast_c, mixed_layer_depth_m,
    surface_temperature_0_10m_c, schmidt_stability_j_m2, antecedent_pc1,
    qualifying_smoke_days, panel_b_qualifying_smoke_days,
    integrated_pm25_aug_dec_ug_d_m3,
    biological_n_dates, biological_span_days,
    biological_temporal_auc_cells_d_m2,
    biological_maximum_cells_m2, biological_coverage_adequate,
    biological_observed_zero
  ) %>%
  left_join(
    smoke_annual %>%
      select(
        year,
        full_year_smoke_days = smoke_days_pm25_ge_9_1,
        full_year_smoke_days_ge_35_5 = `smoke_days_pm25_>=35.5`,
        annual_max_pm25_ug_m3
      ),
    by = "year"
  ) %>%
  mutate(
    panel_b_pm25_definition = paste(
      "August 1-December 31 sum of daily PM2.5 on NOAA HMS smoke days",
      "with composite daily PM2.5 >= 9.1 ug m^-3"
    ),
    biological_response_definition = paste(
      "August 1-December 31 temporal AUC of date-level 5-90 m",
      "depth-integrated Leptolyngbya abundance"
    ),
    inference = "descriptive; no exposure-by-state interaction or causal model"
  ) %>%
  arrange(year)

get_focal_value <- function(target_year, variable) {
  figure5b_focal_years %>%
    filter(year == target_year) %>%
    pull({{ variable }})
}

deposition_2020 <- nutrient_2020_tests %>% filter(source == "Deposition")
deposition_2020_nutrients <- deposition_2020 %>%
  filter(metric %in% c("NO3", "NH4", "TKN", "SRP", "TP")) %>%
  mutate(median_2020_fold_of_month_matched_history = 2^median_2020_log_fold)

figure5b_comparisons <- tribble(
  ~statistic, ~value, ~units_or_detail,
  "2020 full-calendar-year qualifying smoke days", get_focal_value(2020L, full_year_smoke_days), "days; Fig. S1 definition",
  "2020 full-calendar-year smoke days at or above 35.5 ug m^-3", get_focal_value(2020L, full_year_smoke_days_ge_35_5), "days; Fig. S1 definition",
  "2020 annual maximum daily PM2.5", get_focal_value(2020L, annual_max_pm25_ug_m3), "ug m^-3",
  "2011 Figure 5b integrated PM2.5", get_focal_value(2011L, integrated_pm25_aug_dec_ug_d_m3), "ug d m^-3; August-December",
  "2020 Figure 5b integrated PM2.5", get_focal_value(2020L, integrated_pm25_aug_dec_ug_d_m3), "ug d m^-3; August-December",
  "2021 Figure 5b integrated PM2.5", get_focal_value(2021L, integrated_pm25_aug_dec_ug_d_m3), "ug d m^-3; August-December",
  "2021-to-2020 Figure 5b integrated PM2.5 ratio", get_focal_value(2021L, integrated_pm25_aug_dec_ug_d_m3) / get_focal_value(2020L, integrated_pm25_aug_dec_ug_d_m3), "ratio",
  "2011 antecedent thermal-state PC1", get_focal_value(2011L, antecedent_pc1), "dimensionless standardized PCA score",
  "2020 antecedent thermal-state PC1", get_focal_value(2020L, antecedent_pc1), "dimensionless standardized PCA score",
  "2021 antecedent thermal-state PC1", get_focal_value(2021L, antecedent_pc1), "dimensionless standardized PCA score",
  "2011 Leptolyngbya temporal AUC", get_focal_value(2011L, biological_temporal_auc_cells_d_m2), "cells d m^-2; August-December",
  "2020 Leptolyngbya temporal AUC", get_focal_value(2020L, biological_temporal_auc_cells_d_m2), "cells d m^-2; adequately sampled observed zero",
  "2021 Leptolyngbya temporal AUC", get_focal_value(2021L, biological_temporal_auc_cells_d_m2), "cells d m^-2; August-December",
  "2021-to-2011 Leptolyngbya temporal AUC ratio", get_focal_value(2021L, biological_temporal_auc_cells_d_m2) / get_focal_value(2011L, biological_temporal_auc_cells_d_m2), "ratio",
  "Minimum 2020 deposition-nutrient median fold of month-matched history", min(deposition_2020_nutrients$median_2020_fold_of_month_matched_history), "fold; SRP",
  "Maximum 2020 deposition-nutrient median fold of month-matched history", max(deposition_2020_nutrients$median_2020_fold_of_month_matched_history), "fold; TP",
  "Minimum BH-adjusted p among 2020 deposition comparisons", min(deposition_2020$p_adjusted_bh), "q; six deposition metrics",
  "BH-significant 2020 nutrient comparisons", sum(nutrient_2020_tests$p_adjusted_bh < 0.05), "count among 12 deposition and in-lake comparisons"
)

post_spring <- sample_totals %>%
  filter(date > as.Date("2022-05-31"), lepto_abundance > 0)
post_spring_max_abundance_pct <- post_spring %>%
  slice_max(abundance_pct, n = 1, with_ties = FALSE)
post_spring_max_biovolume_pct <- post_spring %>%
  slice_max(biovolume_pct, n = 1, with_ties = FALSE)
post_spring_upper <- upper_date_totals %>%
  filter(date > as.Date("2022-05-31"), lepto_abundance > 0)
post_spring_upper_abundance <- post_spring_upper %>%
  slice_max(abundance_pct, n = 1, with_ties = FALSE)
post_spring_upper_biovolume <- post_spring_upper %>%
  slice_max(biovolume_pct, n = 1, with_ties = FALSE)

# ---- Schmidt stability ------------------------------------------------------
stability <- read_csv(STABILITY_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>% filter(year(date) == 2021L)
oct_stability <- stability %>% filter(month(date) == 10L) %>% slice_min(abs(day(date) - 1L), n = 1)
dec_stability <- stability %>% filter(month(date) == 12L) %>% slice_min(abs(day(date) - 15L), n = 1)
stability_decline_pct <- 100 * (oct_stability$schmidt_stability_j_m2 -
                                  dec_stability$schmidt_stability_j_m2) /
  oct_stability$schmidt_stability_j_m2

stats <- tribble(
  ~statistic, ~value, ~units_or_detail,
  "Deep October 2021 chlorophyll anomaly", as.character(deep_oct$sd_anomaly), "SD",
  "Deep October 2021 chlorophyll t", as.character(deep_oct$t_simple), paste0("df = ", deep_oct$clim_n - 1),
  "Deep October 2021 chlorophyll p", as.character(deep_oct$p_simple), "two-sided comparison with historical October values",
  "Deep October 2021 chlorophyll FDR q", as.character(deep_oct$q_fdr), "Benjamini-Hochberg across 12 Figure 4 month-by-depth comparisons",
  "February 2022 integrated chlorophyll anomaly", as.character(integrated_tests$anomaly_sd[integrated_tests$month == 2L]), "SD",
  "February 2022 integrated chlorophyll t", as.character(integrated_tests$t_simple[integrated_tests$month == 2L]), paste0("df = ", integrated_tests$hist_n[integrated_tests$month == 2L] - 1),
  "February 2022 integrated chlorophyll p", as.character(integrated_tests$p_simple[integrated_tests$month == 2L]), "two-sided comparison with historical February values",
  "February 2022 integrated chlorophyll FDR q", as.character(integrated_tests$q_fdr[integrated_tests$month == 2L]), "Benjamini-Hochberg across 7 months",
  "Integrated chlorophyll Aug 2021-Feb 2022 paired t", as.character(unname(paired_window$statistic)), "df = 6; one-sided",
  "Integrated chlorophyll Aug 2021-Feb 2022 paired p", as.character(paired_window$p.value), "one-sided paired t-test",
  "Integrated chlorophyll Aug 2021-Feb 2022 paired two-sided p", as.character(paired_window_two_sided$p.value), "two-sided paired t-test",
  "Integrated chlorophyll Aug 2021-Feb 2022 mean difference", as.character(unname(paired_window_two_sided$estimate)), "mg m^-2",
  "Integrated chlorophyll Aug 2021-Feb 2022 difference 95% CI lower", as.character(paired_window_two_sided$conf.int[1]), "mg m^-2",
  "Integrated chlorophyll Aug 2021-Feb 2022 difference 95% CI upper", as.character(paired_window_two_sided$conf.int[2]), "mg m^-2",
  "Peak upper-water Leptolyngbya abundance", as.character(peak_upper$lepto_abundance), "cells L^-1",
  "Peak upper-water Leptolyngbya date", format(peak_upper$date, "%Y-%m-%d"), paste0("depth = ", peak_upper$depth, " m"),
  "Leptolyngbya biovolume share at peak sample", as.character(peak_upper$biovolume_pct), "%",
  "Leptolyngbya biovolume share in upper zone on peak date", as.character(peak_upper_date$biovolume_pct), "%",
  "2011 shallow Leptolyngbya maximum date", format(historical_2011_shallow$date, "%Y-%m-%d"), "0-40 m depth-zone sum",
  "2011 shallow Leptolyngbya maximum abundance", as.character(historical_2011_shallow$lepto_abundance_depth_sum_cells_l), "cells L^-1; sum across observed depths",
  "2011 shallow Leptolyngbya biovolume contribution", as.character(historical_2011_shallow$lepto_biovolume_percent), "% of total depth-zone phytoplankton biovolume on that date",
  "2011 deep Leptolyngbya maximum date", format(historical_2011_deep$date, "%Y-%m-%d"), "60-105 m depth-zone sum",
  "2011 deep Leptolyngbya maximum abundance", as.character(historical_2011_deep$lepto_abundance_depth_sum_cells_l), "cells L^-1; sum across observed depths",
  "2011 deep Leptolyngbya biovolume contribution", as.character(historical_2011_deep$lepto_biovolume_percent), "% of total depth-zone phytoplankton biovolume on that date",
  "October 2021 Schmidt stability", as.character(oct_stability$schmidt_stability_j_m2), "J m^-2",
  "December 2021 Schmidt stability", as.character(dec_stability$schmidt_stability_j_m2), "J m^-2",
  "October-December Schmidt stability decline", as.character(stability_decline_pct), "%",
  "Maximum post-spring Leptolyngbya abundance share", as.character(post_spring_max_abundance_pct$abundance_pct), paste0("%; ", post_spring_max_abundance_pct$date, ", ", post_spring_max_abundance_pct$depth, " m"),
  "Maximum post-spring Leptolyngbya biovolume share", as.character(post_spring_max_biovolume_pct$biovolume_pct), paste0("%; ", post_spring_max_biovolume_pct$date, ", ", post_spring_max_biovolume_pct$depth, " m"),
  "Maximum post-spring upper-zone Leptolyngbya abundance share", as.character(post_spring_upper_abundance$abundance_pct), paste0("%; ", post_spring_upper_abundance$date),
  "Maximum post-spring upper-zone Leptolyngbya biovolume share", as.character(post_spring_upper_biovolume$biovolume_pct), paste0("%; ", post_spring_upper_biovolume$date)
)

if (WRITE_EXISTING_OUTPUTS) {
  write_csv(stats, file.path(OUT_DIR, "requested_manuscript_statistics.csv"))
  write_csv(integrated_tests,
            file.path(OUT_DIR, "integrated_chlorophyll_aug2021_feb2022_tests.csv"))
}
write_csv(historical_2011_maxima,
          file.path(OUT_DIR, "leptolyngbya_2011_depth_zone_maxima.csv"))
write_csv(figure5b_focal_years,
          file.path(OUT_DIR, "figure_5b_focal_year_contrasts.csv"))
write_csv(figure5b_comparisons,
          file.path(OUT_DIR, "figure_5b_manuscript_statistics.csv"))

sentences <- c(
  sprintf("Deep 60-105 m chlorophyll, October 2021: %.2f SD; two-sided t test against the historical October mean, t(%d) = %.2f, Benjamini-Hochberg-adjusted p = %.3g across 12 Figure 4 month-by-depth comparisons.",
          deep_oct$sd_anomaly, deep_oct$clim_n - 1L, deep_oct$t_simple, deep_oct$q_fdr),
  sprintf("February 2022 integrated chlorophyll: %+.2f SD; two-sided t test against the historical February mean, t(%d) = %.2f, Benjamini-Hochberg-adjusted p = %.3g across seven monthly comparisons.",
          integrated_tests$anomaly_sd[integrated_tests$month == 2L],
          integrated_tests$hist_n[integrated_tests$month == 2L] - 1L,
          integrated_tests$t_simple[integrated_tests$month == 2L],
          integrated_tests$q_fdr[integrated_tests$month == 2L]),
  sprintf("Peak upper-water Leptolyngbya: %s cells L^-1 at %.0f m in %s; %.1f%% of sample biovolume and %.1f%% of upper-zone biovolume that date.",
          format(round(peak_upper$lepto_abundance), big.mark = ","), peak_upper$depth,
          format(peak_upper$date, "%B %Y"), peak_upper$biovolume_pct,
          peak_upper_date$biovolume_pct),
  sprintf(paste0("Historical 2011 Leptolyngbya depth-zone maxima: %s cells L^-1 on %s in 0-40 m ",
                 "(%.2f%% of community biovolume) and %s cells L^-1 on %s in 60-105 m ",
                 "(%.3f%% of community biovolume)."),
          format(round(historical_2011_shallow$lepto_abundance_depth_sum_cells_l), big.mark = ","),
          format(historical_2011_shallow$date, "%d %B %Y"),
          historical_2011_shallow$lepto_biovolume_percent,
          format(round(historical_2011_deep$lepto_abundance_depth_sum_cells_l), big.mark = ","),
          format(historical_2011_deep$date, "%d %B %Y"),
          historical_2011_deep$lepto_biovolume_percent),
  sprintf(paste0("Figure 5b contrast: August-December integrated PM2.5 was %.1f ug d m^-3 in 2011, ",
                 "%.1f in 2020, and %.1f in 2021; Leptolyngbya temporal AUC was %.3g, zero, ",
                 "and %.3g cells d m^-2, respectively, making the 2021 response %.2f-fold the 2011 response."),
          get_focal_value(2011L, integrated_pm25_aug_dec_ug_d_m3),
          get_focal_value(2020L, integrated_pm25_aug_dec_ug_d_m3),
          get_focal_value(2021L, integrated_pm25_aug_dec_ug_d_m3),
          get_focal_value(2011L, biological_temporal_auc_cells_d_m2),
          get_focal_value(2021L, biological_temporal_auc_cells_d_m2),
          get_focal_value(2021L, biological_temporal_auc_cells_d_m2) /
            get_focal_value(2011L, biological_temporal_auc_cells_d_m2)),
  sprintf("Schmidt stability declined from %s to %s J m^-2 between October and December 2021 (%.1f%% decline).",
          format(round(oct_stability$schmidt_stability_j_m2), big.mark = ","),
          format(round(dec_stability$schmidt_stability_j_m2), big.mark = ","),
          stability_decline_pct),
  sprintf("After spring 2022, the maximum sample-level Leptolyngbya share was %.1f%% of abundance and %.1f%% of biovolume; upper-zone maxima were %.1f%% and %.1f%%, respectively.",
          post_spring_max_abundance_pct$abundance_pct,
          post_spring_max_biovolume_pct$biovolume_pct,
          post_spring_upper_abundance$abundance_pct,
          post_spring_upper_biovolume$biovolume_pct)
)
if (WRITE_EXISTING_OUTPUTS) {
  writeLines(sentences, file.path(OUT_DIR, "requested_manuscript_sentences.txt"))
}
cat(paste(sentences, collapse = "\n"), "\n")
