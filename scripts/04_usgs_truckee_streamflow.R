# =============================================================================
# 04_usgs_truckee_streamflow.R
#
# Two-site Lake Tahoe outlet/inflow discharge and Upper Truckee water quality.
# Sites:
#   10337500 - Truckee River at Tahoe City (Lake Tahoe outlet)
#   10336610 - Upper Truckee River at South Lake Tahoe (major inflow)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(dataRetrieval)
  library(patchwork)
  library(scales)
  library(glue)
  library(fs)
  library(ragg)
})

# ---- Configuration ----------------------------------------------------------
SITES <- tribble(
  ~site_no,    ~site_name,
  "10337500", "Truckee River at Tahoe City",
  "10336610", "Upper Truckee River at South Lake Tahoe"
)
PARAM_CD <- "00060"
STAT_CD <- "00003"
CONTINUOUS_START <- as.POSIXct("2021-09-01 00:00:00", tz = "UTC")
CONTINUOUS_END <- as.POSIXct("2021-11-30 23:59:59", tz = "UTC")
EVENT_START <- as.Date("2021-10-21")
EVENT_END <- as.Date("2021-10-28")
FIRE_START <- as.Date("2021-08-14")
FIRE_END <- as.Date("2021-10-21")
WQ_START <- as.Date("2000-01-01")
WQ_END <- as.Date("2021-12-31")

PROC_DIR <- "data/processed"
FLOW_FIG_DIR <- "figures/streamflow"
WQ_FIG_DIR <- "figures/water_quality"
dir_create(PROC_DIR, recurse = TRUE)
dir_create(FLOW_FIG_DIR, recurse = TRUE)
dir_create(WQ_FIG_DIR, recurse = TRUE)

OUT <- c(
  tahoe_daily = file.path(PROC_DIR, "tahoe_city_discharge_daily.csv"),
  upper_daily = file.path(PROC_DIR, "upper_truckee_discharge_daily.csv"),
  continuous = file.path(PROC_DIR, "truckee_sites_discharge_2021.csv"),
  water_quality = file.path(PROC_DIR, "upper_truckee_water_quality.csv"),
  event_samples = file.path(PROC_DIR, "upper_truckee_event_samples.csv"),
  event_loads = file.path(PROC_DIR, "upper_truckee_event_loads.csv"),
  inventory = file.path(PROC_DIR, "upper_truckee_parameter_inventory.csv"),
  comparison = file.path(WQ_FIG_DIR, "upper_truckee_event_vs_baseline.png"),
  hydrograph = file.path(FLOW_FIG_DIR, "upper_truckee_first_runoff_2021.png"),
  site_comparison = file.path(FLOW_FIG_DIR, "tahoe_city_upper_truckee_2021.png")
)

cat("=== Two-site Truckee discharge and Upper Truckee water quality ===\n")
cat("dataRetrieval version:", as.character(packageVersion("dataRetrieval")), "\n")

# ---- Helpers ----------------------------------------------------------------
site_label <- function(site_no) SITES$site_name[match(site_no, SITES$site_no)]

download_daily <- function(site_no) {
  site_id <- paste0("USGS-", site_no)
  modern <- tryCatch(
    read_waterdata_daily(
      monitoring_location_id = site_id,
      parameter_code = PARAM_CD,
      statistic_id = STAT_CD,
      time = c("1850-01-01", as.character(Sys.Date()))
    ),
    error = function(e) {
      warning("Modern daily API failed for ", site_no, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(modern) && nrow(modern) > 0) {
    return(as_tibble(modern) %>%
      transmute(
        site_no = str_remove(.data[["monitoring_location_id"]], "^USGS-"),
        site_name = site_label(site_no),
        date = as.Date(.data[["time"]]),
        discharge_cfs = as.numeric(.data[["value"]]),
        discharge_cms = discharge_cfs * 0.0283168466,
        discharge_units = .data[["unit_of_measure"]],
        discharge_qualifier = .data[["qualifier"]],
        discharge_approval = .data[["approval_status"]],
        statistic_code = .data[["statistic_id"]],
        time_series_id = .data[["time_series_id"]],
        source_service = "USGS Water Data daily API"
      ))
  }

  legacy <- readNWISdv(
    siteNumbers = site_no,
    parameterCd = PARAM_CD,
    statCd = STAT_CD
  )
  value_col <- grep("00060_00003$", names(legacy), value = TRUE)[1]
  qualifier_col <- grep("00060_00003_cd$", names(legacy), value = TRUE)[1]
  as_tibble(legacy) %>%
    transmute(
      site_no = .data[["site_no"]],
      site_name = site_label(site_no),
      date = as.Date(.data[["Date"]]),
      discharge_cfs = as.numeric(.data[[value_col]]),
      discharge_cms = discharge_cfs * 0.0283168466,
      discharge_units = "ft^3/s",
      discharge_qualifier = .data[[qualifier_col]],
      discharge_approval = NA_character_,
      statistic_code = STAT_CD,
      time_series_id = NA_character_,
      source_service = "NWIS daily values fallback"
    )
}

download_continuous <- function(site_no) {
  site_id <- paste0("USGS-", site_no)
  modern <- tryCatch(
    read_waterdata_continuous(
      monitoring_location_id = site_id,
      parameter_code = PARAM_CD,
      time = c("2021-09-01T00:00:00Z", "2021-11-30T23:59:59Z")
    ),
    error = function(e) {
      warning("Modern continuous API failed for ", site_no, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(modern) && nrow(modern) > 0) {
    return(as_tibble(modern) %>%
      transmute(
        site_no = str_remove(.data[["monitoring_location_id"]], "^USGS-"),
        site_name = site_label(site_no),
        time_utc = as.POSIXct(.data[["time"]], tz = "UTC"),
        time_local = format(with_tz(time_utc, "America/Los_Angeles"),
                            "%Y-%m-%d %H:%M:%S %Z"),
        discharge_cfs = as.numeric(.data[["value"]]),
        discharge_cms = discharge_cfs * 0.0283168466,
        discharge_units = .data[["unit_of_measure"]],
        discharge_qualifier = .data[["qualifier"]],
        discharge_approval = .data[["approval_status"]],
        statistic_code = .data[["statistic_id"]],
        time_series_id = .data[["time_series_id"]],
        source_service = "USGS Water Data continuous API"
      ))
  }

  legacy <- readNWISuv(
    siteNumbers = site_no,
    parameterCd = PARAM_CD,
    startDate = "2021-09-01",
    endDate = "2021-11-30"
  )
  legacy <- renameNWISColumns(legacy)
  as_tibble(legacy) %>%
    transmute(
      site_no = .data[["site_no"]],
      site_name = site_label(site_no),
      time_utc = with_tz(.data[["dateTime"]], "UTC"),
      time_local = format(with_tz(time_utc, "America/Los_Angeles"),
                          "%Y-%m-%d %H:%M:%S %Z"),
      discharge_cfs = as.numeric(.data[["Flow_Inst"]]),
      discharge_cms = discharge_cfs * 0.0283168466,
      discharge_units = "ft^3/s",
      discharge_qualifier = .data[["Flow_Inst_cd"]],
      discharge_approval = NA_character_,
      statistic_code = "00011",
      time_series_id = NA_character_,
      source_service = "NWIS instantaneous values fallback"
    )
}

nearest_index <- function(target, reference) {
  reference_num <- as.numeric(reference)
  target_num <- as.numeric(target)
  ord <- order(reference_num)
  reference_num <- reference_num[ord]
  pos <- findInterval(target_num, reference_num)
  before <- pmax(pos, 1L)
  after <- pmin(pos + 1L, length(reference_num))
  use_after <- abs(reference_num[after] - target_num) <
    abs(reference_num[before] - target_num)
  ord[ifelse(use_after, after, before)]
}

classify_requested_variable <- function(characteristic, user_name) {
  txt <- str_to_lower(str_c(characteristic, user_name, sep = " | "))
  case_when(
    str_detect(txt, "suspended sediment|suspended solids") ~ "Suspended sediment",
    str_detect(txt, "turbid") ~ "Turbidity",
    str_detect(txt, "nitrate and nitrite|nitrate plus nitrite") ~ "Nitrate + nitrite",
    str_detect(txt, "total nitrogen|nitrogen, mixed forms") ~ "Total nitrogen",
    str_detect(txt, "kjeldahl|organic nitrogen|ammonia.*organic nitrogen") ~
      "Kjeldahl/organic nitrogen",
    str_detect(txt, "ammonia|ammonium") ~ "Ammonium/ammonia",
    str_detect(txt, "orthophosphate") ~ "Orthophosphate",
    str_detect(txt, "total phosphorus|phosphorus.*unfiltered") ~ "Total phosphorus",
    str_detect(txt, "dissolved phosphorus|phosphorus.*filtered") ~ "Dissolved phosphorus",
    str_detect(txt, "dissolved organic carbon|organic carbon.*filtered") ~
      "Dissolved organic carbon",
    TRUE ~ NA_character_
  )
}

# ---- Discharge downloads ----------------------------------------------------
daily_list <- map(SITES$site_no, download_daily)
names(daily_list) <- SITES$site_no
tahoe_daily <- daily_list[["10337500"]]
upper_daily <- daily_list[["10336610"]]
write_csv(tahoe_daily, OUT[["tahoe_daily"]], na = "")
write_csv(upper_daily, OUT[["upper_daily"]], na = "")

continuous <- map_dfr(SITES$site_no, download_continuous) %>%
  arrange(site_no, time_utc)
write_csv(continuous, OUT[["continuous"]], na = "")
upper_continuous <- continuous %>% filter(site_no == "10336610")

# ---- Detect the first major post-Caldor runoff rise and peak ----------------
upper_hourly <- upper_continuous %>%
  mutate(hour_utc = floor_date(time_utc, "hour")) %>%
  group_by(hour_utc) %>%
  summarise(discharge_cfs = mean(discharge_cfs, na.rm = TRUE), .groups = "drop")

event_hourly <- upper_hourly %>%
  filter(as.Date(hour_utc) >= EVENT_START, as.Date(hour_utc) <= EVENT_END)
pre_event_baseline <- upper_hourly %>%
  filter(hour_utc >= as.POSIXct("2021-10-21 00:00:00", tz = "UTC"),
         hour_utc < as.POSIXct("2021-10-21 12:00:00", tz = "UTC")) %>%
  summarise(value = median(discharge_cfs, na.rm = TRUE)) %>%
  pull(value)
event_peak <- event_hourly %>% slice_max(discharge_cfs, n = 1, with_ties = FALSE)
rise_threshold <- pre_event_baseline + max(
  2,
  0.10 * (event_peak$discharge_cfs - pre_event_baseline)
)
event_rise <- event_hourly %>%
  filter(discharge_cfs >= rise_threshold) %>%
  slice_head(n = 1)

# ---- Inspect inventory, then download all 2000-2021 sample results ----------
UPPER_SITE_ID <- "USGS-10336610"
sample_summary <- as_tibble(
  summarize_waterdata_samples(monitoringLocationIdentifier = UPPER_SITE_ID)
) %>%
  mutate(
    requested_variable = classify_requested_variable(
      characteristic, characteristicUserSupplied
    ),
    requested = !is.na(requested_variable)
  )

cat("\nUSGS sample inventory rows:", nrow(sample_summary), "\n")
cat("Inventory rows matching requested constituent families:",
    sum(sample_summary$requested), "\n")
print(
  sample_summary %>%
    filter(requested) %>%
    distinct(requested_variable, characteristic, characteristicUserSupplied,
             firstActivity, mostRecentActivity) %>%
    arrange(requested_variable, characteristicUserSupplied),
  n = Inf
)

wq_raw <- as_tibble(
  read_waterdata_samples(
    monitoringLocationIdentifier = UPPER_SITE_ID,
    activityStartDateLower = as.character(WQ_START),
    activityStartDateUpper = as.character(WQ_END),
    dataProfile = "fullphyschem"
  )
)

wq <- wq_raw %>%
  transmute(
    site_id = .data[["Location_Identifier"]],
    site_name = .data[["Location_Name"]],
    sample_id = .data[["Activity_ActivityIdentifier"]],
    sample_id_user = .data[["Activity_ActivityIdentifierUserSupplied"]],
    sample_date = as.Date(.data[["Activity_StartDate"]]),
    sample_datetime_utc = as.POSIXct(.data[["Activity_StartDateTime"]], tz = "UTC"),
    activity_type = .data[["Activity_TypeCode"]],
    activity_media = .data[["Activity_Media"]],
    hydrologic_condition = .data[["Activity_HydrologicCondition"]],
    hydrologic_event = .data[["Activity_HydrologicEvent"]],
    characteristic = .data[["Result_Characteristic"]],
    characteristic_user = .data[["Result_CharacteristicUserSupplied"]],
    characteristic_group = .data[["Result_CharacteristicGroup"]],
    requested_variable = classify_requested_variable(characteristic, characteristic_user),
    usgs_parameter_code = .data[["USGSpcode"]],
    sample_fraction = .data[["Result_SampleFraction"]],
    result_value = suppressWarnings(as.numeric(.data[["Result_Measure"]])),
    result_units = .data[["Result_MeasureUnit"]],
    result_qualifier = .data[["Result_MeasureQualifierCode"]],
    result_status = .data[["Result_MeasureStatusIdentifier"]],
    result_measure_id = .data[["Result_MeasureIdentifier"]],
    detection_condition = .data[["Result_ResultDetectionCondition"]],
    detection_limit_type_a = .data[["DetectionLimit_TypeA"]],
    detection_limit_a = suppressWarnings(as.numeric(.data[["DetectionLimit_MeasureA"]])),
    detection_limit_units_a = .data[["DetectionLimit_MeasureUnitA"]],
    detection_limit_comment_a = .data[["DetectionLimit_CommentA"]],
    detection_limit_type_b = .data[["DetectionLimit_TypeB"]],
    detection_limit_b = suppressWarnings(as.numeric(.data[["DetectionLimit_MeasureB"]])),
    detection_limit_units_b = .data[["DetectionLimit_MeasureUnitB"]],
    detection_limit_comment_b = .data[["DetectionLimit_CommentB"]],
    analytical_method_id = .data[["ResultAnalyticalMethod_Identifier"]],
    analytical_method_context = .data[["ResultAnalyticalMethod_IdentifierContext"]],
    analytical_method_name = .data[["ResultAnalyticalMethod_Name"]],
    analytical_method_description = .data[["ResultAnalyticalMethod_Description"]],
    laboratory = .data[["LabInfo_Name"]],
    result_comment = .data[["DataQuality_ResultComment"]],
    activity_comment = .data[["Activity_Comment"]]
  ) %>%
  filter(!is.na(requested_variable)) %>%
  mutate(
    sample_datetime_utc = coalesce(
      sample_datetime_utc,
      as.POSIXct(sample_date, tz = "UTC") + hours(12)
    ),
    nondetect_flag = str_detect(
      str_to_lower(coalesce(detection_condition, "")),
      "not detected|non-detect|below|less than"
    ) | str_detect(coalesce(result_qualifier, ""), "<")
  ) %>%
  mutate(row_id = row_number())

# ---- Match every sample to nearest available discharge ----------------------
daily_reference_time <- as.POSIXct(upper_daily$date, tz = "UTC") + hours(12)
daily_idx <- nearest_index(wq$sample_datetime_utc, daily_reference_time)

wq <- wq %>%
  mutate(
    daily_discharge_date = upper_daily$date[daily_idx],
    daily_discharge_cfs = upper_daily$discharge_cfs[daily_idx],
    daily_discharge_qualifier = upper_daily$discharge_qualifier[daily_idx],
    daily_discharge_approval = upper_daily$discharge_approval[daily_idx],
    daily_match_difference_hours = abs(
      as.numeric(difftime(sample_datetime_utc, daily_reference_time[daily_idx],
                          units = "hours"))
    )
  )

continuous_candidates <- wq %>%
  filter(sample_datetime_utc >= min(upper_continuous$time_utc),
         sample_datetime_utc <= max(upper_continuous$time_utc))
if (nrow(continuous_candidates) > 0) {
  continuous_idx <- nearest_index(
    continuous_candidates$sample_datetime_utc,
    upper_continuous$time_utc
  )
  continuous_matches <- continuous_candidates %>%
    transmute(
      row_id,
      continuous_discharge_time_utc = upper_continuous$time_utc[continuous_idx],
      continuous_discharge_cfs = upper_continuous$discharge_cfs[continuous_idx],
      continuous_discharge_qualifier = upper_continuous$discharge_qualifier[continuous_idx],
      continuous_discharge_approval = upper_continuous$discharge_approval[continuous_idx],
      continuous_match_difference_minutes = abs(
        as.numeric(difftime(
          sample_datetime_utc,
          upper_continuous$time_utc[continuous_idx],
          units = "mins"
        ))
      )
    )
  wq <- wq %>% left_join(continuous_matches, by = "row_id")
} else {
  wq <- wq %>%
    mutate(
      continuous_discharge_time_utc = as.POSIXct(NA, tz = "UTC"),
      continuous_discharge_cfs = NA_real_,
      continuous_discharge_qualifier = NA_character_,
      continuous_discharge_approval = NA_character_,
      continuous_match_difference_minutes = NA_real_
    )
}

wq <- wq %>%
  mutate(
    matched_discharge_source = if_else(
      !is.na(continuous_discharge_cfs),
      "Nearest 15-minute observation",
      "Nearest daily mean"
    ),
    matched_discharge_time = coalesce(
      continuous_discharge_time_utc,
      as.POSIXct(daily_discharge_date, tz = "UTC") + hours(12)
    ),
    matched_discharge_cfs = coalesce(continuous_discharge_cfs, daily_discharge_cfs),
    matched_discharge_qualifier = coalesce(
      continuous_discharge_qualifier,
      daily_discharge_qualifier
    ),
    matched_discharge_approval = coalesce(
      continuous_discharge_approval,
      daily_discharge_approval
    )
  ) %>%
  select(-row_id)
write_csv(wq, OUT[["water_quality"]], na = "")

# ---- Parameter inventory and availability ----------------------------------
requested_levels <- c(
  "Suspended sediment", "Turbidity", "Nitrate + nitrite",
  "Ammonium/ammonia", "Total nitrogen", "Kjeldahl/organic nitrogen",
  "Orthophosphate", "Total phosphorus", "Dissolved phosphorus",
  "Dissolved organic carbon"
)

parameter_inventory <- wq %>%
  group_by(
    requested_variable, characteristic, characteristic_user,
    usgs_parameter_code, sample_fraction, result_units
  ) %>%
  summarise(
    result_count_2000_2021 = n(),
    numeric_result_count = sum(!is.na(result_value)),
    nondetect_count = sum(nondetect_flag, na.rm = TRUE),
    first_sample = min(sample_date),
    last_sample = max(sample_date),
    event_window_results = sum(sample_date >= EVENT_START & sample_date <= EVENT_END),
    prefire_oct_nov_results = sum(
      year(sample_date) < 2021 & month(sample_date) %in% c(10, 11)
    ),
    .groups = "drop"
  ) %>%
  mutate(available = TRUE)

missing_requested <- tibble(requested_variable = requested_levels) %>%
  anti_join(parameter_inventory %>% distinct(requested_variable),
            by = "requested_variable") %>%
  mutate(
    characteristic = NA_character_, characteristic_user = NA_character_,
    usgs_parameter_code = NA_character_, sample_fraction = NA_character_,
    result_units = NA_character_, result_count_2000_2021 = 0L,
    numeric_result_count = 0L, nondetect_count = 0L,
    first_sample = as.Date(NA), last_sample = as.Date(NA),
    event_window_results = 0L, prefire_oct_nov_results = 0L,
    available = FALSE
  )
parameter_inventory <- bind_rows(parameter_inventory, missing_requested) %>%
  arrange(factor(requested_variable, levels = requested_levels), usgs_parameter_code)
write_csv(parameter_inventory, OUT[["inventory"]], na = "")

# ---- Event samples and sample-date load estimates ---------------------------
event_samples <- wq %>%
  filter(sample_date >= EVENT_START, sample_date <= EVENT_END) %>%
  arrange(sample_datetime_utc, requested_variable, usgs_parameter_code)
write_csv(event_samples, OUT[["event_samples"]], na = "")

event_loads <- event_samples %>%
  filter(str_to_lower(result_units) == "mg/l") %>%
  mutate(
    load_kg_day = if_else(
      !nondetect_flag & is.finite(result_value) & is.finite(matched_discharge_cfs),
      result_value * matched_discharge_cfs * 2.44658,
      NA_real_
    ),
    load_type = "Sample-date load estimate (not event-integrated)",
    load_equation = "concentration_mg_l * discharge_cfs * 2.44658"
  ) %>%
  rename(concentration_mg_l = result_value)
write_csv(event_loads, OUT[["event_loads"]], na = "")

# ---- Canonical comparable constituents for event-vs-baseline figure --------
canonical <- tribble(
  ~requested_variable, ~usgs_parameter_code, ~sample_fraction, ~result_units, ~plot_label,
  "Suspended sediment", "80154", "Unfiltered", "mg/L", "Suspended sediment",
  "Turbidity", "63675", "Unfiltered", "NTU", "Turbidity",
  "Nitrate + nitrite", "00631", "Filtered field and/or lab", "mg/L", "NO₃ + NO₂ as N",
  "Ammonium/ammonia", "00608", "Filtered field and/or lab", "mg/L", "NH₄/NH₃ as N",
  "Total nitrogen", "00600", "Unfiltered", "mg/L", "Total nitrogen",
  "Kjeldahl/organic nitrogen", "00625", "Unfiltered", "mg/L", "TKN/organic N",
  "Orthophosphate", "00671", "Filtered field and/or lab", "mg/L", "Orthophosphate as P",
  "Total phosphorus", "00665", "Unfiltered", "mg/L", "Total phosphorus",
  "Dissolved phosphorus", "00666", "Filtered field and/or lab", "mg/L", "Dissolved phosphorus",
  "Dissolved organic carbon", "00681", "Filtered field and/or lab", "mg/L", "Dissolved organic carbon"
)

comparison_data <- wq %>%
  inner_join(
    canonical,
    by = c("requested_variable", "usgs_parameter_code", "sample_fraction", "result_units")
  ) %>%
  mutate(
    period = case_when(
      sample_date >= EVENT_START & sample_date <= EVENT_END ~ "Oct 2021 runoff",
      year(sample_date) < 2021 & month(sample_date) %in% c(10, 11) ~
        "Historical Oct-Nov",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Historical Oct-Nov", "Oct 2021 runoff"))
  ) %>%
  filter(!is.na(period), !nondetect_flag, is.finite(result_value))

usable_baseline <- comparison_data %>%
  group_by(requested_variable, plot_label, usgs_parameter_code,
           sample_fraction, result_units) %>%
  summarise(
    n_baseline = sum(period == "Historical Oct-Nov"),
    n_event = sum(period == "Oct 2021 runoff"),
    .groups = "drop"
  ) %>%
  mutate(usable = n_baseline >= 3 & n_event >= 1)

comparison_plot_data <- comparison_data %>%
  semi_join(usable_baseline %>% filter(usable),
            by = c("requested_variable", "plot_label", "usgs_parameter_code",
                   "sample_fraction", "result_units")) %>%
  mutate(
    facet_label = paste0(plot_label, "\n", usgs_parameter_code, "; ",
                         sample_fraction, "; ", result_units)
  )

# ---- Figures ----------------------------------------------------------------
flow_2021 <- bind_rows(tahoe_daily, upper_daily) %>%
  filter(year(date) == 2021)

p_sites <- ggplot(flow_2021, aes(date, discharge_cfs)) +
  annotate("rect", xmin = FIRE_START, xmax = FIRE_END,
           ymin = -Inf, ymax = Inf, fill = "#F18A00", alpha = 0.10) +
  geom_line(colour = "#1A3A6B", linewidth = 0.45) +
  facet_wrap(vars(site_name), ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = expression("Daily mean discharge (ft"^3*" s"^-1*")")) +
  theme_classic(base_size = 8) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 7)
  )
ggsave(OUT[["site_comparison"]], p_sites, width = 18, height = 14,
       units = "cm", dpi = 300, device = ragg::agg_png)

sample_marks <- event_samples %>%
  distinct(sample_datetime_utc) %>%
  mutate(idx = nearest_index(sample_datetime_utc, upper_continuous$time_utc),
         discharge_cfs = upper_continuous$discharge_cfs[idx])

p_event <- ggplot(upper_continuous, aes(time_utc, discharge_cfs)) +
  annotate("rect",
           xmin = as.POSIXct(EVENT_START, tz = "UTC"),
           xmax = as.POSIXct(EVENT_END + 1, tz = "UTC"),
           ymin = -Inf, ymax = Inf, fill = "#F18A00", alpha = 0.10) +
  geom_line(colour = "#1A3A6B", linewidth = 0.45) +
  geom_vline(data = sample_marks, aes(xintercept = sample_datetime_utc),
             colour = "#009E73", linetype = "dotted", linewidth = 0.45) +
  geom_point(data = sample_marks, aes(sample_datetime_utc, discharge_cfs),
             inherit.aes = FALSE, colour = "#009E73", size = 2) +
  geom_point(data = event_rise, aes(hour_utc, discharge_cfs),
             inherit.aes = FALSE, colour = "#D55E00", shape = 17, size = 2.4) +
  geom_point(data = event_peak, aes(hour_utc, discharge_cfs),
             inherit.aes = FALSE, colour = "#B71C1C", shape = 17, size = 2.4) +
  annotate("text", x = event_rise$hour_utc, y = event_rise$discharge_cfs,
           label = "Detected rise", hjust = 1.05, vjust = -0.8,
           size = 2.5, colour = "#D55E00") +
  annotate("text", x = event_peak$hour_utc, y = event_peak$discharge_cfs,
           label = "Peak", hjust = -0.1, vjust = -0.8,
           size = 2.5, colour = "#B71C1C") +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b %d",
                   expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL,
    y = expression("Discharge (ft"^3*" s"^-1*")"),
    caption = "Green markers are discrete water-quality sampling times; orange shading is Oct 21-28."
  ) +
  theme_classic(base_size = 8) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25))
ggsave(OUT[["hydrograph"]], p_event, width = 18, height = 10,
       units = "cm", dpi = 300, device = ragg::agg_png)

if (nrow(comparison_plot_data) > 0) {
  p_wq <- ggplot(comparison_plot_data, aes(period, result_value, colour = period)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
    geom_jitter(width = 0.10, height = 0, size = 1.25, alpha = 0.70) +
    facet_wrap(vars(facet_label), scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(
      "Historical Oct-Nov" = "grey45",
      "Oct 2021 runoff" = "#D55E00"
    )) +
    labs(x = NULL, y = "Reported concentration or turbidity", colour = NULL,
         caption = "Only identical USGS parameter code, sample fraction, and unit combinations are compared.") +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 6.5),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25)
    )
} else {
  p_wq <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "No requested variables had both event observations and a usable historical baseline") +
    theme_void()
}
ggsave(OUT[["comparison"]], p_wq, width = 20, height = 15,
       units = "cm", dpi = 300, device = ragg::agg_png)

# ---- Final console summary --------------------------------------------------
cat("\n=== FINAL SUMMARY ===\n")
for (site in SITES$site_no) {
  dat <- daily_list[[site]]
  cat(site, "-", site_label(site), ":",
      as.character(min(dat$date, na.rm = TRUE)), "to",
      as.character(max(dat$date, na.rm = TRUE)), "(", nrow(dat), "days )\n")
}
cat("Detected runoff rise:", format(event_rise$hour_utc, tz = "UTC"),
    sprintf("(%.1f cfs; threshold %.1f cfs)\n", event_rise$discharge_cfs, rise_threshold))
cat("Detected runoff peak:", format(event_peak$hour_utc, tz = "UTC"),
    sprintf("(%.1f cfs)\n", event_peak$discharge_cfs))

cat("Upper Truckee samples during Oct 21-28 runoff:\n")
print(event_samples %>%
        distinct(sample_date, sample_datetime_utc) %>%
        arrange(sample_datetime_utc), n = Inf)

cat("Requested nutrient/sediment families available:\n")
print(parameter_inventory %>%
        filter(available) %>%
        distinct(requested_variable) %>%
        arrange(requested_variable), n = Inf)

cat("Requested families not found:\n")
print(parameter_inventory %>%
        filter(!available) %>%
        select(requested_variable), n = Inf)

cat("Variables with usable pre-fire Oct-Nov baselines and event data:\n")
print(usable_baseline %>% filter(usable), n = Inf)

cat("Output files:\n")
cat(paste0("  ", unname(OUT)), sep = "\n")
cat("\nSample-date loads use concentration_mg_l * discharge_cfs * 2.44658; ",
    "they are not event-integrated loads.\n", sep = "")
