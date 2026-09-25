# =============================================================================
# download_tahoe_epa_pm_daily.R
#
# Purpose : Download daily monitor-level PM data from the EPA AQS API
#           (dailyData/byBox endpoint) for the Lake Tahoe airshed,
#           parse all responses, and write analysis-ready outputs.
#
# Params  : 88101 = PM2.5 FRM/FEM Mass
#           88502 = PM2.5 non-FRM/FEM Mass
#           81102 = PM10 Mass
#
# Requires: AQS_EMAIL and AQS_KEY environment variables
#
# Usage   : Run from the project root.
#           Rscript epa_aqs_tahoe/R/download_tahoe_epa_pm_daily.R
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(httr2)
  library(jsonlite)
  library(lubridate)
  library(fs)
  library(glue)
  library(readr)
  library(arrow)
  library(purrr)
  library(stringr)
})

# =============================================================================
# SECTION 0: Configuration
# =============================================================================

# -- Credentials --------------------------------------------------------------
AQS_EMAIL <- Sys.getenv("AQS_EMAIL")
AQS_KEY   <- Sys.getenv("AQS_KEY")

if (nchar(AQS_EMAIL) == 0 || nchar(AQS_KEY) == 0) {
  stop(
    "\033[31m",
    "AQS credentials not found in environment variables.\n",
    "Please set AQS_EMAIL and AQS_KEY before running.\n",
    "See README.md for instructions on obtaining an API key.\033[0m"
  )
}

# -- Bounding box: Lake Tahoe region ------------------------------------------
BBOX <- list(
  minlat =  38.80,
  maxlat =  39.35,
  minlon = -120.25,
  maxlon = -119.80
)

# -- Parameters ----------------------------------------------------------------
PARAMS <- c(
  "88101" = "PM2.5 FRM/FEM",
  "88502" = "PM2.5 non-FRM/FEM",
  "81102" = "PM10"
)

# -- Year range ----------------------------------------------------------------
YEAR_START <- 1990L
YEAR_END   <- year(today()) - 1L   # last fully completed calendar year

# -- Optional wildfire-season date filter --------------------------------------
# Set SEASON_FILTER to NULL to retain all dates.
# To restrict to wildfire season (Jul 1 – Oct 31), uncomment the second line.
SEASON_FILTER <- NULL
# SEASON_FILTER <- list(start_md = "07-01", end_md = "10-31")

# -- Rate limit (EPA asks for at least 5 s between requests) ------------------
RATE_LIMIT_SEC <- 5L

# -- API endpoint --------------------------------------------------------------
BASE_URL <- "https://aqs.epa.gov/data/api/dailyData/byBox"

# -- Project paths (relative to working directory = project root) -------------
RAW_DIR  <- "data/raw/epa_aqs_tahoe_pm/json"
PROC_DIR <- "data/processed"
FIG_DIR  <- "figures"

fs::dir_create(RAW_DIR,  recurse = TRUE)
fs::dir_create(PROC_DIR, recurse = TRUE)
fs::dir_create(FIG_DIR,  recurse = TRUE)

# -- Columns to retain from API response --------------------------------------
KEEP_COLS <- c(
  "state_code", "county_code", "site_number", "poc",
  "latitude", "longitude", "datum",
  "parameter_code", "parameter",
  "date_local",
  "arithmetic_mean", "first_max_value", "aqi",
  "units_of_measure", "sample_duration",
  "observation_count", "percent_complete",
  "method_code", "method",
  "local_site_name", "site_address",
  "county", "city", "state"
)

# =============================================================================
# SECTION 1: Download raw JSON responses
# =============================================================================

cat("\n=== EPA AQS Daily PM Download: Lake Tahoe ===\n")
cat("Year range   :", YEAR_START, "–", YEAR_END, "\n")
cat("Parameters   :", paste(names(PARAMS), PARAMS, sep = " = ", collapse = "; "), "\n")
cat("Bounding box : lat [", BBOX$minlat, ",", BBOX$maxlat,
    "] lon [", BBOX$minlon, ",", BBOX$maxlon, "]\n")
cat("Season filter:", if (is.null(SEASON_FILTER)) "none (all dates)" else
      glue("{SEASON_FILTER$start_md} – {SEASON_FILTER$end_md}"), "\n\n")

years    <- seq(YEAR_START, YEAR_END)
p_codes  <- names(PARAMS)
n_total  <- length(years) * length(p_codes)
n_done   <- 0L

# Initialise API call log
api_log <- tibble(
  p_code    = character(),
  year      = integer(),
  status    = character(),
  n_rows    = integer(),
  json_file = character()
)

for (pc in p_codes) {
  for (yr in years) {
    n_done    <- n_done + 1L
    json_file <- file.path(RAW_DIR, glue("aqs_{pc}_{yr}.json"))

    # -- Cache: skip if already saved -----------------------------------------
    if (fs::file_exists(json_file)) {
      cat(glue("[{n_done}/{n_total}] CACHED  param={pc}  year={yr}\n"))
      api_log <- add_row(api_log,
        p_code = pc, year = yr, status = "cached",
        n_rows = NA_integer_, json_file = json_file)
      next
    }

    # -- Build URL -------------------------------------------------------------
    url <- glue(
      "{BASE_URL}",
      "?email={AQS_EMAIL}&key={AQS_KEY}",
      "&param={pc}",
      "&bdate={yr}0101&edate={yr}1231",
      "&minlat={BBOX$minlat}&maxlat={BBOX$maxlat}",
      "&minlon={BBOX$minlon}&maxlon={BBOX$maxlon}"
    )

    cat(glue("[{n_done}/{n_total}] Fetching param={pc}  year={yr} ... "))

    # -- HTTP request with retry -----------------------------------------------
    resp <- tryCatch(
      request(url) |>
        req_timeout(90) |>
        req_retry(max_tries = 3, backoff = ~ 20) |>
        req_perform(),
      error = function(e) {
        cat("\033[31mHTTP ERROR:", conditionMessage(e), "\033[0m\n")
        NULL
      }
    )

    if (is.null(resp)) {
      api_log <- add_row(api_log,
        p_code = pc, year = yr, status = "http_error",
        n_rows = NA_integer_, json_file = NA_character_)
      Sys.sleep(RATE_LIMIT_SEC)
      next
    }

    body   <- tryCatch(resp_body_string(resp), error = function(e) NULL)
    parsed <- tryCatch(fromJSON(body, simplifyVector = TRUE), error = function(e) NULL)

    if (is.null(parsed)) {
      cat("\033[31mPARSE ERROR\033[0m\n")
      api_log <- add_row(api_log,
        p_code = pc, year = yr, status = "parse_error",
        n_rows = NA_integer_, json_file = NA_character_)
      Sys.sleep(RATE_LIMIT_SEC)
      next
    }

    api_status <- tryCatch(parsed$Header$status[1], error = function(e) "unknown")
    n_rows     <- if (!is.null(parsed$Data) && is.data.frame(parsed$Data))
                    nrow(parsed$Data) else 0L

    cat(glue("status={api_status}  rows={n_rows}\n"))

    # Save raw JSON (always, even for no-data responses)
    writeLines(body, json_file)

    api_log <- add_row(api_log,
      p_code = pc, year = yr, status = api_status,
      n_rows = as.integer(n_rows), json_file = json_file)

    Sys.sleep(RATE_LIMIT_SEC)
  }
}

# =============================================================================
# SECTION 2: Parse all JSON files into a combined tibble
# =============================================================================

cat("\n--- Parsing JSON responses ---\n")

parse_one_json <- function(json_path) {
  # Extract param code and year from filename: aqs_{param}_{year}.json
  fname  <- fs::path_file(json_path)
  parts  <- str_match(fname, "aqs_(\\d+)_(\\d{4})\\.json")
  pc     <- parts[1, 2]
  yr     <- as.integer(parts[1, 3])

  body   <- paste(readLines(json_path, warn = FALSE), collapse = "\n")
  parsed <- tryCatch(fromJSON(body, simplifyVector = TRUE), error = function(e) NULL)

  if (is.null(parsed) ||
      is.null(parsed$Data) ||
      !is.data.frame(parsed$Data) ||
      nrow(parsed$Data) == 0)
    return(NULL)

  df <- as_tibble(parsed$Data)

  # Ensure all KEEP_COLS exist; fill any absent ones with NA
  for (col in setdiff(KEEP_COLS, names(df))) df[[col]] <- NA_character_
  df <- df[, KEEP_COLS]

  # Type coercions and derived metadata columns
  df %>%
    mutate(
      parameter_code    = pc,
      query_year        = yr,
      pollutant_group   = PARAMS[pc],
      tahoe_bbox        = TRUE,
      downloaded_at     = Sys.time(),
      date_local        = as.Date(date_local),
      arithmetic_mean   = suppressWarnings(as.numeric(arithmetic_mean)),
      first_max_value   = suppressWarnings(as.numeric(first_max_value)),
      aqi               = suppressWarnings(as.integer(aqi)),
      latitude          = suppressWarnings(as.numeric(latitude)),
      longitude         = suppressWarnings(as.numeric(longitude)),
      observation_count = suppressWarnings(as.integer(observation_count)),
      percent_complete  = suppressWarnings(as.numeric(percent_complete))
    )
}

# Gather all JSON files that exist and had data
json_paths <- api_log %>%
  filter(!is.na(json_file), fs::file_exists(json_file)) %>%
  pull(json_file) %>%
  unique()

cat("  Parsing", length(json_paths), "JSON files...\n")

pm_all <- map(json_paths, parse_one_json, .progress = TRUE) %>%
  compact() %>%
  bind_rows()

cat("  Combined table:", nrow(pm_all), "rows ×", ncol(pm_all), "columns\n")

if (nrow(pm_all) == 0) {
  stop("\033[31mNo data was successfully parsed. Check credentials and API log.\033[0m")
}

# -- Optional seasonal filter -------------------------------------------------
if (!is.null(SEASON_FILTER)) {
  pm_all <- pm_all %>%
    filter(
      format(date_local, "%m-%d") >= SEASON_FILTER$start_md,
      format(date_local, "%m-%d") <= SEASON_FILTER$end_md
    )
  cat(glue("  After season filter ({SEASON_FILTER$start_md} to {SEASON_FILTER$end_md}): {nrow(pm_all)} rows\n"))
}

# =============================================================================
# SECTION 3: QA checks
# =============================================================================

cat("\n--- QA Checks ---\n")

# Records by year and pollutant
cat("\nRecords per year per pollutant:\n")
records_by_year <- pm_all %>%
  count(query_year, pollutant_group) %>%
  pivot_wider(names_from = pollutant_group, values_from = n, values_fill = 0L)
print(records_by_year, n = Inf)

# Unique monitors by pollutant
cat("\nUnique monitors by pollutant:\n")
monitor_counts <- pm_all %>%
  mutate(monitor_id = paste(state_code, county_code, site_number, poc, sep = "-")) %>%
  group_by(pollutant_group) %>%
  summarise(n_monitors = n_distinct(monitor_id), .groups = "drop")
print(monitor_counts)

# Missing years per pollutant
cat("\nMissing years (no data returned):\n")
for (pc in p_codes) {
  years_with_data <- pm_all %>%
    filter(parameter_code == pc) %>%
    pull(query_year) %>%
    unique()
  missing <- setdiff(seq(YEAR_START, YEAR_END), years_with_data)
  if (length(missing) == 0) {
    cat(" ", PARAMS[pc], ": none\n")
  } else {
    cat("\033[33m ", PARAMS[pc], ": missing", length(missing), "years (",
        paste(range(missing), collapse = "–"), ")\033[0m\n")
  }
}

# API failures
n_failed <- api_log %>%
  filter(!status %in% c("Success", "No data matched your selection.", "cached")) %>%
  nrow()
if (n_failed > 0) {
  cat("\033[33mWARNING:", n_failed, "API calls failed (not Success/NoData/cached).\033[0m\n")
  print(filter(api_log, !status %in% c("Success", "No data matched your selection.", "cached")))
} else {
  cat("  All API calls succeeded or returned no-data.\n")
}

# PM10 specific check
pm10_rows <- filter(pm_all, parameter_code == "81102")
if (nrow(pm10_rows) == 0) {
  cat("\033[33mWARNING: No PM10 (81102) data found in the Tahoe bounding box.\033[0m\n")
} else {
  cat("  PM10 records found:", nrow(pm10_rows), "\n")
}

# =============================================================================
# SECTION 4: Write outputs
# =============================================================================

cat("\n--- Writing outputs ---\n")

# 4a. Full data — Parquet (efficient, preserves types)
write_parquet(pm_all,
              file.path(PROC_DIR, "tahoe_epa_pm_daily_all_years.parquet"))
cat("  Saved: tahoe_epa_pm_daily_all_years.parquet\n")

# 4b. Full data — CSV
write_csv(pm_all,
          file.path(PROC_DIR, "tahoe_epa_pm_daily_all_years.csv"),
          na = "")
cat("  Saved: tahoe_epa_pm_daily_all_years.csv\n")

# 4c. Inventory: records per year × pollutant × county
inventory <- pm_all %>%
  mutate(monitor_id = paste(state_code, county_code, site_number, poc, sep = "-")) %>%
  group_by(parameter_code, pollutant_group, query_year, county, state) %>%
  summarise(
    n_records        = n(),
    n_monitors       = n_distinct(monitor_id),
    n_days           = n_distinct(date_local),
    mean_pm          = round(mean(arithmetic_mean, na.rm = TRUE), 3),
    max_pm           = round(max(first_max_value,  na.rm = TRUE), 3),
    pct_complete_avg = round(mean(percent_complete, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(parameter_code, query_year, county)

write_csv(inventory,
          file.path(PROC_DIR, "tahoe_epa_pm_daily_inventory.csv"),
          na = "")
cat("  Saved: tahoe_epa_pm_daily_inventory.csv\n")

# 4d. Monitor locations: unique monitors with site metadata
monitor_locs <- pm_all %>%
  mutate(monitor_id = paste(state_code, county_code, site_number, poc, sep = "-")) %>%
  group_by(monitor_id, parameter_code, pollutant_group,
           state_code, county_code, site_number, poc,
           latitude, longitude, datum,
           local_site_name, site_address, county, city, state) %>%
  summarise(
    first_date = min(date_local, na.rm = TRUE),
    last_date  = max(date_local, na.rm = TRUE),
    n_obs      = n(),
    .groups = "drop"
  ) %>%
  arrange(parameter_code, state, county, monitor_id)

write_csv(monitor_locs,
          file.path(PROC_DIR, "tahoe_epa_pm_monitor_locations.csv"),
          na = "")
cat("  Saved: tahoe_epa_pm_monitor_locations.csv\n")

# =============================================================================
# SECTION 5: Session summary
# =============================================================================

cat("\n=== Summary ===\n")
cat("Total records downloaded  :", nrow(pm_all), "\n")
cat("Date range in data        :",
    as.character(min(pm_all$date_local)), "–",
    as.character(max(pm_all$date_local)), "\n")
cat("Unique monitors           :",
    n_distinct(paste(pm_all$state_code, pm_all$county_code,
                     pm_all$site_number, pm_all$poc)), "\n")
cat("Outputs written to        :", PROC_DIR, "\n")
cat("Figures written to        :", FIG_DIR, "\n")
cat("API log (failed calls)    :", n_failed, "\n")

# Write API log for reference
write_csv(api_log,
          file.path(PROC_DIR, "tahoe_epa_pm_api_log.csv"),
          na = "")
cat("API log saved to          : data/processed/tahoe_epa_pm_api_log.csv\n")
cat("\nDone.\n")
