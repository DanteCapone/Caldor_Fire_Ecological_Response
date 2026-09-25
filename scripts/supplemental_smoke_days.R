# =============================================================================
# supplemental_smoke_days.R
#
# Recreate the State of the Lake annual smoke-day plot for the supplement.
# Bars are days with NOAA HMS smoke over Lake Tahoe and daily composite PM2.5
# in either of the two upper concentration categories (>=9.1 ug/m3). The
# secondary-axis line is the maximum daily composite PM2.5 in each calendar
# year (all monitored days, not only HMS smoke days).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ragg)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
} else {
  project_root <- normalizePath(".", winslash = "/")
}
setwd(project_root)
source("scripts/figure_aesthetics.R")

daily_file <- file.path("data", "processed", "tahoe_hms_pm25_daily.csv")
master_pm25_file <- file.path("data", "processed", "tahoe_pm25_compiled.csv")
supplement_dir <- file.path("figures", "supplemental")
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(daily_file)) stop("Missing input: ", daily_file)
if (!file.exists(master_pm25_file)) stop("Missing input: ", master_pm25_file)

years <- 2006:2025
category_levels <- c(
  "9.1-35.4",
  ">=35.5"
)
category_labels <- c(
  "9.1-35.4" = "9.1-35.4 µg m⁻³",
  ">=35.5" = "≥35.5 µg m⁻³"
)
category_colours <- c(
  "9.1-35.4" = "#9E9E9E",
  ">=35.5" = "#1A1A1A"
)

daily <- read_csv(
  daily_file,
  col_types = cols(
    date = col_date(), year = col_integer(), month = col_integer(),
    smoke_day = col_logical(), pm25_max = col_double(),
    smoke_confirmed = col_logical(), .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
  filter(year %in% years)

master_pm25 <- read_csv(
  master_pm25_file,
  col_types = cols(
    date = col_date(), PM25 = col_double(), Source = col_character(),
    interpolated = col_logical(), .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
  transmute(
    date,
    year = as.integer(format(date, "%Y")),
    pm25_composite = PM25,
    pm25_source = Source,
    interpolated
  ) %>%
  filter(year %in% years)

daily <- daily %>%
  left_join(master_pm25, by = c("date", "year")) %>%
  mutate(pm25_for_plot = pm25_composite)

smoke_days <- daily %>%
  filter(smoke_day, is.finite(pm25_for_plot), pm25_for_plot >= 9.1) %>%
  mutate(
    pm25_category = case_when(
      pm25_for_plot >= 35.5 ~ ">=35.5",
      TRUE ~ "9.1-35.4"
    ),
    pm25_category = factor(pm25_category, levels = category_levels)
  )

bar_data <- smoke_days %>%
  count(year, pm25_category, name = "smoke_days") %>%
  complete(
    year = years,
    pm25_category = factor(category_levels, levels = category_levels),
    fill = list(smoke_days = 0L)
  ) %>%
  mutate(pm25_category = factor(pm25_category, levels = category_levels))

annual_max <- master_pm25 %>%
  group_by(year) %>%
  summarise(
    pm25_days_available = sum(is.finite(pm25_composite)),
    annual_max_pm25_ug_m3 = if_else(
      pm25_days_available > 0L,
      max(pm25_composite, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

annual_provenance <- master_pm25 %>%
  filter(is.finite(pm25_composite)) %>%
  mutate(pm25_source = coalesce(pm25_source, "Source not recorded")) %>%
  count(year, pm25_source, name = "n_days") %>%
  group_by(year) %>%
  mutate(percent_of_available_days = 100 * n_days / sum(n_days)) %>%
  ungroup() %>%
  arrange(year, desc(n_days), pm25_source)

annual_summary <- bar_data %>%
  pivot_wider(
    names_from = pm25_category, values_from = smoke_days,
    names_prefix = "smoke_days_pm25_"
  ) %>%
  left_join(
    bar_data %>%
      group_by(year) %>%
      summarise(smoke_days_pm25_ge_9_1 = sum(smoke_days), .groups = "drop"),
    by = "year"
  ) %>%
  left_join(annual_max, by = "year") %>%
  mutate(
    smoke_day_definition = paste(
      "NOAA HMS smoke day and composite daily PM2.5 >=9.1 ug/m3;",
      "bars separate 9.1-35.4 and >=35.5 ug/m3"
    ),
    bar_period = "Full calendar year",
    line_definition = paste(
      "Maximum daily composite PM2.5 from tahoe_pm25_compiled.csv",
      "across all monitored days in calendar year"
    )
  )

summary_file <- file.path(
  "data", "processed", "supplemental_smoke_days_annual_summary.csv"
)
figure_only <- identical(Sys.getenv("CALDOR_S1_FIGURE_ONLY"), "1")
if (!figure_only) write_csv(annual_summary, summary_file, na = "")
provenance_file <- file.path(
  "data", "processed", "supplemental_smoke_days_pm25_provenance_by_year.csv"
)
if (!figure_only) write_csv(annual_provenance, provenance_file, na = "")

max_smoke_days <- max(annual_summary$smoke_days_pm25_ge_9_1, na.rm = TRUE)
max_pm25 <- max(annual_summary$annual_max_pm25_ug_m3, na.rm = TRUE)
left_axis_max <- max(10, ceiling((max_smoke_days * 1.18) / 10) * 10)
right_axis_max <- max(50, ceiling((max_pm25 * 1.06) / 50) * 50)
axis_scale <- left_axis_max / right_axis_max

line_data <- annual_max %>%
  mutate(
    year_factor = factor(year, levels = years),
    pm25_scaled = annual_max_pm25_ug_m3 * axis_scale
  )

p_smoke <- ggplot(
  bar_data,
  aes(x = factor(year, levels = years), y = smoke_days, fill = pm25_category)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_line(
    data = line_data,
    aes(x = year_factor, y = pm25_scaled, group = 1),
    inherit.aes = FALSE, colour = "black", alpha = 0.45, linewidth = 0.65
  ) +
  scale_fill_manual(
    values = category_colours,
    labels = category_labels,
    name = expression("PM"[2.5]~"category")
  ) +
  scale_y_continuous(
    limits = c(0, left_axis_max),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0)),
    name = "Smoke days",
    sec.axis = sec_axis(
      ~ . / axis_scale,
      name = expression("Annual maximum PM"[2.5]~"("*mu*"g m"^{-3}*")"),
      breaks = pretty_breaks(n = 6)
    )
  ) +
  scale_x_discrete(
    name = NULL, drop = FALSE,
    labels = function(x) ifelse(as.integer(as.character(x)) %% 2 == 0, x, "")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    colour = "none"
  ) +
  labs(caption = NULL) +
  lo_theme(base_size = 8, family = LO_FONT) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 8.5),
    axis.title.y.left = element_text(size = 8.5),
    axis.title.y.right = element_text(size = 8.5, colour = "grey35"),
    axis.text.y.right = element_text(size = 8.5, colour = "grey35"),
    axis.line.y.right = element_line(colour = "grey35", linewidth = LO_LW_THIN),
    axis.ticks.y.right = element_line(colour = "grey35", linewidth = LO_LW_THIN),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8.5),
    legend.key.width = unit(4.5, "mm"),
    plot.caption = element_blank(),
    panel.grid.minor = element_blank()
  )

figure_file <- file.path(
  supplement_dir, "supplemental_smoke_days_hms_pm25_annual.png"
)
S1_WIDTH_CM <- 7.62
S1_HEIGHT_CM <- S1_WIDTH_CM * 10.5 / 12.7
save_lo_fig(p_smoke, figure_file, width_type = S1_WIDTH_CM, height_cm = S1_HEIGHT_CM)
ggsave(
  sub("\\.png$", ".pdf", figure_file), p_smoke,
  width = S1_WIDTH_CM, height = S1_HEIGHT_CM, units = "cm", device = cairo_pdf,
  bg = "white"
)

cat("\n=== Annual smoke-day summary ===\n")
print(
  annual_summary %>%
    select(year, smoke_days_pm25_ge_9_1, annual_max_pm25_ug_m3),
  n = Inf
)
cat("\nOutputs:\n  ", figure_file, "\n  ", summary_file, "\n", sep = "")
