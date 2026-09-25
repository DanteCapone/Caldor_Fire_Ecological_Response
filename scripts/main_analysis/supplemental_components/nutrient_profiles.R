# =============================================================================
# supplemental_nutrient_profiles.R
# Caldor Fire Ecosystem Response Project
#
# Purpose: Vertical (depth) distribution of nutrients at the MLTP mid-lake
#          station, Jul 2021 - Apr 2022, one panel per sampling month with
#          each nutrient (NO3, NH4, TKN, TRP, THP) drawn as a differently
#          colored depth profile. Mirrors the layout convention used for
#          figures/supplemental/ctd_profiles/.
#
# Data: data/lake_environmental_data/nutrients/Tahoe_MLTP_Nutrient.csv
#       (full-depth-profile nutrient casts: 0-450 m)
#
# Output: figures/supplemental/nutrient_profiles/
#   nutrient_profiles_mltp_jul2021_apr2022.png (combined reference)
#   nutrient_profile_mltp_<nutrient>_jul2021_apr2022.png (one figure per nutrient)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(ragg)
})
source("scripts/main_analysis/shared_aesthetics.R")

NUTRI_MLTP <- "data/lake_environmental_data/nutrients/Tahoe_MLTP_Nutrient.csv"
OUT_DIR <- "figures/supplemental/nutrient_profiles"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BASE_FAMILY <- "Times New Roman"
PROFILE_START <- as.Date("2021-07-01")
PROFILE_END   <- as.Date("2022-04-30")
CALDOR_START  <- as.Date("2021-08-14")
CALDOR_END    <- as.Date("2021-10-21")

nut_cols <- c(NO3 = "NO3", NH4 = "NH4", TKN = "TKN", TRP = "TRP", THP = "THP")
nut_levels <- names(nut_cols)
nut_display <- c(NO3 = "NO\u2083\u207B", NH4 = "NH\u2084\u207A", TKN = "TKN",
                 TRP = "TRP", THP = "THP")
nut_cols_display <- c(NO3 = "#1B9E77", NH4 = "#D95F02", TKN = "#7570B3",
                      TRP = "#E7298A", THP = "#66A61E")

raw <- read_csv(NUTRI_MLTP, show_col_types = FALSE) %>%
  mutate(date = as.Date(Date), depth = as.numeric(Depth)) %>%
  filter(date >= PROFILE_START, date <= PROFILE_END)

if (nrow(raw) == 0) stop("No MLTP nutrient casts found in ", PROFILE_START, " to ", PROFILE_END)

profiles <- raw %>%
  select(date, depth, all_of(unname(nut_cols))) %>%
  rename(!!!setNames(unname(nut_cols), names(nut_cols))) %>%
  pivot_longer(cols = all_of(nut_levels), names_to = "nutrient", values_to = "value") %>%
  mutate(
    value = suppressWarnings(as.numeric(value)),
    nutrient = factor(nutrient, levels = nut_levels),
    month_label = factor(format(date, "%b %Y"),
                         levels = format(sort(unique(date)), "%b %Y"))
  ) %>%
  filter(!is.na(value), !is.na(depth))

write_csv(
  profiles %>% arrange(date, nutrient, depth),
  file.path(OUT_DIR, "nutrient_profiles_mltp_jul2021_apr2022_data.csv")
)

n_months <- n_distinct(profiles$month_label)
message("Sampling months: ", paste(levels(profiles$month_label), collapse = ", "))

nutrient_theme <- theme_classic(base_size = 7, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY),
    axis.title = element_text(size = 7, family = BASE_FAMILY),
    axis.text = element_text(size = 6, family = BASE_FAMILY),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7, family = BASE_FAMILY),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 6.5, family = BASE_FAMILY),
    legend.text = element_text(size = 6.5, family = BASE_FAMILY),
    legend.key.size = unit(2.5, "mm"),
    panel.grid.major.y = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN)
  )

fire_month_flag <- profiles %>%
  distinct(date, month_label) %>%
  mutate(in_fire_window = date >= floor_date(CALDOR_START, "month") &
           date <= ceiling_date(CALDOR_END, "month") - days(1))

p_full <- ggplot(profiles, aes(x = value, y = depth, colour = nutrient)) +
  geom_path(linewidth = LO_LW_THICK, na.rm = TRUE) +
  geom_point(size = 1.1, na.rm = TRUE) +
  facet_wrap(~ month_label, nrow = 2,
             labeller = labeller(month_label = label_value)) +
  scale_y_reverse(breaks = c(0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450)) +
  scale_colour_manual(values = nut_cols_display, breaks = nut_levels,
                      labels = nut_display, name = "Nutrient") +
  labs(x = expression("Concentration (\u00b5g L"^-1*")"), y = "Depth (m)",
       title = "MLTP nutrient depth profiles: Jul 2021 - Apr 2022") +
  nutrient_theme +
  theme(plot.title = element_text(face = "bold", size = 8, family = BASE_FAMILY))

out_path_full <- file.path(OUT_DIR, "nutrient_profiles_mltp_jul2021_apr2022.png")
ggsave(out_path_full, p_full, width = 26, height = 14, units = "cm",
       dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
message("Saved: ", out_path_full)

# ---- Same figure restricted to the upper 100 m (bloom-relevant depths) ----
p_upper <- p_full %+% (profiles %>% filter(depth <= 100)) +
  scale_y_reverse(breaks = c(0, 10, 25, 50, 75, 100)) +
  labs(title = "MLTP nutrient depth profiles, 0-100 m: Jul 2021 - Apr 2022")

out_path_upper <- file.path(OUT_DIR, "nutrient_profiles_mltp_jul2021_apr2022_0_100m.png")
ggsave(out_path_upper, p_upper, width = 26, height = 14, units = "cm",
       dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
message("Saved: ", out_path_upper)

# ---- Separate figures for each nutrient -------------------------------------
# Each nutrient gets its own figure so the x-axis can be scaled to that
# nutrient's concentration range rather than sharing one axis across nutrients.
nutrient_profile_plot <- function(nutrient_id, depth_limit = NULL) {
  dat <- profiles %>% filter(nutrient == nutrient_id)
  if (!is.null(depth_limit)) dat <- dat %>% filter(depth <= depth_limit)

  y_breaks <- if (is.null(depth_limit)) {
    c(0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450)
  } else {
    c(0, 10, 25, 50, 75, 100)
  }
  depth_label <- if (is.null(depth_limit)) "" else ", 0-100 m"

  ggplot(dat, aes(x = value, y = depth)) +
    geom_path(colour = nut_cols_display[[nutrient_id]], linewidth = LO_LW_THICK,
              na.rm = TRUE) +
    geom_point(colour = nut_cols_display[[nutrient_id]], size = 1.1, na.rm = TRUE) +
    facet_wrap(~ month_label, nrow = 2, labeller = labeller(month_label = label_value)) +
    scale_x_continuous(labels = label_scientific(digits = 2),
                       expand = expansion(mult = c(0.03, 0.10))) +
    scale_y_reverse(breaks = y_breaks, limits = if (is.null(depth_limit)) NULL else c(depth_limit, 0)) +
    labs(x = expression("Concentration (\u00b5g L"^-1*")"), y = "Depth (m)",
         title = paste0("MLTP ", nut_display[[nutrient_id]],
                        " depth profiles: Jul 2021 - Apr 2022", depth_label)) +
    nutrient_theme +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 8, family = BASE_FAMILY))
}

walk(nut_levels, function(nutrient_id) {
  nutrient_name <- tolower(nutrient_id)
  out_path <- file.path(OUT_DIR, paste0(
    "nutrient_profile_mltp_", nutrient_name, "_jul2021_apr2022.png"))
  ggsave(out_path, nutrient_profile_plot(nutrient_id), width = 26, height = 14,
         units = "cm", dpi = LO_DPI, device = ragg::agg_png, limitsize = FALSE)
  message("Saved: ", out_path)

  out_path_upper <- file.path(OUT_DIR, paste0(
    "nutrient_profile_mltp_", nutrient_name, "_jul2021_apr2022_0_100m.png"))
  ggsave(out_path_upper, nutrient_profile_plot(nutrient_id, depth_limit = 100),
         width = 26, height = 14, units = "cm", dpi = LO_DPI,
         device = ragg::agg_png, limitsize = FALSE)
  message("Saved: ", out_path_upper)
})
