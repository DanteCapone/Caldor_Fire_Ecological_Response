# Two-column Caldor forcing figure using Figure 2 climatology aesthetics.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(ragg)
})
source("scripts/main_analysis/figure_2.R")

FORCING_OUT <- Sys.getenv("CALDOR_FIG1_OUT_DIR", unset = "figures_v2")
dir.create(FORCING_OUT, recursive = TRUE, showWarnings = FALSE)

RAD_BACK <- alpha("#F7CACA", 0.32)
NUT_BACK <- alpha("#CDECCF", 0.34)
RUN_BACK <- alpha("#CCE5F5", 0.38)
DEP_COLS <- c(DIN = "#0072B2", TN = "#D55E00", SRP = "#009E73", TP = "#CC79A7")

pathway_theme <- function(fill) {
  theme(
    plot.background = element_rect(fill = fill, colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5,
                              margin = margin(b = 1.5)),
    plot.subtitle = element_blank()
  )
}

dep_force <- read_csv(
  file.path("data", "lake_environmental_data", "deposition",
            "WY2025_TERC_Midlake_AD_Data_FINAL_edit3.csv"),
  show_col_types = FALSE
) %>%
  filter(QA_Code == "OK") %>%
  mutate(
    start = as.POSIXct(Start_Datetime), end = as.POSIXct(End_Datetime),
    date = as.Date(start + (end - start) / 2),
    year = year(date), month = month(date),
    DIN = DIN_Daily_Load, TN = TN_Daily_Load,
    SRP = SRP_Daily_Load, TP = TP_Daily_Load
  )

dep_long <- dep_force %>%
  select(date, year, month, DIN, TN, SRP, TP) %>%
  pivot_longer(c(DIN, TN, SRP, TP), names_to = "nutrient", values_to = "value") %>%
  filter(is.finite(value))

dep_clim <- dep_long %>%
  filter(year != 2021) %>%
  group_by(nutrient, month) %>%
  summarise(n = n(), clim_mean = mean(value), clim_sd = sd(value), .groups = "drop") %>%
  mutate(
    clim_ci = qt(.975, pmax(n - 1, 1)) * clim_sd / sqrt(n),
    clim_lo = pmax(clim_mean - clim_ci, 0),
    clim_hi = clim_mean + clim_ci,
    plot_date = as.Date(sprintf("2021-%02d-15", month))
  )
dep_2021 <- dep_long %>% filter(year == 2021)

make_dep_panel <- function(vars, labels, y_lab, title = NULL) {
  clim <- dep_clim %>% filter(nutrient %in% vars)
  obs <- dep_2021 %>% filter(nutrient %in% vars)
  ggplot() +
    fire_layers +
    geom_ribbon(
      data = clim,
      aes(plot_date, ymin = clim_lo, ymax = clim_hi, fill = nutrient,
          group = nutrient), alpha = .13, colour = NA
    ) +
    geom_line(data = clim, aes(plot_date, clim_mean, colour = nutrient),
              linewidth = .48, alpha = .55) +
    geom_line(data = obs, aes(date, value, colour = nutrient),
              linewidth = .62, alpha = .82) +
    geom_point(data = obs, aes(date, value, colour = nutrient),
               size = 1.25, alpha = .90) +
    scale_colour_manual(values = DEP_COLS[vars], labels = labels, name = NULL) +
    scale_fill_manual(values = DEP_COLS[vars], labels = labels, name = NULL) +
    year_axis +
    labs(y = y_lab, title = title) +
    guides(fill = "none", colour = guide_legend(ncol = 1, byrow = TRUE)) +
    base_theme +
    theme(
      legend.position = c(.09, .84),
      legend.justification = c(0, 0.5),
      legend.background = element_rect(fill = alpha("white", .82), colour = "grey80"),
      legend.text = element_text(size = 5.3),
      legend.key.width = unit(3.6, "mm"),
      legend.key.height = unit(2.2, "mm")
    ) +
    pathway_theme(NUT_BACK)
}

ratio_force <- dep_force %>%
  filter(is.finite(DIN), is.finite(SRP), SRP > 0) %>%
  transmute(date, year, month,
            value = (DIN / 14.0067) / (SRP / 30.973762))
ratio_clim <- ratio_force %>%
  filter(year != 2021) %>%
  group_by(month) %>%
  summarise(n = n(), clim_mean = mean(value), clim_sd = sd(value), .groups = "drop") %>%
  mutate(
    clim_ci = qt(.975, pmax(n - 1, 1)) * clim_sd / sqrt(n),
    clim_lo = pmax(clim_mean - clim_ci, 0),
    clim_hi = clim_mean + clim_ci,
    plot_date = as.Date(sprintf("2021-%02d-15", month))
  )

p_ratio <- ggplot() +
  fire_layers +
  geom_ribbon(data = ratio_clim, aes(plot_date, ymin = clim_lo, ymax = clim_hi),
              fill = COL_CLIM, alpha = .24) +
  geom_line(data = ratio_clim, aes(plot_date, clim_mean),
            colour = COL_CLIM, linewidth = .55) +
  geom_line(data = ratio_force %>% filter(year == 2021), aes(date, value),
            colour = "#6A3D9A", linewidth = .62) +
  geom_point(data = ratio_force %>% filter(year == 2021), aes(date, value),
             colour = "#6A3D9A", size = 1.3) +
  year_axis +
  labs(y = "Atmospheric N:P") +
  base_theme + pathway_theme(NUT_BACK)

flow_raw <- read_csv(
  file.path("data", "processed", "usgs_tahoe_tributary_streamflow_daily_2015_2025.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date), year = year(date), month = month(date))

# The historical reference is the monthly mean of the combined daily flows,
# while 2021 retains the native daily observations for each tributary.
flow_daily <- flow_raw %>%
  group_by(date, year, month) %>%
  summarise(value = sum(flow_cms, na.rm = TRUE), n_sites = n_distinct(site_no),
            .groups = "drop") %>%
  filter(n_sites >= 2)
flow_month <- flow_daily %>%
  group_by(year, month) %>%
  summarise(value = mean(value), .groups = "drop")
flow_clim <- monthly_climatology(flow_month, "value", 2015L)
flow_2021 <- flow_raw %>%
  filter(year == 2021) %>%
  transmute(date, stream = site, value = flow_cms)

FLOW_COLS <- c("Upper Truckee River" = "#1F78B4", "Blackwood Creek" = "#E66101")

p_runoff <- ggplot() +
  fire_layers +
  geom_ribbon(data = flow_clim, aes(plot_date, ymin = clim_lo, ymax = clim_hi),
              fill = COL_CLIM, alpha = 0.24, colour = NA) +
  geom_line(data = flow_clim, aes(plot_date, clim_mean),
            colour = COL_CLIM, linewidth = 0.65) +
  geom_point(data = flow_clim, aes(plot_date, clim_mean),
             colour = COL_CLIM, size = 1.15) +
  geom_line(data = flow_2021, aes(date, value, colour = stream),
            linewidth = 0.48, alpha = 0.9) +
  scale_colour_manual(values = FLOW_COLS, name = NULL) +
  year_axis +
  labs(y = "Discharge", title = "Runoff pathway", tag = "h") +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE)) +
  base_theme +
  theme(
    legend.position = c(0.03, 0.84),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = alpha("white", 0.82), colour = "grey80"),
    legend.text = element_text(size = 5.3),
    legend.key.width = unit(4.2, "mm"),
    legend.key.height = unit(2.2, "mm")
  ) +
  pathway_theme(RUN_BACK)

# The source figure uses hollow circles to flag non-QC observations. Figure 1
# retains those observations; selected panels preserve that explicit QC flag.
solidify_open_circle_layer <- function(plot) {
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (inherits(layer$geom, "GeomPoint") && identical(layer$aes_params$shape, 1)) {
      layer$aes_params$shape <- 16
      plot$layers[[i]] <- layer
    }
  }
  plot
}

# Fold flagged observations into the original 2021 line layer, avoiding a
# duplicate trajectory while keeping every observation connected.
connect_flagged_points <- function(plot) {
  flagged_idx <- which(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomPoint") && identical(layer$aes_params$shape, 1),
    logical(1)
  ))
  line_idx <- tail(which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))), 1)
  if (!length(flagged_idx) || !length(line_idx)) return(plot)

  observed <- plot$layers[[line_idx]]$data
  flagged <- plot$layers[[flagged_idx[[1]]]]$data
  if (!all(c("date", "value") %in% names(observed)) ||
      !all(c("date", "value") %in% names(flagged))) return(plot)

  plot$layers[[line_idx]]$data <- bind_rows(observed, flagged) %>% arrange(date)
  plot
}

p_rad_par <- connect_flagged_points(p_a) + labs(tag = "a", title = "Radiative pathway") + pathway_theme(RAD_BACK)
p_rad_uv <- connect_flagged_points(p_b) + labs(tag = "c") + pathway_theme(RAD_BACK)
p_rad_secchi <- solidify_open_circle_layer(connect_flagged_points(p_c)) + labs(tag = "e") + pathway_theme(RAD_BACK)
p_rad_temp <- solidify_open_circle_layer(connect_flagged_points(p_d)) +
  labs(tag = "g", y = expression(bold("Mean 0-10 m temperature ("*degree*"C)"))) +
  pathway_theme(RAD_BACK)

p_n <- make_dep_panel(
  c("DIN", "TN"), c(DIN = "DIN", TN = "TN"),
  "N deposition\n(mg N m⁻² d⁻¹)",
  "Nutrient pathway"
) + labs(tag = "b")
p_p <- make_dep_panel(
  c("SRP", "TP"), c(SRP = "SRP", TP = "TP"),
  "P deposition\n(mg P m⁻² d⁻¹)"
) + labs(tag = "d")
p_ratio <- p_ratio + labs(tag = "f")

# Only the bottom row displays month labels; every panel retains a unit-bearing
# y-axis as in the original Figure 2.
blank_upper_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_rad_par <- p_rad_par + blank_upper_x
p_rad_uv <- p_rad_uv + blank_upper_x
p_rad_secchi <- p_rad_secchi + blank_upper_x
p_n <- p_n + blank_upper_x
p_p <- p_p + blank_upper_x
p_ratio <- p_ratio + blank_upper_x

# Keep all Figure 1 panel labels consistently legible at manuscript scale.
figure1_label_theme <- theme(
  axis.title = element_text(size = 8.5, face = "bold"),
  axis.title.y = element_text(size = 7.5, face = "bold"),
  axis.text = element_text(size = 8.5, face = "bold"),
  axis.text.y = element_text(size = 7.5, face = "bold"),
  legend.title = element_text(size = 8.5, face = "bold"),
  legend.text = element_text(size = 8.5, face = "bold"),
  plot.title = element_text(size = 9, face = "bold"),
  plot.tag = element_text(size = 9, face = "bold", hjust = 1, vjust = 1),
  plot.tag.position = c(.985, .985)
)
p_rad_par <- p_rad_par + figure1_label_theme
p_rad_uv <- p_rad_uv + figure1_label_theme
p_rad_secchi <- p_rad_secchi + figure1_label_theme
p_rad_temp <- p_rad_temp + figure1_label_theme
p_n <- p_n + figure1_label_theme
p_p <- p_p + figure1_label_theme
p_ratio <- p_ratio + figure1_label_theme
p_runoff <- p_runoff + figure1_label_theme

# Figure 1 uses line trajectories only; point layers are removed after the
# source plots have supplied their lines and reference envelopes.
remove_point_layers <- function(plot, preserve_open = FALSE) {
  keep <- vapply(plot$layers, function(layer) {
    if (!inherits(layer$geom, "GeomPoint")) return(TRUE)
    preserve_open && identical(layer$aes_params$shape, 1)
  }, logical(1))
  plot$layers <- plot$layers[keep]
  plot
}
p_rad_par <- remove_point_layers(p_rad_par, preserve_open = TRUE)
p_rad_uv <- remove_point_layers(p_rad_uv, preserve_open = TRUE)
p_rad_secchi <- remove_point_layers(p_rad_secchi)
p_rad_temp <- remove_point_layers(p_rad_temp)
p_n <- remove_point_layers(p_n)
p_p <- remove_point_layers(p_p)
p_ratio <- remove_point_layers(p_ratio)
p_runoff <- remove_point_layers(p_runoff)

caldor_forcing_2col <-
  (p_rad_par | p_n) /
  (p_rad_uv | p_p) /
  (p_rad_secchi | p_ratio) /
  (p_rad_temp | p_runoff) +
  plot_layout(widths = c(1, 1), heights = c(1, 1, 1, 1))

ggsave(file.path(FORCING_OUT, "fig01_caldor_forcing.png"), caldor_forcing_2col,
       width = 12.7, height = 15.24, units = "cm", dpi = 1200,
       device = ragg::agg_png, bg = "white")
ggsave(file.path(FORCING_OUT, "fig01_caldor_forcing.pdf"), caldor_forcing_2col,
       width = 12.7, height = 15.24, units = "cm", device = cairo_pdf, bg = "white")

# Poster-ready pathway extracts use the same panels and encodings as Figure 1.
POSTER_OUT <- file.path("figures", "for_poster")
dir.create(POSTER_OUT, recursive = TRUE, showWarnings = FALSE)
poster_panel_theme <- theme(
  plot.background = element_rect(fill = "white", colour = NA),
  axis.text.x = element_text(size = 13, face = "bold", colour = "black"),
  axis.text.y = element_text(size = 13, face = "bold", colour = "black"),
  axis.ticks.x = element_line(),
  axis.title.y = element_text(size = 16, face = "bold"),
  axis.title.x = element_text(size = 16, face = "bold"),
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 13, hjust = 0.5)
)
poster_month_axis <- scale_x_date(
  limits = c(as.Date("2021-01-01"), as.Date("2021-12-31")),
  date_breaks = "2 months", date_labels = "%b",
  expand = expansion(mult = 0.01)
)

poster_radiative <-
  solidify_open_circle_layer(connect_flagged_points(p_a)) +
  poster_month_axis +
  labs(x = NULL, tag = NULL, title = NULL) +
  poster_panel_theme +
  plot_annotation(
    title = "Radiative pathways",
    theme = theme(plot.title = element_text(size = 22, hjust = 0.5,
                                            face = "bold"))
  )

lake_nutrients <- read_csv(
  file.path("data", "lake_environmental_data", "nutrients", "Tahoe_MLTP_Nutrient.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(Date), year = year(date), month = month(date)) %>%
  filter(Depth <= 10) %>%
  group_by(date, year, month) %>%
  summarise(
    DIN = mean(as.numeric(NO3) + as.numeric(NH4), na.rm = TRUE),
    SRP = mean(as.numeric(TRP), na.rm = TRUE), .groups = "drop"
  ) %>%
  filter(is.finite(DIN), is.finite(SRP))
lake_force <- lake_nutrients %>% filter(year == 2021) %>% select(date, DIN, SRP)

dep_nutrients <- dep_force %>%
  filter(is.finite(DIN), is.finite(SRP)) %>%
  transmute(date, year, month, DIN, SRP)
dep_force_2021 <- dep_nutrients %>% filter(year == 2021) %>% select(date, DIN, SRP)
make_nutrient_climatology <- function(data) {
  data %>%
    filter(year != 2021) %>%
    pivot_longer(c(DIN, SRP), names_to = "nutrient", values_to = "value") %>%
    group_by(nutrient, month) %>%
    summarise(
      n = sum(is.finite(value)),
      clim_mean = mean(value, na.rm = TRUE),
      clim_sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      clim_ci = qt(0.975, pmax(n - 1, 1)) * clim_sd / sqrt(n),
      clim_lo = pmax(clim_mean - clim_ci, 0),
      clim_hi = clim_mean + clim_ci,
      date = as.Date(sprintf("2021-%02d-15", month))
    )
}

lake_clim <- make_nutrient_climatology(lake_nutrients)
dep_clim <- make_nutrient_climatology(dep_nutrients)

make_nutrient_panel <- function(data, clim, nutrient, colour, y_label, title = NULL,
                                y_position = "left") {
  observed <- data %>% transmute(date, value = .data[[nutrient]])
  reference <- clim %>% filter(nutrient == !!nutrient)
  ggplot() +
    fire_layers +
    geom_ribbon(data = reference, aes(date, ymin = clim_lo, ymax = clim_hi),
                fill = colour, alpha = 0.16, colour = NA) +
    geom_line(data = reference, aes(date, clim_mean),
              colour = colour, alpha = 0.38, linewidth = 0.72) +
    geom_point(data = reference, aes(date, clim_mean),
               colour = colour, alpha = 0.38, size = 1.35) +
    geom_line(data = observed, aes(date, value), colour = colour, linewidth = 0.95) +
    geom_point(data = observed, aes(date, value), colour = colour, size = 2.0) +
    poster_month_axis +
    scale_y_continuous(name = y_label, position = y_position) +
    labs(x = NULL, y = y_label, title = title) +
    base_theme + poster_panel_theme +
    theme(
      axis.title.y.left = element_text(size = 14, face = "bold", colour = "black"),
      axis.title.y.right = element_text(size = 14, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 13, face = "bold", colour = "black"),
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      plot.margin = margin(3, 5, 2, 5, "mm")
    )
}

p_dep_din <- make_nutrient_panel(
  dep_force_2021, dep_clim, "DIN", "#D55E00",
  expression("Deposition (mg N m"^{-2}~d^{-1}*")"), "Deposition"
)
p_lake_din <- make_nutrient_panel(
  lake_force, lake_clim, "DIN", "#0072B2",
  expression("Concentration ("*mu*"g N L"^{-1}*")"), "In-lake (0-10 m)", "right"
)
p_dep_srp <- make_nutrient_panel(
  dep_force_2021, dep_clim, "SRP", "#F4A261",
  expression("Deposition (mg P m"^{-2}~d^{-1}*")")
)
p_lake_srp <- make_nutrient_panel(
  lake_force, lake_clim, "SRP", "#56B4E9",
  expression("Concentration ("*mu*"g P L"^{-1}*")"), NULL, "right"
)

make_nutrient_row_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90,
             size = 6.0, fontface = "bold", colour = "black") +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 0, "mm"))
}

poster_nutrient <- (
  make_nutrient_row_label("Dissolved inorganic nitrogen") | p_dep_din | p_lake_din
) + plot_layout(widths = c(0.055, 1, 1)) &
  theme(plot.background = element_rect(fill = NUT_BACK, colour = NA))
poster_nutrient <- poster_nutrient + plot_annotation(
  title = "Nutrient pathways",
  theme = theme(plot.title = element_text(size = 22, hjust = 0.5,
                                          face = "bold"))
)

poster_runoff <- p_runoff + labs(x = NULL, tag = NULL, title = "Runoff pathway") +
  poster_panel_theme +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 2, unit = "mm")
  )

for (spec in list(
  list(plot = poster_radiative, file = "radiative_pathways.png", width = 15.0, height = 11.0, dpi = 300),
  list(plot = poster_nutrient, file = "nutrient_pathways.png", width = 25.4, height = 11.0, dpi = 300),
  list(plot = poster_runoff, file = "runoff_pathway.png", width = 10.0, height = 8.2, dpi = 600)
)) {
  ggsave(file.path(POSTER_OUT, spec$file), spec$plot,
         width = spec$width, height = spec$height, units = "cm", dpi = spec$dpi,
         device = ragg::agg_png, bg = "white")
}

message("Two-column Caldor forcing figure written without subtitle or timeline strip.")
