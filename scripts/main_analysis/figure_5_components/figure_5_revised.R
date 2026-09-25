# =============================================================================
# figure_5_revised.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Revised Figure 5 - SIMPER drivers + distance to season-matched
#          pre-fire centroids using corrected Bray-Curtis PCoA coordinates.
#
# Outputs: figures/figure_5_nmds_braycurtis/revised/
#   figure_5_ltp_<metric>.png  (abundance and biovolume only)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
  library(scales)
  library(patchwork)
  library(ragg)
})

source("scripts/main_analysis/shared_aesthetics.R")

# ---- Paths -------------------------------------------------------------------
PHYTO_FILE <- file.path(
  "data", "lake_environmental_data", "phytoplankton",
  "caldor_fire_phytoplankton_LTP_2005_2025_selected_depths.csv"
)
OUT_DIR    <- file.path("figures", "figure_5_nmds_braycurtis", "archived",
                        "revised_analysis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FAMILY  <- "Times New Roman"
CALDOR_START <- as.Date("2021-08-14")
CALDOR_END   <- as.Date("2021-10-21")
FIRE_YEAR    <- 2021
MIN_PREFIRE_REFS <- 2L
MIN_SEASON_N <- 4L  # min pre-fire samples in a season before using a season-specific 95% band
SIMPER_RECOVERY_END <- c(Surface = "2022-03-08", Deep = "2022-04-06")
LEPTO_TAXON  <- "Leptolyngbya spp."
LEPTO_COL    <- "#009E73"  # matches Leptolyngbya colour in figure_4_phytoplankton_by_depth.R
DIST_Y_MIN   <- 0.25
DIST_Y_MAX   <- 0.9

# Meteorological season (matches pre-fire reference samples by season
# instead of exact calendar month, giving a larger, less gap-prone
# reference pool per sample).
get_season <- function(dates) {
  m <- month(dates)
  case_when(
    m %in% c(12, 1, 2) ~ "Winter",
    m %in% c(3, 4, 5)  ~ "Spring",
    m %in% c(6, 7, 8)  ~ "Summer",
    m %in% c(9, 10, 11) ~ "Fall"
  )
}

# ---- Themes / palettes -------------------------------------------------------
base_theme <- theme_bw(base_size = 8, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
    axis.text        = element_text(size = LO_FS_AXIS_TEXT),
    axis.title       = element_text(size = LO_FS_AXIS_TITLE),
    plot.title       = element_text(size = LO_FS_TITLE, face = "bold"),
    plot.subtitle    = element_text(size = LO_FS_CAPTION, colour = "grey40"),
    legend.text      = element_text(size = 5.8 * LO_FIG_SCALE),
    legend.title     = element_text(size = 6.5 * LO_FIG_SCALE, face = "bold"),
    legend.key.size  = unit(2.0 * LO_FIG_SCALE, "mm"),
    legend.margin    = margin(0, 0, 0, 0, "mm")
  )

period_pal <- c("Pre-fire" = "grey55",
                "Fire"     = "firebrick",
                "Recovery" = "#E69F00",
                "Post"     = "#009E73")

# ---- Date parser -------------------------------------------------------------
parse_excel_date <- function(x) {
  result <- suppressWarnings(as.Date(as.numeric(x), origin = "1899-12-30"))
  iso    <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  as.Date(ifelse(!is.na(iso), as.character(iso), as.character(result)))
}

# ---- Load the curated full 2005-2025 LTP community series -------------------
# The curated table stores biovolume as mm3 m-3, numerically equivalent to
# 10^-6 um3 L-1; convert back to the plotting unit used by the original figure.
clean_ltp <- read_csv(PHYTO_FILE, show_col_types = FALSE) %>%
  transmute(
    station_id = as.character(station_id),
    event_id = as.character(event_id),
    date = as.Date(date),
    depth_m = as.numeric(depth_num),
    taxon = str_squish(as.character(taxon)),
    abundance = replace_na(as.numeric(abundance), 0),
    biovolume = replace_na(as.numeric(biovolume), 0) * 1e6
  ) %>%
  filter(
    station_id == "LTP", !is.na(date), !is.na(taxon), taxon != "",
    !is.na(depth_m), year(date) >= 2005L, year(date) <= 2025L
  )

add_site_label <- function(dat, site_label) {
  dat %>%
    mutate(site = site_label,
           depth_raw = as.character(depth_m),
           depth_num = parse_number(as.character(depth_m))) %>%
    select(site, date, taxon, abundance, biovolume, depth_raw, depth_num)
}

add_depth_bin <- function(dat) {
  dat %>%
    mutate(depth_bin = case_when(
      site == "LTP"  & depth_num >= 0  & depth_num <= 40  ~ "Surface",
      site == "LTP"  & depth_num >= 60 & depth_num <= 105 ~ "Deep",
      site == "MLTP" & depth_raw == "0-100"              ~ "Surface",
      site == "MLTP" & depth_raw == "150-450"            ~ "Deep",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(depth_bin))
}

ltp_site  <- add_site_label(clean_ltp,  "LTP")

build_simper_table <- function(comm, group_results, top_n = 8L) {
  map_dfr(split(comm, comm$line_group), function(grp_comm) {
    group_name <- unique(grp_comm$line_group)
    meta <- grp_comm %>% select(site, depth_bin, date, line_group)
    taxa <- setdiff(names(grp_comm), names(meta))
    M <- grp_comm %>% select(all_of(taxa)) %>% as.matrix()
    keep_taxa <- colSums(M, na.rm = TRUE) > 0
    M <- M[, keep_taxa, drop = FALSE]
    if (ncol(M) == 0) return(tibble())

    depth_name <- unique(meta$depth_bin)
    post_end <- as.Date(unname(SIMPER_RECOVERY_END[depth_name]))
    post_dates <- meta$date[meta$date >= CALDOR_START & meta$date <= post_end]
    if (length(post_dates) < 2) return(tibble())

    post_months <- sort(unique(month(post_dates)))
    sample_status <- meta %>%
      mutate(comparison_group = case_when(
        date < CALDOR_START & month(date) %in% post_months ~ "Month-matched pre-fire",
        date >= CALDOR_START & date <= post_end ~ "Post-fire through recovery",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(comparison_group))

    idx <- match(sample_status$date, meta$date)
    ok <- !is.na(idx)
    sample_status <- sample_status[ok, , drop = FALSE]
    idx <- idx[ok]
    group_fac <- factor(sample_status$comparison_group,
                        levels = c("Month-matched pre-fire", "Post-fire through recovery"))
    counts <- table(group_fac)
    if (length(counts) < 2 || any(counts < 2)) return(tibble())

    H <- decostand(M, method = "hellinger")[idx, , drop = FALSE]
    set.seed(42)
    sim <- vegan::simper(H, group_fac, permutations = 999)
    if (length(sim) == 0) return(tibble())
    sim_df <- as.data.frame(sim[[1]])
    if (nrow(sim_df) == 0) return(tibble())

    taxa_out <- rownames(sim_df)
    pre_matrix <- H[group_fac == "Month-matched pre-fire", , drop = FALSE]
    pre_mean <- colMeans(pre_matrix)
    pre_sd <- apply(pre_matrix, 2, sd, na.rm = TRUE)
    post_matrix <- H[group_fac == "Post-fire through recovery", , drop = FALSE]
    post_mean <- colMeans(post_matrix)
    post_sd <- apply(post_matrix, 2, sd, na.rm = TRUE)
    total_average <- sum(sim_df$average, na.rm = TRUE)

    tibble(
      site = unique(meta$site),
      depth_bin = unique(meta$depth_bin),
      line_group = group_name,
      taxon = taxa_out,
      taxon_label = lo_abbreviate_taxon(taxa_out),
      average_contribution = sim_df$average,
      contribution_pct = 100 * sim_df$average / total_average,
      simper_p = sim_df$p,
      prefire_mean = unname(pre_mean[taxa_out]),
      prefire_sd = unname(pre_sd[taxa_out]),
      postfire_mean = unname(post_mean[taxa_out]),
      postfire_sd = unname(post_sd[taxa_out]),
      higher_in = if_else(postfire_mean >= prefire_mean,
                          "Higher post-fire", "Higher pre-fire"),
      n_prefire = unname(counts["Month-matched pre-fire"]),
      n_postfire = unname(counts["Post-fire through recovery"]),
      postfire_start = min(post_dates),
      postfire_end = max(post_dates),
      postfire_duration_months = ceiling(as.numeric(max(post_dates) - CALDOR_START) / 30.4375),
      matched_months = paste(month.abb[post_months], collapse = ", ")
    ) %>%
      arrange(desc(average_contribution)) %>%
      mutate(rank = row_number(),
             cumulative_pct = cumsum(contribution_pct)) %>%
      slice_head(n = top_n)
  })
}

find_recovery <- function(dat, threshold) {
  dat <- dat %>%
    filter(date >= CALDOR_START, !is.na(dist), n_ref >= MIN_PREFIRE_REFS) %>%
    arrange(date)
  if (nrow(dat) < 3 || is.na(threshold)) {
    return(tibble(recovery_date = as.Date(NA), recovery_months = NA_real_))
  }

  peak_idx <- which.max(dat$dist)
  after_peak <- dat %>% slice(peak_idx:n())
  inside <- after_peak$dist <= threshold
  consecutive <- which(inside & lead(inside, default = FALSE))
  if (length(consecutive) == 0) {
    return(tibble(recovery_date = as.Date(NA), recovery_months = NA_real_))
  }

  recovery_date <- after_peak$date[min(consecutive)]
  tibble(
    recovery_date = recovery_date,
    recovery_months = ceiling(as.numeric(recovery_date - CALDOR_START) / 30.4375)
  )
}

make_simper_plot <- function(simper_results) {
  if (nrow(simper_results) == 0 ||
      !all(c("Surface", "Deep") %in% simper_results$depth_bin)) {
    return(NULL)
  }

  plot_dat <- simper_results %>%
    filter(depth_bin %in% c("Surface", "Deep")) %>%
    group_by(depth_bin) %>%
    slice_min(rank, n = 8, with_ties = FALSE) %>%
    arrange(depth_bin, contribution_pct) %>%
    mutate(
      y = row_number(),
      signed_contribution = if_else(
        higher_in == "Higher pre-fire",
        -contribution_pct,
        contribution_pct
      )
    ) %>%
    ungroup()

  max_x <- max(abs(plot_dat$signed_contribution), na.rm = TRUE) * 1.12

  make_depth_plot <- function(depth_name) {
    depth_dat <- plot_dat %>% filter(depth_bin == depth_name)
    depth_label <- if_else(depth_name == "Surface",
                           "Upper (0-40 m)", "Lower (60-105 m)")
    panel_fill <- if (depth_name == "Deep") "grey92" else "white"

    ggplot(depth_dat, aes(x = signed_contribution, y = y, fill = higher_in)) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = min(depth_dat$y) - 0.5,
               ymax = max(depth_dat$y) + 0.5,
               fill = panel_fill, colour = NA) +
      geom_vline(xintercept = 0, colour = "grey25",
                 linewidth = 0.30 * LO_FIG_SCALE) +
      geom_col(width = 0.72, colour = "grey25",
               linewidth = 0.18 * LO_FIG_SCALE, orientation = "y") +
      annotate("text", x = max_x, y = min(depth_dat$y), label = depth_label,
               hjust = 1, vjust = 0.5, size = 2.05, fontface = "bold",
               family = BASE_FAMILY) +
      scale_y_continuous(breaks = depth_dat$y, labels = depth_dat$taxon_label,
                         expand = expansion(mult = c(0.01, 0.02))) +
      scale_x_continuous(
                         labels = function(x) label_number(accuracy = 1)(abs(x)),
                         limits = c(-max_x, max_x),
                         expand = expansion(mult = c(0, 0))) +
      scale_fill_manual(values = c("Higher post-fire" = "firebrick",
                                   "Higher pre-fire" = "#1B4F72"), name = NULL) +
      labs(x = "Contribution to dissimilarity (%)", y = NULL,
           title = NULL, subtitle = NULL) +
      base_theme +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = ggplot2::element_text(size = 5.5, face = "italic"),
        axis.text.x = ggplot2::element_text(size = 5.5),
        axis.title.x = ggplot2::element_text(size = 6),
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 5.5),
        legend.key.height = unit(1.8, "mm"),
        plot.margin = margin(1, 2, 1, 1, "mm")
      )
  }

  (make_depth_plot("Surface") / make_depth_plot("Deep")) +
    plot_layout(heights = c(1, 1)) &
    theme(legend.position = "none")
}

make_distance_subplot <- function(ts_dat, show_y = TRUE, lepto_dat = NULL,
                                  date_start = as.Date("2019-01-01"),
                                  date_end = as.Date("2025-12-31"),
                                  recovery_label_override = NULL,
                                  date_breaks_override = NULL,
                                  log_lepto = FALSE,
                                  show_bc_points = TRUE,
                                  y_lower = 0.1,
                                  recovery_position = "upper_right",
                                  show_fire_label = TRUE,
                                  recovery_text_size = 3.0,
                                  recovery_fontface = "bold") {
  depth_name <- unique(ts_dat$depth_bin)
  depth_text <- if_else(depth_name == "Surface", "Upper (0-40 m)", "Lower (60-105 m)")
  pt_col <- if (depth_name == "Surface") "#7FB3D5" else "#1B4F72"
  recovery_date <- first(ts_dat$recovery_date)
  recovery_months <- first(ts_dat$recovery_months)
  recovery_label <- if (is.na(recovery_months)) {
    "Recovery: not observed"
  } else {
    paste0("Recovery: ", format(recovery_months, trim = TRUE), " months")
  }
  if (!is.null(recovery_label_override)) recovery_label <- recovery_label_override

  date_range <- c(as.Date(date_start), as.Date(date_end))
  date_break_interval <- if (as.numeric(diff(date_range)) > 15 * 365.25) "2 years" else "1 year"
  if (!is.null(date_breaks_override)) date_break_interval <- date_breaks_override
  visible_bc <- ts_dat %>%
    filter(date >= date_range[1], date <= date_range[2], is.finite(dist)) %>%
    pull(dist)
  if (length(visible_bc) == 0L) visible_bc <- ts_dat$dist[is.finite(ts_dat$dist)]
  bc_span <- diff(range(visible_bc, na.rm = TRUE))
  if (!is.finite(bc_span) || bc_span <= 0) bc_span <- 0.1
  # Focus the displayed Bray-Curtis range on observed variation while keeping
  # the constant historical reference and trajectory interpretable.
  plot_y_min <- y_lower
  y_max <- max(visible_bc, na.rm = TRUE) + 0.10 * bc_span

  pre05 <- first(ts_dat$pre_05)
  premed <- first(ts_dat$pre_median)
  pre95 <- first(ts_dat$pre_95)

  lepto_join <- NULL
  lepto_curve <- NULL
  k <- 1
  if (!is.null(lepto_dat) && nrow(lepto_dat) > 0) {
    lepto_join <- lepto_dat %>%
      filter(date >= date_range[1], date <= date_range[2],
             !is.na(lepto_abundance)) %>%
      arrange(date)
    if (log_lepto) {
      lepto_join <- lepto_join %>%
        # Zero abundance represents absence rather than a finite logarithmic
        # value. Retaining these rows as NA below breaks the curve across the
        # 2013-2021 zero interval instead of drawing an artificial connection.
        mutate(lepto_display = if_else(lepto_abundance > 0,
                                       log(lepto_abundance), NA_real_))
    } else {
      lepto_join <- lepto_join %>% mutate(lepto_display = lepto_abundance)
    }
    max_lepto <- max(lepto_join$lepto_display, na.rm = TRUE)
    if (is.finite(max_lepto) && max_lepto > 0) {
      # Scale the contextual Leptolyngbya series so its observed maximum is
      # approximately halfway up the panel while retaining Bray-Curtis as the
      # primary scale.
      lepto_room <- 0.50 * (y_max - plot_y_min)
      if (lepto_room > 0) {
        lepto_min <- min(lepto_join$lepto_display, na.rm = TRUE)
        lepto_span <- max_lepto - lepto_min
        if (!is.finite(lepto_span) || lepto_span <= 0) lepto_span <- 1
        k <- lepto_room / lepto_span
        lepto_curve <- lepto_join %>%
          mutate(lepto_scaled = if_else(
            is.finite(lepto_display),
            plot_y_min + (lepto_display - lepto_min) * k,
            NA_real_
          ),
          # A long interval without a nonzero observation is not a trajectory.
          # Keep it as a separate group so neither line nor ribbon bridges the
          # 2013-2021 Leptolyngbya absence.
          curve_group = cumsum(replace_na(as.numeric(date - lag(date)) > 180, FALSE)))
      } else {
        lepto_join <- NULL
      }
    } else {
      lepto_join <- NULL
    }
  }

  p <- ggplot(ts_dat, aes(x = date, y = dist)) +
    geom_ribbon(aes(ymin = pre_05, ymax = pre_95, group = 1),
                fill = "grey80", alpha = 0.11, colour = NA) +
    geom_line(aes(y = pre_05, group = 1), linetype = "solid", colour = "grey35",
              linewidth = 0.22 * LO_FIG_SCALE) +
    geom_line(aes(y = pre_median, group = 1), linetype = "solid", colour = "grey45",
              linewidth = 0.20 * LO_FIG_SCALE) +
    geom_line(aes(y = pre_95, group = 1), linetype = "solid", colour = "grey35",
              linewidth = 0.22 * LO_FIG_SCALE)

  if (!is.null(lepto_curve) && nrow(lepto_curve) > 0) {
    p <- p +
      geom_ribbon(
        data = lepto_curve,
        aes(x = date, ymin = plot_y_min, ymax = lepto_scaled, group = curve_group),
        inherit.aes = FALSE, fill = "#009E73", alpha = 0.20,
        colour = NA, na.rm = TRUE
      ) +
      geom_line(
        data = lepto_curve, aes(x = date, y = lepto_scaled, group = curve_group),
        inherit.aes = FALSE, colour = "#006B4F", alpha = 0.70,
        linewidth = 0.42 * LO_FIG_SCALE, na.rm = TRUE
      )
  }

  p <- p +
    # Draw the fire window above both background layers so its red shading is
    # visible over the full y range, including the Leptolyngbya strip.
    annotate("rect", xmin = CALDOR_START, xmax = CALDOR_END,
             ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.08) +
    geom_line(aes(group = 1), colour = pt_col,
              linewidth = 0.55 * LO_FIG_SCALE, na.rm = TRUE) +
    geom_vline(xintercept = CALDOR_START, linetype = "dashed",
               colour = "firebrick", linewidth = 0.40 * LO_FIG_SCALE) +
    annotate("text", x = CALDOR_START - 20, y = y_max,
             label = if (show_fire_label) "Caldor fire" else "",
             angle = 90, hjust = 1, vjust = 0.5,
             size = 3.0, colour = "firebrick", family = BASE_FAMILY) +
    annotate("text",
             x = if (recovery_position == "upper_left") date_range[1] else date_range[2],
             y = y_max, label = recovery_label,
             hjust = if (recovery_position == "upper_left") -0.05 else 1.05,
             vjust = 1.25, size = recovery_text_size,
             fontface = recovery_fontface,
             family = BASE_FAMILY) +
    scale_x_date(date_breaks = date_break_interval, date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                       sec.axis = if (!is.null(lepto_curve) && nrow(lepto_curve) > 0) {
                         sec_axis(~ (. - plot_y_min) / k + lepto_min,
                                  name = NULL,
                                  labels = function(x) {
                                    if (log_lepto) {
                                      formatC(x, format = "f", digits = 1)
                                    } else {
                                      if_else(x == 0, "0", formatC(x, format = "e", digits = 1))
                                    }
                                  })
                       } else waiver()) +
    coord_cartesian(xlim = date_range, ylim = c(plot_y_min, y_max)) +
    labs(title = depth_text, x = NULL, y = NULL) +
    base_theme +
    theme(
      plot.title = ggplot2::element_text(size = LO_FONT_MIN_PT, face = "bold"),
      axis.text = ggplot2::element_text(size = LO_FONT_MIN_PT),
      legend.position = "none",
      plot.margin = margin(1, 1.5, 1, 1.5, "mm"),
      axis.text.y.right = element_text(colour = LEPTO_COL, size = LO_FONT_MIN_PT),
      axis.title.y.right = element_blank(),
      axis.ticks.y.right = element_line(colour = LEPTO_COL)
    )

  if (show_bc_points) {
    p <- p + geom_point(colour = pt_col, size = 1.3 * LO_FIG_SCALE, na.rm = TRUE)
  }

  if (!is.na(recovery_date)) {
    p <- p + geom_vline(xintercept = recovery_date, linetype = "dashed",
                        colour = "black", linewidth = 0.40 * LO_FIG_SCALE)
  }

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

# ---- Core analysis function --------------------------------------------------
make_figure5_rev <- function(dat, value_col, metric) {
  cat("\n=== Processing: LTP |", metric, "===\n")

  comm <- dat %>%
    add_depth_bin() %>%
    group_by(site, depth_bin, date, taxon) %>%
    summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = val, values_fill = 0) %>%
    arrange(depth_bin, date) %>%
    mutate(line_group = paste(site, depth_bin, sep = " | "))

  group_results <- map_dfr(split(comm, comm$depth_bin), function(grp_comm) {
    site_name <- unique(grp_comm$site)
    bin_name <- unique(grp_comm$depth_bin)
    line_group <- unique(grp_comm$line_group)

    M <- grp_comm %>% select(-site, -depth_bin, -date, -line_group) %>% as.matrix()
    keep_rows <- rowSums(M) > 0
    M <- M[keep_rows, , drop = FALSE]
    dates <- grp_comm$date[keep_rows]
    keep_taxa <- colSums(M, na.rm = TRUE) > 0
    M <- M[, keep_taxa, drop = FALSE]

    if (nrow(M) < 6 || ncol(M) < 2) {
      warning("Insufficient community data for ", line_group)
      return(tibble())
    }

    H <- decostand(M, method = "hellinger")
    D <- vegdist(H, method = "bray")
    pc <- cmdscale(D, k = nrow(M) - 1L, eig = TRUE, add = TRUE)
    eig_k <- pc$eig[seq_len(ncol(pc$points))]
    positive <- eig_k > 1e-12
    coords <- pc$points[, positive, drop = FALSE]
    if (ncol(coords) == 0) return(tibble())

    pre_mask <- dates < CALDOR_START
    sample_season <- get_season(dates)
    seasonal <- map_dfr(seq_along(dates), function(i) {
      ref <- which(pre_mask & sample_season == sample_season[i])
      if (pre_mask[i]) ref <- setdiff(ref, i)
      n_ref <- length(ref)
      if (n_ref < MIN_PREFIRE_REFS) {
        return(tibble(date = dates[i], n_ref = n_ref, dist = NA_real_))
      }
      centroid <- colMeans(coords[ref, , drop = FALSE])
      tibble(date = dates[i], n_ref = n_ref,
             dist = sqrt(sum((coords[i, ] - centroid)^2)))
    })

    pre_d <- seasonal$dist[pre_mask & seasonal$n_ref >= MIN_PREFIRE_REFS]
    if (length(pre_d) < MIN_PREFIRE_REFS) return(tibble())
    pre_median <- median(pre_d, na.rm = TRUE)
    pre_95 <- as.numeric(quantile(pre_d, 0.95, na.rm = TRUE, type = 8))
    pre_99 <- as.numeric(quantile(pre_d, 0.99, na.rm = TRUE, type = 8))
    pre_05 <- as.numeric(quantile(pre_d, 0.05, na.rm = TRUE, type = 8))

    # Season-specific 95% range: each date is compared against the pre-fire
    # distribution of distances for *its own season* (falls back to the
    # pooled pre-fire distribution when a season has too few pre-fire samples).
    season_lookup <- tibble(date = dates, season = sample_season,
                            dist = seasonal$dist, n_ref = seasonal$n_ref) %>%
      filter(pre_mask, n_ref >= MIN_PREFIRE_REFS, !is.na(dist)) %>%
      group_by(season) %>%
      summarise(n_season = n(),
                s_premed = median(dist, na.rm = TRUE),
                s_pre95 = as.numeric(quantile(dist, 0.95, na.rm = TRUE, type = 8)),
                s_pre05 = as.numeric(quantile(dist, 0.05, na.rm = TRUE, type = 8)),
                .groups = "drop") %>%
      right_join(tibble(season = c("Winter", "Spring", "Summer", "Fall")), by = "season") %>%
      mutate(
        n_season = replace_na(n_season, 0L),
        season_premed = if_else(n_season >= MIN_SEASON_N, s_premed, pre_median),
        season_pre95  = if_else(n_season >= MIN_SEASON_N, s_pre95, pre_95),
        season_pre05  = if_else(n_season >= MIN_SEASON_N, s_pre05, pre_05)
      ) %>%
      select(season, season_premed, season_pre95, season_pre05)

    seasonal <- seasonal %>%
      mutate(season = sample_season) %>%
      left_join(season_lookup, by = "season")

    post_d <- seasonal %>%
      filter(date > CALDOR_END, n_ref >= MIN_PREFIRE_REFS, !is.na(dist)) %>%
      pull(dist)
    post_95 <- if (length(post_d) >= MIN_PREFIRE_REFS) {
      as.numeric(quantile(post_d, 0.95, na.rm = TRUE, type = 8))
    } else {
      Inf
    }
    extreme_threshold <- max(pre_99, post_95, na.rm = TRUE)
    recovery <- find_recovery(seasonal, pre_95)

    seasonal %>%
      mutate(
        site = site_name,
        depth_bin = bin_name,
        line_group = line_group,
        pre_median = pre_median,
        pre_05 = pre_05,
        pre_95 = pre_95,
        pre_99 = pre_99,
        post_95 = post_95,
        extreme_threshold = extreme_threshold,
        extreme = date > CALDOR_END & n_ref >= MIN_PREFIRE_REFS & !is.na(dist) &
          dist > extreme_threshold,
        recovery_date = recovery$recovery_date,
        recovery_months = recovery$recovery_months
      )
  })

  # Figure 5a and its SIMPER comparison window use the same primary recovery
  # estimate as the resistance/resilience synthesis (season-specific centroid,
  # historical 95th percentile, two consecutive observations). This removes
  # the former one-month discrepancy caused by rounding a separate calculation.
  shared_marker_file <- file.path(
    "figures", "supplemental", "resistance_and_resilience",
    "pathway_resistance_resilience_markers.csv"
  )
  shared_series_file <- file.path(
    "figures", "supplemental", "resistance_and_resilience",
    "pathway_resistance_resilience_timeseries.csv"
  )
  if (metric == "abundance" && file.exists(shared_series_file)) {
    shared_series <- read_csv(shared_series_file, show_col_types = FALSE) %>%
      filter(str_detect(Variable, fixed("Community composition (Bray-Curtis;"))) %>%
      transmute(
        depth_bin = if_else(str_detect(Variable, fixed("0-40 m")), "Surface", "Deep"),
        date = as.Date(date), shared_dist = as.numeric(value),
        shared_expected = as.numeric(expected),
        shared_reference_lo = as.numeric(reference_lo),
        shared_reference_hi = as.numeric(reference_hi),
        shared_n_ref = as.integer(n_hist)
      )
    group_results <- group_results %>%
      inner_join(shared_series, by = c("depth_bin", "date")) %>%
      mutate(
        dist = shared_dist,
        pre_median = shared_expected,
        pre_05 = shared_reference_lo,
        pre_95 = shared_reference_hi,
        n_ref = shared_n_ref
      ) %>%
      select(-starts_with("shared_"))
  }
  if (metric == "abundance" && file.exists(shared_marker_file)) {
    shared_recovery <- read_csv(shared_marker_file, show_col_types = FALSE) %>%
      filter(str_detect(Variable, fixed("Community composition (Bray-Curtis;"))) %>%
      transmute(
        depth_bin = if_else(str_detect(Variable, fixed("0-40 m")), "Surface", "Deep"),
        shared_recovery_date = as.Date(recovery_date),
        shared_recovery_months = as.numeric(recovery_months)
      )
    group_results <- group_results %>%
      select(-recovery_date, -recovery_months) %>%
      left_join(shared_recovery, by = "depth_bin") %>%
      rename(recovery_date = shared_recovery_date,
             recovery_months = shared_recovery_months)
  }

  if (nrow(group_results) == 0) {
    stop("No usable groups for LTP ", metric)
  }

  diagnostics <- group_results %>%
    group_by(depth_bin) %>%
    summarise(
      metric = metric,
      n_prefire_reference = sum(date < CALDOR_START & n_ref >= MIN_PREFIRE_REFS & !is.na(dist)),
      n_postfire = sum(date > CALDOR_END & n_ref >= MIN_PREFIRE_REFS & !is.na(dist)),
      n_above_prefire_95 = sum(date > CALDOR_END & dist > pre_95, na.rm = TRUE),
      n_extreme = sum(extreme, na.rm = TRUE),
      prefire_median = first(pre_median),
      prefire_95 = first(pre_95),
      prefire_99 = first(pre_99),
      postfire_95 = first(post_95),
      extreme_threshold = first(extreme_threshold),
      recovery_date = first(recovery_date),
      recovery_months = first(recovery_months),
      .groups = "drop"
    )
  print(diagnostics, n = Inf)

  simper_results <- build_simper_table(comm, group_results, top_n = 8L)
  p_simper <- make_simper_plot(simper_results)
  if (is.null(p_simper)) {
    stop("SIMPER did not return both LTP depth sections for ", metric)
  }

  lepto_ts <- dat %>%
    add_depth_bin() %>%
    filter(taxon == LEPTO_TAXON) %>%
    group_by(depth_bin, date) %>%
    summarise(lepto_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

  ts_base <- group_results %>%
    filter(n_ref >= MIN_PREFIRE_REFS, !is.na(dist)) %>%
    arrange(depth_bin, date)
  p_surface <- make_distance_subplot(
    filter(ts_base, depth_bin == "Surface"), show_y = TRUE,
    lepto_dat = filter(lepto_ts, depth_bin == "Surface"))
  p_deep <- make_distance_subplot(
    filter(ts_base, depth_bin == "Deep"), show_y = TRUE,
    lepto_dat = filter(lepto_ts, depth_bin == "Deep"))

  y_title <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "Bray-Curtis community dissimilarity to pre-fire centroid",
             angle = 90, size = 2.1, family = BASE_FAMILY) +
    theme_void()
  lepto_title <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = expression(italic("Leptolyngbya")~"abundance (cells L"^-1*")"),
             angle = -90, size = 2.1, colour = LEPTO_COL,
             family = BASE_FAMILY) +
    theme_void()
  p_ts <- y_title + (p_surface / p_deep) + lepto_title +
    plot_layout(widths = c(0.055, 1, 0.055))

  fig5 <- wrap_elements(full = p_ts) + wrap_elements(full = p_simper) +
    plot_layout(widths = c(0.54, 0.46)) +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 8,
                                       family = BASE_FAMILY),
      plot.tag.position = c(0.03, 0.99)
    )

  stub <- paste0("figure_5_ltp_", metric)
  write_csv(group_results, file.path(OUT_DIR, paste0(stub, "_distance.csv")))
  write_csv(diagnostics, file.path(OUT_DIR, paste0(stub, "_diagnostics.csv")))
  write_csv(simper_results, file.path(OUT_DIR, paste0(stub, "_simper.csv")))

  ragg::agg_png(file.path(OUT_DIR, paste0(stub, ".png")),
                width = 17.8, height = 13.5, units = "cm", res = LO_DPI)
  print(fig5)
  invisible(dev.off())
  cat("  Saved:", paste0(stub, ".png"), "\n")

  invisible(list(fig = fig5, distance = group_results,
                 diagnostics = diagnostics, simper = simper_results))
}

results <- list(
  abundance = make_figure5_rev(ltp_site, "abundance", "abundance"),
  biovolume = make_figure5_rev(ltp_site, "biovolume", "biovolume")
)

message("\nDone. LTP abundance and biovolume outputs in: ", OUT_DIR)
