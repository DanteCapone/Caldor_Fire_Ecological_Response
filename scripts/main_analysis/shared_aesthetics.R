# =============================================================================
# figure_aesthetics.R
# Caldor Fire Ecosystem Response Project
# Author: Dante A. Capone  |  Date: 2026
#
# Purpose: Shared publication aesthetics for all figures.
#   Standards: Limnology and Oceanography (Wiley) author guidelines
#
# Usage: source("scripts/main_analysis/shared_aesthetics.R")
# =============================================================================

# =============================================================================
# LIMNOLOGY AND OCEANOGRAPHY FIGURE SPECIFICATIONS (Wiley)
# =============================================================================
# Column widths (cm):
LO_WIDTH_SINGLE <- 8.5    # single-column figure
LO_WIDTH_1HALF  <- 11.4   # 1.5-column figure
LO_WIDTH_DOUBLE <- 12.7   # maximum figure width (5 in)
LO_HEIGHT_MAX   <- 15.24  # maximum figure height (6 in)
LO_WIDTH_MAX    <- 12.7   # hard 5-in cap for every helper-based export
LO_DPI_LINE        <- 1200 # L&O minimum for black-and-white or color line art
LO_DPI_COMBINATION <- 600  # L&O minimum for line art combined with halftones
LO_DPI_HALFTONE    <- 300  # L&O minimum for photographic/halftone material
LO_DPI             <- LO_DPI_LINE
LO_DPI_BW          <- LO_DPI_LINE

# Geometry scale for compact preview canvases. Do not use this for font sizes.
LO_FIG_SCALE     <- 0.6

# Font sizes (pt) are evaluated at final exported size. L&O requires text
# greater than 8 pt; dense figures must be split rather than shrinking below it.
LO_FONT_MIN_PT   <- 8.5
LO_FONT_MAX_PT   <- 12
LO_PT_PER_MM     <- 72.27 / 25.4
lo_clamp_pt <- function(size_pt) {
  if (!is.numeric(size_pt)) return(size_pt)
  pmin(LO_FONT_MAX_PT, pmax(LO_FONT_MIN_PT, size_pt))
}
lo_text_mm <- function(size_mm) {
  if (!is.numeric(size_mm)) return(size_mm)
  lo_clamp_pt(size_mm * LO_PT_PER_MM) / LO_PT_PER_MM
}
lo_geom_text_size <- function(size_pt = LO_FONT_MIN_PT) {
  lo_clamp_pt(size_pt) / LO_PT_PER_MM
}

lo_abbreviate_taxon <- function(x) {
  dplyr::recode(
    as.character(x),
    "Tetraedron minimum var. tetralobulatum" = "Tetraedron minimum",
    "Dinobryon sociale var. americanum" = "Dinobryon sociale",
    "Synedra acus var. radians" = "Synedra acus",
    .default = as.character(x)
  )
}

LO_WORD_TARGET_WIDTH_CM <- 10.16  # 4 in; keeps pasted figures compact in Word
LO_WORD_MAX_WIDTH_CM    <- 12.70  # 5 in
LO_WORD_MAX_HEIGHT_CM   <- 15.24  # 6 in
lo_word_width <- function(width_cm) {
  min(LO_WORD_MAX_WIDTH_CM, width_cm * LO_WORD_TARGET_WIDTH_CM / 7.62)
}
lo_word_height <- function(width_cm, height_cm) {
  min(LO_WORD_MAX_HEIGHT_CM, height_cm * lo_word_width(width_cm) / width_cm)
}

LO_FS_AXIS_TEXT  <- 8.0    # tick labels
LO_FS_AXIS_TITLE <- 8.0    # axis titles
LO_FS_STRIP      <- 8.0    # facet strip text
LO_FS_TITLE      <- 9.0    # panel title / plot title
LO_FS_TAG        <- 9.0    # panel tags (a, b, c...)
LO_FS_CAPTION    <- 8.0    # figure caption (in-plot annotations)
LO_FS_LEGEND     <- 8.0    # legend text
LO_FS_LEGEND_TTL <- 8.0    # legend title

# Subpanel tag style — use in plot_annotation(tag_levels = LO_TAG_LEVEL)
# "a" = lowercase  |  "A" = uppercase  |  "1" = numeric
LO_TAG_LEVEL    <- "a"    # L&O convention: lowercase (a, b, c...)

# Line widths (mm → ggplot linewidth units ≈ mm / 0.75)
# Scaled by LO_FIG_SCALE for compact linework on the preview canvas.
LO_LW_THIN       <- 0.3 * LO_FIG_SCALE   # gridlines, ticks
LO_LW_MID        <- 0.6 * LO_FIG_SCALE   # secondary data lines
LO_LW_THICK      <- 0.9 * LO_FIG_SCALE   # primary data lines

# Default font family for all manuscript and supplemental figures.
LO_FONT <- "Times New Roman"

# =============================================================================
# SHARED BASE THEME (L&O-compliant)
# =============================================================================
library(ggplot2)

theme_bw <- function(base_size = 11, base_family = LO_FONT, ...) {
  ggplot2::theme_bw(base_size = lo_clamp_pt(base_size),
                    base_family = base_family, ...)
}

theme_classic <- function(base_size = 11, base_family = LO_FONT, ...) {
  ggplot2::theme_classic(base_size = lo_clamp_pt(base_size),
                         base_family = base_family, ...)
}

theme_minimal <- function(base_size = 11, base_family = LO_FONT, ...) {
  ggplot2::theme_minimal(base_size = lo_clamp_pt(base_size),
                         base_family = base_family, ...)
}

theme_void <- function(base_size = 11, base_family = LO_FONT, ...) {
  ggplot2::theme_void(base_size = lo_clamp_pt(base_size),
                      base_family = base_family, ...)
}

element_text <- function(...) {
  args <- list(...)
  if (!is.null(args$size)) args$size <- lo_clamp_pt(args$size)
  if (is.null(args$family)) args$family <- LO_FONT
  do.call(ggplot2::element_text, args)
}

geom_text <- function(..., size = NULL) {
  args <- list(...)
  if (is.null(args$family)) args$family <- LO_FONT
  if (is.null(size)) {
    do.call(ggplot2::geom_text, c(args, list(size = lo_geom_text_size())))
  } else {
    do.call(ggplot2::geom_text, c(args, list(size = lo_text_mm(size))))
  }
}

geom_label <- function(..., size = NULL) {
  args <- list(...)
  if (is.null(args$family)) args$family <- LO_FONT
  if (is.null(size)) {
    do.call(ggplot2::geom_label, c(args, list(size = lo_geom_text_size())))
  } else {
    do.call(ggplot2::geom_label, c(args, list(size = lo_text_mm(size))))
  }
}

annotate <- function(geom, ...) {
  args <- list(...)
  if (identical(geom, "text") || identical(geom, "label")) {
    if (is.null(args$family)) args$family <- LO_FONT
    if (is.null(args$size)) {
      args$size <- lo_geom_text_size()
    } else {
      args$size <- lo_text_mm(args$size)
    }
  }
  do.call(ggplot2::annotate, c(list(geom = geom), args))
}

lo_theme <- function(base_size = 8, family = LO_FONT) {
  theme_bw(base_size = lo_clamp_pt(base_size), base_family = family) +
    theme(
      text               = element_text(family = family),
      axis.text          = element_text(size = LO_FS_AXIS_TEXT,  family = family),
      axis.title         = element_text(size = LO_FS_AXIS_TITLE, family = family),
      plot.title         = element_text(size = LO_FS_TITLE, face = "bold",
                                        family = family),
      plot.subtitle      = element_text(size = LO_FS_CAPTION, colour = "grey40",
                                        family = family),
      plot.caption       = element_text(size = LO_FS_CAPTION, colour = "grey40",
                                        hjust = 0, family = family),
      legend.title       = element_text(size = LO_FS_LEGEND_TTL, face = "bold",
                                        family = family),
      legend.text        = element_text(size = LO_FS_LEGEND, family = family),
      legend.key.size    = unit(3.5 * LO_FIG_SCALE, "mm"),
      strip.text         = element_text(size = LO_FS_STRIP, face = "bold",
                                        family = family),
      strip.background   = element_rect(fill = "#F0F0F0", colour = "grey70"),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "#EBEBEB", linewidth = LO_LW_THIN),
      plot.margin        = margin(1.5 * LO_FIG_SCALE, 2 * LO_FIG_SCALE,
                                  1.5 * LO_FIG_SCALE, 2 * LO_FIG_SCALE, "mm")
    )
}

# =============================================================================
# SAVE HELPER — enforces L&O dimensions and DPI
# =============================================================================
# width_type: "single" | "1half" | "double" | numeric (cm)
# height_cm : figure height in cm (must be <= LO_HEIGHT_MAX)
# bw        : set TRUE for line art (uses LO_DPI_BW)

# Rule: helper-based manuscript figures must never exceed 5 in wide by 6 in
# high. Numeric widths are therefore clamped to LO_WIDTH_MAX as well.
save_lo_fig <- function(plot, filename, width_type = "double",
                        height_cm = 12, bw = FALSE, ...) {
  requested_w <- switch(as.character(width_type),
                        single  = LO_WIDTH_SINGLE,
                        `1half` = LO_WIDTH_1HALF,
                        double  = LO_WIDTH_DOUBLE,
                        as.numeric(width_type))
  w <- min(requested_w, LO_WIDTH_MAX)
  h <- min(height_cm, LO_HEIGHT_MAX)
  dpi <- if (bw) LO_DPI_BW else LO_DPI

  if (h < height_cm)
    warning("Height clamped from ", height_cm, " to ", LO_HEIGHT_MAX,
            " cm (L&O 6-in maximum).", call. = FALSE)
  if (w < requested_w)
    warning("Width clamped from ", requested_w, " to ", LO_WIDTH_MAX,
            " cm (L&O 5-in maximum).", call. = FALSE)

  ragg::agg_png(filename, width = w, height = h, units = "cm", res = dpi, ...)
  print(plot)
  invisible(dev.off())
  message("Saved (", round(w, 1), " x ", round(h, 1), " cm, ",
          dpi, " dpi): ", basename(filename))
}
