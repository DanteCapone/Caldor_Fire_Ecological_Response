# Build the former Figure 5 abundance-SIMPER panel as manuscript Figure S6.
# The upstream SIMPER calculation remains in figure_5_components/figure_5_revised.R;
# this script only renders its authoritative machine-readable export.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ragg)
})

source("scripts/main_analysis/shared_aesthetics.R")

project_root <- normalizePath(".", winslash = "/")
input_file <- file.path(
  project_root, "figures", "figure_5_nmds_braycurtis", "archived",
  "revised_analysis", "figure_5_ltp_abundance_simper.csv"
)
# Stage beside the authoritative SIMPER export. This existing analysis-output
# directory is also preserved with the former Figure 5 products.
output_dir <- file.path(
  project_root, "figures", "figure_5_nmds_braycurtis", "archived", "revised_analysis"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(
    "Authoritative abundance-SIMPER export is missing: ", input_file,
    "\nRegenerate it with scripts/main_analysis/figure_5_components/figure_5_revised.R"
  )
}

simper <- read_csv(input_file, show_col_types = FALSE)
required_columns <- c(
  "depth_bin", "taxon", "taxon_label", "contribution_pct", "higher_in",
  "postfire_start", "postfire_end", "matched_months", "rank"
)
missing_columns <- setdiff(required_columns, names(simper))
if (length(missing_columns) > 0L) {
  stop("SIMPER export lacks required columns: ", paste(missing_columns, collapse = ", "))
}
if (!setequal(unique(simper$depth_bin), c("Surface", "Deep"))) {
  stop("Expected exactly the Surface and Deep SIMPER strata.")
}
if (any(!simper$higher_in %in% c("Higher pre-fire", "Higher post-fire"))) {
  stop("SIMPER direction is missing or unrecognized for at least one taxon.")
}

plot_data <- simper %>%
  group_by(depth_bin) %>%
  slice_min(rank, n = 8L, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stratum = recode(
      depth_bin,
      Surface = "Upper water column (0\u201340 m)",
      Deep = "Lower water column (60\u2013105 m)"
    ),
    signed_contribution_pct = if_else(
      higher_in == "Higher pre-fire", -contribution_pct, contribution_pct
    )
  )

if (any(!is.finite(plot_data$signed_contribution_pct))) {
  stop("Non-finite contribution in the plotted SIMPER data.")
}

# Use a common signed scale so contribution magnitudes are comparable between panels.
x_limit <- ceiling(max(abs(plot_data$signed_contribution_pct)) / 2) * 2
direction_colors <- c(
  "Higher pre-fire" = "#0072B2",
  "Higher post-fire" = "#D55E00"
)

make_panel <- function(depth_name, show_x = TRUE) {
  depth_data <- plot_data %>%
    filter(depth_bin == depth_name) %>%
    arrange(signed_contribution_pct) %>%
    mutate(taxon_axis = factor(taxon_label, levels = taxon_label))

  ggplot(depth_data, aes(signed_contribution_pct, taxon_axis, fill = higher_in)) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.35) +
    geom_col(width = 0.70, colour = "grey25", linewidth = 0.18) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = seq(-x_limit, x_limit, by = 4),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = direction_colors, name = NULL) +
    labs(
      title = unique(depth_data$stratum),
      x = if (show_x) "Contribution to pre\u2013post-fire dissimilarity (%)" else NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9, base_family = "Times New Roman") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = ggplot2::element_text(size = 8.5, face = "italic", family = "Times New Roman"),
      axis.text.x = ggplot2::element_text(size = 8.5, family = "Times New Roman"),
      axis.title.x = ggplot2::element_text(size = 9, family = "Times New Roman",
                                          margin = margin(t = 5, unit = "pt")),
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5,
                                        family = "Times New Roman", margin = margin(b = 4, unit = "pt")),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 8.5, family = "Times New Roman"),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(6, "mm"),
      plot.margin = margin(4, 3, 3, 8, unit = "mm")
    )
}

figure_s6 <- make_panel("Surface", show_x = FALSE) /
  make_panel("Deep", show_x = TRUE) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = ggplot2::element_text(face = "bold", size = 9, family = "Times New Roman"),
    plot.tag.position = c(0.005, 0.995)
  )

data_output <- file.path(output_dir, "figure_S6_simper_taxon_contributions_data.csv")
png_output <- file.path(output_dir, "figure_S6_simper_taxon_contributions.png")
pdf_output <- file.path(output_dir, "figure_S6_simper_taxon_contributions.pdf")
caption_output <- file.path(output_dir, "figure_S6_simper_taxon_contributions_caption.txt")

write_csv(
  plot_data %>%
    arrange(factor(depth_bin, levels = c("Surface", "Deep")), rank),
  data_output
)

caption <- paste0(
  "Figure S6. Taxon contributions to pre- versus post-fire phytoplankton community ",
  "dissimilarity in (a) the upper water column (0\u201340 m) and (b) the lower water ",
  "column (60\u2013105 m). SIMPER was applied to Bray\u2013Curtis dissimilarities after ",
  "Hellinger transformation of abundance. Post-fire observations were compared with ",
  "month-matched pre-fire observations through the first observed community return ",
  "for each stratum (8 March 2022 for 0\u201340 m; 6 April 2022 for 60\u2013105 m). ",
  "Bars show the eight taxa with the largest average contributions; direction indicates ",
  "whether transformed abundance was greater before or after the fire. Contributions ",
  "are descriptive decompositions of dissimilarity and do not identify causal drivers."
)
writeLines(caption, caption_output, useBytes = TRUE)

ragg::agg_png(
  png_output, width = LO_WIDTH_MAX, height = 13.4, units = "cm", res = LO_DPI_LINE
)
print(figure_s6)
invisible(dev.off())
ggsave(
  pdf_output, figure_s6, width = LO_WIDTH_MAX, height = 13.4,
  units = "cm", device = cairo_pdf
)

message("Saved Figure S6 source products to: ", output_dir)
