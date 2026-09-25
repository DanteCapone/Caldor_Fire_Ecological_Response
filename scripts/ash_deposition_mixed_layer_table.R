# Recalculate the ash-deposition concentration equivalents for the September
# 2021 mixed layer rather than the full Lake Tahoe volume.

suppressPackageStartupMessages({
  library(tidyverse)
})

OUT_DIR <- file.path("figures", "supplemental_final")
PROC_DIR <- file.path("data", "processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COLLECTION_START <- as.Date("2021-08-26")
COLLECTION_END <- as.Date("2021-09-10")
COLLECTION_DAYS <- 16
LAKE_AREA_M2 <- 490e6
MLD_DATE <- as.Date("2021-09-01")

mld_source <- read_csv(
  file.path(PROC_DIR, "mltp_mixed_layer_depth_density_2005_2025.csv"),
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  filter(date == MLD_DATE)

if (nrow(mld_source) != 1 || !is.finite(mld_source$mld_m[[1]])) {
  stop("Expected one finite MLTP mixed-layer depth on 2021-09-01.")
}

MLD_M <- mld_source$mld_m[[1]]

# Mean rates were transcribed from the manuscript table supplied for this task.
# Totals are recomputed rather than copied so the table is formula-derived.
ash_rates <- tribble(
  ~element, ~mean_rate_ug_m2_d,
  "TP", 3090.4,
  "Mg", 4613.4,
  "K", 8676.1,
  "Ca", 36389.1,
  "Na", 3362.6,
  "Fe", 6197.7,
  "Zn", 97.1,
  "Cu", 19.5,
  "Co", 3.5,
  "Mn", 2187.8,
  "As", 0.8,
  "Cd", 0.7,
  "Cr", 7.0,
  "Ni", 5.8,
  "Pb", 8.3
)

ash_mixed_layer <- ash_rates %>%
  mutate(
    collection_days = COLLECTION_DAYS,
    lake_area_m2 = LAKE_AREA_M2,
    mld_date = MLD_DATE,
    mixed_layer_depth_m = MLD_M,
    areal_load_ug_m2 = mean_rate_ug_m2_d * collection_days,
    total_deposition_tons = areal_load_ug_m2 * lake_area_m2 / 1e12,
    mixed_layer_concentration_ug_L = areal_load_ug_m2 /
      (mixed_layer_depth_m * 1000)
  )

write_csv(
  ash_mixed_layer,
  file.path(PROC_DIR, "ash_deposition_mixed_layer_equivalent_2021_08_26_09_10.csv")
)

table_display <- ash_mixed_layer %>%
  transmute(
    Elements = element,
    `Mean deposition rate (µg m⁻² d⁻¹)` = formatC(
      mean_rate_ug_m2_d, format = "f", digits = 1
    ),
    `Total 16-day deposition (tons)` = formatC(
      total_deposition_tons, format = "f", digits = 3
    ),
    `13-m mixed-layer concentration equivalent (µg L⁻¹)` = formatC(
      mixed_layer_concentration_ug_L, format = "f", digits = 4
    )
  )

write_tsv(
  table_display,
  file.path(OUT_DIR, "Table_S1_ash_deposition_mixed_layer.tsv")
)

header <- paste(names(table_display), collapse = " | ")
separator <- paste(rep("---", ncol(table_display)), collapse = " | ")
rows <- apply(table_display, 1, paste, collapse = " | ")
caption <- paste(
  "Table S1. Mean elemental deposition rates across 10 sampling sites,",
  "estimated total loading to Lake Tahoe during the 16-day collection period,",
  "and corresponding concentration equivalents for the 13-m mixed layer",
  "observed at MLTP on 1 September 2021."
)
note <- paste(
  "Note: Estimates assume a lake surface area of 490 km².",
  "The concentration equivalent divides the 16-day areal load by the observed",
  "13-m mixed-layer depth and assumes lake-wide deposition, complete retention,",
  "and uniform mixing within that layer. It does not account for solubility,",
  "bioavailability, spatial heterogeneity, or removal."
)

md_path <- file.path(OUT_DIR, "Table_S1_ash_deposition_mixed_layer.md")
writeLines(
  c(caption, "", paste0("| ", header, " |"),
    paste0("| ", separator, " |"), paste0("| ", rows, " |"), "", note),
  md_path,
  useBytes = TRUE
)

docx_path <- file.path(OUT_DIR, "Table_S1_ash_deposition_mixed_layer.docx")
pandoc_status <- system2(
  "quarto",
  c("pandoc", shQuote(md_path), "--output", shQuote(docx_path))
)
if (!identical(pandoc_status, 0L) || !file.exists(docx_path)) {
  stop("Could not create the editable Word version of Table S1 with quarto pandoc.")
}

method_path <- file.path(OUT_DIR, "Table_S1_ash_deposition_mixed_layer_method.txt")
writeLines(
  c(
    paste0("Collection interval: ", COLLECTION_START, " to ", COLLECTION_END,
           " (", COLLECTION_DAYS, " inclusive days)"),
    paste0("Lake area: ", format(LAKE_AREA_M2, scientific = FALSE), " m2"),
    paste0("MLD: ", MLD_M, " m on ", MLD_DATE,
           "; density threshold = +0.1 kg m-3 from the near-surface reference"),
    "Total deposition (tons) = mean rate (ug m-2 d-1) * 16 d * 490e6 m2 / 1e12 ug ton-1",
    "Mixed-layer concentration (ug L-1) = mean rate * 16 d / (13 m * 1000 L m-3)"
  ),
  method_path,
  useBytes = TRUE
)

message("Mixed-layer ash-deposition Table S1 written to ", OUT_DIR)
