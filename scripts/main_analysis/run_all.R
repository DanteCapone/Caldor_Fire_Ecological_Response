# Rebuild final manuscript Figures 2-5 from the repository root.

prerequisite_scripts <- file.path(
  "scripts", "main_analysis", "supplemental_components",
  "ctd_profiles_stability_pkl.R"
)

lake_tools_plot_script <- file.path(
  "scripts", "main_analysis", "supplemental_components",
  "lake_tools_stability_plots.R"
)
lake_tools_python_script <- file.path(
  "scripts", "main_analysis", "calculate_lake_tools_stability.py"
)

main_figure_scripts <- file.path(
  "scripts", "main_analysis",
  c("figure_2_2021_final.R", "figure_3_v2.R", "figure_4.R", "figure_5.R")
)

figure4_variants_script <- file.path(
  "scripts", "figure_4_phytoplankton_by_depth.R"
)

missing_scripts <- c(
  prerequisite_scripts, lake_tools_plot_script, lake_tools_python_script,
  main_figure_scripts, figure4_variants_script
)
missing_scripts <- missing_scripts[!file.exists(missing_scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing main-analysis scripts: ", paste(missing_scripts, collapse = ", "))
}

for (script in prerequisite_scripts) {
  message("Prerequisite: ", script)
  sys.source(script, envir = new.env(parent = globalenv()))
  if ("package:data.table" %in% search()) {
    detach("package:data.table", unload = FALSE, character.only = TRUE)
  }
}

python_candidates <- unique(c(
  Sys.getenv("CTD_PYTHON", unset = NA_character_),
  file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "geomap", "python.exe"),
  Sys.which("python")
))
python_candidates <- python_candidates[
  !is.na(python_candidates) & nzchar(python_candidates) & file.exists(python_candidates)
]
python_ok <- vapply(python_candidates, function(candidate) {
  identical(
    suppressWarnings(system2(candidate, c("-c", shQuote("import pandas, numpy")),
                             stdout = FALSE, stderr = FALSE)),
    0L
  )
}, logical(1))
if (!any(python_ok)) stop("No Python environment with pandas and numpy was found.")
python <- python_candidates[which(python_ok)[1]]
message("Lake-Tools metrics: ", lake_tools_python_script)
status <- system2(
  python,
  c(shQuote(lake_tools_python_script), "--project-root",
    shQuote(normalizePath(".", winslash = "/")))
)
if (!identical(status, 0L)) stop("Lake-Tools calculation failed with status ", status)
sys.source(lake_tools_plot_script, envir = new.env(parent = globalenv()))

for (script in main_figure_scripts) {
  message("\n=== Running ", script, " ===")
  sys.source(script, envir = new.env(parent = globalenv()))
  if ("package:data.table" %in% search()) {
    detach("package:data.table", unload = FALSE, character.only = TRUE)
  }
}

message("\n=== Running ", figure4_variants_script, " ===")
sys.source(figure4_variants_script, envir = new.env(parent = globalenv()))

message("\nCompleted final Figures 2-5.")
