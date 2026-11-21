# params.R
#
# Working parameter set of the default species (config/run_defaults.yaml,
# 'species': M. nodosa): the default for single-site runs, the Pegg et al.
# validation and the base for range-wide runs (which override fire and rust
# pressure per grid cell; runner/site_params.R).
#
# The values no longer live here (PLAN.md Phase 2; methodological notes T7, T8).
# They are read and validated by R/species_config.R from:
#   species/<sp>.yaml         biology and the rust anchor; every value cites
#                             its CSV row (data/m_nodosa_params.csv) or is a prior
#   config/landscape.yaml     shared fire defaults (single-site runs)
#   config/run_defaults.yaml  the default species; N0, n_years, seed, fire_years,
#                             rust_start_year
#   config/schema.yaml        field types, bounds and list order
# To change a value, edit the species YAML and its CSV row together; the loader
# refuses values that disagree. Assumptions behind the values:
# docs/methodological_notes.md.
#
# This file keeps the old interface: get_default_params() and RUST_ANCHOR.
# The configuration is loaded once, when this file is sourced, so any error in
# it shows up at that point. Sourcing from any working directory works (paths
# are found from this file's location); MNSIM_ROOT overrides the repo root.
#
# RUST_ANCHOR: the single place where the absolute level of rust pressure is
# set, shared by all species (config/landscape.yaml, rust_anchor; notes C2,
# T18): Yarringully NR, mean annual pathogen pressure 0.70 in the anchor cell
# (cell_lon, cell_lat). With mode sdm (the default; a species file may choose
# another) each grid cell's mean pressure is
#   mu * (SDM_cell / SDM_anchor_cell), capped at 'cap'
# (spatial/rust_pressure.R); 'conc' = alpha + beta of the annual Beta draw.

.mnsim_params_dir <- local({
  d <- NA_character_
  for (i in rev(seq_len(sys.nframe()))) {
    fn <- sys.function(i)
    fr <- sys.frame(i)
    if (identical(fn, base::source) && !is.null(fr$ofile)) {
      d <- if (isTRUE(fr$chdir)) getwd() else dirname(normalizePath(fr$ofile))
      break
    }
    if (identical(fn, base::sys.source) && !is.null(fr$file)) {
      d <- if (isTRUE(fr$chdir)) getwd() else dirname(normalizePath(fr$file))
      break
    }
  }
  if (is.na(d)) d <- normalizePath(".")   # not sourced from a file: assume plant_ibm/
  d
})

source(file.path(.mnsim_params_dir, "R", "species_config.R"), local = TRUE)

.mnsim_config <- species_config(root = Sys.getenv("MNSIM_ROOT", dirname(.mnsim_params_dir)))

RUST_ANCHOR <- .mnsim_config$rust_anchor

get_default_params <- function() .mnsim_config$params
