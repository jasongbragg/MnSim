# spatial/fire/R/regime.R -- from cell-year tables to the per-cell fire regime
# that 02 writes into mn_site_params.csv (docs/methodological_notes.md D2, D5).
#
# Needs: fire_config.R, R/estimate.R; cell-year tables in FIRE_DATA_DIR
# (spatial/fire/02_fire_cell_years.R); cells with cell_id and state
# (load_fire_cells()).

read_cy <- function(name, dir = FIRE_DATA_DIR)
  utils::read.csv(file.path(dir, paste0("cell_years_", name, ".csv")))

# Cells with state and land area, cached (geometry work takes a minute).
load_fire_cells <- function(cache = file.path(FIRE_DATA_DIR, "cells.rds"), grid_csv = MN_GRID_CSV) {
  if (file.exists(cache)) return(readRDS(cache))
  cells <- load_cells(grid_csv)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cells, cache)
  cells
}

# One scenario (fire_scenario() in fire_config.R): the estimate_fire_regime()
# table for every cell, plus state and record length.
run_fire_scenario <- function(sc, cells, dir = FIRE_DATA_DIR) {
  nsw <- read_cy(cy_table("NSW", sc$fire, sc$area_base), dir)
  qld <- read_cy(cy_table("Qld", sc$fire, sc$area_base), dir)
  # each state's record applies to that state's cells only
  nsw <- nsw[nsw$year %in% sc$nsw & nsw$cell_id %in% cells$cell_id[cells$state == "NSW"], ]
  qld <- qld[qld$year %in% sc$qld & qld$cell_id %in% cells$cell_id[cells$state == "Qld"], ]
  cy  <- rbind(nsw[, c("cell_id", "year", "frac_burned")], qld[, c("cell_id", "year", "frac_burned")])
  T_obs <- stats::setNames(ifelse(cells$state == "NSW", length(sc$nsw), length(sc$qld)), cells$cell_id)
  grp <- if (sc$prior_by == "state") cells$state else NULL
  reg <- estimate_fire_regime(cy, cells, T_obs, prior = sc$prior, prior_group = grp,
                              fire_min = sc$fire_min)
  reg$state <- cells$state[match(reg$cell_id, cells$cell_id)]
  attr(reg, "cy") <- cy
  reg
}

# Per-cell table for mn_site_params.csv:
#   wildfire (the default regime) under the standard names that the model
#   reads (fire_prob_annual, fire_p_fimp, sd_frac_burned, ...), plus the
#   record behind it (fire_state, fire_T_obs, fire_record);
#   all fire and planned burns with suffixes _all and _planned.
site_fire_regime <- function(cells, regimes = FIRE_REGIMES, dir = FIRE_DATA_DIR) {
  regs <- lapply(regimes, run_fire_scenario, cells = cells, dir = dir)
  w <- regs$wild
  out <- data.frame(
    cell_id = w$cell_id, fire_state = w$state, fire_T_obs = w$T_obs,
    fire_record = ifelse(w$state == "NSW",
      sprintf("NSW NPWS wildfire %d-%d", min(regimes$wild$nsw), max(regimes$wild$nsw)),
      sprintf("Qld fire scars minus QPWS planned burns %d-%d", min(regimes$wild$qld), max(regimes$wild$qld))),
    n_wildfire_years_all = w$n_wildfire_years_all,
    n_wildfire_years_sig = w$n_wildfire_years_sig,
    fire_prob_annual = w$fire_prob_annual,
    mean_frac_burned = w$mean_frac_burned,
    sd_frac_burned   = w$sd_frac_burned,
    var_frac_burned  = w$var_frac_burned,
    fire_p_fimp      = w$fire_p_fimp,
    FRI_years        = w$FRI_years)
  for (nm in setdiff(names(regs), "wild")) {
    r <- regs[[nm]][match(out$cell_id, regs[[nm]]$cell_id), ]
    out[[paste0("n_fire_years_", nm)]]     <- r$n_wildfire_years_sig
    out[[paste0("fire_prob_annual_", nm)]] <- r$fire_prob_annual
    out[[paste0("mean_frac_burned_", nm)]] <- r$mean_frac_burned
    out[[paste0("sd_frac_burned_", nm)]]   <- r$sd_frac_burned
    out[[paste0("fire_p_fimp_", nm)]]      <- r$fire_p_fimp
  }
  attr(out, "priors") <- do.call(rbind, lapply(names(regs), function(nm) {
    p <- attr(regs[[nm]], "priors"); p$regime <- nm; p }))
  out
}

# Columns site_fire_regime() adds to the site table (02 drops these before merging).
SITE_FIRE_COLUMNS <- function(regimes = FIRE_REGIMES)
  c("fire_state", "fire_T_obs", "fire_record", "n_wildfire_years_all", "n_wildfire_years_sig",
    "n_wildfire_events", "fire_prob_annual", "mean_frac_burned", "sd_frac_burned",
    "var_frac_burned", "fire_p_fimp", "FRI_years",
    unlist(lapply(setdiff(names(regimes), "wild"), function(nm)
      paste0(c("n_fire_years_", "fire_prob_annual_", "mean_frac_burned_", "sd_frac_burned_",
               "fire_p_fimp_"), nm))))
