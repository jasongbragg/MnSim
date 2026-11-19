# tests/regression/scenarios.R
#
# Fixed-seed M. nodosa scenarios for the Phase 1 regression harness
# (PLAN.md Phase 1). Their outputs are frozen in ../tests/reference/ (repo
# root) by make_reference.R; test-regression.R re-runs them and compares.
#
# WHAT IS FROZEN
# Each reference file holds the exact params list each scenario ran with, as
# well as its output. That gives two separate checks:
#   code        re-running the FROZEN params must reproduce the frozen output
#               exactly (same R version, same seed).
#   calibration building the params afresh from today's params.R must give
#               the frozen params. If it fails, params.R has changed (e.g. a
#               ROADMAP recalibration), not the model code.
#
# TEST SCALE
# params.R settles at ~13,000 plants; CLAUDE.md caps test runs at 500.
# regression_base() multiplies the density-dependent parameters by
# REGRESSION_SCALE: R_max and K_half (counts of flowering plants) and the
# half_sat of juvenile decline and of shading (counts of canopy plants).
# Per-capita dynamics stay roughly the same and the population settles near
# 400. Every other value comes from params.R unchanged.
#
# Run with wd = plant_ibm/, after sourcing R/*.R, runner/params_utils.R,
# runner/iucn_a4.R and runner/pegg_validation.R (the test helper does this).

REGRESSION_SCALE <- 0.03
REGRESSION_N0    <- 400L

# Low-rust cell from spatial/data/mn_site_params.csv (Nov 2026): cell 7006,
# mean annual pressure 0.148 (anchor: 0.70). Values copied here so the
# scenario does not change when the site table is rebuilt.
REGRESSION_LOW_RUST <- list(cell_id = 7006L, alpha = 2.9527, beta = 17.0473)

# Working M. nodosa params from params.R, read into their own environment
# (the test helper defines a different get_default_params()), at test scale.
regression_base <- function(params_file = "params.R") {
  e <- new.env()
  sys.source(params_file, envir = e)
  p <- e$get_default_params()
  s <- REGRESSION_SCALE
  p$N0     <- REGRESSION_N0
  p$R_max  <- p$R_max  * s
  p$K_half <- p$K_half * s
  p$juv_decline_dose_response$half_sat <- p$juv_decline_dose_response$half_sat * s
  p$shade_dose_response$half_sat       <- p$shade_dose_response$half_sat * s
  p
}

# Each scenario: build(base) -> params; run(params) -> output (a list).
# run() must use nothing but its params argument, so that the frozen params
# alone reproduce the output.
REGRESSION_SCENARIOS <- list(

  no_fire = list(
    description = "No fire, no rust: background, senescence, juvenile decline, recruitment, genetics. 1850-1999.",
    build = function(p) {
      p$seed <- 101L; p$n_years <- 150L
      p$fire_prob_annual <- 0; p$fire_years <- integer(0)
      p$rust_start_year  <- Inf
      p
    },
    run = function(p) run_simulation(p, year0 = 1850L)
  ),

  fire = list(
    description = "Fixed fire years 1900/1930/1960/1990, stand mode (params.R fire settings), no rust. 1850-1999.",
    build = function(p) {
      p$seed <- 102L; p$n_years <- 150L
      p$fire_prob_annual <- 0; p$fire_years <- c(1900L, 1930L, 1960L, 1990L)
      p$rust_start_year  <- Inf
      p
    },
    run = function(p) run_simulation(p, year0 = 1850L)
  ),

  rust_high = list(
    description = "Rust from 2010 at the anchor pressure, Beta(14, 6); params.R stochastic fire (1/50 yr). 1900-2060.",
    build = function(p) {
      p$seed <- 103L; p$n_years <- 161L
      p$rust_start_year <- 2010L
      p
    },
    run = function(p) run_simulation(p, year0 = 1900L)
  ),

  rust_low = list(
    description = "As rust_high (same seed) at a low-rust cell (7006, mean pressure 0.148). 1900-2060.",
    build = function(p) {
      p$seed <- 103L; p$n_years <- 161L
      p$rust_start_year  <- 2010L
      p$annual_env_alpha <- REGRESSION_LOW_RUST$alpha
      p$annual_env_beta  <- REGRESSION_LOW_RUST$beta
      p
    },
    run = function(p) run_simulation(p, year0 = 1900L)
  ),

  pegg = list(
    description = "Pegg et al. 2025 scenario: spin-up 1910, rust 2010, fire burning every plant in 2020, survey Jan 2022.",
    build = function(p) { p$seed <- 105L; p },
    run = function(p) run_pegg_validation(p, site = NULL, plot = FALSE, verbose = FALSE)
  ),

  a4 = list(
    description = "Shortened IUCN A4 run with counterfactual: spin-up 1700, gen time 1800-1909, rust 2010, end 2126.",
    build = function(p) { p$seed <- 106L; p },
    # Shortened on purpose: its 100-yr spin-up is ~4.4 T_g, under the 5 T_g that
    # check_tg_windows() asks for, so that one warning is muffled (notes T11).
    run = function(p) withCallingHandlers(
      suppressMessages(run_iucn_a4(
        p, year_spinup = 1700L, year_gentime_start = 1800L, year_gentime_end = 1909L,
        year_rust = 2010L, year_assess = 2026L, year_end = 2126L,
        run_counterfactual = TRUE, plot = FALSE, verbose = FALSE)),
      mnsim_tg_window_warning = function(w) invokeRestart("muffleWarning"))
  )
)

regression_params <- function(name, base = regression_base()) {
  REGRESSION_SCENARIOS[[name]]$build(base)
}

regression_run <- function(name, params) {
  REGRESSION_SCENARIOS[[name]]$run(params)
}

regression_reference_path <- function(name, ref_dir) {
  file.path(ref_dir, paste0(name, ".rds"))
}
