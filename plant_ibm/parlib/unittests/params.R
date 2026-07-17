# parlib/unittests/params.R
#
# Conservative parameter baseline for the unit test suite.
# Sourced by tests/testthat/helper-fixtures.R INSTEAD of params.R.
#
# PURPOSE: the test suite needs get_default_params() to return a
# safe, predictable baseline with all calibrated mechanisms OFF.
# The main params.R now holds the working M. nodosa calibration
# (mechanisms active, site-specific values), which is correct for
# real runs but breaks tests that were written assuming a conservative
# "everything off" baseline -- e.g. the matrix-equivalence validation,
# which builds its survival schedule from weibull_hazard() alone and
# would be confounded by active senescence, juv_decline or shade terms
# inherited through make_geometric_params() / make_lefkovitch_params().
#
# SCHEMA: identical to params.R (same parameter names, same structure).
# Only VALUES differ -- conservative / off where the working calibration
# has active mechanisms.
#
# DO NOT calibrate this file. It is test infrastructure, not biology.
# Calibrated scenarios belong in parlib/diags/ or similar.

get_default_params <- function() {
  list(

    # ---- Founding population / run control --------------------------------
    N0      = 10000,
    n_years = 200,
    seed    = 42,

    # ---- Background mortality (Weibull, k=2: senescence via Weibull) ------
    # Original default -- k=2 gives a rising hazard that provides baseline
    # senescence for the matrix-equivalence tests without needing the Hill
    # senescence term active.
    weibull_k      = 2.0,
    weibull_lambda = 30,

    # ---- Senescence hazard (Hill function) -- OFF -------------------------
    senescence_dose_response = list(
      max_effect = 0,
      form        = "sigmoid",
      half_sat    = 30,
      hill         = 3
    ),

    # ---- Juvenile-decline hazard -- OFF -----------------------------------
    juv_decline_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 500,
      hill         = 1
    ),
    juv_decline_age_half_sat = 2,
    juv_decline_age_hill     = 1,

    # ---- Recruitment (Beverton-Holt) --------------------------------------
    R_max  = 4000,
    K_half = 5000,

    # ---- Canopy shading suppression of establishment -- OFF ---------------
    shade_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 5000,
      hill         = 1
    ),

    # ---- Flowering --------------------------------------------------------
    age_first_flower_mean = 5,
    age_first_flower_sd   = 1.5,

    # ---- Mating system ----------------------------------------------------
    selfing_rate = 0.10,

    # ---- Fire -------------------------------------------------------------
    fire_years            = c(40, 80, 120),
    fire_prob_annual      = 0,
    resprout_yrs_base     = 3,
    resprout_recovery     = "countdown",
    resprout_prob_recovery = 1 / 3,

    fire_kill_prob      = 0.35,
    fire_p_fimp         = 0,
    fire_kill_scalar    = 0.30,
    fire_kill_half_sat  = 5,
    fire_kill_hill      = 2,

    # ---- Myrtle rust -- OFF (rust_start_year = Inf) ----------------------
    rust_start_year    = Inf,
    rust_pressure      = 1.0,
    microclim_alpha    = 1,
    microclim_beta     = 1,
    annual_env_alpha   = 1,
    annual_env_beta    = 1,

    rust_dose_response = list(
      age_peak     = 0.20,
      age_floor    = 0,
      age_half_sat = 3,
      age_hill     = 2,
      resprout = list(
        max_effect = 0.30,
        form        = "saturating",
        half_sat    = 0.5,
        hill         = 1
      ),
      delay = list(
        max_effect = 2,
        form        = "linear",
        half_sat    = 1,
        hill         = 1
      )
    ),

    # Rust flower- and fecundity-suppression -- OFF
    rust_flower_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 2
    ),
    rust_recruit_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 2
    ),

    # ---- Genetic architecture ---------------------------------------------
    resist_locus_effect = c(1.0),
    resist_dominance     = c(1.0),
    resist_freq0         = c(0.05)
  )
}
