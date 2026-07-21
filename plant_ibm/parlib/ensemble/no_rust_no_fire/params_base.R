# parlib/ensemble/no_rust_no_fire/params_base.R
#
# Fixed (non-varying) parameters for the no-rust / no-fire ensemble.
# Used alongside param_ranges.R: this file defines every parameter;
# param_ranges.R specifies which subset to vary and over what range.
# vec_to_params() overwrites only the varying fields; everything else
# stays as defined here.
#
# PURPOSE: establish a demographically plausible baseline without
# disturbance, suitable for Morris screening and ABC calibration of
# intrinsic population structure (adult fraction, age profile,
# recruitment rate, stability).

get_default_params <- function() {
  list(

    # ---- Run control -------------------------------------------------------
    N0      = 5000L,
    n_years = 200L,
    seed    = 42L,

    # ---- Background mortality (Weibull, k=1: flat hazard) -----------------
    # k=1 paired with senescence_dose_response for the rising component.
    # weibull_lambda is a VARYING parameter in param_ranges.R.
    weibull_k      = 1.0,
    weibull_lambda = 60,

    # ---- Senescence hazard (Hill, rising with age) ------------------------
    # max_effect and half_sat VARY; hill is fixed (shape constraint).
    senescence_dose_response = list(
      max_effect = 0.30,
      form       = "sigmoid",
      half_sat   = 30,
      hill       = 4
    ),

    # ---- Juvenile decline (canopy-dependent height, Hill age-decay) -------
    # max_effect and age_half_sat VARY; half_sat (canopy) and hill fixed.
    juv_decline_dose_response = list(
      max_effect = 0.90,
      form       = "saturating",
      half_sat   = 7000,
      hill       = 4
    ),
    juv_decline_age_half_sat = 3,
    juv_decline_age_hill     = 8,

    # ---- Recruitment -------------------------------------------------------
    # R_max and K_half VARY.
    R_max  = 500,
    K_half = 1000,

    # ---- Canopy shading suppression of establishment ----------------------
    # max_effect and half_sat VARY; hill fixed.
    shade_dose_response = list(
      max_effect = 0.70,
      form       = "saturating",
      half_sat   = 3000,
      hill       = 4
    ),

    # ---- Flowering ---------------------------------------------------------
    # age_first_flower_mean VARIES; sd fixed.
    age_first_flower_mean = 5,
    age_first_flower_sd   = 1.5,

    # ---- Mating system -----------------------------------------------------
    selfing_rate = 0.10,

    # ---- Fire: OFF for this ensemble ---------------------------------------
    fire_years            = integer(0),
    fire_prob_annual      = 0,
    resprout_yrs_base     = 3,
    resprout_recovery     = "countdown",
    resprout_prob_recovery = 1 / 3,
    fire_kill_prob        = 0.35,
    fire_p_fimp           = 0,
    fire_kill_scalar      = 0.30,
    fire_kill_half_sat    = 5,
    fire_kill_hill        = 2,

    # ---- Rust: OFF for this ensemble ---------------------------------------
    rust_start_year    = Inf,
    rust_pressure      = 1.0,
    microclim_alpha    = 1,
    microclim_beta     = 1,
    annual_env_alpha   = 1,
    annual_env_beta    = 1,

    rust_dose_response = list(
      age_peak     = 0.60,
      age_floor    = 0.05,
      age_half_sat = 3,
      age_hill     = 8,
      resprout = list(
        max_effect = 0.50,
        form       = "saturating",
        half_sat   = 0.5,
        hill       = 1
      ),
      delay = list(
        max_effect = 2,
        form       = "linear",
        half_sat   = 1,
        hill       = 1
      )
    ),

    rust_flower_dose_response = list(
      max_effect = 0, form = "saturating", half_sat = 0.5, hill = 2
    ),
    rust_recruit_dose_response = list(
      max_effect = 0, form = "saturating", half_sat = 0.5, hill = 2
    ),

    # ---- Genetics ----------------------------------------------------------
    resist_locus_effect = c(1.0),
    resist_dominance    = c(1.0),
    resist_freq0        = c(0.0)
  )
}
