# parlib/diags/rust_age/params.R
#
# Parameters for diag_rust_age.R.
# Focus: rust age-hazard shape and the resprout bonus.
# Fire and rust are both active so the Panel C trajectories are meaningful.
# The key parameters for this diagnostic are annotated.

get_default_params <- function() {
  list(

    N0      = 10000,
    n_years = 200,
    seed    = 42,

    weibull_k      = 1.0,
    weibull_lambda = 60,

    senescence_dose_response = list(
      max_effect = 0.2, form = "sigmoid", half_sat = 40, hill = 4
    ),

    juv_decline_dose_response = list(
      max_effect = 0.5, form = "saturating", half_sat = 5000, hill = 2
    ),
    juv_decline_age_half_sat = 2,
    juv_decline_age_hill     = 2,

    R_max  = 500,
    K_half = 1000,

    shade_dose_response = list(
      max_effect = 0, form = "saturating", half_sat = 2500, hill = 4
    ),

    age_first_flower_mean = 5,
    age_first_flower_sd   = 1.5,
    selfing_rate          = 0.10,

    fire_years            = c(40, 80, 120),
    fire_prob_annual      = 0,
    resprout_yrs_base     = 3,     # Panel C window length -- change to explore
    resprout_recovery     = "countdown",
    resprout_prob_recovery = 1 / 3,
    fire_kill_prob        = 0.35,
    fire_p_fimp           = 0.8,
    fire_kill_scalar      = 0.30,
    fire_kill_half_sat    = 5,
    fire_kill_hill        = 2,

    # ---- Rust: PRIMARY CALIBRATION TARGETS for this diagnostic ------------
    rust_start_year = 1,
    rust_pressure   = 1.0,

    rust_dose_response = list(
      # age_peak: maximum rust hazard at age 0 for fully susceptible plant
      # at rust_pressure = 1. Panel A peak.
      age_peak     = 0.025,

      # age_floor: residual hazard at very old age. Panel A asymptote.
      # Set > 0 to allow occasional rust-driven adult death.
      age_floor    = 0.025,

      # age_half_sat: age at which hazard halves toward floor.
      # The Panel A inflection point.
      age_half_sat = 2,

      # age_hill: steepness of the Hill decay. Large values = sharper
      # transition; effectively a gate at age_half_sat.
      age_hill     = 4,

      # resprout: the POST-FIRE rust bonus. max_effect controls the
      # magnitude of the step-up in Panel C. Calibrate against the
      # post-fire survival data from the new literature.
      resprout = list(
        max_effect = 0.5,
        form        = "saturating",
        half_sat    = 0.5,
        hill         = 1
      ),

      delay = list(
        max_effect = 2, form = "linear", half_sat = 1, hill = 1
      )
    ),

    microclim_alpha = 20,   # env_susc: mean=0.80, tight
    microclim_beta  = 5,
    annual_env_alpha = 14,  # annual pressure: mean=0.70
    annual_env_beta  = 6,

    rust_flower_dose_response = list(
      max_effect = 0, form = "saturating", half_sat = 0.5, hill = 2
    ),
    rust_recruit_dose_response = list(
      max_effect = 0, form = "saturating", half_sat = 0.5, hill = 2
    ),

    resist_locus_effect = c(1.0),
    resist_dominance    = c(1.0),
    resist_freq0        = c(0.0)
  )
}
