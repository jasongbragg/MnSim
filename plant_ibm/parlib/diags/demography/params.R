# parlib/diags/demography/params.R
#
# Base parameterisation for the demography diagnostic (diag_demography.R).
# Starting point: the working M. nodosa calibration from plant_ibm/params.R.
#
# The diagnostic runs three phases automatically (no need to edit here):
#   1. 500-year fire-free spin-up (rust off regardless of rust_start_year)
#   2. Post-fire rust-ON scenario  (annual_env_t ≈ 1.0, max constant pressure)
#   3. Post-fire rust-OFF scenario (rust_start_year forced to Inf)
#
# Parameters to consider adjusting here:
#   rust_dose_response  -- determines how hard rust hits each age class
#   resist_freq0        -- allele frequency of resistance (0 = fully susceptible,
#                          shows absolute worst-case rust impact)
#   fire_p_fimp / fire_kill_scalar / fire_kill_half_sat -- fire severity
#   resprout_yrs_base   -- recovery window length

get_default_params <- function() {
  list(

    N0      = 10000,
    n_years = 200,    # overridden internally by diag_demography()
    seed    = 42,

    # ---- Background mortality ----------------------------------------------
    weibull_k      = 1.0,
    weibull_lambda = 60,

    senescence_dose_response = list(
      max_effect = 0.3,
      form        = "sigmoid",
      half_sat    = 30,
      hill         = 4
    ),

    juv_decline_dose_response = list(
      max_effect = 0.99,
      form        = "saturating",
      half_sat    = 7000,
      hill         = 4
    ),
    juv_decline_age_half_sat = 3,
    juv_decline_age_hill     = 8,

    # ---- Recruitment -------------------------------------------------------
    R_max  = 500,
    K_half = 1000,

    shade_dose_response = list(
      max_effect = 0.7,
      form        = "saturating",
      half_sat    = 3000,
      hill         = 4
    ),

    # ---- Flowering ---------------------------------------------------------
    age_first_flower_mean = 5,
    age_first_flower_sd   = 1.5,

    # ---- Mating system -----------------------------------------------------
    selfing_rate = 0.10,

    # ---- Fire --------------------------------------------------------------
    # fire_years / fire_p_fimp overridden internally for the single fire event.
    # Edit fire_kill_scalar, fire_kill_half_sat, fire_kill_hill to tune severity.
    fire_years            = c(40, 80, 120),  # ignored during diagnostic runs
    fire_prob_annual      = 0,
    resprout_yrs_base     = 3,
    resprout_recovery     = "countdown",
    resprout_prob_recovery = 1 / 3,
    fire_kill_prob        = 0.35,
    fire_p_fimp           = 0.5,
    fire_kill_scalar      = 0.30,
    fire_kill_half_sat    = 5,
    fire_kill_hill        = 2,

    # ---- Myrtle rust -------------------------------------------------------
    # rust_start_year and annual_env_alpha/beta are overridden by the
    # diagnostic for each scenario. Edit rust_dose_response to calibrate
    # how rust affects each age class and the resprout bonus.
    rust_start_year    = 1,
    rust_pressure      = 1.0,
    microclim_alpha    = 1,
    microclim_beta     = 1,
    annual_env_alpha   = 1,   # overridden to ~1e4 in rust-ON scenario
    annual_env_beta    = 1,   # overridden to 1   in rust-ON scenario

    rust_dose_response = list(
      age_peak     = 0.60,
      age_floor    = 0.05,
      age_half_sat = 3,
      age_hill     = 8,
      resprout = list(
        max_effect = 0.50,
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

    # ---- Genetic architecture ----------------------------------------------
    # resist_freq0 = c(0): fully susceptible population -- shows the
    # absolute worst-case rust impact (upper bound on demographic damage).
    # Change to c(0.05) or higher to add resistance buffering.
    resist_locus_effect = c(1.0),
    resist_dominance     = c(1.0),
    resist_freq0         = c(0)
  )
}
