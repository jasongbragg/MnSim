# parlib/diags/recruitment/params.R
#
# Parameters for the recruitment diagnostic (diag_recruitment.R).
# Based on the working M. nodosa calibration; key recruitment-relevant
# parameters are annotated with calibration guidance below.
#
# RECRUITMENT MECHANISMS IN THIS PARAMETERISATION
# -----------------------------------------------
# 1. Beverton-Holt seed pool: expected_rec = R_max * n_flowering / (K_half + n_flowering)
#    R_max   = 500  -- hard ceiling on annual recruits
#    K_half  = 1000 -- n_flowering at which recruitment reaches half R_max
#    With these values, the "linear regime" (approximately linear in
#    n_flowering) extends to roughly 200-300 flowering adults; above that
#    the curve begins to saturate. The diagnostic will show whether your
#    equilibrium n_flowering is in the linear or saturated regime.
#
# 2. Canopy shade suppression: scales bevholt down by (1 - shade_suppression)
#    shade_dose_response$half_sat = 3000 -- N_canopy at 50% suppression
#    max_effect = 0.7 -- maximum suppression at infinite canopy density
#    hill = 4 -- steep sigmoid; suppression stays low until canopy is
#    substantial, then rises sharply. Check against your equilibrium
#    N_canopy in the diagnostic: if the equilibrium point falls well to
#    the left of half_sat, suppression is mild; well to the right, it is
#    near-maximal.
#
# 3. Rust fecundity suppression: currently OFF (max_effect = 0)
#    Each flowering adult's contribution is weighted by
#    (1 - dose_response(eff_susc, cfg)). Turn on once susceptibility data
#    arrive.

get_default_params <- function() {
  list(

    # ---- Run control -------------------------------------------------------
    N0      = 10000,
    n_years = 200,
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

    # ---- Recruitment (Beverton-Holt) ---------------------------------------
    # PRIMARY CALIBRATION TARGETS for this diagnostic.
    # R_max: set relative to measured or estimated annual seed production.
    # K_half: set relative to typical standing adult density. Lower K_half
    # means recruitment saturates earlier (each additional adult adds less
    # at equilibrium). The diagnostic marks the equilibrium n_flowering on
    # the BH curve so you can see directly whether you are in the linear
    # or saturated regime.
    R_max  = 4000,
    K_half = 3000,

    # ---- Canopy shading suppression of establishment -----------------------
    # SECONDARY CALIBRATION TARGET.
    # half_sat is in raw N_canopy units from canopy_density() -- set it
    # relative to the N_canopy your equilibrium population actually reaches
    # (the diagnostic will mark this). At the current half_sat = 3000,
    # a canopy of 3000 gives 50% suppression at max_effect = 0.7.
    shade_dose_response = list(
      max_effect = 1.0,
      form        = "saturating",
      half_sat    = 5000,
      hill         = 1
    ),

    # ---- Flowering ---------------------------------------------------------
    age_first_flower_mean = 5,
    age_first_flower_sd   = 1.5,

    # ---- Mating system -----------------------------------------------------
    selfing_rate = 0.10,

    # ---- Fire --------------------------------------------------------------
    fire_years            = c(40, 80, 120),
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
    rust_start_year    = 1,
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

    # ---- Rust fecundity suppression ----------------------------------------
    # OFF by default. When turned on, each flowering adult's contribution
    # to the expected recruit count is weighted by
    # (1 - dose_response(eff_susc, cfg)). Panel D of the diagnostic shows
    # this curve and labels it "off" while max_effect = 0.
    rust_recruit_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 2
    ),

    # ---- Genetics ----------------------------------------------------------
    resist_locus_effect = c(1.0),
    resist_dominance     = c(1.0),
    resist_freq0         = c(0.05)
  )
}
