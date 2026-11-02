# params.R
#
# Working M. nodosa calibration: the default parameter set for single-site
# runs, the Pegg et al. validation and the base for range-wide runs (which
# override fire and rust pressure per grid cell; see runner/site_params.R).
# Assumptions behind these values: docs/methodological_notes.md.
#
# Parameters to consider adjusting here:
#   rust_dose_response  -- determines how hard rust hits each age class
#   resist_freq0        -- allele frequency of resistance (0 = fully susceptible,
#                          shows absolute worst-case rust impact)
#   fire_p_fimp / fire_kill_scalar / fire_kill_half_sat -- fire severity
#   resprout_yrs_base   -- recovery window length

# ---- Rust-pressure anchor ---------------------------------------------------
# The single place where the absolute level of rust pressure is set.
# Yarringully NR (Pegg et al. 2025 site): field incidence gives a mean annual
# pathogen pressure of 0.70. Range-wide, each grid cell's mean pressure is
#   mu * (SDM_cell / SDM_anchor_cell), capped at 'cap'
# where the anchor cell is the occupied 0.1-degree cell nearest (lon, lat),
# so the anchor cell always gets exactly 'mu' whatever SDM is used
# (spatial/rust_pressure.R). 'conc' = alpha + beta of the annual Beta draw
# (year-to-year variability; not yet constrained by data).
RUST_ANCHOR <- list(
  site = "Yarringully NR",
  lon  = 153.21967,
  lat  = -29.09357,
  mu   = 0.70,
  conc = 20,
  cap  = 0.95
)

get_default_params <- function() {
  list(

    N0      = 10000,
    n_years = 200,    # overridden internally by diag_demography()
    seed    = 42,

    # ---- Background mortality ----------------------------------------------
    weibull_k      = 1.0,
    weibull_lambda = 60,

    senescence_dose_response = list(
      max_effect = 0.1,
      form        = "sigmoid",
      half_sat    = 40,
      hill         = 8
    ),

    juv_decline_dose_response = list(
      max_effect = 0.5,
      form        = "saturating",
      half_sat    = 5000,
      hill         = 2
    ),
    juv_decline_age_half_sat = 2,
    juv_decline_age_hill     = 2,

    # ---- Recruitment -------------------------------------------------------
    R_max  = 1000,
    K_half = 4000,

    shade_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 1000,
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
    fire_years            = integer(0),  # ignored during diagnostic runs
    fire_prob_annual      = 1/50,
    resprout_yrs_base     = 3,
    resprout_recovery     = "countdown",
    resprout_prob_recovery = 1 / 3,
    fire_kill_prob        = 0.35,
    fire_p_fimp           = 0.6,
    fire_kill_scalar      = 0.30,
    fire_kill_half_sat    = 2,
    fire_kill_hill        = 2,

    # ---- Myrtle rust -------------------------------------------------------
    # rust_start_year and annual_env_alpha/beta are overridden by the
    # diagnostic for each scenario. Edit rust_dose_response to calibrate
    # how rust affects each age class and the resprout bonus.
    rust_start_year    = Inf,
    rust_pressure      = 1,
    microclim_alpha    = 20,
    microclim_beta     = 5,
    # Annual pathogen draw at the anchor site: Beta(mu*conc, (1-mu)*conc),
    # i.e. Beta(14, 6). Range-wide runs replace these per cell.
    annual_env_alpha   = RUST_ANCHOR$mu * RUST_ANCHOR$conc,
    annual_env_beta    = (1 - RUST_ANCHOR$mu) * RUST_ANCHOR$conc,

    rust_dose_response = list(
      age_peak     = 0.3,
      age_floor    = 0.025,
      age_half_sat = 2,
      age_hill     = 4,
      resprout = list(
        # 0.53: Pegg et al. median at the Yarringully cell (Nov 2026 sweep;
        # docs/methodological_notes.md V4). Plausible range 0.32-0.78.
        max_effect = 0.53,
        form        = "saturating",
        half_sat    = 0.2,
        hill         = 1
      ),
      delay = list(
        max_effect = 1,
        form        = "linear",
        half_sat    = 1,
        hill         = 1
      )
    ),

    rust_flower_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 1
    ),

    rust_recruit_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 1
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
