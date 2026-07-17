# params.R
#
# Single source of truth for all model parameters. Returns a named list;
# nothing here is hard-coded into the simulation functions, so scenarios
# are built by copying this list and overriding individual elements.
#
# Genetic architecture of myrtle rust resistance (resist_locus_effect,
# resist_dominance, resist_freq0) are equal-length vectors, one per
# locus. n_loci = length(params$resist_locus_effect) -- never set
# directly.
#
# DISEASE TRIANGLE (see rust_* params below)
# ------------------------------------------
# Rust effects are structured around the classic disease-triangle
# formulation: infection requires all three of host susceptibility,
# a permissive environment, and inoculum pressure. Any vertex at zero
# → no infection, regardless of the other two.
#
#   effective_susceptibility_i_t
#       = (1 - resist_score_i)   # host vertex (genetic, heritable)
#       × env_susceptibility_i   # environment vertex (microsite, fixed at birth)
#       × annual_env_t           # pathogen vertex  (site-wide, drawn each year)
#
# This product is naturally bounded in (0,1) -- no clamping needed.
# A resistant individual (resist_score → 1) is protected regardless
# of microsite or annual inoculum. A sheltered microsite (env_susc → 0)
# protects even a susceptible host in a bad year. A good year for rust
# (annual_env_t → 1) only harms susceptible individuals in exposed
# microsites.

get_default_params <- function() {
  list(

    # ---- Founding population / run control ---------------------------------
    N0      = 10000,
    n_years = 200,
    seed    = 42,

    # ---- Age-dependent background mortality (Weibull hazard) ---------------
    # k > 1: senescence via Weibull (hazard rises with age, unbounded,
    #         shape set jointly by k and lambda -- less interpretable to
    #         calibrate than the Hill-based alternative below).
    # k = 1: hazard = (k/lambda)*(age/lambda)^(k-1) reduces to a CONSTANT
    #         1/lambda at every age -- NOT zero. At lambda=30 (default)
    #         that's a flat ~3.3% annual hazard for every individual,
    #         every age. RECOMMENDED pairing with senescence_dose_response
    #         below if you want a deliberate non-senescent background rate
    #         alongside age-driven senescence; if you want
    #         senescence_dose_response to be the ONLY source of
    #         age-related mortality, also raise weibull_lambda to
    #         something very large (e.g. 1e9) -- k=1 alone does not zero
    #         this term out.
    # k < 1: front-loaded (hazard falls with age) -- rarely useful alone;
    #         juv_decline_dose_response is the better tool for front-loaded
    #         juvenile mortality, since it's also canopy-dependent.
    weibull_k      = 1.0,
    weibull_lambda = 60,       # characteristic lifespan, years

    # ---- Senescence hazard: rising Hill function of age --------------------
    # A more interpretable alternative to Weibull (k > 1): bounded,
    # asymptotic ceiling at max_effect, with a single clearly-named
    # parameter (half_sat) for "the age at which senescence reaches half
    # its eventual severity" -- e.g. half_sat = 30 means hazard is half
    # of max_effect by age 30. See senescence_hazard() in mortality.R for
    # the exact formula (reuses dose_response()'s "sigmoid" form with age
    # as the pressure axis).
    #
    # max_effect = 0 by default: off unless explicitly parameterised, so
    # the matrix-equivalence validation (which builds its survival
    # schedule from weibull_hazard only) is unaffected.
    #
    # TO ACTUALLY REPLACE Weibull senescence in a scenario: set
    # weibull_k = 1 AND weibull_lambda large (e.g. 1e9, NOT the default 30
    # -- see the weibull_k=1 note above, k=1 alone leaves a substantial
    # flat hazard behind) AND set senescence_dose_response$max_effect > 0
    # here. Both terms remain independently addressable (the union
    # formula in nominate_deaths() combines them like every other hazard
    # source), but running both weibull_k > 1 AND a non-zero
    # senescence_dose_response at once double-counts age-related
    # mortality for most use cases.
    senescence_dose_response = list(
      max_effect = 0.3,
      form        = "sigmoid",
      half_sat    = 30,
      hill         = 4
    ),

    # ---- Juvenile decline hazard (bathtub component) -----------------------
    # Adds a SEPARATE declining-with-age excess hazard for true juveniles
    # (age < AFR, not resprouting). Height scales with canopy density via
    # dose_response(N_canopy, juv_decline_dose_response). Age-decay uses
    # a Hill function: hill_weight(age, age_half_sat, age_hill) = same
    # shape as the rust age-decay but with independent parameters.
    #
    # Complementary to shade_dose_response in recruit(): that asks "does
    # a seed establish at all"; this asks "given a seedling, how fast
    # does it thin in its first couple of years under canopy."
    #
    # max_effect = 0 by default: off unless explicitly parameterised,
    # so the matrix-equivalence validation is unaffected.
    juv_decline_dose_response = list(
      max_effect = 0.99,
      form        = "saturating",
      half_sat    = 7000,
      hill         = 4
    ),
    juv_decline_age_half_sat = 3,  # age at which juv-decline hazard halves
    juv_decline_age_hill     = 8,  # Hill coeff (large = sharper early peak)

    # ---- Density-dependent mortality, with Hill age-weighting --------------
    K        = 5000,           # carrying capacity (total individuals)
    dd_alpha = 0.01,             # strength of crowding hazard at N = K
    #
    # dd_age_half_sat / dd_age_hill: Hill weight applied to the dd
    # scalar, so crowding hits the youngest/most-shaded individuals
    # hardest. hill_weight(age, half_sat, hill) = 1 at age 0 (full dd),
    # declines to 0.5 at age = half_sat, approaches 0 at old age.
    # Set dd_age_half_sat = Inf (default) to recover completely flat
    # density dependence -- equivalent to the old single-scalar dd_hazard.
    # These share the same functional shape as juv_decline but are
    # SEPARATELY parameterised: different biological claims, can be turned
    # on or off independently.
    dd_age_half_sat = Inf,
    dd_age_hill     = 8,

    # ---- Recruitment (Beverton-Holt) ---------------------------------------
    R_max  =  500,              # asymptotic max recruits/year
    K_half = 1000,              # flowering-adult count giving half R_max

    # ---- Canopy shading suppression of establishment -----------------------
    # Reduces how many seeds from bevholt() actually establish, as a
    # function of canopy_density() (see recruitment.R). max_effect = 0
    # by default (off). half_sat is in raw canopy-count units -- calibrate
    # to the actual canopy density your population reaches.
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
    fire_years       = c(40, 80, 120),
    fire_prob_annual = 0,          # if > 0 replaces fire_years (stochastic)
    resprout_yrs_base = 3,
    resprout_recovery      = "countdown",
    resprout_prob_recovery = 1 / 3,  # geometric mode only

    # Fire intensity model: two-stage, age-dependent
    # --------------------------------------------------
    # When fire_p_fimp = 0 (default), the LEGACY model is used:
    # every alive individual faces flat P(die) = fire_kill_prob,
    # survivors always resprout. This preserves backward compatibility
    # with existing scenarios and the Lefkovitch matrix-equivalence test.
    #
    # When fire_p_fimp > 0, the TWO-STAGE age-dependent model runs:
    #   Stage 1 -- impact: each plant impacted with P = fire_p_fimp.
    #              fire_p_fimp is an index of fire destructiveness (larger
    #              fires impact more plants and proportionally kill more).
    #   Stage 2 -- kill vs resprout given impact:
    #     p_kill(age) = c*p_fimp + (1 - c*p_fimp) * hill_weight(age, half_sat, hill)
    #     At age 0: p_kill = 1 (seedlings always killed if impacted).
    #     At age → ∞: p_kill → c*fire_p_fimp (ceiling for old plants).
    #     fire_kill_scalar (c): sets the ceiling for old-plant kill probability
    #              as a fraction of fire intensity -- hotter fires kill more
    #              old plants too, not just more total plants.
    #
    # The legacy fire_kill_prob is kept and used when fire_p_fimp = 0.
    fire_kill_prob      = 0.35,    # legacy: flat P(die | fire year)
    fire_p_fimp         = 0,       # 0 = legacy mode; > 0 = two-stage model
    fire_kill_scalar    = 0.30,    # c: old-plant kill ceiling = c * fire_p_fimp
    fire_kill_half_sat  = 5,       # age at which kill/resprout split is halfway
    fire_kill_hill      = 2,       # Hill coeff (large = sharper seedling selectivity)

    # ---- Myrtle rust: disease-triangle architecture ------------------------
    # See file header for the full effective_susceptibility formula.
    #
    # rust_start_year = Inf → rust never arrives (pre-rust calibration).
    # rust_pressure: site-level global scalar (0 = none, 1 = reference,
    # > 1 = more severe). Scales the age-decay peak and resprout bonus.
    #
    # Environment vertex of the disease triangle:
    #   env_susceptibility_i ~ Beta(microclim_alpha, microclim_beta)
    #     Drawn once per individual at birth; not heritable. Represents
    #     microsite variation in rust expression (aspect, local humidity,
    #     proximity to inoculum sources, etc). Beta(1,1) = Uniform(0,1)
    #     = maximum uncertainty / no microsite structure.
    #   annual_env_t ~ Beta(annual_env_alpha, annual_env_beta)
    #     One draw per timestep, shared by ALL individuals that year.
    #     Represents interannual variation in inoculum pressure / weather
    #     favourability for rust. Beta(1,1) = maximum year-to-year noise.
    #     Set annual_env_alpha >> annual_env_beta for consistently high
    #     pressure sites; reverse for low-pressure sites. Setting both
    #     large with alpha ≈ beta gives a near-constant site pressure.
    rust_start_year    = 1,
    rust_pressure      = 2.0,
    microclim_alpha    = 1,        # Beta shape: env_susceptibility per individual
    microclim_beta     = 1,
    annual_env_alpha   = 1,        # Beta shape: annual pathogen pressure
    annual_env_beta    = 1,

    # Rust hazard: continuous Hill-shaped age-decay with floor
    # ---------------------------------------------------------
    # Replaces the old discrete juvenile/adult gate in nominate_deaths().
    # Hazard = rust_pressure * age_peak * hill_weight(age, age_half_sat, age_hill)
    #          + age_floor
    # ... multiplied per-individual by effective_susceptibility.
    #
    # age_peak:    max hazard at age 0 (fully susceptible, at rust_pressure=1)
    # age_floor:   residual hazard at very old age. Default 0 (no adult mortality).
    #              Set > 0 to allow occasional rust-driven adult death -- rarely
    #              needed for M. nodosa but makes the model applicable to species
    #              where adult mortality from rust is common.
    # age_half_sat: age at which hazard halves toward floor
    # age_hill:    Hill coefficient (large = more gate-like; 1 = smooth decay)
    #
    # Resprout extra: additional hazard for resprouting individuals
    # (state-based, in addition to the age curve). A resprouting tree is
    # fire-stressed and its new growth is highly susceptible.
    rust_dose_response = list(
      age_peak     = 0.40,
      age_floor    = 0,
      age_half_sat = 3,
      age_hill     = 8,
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

    # Rust-suppressed flowering
    # -------------------------
    # Per-year Bernoulli: P(not flower | eligible) = dose_response(eff_susc, cfg).
    # Applied after the age/resprout gate: only eligible plants can have
    # flowering suppressed. max_effect = 0 by default (off).
    # NOTE: This is a per-year draw, so flowering can be recovered in a
    # good year (low annual_env_t). Revisit if data show persistent
    # multi-year flowering suppression once infection is established.
    rust_flower_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 2
    ),

    # Rust-suppressed fecundity (per-parent contribution to recruit pool)
    # -------------------------------------------------------------------
    # Each flowering adult's contribution to the expected recruit count
    # is discounted by dose_response(eff_susc, cfg). max_effect = 0 (off).
    rust_recruit_dose_response = list(
      max_effect = 0,
      form        = "saturating",
      half_sat    = 0.5,
      hill         = 2
    ),

    # ---- Genetic architecture of rust resistance ---------------------------
    resist_locus_effect = c(1.0),
    resist_dominance     = c(1.0),
    resist_freq0          = c(0.05)
  )
}
