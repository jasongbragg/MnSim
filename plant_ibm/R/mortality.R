# R/mortality.R
#
# Mortality is fully vectorised: every alive individual gets a single
# p_death computed in one pass, then one rbinom() draw nominates deaths.
# There is no per-individual for-loop.
#
# COMBINATION RULE
# ----------------
# All hazard components are combined via the proper probability-union
# formula rather than the additive approximation + clip:
#
#   p_death = 1 - (1-p_age) * (1-p_dd) * (1-p_jd) * (1-p_rust)
#
# This is the exact probability of "at least one of these independent
# causes kills this individual." It never exceeds 1 by construction, so
# the old pmin(...,1) cap is no longer needed (retained as a defensive
# assertion only). Unlike the additive sum, it remains correct even when
# individual hazard terms are large.

#' Age-dependent hazard, Weibull form: (k/lambda) * (age/lambda)^(k-1)
weibull_hazard <- function(age, params) {
  k        <- params$weibull_k
  lambda   <- params$weibull_lambda
  age_safe <- pmax(age, 0.01)   # floor: avoids 0^(k-1) blowing up for k<1
  (k / lambda) * (age_safe / lambda) ^ (k - 1)
}

#' Density-dependent hazard scalar (flat across ages; age-weighting via
#' dd_age_weight() is applied separately in nominate_deaths()).
dd_hazard <- function(N, params) {
  params$dd_alpha * (N / params$K)
}

#' Declining Hill-function weight: 1 at x = 0, 0.5 at x = half_sat,
#' approaching 0 as x → ∞. Used for both the juvenile-decline age-decay
#' and the density-dependence age-weighting (with independent parameters
#' for each -- see params.R for why they're kept separate).
#'
#' Setting half_sat = Inf recovers a flat weight of 1 everywhere (the
#' default for dd_age_half_sat, which preserves old flat dd behaviour).
#' Large hill values produce a sharper, more gate-like transition.
hill_weight <- function(x, half_sat, hill) {
  x <- pmax(x, 0)
  1 / (1 + (x / half_sat) ^ hill)
}

#' Age-weighting for the density-dependent hazard. Multiplied onto
#' dd_hazard() scalar in nominate_deaths() so crowding hits younger,
#' shorter, more-shaded individuals disproportionately.
#'
#' Default parameters (dd_age_half_sat = Inf) give weight = 1 at all
#' ages -- identical to the old unweighted dd_hazard.
dd_age_weight <- function(age, params) {
  hill_weight(age, params$dd_age_half_sat, params$dd_age_hill)
}

#' Generic dose-response transform, mapping a non-negative pressure
#' value to a realised effect size. `cfg` needs: max_effect, form,
#' half_sat, hill. Forms: "linear", "power", "saturating", "sigmoid",
#' "threshold".
dose_response <- function(pressure, cfg) {
  pressure <- pmax(pressure, 0)
  max_eff  <- cfg$max_effect
  hs       <- cfg$half_sat
  hill     <- cfg$hill

  out <- switch(cfg$form,
    linear     = max_eff * pressure,
    power      = max_eff * pressure ^ hill,
    saturating = max_eff * pressure / (pressure + hs),
    sigmoid    = max_eff * pressure ^ hill / (pressure ^ hill + hs ^ hill),
    threshold  = ifelse(pressure >= hs, max_eff, 0),
    stop("Unknown dose_response form: '", cfg$form, "'. Must be one of: ",
         "linear, power, saturating, sigmoid, threshold.")
  )
  pmax(out, 0)
}

#' Returns the resprout-delay modifier rust imposes this year.
#' Returns delay_extra = 0 before rust_start_year (resistance locus is
#' then selectively neutral -- freely drifting).
#'
#' Note: the age-decay and resprout-bonus rust mortality components are
#' now computed per-individual directly in nominate_deaths() via the
#' disease-triangle effective_susceptibility, NOT via this function.
#' rust_modifiers() is kept only for the resprout delay used in
#' apply_fire() and by the diagnostic functions in
#' diag_rust_pressure.R.
rust_modifiers <- function(params, t) {
  if (t < params$rust_start_year) {
    return(list(delay_extra = 0))
  }
  list(
    delay_extra = dose_response(params$rust_pressure, params$rust_dose_response$delay)
  )
}

#' Juvenile-decline excess hazard: declining Hill in age, with height
#' scaled by canopy density (dose_response on N_canopy). Applies only to
#' true juveniles (age < AFR, !resprout) in nominate_deaths().
#'
#' hill_weight replaces the old exp(-age/tau) decay: same biological
#' intent (rapid early decline) but parameterised with half_sat and hill,
#' which subsumes the exponential as a special case (hill = 1 is similar
#' in shape, large hill recovers a near-gate). juv_decline_tau is removed.
juv_decline_hazard <- function(age, N_canopy, params) {
  height <- dose_response(N_canopy, params$juv_decline_dose_response)
  height * hill_weight(age, params$juv_decline_age_half_sat,
                        params$juv_decline_age_hill)
}

#' Senescence hazard: a RISING Hill function of age, replacing Weibull
#' (k > 1) as a more interpretable alternative. Reuses dose_response()'s
#' "sigmoid" form directly, with age as the pressure axis:
#'
#'   hazard(age) = max_effect * age^hill / (age^hill + half_sat^hill)
#'
#' hazard(0) = 0. hazard(half_sat) = max_effect / 2 EXACTLY -- half_sat
#' is the age at which senescence reaches half its eventual ceiling, a
#' single interpretable number (e.g. half_sat = 30 means "by age 30,
#' senescence-driven mortality is halfway to its max"). hazard -> max_effect
#' as age -> Inf, a hard asymptotic ceiling (unlike Weibull, which is
#' unbounded for any k > 1). hill controls steepness of the transition --
#' larger hill = senescence stays negligible for longer, then rises more
#' sharply around half_sat.
#'
#' max_effect = 0 by default: off unless explicitly parameterised, so
#' existing scenarios/tests (including matrix-equivalence, which builds
#' its survival schedule from weibull_hazard only) are unaffected.
#'
#' RECOMMENDED PAIRING: this is meant to REPLACE weibull_k > 1 senescence,
#' not stack with it. To switch a scenario over: set weibull_k = 1 (flat
#' background hazard only) and senescence_dose_response$max_effect > 0.
#' Both terms remain independently addressable in the union formula in
#' nominate_deaths() if you have a specific reason to want both shapes
#' active simultaneously, but that would double-count age-related
#' mortality for most use cases.
senescence_hazard <- function(age, params) {
  dose_response(age, params$senescence_dose_response)
}

#' Compute p_death for every alive individual and draw deaths.
#' Returns a vector of ROW INDICES (into pop) nominated to die this year.
#'
#' Components, combined via the probability-union formula:
#'   p_age  = Weibull background hazard (flat when weibull_k = 1)
#'   p_sen  = senescence hazard, rising Hill function of age (off by
#'            default; see senescence_hazard() -- intended to REPLACE
#'            weibull_k > 1, paired with weibull_k = 1 above)
#'   p_dd   = density hazard × Hill age-weight (youngest bear most)
#'   p_jd   = juvenile-decline hazard (canopy-dependent height, Hill age-decay;
#'            juveniles only)
#'   p_rust = rust hazard via disease triangle (continuous Hill age-decay
#'            with floor; resprout individuals get an additional state bonus)
#'
#' annual_env_t is the pathogen vertex of the disease triangle (one Beta
#' draw per year, shared by all individuals, passed in from simulate.R).
#' Default 0 → no rust effect (used by tests that call nominate_deaths()
#' directly without simulating an annual loop).
nominate_deaths <- function(pop, t, params, annual_env_t = 0) {
  alive_idx <- which(pop$alive)
  n_alive   <- length(alive_idx)
  if (n_alive == 0) return(integer(0))

  N        <- n_alive
  N_canopy <- canopy_density(pop)

  age_i      <- pop$age[alive_idx]
  resprout_i <- pop$resprout[alive_idx]
  juv_i      <- (age_i < pop$age_first_flower[alive_idx]) & !resprout_i
  resist_i   <- pop$resist_score[alive_idx]
  env_susc_i <- pop$env_susceptibility[alive_idx]

  # --- Background Weibull hazard ------------------------------------------
  # With weibull_k = 1 (recommended pairing for senescence_hazard below),
  # this is a flat constant background, not a senescence term.
  p_age <- weibull_hazard(age_i, params)

  # --- Senescence hazard: rising Hill function of age ----------------------
  # Off by default (max_effect = 0). See senescence_hazard() docstring --
  # intended to replace weibull_k > 1, not stack with it.
  p_sen <- senescence_hazard(age_i, params)

  # --- Density-dependent hazard, Hill-weighted by age ---------------------
  p_dd <- dd_hazard(N, params) * dd_age_weight(age_i, params)

  # --- Juvenile-decline hazard (juveniles only) ---------------------------
  p_jd <- numeric(n_alive)
  if (any(juv_i)) {
    p_jd[juv_i] <- juv_decline_hazard(age_i[juv_i], N_canopy, params)
  }

  # --- Rust hazard via disease triangle ------------------------------------
  # effective_susceptibility = host × environment × pathogen
  #   host        = (1 - resist_score)         genetic, heritable
  #   environment = env_susceptibility          microsite, fixed at birth
  #   pathogen    = annual_env_t               site-wide draw each year
  p_rust <- numeric(n_alive)
  if (t >= params$rust_start_year && annual_env_t > 0) {
    eff_susc <- (1 - resist_i) * env_susc_i * annual_env_t

    # Continuous Hill age-decay with floor, scaled by rust_pressure
    rust_age_decay <- hill_weight(age_i,
                                   params$rust_dose_response$age_half_sat,
                                   params$rust_dose_response$age_hill)
    peak <- params$rust_pressure * params$rust_dose_response$age_peak
    rust_age_h <- params$rust_dose_response$age_floor +
                  (peak - params$rust_dose_response$age_floor) * rust_age_decay

    # State-based resprout bonus (post-fire stressed tissue)
    rust_resprout_h <- ifelse(resprout_i,
      dose_response(params$rust_pressure, params$rust_dose_response$resprout), 0)

    p_rust <- eff_susc * (rust_age_h + rust_resprout_h)
  }

  # --- Probability union ---------------------------------------------------
  p_death <- 1 - (1 - p_age) * (1 - p_sen) * (1 - p_dd) * (1 - p_jd) * (1 - p_rust)
  p_death <- pmin(p_death, 1)   # defensive -- union formula should not need this

  dead_draw <- rbinom(n_alive, 1, p_death) == 1
  alive_idx[dead_draw]
}
