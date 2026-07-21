# parlib/ensemble/no_rust_no_fire/param_ranges.R
#
# Defines which parameters vary and their plausible ranges for the
# no-rust / no-fire Morris screening and ABC ensemble.
#
# STRUCTURE
# ---------
# keys:  character vector of __ names matching the nested params list
# mins:  numeric vector of lower bounds (same order as keys)
# maxs:  numeric vector of upper bounds (same order as keys)
# Sigma: NULL for independent draws; an n_keys × n_keys correlation
#        matrix for correlated draws via the Gaussian copula
#        (used with ensemble.R generate_params_matrix(method="mvn_copula"))
#
# RANGE RATIONALE
# ---------------
# weibull_lambda:          [20, 120]
#   Flat background hazard = 1/lambda. Range spans 0.8% to 5% annual
#   background mortality -- covers from very long-lived to moderate.
#
# senescence_dose_response__max_effect: [0, 0.60]
#   Maximum senescence hazard at old age. Upper bound at 0.60 means even
#   very old individuals have at most 60% annual mortality from senescence
#   alone (before background or juvenile terms).
#
# senescence_dose_response__half_sat:   [15, 60]
#   Age at which senescence reaches half its ceiling. Range covers
#   "starts early at 15yr" to "barely noticeable until 60yr".
#
# juv_decline_dose_response__max_effect: [0.30, 0.99]
#   Height of the juvenile-decline hazard at full canopy. Lower bound
#   at 0.30 (still substantial early mortality); upper at 0.99 (near-
#   certain death for age-0 seedlings under full canopy).
#
# juv_decline_age_half_sat:             [1, 8]
#   Age at which juvenile-decline hazard halves. Range covers "falls off
#   in 1 year" to "still declining at age 8".
#
# shade_dose_response__max_effect:      [0.30, 0.95]
#   Maximum fraction of recruits suppressed by full canopy. Lower bound
#   means shade always suppresses at least some establishment even at
#   moderate canopy; upper bound means near-complete suppression.
#
# shade_dose_response__half_sat:        [500, 8000]
#   Canopy density at 50% suppression. Wide range -- this parameter is
#   in raw canopy_density() units and is hard to constrain from theory.
#
# R_max:                                [100, 2000]
#   Asymptotic annual recruit ceiling. Wide range because this parameter
#   is tightly coupled with juv_decline (high R_max + high juv mortality
#   can produce the same equilibrium as low R_max + low mortality).
#
# K_half:                               [200, 4000]
#   Flowering adults giving half R_max. Range should span from "saturates
#   quickly" to "still nearly linear at observed adult densities".
#
# age_first_flower_mean:                [3, 10]
#   Mean age at first reproduction. Range spans shrub to small-tree life
#   histories.

get_param_ranges <- function() {
  list(
    keys = c(
      "weibull_lambda",
      "senescence_dose_response__max_effect",
      "senescence_dose_response__half_sat",
      "juv_decline_dose_response__max_effect",
      "juv_decline_age_half_sat",
      "shade_dose_response__max_effect",
      "shade_dose_response__half_sat",
      "R_max",
      "K_half",
      "age_first_flower_mean"
    ),
    mins = c( 20,  0.00, 15,  0.30, 1,  0.30,  500,  100,  200,  3),
    maxs = c(120,  0.60, 60,  0.99, 8,  0.95, 8000, 2000, 4000, 10),
    Sigma = NULL   # NULL = independent marginals
                   # Replace with a 10×10 correlation matrix to impose
                   # covariance structure (e.g. if shade and juv_decline
                   # are expected to covary across sites)
  )
}
