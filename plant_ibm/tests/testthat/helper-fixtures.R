# tests/testthat/helper-fixtures.R
#
# Auto-sourced by testthat before any test file runs. Sources all model
# code and defines small fixture helpers shared across test files.

# --- Resolve project root (testthat sources helpers with wd = tests/testthat/)
proj_root <- normalizePath(file.path("..", ".."))
setwd(proj_root)

# --- Source all model code -----------------------------------------------
# parlib/unittests/params.R defines get_default_params() with conservative
# "everything off" values safe for unit tests. The main params.R holds the
# working M. nodosa calibration and is NOT sourced here -- active mechanisms
# (senescence, juv_decline, shade, rust) in that file would leak into the
# matrix-equivalence fixtures and confound model-mechanic validation.
source("parlib/unittests/params.R")
source("R/individuals.R")
source("R/genetics.R")
source("R/mortality.R")
source("R/fire.R")
source("R/recruitment.R")
source("R/census.R")
source("R/simulate.R")
source("R/matrix_model.R")

# --- Minimal-population params for unit tests (fast, small) ---------------
make_toy_params <- function(...) {
  p <- get_default_params()
  p$N0               <- 200L
  p$K_half           <- 100L
  p$R_max            <- 60L
  p$n_years          <- 10L
  p$fire_years       <- integer(0)
  p$fire_prob_annual <- 0
  p$rust_start_year  <- Inf
  p$seed             <- 1L
  # Override with any named arguments
  modifyList(p, list(...))
}

# --- Small population for quick smoke runs --------------------------------
# IMPORTANT: overrides params$N0 <- n so the requested sample size is
# actually honoured -- statistical tests rely on this for their tolerance
# to be calibrated against the right sampling noise.
make_toy_pop <- function(n = 100, params = make_toy_params()) {
  params$N0 <- n
  set.seed(params$seed)
  init <- create_population(params)
  init$pop
}

make_toy_gt <- function(n = 100, params = make_toy_params()) {
  params$N0 <- n
  set.seed(params$seed)
  init <- create_population(params)
  init$resist_gt
}

# --- Level-1 params: geometric growth (constant hazard, AFR=1) ------------
# Expected lambda = s*(1+F) = 0.9 * 1.15 = 1.035 (closed form).
# K_half is set far beyond any population size reached during the run
# (N0=1000 growing at lambda~1.035 over 80 years reaches ~17000), so the
# Beverton-Holt nonlinearity stays negligible and recruitment tracks the
# matrix's linear approximation R = F*N essentially exactly. Validated
# numerically: |lambda_ibm - lambda_matrix| < 0.002 (tolerance used: 0.015).
make_geometric_params <- function() {
  p <- get_default_params()
  p$N0                    <- 1000L
  p$n_years               <- 80L
  p$seed                  <- 7L
  p$weibull_k             <- 1.0      # constant hazard: h=0.1, s=0.9
  p$weibull_lambda        <- 10.0
  p$K_half                <- 1e9      # >> reachable N: keeps bevholt linear
  p$R_max                 <- 0.15 * 1e9   # F = R_max/K_half = 0.15
  p$age_first_flower_mean <- 1
  p$age_first_flower_sd   <- 0
  p$fire_years            <- integer(0)
  p$fire_prob_annual      <- 0
  p$rust_start_year       <- Inf
  # Fully resistant population: no genetic source of mortality variation
  p$resist_locus_effect   <- c(1.0)
  p$resist_dominance      <- c(1.0)
  p$resist_freq0          <- c(1.0)
  p
}

# --- Level-2 params: age-structured Leslie (Weibull k=2, AFR=5) -----------
# Lambda has no closed form here -- compared against the dominant
# eigenvalue of build_leslie(). Validated: |lambda_ibm - lambda_matrix|
# < 0.002 (tolerance used: 0.020).
make_leslie_params <- function() {
  p <- get_default_params()
  p$N0                    <- 1000L
  p$n_years               <- 80L
  p$seed                  <- 13L
  p$weibull_k             <- 2.0
  p$weibull_lambda        <- 25.0
  p$K_half                <- 1e9
  p$R_max                 <- 0.10 * 1e9   # F = 0.10 per surviving adult
  p$age_first_flower_mean <- 5
  p$age_first_flower_sd   <- 0
  p$fire_years            <- integer(0)
  p$fire_prob_annual      <- 0
  p$rust_start_year       <- Inf
  p$resist_locus_effect   <- c(1.0)
  p$resist_dominance      <- c(1.0)
  p$resist_freq0          <- c(1.0)
  p
}

# --- Level-3 params: Lefkovitch (constant hazard + FIXED fire years + ------
#     geometric resprout recovery) -----------------------------------------
#
# Uses a fixed, deterministic fire_years schedule rather than a stochastic
# fire_prob_annual. This is a deliberate test-design choice, not a
# simplification of convenience: with random annual fire, the realised
# stochastic growth rate of a Markov-switching matrix model is NOT in
# general equal to the dominant eigenvalue of the probability-weighted
# mean matrix (1-p)*M_nofire + p*M_fire (a well-known result in stochastic
# demography -- the "extra" effect of variance in switching reduces
# realised growth relative to the deterministic mean-matrix eigenvalue).
# An asymptotic-lambda comparison against the mean-matrix eigenvalue
# would therefore not be a clean implementation test: it conflates "did I
# implement the per-year transition mechanics correctly" with "does this
# specific draw of fire timing happen to land close to its long-run
# expectation," and was empirically far too noisy across seeds (final-
# year ratios IBM/matrix ranging ~0.85-1.03 in exploration).
#
# Fixing fire_years and projecting the matrix through the IDENTICAL
# year-by-year fire/no-fire sequence (see project_lefkovitch_seq) removes
# that confound entirely: both trajectories experience exactly the same
# fire events in the same years, so any remaining mismatch reflects
# genuine demographic sampling noise in the IBM, not switching-process
# noise. With N0=20000 over 40 years (two fixed fires at years 15 and 30),
# final-year IBM/matrix ratios were validated across 7 seeds to fall in
# [0.965, 1.033] -- comfortably inside the 0.08 tolerance used by the test.
make_lefkovitch_params <- function() {
  p <- get_default_params()
  p$N0                      <- 20000L
  p$n_years                 <- 40L
  p$seed                    <- 42L
  p$weibull_k               <- 1.0     # constant hazard, same s as Level 1
  p$weibull_lambda          <- 10.0
  p$K_half                  <- 1e9
  p$R_max                   <- 0.18 * 1e9   # F = 0.18
  p$age_first_flower_mean   <- 1
  p$age_first_flower_sd     <- 0
  p$fire_years              <- c(15L, 30L)  # fixed, deterministic schedule
  p$fire_prob_annual        <- 0
  p$fire_kill_prob          <- 0.30
  p$resprout_recovery       <- "geometric"  # required for matrix equivalence
  p$resprout_prob_recovery  <- 1 / 3
  p$rust_start_year         <- Inf
  p$resist_locus_effect     <- c(1.0)
  p$resist_dominance        <- c(1.0)
  p$resist_freq0            <- c(1.0)
  p
}

# --- Build the IBM's initial founding population + the matching matrix ----
#     stage vector n0 = c(n_Juvenile, n_Adult, n_Resprout). Used by the
#     Level-3 trajectory comparison so both the IBM and the matrix start
#     from EXACTLY the same initial stage distribution.
make_lefkovitch_init <- function(params) {
  set.seed(params$seed)
  init <- create_population(params)
  pop0 <- init$pop
  AFR  <- as.integer(round(params$age_first_flower_mean))
  n0 <- c(
    sum(pop0$alive & !pop0$resprout & pop0$age <  AFR),
    sum(pop0$alive & !pop0$resprout & pop0$age >= AFR),
    sum(pop0$alive &  pop0$resprout)
  )
  list(pop0 = pop0, resist_gt0 = init$resist_gt, n0 = n0)
}
