# tests/testthat/test-matrix-equivalence.R
#
# Validation tests: contrive IBM parameterisations that degenerate to a
# known matrix population model, run both, and confirm the IBM's realised
# dynamics match the matrix's to within sampling tolerance. A mismatch here
# is a genuine implementation error, not a calibration question.
#
# Levels of equivalence:
#
#   Level 1 -- Geometric (single-class, constant hazard, AFR=1)
#     Matrix: degenerate 1-class model. Lambda = s*(1+F) in closed form.
#     Comparison: asymptotic IBM growth rate (post burn-in) vs closed form
#     and vs the Leslie matrix's dominant eigenvalue.
#
#   Level 2 -- Leslie (age-structured Weibull, delayed flowering)
#     Matrix: full Leslie matrix built from the same Weibull params.
#     Comparison: asymptotic IBM growth rate vs Leslie dominant eigenvalue.
#
#   Level 3 -- Lefkovitch (Leslie + Resprout stage, fire-driven)
#     Matrix: 3-stage Lefkovitch (Juvenile, Adult, Resprout), projected
#     year-by-year through a FIXED fire schedule identical to the IBM's.
#     Comparison: year-by-year trajectory, not asymptotic lambda -- see
#     make_lefkovitch_params() in helper-fixtures.R for why a stochastic
#     fire process would make an asymptotic-lambda comparison unreliable.
#
# What the IBM can do that even the Lefkovitch level cannot represent, and
# is therefore NOT covered by equivalence testing here:
#   - individual genetic resistance variation (genetics tested separately
#     in test-genetics.R)
#   - density dependence in mortality or recruitment
#   - individual variation in age-at-first-flower
#   - resistance-scaled resprout delay (fire x rust compounding)
# Those are the domain of diagnostics/diag_rust_pressure.R, which checks
# biological plausibility rather than exact equivalence.

# ===========================================================================
# Level 1: Geometric (constant hazard, AFR=1) -- closed-form benchmark
# ===========================================================================

test_that("L1 geometric: IBM asymptotic lambda matches closed form s*(1+F)", {
  p <- make_geometric_params()
  s <- 1 - weibull_hazard(1, p)        # constant survival (k=1)
  F <- p$R_max / p$K_half               # per-adult fecundity, linear BH limit
  lambda_theory <- s * (1 + F)

  res        <- run_simulation(p, year0 = 1L, verbose = FALSE)
  lambda_ibm <- estimate_lambda_ibm(res$census, burn = 30L)

  expect_equal(lambda_ibm, lambda_theory, tolerance = 0.015,
    label = sprintf("IBM lambda=%.4f vs theory=%.4f", lambda_ibm, lambda_theory))
})

test_that("L1 geometric: IBM asymptotic lambda matches Leslie matrix eigenvalue", {
  p             <- make_geometric_params()
  lambda_matrix <- dominant_lambda(build_leslie(p, max_age = 80L))

  res        <- run_simulation(p, year0 = 1L, verbose = FALSE)
  lambda_ibm <- estimate_lambda_ibm(res$census, burn = 30L)

  expect_equal(lambda_ibm, lambda_matrix, tolerance = 0.015,
    label = sprintf("IBM=%.4f vs Leslie=%.4f", lambda_ibm, lambda_matrix))
})

test_that("L1 geometric: Leslie matrix eigenvalue itself matches the closed form", {
  # Purely deterministic internal-consistency check -- no IBM involved.
  # Confirms build_leslie()'s degenerate single-class case is correct
  # before trusting it as a benchmark for the IBM.
  p             <- make_geometric_params()
  s             <- 1 - weibull_hazard(1, p)
  F             <- p$R_max / p$K_half
  lambda_theory <- s * (1 + F)
  lambda_matrix <- dominant_lambda(build_leslie(p, max_age = 80L))
  expect_equal(lambda_matrix, lambda_theory, tolerance = 0.001)
})

# ===========================================================================
# Level 2: Leslie (age-structured, Weibull k=2, AFR=5)
# ===========================================================================

test_that("L2 Leslie: IBM asymptotic lambda matches Leslie matrix eigenvalue", {
  p             <- make_leslie_params()
  lambda_matrix <- dominant_lambda(build_leslie(p, max_age = 100L))

  res        <- run_simulation(p, year0 = 1L, verbose = FALSE)
  lambda_ibm <- estimate_lambda_ibm(res$census, burn = 30L)

  expect_equal(lambda_ibm, lambda_matrix, tolerance = 0.020,
    label = sprintf("IBM=%.4f vs Leslie=%.4f", lambda_ibm, lambda_matrix))
})

test_that("L2 Leslie: IBM and matrix agree on direction under sub-replacement fecundity", {
  # A second, independent check at different params: confirms the
  # agreement isn't an artefact of one specific (F, AFR) combination.
  p        <- make_leslie_params()
  p$R_max  <- 0.05 * p$K_half   # F=0.05, well below replacement here
  p$seed   <- 99L

  lambda_matrix <- dominant_lambda(build_leslie(p, max_age = 100L))
  res        <- run_simulation(p, year0 = 1L, verbose = FALSE)
  lambda_ibm <- estimate_lambda_ibm(res$census, burn = 30L)

  expect_true(lambda_matrix < 1)
  expect_true(lambda_ibm    < 1)
  expect_equal(lambda_ibm, lambda_matrix, tolerance = 0.02)
})

# ===========================================================================
# Level 3: Lefkovitch (constant hazard + fixed fire schedule + geometric
#          resprout) -- year-by-year trajectory comparison
# ===========================================================================

test_that("L3 Lefkovitch: IBM trajectory tracks the matrix through fixed fire years", {
  p    <- make_lefkovitch_params()
  init <- make_lefkovitch_init(p)   # exact same starting pop for both

  res <- run_simulation(p,
                         init_state = list(individuals = init$pop0,
                                            resist_gt   = init$resist_gt0),
                         year0 = 1L, verbose = FALSE)

  M_nf   <- build_lefkovitch_nofire(p)
  M_fire <- build_lefkovitch_fire(p)
  N_mat  <- project_lefkovitch_seq(M_nf, M_fire, p$fire_years, init$n0, p$n_years)
  N_ibm  <- c(sum(init$n0), res$census$N_alive)

  # Compare at checkpoints spanning both fire events, not just the end --
  # this catches a mismatch that happens to cancel out by the final year.
  checkpoints <- c(10L, 15L, 25L, 30L, 40L) + 1L   # +1 for the year-0 offset
  ratios <- N_ibm[checkpoints] / N_mat[checkpoints]

  expect_true(all(ratios > 0.92 & ratios < 1.08),
    label = paste("ratios:", paste(sprintf("%.4f", ratios), collapse = ", ")))
})

test_that("L3 Lefkovitch: fixed-fire matrix trajectory shows the expected dips at fire years", {
  # Deterministic, IBM-independent sanity check on the matrix construction
  # itself: N should be lower immediately after a fire year than a
  # no-fire counterfactual trajectory from the same starting point.
  p      <- make_lefkovitch_params()
  init   <- make_lefkovitch_init(p)
  M_nf   <- build_lefkovitch_nofire(p)
  M_fire <- build_lefkovitch_fire(p)

  N_fire   <- project_lefkovitch_seq(M_nf, M_fire, p$fire_years, init$n0, p$n_years)
  N_nofire <- project_lefkovitch_seq(M_nf, M_fire, integer(0),   init$n0, p$n_years)

  # At year 16 (just after the first fire at year 15), the fire trajectory
  # should be measurably behind the no-fire counterfactual
  expect_lt(N_fire[17L], N_nofire[17L])
})

test_that("L3 Lefkovitch: zero fire_kill_prob makes fire-year survival equal no-fire survival", {
  # Internal consistency check on build_lefkovitch_fire(): with no killing,
  # the fire matrix's survival terms should equal the no-fire matrix's
  # baseline survival s (fire then only suppresses flowering that year,
  # which the M[1,] / M[2,] zero rows already capture).
  p   <- make_lefkovitch_params()
  p$fire_kill_prob <- 0
  s   <- max(1 - weibull_hazard(2, p), 0)
  M_f <- build_lefkovitch_fire(p)
  expect_equal(unique(M_f[3L, ]), s, tolerance = 1e-10)
})

# ===========================================================================
# Matrix-builder internal sanity checks (deterministic, no IBM needed)
# ===========================================================================

test_that("build_leslie: fecundity is zero for ages below age_first_flower", {
  p <- make_leslie_params()   # AFR = 5
  L <- build_leslie(p, max_age = 50L)
  for (age_j in 0:3) {
    expect_equal(L[1L, age_j + 1L], 0,
      label = paste("fecundity at age", age_j, "should be 0 (AFR=5)"))
  }
})

test_that("build_leslie: subdiagonal entries equal Weibull survival probabilities", {
  p <- make_leslie_params()
  L <- build_leslie(p, max_age = 30L)
  for (j in 0:29) {
    expected_s <- max(1 - weibull_hazard(j + 1, p), 0)
    expect_equal(L[j + 2L, j + 1L], expected_s, tolerance = 1e-10,
      label = paste("survival age", j, "to", j + 1))
  }
})

test_that("build_lefkovitch_nofire: with AFR=1, Juvenile and Adult columns are identical", {
  # Documents the key timing subtlety: age-0 individuals age to 1 within
  # the same simulation year, which already meets AFR=1, so they flower
  # exactly like existing adults. A naive implementation that zeroes the
  # J column's fecundity entry would be a timing bug, not a simplification.
  p    <- make_lefkovitch_params()
  M_nf <- build_lefkovitch_nofire(p)
  expect_equal(M_nf[, 1L], M_nf[, 2L])
})

test_that("project_matrix: long-run growth rate converges to the dominant eigenvalue", {
  p             <- make_leslie_params()
  L             <- build_leslie(p, max_age = 50L)
  lambda_matrix <- dominant_lambda(L)
  n0            <- rep(100, ncol(L))
  proj          <- project_matrix(L, n0, n_years = 200L)
  N             <- rowSums(proj[, -1L])
  ratios        <- N[152:200] / N[151:199]
  expect_equal(mean(ratios), lambda_matrix, tolerance = 0.001)
})
