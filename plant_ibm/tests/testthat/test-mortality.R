# tests/testthat/test-mortality.R

# Helper: minimal pop row for nominate_deaths() that includes all new columns
make_pop_rows <- function(n, age, afr = 5, resprout = FALSE, resist = 0,
                           ryr = 0L, rtot = 0L, env_susc = 0.5) {
  data.frame(alive = TRUE, age = age, age_first_flower = afr,
             resprout = resprout, resprout_yrs_remain = ryr,
             resprout_total_yrs = rtot, resist_score = resist,
             env_susceptibility = env_susc)[rep(1, n), ]
}

# --- weibull_hazard -------------------------------------------------------

test_that("weibull_hazard is negligible (not exactly zero) at age 0 when k > 1", {
  p <- make_toy_params(weibull_k = 2.0, weibull_lambda = 20.0)
  h0 <- weibull_hazard(0, p)
  expect_lt(h0, 1e-3)
  expect_gt(h0, 0)
})

test_that("weibull_hazard is constant for k = 1 (memoryless)", {
  p <- make_toy_params(weibull_k = 1.0, weibull_lambda = 10.0)
  h <- weibull_hazard(1:50, p)
  expect_true(all(abs(h - h[1]) < 1e-10))
})

test_that("weibull_hazard is strictly increasing for k > 1 (senescence)", {
  p <- make_toy_params(weibull_k = 2.0, weibull_lambda = 20.0)
  h <- weibull_hazard(1:80, p)
  expect_true(all(diff(h) > 0))
})

test_that("weibull_hazard equals 1/lambda at every age when k = 1", {
  p <- make_toy_params(weibull_k = 1.0, weibull_lambda = 12.0)
  expect_equal(weibull_hazard(5, p), 1 / 12, tolerance = 1e-10)
})

test_that("weibull_k = 1 does NOT turn weibull_hazard off -- it remains a substantial nonzero constant unless lambda is also large", {
  # This is exactly the fact that caused real confusion when documenting
  # how to pair senescence_hazard() with weibull_k=1 (see README.md /
  # params.R) -- k=1 removes the AGE-DEPENDENT shape, not the hazard
  # itself. At the project default lambda=30, that's a flat ~3.3%
  # annual hazard, not zero.
  p_default_lambda <- make_toy_params(weibull_k = 1.0, weibull_lambda = 30)
  expect_equal(weibull_hazard(0, p_default_lambda), 1/30, tolerance = 1e-10)
  expect_gt(weibull_hazard(0, p_default_lambda), 0.03)   # NOT negligible

  # Only a large lambda actually makes the k=1 hazard negligible.
  p_large_lambda <- make_toy_params(weibull_k = 1.0, weibull_lambda = 1e9)
  expect_lt(weibull_hazard(0, p_large_lambda), 1e-8)
})

test_that("weibull_hazard is vectorised over age", {
  p    <- make_toy_params(weibull_k = 1.5, weibull_lambda = 15.0)
  ages <- 1:10
  h    <- weibull_hazard(ages, p)
  h_scalar <- sapply(ages, function(a) weibull_hazard(a, p))
  expect_equal(h, h_scalar)
})

# --- dd_hazard ------------------------------------------------------------

test_that("dd_hazard is zero when dd_alpha = 0", {
  p <- make_toy_params(dd_alpha = 0)
  expect_equal(dd_hazard(500, p), 0)
})

test_that("dd_hazard equals dd_alpha exactly when N = K", {
  p <- make_toy_params(K = 1000L, dd_alpha = 0.4)
  expect_equal(dd_hazard(1000, p), 0.4)
})

test_that("dd_hazard scales linearly with N/K", {
  p <- make_toy_params(K = 1000L, dd_alpha = 0.4)
  expect_equal(dd_hazard(500, p), 0.2)
  expect_equal(dd_hazard(250, p), 0.1)
})

# --- hill_weight ----------------------------------------------------------

test_that("hill_weight is exactly 1 at x = 0", {
  expect_equal(hill_weight(0,  500, 1), 1)
  expect_equal(hill_weight(0,  500, 3), 1)
  expect_equal(hill_weight(0, Inf,  1), 1)
})

test_that("hill_weight is exactly 0.5 at x = half_sat (for any hill)", {
  expect_equal(hill_weight(5,   5, 1), 0.5)
  expect_equal(hill_weight(5,   5, 3), 0.5)
  expect_equal(hill_weight(100, 100, 2), 0.5)
})

test_that("hill_weight is 1 everywhere when half_sat = Inf (turns off age-weighting)", {
  expect_equal(hill_weight(0,   Inf, 1), 1)
  expect_equal(hill_weight(100, Inf, 1), 1)
  expect_equal(hill_weight(1e6, Inf, 2), 1)
})

test_that("hill_weight is strictly decreasing in x for any positive half_sat", {
  w <- hill_weight(0:20, half_sat = 5, hill = 2)
  expect_true(all(diff(w) < 0))
})

test_that("hill_weight with large hill is more gate-like than hill = 1", {
  # At x = 3, half_sat = 5 (below threshold):
  # hill=1 gives a smooth intermediate value; hill=10 is closer to 1 (less suppressed)
  # because the steep transition hasn't hit yet (x < half_sat)
  expect_lt(hill_weight(3, 5, 1), hill_weight(3, 5, 10))
  # Above half_sat the reverse holds: large hill gives sharper suppression
  expect_gt(hill_weight(7, 5, 1), hill_weight(7, 5, 10))
})

# --- dd_age_weight --------------------------------------------------------

test_that("dd_age_weight = 1 at all ages when dd_age_half_sat = Inf (flat dd, backward compat)", {
  p <- get_default_params()
  expect_equal(dd_age_weight(0,  p), 1)
  expect_equal(dd_age_weight(50, p), 1)
})

test_that("dd_age_weight decreases with age when dd_age_half_sat is finite", {
  p <- make_toy_params(dd_age_half_sat = 5, dd_age_hill = 1)
  w <- dd_age_weight(0:20, p)
  expect_true(all(diff(w) < 0))
})

# --- rust_modifiers -------------------------------------------------------

test_that("rust_modifiers returns delay_extra = 0 before rust_start_year", {
  p  <- make_toy_params(rust_start_year = 50L)
  rm <- rust_modifiers(p, t = 49L)
  expect_equal(rm$delay_extra, 0)
})

test_that("rust_modifiers returns delay_extra = 0 when rust_pressure = 0", {
  p  <- make_toy_params(rust_start_year = 1L, rust_pressure = 0)
  rm <- rust_modifiers(p, t = 5L)
  expect_equal(rm$delay_extra, 0)
})

test_that("rust_modifiers delay_extra scales correctly with rust_pressure (linear form)", {
  p <- make_toy_params(rust_start_year = 1L, rust_pressure = 2.0)
  p$rust_dose_response$delay$max_effect <- 2.0
  p$rust_dose_response$delay$form       <- "linear"
  rm <- rust_modifiers(p, t = 5L)
  expect_equal(rm$delay_extra, 4.00, tolerance = 1e-10)
})

# --- dose_response --------------------------------------------------------

test_that("dose_response 'linear' returns 0 at pressure=0", {
  cfg <- list(max_effect = 1, form = "linear", half_sat = 1, hill = 1)
  expect_equal(dose_response(0, cfg), 0)
})

test_that("dose_response 'saturating' returns max_effect/2 at pressure=half_sat", {
  cfg <- list(max_effect = 1, form = "saturating", half_sat = 0.5, hill = 1)
  expect_equal(dose_response(0.5, cfg), 0.5, tolerance = 1e-10)
})

test_that("dose_response 'threshold' is 0 below threshold and max_effect above", {
  cfg <- list(max_effect = 0.8, form = "threshold", half_sat = 1.0, hill = 1)
  expect_equal(dose_response(0.9, cfg), 0)
  expect_equal(dose_response(1.0, cfg), 0.8)
  expect_equal(dose_response(1.5, cfg), 0.8)
})

test_that("dose_response raises error for unknown form", {
  cfg <- list(max_effect = 1, form = "banana", half_sat = 1, hill = 1)
  expect_error(dose_response(1, cfg))
})

# --- juv_decline_hazard ---------------------------------------------------

test_that("juv_decline_hazard is zero at the default max_effect=0", {
  p <- get_default_params()
  expect_equal(juv_decline_hazard(0,  0,   p), 0)
  expect_equal(juv_decline_hazard(0,  1e6, p), 0)
  expect_equal(juv_decline_hazard(10, 1e6, p), 0)
})

test_that("juv_decline_hazard is zero at any age when N_canopy=0", {
  p <- make_toy_params(juv_decline_dose_response = list(max_effect = 0.5, form = "saturating",
                                                          half_sat = 500, hill = 1),
                        juv_decline_age_half_sat = 2, juv_decline_age_hill = 1)
  expect_equal(juv_decline_hazard(0, 0, p), 0)
  expect_equal(juv_decline_hazard(5, 0, p), 0)
})

test_that("juv_decline_hazard decays via hill_weight (not exp) and matches formula exactly", {
  p <- make_toy_params(juv_decline_dose_response = list(max_effect = 0.4, form = "linear",
                                                          half_sat = 1, hill = 1),
                        juv_decline_age_half_sat = 2, juv_decline_age_hill = 1)
  N_canopy <- 1.5   # linear: height = 0.4 * 1.5 = 0.6
  for (age in c(0, 1, 2, 4)) {
    expected <- 0.6 * hill_weight(age, 2, 1)
    expect_equal(juv_decline_hazard(age, N_canopy, p), expected, tolerance = 1e-10)
  }
})

test_that("juv_decline_hazard is strictly decreasing in age", {
  p <- make_toy_params(juv_decline_dose_response = list(max_effect = 0.4, form = "saturating",
                                                          half_sat = 500, hill = 1),
                        juv_decline_age_half_sat = 2, juv_decline_age_hill = 1)
  h <- sapply(0:10, juv_decline_hazard, N_canopy = 2000, params = p)
  expect_true(all(diff(h) < 0))
})

test_that("juv_decline_hazard matches dose_response() for its height at age 0", {
  p <- make_toy_params(juv_decline_dose_response = list(max_effect = 0.5, form = "saturating",
                                                          half_sat = 200, hill = 1),
                        juv_decline_age_half_sat = 2, juv_decline_age_hill = 1)
  # At age 0: hill_weight(0, ...) = 1, so hazard = height * 1 = height
  expect_equal(juv_decline_hazard(0, 200, p), dose_response(200, p$juv_decline_dose_response))
})

# --- senescence_hazard -----------------------------------------------------

test_that("senescence_hazard is zero at the default max_effect=0", {
  p <- get_default_params()
  expect_equal(senescence_hazard(0,   p), 0)
  expect_equal(senescence_hazard(30,  p), 0)
  expect_equal(senescence_hazard(200, p), 0)
})

test_that("senescence_hazard is exactly 0 at age 0 regardless of max_effect", {
  p <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                          half_sat = 30, hill = 3))
  expect_equal(senescence_hazard(0, p), 0)
})

test_that("senescence_hazard is exactly max_effect/2 at age = half_sat", {
  p <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                          half_sat = 30, hill = 3))
  expect_equal(senescence_hazard(30, p), 0.3, tolerance = 1e-10)
})

test_that("senescence_hazard approaches max_effect as age grows large", {
  p <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                          half_sat = 30, hill = 3))
  expect_equal(senescence_hazard(10000, p), 0.6, tolerance = 1e-6)
})

test_that("senescence_hazard is strictly increasing in age", {
  p <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                          half_sat = 30, hill = 3))
  h <- sapply(0:100, senescence_hazard, params = p)
  expect_true(all(diff(h) > 0))
})

test_that("senescence_hazard never exceeds max_effect at any age", {
  p <- make_toy_params(senescence_dose_response = list(max_effect = 0.45, form = "sigmoid",
                                                          half_sat = 30, hill = 3))
  h <- sapply(c(0, 1, 10, 30, 100, 1000, 1e6), senescence_hazard, params = p)
  expect_true(all(h <= 0.45 + 1e-10))
})

test_that("larger hill makes senescence stay lower below half_sat, and closer to max_effect above it", {
  p_low_hill  <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                                    half_sat = 30, hill = 1))
  p_high_hill <- make_toy_params(senescence_dose_response = list(max_effect = 0.6, form = "sigmoid",
                                                                    half_sat = 30, hill = 8))
  # Below half_sat: steeper hill means LESS senescence so far (sharper later rise)
  expect_lt(senescence_hazard(10, p_high_hill), senescence_hazard(10, p_low_hill))
  # Above half_sat: steeper hill means closer to the ceiling already
  expect_gt(senescence_hazard(60, p_high_hill), senescence_hazard(60, p_low_hill))
})

# --- nominate_deaths: senescence integration -------------------------------

test_that("nominate_deaths(): senescence_hazard off by default leaves mortality unaffected by age beyond weibull", {
  set.seed(41L)
  n <- 5000
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9)  # near-zero background
  pop_young <- make_pop_rows(n, age = 1,   afr = 5)
  pop_old   <- make_pop_rows(n, age = 200, afr = 5)
  rate_young <- length(nominate_deaths(pop_young, t = 1L, p)) / n
  rate_old   <- length(nominate_deaths(pop_old,   t = 1L, p)) / n
  expect_equal(rate_young, rate_old, tolerance = 0.01)
})

test_that("nominate_deaths(): turning on senescence_dose_response elevates old-age mortality, paired with weibull_k=1", {
  set.seed(43L)
  n <- 5000
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9)
  p$senescence_dose_response <- list(max_effect = 0.7, form = "sigmoid", half_sat = 30, hill = 3)

  pop_young <- make_pop_rows(n, age = 1,   afr = 5)
  pop_old   <- make_pop_rows(n, age = 100, afr = 5)
  rate_young <- length(nominate_deaths(pop_young, t = 1L, p)) / n
  rate_old   <- length(nominate_deaths(pop_old,   t = 1L, p)) / n
  expect_gt(rate_old, rate_young)
})

test_that("nominate_deaths(): senescence_hazard combines with other hazard terms via the union formula, not additively", {
  # Verify the union formula directly: with weibull off (k=1, huge lambda),
  # dd off, juv_decline off, rust off, p_death should equal senescence_hazard exactly.
  set.seed(47L)
  n <- 5000
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9)
  p$senescence_dose_response <- list(max_effect = 0.5, form = "sigmoid", half_sat = 30, hill = 3)
  pop <- make_pop_rows(n, age = 30, afr = 5)   # age = half_sat -> expected p = 0.25
  rate <- length(nominate_deaths(pop, t = 1L, p)) / n
  expect_equal(rate, 0.25, tolerance = 0.02)
})

# --- nominate_deaths ------------------------------------------------------

test_that("nominate_deaths returns empty vector when no alive individuals", {
  p   <- make_toy_params()
  pop <- make_toy_pop(100, p)
  pop$alive[] <- FALSE
  dead <- nominate_deaths(pop, t = 1L, p)
  expect_length(dead, 0L)
})

test_that("nominate_deaths returns only alive row indices", {
  set.seed(42L)
  p   <- make_toy_params()
  pop <- make_toy_pop(200, p)
  dead <- nominate_deaths(pop, t = 1L, p)
  expect_true(all(pop$alive[dead]))
})

test_that("nominate_deaths: hazard = 1 exactly kills everyone (deterministic)", {
  set.seed(42L)
  p <- make_toy_params(weibull_k = 1, weibull_lambda = 1.0, dd_alpha = 0)
  pop  <- make_toy_pop(1000, p)
  dead <- nominate_deaths(pop, t = 1L, p)
  expect_equal(length(dead), sum(pop$alive))
})

test_that("nominate_deaths: zero hazard kills nobody", {
  set.seed(42L)
  p <- make_toy_params(weibull_k = 1, weibull_lambda = 1e9, dd_alpha = 0)
  pop  <- make_toy_pop(500, p)
  dead <- nominate_deaths(pop, t = 1L, p)
  expect_length(dead, 0L)
})

# --- disease triangle in nominate_deaths ----------------------------------

test_that("disease triangle: resistant plants (resist_score=1) have zero rust hazard regardless of env", {
  set.seed(7L)
  n <- 5000
  pop <- make_pop_rows(n, age = 1, afr = 5, resist = 1, env_susc = 1)
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                        rust_start_year = 1L, rust_pressure = 1.0)
  p$rust_dose_response$age_peak <- 0.9
  # annual_env_t = 1: maximum pathogen; env_susc = 1; but resist = 1 → host=0 → eff_susc=0
  rate <- length(nominate_deaths(pop, t = 1L, p, annual_env_t = 1)) / n
  expect_equal(rate, 0)
})

test_that("disease triangle: zero annual_env_t means no rust regardless of genetics or microsite", {
  set.seed(9L)
  n <- 5000
  pop <- make_pop_rows(n, age = 1, afr = 5, resist = 0, env_susc = 1)
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                        rust_start_year = 1L, rust_pressure = 1.0)
  p$rust_dose_response$age_peak <- 0.9
  rate <- length(nominate_deaths(pop, t = 1L, p, annual_env_t = 0)) / n
  expect_equal(rate, 0)
})

test_that("disease triangle: zero env_susceptibility protects even susceptible plants in bad years", {
  set.seed(11L)
  n <- 5000
  pop <- make_pop_rows(n, age = 1, afr = 5, resist = 0, env_susc = 0)
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                        rust_start_year = 1L, rust_pressure = 1.0)
  p$rust_dose_response$age_peak <- 0.9
  rate <- length(nominate_deaths(pop, t = 1L, p, annual_env_t = 1)) / n
  expect_equal(rate, 0)
})

test_that("disease triangle: mortality increases with effective_susceptibility", {
  set.seed(13L)
  n <- 8000
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                        rust_start_year = 1L, rust_pressure = 1.0)
  p$rust_dose_response$age_peak <- 0.8

  pop_low  <- make_pop_rows(n, age = 1, afr = 5, resist = 0.9, env_susc = 0.1)
  pop_high <- make_pop_rows(n, age = 1, afr = 5, resist = 0,   env_susc = 1.0)

  rate_low  <- length(nominate_deaths(pop_low,  t = 1L, p, annual_env_t = 1)) / n
  rate_high <- length(nominate_deaths(pop_high, t = 1L, p, annual_env_t = 1)) / n
  expect_gt(rate_high, rate_low)
})

test_that("rust continuous age-decay: young individuals face higher rust hazard than old ones", {
  set.seed(17L)
  n <- 5000
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                        rust_start_year = 1L, rust_pressure = 1.0)
  p$rust_dose_response$age_peak     <- 0.7
  p$rust_dose_response$age_floor    <- 0
  p$rust_dose_response$age_half_sat <- 3
  p$rust_dose_response$age_hill     <- 2

  pop_young <- make_pop_rows(n, age = 1,  afr = 5, resist = 0, env_susc = 1)
  pop_old   <- make_pop_rows(n, age = 30, afr = 5, resist = 0, env_susc = 1)

  rate_young <- length(nominate_deaths(pop_young, t = 1L, p, annual_env_t = 1)) / n
  rate_old   <- length(nominate_deaths(pop_old,   t = 1L, p, annual_env_t = 1)) / n
  expect_gt(rate_young, rate_old)
})

test_that("rust age_floor gives old plants a non-zero rust hazard", {
  set.seed(19L)
  n <- 5000
  p_no_floor  <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 1e9,
                                   rust_start_year = 1L, rust_pressure = 1.0)
  p_no_floor$rust_dose_response$age_peak  <- 0.5
  p_no_floor$rust_dose_response$age_floor <- 0

  p_floor <- p_no_floor
  p_floor$rust_dose_response$age_floor <- 0.2  # residual adult mortality

  pop <- make_pop_rows(n, age = 100, afr = 5, resist = 0, env_susc = 1)

  rate_no_floor <- length(nominate_deaths(pop, t = 1L, p_no_floor, annual_env_t = 1)) / n
  rate_floor    <- length(nominate_deaths(pop, t = 1L, p_floor,    annual_env_t = 1)) / n
  expect_gt(rate_floor, rate_no_floor)
})

# --- nominate_deaths: juv_decline integration ----------------------------

test_that("nominate_deaths(): juv_decline contributes nothing when N_canopy = 0", {
  set.seed(21L)
  n   <- 5000
  pop <- make_pop_rows(n, age = 1, afr = 5)   # all juveniles, canopy = 0

  p_off <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 30)
  p_on  <- p_off
  p_on$juv_decline_dose_response <- list(max_effect = 0.5, form = "saturating",
                                          half_sat = 500, hill = 1)
  expect_equal(length(nominate_deaths(pop, t = 1L, p_off)) / n,
               length(nominate_deaths(pop, t = 1L, p_on))  / n,
               tolerance = 0.02)
})

test_that("nominate_deaths(): juveniles elevated under mature canopy when juv_decline active", {
  set.seed(23L)
  n_juv <- 5000; n_adult <- 2000
  pop <- rbind(make_pop_rows(n_juv,  age = 1,  afr = 5),
               make_pop_rows(n_adult, age = 10, afr = 5))
  p_off <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 30)
  p_on  <- p_off
  p_on$juv_decline_dose_response <- list(max_effect = 0.5, form = "saturating",
                                          half_sat = 500, hill = 1)
  juv_rows <- 1:n_juv
  expect_gt(mean(juv_rows %in% nominate_deaths(pop, t = 1L, p_on)),
            mean(juv_rows %in% nominate_deaths(pop, t = 1L, p_off)))
})

test_that("nominate_deaths(): juv_decline does not affect adults or resprouters", {
  set.seed(29L)
  n <- 5000
  p_off <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 30)
  p_on  <- p_off
  p_on$juv_decline_dose_response <- list(max_effect = 0.9, form = "saturating",
                                          half_sat = 1, hill = 1)
  for (pop in list(make_pop_rows(n, age = 10, afr = 5, resprout = FALSE),
                   make_pop_rows(n, age = 10, afr = 5, resprout = TRUE, ryr = 2L, rtot = 4L))) {
    expect_equal(length(nominate_deaths(pop, t = 1L, p_off)) / n,
                 length(nominate_deaths(pop, t = 1L, p_on))  / n,
                 tolerance = 0.03)
  }
})

test_that("nominate_deaths(): younger juveniles face higher mortality under mature canopy", {
  set.seed(31L)
  n_each  <- 5000; n_adult <- 2000
  pop <- rbind(make_pop_rows(n_each, age = 0, afr = 5),
               make_pop_rows(n_each, age = 4, afr = 5),
               make_pop_rows(n_adult, age = 10, afr = 5))
  p <- make_toy_params(dd_alpha = 0, weibull_k = 1, weibull_lambda = 30,
                        juv_decline_age_half_sat = 2, juv_decline_age_hill = 1)
  p$juv_decline_dose_response <- list(max_effect = 0.5, form = "saturating",
                                       half_sat = 500, hill = 1)
  dead <- nominate_deaths(pop, t = 1L, p)
  expect_gt(mean(1:n_each %in% dead), mean((n_each+1):(2*n_each) %in% dead))
})
