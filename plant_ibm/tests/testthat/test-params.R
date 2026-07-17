# tests/testthat/test-params.R

test_that("get_default_params returns a complete parameter set", {
  p <- get_default_params()
  required <- c("N0", "n_years", "seed",
                "weibull_k", "weibull_lambda",
                "senescence_dose_response",
                "juv_decline_dose_response", "juv_decline_age_half_sat", "juv_decline_age_hill",
                "R_max", "K_half", "shade_dose_response",
                "age_first_flower_mean", "age_first_flower_sd",
                "selfing_rate",
                "fire_years", "fire_prob_annual",
                "fire_kill_prob", "fire_p_fimp", "fire_kill_scalar",
                "fire_kill_half_sat", "fire_kill_hill",
                "resprout_yrs_base", "resprout_recovery", "resprout_prob_recovery",
                "rust_start_year", "rust_pressure",
                "microclim_alpha", "microclim_beta",
                "annual_env_alpha", "annual_env_beta",
                "rust_dose_response",
                "rust_flower_dose_response", "rust_recruit_dose_response",
                "resist_locus_effect", "resist_dominance", "resist_freq0")
  missing <- required[!required %in% names(p)]
  expect_true(length(missing) == 0,
    label = paste("Missing params:", paste(missing, collapse = ", ")))
})

test_that("genetic architecture vectors all have matching length", {
  p      <- get_default_params()
  n_loci <- length(p$resist_locus_effect)
  expect_equal(length(p$resist_dominance), n_loci)
  expect_equal(length(p$resist_freq0), n_loci)
})

test_that("rust_dose_response has age, resprout, and delay components", {
  p <- get_default_params()
  expect_true(all(c("age_peak","age_floor","age_half_sat","age_hill") %in%
                    names(p$rust_dose_response)))
  expect_true(all(c("max_effect","form","half_sat","hill") %in%
                    names(p$rust_dose_response$resprout)))
  expect_true(all(c("max_effect","form","half_sat","hill") %in%
                    names(p$rust_dose_response$delay)))
})

test_that("resist_freq0 values are valid allele frequencies in [0, 1]", {
  p <- get_default_params()
  expect_true(all(p$resist_freq0 >= 0 & p$resist_freq0 <= 1))
})

test_that("default resprout_recovery is countdown (geometric mode is opt-in)", {
  expect_equal(get_default_params()$resprout_recovery, "countdown")
})

test_that("default rust_start_year is Inf (rust off by default)", {
  expect_equal(get_default_params()$rust_start_year, Inf)
})

test_that("disease-triangle Beta params default to Beta(1,1) = uniform", {
  p <- get_default_params()
  expect_equal(p$microclim_alpha, 1); expect_equal(p$microclim_beta, 1)
  expect_equal(p$annual_env_alpha, 1); expect_equal(p$annual_env_beta, 1)
})

test_that("new rust suppression mechanisms are off by default (max_effect = 0 for flower and recruit)", {
  p <- get_default_params()
  expect_equal(p$rust_flower_dose_response$max_effect, 0)
  expect_equal(p$rust_recruit_dose_response$max_effect, 0)
  expect_equal(p$rust_dose_response$age_floor, 0)
  expect_equal(p$rust_start_year, Inf)
})

test_that("shade_dose_response and juv_decline are off by default", {
  p <- get_default_params()
  expect_equal(p$shade_dose_response$max_effect, 0)
  expect_equal(p$juv_decline_dose_response$max_effect, 0)
})

test_that("senescence_dose_response is off by default, with half_sat=30 ready to use", {
  p <- get_default_params()
  expect_equal(p$senescence_dose_response$max_effect, 0)
  expect_equal(p$senescence_dose_response$half_sat, 30)
})

test_that("fire two-stage model is off by default (legacy mode)", {
  p <- get_default_params()
  expect_equal(p$fire_p_fimp, 0)
})


test_that("test suite is using parlib/unittests/params.R, not the working M. nodosa params.R", {
  # The test suite sources parlib/unittests/params.R (conservative baseline,
  # everything off) NOT params.R (working M. nodosa calibration, active
  # mechanisms). If someone accidentally swaps them back, this test catches
  # it immediately rather than having the matrix-equivalence tests fail with
  # a confusing mismatch.
  p <- get_default_params()
  expect_equal(p$rust_start_year, Inf,
    label = "test params must have rust off (rust_start_year=Inf)")
  expect_equal(p$senescence_dose_response$max_effect, 0,
    label = "test params must have senescence off")
  expect_equal(p$juv_decline_dose_response$max_effect, 0,
    label = "test params must have juv_decline off")
  expect_equal(p$shade_dose_response$max_effect, 0,
    label = "test params must have shade suppression off")
  expect_false("dd_alpha" %in% names(p),
    label = "test params must NOT contain dd_alpha (dd mortality removed)")
})
