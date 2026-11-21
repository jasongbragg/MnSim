# tests/testthat/test-olney.R
#
# Olney exclusion-trial validation (docs/design_olney_validation.md; notes T20):
# rust protection, fixed yearly pressures, crown transparency, and the runner.

source(file.path(proj_root, "runner", "olney_validation.R"), local = TRUE)

dr <- function(max_effect, form, half_sat, hill) list(max_effect = max_effect, form = form,
                                                      half_sat = half_sat, hill = hill)

test_that("protected share: overlap of the year with the protected periods", {
  set.seed(1)
  p <- list(rust_protection = c(2011, 2014 + 2 / 12))
  expect_equal(vapply(2010:2015, rust_protected_share, 1, params = p), c(0, 1, 1, 1, 2 / 12, 0))
  expect_equal(rust_protected_share(list(rust_protection = c(2011.5, 2011.75, 2011.75, 2013)), 2011), 0.5)
  expect_equal(rust_protected_share(list(), 2011), 0)
  expect_error(rust_protected_share(list(rust_protection = 2011), 2011), "pairs")
})

test_that("protection scales the year's pressure without extra draws", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 1, annual_env_alpha = 14, annual_env_beta = 6)
  set.seed(3); a <- draw_annual_env(p, 5); ua <- runif(1)
  set.seed(3); b <- draw_annual_env(modifyList(p, list(rust_protection = c(5, 5.25))), 5); ub <- runif(1)
  expect_equal(b, a * 0.75); expect_identical(ua, ub)
  expect_identical(draw_annual_env(modifyList(p, list(rust_protection = c(0, 100))), 5), 0)
})

test_that("a fixed pressure series replaces the draw and makes no random draws", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 3, annual_env_series = c(`3` = 0.2, `4` = 0.9))
  set.seed(5); u0 <- runif(1)
  set.seed(5); x <- draw_annual_env(p, 4); u1 <- runif(1)
  expect_identical(x, 0.9); expect_identical(u0, u1)
  expect_identical(draw_annual_env(p, 2), 0)
  expect_error(draw_annual_env(p, 7), "no value for year 7")
})

test_that("full protection: no damage gain and no rust deaths", {
  set.seed(1)
  p <- make_toy_params(n_years = 15L, rust_start_year = 2, annual_env_alpha = 14, annual_env_beta = 6,
                       rust_hazard_mode = "damage", damage_gain_dose_response = dr(0.5, "linear", 0, 1),
                       damage_recovery = 0.3, damage_mortality_dose_response = dr(0.6, "sigmoid", 0.7, 6),
                       rust_protection = c(0, 1e4))
  p$rust_pressure <- 1
  r <- run_simulation(p)
  expect_true(all(r$individuals$rust_damage == 0))
  # With a fixed series (no Beta draws), full protection is exactly no rust
  p$annual_env_series <- stats::setNames(rep(0.9, 15), 1:15)
  q <- modifyList(p, list(rust_start_year = Inf, rust_protection = NULL))
  expect_identical(run_simulation(p)$census$N_alive, run_simulation(q)$census$N_alive)
})

test_that("crown transparency maps damage above the healthy baseline; dead trees count as 1", {
  set.seed(1)
  expect_equal(crown_transparency(c(TRUE, TRUE, FALSE), c(0, 0.5, 0.2), 0.3), c(0.3, 0.65, 1))
})

test_that("the runner: arms share their start and years; treated recovers, untreated does not", {
  set.seed(1)
  base <- make_toy_params(N0 = 300L, rust_hazard_mode = "damage",
                          annual_env_alpha = 14, annual_env_beta = 6,
                          damage_gain_dose_response = dr(1, "linear", 0, 1),
                          damage_recovery = 0.6, weibull_lambda = 1e9,
                          damage_mortality_dose_response = dr(0.5, "sigmoid", 0.9, 8))
  base$rust_pressure <- 1
  r <- run_olney_replicate(base, seed = 11L, tree_age = c(1, 200), baseline = 0.3, spinup_years = 20L)
  tr <- r$trees
  first <- tr[tr$model_year == 2010, ]
  expect_length(unique(first$mean_transparency[first$arm != "partially_treated"]), 2L)  # different trees
  expect_setequal(unique(tr$arm), names(OLNEY_ARMS))
  get <- function(arm, yr) tr$mean_transparency[tr$arm == arm & tr$model_year == yr]
  expect_lt(get("treated", 2012), get("treated", 2010))
  expect_lt(get("treated", 2012), 0.36)
  expect_gt(get("untreated", 2012), get("treated", 2012) + 0.2)
  expect_gt(get("partially_treated", 2013), get("partially_treated", 2012))
  expect_true(r$stand$stand_mortality >= 0 && r$stand$stand_mortality <= 1)
  expect_error(run_olney_replicate(make_toy_params(), 1L, c(1, 2), 0.3), "damage mode")
})

test_that("the observed August values load with their model years", {
  set.seed(1)
  o <- olney_observed()
  expect_setequal(unique(o$treatment), c("untreated", "treated", "partially_treated"))
  expect_equal(o$mean_transparency[o$treatment == "untreated" & o$model_year == 2010], 0.662)
  expect_true(all(o$model_year %in% 2010:2013))
})
