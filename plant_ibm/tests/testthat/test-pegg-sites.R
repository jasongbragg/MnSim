# tests/testthat/test-pegg-sites.R
#
# Phase 2 pieces (docs/methodological_notes.md D3, C3, V6): the fixed
# fire-only resprout window, the persistent site effect on rust pressure,
# interval mortality from annual survival, and the multi-plot likelihood.

test_that("fixed resprout window: every fire survivor resprouts exactly resprout_yrs_base years without rust", {
  p <- make_toy_params(fire_kill_prob = 0, fire_p_fimp = 0,
                       resprout_yrs_base = 2L, resprout_delay_dist = "fixed")
  pop <- make_toy_pop(300, p)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_true(all(pop2$resprout_yrs_remain[pop2$alive] == 2L))
  expect_true(all(!pop2$flowering[pop2$alive]))
  # the original Poisson rule is still the default when the option is absent
  q <- make_toy_params(fire_kill_prob = 0, fire_p_fimp = 0, resprout_yrs_base = 2L)
  set.seed(5L)
  pop3 <- apply_fire(make_toy_pop(300, q), q, t = 1L)
  expect_gt(length(unique(pop3$resprout_yrs_remain[pop3$alive])), 1L)
})

test_that("fixed window: rust adds its delay on top of the two years", {
  p <- make_toy_params(fire_kill_prob = 0, fire_p_fimp = 0, resprout_yrs_base = 2L,
                       resprout_delay_dist = "fixed", rust_start_year = 1L, rust_pressure = 1)
  p$rust_dose_response$delay$max_effect <- 3
  pop <- make_toy_pop(500, p)
  pop2 <- apply_fire(pop, p, t = 5L)
  r <- pop2$resprout_yrs_remain[pop2$alive]
  expect_true(all(r >= 2L))
  expect_gt(mean(r), 2.5)
})

test_that("params.R uses the fixed two-year window", {
  e <- new.env(); sys.source(file.path(proj_root, "params.R"), e)
  p <- e$get_default_params()
  expect_equal(p$resprout_yrs_base, 2)
  expect_identical(p$resprout_delay_dist, "fixed")
  expect_equal(p$rust_site_sd, 0)
})

test_that("site effect: none drawn when rust_site_sd = 0 (random stream untouched)", {
  p <- make_toy_params(rust_site_sd = 0)
  set.seed(11); before <- .Random.seed
  z <- draw_rust_site_effect(p)
  expect_identical(.Random.seed, before)
  expect_equal(z, 0)
  q <- apply_rust_site_effect(p, z)
  expect_equal(q$annual_env_alpha, p$annual_env_alpha)
  expect_equal(q$rust_site_z, 0)
})

test_that("site effect shifts the mean on the logit scale and keeps the concentration", {
  p <- make_toy_params(annual_env_alpha = 14, annual_env_beta = 6)
  q <- apply_rust_site_effect(p, -1)
  expect_equal(q$annual_env_alpha + q$annual_env_beta, 20)
  expect_equal(q$annual_env_alpha / 20, plogis(qlogis(0.7) - 1))
  r <- set_rust_concentration(p, 5)
  expect_equal(r$annual_env_alpha / (r$annual_env_alpha + r$annual_env_beta), 0.7)
  expect_equal(r$annual_env_alpha + r$annual_env_beta, 5)
})

test_that("interval mortality assumes a constant death rate within each post-fire year", {
  S <- c(0.8, 0.64, 0.512)                      # 20% per year
  expect_equal(interval_mortality(S, 0, 12), 0.2)
  expect_equal(interval_mortality(S, 12, 24), 0.2)
  expect_equal(interval_mortality(S, 6, 18), 0.2)
  expect_equal(interval_mortality(S, 3, 12), 1 - 0.8^(9 / 12))
  expect_error(interval_mortality(S, 12, 40), "n_post")
})

test_that("the mixture likelihood rewards between-site spread when plots include zero deaths", {
  plots <- data.frame(plot = c("A", "B"), lon = c(150, 150), lat = c(-30, -30),
                      t1_months = 0, t2_months = 24, n_at_risk = 50, deaths = c(0, 25))
  sites <- data.frame(cell_id = 1L, lon_centre = 150, lat_centre = -30)
  # same mean two-year mortality (~0.25 per year on average), different spread
  narrow <- data.frame(conc = 20, site_sd = 0, cell_id = 1L, seed = 1:4, z = 0, max_effect = 0.5,
                       S12 = rep(0.75, 4), S24 = rep(0.5625, 4), S36 = rep(0.42, 4))
  wide   <- data.frame(conc = 2, site_sd = 0, cell_id = 1L, seed = 1:4, z = 0, max_effect = 0.5,
                       S12 = c(0.999, 0.999, 0.5, 0.5), S24 = c(0.998, 0.998, 0.25, 0.25),
                       S36 = c(0.997, 0.997, 0.12, 0.12))
  fit <- pegg_plot_loglik(rbind(narrow, wide), plots, sites)
  expect_equal(fit$loglik$conc[1], 2)                 # the wide set fits better
  pz <- fit$per_plot[fit$per_plot$plot == "A", ]
  expect_lt(pz$p_zero[pz$conc == 20], 0.01)
  expect_gt(pz$p_zero[pz$conc == 2], 0.3)
})
