# tests/testthat/test-range-ensemble.R
#
# Range-wide ensemble design (docs/methodological_notes.md B1-B3): per-run
# random streams, plan contents, parameter draws, category summaries.
# (Runs themselves are too slow for unit tests.)

.sites <- function() load_range_sites(SITE_PARAMS_PATH)

test_that("each run's plan depends only on its own index and the master seed", {
  s <- .sites(); b <- make_toy_params()
  a1 <- plan_range_run(7, s, b)$row
  a2 <- plan_range_run(7, s, b)$row
  expect_identical(a1, a2)
  # planning other runs first (or more runs) does not change run 7
  invisible(lapply(1:20, plan_range_run, sites = s, base = b))
  expect_identical(plan_range_run(7, s, b)$row, a1)
  expect_false(identical(plan_range_run(7, s, b, master_seed = 1L)$row$run_id, a1$run_id))
})

test_that("plans carry the cell's values and the run's own seed", {
  s <- .sites(); b <- make_toy_params()
  pl <- plan_range_run(3, s, b)
  k  <- pl$row$cell_row
  expect_equal(pl$params$annual_env_alpha, s$annual_env_alpha[k])
  expect_equal(pl$params$fire_prob_annual, s$fire_prob_annual[k])
  expect_equal(pl$params$seed, pl$row$sim_seed)
  expect_equal(pl$row$region, s$region[k])
  expect_true(all(s$region %in% c("NSW", "Qld")))
  expect_true(all(c("site_fire_prob", "site_rust", "resprout_rust_mortality") %in% names(pl$row)))
})

test_that("parameter draws land in the params list within their range", {
  b  <- make_toy_params()
  vals <- vapply(1:200, function(i) {
    set.seed(i); dr <- draw_params(b)
    expect_equal(dr$params$rust_dose_response$resprout$max_effect, dr$values$resprout_rust_mortality)
    dr$values$resprout_rust_mortality
  }, numeric(1))
  d <- DEFAULT_PARAM_DRAWS$resprout_rust_mortality
  expect_true(all(vals >= d$min & vals <= d$max))
  expect_identical(draw_params(b, list())$params, b)
})

test_that("cells are drawn roughly in proportion to abundance weight", {
  s <- .sites(); b <- make_toy_params()
  k <- vapply(1:3000, function(i) plan_range_run(i, s, b, param_draws = list())$row$cell_row, integer(1))
  top <- which.max(s$weight_n_raw)
  expect_equal(mean(k == top), s$weight_n_raw[top] / sum(s$weight_n_raw), tolerance = 0.02)
})

test_that("range_category_shares counts categories overall and by region", {
  df <- data.frame(pct_decline = c(10, 35, 55, 85, 40, 60),
                   region = c("NSW", "NSW", "NSW", "NSW", "Qld", "Qld"))
  sh <- range_category_shares(df)
  expect_equal(unname(sh["all", c("LC/NT", "VU", "EN", "CR")]), round(100 * c(1, 2, 2, 1) / 6, 1))
  expect_equal(unname(sh["Qld", "n"]), 2)
  expect_equal(unname(sh["NSW", "median_decline"]), 45)
})
