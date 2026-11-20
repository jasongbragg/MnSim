# tests/testthat/test-fire-change.R
#
# Species fire exposure and the fire-change scenario
# (docs/design_fire_exposure_and_change.md; notes T15).

repo_root <- normalizePath(file.path(proj_root, ".."))

test_that("absent settings leave fire exactly as before", {
  set.seed(1)
  p <- make_toy_params(fire_prob_annual = 0.2, fire_p_fimp = 0.4)
  expect_identical(fire_prob_at(p, 2050L), 0.2)
  expect_identical(fire_params_at(p, 2050L), p)
  m <- load_species("m_nodosa", root = repo_root)
  expect_false(any(c("fire_exposure", "fire_change_year", "fire_change_rate",
                     "fire_change_size", "fire_change_ramp") %in% names(m)))
})

test_that("fire_exposure scales the probability for the whole run", {
  set.seed(1)
  p <- make_toy_params(fire_prob_annual = 0.2, fire_exposure = 0.25)
  expect_equal(fire_prob_at(p, 1700L), 0.05)
  expect_equal(fire_prob_at(modifyList(p, list(fire_exposure = 0)), 2100L), 0)
  # Exposure 0 means no fire at all in a run.
  q <- make_toy_params(fire_prob_annual = 0.5, fire_exposure = 0, n_years = 30L)
  expect_identical(sum(run_simulation(q)$census$fire_in_cell), 0L)
})

test_that("the change applies from its year, as a step or a ramp", {
  set.seed(1)
  p <- set_fire_change(make_toy_params(fire_prob_annual = 0.1, fire_p_fimp = 0.4),
                       list(year = 2026, rate = 3, size = 1.5))
  expect_equal(fire_prob_at(p, 2025L), 0.1)
  expect_equal(fire_prob_at(p, 2026L), 0.3)
  expect_equal(fire_params_at(p, 2025L)$fire_p_fimp, 0.4)
  expect_equal(fire_params_at(p, 2030L)$fire_p_fimp, 0.6)
  r <- set_fire_change(p, list(year = 2026, rate = 3, ramp = 4))
  expect_equal(vapply(2025:2031, function(t) fire_prob_at(r, t), 0),
               c(0.1, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3))
  big <- set_fire_change(make_toy_params(fire_prob_annual = 0.6, fire_p_fimp = 0.8), list(year = 2000, rate = 5, size = 5))
  expect_identical(fire_prob_at(big, 2001L), 1)
  expect_identical(fire_params_at(big, 2001L)$fire_p_fimp, 0.99)
})

test_that("a run is identical to the baseline, draw for draw, until the change year", {
  set.seed(1)
  p  <- make_toy_params(fire_prob_annual = 0.3, fire_p_fimp = 0.5, fire_extent_mode = "stand",
                        fire_size_sd = 0.2, n_years = 40L, seed = 7L)
  pc <- set_fire_change(p, list(year = 2020, rate = 3, size = 1.5))
  b  <- run_simulation(p,  year0 = 2000L)$census
  c  <- run_simulation(pc, year0 = 2000L)$census
  before <- b$year < 2020
  expect_identical(c[before, ], b[before, ])
  expect_gt(sum(c$fire_in_cell[c$year >= 2020]), sum(b$fire_in_cell[b$year >= 2020]))
})

test_that("the scenario changes the run ID, and must start after the gen-time window", {
  set.seed(1)
  p  <- make_toy_params(fire_prob_annual = 0.1)
  pc <- set_fire_change(p, list(year = 2026, rate = 2))
  expect_false(identical(make_run_id(p), make_run_id(pc)))
  expect_error(run_iucn_a4(set_fire_change(p, list(year = 1900, rate = 2)), plot = FALSE, verbose = FALSE),
               "inside or before the gen-time window")
  expect_error(set_fire_change(p, list(year = 2026, rate = -1)))
  expect_identical(fire_change_label(NULL), "none")
  expect_identical(fire_change_label(list(year = 2026, rate = 2, size = 1.5)),
                   "from 2026: rate x2, size x1.5, ramp 0 yr")
})

test_that("range plans pair up across fire-change scenarios", {
  set.seed(1)
  sites <- load_range_sites(SITE_PARAMS_PATH)
  base  <- get_default_params()
  a <- plan_range_run(3L, sites, base)
  b <- plan_range_run(3L, sites, base, fire_change = list(year = 2026, rate = 2, size = 1.5))
  keep <- c("cell_row", "cell_id", "sim_seed", "site_fire_prob", "site_fire_p_fimp", "site_rust")
  expect_identical(a$row[keep], b$row[keep])
  expect_identical(b$params$fire_change_rate, 2)
  expect_false(identical(a$row$run_id, b$row$run_id))
})
