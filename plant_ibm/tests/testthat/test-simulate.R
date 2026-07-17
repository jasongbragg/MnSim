# tests/testthat/test-simulate.R

test_that("run_simulation runs end-to-end and keeps pop/resist_gt row-aligned", {
  p   <- make_toy_params(n_years = 10L)
  res <- run_simulation(p, year0 = 1L)
  expect_equal(nrow(res$individuals), nrow(res$resist_gt))
  expect_true(all(c("individuals", "resist_gt", "census") %in% names(res)))
  expect_equal(nrow(res$census), p$n_years)
})

test_that("identical seed and params give byte-identical results", {
  p  <- make_toy_params()
  r1 <- run_simulation(p, year0 = 1L)
  r2 <- run_simulation(p, year0 = 1L)
  expect_identical(r1$census, r2$census)
})

test_that("different seeds give different trajectories", {
  p1 <- make_toy_params(seed = 1L)
  p2 <- make_toy_params(seed = 2L)
  r1 <- run_simulation(p1, year0 = 1L)
  r2 <- run_simulation(p2, year0 = 1L)
  expect_false(isTRUE(all.equal(r1$census$N_alive, r2$census$N_alive)))
})

test_that("extinction stops the simulation loop early and is flagged in the final census row", {
  p   <- make_toy_params(n_years = 50L, fire_years = 1:50, fire_kill_prob = 1)
  res <- run_simulation(p, year0 = 1L)
  expect_equal(nrow(res$census), 1L)
  expect_equal(tail(res$census$N_alive, 1), 0L)
  expect_true(tail(res$census$extinct, 1))
})

test_that("no alive individuals have NA resist_score after a multi-year run", {
  p   <- make_toy_params(n_years = 20L)
  res <- run_simulation(p, year0 = 1L)
  expect_false(any(is.na(res$individuals$resist_score[res$individuals$alive])))
})

test_that("all dead individuals have a recorded, non-NA death_year", {
  p    <- make_toy_params(n_years = 20L)
  res  <- run_simulation(p, year0 = 1L)
  dead <- !res$individuals$alive
  expect_false(any(is.na(res$individuals$death_year[dead])))
})

test_that("rust and counterfactual runs are identical up to rust onset (common random numbers)", {
  # Documents and locks in a deliberate design property: run_simulation()
  # always calls set.seed(params$seed) first, so two runs that only differ
  # in rust_start_year consume identical random draws (age structure,
  # genotypes, mortality, recruitment) right up until rust_modifiers()
  # starts returning nonzero values -- this is what lets a rust scenario
  # and its no-rust counterfactual be compared as a clean paired design.
  p_rust  <- make_toy_params(rust_start_year = 5L, n_years = 10L)
  p_cfact <- p_rust
  p_cfact$rust_start_year <- Inf

  res_rust  <- run_simulation(p_rust,  year0 = 1L)
  res_cfact <- run_simulation(p_cfact, year0 = 1L)

  pre_onset <- res_rust$census$year < 5L
  expect_identical(res_rust$census[pre_onset, ], res_cfact$census[pre_onset, ])
  expect_false(isTRUE(all.equal(res_rust$census, res_cfact$census)))
})

test_that("init_state correctly resumes a run from a saved population and genotype matrix", {
  p1   <- make_toy_params(n_years = 5L)
  res1 <- run_simulation(p1, year0 = 1L)
  init_state <- list(individuals = res1$individuals, resist_gt = res1$resist_gt)

  p2   <- make_toy_params(n_years = 5L)
  res2 <- run_simulation(p2, init_state = init_state, year0 = 6L)

  expect_equal(res2$census$year, 6:10)
  expect_true(nrow(res2$individuals) >= nrow(res1$individuals))
})
