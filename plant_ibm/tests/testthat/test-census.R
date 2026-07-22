# tests/testthat/test-census.R

test_that("census_year reports N_IUCN as exactly the alive & flowering count", {
  p   <- make_toy_params()
  pop <- make_toy_pop(50, p)
  pop$flowering <- c(rep(TRUE, 10), rep(FALSE, 40))
  gt  <- matrix(0L, nrow = 50L, ncol = 1L)

  row <- census_year(pop, gt, t = 1L, n_rec = 3L, n_dead = 2L, params = p)
  expect_equal(row$N_IUCN, 10L)
  expect_equal(row$N_alive, sum(pop$alive))
})

test_that("census_year's N_juvenile and N_resprout match pop columns directly", {
  p   <- make_toy_params()
  pop <- make_toy_pop(80, p)
  gt  <- matrix(0L, nrow = 80L, ncol = 1L)

  row <- census_year(pop, gt, t = 1L, n_rec = 0L, n_dead = 0L, params = p)
  expect_equal(row$N_juvenile,
               sum(pop$alive & pop$age < pop$age_first_flower & !pop$resprout))
  expect_equal(row$N_resprout, sum(pop$alive & pop$resprout))
})

test_that("census_year passes through births and deaths counts unchanged", {
  p   <- make_toy_params()
  pop <- make_toy_pop(20, p)
  gt  <- matrix(0L, nrow = 20L, ncol = 1L)
  row <- census_year(pop, gt, t = 1L, n_rec = 7L, n_dead = 4L, params = p)
  expect_equal(row$births, 7L)
  expect_equal(row$deaths, 4L)
})

test_that("census_year flags extinction and gives NA mean_resist_score with nobody alive", {
  p   <- make_toy_params()
  pop <- make_toy_pop(20, p)
  pop$alive[] <- FALSE
  gt  <- matrix(0L, nrow = 20L, ncol = 1L)

  row <- census_year(pop, gt, t = 5L, n_rec = 0L, n_dead = 20L, params = p)
  expect_true(row$extinct)
  expect_equal(row$N_alive, 0L)
  expect_true(is.na(row$mean_resist_score))
})

test_that("census_year's freq_locus columns exactly match allele_freqs", {
  p   <- make_toy_params(resist_locus_effect = c(0.4, 0.35),
                          resist_dominance     = c(0.5, 0.5),
                          resist_freq0         = c(0.2, 0.1))
  pop <- make_toy_pop(60, p)
  gt  <- make_toy_gt(60, p)

  row   <- census_year(pop, gt, t = 1L, n_rec = 0L, n_dead = 0L, params = p)
  freqs <- allele_freqs(gt, pop$alive)
  expect_true(all(c("freq_locus1", "freq_locus2") %in% names(row)))
  expect_equal(row$freq_locus1, unname(freqs[1]))
  expect_equal(row$freq_locus2, unname(freqs[2]))
})

test_that("census_year's fire_event and rust_active flags match params exactly", {
  p   <- make_toy_params(fire_years = c(5L, 10L), rust_start_year = 8L)
  pop <- make_toy_pop(30, p)
  gt  <- make_toy_gt(30, p)

  # fire_this_year is now passed by the caller (simulate.R) not recomputed
  # internally -- so the test must supply it explicitly.
  row5 <- census_year(pop, gt, t = 5L, n_rec = 0L, n_dead = 0L, params = p,
                       fire_this_year = TRUE)
  row7 <- census_year(pop, gt, t = 7L, n_rec = 0L, n_dead = 0L, params = p,
                       fire_this_year = FALSE)
  row8 <- census_year(pop, gt, t = 8L, n_rec = 0L, n_dead = 0L, params = p,
                       fire_this_year = FALSE)

  expect_true(row5$fire_event);  expect_false(row5$rust_active)
  expect_false(row7$fire_event); expect_false(row7$rust_active)
  expect_true(row8$rust_active)
})

test_that("pct_decline computes the exact percentage between two census years", {
  census <- data.frame(year = c(2010, 2020), N_IUCN = c(100, 60))
  expect_equal(pct_decline(census, 2010, 2020), 40)
})

test_that("pct_decline returns NA for a missing year or a zero baseline", {
  census <- data.frame(year = c(2010, 2020), N_IUCN = c(100, 60))
  expect_true(is.na(pct_decline(census, 1999, 2020)))

  census0 <- data.frame(year = c(2010, 2020), N_IUCN = c(0, 0))
  expect_true(is.na(pct_decline(census0, 2010, 2020)))
})
