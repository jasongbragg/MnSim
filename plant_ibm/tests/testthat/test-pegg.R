# tests/testthat/test-pegg.R
#
# compute_pegg_stats() on a hand-built pop (no simulation), plus short
# end-to-end checks that death_cause is recorded, the run stops at the
# survey, and the fates partition the pre-fire cohort.

.pegg_pop <- function() {
  # Founders encoded with birth_year = -(age at spin-up); spin-up 1910.
  # rows: 1 fire-killed | 2 died of mortality in the fire year (post-fire)
  #       3 died 2021 | 4 still resprouting | 5 recovered
  #       6 died before the fire (not in cohort) | 7 recovered (recruit 1950)
  #       8 born after the fire (not in cohort)
  data.frame(
    birth_year  = c(-30L, -30L, -30L, -30L, -30L, -30L, 1950L, 2020L),
    alive       = c(FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE,  TRUE),
    death_year  = c(2020L, 2020L, 2021L, NA,    NA,    1990L, NA,    NA),
    death_cause = c("fire", "mortality", "mortality", NA, NA, "mortality", NA, NA),
    resprout    = c(FALSE, FALSE, FALSE, TRUE,  FALSE, FALSE, FALSE, FALSE)
  )
}
.pegg_census <- data.frame(year = 2018:2021, N_IUCN = c(10, 11, 5, 6))

test_that("compute_pegg_stats counts fire-year mortality deaths as post-fire", {
  s <- compute_pegg_stats(.pegg_pop(), .pegg_census, year_fire = 2020L,
                          year_survey = 2022L, year_spinup = 1910L)
  expect_equal(s$N_prefire, 6L)
  expect_equal(s$n_fire_dead, 1L)
  expect_equal(s$N_resprouted, 5L)
  expect_equal(s$n_postfire_dead, 2L)          # rows 2 and 3
  expect_equal(s$fraction_resprout_dead, 2 / 5)
  expect_equal(s$N_iucn_prefire, 11)
  expect_equal(s$N_iucn_assess, 6)
})

test_that("compute_pegg_stats fates partition the pre-fire cohort", {
  s <- compute_pegg_stats(.pegg_pop(), .pegg_census, 2020L, 2022L, 1910L)
  expect_equal(s$n_fire_dead + s$n_postfire_dead + s$n_resprouting + s$n_recovered,
               s$N_prefire)
  expect_equal(s$fraction_fire + s$fraction_postfire + s$fraction_resprout +
               s$fraction_recover, 1)
})

test_that("compute_pegg_stats refuses runs that continue past the survey", {
  pop <- .pegg_pop(); pop$death_year[3] <- 2023L
  expect_error(compute_pegg_stats(pop, .pegg_census, 2020L, 2022L, 1910L),
               "past the survey")
})

test_that("compute_pegg_stats requires death_cause", {
  pop <- .pegg_pop(); pop$death_cause <- NULL
  expect_error(compute_pegg_stats(pop, .pegg_census, 2020L, 2022L, 1910L),
               "death_cause")
})

test_that("pegg_run_id depends on the scenario years", {
  p  <- get_default_params()
  s1 <- list(year_spinup = 1910L, year_rust = 2010L, year_fire = 2020L, year_survey = 2022L)
  s2 <- modifyList(s1, list(year_survey = 2023L))
  expect_false(identical(pegg_run_id(p, s1), pegg_run_id(p, s2)))
  expect_identical(pegg_run_id(p, s1), pegg_run_id(p, s1))
  expect_false(identical(pegg_run_id(p, s1), make_run_id(p)))
})

test_that("run_simulation records death_cause for fire and mortality deaths", {
  p <- make_toy_params()
  p$N0 <- 200L; p$n_years <- 6L; p$fire_years <- 3L; p$fire_prob_annual <- 0
  p$fire_p_fimp <- 1
  res <- run_simulation(p, year0 = 1L)
  pop <- res$individuals
  expect_true(all(is.na(pop$death_cause[pop$alive])))
  expect_true(all(pop$death_cause[!pop$alive] %in% c("fire", "mortality")))
  expect_true(any(pop$death_cause %in% "fire"))
  expect_true(all(pop$death_year[pop$death_cause %in% "fire"] == 3L))
})

test_that("run_pegg_validation ends at the survey and fates sum to 1", {
  p <- make_toy_params(); p$N0 <- 300L; p$seed <- 3L
  r <- run_pegg_validation(p, year_spinup = 1990L, plot = FALSE, verbose = FALSE)
  expect_equal(max(r$census$year), 2021L)
  expect_equal(r$fraction_fire + r$fraction_postfire + r$fraction_resprout +
               r$fraction_recover, 1)
  expect_true(r$fraction_resprout_dead >= 0 && r$fraction_resprout_dead <= 1)
})
