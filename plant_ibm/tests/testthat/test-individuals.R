# tests/testthat/test-individuals.R

test_that("create_population produces N0 alive individuals with internally consistent flowering", {
  p   <- make_toy_params(N0 = 500L)
  out <- create_population(p)
  expect_equal(nrow(out$pop), 500L)
  expect_true(all(out$pop$alive))
  expect_true(all(out$pop$age_first_flower >= 1L))
  expect_equal(out$pop$flowering,
               out$pop$alive & !out$pop$resprout & out$pop$age >= out$pop$age_first_flower)
})

test_that("create_population's resist_score matches resist_score_from_gt(resist_gt, params)", {
  p   <- make_toy_params(N0 = 300L)
  out <- create_population(p)
  expect_equal(out$pop$resist_score, resist_score_from_gt(out$resist_gt, p))
})

test_that("create_population's resist_gt has dimensions matching N0 and n_loci", {
  p <- make_toy_params(N0 = 100L,
                        resist_locus_effect = c(0.4, 0.3, 0.2),
                        resist_dominance     = c(0.5, 0.5, 0.5),
                        resist_freq0         = c(0.1, 0.1, 0.1))
  out <- create_population(p)
  expect_equal(dim(out$resist_gt), c(100L, 3L))
})

test_that("create_population includes env_susceptibility column in (0,1) for all founders", {
  p   <- make_toy_params(N0 = 500L, microclim_alpha = 2, microclim_beta = 5)
  out <- create_population(p)
  expect_true("env_susceptibility" %in% names(out$pop))
  expect_true(all(out$pop$env_susceptibility >= 0 & out$pop$env_susceptibility <= 1))
})

test_that("env_susceptibility in recruits is drawn from the same distribution as founders", {
  # Statistical: mean of Beta(3,3) should be close to 0.5
  set.seed(77L)
  p   <- make_toy_params(N0 = 5000L, microclim_alpha = 3, microclim_beta = 3)
  out <- create_population(p)
  expect_equal(mean(out$pop$env_susceptibility), 0.5, tolerance = 0.03)
})

test_that("make_recruit_rows assigns sequential IDs and includes env_susceptibility", {
  p        <- make_toy_params()
  pop      <- make_toy_pop(20, p)
  pop$id   <- seq_len(20L) * 10L   # max id = 200
  new_rows <- make_recruit_rows(5L, t = 99L, mother_idx = rep(1L, 5L),
                                 father_idx = rep(2L, 5L), pop, p)
  expect_equal(new_rows$id, 201:205)
  expect_true("env_susceptibility" %in% names(new_rows))
  expect_true(all(new_rows$env_susceptibility >= 0 & new_rows$env_susceptibility <= 1))
  expect_equal(new_rows$age, rep(0L, 5L))
})

test_that("make_recruit_rows returns NULL for zero recruits", {
  p   <- make_toy_params()
  pop <- make_toy_pop(5, p)
  expect_null(make_recruit_rows(0L, t = 1L, integer(0), integer(0), pop, p))
})


test_that("founders get calendar birth years and a founder flag; recruits do not", {
  p <- make_toy_params(); p$N0 <- 200L; p$n_years <- 15L
  res <- run_simulation(p, year0 = 1900L)
  pop <- res$individuals
  f <- pop$founder
  expect_true(any(f) && any(!f))
  expect_true(all(pop$birth_year[f] < 1900L))
  expect_true(all(pop$birth_year[!f] >= 1900L & pop$birth_year[!f] <= 1914L))
  # age = year - birth_year at the end of the last year, founders and recruits alike
  al <- pop$alive
  expect_true(all(pop$age[al] == 1914L - pop$birth_year[al]))
  expect_identical(calendar_birth_year(pop, 1900L), pop$birth_year)
})

test_that("calendar_birth_year and is_founder still read the legacy encoding", {
  old <- data.frame(birth_year = c(-10, 150), alive = TRUE)   # no founder column
  expect_equal(calendar_birth_year(old, 100), c(90, 150))
  expect_equal(is_founder(old, 100), c(TRUE, FALSE))
})
