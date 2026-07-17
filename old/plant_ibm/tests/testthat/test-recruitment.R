# tests/testthat/test-recruitment.R

# --- bevholt ----------------------------------------------------------

test_that("bevholt hits exactly half R_max at N = K_half, and zero at N = 0", {
  p <- make_toy_params(R_max = 1000, K_half = 500)
  expect_equal(bevholt(500, p), 500)
  expect_equal(bevholt(0, p), 0)
})

test_that("bevholt is monotonically increasing in N and bounded below R_max", {
  p    <- make_toy_params(R_max = 1000, K_half = 500)
  Ns   <- seq(0, 20000, by = 200)
  vals <- bevholt(Ns, p)
  expect_true(all(diff(vals) >= 0))
  expect_true(all(vals < p$R_max))
})

# --- sample_parents -----------------------------------------------------

test_that("sample_parents forces selfing when only one flowering individual exists", {
  p   <- make_toy_params()
  pop <- make_toy_pop(10, p)
  pop$flowering <- c(TRUE, rep(FALSE, 9))
  parents <- sample_parents(pop, n_rec = 200L, params = p)
  expect_true(all(parents$mother_idx == 1L))
  expect_true(all(parents$father_idx == 1L))
})

test_that("sample_parents never selfs when selfing_rate = 0 and n_flower > 1", {
  p   <- make_toy_params(selfing_rate = 0)
  pop <- make_toy_pop(30, p)
  pop$flowering <- TRUE
  parents <- sample_parents(pop, n_rec = 3000L, params = p)
  expect_true(all(parents$mother_idx != parents$father_idx))
})

test_that("sample_parents self-fertilisation rate increases with selfing_rate", {
  pop <- make_toy_pop(30, make_toy_params())
  pop$flowering <- TRUE

  p_low  <- make_toy_params(selfing_rate = 0.02)
  p_high <- make_toy_params(selfing_rate = 5)

  par_low  <- sample_parents(pop, n_rec = 5000L, params = p_low)
  par_high <- sample_parents(pop, n_rec = 5000L, params = p_high)

  rate_low  <- mean(par_low$mother_idx  == par_low$father_idx)
  rate_high <- mean(par_high$mother_idx == par_high$father_idx)
  expect_gt(rate_high, rate_low)
})

test_that("sample_parents returns NULL when nobody is flowering", {
  p   <- make_toy_params()
  pop <- make_toy_pop(10, p)
  pop$flowering <- FALSE
  expect_null(sample_parents(pop, n_rec = 10L, params = p))
})

# --- recruit (full step) -------------------------------------------------

test_that("recruit() adds zero rows and leaves pop/resist_gt unchanged when nobody flowers", {
  p   <- make_toy_params()
  pop <- make_toy_pop(20, p)
  pop$flowering <- FALSE
  gt  <- make_toy_gt(20, p)

  out <- recruit(pop, gt, t = 1L, params = p)
  expect_equal(out$n_rec, 0L)
  expect_identical(out$pop, pop)
  expect_identical(out$resist_gt, gt)
})

test_that("recruit() keeps pop and resist_gt row-aligned after adding recruits", {
  set.seed(11L)
  p   <- make_toy_params(R_max = 500, K_half = 20)
  pop <- make_toy_pop(50, p)
  pop$flowering <- TRUE
  gt  <- make_toy_gt(50, p)

  out <- recruit(pop, gt, t = 1L, params = p)
  expect_equal(nrow(out$pop), nrow(out$resist_gt))
  expect_true(nrow(out$pop) >= 50L)   # recruit() only ever adds rows
})

test_that("recruit() gives offspring of two homozygous-resistant parents resist_score = 1", {
  p <- get_default_params()
  p$resist_locus_effect <- c(1.0); p$resist_dominance <- c(1.0)
  p$R_max <- 5000; p$K_half <- 5

  pop <- make_toy_pop(10, p)
  pop$flowering <- TRUE
  gt  <- matrix(2L, nrow = 10L, ncol = 1L)   # everyone RR

  out <- recruit(pop, gt, t = 1L, params = p)
  expect_gt(out$n_rec, 0L)
  new_scores <- tail(out$pop$resist_score, out$n_rec)
  expect_true(all(new_scores == 1))
})

# --- canopy_density -------------------------------------------------------

test_that("canopy_density: full (non-resprouting) flowering-age individuals count as exactly 1", {
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = FALSE, resprout_yrs_remain = 0L, resprout_total_yrs = 0L
  )
  expect_equal(canopy_density(pop), 1)
})

test_that("canopy_density: a freshly-resprouted adult (remain == total) counts as exactly 0", {
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = TRUE, resprout_yrs_remain = 4L, resprout_total_yrs = 4L
  )
  expect_equal(canopy_density(pop), 0)
})

test_that("canopy_density: a resprouting adult halfway through recovery counts as exactly 0.5", {
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = TRUE, resprout_yrs_remain = 2L, resprout_total_yrs = 4L
  )
  expect_equal(canopy_density(pop), 0.5)
})

test_that("canopy_density: a resprouting adult right at the point of recovery (remain=0) counts as exactly 1", {
  # A transient state -- simulate.R's loop would flip resprout to FALSE
  # in the same step remain hits 0 -- but canopy_density() should still
  # treat full recovery progress as full canopy if ever handed this
  # snapshot directly.
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = TRUE, resprout_yrs_remain = 0L, resprout_total_yrs = 4L
  )
  expect_equal(canopy_density(pop), 1)
})

test_that("canopy_density: progress is clipped to [0, 1] even in an inconsistent state", {
  # remain > total shouldn't arise from the normal countdown, but defend
  # against it producing a negative weight regardless.
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = TRUE, resprout_yrs_remain = 6L, resprout_total_yrs = 4L
  )
  expect_equal(canopy_density(pop), 0)
})

test_that("canopy_density: resprout_total_yrs == 0 for a resprouting individual defaults to full weight", {
  # Defensive case: only reachable in practice under geometric
  # resprout_recovery mode (which never decrements resprout_yrs_remain),
  # not under the normal countdown mode (apply_fire() always draws
  # delay >= 1). Should not divide by zero or silently understate canopy.
  pop <- data.frame(
    alive = TRUE, age = 8, age_first_flower = 5,
    resprout = TRUE, resprout_yrs_remain = 3L, resprout_total_yrs = 0L
  )
  expect_equal(canopy_density(pop), 1)
})

test_that("canopy_density counts alive flowering-age individuals, excluding resprouting juveniles", {
  pop <- data.frame(
    alive               = c(TRUE,  TRUE,  TRUE, TRUE,  FALSE, TRUE),
    age                 = c(2,     7,     7,    2,     10,    5),
    age_first_flower    = c(5,     5,     5,    5,     5,     5),
    resprout            = c(FALSE, FALSE, TRUE, TRUE,  FALSE, FALSE),
    resprout_yrs_remain = c(0L,    0L,    1L,   2L,    0L,    0L),
    resprout_total_yrs  = c(0L,    0L,    4L,   4L,    0L,    0L)
  )
  # row1: juvenile, not flowering age -> not counted (0)
  # row2: flowering adult -> full weight (1)
  # row3: resprouting ADULT, 3/4 through recovery -> weight 0.75
  # row4: resprouting JUVENILE (age < AFR) -> must NOT be counted at all,
  #       regardless of its own recovery progress -- this is exactly the
  #       bug found and fixed while building this mechanism: apply_fire()
  #       sends surviving juveniles into resprout too, and a resprouting
  #       former-juvenile never had canopy to lose
  # row5: dead -> not counted regardless of age
  # row6: age exactly equals age_first_flower, not resprouting -> full weight (1)
  expect_equal(canopy_density(pop), 1 + 0.75 + 1)
})

test_that("canopy_density: a resprouting adult counts less than a flowering adult mid-recovery, the same once fully recovered", {
  pop_flowering    <- data.frame(alive = TRUE, age = 8, age_first_flower = 5,
                                  resprout = FALSE, resprout_yrs_remain = 0L, resprout_total_yrs = 0L)
  pop_resprout_mid <- data.frame(alive = TRUE, age = 8, age_first_flower = 5,
                                  resprout = TRUE, resprout_yrs_remain = 2L, resprout_total_yrs = 4L)
  pop_resprout_done <- data.frame(alive = TRUE, age = 8, age_first_flower = 5,
                                   resprout = TRUE, resprout_yrs_remain = 0L, resprout_total_yrs = 4L)

  expect_lt(canopy_density(pop_resprout_mid), canopy_density(pop_flowering))
  expect_equal(canopy_density(pop_resprout_done), canopy_density(pop_flowering))
})

test_that("canopy_density excludes a resprouting individual still below flowering age, regardless of recovery progress", {
  pop <- data.frame(alive = TRUE, age = 2, age_first_flower = 5,
                     resprout = TRUE, resprout_yrs_remain = 0L, resprout_total_yrs = 4L)
  expect_equal(canopy_density(pop), 0L)
})

# --- shade_suppression -----------------------------------------------------

test_that("shade_suppression is zero at the default max_effect = 0, regardless of canopy density", {
  p <- get_default_params()
  expect_equal(shade_suppression(0, p), 0)
  expect_equal(shade_suppression(1e6, p), 0)
})

test_that("shade_suppression matches dose_response() directly", {
  p <- make_toy_params(shade_dose_response = list(max_effect = 0.6, form = "saturating",
                                                   half_sat = 500, hill = 1))
  expect_equal(shade_suppression(500, p), dose_response(500, p$shade_dose_response))
})

test_that("shade_suppression is defensively capped at 1 even with an improper max_effect > 1", {
  p <- make_toy_params(shade_dose_response = list(max_effect = 5, form = "linear",
                                                   half_sat = 1, hill = 1))
  # "linear" is unbounded by construction -- shade_suppression()'s own
  # pmin(...,1) is what has to catch this, not dose_response() itself
  expect_equal(shade_suppression(10, p), 1)
})

# --- recruit() with shading active -----------------------------------------

test_that("recruit(): default params (max_effect=0) match the pre-shading bevholt() expectation", {
  set.seed(7L)
  p   <- make_toy_params(R_max = 3000, K_half = 60)
  pop <- make_toy_pop(100, p)
  pop$flowering <- TRUE
  gt  <- make_toy_gt(100, p)

  n_flowering <- sum(pop$alive & pop$flowering)
  expected    <- bevholt(n_flowering, p)   # shade off -> no discount

  counts <- replicate(200, recruit(pop, gt, t = 1L, params = p)$n_rec)
  expect_equal(mean(counts), expected, tolerance = 0.05)
})

test_that("recruit(): mean recruit count matches bevholt() x (1 - shade_suppression) exactly", {
  set.seed(11L)
  p <- make_toy_params(R_max = 3000, K_half = 60,
                        shade_dose_response = list(max_effect = 0.5, form = "saturating",
                                                    half_sat = 80, hill = 1))
  pop <- make_toy_pop(100, p)
  pop$flowering <- TRUE
  gt  <- make_toy_gt(100, p)

  n_flowering <- sum(pop$alive & pop$flowering)
  N_canopy    <- canopy_density(pop)
  expected    <- bevholt(n_flowering, p) * (1 - shade_suppression(N_canopy, p))

  counts <- replicate(200, recruit(pop, gt, t = 1L, params = p)$n_rec)
  expect_equal(mean(counts), expected, tolerance = 0.05)
})

test_that("recruit(): shading strictly suppresses the recruit count relative to shading off", {
  set.seed(13L)
  p_off <- make_toy_params(R_max = 2000, K_half = 50)
  p_on  <- p_off
  p_on$shade_dose_response <- list(max_effect = 0.7, form = "saturating", half_sat = 50, hill = 1)

  pop <- make_toy_pop(80, p_off)
  pop$flowering <- TRUE
  gt  <- make_toy_gt(80, p_off)

  counts_off <- replicate(200, recruit(pop, gt, t = 1L, params = p_off)$n_rec)
  counts_on  <- replicate(200, recruit(pop, gt, t = 1L, params = p_on)$n_rec)

  expect_lt(mean(counts_on), mean(counts_off))
})
