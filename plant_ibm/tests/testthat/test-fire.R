# tests/testthat/test-fire.R
#
# Tests cover BOTH modes of apply_fire():
#   Legacy mode (fire_p_fimp = 0): flat fire_kill_prob, backward-compatible
#   Two-stage mode (fire_p_fimp > 0): age-dependent kill/resprout split

# ===========================================================================
# Legacy mode tests (fire_p_fimp = 0, the default)
# ===========================================================================

test_that("legacy: fire_kill_prob=0 -- nobody dies, everyone resproutes", {
  p   <- make_toy_params(fire_kill_prob = 0, fire_p_fimp = 0)
  pop <- make_toy_pop(300, p)
  n_alive_before <- sum(pop$alive)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_equal(sum(pop2$alive), n_alive_before)
  expect_true(all(pop2$resprout[pop2$alive]))
  expect_true(all(!pop2$flowering[pop2$alive]))
})

test_that("legacy: fire_kill_prob=1 -- everyone dies", {
  p   <- make_toy_params(fire_kill_prob = 1, fire_p_fimp = 0)
  pop <- make_toy_pop(200, p)
  pop2 <- apply_fire(pop, p, t = 5L)
  expect_equal(sum(pop2$alive), 0L)
  expect_true(all(!pop2$resprout))
})

test_that("legacy: kills approximately fire_kill_prob fraction (statistical)", {
  set.seed(77L)
  fkp <- 0.40
  p   <- make_toy_params(fire_kill_prob = fkp, fire_p_fimp = 0)
  pop <- make_toy_pop(5000, p)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_equal(1 - sum(pop2$alive)/sum(pop$alive), fkp, tolerance = 0.03)
})

test_that("legacy: fire sets death_year for killed individuals", {
  set.seed(3L)
  p   <- make_toy_params(fire_kill_prob = 0.5, fire_p_fimp = 0)
  pop <- make_toy_pop(200, p)
  pop2 <- apply_fire(pop, p, t = 42L)
  killed <- which(!pop2$alive & pop$alive)
  expect_true(all(pop2$death_year[killed] == 42L))
})

test_that("legacy: survivors have resprout_yrs_remain >= 1", {
  set.seed(9L)
  p   <- make_toy_params(fire_kill_prob = 0.3, resprout_yrs_base = 3L, fire_p_fimp = 0)
  pop <- make_toy_pop(500, p)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_true(all(pop2$resprout_yrs_remain[pop2$alive] >= 1L))
})

test_that("legacy: pre-existing dead individuals are unaffected", {
  set.seed(21L)
  p   <- make_toy_params(fire_kill_prob = 0.5, fire_p_fimp = 0)
  pop <- make_toy_pop(200, p)
  pop$alive[1:50] <- FALSE; pop$death_year[1:50] <- 0L
  pop2 <- apply_fire(pop, p, t = 10L)
  expect_true(all(pop2$death_year[1:50] == 0L))
})

test_that("legacy: rust before rust_start_year does not extend resprout delay", {
  set.seed(31L)
  p <- make_toy_params(fire_kill_prob = 0, resprout_yrs_base = 3L,
                        rust_start_year = 100L, rust_pressure = 1.0, fire_p_fimp = 0)
  p$rust_dose_response$delay$max_effect <- 5
  pop <- make_toy_pop(1000, p)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_lt(mean(pop2$resprout_yrs_remain[pop2$alive]), 6)
})

test_that("apply_fire with empty population returns pop unchanged", {
  p   <- make_toy_params()
  pop <- make_toy_pop(100, p)
  pop$alive[] <- FALSE
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_equal(sum(pop2$alive), 0L)
})

# ===========================================================================
# Two-stage mode tests (fire_p_fimp > 0)
# ===========================================================================

make_fire2 <- function(...) {
  make_toy_params(fire_p_fimp = 0.8, fire_kill_scalar = 0.3,
                   fire_kill_half_sat = 3, fire_kill_hill = 2, ...)
}

test_that("two-stage: unimpacted plants are completely unaffected", {
  set.seed(55L)
  p   <- make_toy_params(fire_p_fimp = 0.0001,  # almost nobody impacted
                          fire_kill_scalar = 0.5,
                          fire_kill_half_sat = 3, fire_kill_hill = 2)
  pop <- make_toy_pop(5000, p)
  pop2 <- apply_fire(pop, p, t = 1L)
  # With p_impact = 0.0001, expect ≈ 0.5 impacted: essentially nobody dead
  expect_lt(sum(!pop2$alive), 10)
})

test_that("two-stage: seedlings (age 0) are always killed when impacted", {
  # p_kill(0) = ceiling + (1-ceiling)*1 = 1 by formula
  set.seed(57L)
  p   <- make_fire2()
  pop <- make_toy_pop(1000, p)
  pop$age <- 0L   # force everyone to age 0
  pop2 <- apply_fire(pop, p, t = 1L)
  impacted_survivors <- sum(pop2$alive & pop2$resprout)
  expect_equal(impacted_survivors, 0L)   # no seedlings survive to resprout
})

test_that("two-stage: older plants have lower kill probability than seedlings", {
  set.seed(59L)
  p   <- make_fire2()
  # Force p_fimp=1 (all impacted) to isolate Stage 2
  p$fire_p_fimp <- 1

  pop_young <- make_toy_pop(3000, p); pop_young$age <- 1L
  pop_old   <- make_toy_pop(3000, p); pop_old$age   <- 40L

  pop_y2 <- apply_fire(pop_young, p, t = 1L)
  pop_o2 <- apply_fire(pop_old,   p, t = 1L)

  kill_young <- 1 - sum(pop_y2$alive) / 3000
  kill_old   <- 1 - sum(pop_o2$alive) / 3000
  expect_gt(kill_young, kill_old)
})

test_that("two-stage: old-plant kill ceiling scales with fire_p_fimp (intensity)", {
  # At very old age, p_kill → c * p_fimp. Hotter fire → more old plants killed.
  set.seed(61L)
  n <- 5000
  make_intense_fire <- function(pfimp) {
    p <- make_toy_params(fire_p_fimp = pfimp, fire_kill_scalar = 0.4,
                          fire_kill_half_sat = 1, fire_kill_hill = 10)
    pop <- make_toy_pop(n, p); pop$age <- 100L
    pop2 <- apply_fire(pop, p, t = 1L)
    1 - sum(pop2$alive) / n
  }
  expect_gt(make_intense_fire(0.9), make_intense_fire(0.3))
})

test_that("two-stage: resprout_total_yrs recorded for survivors", {
  set.seed(63L)
  p   <- make_fire2()
  pop <- make_toy_pop(500, p); pop$age <- 20L
  pop2 <- apply_fire(pop, p, t = 1L)
  surv <- pop2$alive & pop2$resprout
  if (any(surv)) {
    expect_true(all(pop2$resprout_total_yrs[surv] >= 1L))
    expect_true(all(pop2$resprout_yrs_remain[surv] == pop2$resprout_total_yrs[surv]))
  }
})

# ============================================================================
# Stand extent mode (fire_extent_mode = "stand"; notes D3)
# ============================================================================
.stand_params <- function(...) {
  p <- make_toy_params(fire_p_fimp = 0.3, fire_kill_scalar = 0.3,
                       fire_kill_half_sat = 2, fire_kill_hill = 2)
  p$fire_extent_mode <- "stand"; p$fire_size_sd <- 0; p$fire_stand_cover <- 1
  modifyList(p, list(...))
}

test_that("stand: a fire either misses the stand entirely or reaches it", {
  p <- .stand_params()
  pop <- make_toy_pop(500, p)
  set.seed(11)
  res <- replicate(400, {
    pop2 <- apply_fire(pop, p, t = 1L)
    c(missed = identical(attr(pop2, "fire_reached"), FALSE),
      touched = sum(pop2$resprout | !pop2$alive))
  })
  missed <- res["missed", ] == 1
  expect_true(all(res["touched", missed] == 0))      # misses touch nobody
  expect_true(all(res["touched", !missed] == 500))   # cover 1: everyone burns
  expect_equal(mean(!missed), 0.3, tolerance = 0.06) # reached with P = f = 0.3
})

test_that("stand: fire_stand_cover is the share burnt when a fire arrives", {
  p <- .stand_params(fire_p_fimp = 1, fire_stand_cover = 0.9)
  pop <- make_toy_pop(4000, p)
  set.seed(3)
  pop2 <- apply_fire(pop, p, t = 1L)
  expect_equal(mean(pop2$resprout | !pop2$alive), 0.9, tolerance = 0.02)
})

test_that("stand: old-plant kill ceiling scales with the drawn fire size", {
  old <- function(p, n = 6000) { pop <- make_toy_pop(n, p); pop$age <- 200L; pop }
  set.seed(5)
  p1 <- .stand_params(fire_p_fimp = 1)
  k1 <- mean(!apply_fire(old(p1), p1, t = 1L)$alive)
  expect_equal(k1, 0.3, tolerance = 0.02)            # f = 1: ceiling = c
  # Small fires that do reach the stand are gentler: kill share ~ c * f.
  p2 <- .stand_params(fire_p_fimp = 0.25)
  pop_old <- old(p2, 2000)              # built once: make_toy_pop() resets the seed
  ks <- replicate(200, {
    pop2 <- apply_fire(pop_old, p2, t = 1L)
    if (identical(attr(pop2, "fire_reached"), FALSE)) NA else mean(!pop2$alive)
  })
  expect_equal(mean(ks, na.rm = TRUE), 0.3 * 0.25, tolerance = 0.015)
})

test_that("stand: fire sizes vary between fires with the cell's sd", {
  set.seed(9)
  f <- replicate(4000, draw_fire_size(0.25, 0.12))
  expect_equal(mean(f), 0.25, tolerance = 0.01)
  expect_equal(sd(f), 0.12, tolerance = 0.01)
  expect_true(all(f > 0 & f < 1))
  expect_equal(draw_fire_size(0.25, 0), 0.25)
  expect_equal(draw_fire_size(1, 0.3), 1)
  big <- replicate(500, draw_fire_size(0.2, 0.9))      # sd above the Beta maximum: capped
  expect_true(all(is.finite(big)))
})

test_that("plant mode is the original model (ceiling = c * fire_p_fimp)", {
  set.seed(5)
  p <- .stand_params(fire_p_fimp = 1); p$fire_extent_mode <- "plant"
  pop <- make_toy_pop(6000, p); pop$age <- 200L
  expect_equal(mean(!apply_fire(pop, p, t = 1L)$alive), 0.3, tolerance = 0.02)
  p$fire_extent_mode <- NULL                          # absent = plant
  expect_null(attr(apply_fire(pop, p, t = 1L), "fire_reached"))
})

test_that("census separates fires in the cell from fires reaching the stand", {
  p <- .stand_params(fire_p_fimp = 0.05)
  p$N0 <- 200L; p$n_years <- 60L; p$fire_years <- integer(0); p$fire_prob_annual <- 0.5
  cen <- run_simulation(p, year0 = 1L)$census
  expect_true(all(cen$fire_event <= cen$fire_in_cell))
  expect_gt(sum(cen$fire_in_cell), sum(cen$fire_event))
})
