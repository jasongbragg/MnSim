# tests/testthat/test-iucn-a3.R
#
# Unit tests for runner/iucn_a3.R using entirely invented small-population
# data. No model simulation is run -- all inputs are hand-crafted so the
# expected outputs can be verified by arithmetic.
#
# compute_gen_time() and summarise_iucn_a3() are sourced via
# helper-fixtures.R (runner/iucn_a3.R), consistent with all other test files.

# ============================================================================
# Helper: build a minimal pop data frame for gen-time tests
# ============================================================================
#
# Scenario (year_spinup = 100, gen-time window = [200, 300]):
#
# id | type    | birth_year_stored | cal_birth | alive | death_year
# ---+---------+-------------------+-----------+-------+-----------
#  1 | FOUNDER | -10               |  90       | FALSE | 250    <- born before sim
#  2 | recruit |  150              | 150       | FALSE | 280    <- born during sim
#  3 | recruit |  180              | 180       | TRUE  | NA     <- still alive
#  4 | recruit |  210              | 210       | FALSE | 290    <- offspring of M1 (founder)
#  5 | recruit |  240              | 240       | FALSE | 295    <- offspring of M1 (founder)
#  6 | recruit |  260              | 260       | FALSE | 310    <- offspring of M2
#  7 | recruit |  270              | 270       | FALSE | 315    <- offspring of M3 (alive)
#  8 | recruit |  280              | 280       | FALSE | 320    <- offspring of M2
#  9 | recruit |  350              | 350       | FALSE | 400    <- outside window
#
# TWO FILTERS applied:
#   (A) Dead-mother filter: exclude id=7 (mother M3 is alive)
#   (B) Born-in-sim filter: exclude id=4, id=5 (mother M1 is a founder,
#       cal_birth=90 < year_spinup=100)
#
# Eligible offspring after both filters: id=6 and id=8 only
#   id=6: mother M2 born 150, offspring born 260, maternal age = 260-150 = 110
#   id=8: mother M2 born 150, offspring born 280, maternal age = 280-150 = 130
#
# Expected T_g = mean(110, 130) = 120
.make_test_pop <- function() {
  data.frame(
    id         = 1:9,
    birth_year = c(-10, 150, 180, 210, 240, 260, 270, 280, 350),
    alive      = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    death_year = c(250L, 280L, NA_integer_, 290L, 295L, 310L, 315L, 320L, 400L),
    mother_id  = c(NA, NA, NA, 1L, 1L, 2L, 3L, 2L, 1L),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# compute_gen_time tests
# ============================================================================

test_that("compute_gen_time returns expected T_g = 120 after both filters", {
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_equal(gt$T_g, 120, tolerance = 1e-10)
})

test_that("compute_gen_time counts exactly 2 eligible offspring after both filters", {
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_equal(gt$n_offspring, 2L)
})

test_that("compute_gen_time returns correct individual maternal ages", {
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_setequal(gt$mother_ages, c(110, 130))
})

test_that("compute_gen_time correctly adjusts founder birth_year using year_spinup", {
  # Verify the adjustment works: M1 stored birth_year=-10, year_spinup=100 →
  # cal_birth=90. But M1 is a founder, so her offspring are excluded by filter B.
  # If the adjustment were WRONG (unadjusted -10 used), maternal ages would be
  # 210-(-10)=220 and 240-(-10)=250. Neither should appear.
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_false(220 %in% gt$mother_ages)
  expect_false(250 %in% gt$mother_ages)
})

test_that("compute_gen_time FILTER B: excludes offspring of founder mothers", {
  # M1 is a founder (cal_birth=90 < year_spinup=100).
  # id=4 and id=5 are M1's offspring; their maternal ages (120 and 150)
  # must NOT appear even though M1 is dead.
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_false(120 %in% gt$mother_ages)  # id=4: 210-90=120, should be excluded
  expect_false(150 %in% gt$mother_ages)  # id=5: 240-90=150, should be excluded
})

test_that("compute_gen_time FILTER A: excludes offspring of living mothers", {
  # M3 is alive (alive=TRUE). id=7 is M3's offspring; maternal age would be
  # 270-180=90, but M3 is alive so it must be excluded.
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_false(90 %in% gt$mother_ages)
})

test_that("compute_gen_time excludes offspring born outside the window", {
  # id=9 is born at 350, outside window [200,300].
  # Mother is M1 (founder, excluded anyway), maternal age would be 350-90=260.
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_false(260 %in% gt$mother_ages)
})

test_that("compute_gen_time includes only the two recruit-mother offspring", {
  # Final positive check: 110 and 130 are present, nothing else.
  pop <- .make_test_pop()
  gt  <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100)
  expect_true(110 %in% gt$mother_ages)
  expect_true(130 %in% gt$mother_ages)
  expect_equal(length(gt$mother_ages), 2L)
})

test_that("compute_gen_time warns and returns NA when no offspring in window", {
  pop <- .make_test_pop()
  expect_warning(
    gt <- compute_gen_time(pop, window_start = 500, window_end = 600, year_spinup = 100),
    "No offspring with known mothers"
  )
  expect_true(is.na(gt$T_g))
  expect_equal(gt$n_offspring, 0L)
})

test_that("compute_gen_time warns and returns NA when no eligible mothers remain", {
  # Force all mothers to be alive so both filters together leave nothing.
  pop <- .make_test_pop()
  pop$alive[] <- TRUE
  expect_warning(
    gt <- compute_gen_time(pop, window_start = 200, window_end = 300, year_spinup = 100),
    "No eligible offspring"
  )
  expect_true(is.na(gt$T_g))
})

test_that("compute_gen_time works correctly with recruit-only mothers (no founders)", {
  # Simple case: year_spinup=0, all individuals have birth_year >= 0 (recruits).
  pop <- data.frame(
    id         = 1:4,
    birth_year = c(10L, 20L, 30L, 40L),
    alive      = c(FALSE, FALSE, FALSE, FALSE),
    death_year = c(50L, 60L, 70L, 80L),
    mother_id  = c(NA_integer_, NA_integer_, 1L, 2L),
    stringsAsFactors = FALSE
  )
  # id=3: offspring of id=1 (born 10), born 30 → maternal age = 20
  # id=4: offspring of id=2 (born 20), born 40 → maternal age = 20
  gt <- compute_gen_time(pop, window_start = 0, window_end = 50, year_spinup = 0)
  expect_equal(gt$T_g, 20, tolerance = 1e-10)
  expect_equal(gt$n_offspring, 2L)
})

# ============================================================================
# summarise_iucn_a3 tests
# ============================================================================

.make_mock_result <- function(pct_decline, T_g = 25, N_ref = 1000) {
  N_compare <- N_ref * (1 - pct_decline / 100)
  t_inf     <- min(100, 3 * T_g)
  list(
    run_id          = "test123abc456",
    T_g             = T_g,
    t_inf           = t_inf,
    N_ref           = N_ref,
    N_compare       = N_compare,
    pct_decline     = pct_decline,
    iucn_cat        = if (pct_decline >= 80) "CR" else if (pct_decline >= 50) "EN"
                      else if (pct_decline >= 30) "VU" else "LC/NT",
    n_offspring_gt  = 42L,
    compare_mid     = 2085L,
    mother_ages     = c(20, 25, 30),
    year_rust       = 2010L,
    year_spinup     = 1700L,
    year_gentime_start = 1810L,
    year_gentime_end   = 1910L,
    year_end        = 2210L,
    ref_start       = 2004L, ref_end     = 2009L,
    compare_start   = 2080L, compare_end = 2090L,
    census          = data.frame(year = integer(0), N_IUCN = integer(0))
  )
}

test_that("summarise_iucn_a3 returns a single-row data.frame", {
  s <- summarise_iucn_a3(.make_mock_result(40))
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
})

test_that("summarise_iucn_a3 preserves numeric types for T_g and pct_decline", {
  s <- summarise_iucn_a3(.make_mock_result(40))
  expect_type(s$T_g,         "double")
  expect_type(s$pct_decline, "double")
})

test_that("summarise_iucn_a3 preserves run_id as character (not coerced numeric)", {
  s <- summarise_iucn_a3(.make_mock_result(40))
  expect_type(s$run_id, "character")
  expect_equal(s$run_id, "test123abc456")
})

test_that("summarise_iucn_a3 VU boundary: exactly 30% triggers VU only", {
  s <- summarise_iucn_a3(.make_mock_result(30))
  expect_equal(s$iucn_vu, 1L)
  expect_equal(s$iucn_en, 0L)
  expect_equal(s$iucn_cr, 0L)
})

test_that("summarise_iucn_a3 below VU: 29.9% gives iucn_vu = 0", {
  s <- summarise_iucn_a3(.make_mock_result(29.9))
  expect_equal(s$iucn_vu, 0L)
})

test_that("summarise_iucn_a3 EN boundary: 50% triggers VU and EN", {
  s <- summarise_iucn_a3(.make_mock_result(50))
  expect_equal(s$iucn_vu, 1L)
  expect_equal(s$iucn_en, 1L)
  expect_equal(s$iucn_cr, 0L)
})

test_that("summarise_iucn_a3 CR boundary: 80% triggers all three thresholds", {
  s <- summarise_iucn_a3(.make_mock_result(80))
  expect_equal(s$iucn_vu, 1L)
  expect_equal(s$iucn_en, 1L)
  expect_equal(s$iucn_cr, 1L)
})

test_that("summarise_iucn_a3 t_inf = min(100, 3*T_g): short generation", {
  # T_g=25: 3*25=75 < 100, so t_inf=75
  s <- summarise_iucn_a3(.make_mock_result(40, T_g = 25))
  expect_equal(s$t_inf, 75)
})

test_that("summarise_iucn_a3 t_inf = min(100, 3*T_g): long generation caps at 100", {
  # T_g=40: 3*40=120 > 100, so t_inf=100
  s <- summarise_iucn_a3(.make_mock_result(40, T_g = 40))
  expect_equal(s$t_inf, 100)
})

test_that("summarise_iucn_a3 N_ref and N_compare are numerically correct", {
  s <- summarise_iucn_a3(.make_mock_result(40, N_ref = 2000))
  expect_equal(s$N_ref,      2000, tolerance = 1e-10)
  expect_equal(s$N_compare,  1200, tolerance = 1e-10)  # 2000*(1-0.40)
  expect_equal(s$pct_decline,  40, tolerance = 1e-10)
})
