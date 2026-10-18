# tests/testthat/test-iucn-a4.R
#
# Unit tests for runner/iucn_a4.R using entirely invented small-population
# data. No model simulation is run -- all inputs are hand-crafted so the
# expected outputs can be verified by arithmetic.
#
# compute_gen_time() and summarise_iucn_a4() and assess_a4() are sourced via
# helper-fixtures.R (runner/iucn_a4.R), consistent with all other test files.

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
# assess_a4 tests (hand-built census trajectories, no simulation)
# ============================================================================
#
# Reference trajectory: N = 1000 up to 2010, then a straight-line decline to
# 400 at 2080, then flat at 400. With T_g = 70/3 the window is L = 70 yr and,
# with no smoothing (smooth_half = 0), the worst window is 2010-2080:
#   decline = (1000 - 400) / 1000 = 60%.
# Windows starting before 2010 end before 2080 (less decline); windows
# starting after 2010 start lower but end at the same 400 floor.
.make_census <- function(years = 1900:2220, extinct_at = NA) {
  N <- ifelse(years <= 2010, 1000,
       ifelse(years >= 2080, 400, 1000 - 600 * (years - 2010) / 70))
  cen <- data.frame(year = years, N_IUCN = N, extinct = FALSE)
  if (!is.na(extinct_at)) {
    cen <- cen[cen$year <= extinct_at, ]
    cen$N_IUCN[cen$year == extinct_at] <- 0
    cen$extinct[cen$year == extinct_at] <- TRUE
  }
  cen
}

test_that("assess_a4 finds the worst window and its decline", {
  a <- assess_a4(.make_census(), T_g = 70 / 3, year_assess = 2026L, smooth_half = 0L)
  expect_equal(a$window_len, 70L)
  expect_equal(a$win_start, 2010L)
  expect_equal(a$win_end,   2080L)
  expect_equal(a$pct_decline, 60, tolerance = 1e-10)
  expect_equal(a$iucn_cat, "EN")
})

test_that("assess_a4 windows all straddle the present and reach <= 100 yr ahead", {
  a  <- assess_a4(.make_census(), T_g = 70 / 3, year_assess = 2026L, smooth_half = 0L)
  dc <- a$decline_curve
  expect_true(all(dc$win_start <  2026L))
  expect_true(all(dc$win_end   >  2026L))
  expect_true(all(dc$win_end   <= 2026L + 100L))
  expect_equal(range(dc$win_start), c(2026L - 70L + 1L, 2025L))
})

test_that("assess_a4 uses a 10-year minimum window for short generations", {
  a <- assess_a4(.make_census(), T_g = 2, year_assess = 2026L, smooth_half = 0L)
  expect_equal(a$window_len, 10L)
  # Steepest 10-yr decline in the linear phase: start as late as possible
  # (2025) because the per-year drop is constant but N(s) keeps falling.
  expect_equal(a$win_start, 2025L)
  N_s <- 1000 - 600 * 15 / 70
  expect_equal(a$pct_decline, 100 * (600 * 10 / 70) / N_s, tolerance = 1e-10)
})

test_that("assess_a4 caps the future at 100 yr when 3*T_g exceeds 100", {
  # T_g = 40 -> L = 120: windows must start in [2026-119, 2026+100-120].
  a  <- assess_a4(.make_census(), T_g = 40, year_assess = 2026L, smooth_half = 0L)
  expect_equal(a$window_len, 120L)
  expect_equal(max(a$decline_curve$win_start), 2006L)
  expect_true(all(a$decline_curve$win_end <= 2126L))
  expect_equal(a$pct_decline, 60, tolerance = 1e-10)
})

test_that("assess_a4 running mean is applied at both window ends", {
  a <- assess_a4(.make_census(), T_g = 70 / 3, year_assess = 2026L, smooth_half = 5L)
  dc <- a$decline_curve
  row <- dc[dc$win_start == 2015L, ]
  N_y <- function(y) ifelse(y <= 2010, 1000, ifelse(y >= 2080, 400, 1000 - 600 * (y - 2010) / 70))
  expect_equal(row$N_start, mean(N_y(2010:2020)), tolerance = 1e-10)
  expect_equal(row$N_end,   mean(N_y(2080:2090)), tolerance = 1e-10)
})

test_that("assess_a4 treats years after extinction as N = 0", {
  a <- assess_a4(.make_census(extinct_at = 2050L), T_g = 70 / 3,
                 year_assess = 2026L, smooth_half = 2L)
  expect_equal(a$pct_decline, 100, tolerance = 1e-10)
  expect_equal(a$iucn_cat, "CR")
})

test_that("assess_a4 warns and returns NA when the census is too short", {
  expect_warning(
    a <- assess_a4(.make_census(years = 1900:2020), T_g = 70 / 3,
                   year_assess = 2026L, smooth_half = 0L),
    "No A4 window")
  expect_true(is.na(a$pct_decline))
  expect_true(is.na(a$iucn_cat))
})

test_that("iucn_category thresholds", {
  expect_equal(iucn_category(29.9), "LC/NT")
  expect_equal(iucn_category(30),   "VU")
  expect_equal(iucn_category(50),   "EN")
  expect_equal(iucn_category(80),   "CR")
  expect_true(is.na(iucn_category(NA_real_)))
})

test_that("iucn_run_id changes when assessment settings change", {
  p  <- get_default_params()
  id1 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2026L))
  id2 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2030L))
  id3 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2026L))
  expect_false(identical(id1, id2))
  expect_identical(id1, id3)
  expect_false(identical(id1, make_run_id(p)))
})

# ============================================================================
# summarise_iucn_a4 tests
# ============================================================================

.make_mock_result <- function(pct_decline, T_g = 25, N_start = 1000) {
  list(
    criterion      = "A4",
    run_id         = "test123abc456",
    T_g            = T_g,
    window_len     = as.integer(round(max(10, 3 * T_g))),
    year_assess    = 2026L,
    win_start      = 2005L,
    win_end        = 2080L,
    N_start        = N_start,
    N_end          = N_start * (1 - pct_decline / 100),
    pct_decline    = pct_decline,
    iucn_cat       = iucn_category(pct_decline),
    n_offspring_gt = 42L
  )
}

test_that("summarise_iucn_a4 returns a single-row data.frame with criterion A4", {
  s <- summarise_iucn_a4(.make_mock_result(40))
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
  expect_equal(s$criterion, "A4")
})

test_that("summarise_iucn_a4 preserves types", {
  s <- summarise_iucn_a4(.make_mock_result(40))
  expect_type(s$T_g,         "double")
  expect_type(s$pct_decline, "double")
  expect_type(s$run_id,      "character")
  expect_equal(s$run_id, "test123abc456")
})

test_that("summarise_iucn_a4 threshold flags at the boundaries", {
  s <- summarise_iucn_a4(.make_mock_result(29.9))
  expect_equal(c(s$iucn_vu, s$iucn_en, s$iucn_cr), c(0L, 0L, 0L))
  s <- summarise_iucn_a4(.make_mock_result(30))
  expect_equal(c(s$iucn_vu, s$iucn_en, s$iucn_cr), c(1L, 0L, 0L))
  s <- summarise_iucn_a4(.make_mock_result(50))
  expect_equal(c(s$iucn_vu, s$iucn_en, s$iucn_cr), c(1L, 1L, 0L))
  s <- summarise_iucn_a4(.make_mock_result(80))
  expect_equal(c(s$iucn_vu, s$iucn_en, s$iucn_cr), c(1L, 1L, 1L))
})

test_that("summarise_iucn_a4 window_len = max(10, 3*T_g), no 100-yr cap", {
  expect_equal(summarise_iucn_a4(.make_mock_result(40, T_g = 25))$window_len, 75L)
  expect_equal(summarise_iucn_a4(.make_mock_result(40, T_g = 40))$window_len, 120L)
  expect_equal(summarise_iucn_a4(.make_mock_result(40, T_g = 2))$window_len,  10L)
})

test_that("summarise_iucn_a4 N_start and N_end are numerically correct", {
  s <- summarise_iucn_a4(.make_mock_result(40, N_start = 2000))
  expect_equal(s$N_start,     2000, tolerance = 1e-10)
  expect_equal(s$N_end,       1200, tolerance = 1e-10)
  expect_equal(s$pct_decline,   40, tolerance = 1e-10)
})

test_that("collect_iucn_results skips non-A4 result files", {
  d <- file.path(tempdir(), "a4_collect_test")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  saveRDS(.make_mock_result(40), file.path(d, "a4.rds"))
  old <- .make_mock_result(55); old$criterion <- NULL; old$run_id <- "oldA3run0000"
  saveRDS(old, file.path(d, "a3.rds"))
  df <- suppressMessages(collect_iucn_results(d))
  expect_equal(nrow(df), 1L)
  expect_equal(df$criterion, "A4")
})

test_that("assess_a4 warns when only some windows can be evaluated", {
  expect_warning(
    a <- assess_a4(.make_census(years = 1900:2030), T_g = 70 / 3,
                   year_assess = 2026L, smooth_half = 0L),
    "could not be evaluated")
  expect_equal(sum(!is.na(a$decline_curve$pct_decline)), 4L)
})
