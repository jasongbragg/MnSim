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
# Fixture: N = 1000 up to 2010, then a straight-line decline to 400 at 2080,
# then flat at 400. Rust arrives in 2010, so the reference window is the 11 yr
# ending 2009 (1999-2009, mean 1000). With T_g = 70/3, 3*T_g = 70 yr and the
# comparison window is the 11 yr ending 2009 + 70 = 2079 (2069-2079). N is
# linear there, so its mean is N at the mid-year 2074:
#   1000 - 600 * 64 / 70 = 451.43  ->  decline 54.86%.
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
.N_at <- function(y) ifelse(y <= 2010, 1000, ifelse(y >= 2080, 400, 1000 - 600 * (y - 2010) / 70))

test_that("assess_a4 compares the pre-rust window with the window 3*T_g later", {
  a <- assess_a4(.make_census(), T_g = 70 / 3, year_rust = 2010L, year_assess = 2026L)
  expect_equal(a$window_len, 70L)
  expect_equal(c(a$ref_start, a$ref_end), c(1999L, 2009L))
  expect_equal(c(a$cmp_start, a$cmp_end), c(2069L, 2079L))
  expect_equal(a$N_ref, 1000)
  expect_equal(a$N_compare, mean(.N_at(2069:2079)))
  expect_equal(a$pct_decline, 100 * (1000 - mean(.N_at(2069:2079))) / 1000)
  expect_equal(a$iucn_cat, "EN")
})

test_that("assess_a4 window is fixed by rust arrival, not by the present", {
  a1 <- assess_a4(.make_census(), T_g = 70 / 3, year_rust = 2010L, year_assess = 2026L)
  a2 <- assess_a4(.make_census(), T_g = 70 / 3, year_rust = 2010L, year_assess = 2040L)
  expect_identical(a1[c("ref_end", "cmp_end", "pct_decline")],
                   a2[c("ref_end", "cmp_end", "pct_decline")])
})

test_that("assess_a4 averages over avg_years at each end", {
  a <- assess_a4(.make_census(), T_g = 70 / 3, year_rust = 2010L,
                 year_assess = 2026L, avg_years = 1L)
  expect_equal(c(a$ref_start, a$ref_end, a$cmp_start, a$cmp_end),
               c(2009L, 2009L, 2079L, 2079L))
  expect_equal(a$N_compare, .N_at(2079))
})

test_that("assess_a4 uses a 10-year minimum window (and warns it is not A4)", {
  expect_warning(
    a <- assess_a4(.make_census(), T_g = 2, year_rust = 2010L, year_assess = 2026L),
    "not an A4 window")
  expect_equal(a$window_len, 10L)
  expect_equal(a$cmp_end, 2019L)
})

test_that("assess_a4 caps the comparison window at 100 yr after the present", {
  expect_warning(
    a <- assess_a4(.make_census(), T_g = 40, year_rust = 2010L, year_assess = 2026L),
    "capping")
  expect_equal(a$cmp_end, 2126L)
  expect_equal(a$window_len, 2126L - 2009L)
  expect_equal(a$pct_decline, 60)
})

test_that("assess_a4 treats years after extinction as N = 0", {
  a <- assess_a4(.make_census(extinct_at = 2050L), T_g = 70 / 3,
                 year_rust = 2010L, year_assess = 2026L)
  expect_equal(a$N_compare, 0)
  expect_equal(a$pct_decline, 100)
  expect_equal(a$iucn_cat, "CR")
})

test_that("assess_a4 warns and returns NA when the census is too short", {
  expect_warning(
    a <- assess_a4(.make_census(1900:2050), T_g = 70 / 3,
                   year_rust = 2010L, year_assess = 2026L),
    "does not cover")
  expect_true(is.na(a$pct_decline))
  expect_true(is.na(a$cmp_end))
})

test_that("iucn_category thresholds", {
  expect_equal(iucn_category(NA),   NA_character_)
  expect_equal(iucn_category(29.9), "LC/NT")
  expect_equal(iucn_category(30),   "VU")
  expect_equal(iucn_category(50),   "EN")
  expect_equal(iucn_category(80),   "CR")
})

test_that("iucn_run_id changes when assessment settings change", {
  p  <- get_default_params()
  id1 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2026L))
  id2 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2030L))
  id3 <- iucn_run_id(p, iucn_assess_spec(year_assess = 2026L))
  id4 <- iucn_run_id(p, iucn_assess_spec(avg_years = 5L))
  expect_false(identical(id1, id2))
  expect_false(identical(id1, id4))
  expect_identical(id1, id3)
  expect_false(identical(id1, make_run_id(p)))
})

# ============================================================================
# summarise_iucn_a4 tests
# ============================================================================

.make_mock_result <- function(pct_decline, T_g = 25, N_ref = 1000) {
  L <- as.integer(round(max(10, 3 * T_g)))
  list(
    criterion      = "A4",
    run_id         = "test123abc456",
    T_g            = T_g,
    window_len     = L,
    year_assess    = 2026L,
    ref_start      = 1999L,
    ref_end        = 2009L,
    cmp_start      = 2009L + L - 10L,
    cmp_end        = 2009L + L,
    N_ref          = N_ref,
    N_compare      = N_ref * (1 - pct_decline / 100),
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

test_that("summarise_iucn_a4 carries window bounds and N values", {
  s <- summarise_iucn_a4(.make_mock_result(40, N_ref = 2000))
  expect_equal(c(s$ref_start, s$ref_end), c(1999L, 2009L))
  expect_equal(s$cmp_end,   2084L)
  expect_equal(s$N_ref,     2000, tolerance = 1e-10)
  expect_equal(s$N_compare, 1200, tolerance = 1e-10)
  expect_equal(s$pct_decline, 40, tolerance = 1e-10)
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

test_that("assess_a4 can assess the flowering-only count over the same window", {
  cen <- .make_census()
  cen$N_flowering <- cen$N_IUCN / 2
  cen$N_flowering[cen$year >= 2069] <- 100
  a  <- assess_a4(cen, T_g = 70 / 3, year_rust = 2010L, year_assess = 2026L)
  af <- assess_a4(cen, T_g = 70 / 3, year_rust = 2010L, year_assess = 2026L,
                  n_col = "N_flowering")
  expect_equal(af$cmp_end, a$cmp_end)
  expect_equal(af$N_ref, 500)
  expect_equal(af$pct_decline, 80)
  expect_error(assess_a4(cen, 20, 2010L, 2026L, n_col = "nope"), "no column")
})
