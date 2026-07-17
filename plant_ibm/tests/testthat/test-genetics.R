# tests/testthat/test-genetics.R

# Deterministic edge cases -----------------------------------------------

test_that("SS x SS cross always produces SS offspring", {
  gt <- matrix(c(0L, 0L), nrow = 2L)
  rec <- inherit_alleles(gt, mother_idx = rep(1L, 500L),
                             father_idx = rep(2L, 500L))
  expect_true(all(rec == 0L))
})

test_that("RR x RR cross always produces RR offspring", {
  gt <- matrix(c(2L, 2L), nrow = 2L)
  rec <- inherit_alleles(gt, mother_idx = rep(1L, 500L),
                             father_idx = rep(2L, 500L))
  expect_true(all(rec == 2L))
})

test_that("resist_score = 0 for SS under any architecture", {
  p1 <- make_toy_params()
  p2 <- make_toy_params(resist_locus_effect = c(0.4, 0.35, 0.25),
                         resist_dominance     = c(0.5, 0.5, 0.5))
  for (p in list(p1, p2)) {
    n_loci <- length(p$resist_locus_effect)
    gt_ss  <- matrix(0L, nrow = 1L, ncol = n_loci)
    expect_equal(resist_score_from_gt(gt_ss, p), 0)
  }
})

test_that("resist_score = 1 for RR under single dominant locus", {
  p  <- make_toy_params(resist_locus_effect = c(1.0), resist_dominance = c(1.0))
  gt <- matrix(2L, nrow = 1L, ncol = 1L)
  expect_equal(resist_score_from_gt(gt, p), 1)
})

test_that("resist_score for RS equals effect*dominance for single locus", {
  eff <- 0.7; dom <- 0.4
  p   <- make_toy_params(resist_locus_effect = c(eff), resist_dominance = c(dom))
  gt  <- matrix(1L, nrow = 1L, ncol = 1L)
  expect_equal(resist_score_from_gt(gt, p), eff * dom)
})

test_that("resist_score is clamped at 1 even when multi-locus effects sum > 1", {
  p  <- make_toy_params(resist_locus_effect = c(0.6, 0.6),
                         resist_dominance     = c(1.0, 1.0))
  gt <- matrix(2L, nrow = 1L, ncol = 2L)   # RR at both loci
  expect_equal(resist_score_from_gt(gt, p), 1)
})

test_that("full dominance: RS score equals RR score", {
  p    <- make_toy_params(resist_locus_effect = c(1.0), resist_dominance = c(1.0))
  gt_h <- matrix(1L, 1L, 1L); gt_r <- matrix(2L, 1L, 1L)
  expect_equal(resist_score_from_gt(gt_h, p), resist_score_from_gt(gt_r, p))
})

test_that("full recessiveness: RS score equals SS score", {
  p    <- make_toy_params(resist_locus_effect = c(1.0), resist_dominance = c(0.0))
  gt_h <- matrix(1L, 1L, 1L); gt_s <- matrix(0L, 1L, 1L)
  expect_equal(resist_score_from_gt(gt_h, p), resist_score_from_gt(gt_s, p))
})

# Statistical checks (stochastic, large N for stability) -----------------

test_that("RS x RS cross gives ~1:2:1 Mendelian ratio", {
  set.seed(99L)
  gt  <- matrix(c(1L, 1L), nrow = 2L)
  rec <- inherit_alleles(gt, mother_idx = rep(1L, 40000L),
                             father_idx = rep(2L, 40000L))
  freq <- table(rec) / length(rec)
  expect_equal(as.numeric(freq["0"]), 0.25, tolerance = 0.02)
  expect_equal(as.numeric(freq["1"]), 0.50, tolerance = 0.02)
  expect_equal(as.numeric(freq["2"]), 0.25, tolerance = 0.02)
})

test_that("init_resist_gt produces allele frequencies near resist_freq0", {
  set.seed(11L)
  freq0 <- 0.20
  p     <- make_toy_params(resist_freq0 = c(freq0))
  gt    <- init_resist_gt(20000L, p)
  # Allele frequency = mean dose / 2
  obs   <- mean(gt[, 1L]) / 2
  expect_equal(obs, freq0, tolerance = 0.015)
})

test_that("allele_freqs returns NA when no alive individuals", {
  gt   <- matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L)
  alive <- c(FALSE, FALSE, FALSE)
  expect_true(is.na(allele_freqs(gt, alive)))
})

# Structural invariants --------------------------------------------------

test_that("inherit_alleles output has same ncol as input resist_gt", {
  set.seed(5L)
  p  <- make_toy_params(resist_locus_effect = c(0.4, 0.35),
                         resist_dominance     = c(0.5, 0.5),
                         resist_freq0         = c(0.1, 0.05))
  gt <- init_resist_gt(50L, p)
  expect_equal(ncol(gt), 2L)
  rec <- inherit_alleles(gt, mother_idx = 1:10, father_idx = 11:20)
  expect_equal(ncol(rec), 2L)
})

test_that("all offspring allele doses are in {0, 1, 2}", {
  set.seed(3L)
  p   <- make_toy_params()
  gt  <- init_resist_gt(100L, p)
  rec <- inherit_alleles(gt, mother_idx = sample(100L, 200L, TRUE),
                              father_idx = sample(100L, 200L, TRUE))
  expect_true(all(rec %in% 0:2))
})
