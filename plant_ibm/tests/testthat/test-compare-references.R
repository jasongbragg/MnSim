# tests/testthat/test-compare-references.R
#
# The reference-comparison functions (tests/regression/compare.R) used in
# step 1 of the reference-update procedure (docs/reference_updates.md).
# Builds altered copies of the real reference set in memory; never writes.

source(file.path(proj_root, "tests", "regression", "compare.R"), local = TRUE)

ref_dir <- Sys.getenv("MNSIM_REFERENCE_DIR",
                      normalizePath(file.path(proj_root, "..", "tests", "reference"), mustWork = FALSE))

test_that("a reference set compared with itself shows nothing to update", {
  set.seed(1)
  skip_if_not(dir.exists(ref_dir), "no reference directory")
  ref <- load_reference_set(ref_dir)
  cmp <- compare_reference_sets(ref, ref)
  expect_equal(nrow(cmp$params), 0)
  expect_true(all(cmp$output_differs == ""))
  expect_true(all(cmp$stats$verdict == "within"))
  expect_match(classify_comparison(cmp), "nothing to update", all = FALSE)
  expect_true(length(report_comparison(cmp)) > 10)
})

test_that("a changed parameter and output are reported as a calibration change", {
  set.seed(1)
  skip_if_not(dir.exists(ref_dir), "no reference directory")
  old <- load_reference_set(ref_dir)
  new <- old
  sc  <- names(new$scenarios)[1]
  new$scenarios[[sc]]$params$weibull_lambda <- 55
  new$scenarios[[sc]]$output$census$N_alive[3] <- new$scenarios[[sc]]$output$census$N_alive[3] + 1L
  cmp <- compare_reference_sets(old, new)
  expect_identical(cmp$params$parameter, "weibull_lambda")
  expect_identical(cmp$params$scenarios, sc)
  expect_match(cmp$output_differs[[sc]], "N_alive: first difference at row 3")
  expect_match(classify_comparison(cmp), "row 1", all = FALSE)
})

test_that("same params, changed output: statistics decide between rows 2 and 3", {
  set.seed(1)
  skip_if_not(dir.exists(ref_dir), "no reference directory")
  old <- load_reference_set(ref_dir)
  skip_if(is.null(old$stats), "no summary_stats.rds")
  new <- old
  sc  <- names(new$scenarios)[1]
  new$scenarios[[sc]]$output$census$deaths[5] <- new$scenarios[[sc]]$output$census$deaths[5] + 1L
  expect_match(classify_comparison(compare_reference_sets(old, new)), "row 3", all = FALSE)

  st <- names(new$stats$stats)[1]
  new$stats$stats[[st]]$mean <- old$stats$stats[[st]]$mean + 10 * old$stats$stats[[st]]$sd
  cmp <- compare_reference_sets(old, new)
  expect_identical(cmp$stats$verdict[cmp$stats$statistic == st], "OUTSIDE")
  expect_match(classify_comparison(cmp), "row 2", all = FALSE)
})
