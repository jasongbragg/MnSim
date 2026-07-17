# tests/testthat.R
#
# Entry point for the unit test suite. Run from the project root:
#   Rscript tests/testthat.R
#
# This project is a plain set of sourced R scripts, not a formal
# package, so there's no testthat::test_check() to hand off to.
# Instead, tests/testthat/helper-fixtures.R (auto-loaded by
# test_dir() before any test-*.R file runs) sources every R/*.R file
# and params.R itself, and defines the small fixture helpers shared
# across test files.

library(testthat)
test_dir("tests/testthat", reporter = "summary")
