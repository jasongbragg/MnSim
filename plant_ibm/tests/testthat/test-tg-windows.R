# tests/testthat/test-tg-windows.R
#
# Warnings when a run's spin-up or generation-time window is short relative to
# that run's own T_g (runner/iucn_a4.R, check_tg_windows(); notes T11).

repo_root <- normalizePath(file.path(proj_root, ".."))
assess    <- iucn_assess_spec()   # M. nodosa: spin-up 1600, gen-time window 1810-1910

test_that("the thresholds come from config/assessment.yaml", {
  set.seed(1)
  expect_identical(load_window_checks("m_nodosa", root = repo_root),
                   list(gentime_window_min_tg = 3, spinup_min_tg = 5))
  expect_identical(IUCN_WINDOW_CHECKS, list(gentime_window_min_tg = 3, spinup_min_tg = 5))
})

test_that("M. nodosa's windows pass at its T_g (about 23 yr)", {
  set.seed(1)
  expect_identical(check_tg_windows(23, assess), character(0))
  expect_identical(check_tg_windows(NA_real_, assess), character(0))
})

test_that("each short window gives one message", {
  set.seed(1)
  # T_g 40: gen-time window 100 yr = 2.5 T_g (< 3); spin-up 210 yr = 5.25 T_g (ok).
  m <- check_tg_windows(40, assess)
  expect_length(m, 1L)
  expect_match(m, "gen-time window (1810-1910) is 100 yr = 2.5 T_g, short of 3 T_g", fixed = TRUE)
  # T_g 50: both short (spin-up 210 yr = 4.2 T_g).
  m <- check_tg_windows(50, assess)
  expect_length(m, 2L)
  expect_match(m[2], "spin-up before it (1600-1810) is 210 yr = 4.2 T_g, short of 5 T_g", fixed = TRUE)
  # Exactly at the threshold is not short.
  edge <- assess
  edge$year_spinup <- 1710L                                      # spin-up 100 yr = 5 * 20
  edge$year_gentime_end <- 1870L                                 # window   60 yr = 3 * 20
  expect_identical(check_tg_windows(20, edge), character(0))
  expect_length(check_tg_windows(20.01, edge), 2L)
})

test_that("the warnings have their own class, so callers can muffle just these", {
  set.seed(1)
  expect_warning(warn_tg_windows(check_tg_windows(50, assess)), class = "mnsim_tg_window_warning")
  expect_silent(withCallingHandlers(warn_tg_windows(check_tg_windows(50, assess)),
                                    mnsim_tg_window_warning = function(w) invokeRestart("muffleWarning")))
})

test_that("summarise_iucn_a4() flags results with a short window", {
  set.seed(1)
  # The shortened regression A4 run: spin-up 1700-1800 = 100 yr, T_g 22.7 -> 4.4 T_g.
  r <- readRDS(file.path(repo_root, "tests", "reference", "a4.rds"))$output
  expect_false(summarise_iucn_a4(r)$tg_windows_ok)
  r$year_spinup <- 1600L   # 200 yr = 8.8 T_g
  expect_true(summarise_iucn_a4(r)$tg_windows_ok)
  expect_message(report_tg_windows(c(TRUE, FALSE, NA)), "1 of 3 runs had a window short")
})

test_that("validation rejects a window check without a source", {
  set.seed(1)
  root <- tempfile("mnsim_cfg_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  f   <- file.path(root, "config", "assessment.yaml")
  txt <- readLines(f)
  i   <- grep("source: \"docs/methodological_notes.md T11: spin-up", txt, fixed = TRUE)
  expect_length(i, 1L)
  txt[i] <- sub("source:", "note:", txt[i], fixed = TRUE)
  writeLines(txt, f)
  expect_error(load_window_checks("m_nodosa", root = root),
               "config/assessment.yaml: window_checks$spinup_min_tg: needs a 'source'", fixed = TRUE)
})
