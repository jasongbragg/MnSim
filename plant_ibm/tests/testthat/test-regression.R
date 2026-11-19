# tests/testthat/test-regression.R
#
# Phase 1 regression harness (PLAN.md): re-runs the fixed-seed M. nodosa
# scenarios in tests/regression/scenarios.R and compares with the frozen
# outputs in ../tests/reference/ (repo root). See docs/methodological_notes.md T1.
#
# Two checks per scenario:
#   code        the frozen params must reproduce the frozen output exactly.
#               A failure means the model code changed its behaviour.
#   calibration today's params.R must build the frozen params. A failure
#               means params.R changed (e.g. a recalibration), not the code.
#
# If a check fails, the code is wrong, not the reference (CLAUDE.md rule 2):
# stop and report the difference. Never edit tests/reference/.

source(file.path(proj_root, "tests", "regression", "scenarios.R"), local = TRUE)

source(file.path(proj_root, "tests", "regression", "stats.R"), local = TRUE)

source(file.path(proj_root, "tests", "regression", "compare.R"), local = TRUE)   # first_difference()

# MNSIM_REFERENCE_DIR points the tests at another copy of the references
# (e.g. a dry run); default: tests/reference at the repo root.
ref_dir     <- Sys.getenv("MNSIM_REFERENCE_DIR",
                          normalizePath(file.path(proj_root, "..", "tests", "reference"), mustWork = FALSE))
params_file <- file.path(proj_root, "params.R")

test_that("every regression scenario has a reference file, made with this R version", {
  set.seed(1)
  files <- regression_reference_path(names(REGRESSION_SCENARIOS), ref_dir)
  expect_true(all(file.exists(files)),
              info = paste("missing:", paste(basename(files[!file.exists(files)]), collapse = ", "),
                           "-- generate with tests/regression/make_reference.R"))
  skip_if_not(all(file.exists(files)))
  for (f in files) {
    meta <- readRDS(f)$meta
    expect_identical(R.version.string, meta$R_version,
                     info = paste(basename(f), "was made with", meta$R_version,
                                  "-- exact matching needs the same R version (CLAUDE.md rule 3)"))
  }
})

for (nm in names(REGRESSION_SCENARIOS)) local({
  nm   <- nm
  file <- regression_reference_path(nm, ref_dir)

  test_that(sprintf("regression [%s]: frozen params reproduce the frozen output exactly", nm), {
    set.seed(1)   # convention; every scenario re-seeds from params$seed
    skip_if_not(file.exists(file), paste("no reference file", file))
    ref <- readRDS(file)
    out <- regression_run(nm, ref$params)
    d   <- first_difference(out, ref$output)
    expect(is.null(d), paste0("Model output differs from tests/reference/", nm, ".rds -- ",
                              d, ". The code changed its behaviour: stop and report (CLAUDE.md rule 2)."))
  })

  test_that(sprintf("regression [%s]: params.R still builds the frozen params", nm), {
    set.seed(1)
    skip_if_not(file.exists(file), paste("no reference file", file))
    ref <- readRDS(file)
    p   <- regression_params(nm, regression_base(params_file))
    d   <- first_difference(p, ref$params, path = "params")
    expect(is.null(d), paste0("params.R no longer gives the frozen params of '", nm, "' -- ", d,
                              ". params.R (the calibration) changed, not the code; references are",
                              " updated only by the project owner (PLAN.md Phase 1)."))
  })
})

# ---- Second check: summary statistics over 20 seeds (tests/regression/stats.R)
stats_file <- file.path(ref_dir, REGRESSION_STATS_FILE)

for (nm in names(REGRESSION_STATS)) local({
  nm <- nm
  test_that(sprintf("regression stats [%s]: 20-seed mean within tolerance of the reference", nm), {
    set.seed(1)
    skip_if_not(file.exists(stats_file), paste("no", stats_file,
                "-- generate with tests/regression/make_reference_stats.R"))
    ref <- readRDS(stats_file)$stats[[nm]]
    expect_false(is.null(ref), info = paste("no reference values for", nm))
    v   <- regression_stats_values(nm, ref$params)
    tol <- REGRESSION_STATS_Z * sqrt(2) * ref$sd / sqrt(length(ref$values))
    expect(abs(mean(v) - ref$mean) <= tol,
           sprintf(paste0("%s: 20-seed mean %.4f, reference %.4f (tolerance +/- %.4f). ",
                          "Model behaviour changed beyond sampling error: stop and report (CLAUDE.md rule 2)."),
                   nm, mean(v), ref$mean, tol))
  })
})
