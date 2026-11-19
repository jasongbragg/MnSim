# tests/testthat/test-species-config.R
#
# The species configuration loader (R/species_config.R; notes T7):
# species/m_nodosa.yaml + config/ must rebuild params.R's get_default_params()
# exactly, and the validation must reject broken files.

source(file.path(proj_root, "R", "species_config.R"), local = TRUE)

repo_root <- normalizePath(file.path(proj_root, ".."))

working_params_env <- function() {
  env <- new.env()
  sys.source(file.path(proj_root, "params.R"), envir = env)
  env
}

# params.R now only wraps the loader (notes T8). Its frozen pre-switch output is
# guarded by test-regression.R ("params.R still builds the frozen params"); these
# two tests check the wiring.
test_that("params.R's get_default_params() is load_species('m_nodosa')", {
  set.seed(1)
  env <- working_params_env()
  p   <- load_species("m_nodosa", root = repo_root)
  expect_identical(p, env$get_default_params())
  # Same hash, so cached results stay valid (PLAN.md Phase 2 note on YAML types).
  expect_identical(make_run_id(p), make_run_id(env$get_default_params()))
})

test_that("load_rust_anchor('m_nodosa') is identical to RUST_ANCHOR", {
  set.seed(1)
  env <- working_params_env()
  expect_identical(load_rust_anchor("m_nodosa", root = repo_root), env$RUST_ANCHOR)
})

test_that("every M. nodosa CSV parameter row is used by the species or landscape file", {
  set.seed(1)
  csv  <- utils::read.csv(file.path(repo_root, "data", "m_nodosa_params.csv"), colClasses = "character")
  rows <- csv$parameter[csv$role == "parameter"]
  sp   <- read_config_yaml(file.path(repo_root, "species", "m_nodosa.yaml"))
  land <- read_config_yaml(file.path(repo_root, "config", "landscape.yaml"))
  key  <- function(x, pre) vapply(config_paths(x, is_entry), function(p) paste(c(pre, p), collapse = "__"), "")
  used <- c(key(sp$params, character()), key(sp$rust_anchor, "rust_anchor"), key(land$params, character()))
  expect_identical(sort(setdiff(rows, used)), character(0))
  expect_identical(sort(setdiff(used, rows)), character(0))
})

test_that("species_priors() lists the priors with their ranges", {
  set.seed(1)
  pr <- species_priors("m_nodosa", root = repo_root)
  expect_true(all(c("key", "layer", "min", "max", "choices", "reason") %in% names(pr)))
  sen <- pr[pr$key == "senescence_dose_response__max_effect", ]
  expect_equal(c(sen$min, sen$max), c(0, 0.6))
  expect_true("fire_prob_annual" %in% pr$key[pr$layer == "landscape"])
})

# ---- Validation: copy the config to a temporary root, break one thing ----------
broken_root <- function(file, from, to) {
  root <- tempfile("mnsim_cfg_")
  dir.create(root)
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  path <- file.path(root, file)
  txt  <- readLines(path)
  new  <- sub(from, to, txt, fixed = TRUE)
  if (identical(new, txt)) stop("test setup: '", from, "' not found in ", file)
  writeLines(new, path)
  root
}
expect_config_error <- function(file, from, to, pattern) {
  root <- broken_root(file, from, to)
  on.exit(unlink(root, recursive = TRUE))
  expect_error(load_species("m_nodosa", root = root), pattern, fixed = TRUE)
}

test_that("validation rejects a missing field", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "  weibull_k: {value: 1, cite: m_nodosa_params.csv}", "",
                      "weibull_k: missing")
})

test_that("validation rejects a field the schema does not know", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "  weibull_k: {value: 1, cite: m_nodosa_params.csv}",
                      "  weibull_k: {value: 1, cite: m_nodosa_params.csv}\n  weibul_lambda: {value: 60, cite: m_nodosa_params.csv}",
                      "weibul_lambda: not a field of this layer")
})

test_that("validation rejects a value of the wrong type", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "weibull_lambda: {value: 60,", "weibull_lambda: {value: sixty,",
                      "expected a number, got 'sixty'")
})

test_that("validation rejects a value outside the schema bounds", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "selfing_rate: {value: 0.10,", "selfing_rate: {value: 1.5,",
                      "1.5 is above the maximum 1")
})

test_that("validation rejects a value not among the allowed choices", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "fire_extent_mode: {value: \"stand\",",
                      "fire_extent_mode: {value: \"everywhere\",", "'everywhere' is not one of: stand, plant")
})

test_that("validation rejects a species value with neither citation nor prior (rule 5)", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "weibull_lambda: {value: 60, cite: m_nodosa_params.csv}",
                      "weibull_lambda: {value: 60}", "no citation")
})

test_that("validation rejects a value that differs from its cited CSV row", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "weibull_lambda: {value: 60,", "weibull_lambda: {value: 61,",
                      "value 61 differs from data/m_nodosa_params.csv (60)")
})

test_that("validation rejects a prior value outside its own range", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "      max: 0.6", "      max: 0.05",
                      "senescence_dose_response$max_effect: value outside its prior [0, 0.05]")
})

test_that("validation rejects a landscape value with neither source nor prior", {
  set.seed(1)
  expect_config_error("config/landscape.yaml", "    source: \"docs/methodological_notes.md D3",
                      "    note: \"docs/methodological_notes.md D3", "needs a 'source' or must be a prior")
})

# ---- Settings outside the params list (notes T9) --------------------------------
# The values the runners had hard-coded before the move, which the config must
# reproduce exactly (they are hashed into run IDs and stored in results).
test_that("load_assessment() gives the old iucn_assess_spec() defaults", {
  set.seed(1)
  old <- list(criterion = "A4", year_spinup = 1600L, year_gentime_start = 1810L,
              year_gentime_end = 1910L, year_rust = 2010L, year_assess = 2026L,
              year_end = 2210L, avg_years = 11L, max_future = 100L)
  expect_identical(load_assessment("m_nodosa", root = repo_root), old)
  expect_identical(iucn_assess_spec(), old)   # runner/iucn_a4.R, sourced by the helper
})

test_that("load_categories() gives the IUCN thresholds 30/50/80", {
  set.seed(1)
  expect_identical(load_categories("m_nodosa", root = repo_root), c(VU = 30, EN = 50, CR = 80))
  expect_identical(vapply(c(29.9, 30, 50, 80), iucn_category, ""), c("LC/NT", "VU", "EN", "CR"))
})

test_that("load_param_draws() gives the old DEFAULT_PARAM_DRAWS", {
  set.seed(1)
  old <- list(resprout_rust_mortality = list(
    path  = c("rust_dose_response", "resprout", "max_effect"),
    min   = 0.32, max = 0.78,
    label = "Extra rust mortality of resprouting plants (max_effect)",
    note  = "values whose Pegg-scenario median lies in the 95% CI of 24/50 dead (0.34-0.63), at the Yarringully cell"))
  expect_identical(load_param_draws("m_nodosa", root = repo_root), old)
  expect_identical(DEFAULT_PARAM_DRAWS, old)
})

test_that("validation_config() gives the old Pegg settings", {
  set.seed(1)
  pegg <- validation_config("m_nodosa", root = repo_root)$pegg
  expect_identical(pegg[c("year_spinup", "year_rust", "year_fire", "year_survey")],
                   list(year_spinup = 1910L, year_rust = 2010L, year_fire = 2020L, year_survey = 2022L))
  expect_identical(pegg$observed, 0.48)
  expect_identical(pegg$observed_n, 50L)
  expect_identical(pegg$site, "anchor")
  expect_identical(pegg$n_post, 3L)
  expect_identical(pegg$plots, "runner/data/pegg2025_mnodosa_plots.csv")
  expect_identical(pegg$scenario, list(fire_p_fimp = 1, fire_size_sd = 0, fire_stand_cover = 1))
  expect_true(file.exists(file.path(proj_root, pegg$plots)))
  expect_null(validation_config("no_such_species", root = repo_root))
})

test_that("spatial settings: species range and shared landscape values", {
  set.seed(1)
  sp <- load_species_spatial("m_nodosa", root = repo_root)
  expect_identical(sp$name, "Melaleuca nodosa")
  expect_identical(unlist(sp$range_box), c(lon_min = 148, lon_max = 155, lat_min = -38, lat_max = -24))
  expect_identical(sp$data_dir, "spatial/data/m_nodosa")
  expect_true(file.exists(file.path(proj_root, sp$data_dir, "site_params.csv")))
  land <- load_landscape_spatial(root = repo_root)
  expect_identical(land$grid_deg, 0.1)
  expect_identical(land$qld_border_lat, -28.16)
  expect_identical(default_species(repo_root), "m_nodosa")
})

test_that("validation rejects a setting without a source", {
  set.seed(1)
  expect_config_error("config/assessment.yaml", '    source: "docs/methodological_notes.md A2', '    note: "docs/methodological_notes.md A2',
                      "config/assessment.yaml: assessment$avg_years: needs a 'source'")
})

test_that("validation rejects assessment years out of order", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "    value: 1810", "    value: 1500",
                      "years must run spin-up < gen-time start")
})

test_that("validation rejects an unknown setting", {
  set.seed(1)
  expect_config_error("config/assessment.yaml", "  avg_years:", "  avg_yrs:\n    value: 11\n    source: x\n  avg_years:",
                      "assessment$avg_yrs: not a setting of this file")
})

test_that("validation rejects a parameter draw whose range differs from the CSV", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "    min: 0.32", "    min: 0.30",
                      "range differs from 'RANGE [min, max]'")
})

test_that("validation rejects a parameter draw outside the parameter's own value", {
  set.seed(1)
  expect_config_error("species/m_nodosa.yaml", "    max: 0.78", "    max: 0.5",
                      "the parameter's value 0.53 is outside [0.32, 0.5]")
})

test_that("validation of the validation file: species name and cited row", {
  set.seed(1)
  root <- broken_root("species/m_nodosa_validation.yaml", 'species: "Melaleuca nodosa"', 'species: "Melaleuca nodos"')
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_error(validation_config("m_nodosa", root = root), "'species' must be \"Melaleuca nodosa\"", fixed = TRUE)
  root2 <- broken_root("species/m_nodosa_validation.yaml", "row: pegg_yarringully_resprout_dead", "row: pegg_no_such_row")
  on.exit(unlink(root2, recursive = TRUE), add = TRUE)
  expect_error(validation_config("m_nodosa", root = root2), "0 rows in data/m_nodosa_params.csv", fixed = TRUE)
})
