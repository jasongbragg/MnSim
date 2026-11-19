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

test_that("load_species('m_nodosa') is identical to params.R's get_default_params()", {
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
