# tests/testthat/test-csv-to-yaml.R
#
# scripts/csv_to_yaml.R (PLAN.md Phase 4; notes T22): the CSV-to-YAML contract
# reproduces M. nodosa's params exactly, and the gap logic follows the modes.

repo_root <- normalizePath(file.path(proj_root, ".."))
c2y <- new.env()
sys.source(file.path(repo_root, "scripts", "csv_to_yaml.R"), envir = c2y)

test_that("M. nodosa's params rebuild from its CSV to the same parameter list and priors", {
  set.seed(1)
  r <- c2y$build_species_yaml("m_nodosa", repo_root, blocks = c("rust_anchor", "params"))
  expect_true(r$ok)
  root <- tempfile("mnsim_c2y_"); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  old <- read_config_yaml(file.path(repo_root, "species", "m_nodosa.yaml"))
  new <- yaml::yaml.load(paste(r$text, collapse = "\n"))
  for (b in c("assessment", "spatial", "ensemble_draws")) new[[b]] <- old[[b]]
  yaml::write_yaml(new, file.path(root, "species", "m_nodosa.yaml"))
  a <- species_config("m_nodosa", root = repo_root)
  b <- species_config("m_nodosa", root = root)
  expect_identical(b$params, a$params)
  expect_identical(b$rust_anchor, a$rust_anchor)
  expect_identical(b$priors[, c("key", "min", "max", "choices")], a$priors[, c("key", "min", "max", "choices")])
})

test_that("field status follows the modes: required_when fields only under their mode", {
  set.seed(1)
  schema <- read_config_yaml(file.path(repo_root, "config", "schema.yaml"))
  row <- function(p, v) list(parameter = p, value = v, notes = "", evidence_type = "model_choice")
  st0 <- c2y$field_status(schema, list())
  expect_identical(st0$status[st0$key == "weibull_k"], "gap")
  expect_identical(st0$status[st0$key == "clonal_mode"], "default")
  expect_identical(st0$status[st0$key == "sucker_rate"], "inactive")
  expect_identical(st0$status[st0$key == "rust_anchor__mode"], "default")
  st1 <- c2y$field_status(schema, list(clonal_mode = row("clonal_mode", "genet"),
                                       rust_hazard_mode = row("rust_hazard_mode", "damage")))
  expect_identical(st1$status[st1$key == "sucker_rate"], "gap")
  expect_identical(st1$status[st1$key == "damage_recovery"], "gap")
  expect_identical(st1$status[st1$key == "clonal_mode"], "row")
})

test_that("a prior row becomes a prior entry with its range and reason", {
  set.seed(1)
  spec <- list(type = "double")
  e <- c2y$param_entry(list(value = "40", notes = "PRIOR [15, 60]: a reason", evidence_type = "model_assumption"),
                       spec, "rust_species_params.csv")
  expect_identical(e, list(value = 40, prior = TRUE, min = 15, max = 60, reason = "a reason",
                           cite = "rust_species_params.csv"))
  e <- c2y$param_entry(list(value = "damage", notes = "PRIOR {instant, damage}: why", evidence_type = "model_choice"),
                       list(type = "character"), "f.csv")
  expect_identical(e$choices, list("instant", "damage"))
})
