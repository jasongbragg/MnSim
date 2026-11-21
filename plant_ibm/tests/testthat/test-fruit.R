# tests/testthat/test-fruit.R
#
# Rust suppression of flowering and fruit (docs/design_fruit_infection.md;
# notes T21): site pressure in both hooks, mothers weighted by fruit success,
# flowering driven by crown damage, and the census columns.

repo_root <- normalizePath(file.path(proj_root, ".."))

dr <- function(max_effect, form, half_sat, hill) list(max_effect = max_effect, form = form,
                                                      half_sat = half_sat, hill = hill)

flowering_pop <- function(n, p) {
  pop <- make_toy_pop(n, p)
  pop$alive <- TRUE; pop$resprout <- FALSE
  pop$age <- pmax(pop$age, pop$age_first_flower)
  pop$flowering <- TRUE
  pop
}

test_that("site rust pressure scales fruit success; zero pressure means no loss", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 1, rust_recruit_dose_response = dr(0.8, "linear", 0, 1))
  pop <- flowering_pop(50, p); pop$resist_score <- 0; pop$env_susceptibility <- 1
  p$rust_pressure <- 1
  expect_equal(fruit_success(pop, 1, p, 0.5), rep(1 - 0.8 * 0.5, 50))
  p$rust_pressure <- 0.5
  expect_equal(fruit_success(pop, 1, p, 0.5), rep(1 - 0.8 * 0.25, 50))
  p$rust_pressure <- 0
  expect_equal(fruit_success(pop, 1, p, 0.5), rep(1, 50))
  expect_null(fruit_success(pop, 0, p, 0.5))                       # before rust
  expect_null(fruit_success(pop, 1, modifyList(p, list(rust_recruit_dose_response = dr(0, "linear", 0, 1))), 0.5))
})

test_that("site rust pressure scales flowering suppression", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 1, rust_flower_dose_response = dr(1, "linear", 0, 1))
  pop <- flowering_pop(400, p); pop$resist_score <- 0; pop$env_susceptibility <- 1
  share <- function(rp) { p$rust_pressure <- rp; mean(!update_flowering(pop, 1, p, 0.8)$flowering) }
  expect_identical(share(0), 0)
  expect_true(share(1) > 0.7 && share(1) < 0.9)
  expect_true(share(0.5) > 0.3 && share(0.5) < 0.5)
})

test_that("mothers follow fruit success; fathers do not; no weights draws as before", {
  set.seed(1)
  p <- make_toy_params(selfing_rate = 0)
  pop <- flowering_pop(20, p)
  w <- c(0, rep(1, 19))
  par <- sample_parents(pop, 4000L, p, mother_w = w)
  expect_false(any(par$mother_idx == 1L))
  expect_true(any(par$father_idx == 1L))
  set.seed(7); a <- sample_parents(pop, 50L, p)
  set.seed(7); b <- sample_parents(pop, 50L, p, mother_w = NULL)
  expect_identical(a, b)
  expect_null(sample_parents(pop, 5L, p, mother_w = rep(0, 20)))
  expect_error(sample_parents(pop, 5L, p, mother_w = 1), "one weight per flowering plant")
})

test_that("recruitment selects for resistance through fruit success", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 1, R_max = 2000, K_half = 10,
                       rust_recruit_dose_response = dr(1, "linear", 0, 1))
  p$rust_pressure <- 1
  pop <- flowering_pop(200, p); pop$env_susceptibility <- 1
  pop$resist_score <- rep(c(0, 0.9), each = 100)
  gt <- make_toy_gt(200, p)
  r <- recruit(pop, gt, 1, p, annual_env_t = 1)
  new <- r$pop[(nrow(pop) + 1):nrow(r$pop), ]
  m <- match(new$mother_id, pop$id)
  expect_gt(mean(m > 100), 0.85)                        # resistant mothers (weight 0.9 vs 0)
  expect_equal(r$fruit_success, mean(c(rep(0, 100), rep(0.9, 100))))
})

test_that("damage-driven flowering follows crown damage and persists in light years", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 1, rust_flower_dose_response = dr(1, "threshold", 0.5, 1),
                       rust_flower_driver = "damage", rust_hazard_mode = "damage")
  pop <- flowering_pop(100, p)
  pop$rust_damage <- rep(c(0.2, 0.9), each = 50)
  f <- update_flowering(pop, 2, p, annual_env_t = 0)$flowering   # no pressure this year
  expect_true(all(f[1:50])); expect_false(any(f[51:100]))
  expect_true(all(update_flowering(pop, 0, p, 0.9)$flowering))  # before rust
  pop$rust_damage <- NULL
  expect_error(update_flowering(pop, 2, p, 0.5), "needs rust_hazard_mode: damage")
  expect_error(rust_flower_driver(list(rust_flower_driver = "fruit")), "Unknown rust_flower_driver")
})

test_that("census columns appear only with a reproductive hook on; runs otherwise unchanged", {
  set.seed(1)
  p <- make_toy_params(n_years = 12L, rust_start_year = 3, annual_env_alpha = 14, annual_env_beta = 6)
  p$rust_pressure <- 1
  r0 <- run_simulation(p)
  expect_false(any(c("share_flowering_suppressed", "fruit_success") %in% names(r0$census)))
  q <- modifyList(p, list(rust_flower_dose_response = dr(0.9, "linear", 0, 1),
                          rust_recruit_dose_response = dr(0.9, "linear", 0, 1)))
  r1 <- run_simulation(q)
  cs <- r1$census
  expect_true(all(c("share_flowering_suppressed", "fruit_success") %in% names(cs)))
  expect_true(all(cs$share_flowering_suppressed[cs$year < 3] == 0))
  expect_gt(mean(cs$share_flowering_suppressed[cs$year >= 3]), 0.1)
  expect_true(all(cs$fruit_success[cs$year >= 3] < 1, na.rm = TRUE))
  expect_identical(cs$N_IUCN[cs$year < 3], r0$census$N_IUCN[r0$census$year < 3])
})

test_that("rust_flower_driver validates through the schema", {
  set.seed(1)
  schema <- read_config_yaml(file.path(repo_root, "config", "schema.yaml"))
  expect_setequal(RUST_FLOWER_DRIVERS, unlist(schema$params$rust_flower_driver$choices))
  expect_null(load_species("m_nodosa", root = repo_root)$rust_flower_driver)
  root <- tempfile("mnsim_cfg_"); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  f <- file.path(root, "species", "m_nodosa.yaml"); txt <- readLines(f)
  writeLines(append(txt, '  rust_flower_driver: {value: "damage", prior: true, choices: [pressure, damage], reason: "test"}',
                    after = match("params:", txt)), f)
  expect_error(load_species("m_nodosa", root = root, check_sources = FALSE), "needs rust_hazard_mode: damage")
})
