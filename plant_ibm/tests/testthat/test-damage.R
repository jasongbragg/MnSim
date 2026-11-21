# tests/testthat/test-damage.R
#
# Accumulating rust damage (docs/design_accumulating_damage.md; notes T19).

repo_root <- normalizePath(file.path(proj_root, ".."))

dr <- function(max_effect, form, half_sat, hill) list(max_effect = max_effect, form = form,
                                                      half_sat = half_sat, hill = hill)
damage_params <- function(...) {
  base <- list(rust_hazard_mode = "damage",
               damage_gain_dose_response      = dr(0.5, "linear", 0, 1),
               damage_recovery                = 0.2,
               damage_mortality_dose_response = dr(0.6, "sigmoid", 0.8, 8))
  do.call(make_toy_params, modifyList(base, list(...)))
}
# A small population with given susceptibility, all mature.
toy_pop <- function(p, n = 50L) {
  pop <- create_population(modifyList(p, list(N0 = n)))$pop
  pop$age <- 50L; pop$resist_score <- 0; pop$env_susceptibility <- 1
  pop$rust_damage <- 0
  pop
}

test_that("instant mode leaves runs and outputs exactly as before", {
  set.seed(1)
  p <- make_toy_params(n_years = 30L, rust_start_year = 5, fire_prob_annual = 0.2, fire_p_fimp = 0.5)
  p$rust_pressure <- 1
  a <- run_simulation(p)
  b <- run_simulation(modifyList(p, list(rust_hazard_mode = "instant")))
  expect_identical(a, b)
  expect_false("rust_damage" %in% names(a$individuals))
  expect_false(any(c("mean_damage", "share_damaged") %in% names(a$census)))
  expect_null(load_species("m_nodosa", root = repo_root)$rust_hazard_mode)
})

test_that("the damage update makes no random draws", {
  set.seed(1)
  p <- damage_params(rust_start_year = 1)
  pop <- toy_pop(p)
  set.seed(7); u0 <- runif(1)
  set.seed(7); pop2 <- update_rust_damage(pop, 5, p, 0.7); u1 <- runif(1)
  expect_identical(u0, u1)
  expect_true(all(pop2$rust_damage > 0))
})

test_that("damage rises towards 1 under constant pressure, recovers without it, stays in [0, 1]", {
  set.seed(1)
  p <- damage_params(rust_start_year = 1)
  pop <- toy_pop(p, 5L)
  Ds <- numeric(0)
  for (t in 1:40) { pop <- update_rust_damage(pop, t, p, 0.9); Ds <- c(Ds, pop$rust_damage[1]) }
  expect_true(all(diff(Ds) >= -1e-12))
  # Equilibrium of D = D(1-r) + (1 - D(1-r)) g, with g = 0.45, r = 0.2
  g <- 0.45; r <- 0.2
  expect_equal(tail(Ds, 1), g / (1 - (1 - r) * (1 - g)), tolerance = 1e-6)
  D0 <- pop$rust_damage[1]
  pop <- update_rust_damage(pop, 41, p, 0)
  expect_equal(pop$rust_damage[1], D0 * (1 - r))
  # Before rust: no gain
  q <- modifyList(p, list(rust_start_year = 100))
  pop0 <- toy_pop(q, 5L)
  expect_identical(update_rust_damage(pop0, 10, q, 0.9)$rust_damage, rep(0, 5))
  # A huge gain is capped at full damage
  h <- modifyList(p, list(damage_gain_dose_response = dr(1, "linear", 0, 1), rust_pressure = 50))
  expect_identical(update_rust_damage(toy_pop(h, 3L), 5, h, 1)$rust_damage, rep(1, 3))
})

test_that("new plants start undamaged, and damage columns appear only in damage mode", {
  set.seed(1)
  p <- damage_params(n_years = 25L, rust_start_year = 5, annual_env_alpha = 14, annual_env_beta = 6)
  p$rust_pressure <- 1
  r <- run_simulation(p)
  ind <- r$individuals
  expect_true("rust_damage" %in% names(ind))
  expect_true(all(ind$rust_damage >= 0 & ind$rust_damage <= 1))
  born_last <- ind$birth_year == max(r$census$year)
  expect_true(all(ind$rust_damage[born_last] == 0))
  expect_true(all(c("mean_damage", "share_damaged") %in% names(r$census)))
  expect_true(all(r$census$mean_damage[r$census$year < 5] == 0))
  expect_gt(max(r$census$mean_damage, na.rm = TRUE), 0)
  # Clonal suckers too
  cp <- modifyList(p, list(clonal_mode = "genet", clonal_founder_genets = 20L, sucker_rate = 0.1,
                           sucker_rate_rootstock = 0.5, sucker_pulse = 1,
                           rootstock_fire_survival = 0.8, rootstock_hazard_no_canopy = 0))
  rc <- run_simulation(cp)
  s  <- rc$individuals[rc$individuals$origin == "sucker", ]
  expect_gt(nrow(s), 0)
  expect_true(all(s$rust_damage[s$birth_year == max(rc$census$year)] == 0))
})

test_that("fire top-kill resets damage when the reset is on, keeps it when off", {
  set.seed(1)
  p <- damage_params(fire_p_fimp = 1, fire_extent_mode = "plant", fire_size_sd = 0,
                     fire_kill_scalar = 0, fire_kill_prob = 0)
  pop <- toy_pop(p, 40L); pop$rust_damage <- 0.7
  a <- apply_fire(pop, p, 10)
  rs <- which(a$alive & a$resprout)
  expect_gt(length(rs), 0)
  expect_true(all(a$rust_damage[rs] == 0))
  set.seed(1)
  b <- apply_fire(pop, modifyList(p, list(damage_fire_reset = "keep")), 10)
  expect_true(all(b$rust_damage[b$alive & b$resprout] == 0.7))
})

test_that("with constant pressure, survivors' death rate rises in damage mode but not in instant mode", {
  set.seed(1)
  # Background and age-driven rust risk switched off, so rust is the only cause.
  base <- make_toy_params(n_years = 12L, N0 = 500L, R_max = 0, rust_start_year = 1,
                          annual_env_alpha = 1e4, annual_env_beta = 1e4 * 0.3 / 0.7,
                          microclim_alpha = 1e4, microclim_beta = 1, weibull_lambda = 1e9,
                          age_first_flower_mean = 1, age_first_flower_sd = 0)
  base$rust_pressure <- 1; base$resist_freq0[] <- 0
  base$rust_dose_response$age_peak  <- 0.15
  base$rust_dose_response$age_floor <- 0.15
  rate <- function(cen) { n <- cen$N_alive; (head(n, -1) - tail(n, -1)) / head(n, -1) }
  inst <- rate(run_simulation(base)$census)
  dmg  <- rate(run_simulation(modifyList(base, list(
    rust_hazard_mode = "damage",
    damage_gain_dose_response      = dr(0.4, "linear", 0, 1),
    damage_recovery                = 0,
    damage_mortality_dose_response = dr(0.6, "sigmoid", 0.7, 6))))$census)
  expect_lt(mean(inst[7:11]) - mean(inst[1:3]), 0.03)        # flat (selection can only lower it)
  expect_gt(mean(dmg[7:11]), mean(dmg[1:3]) + 0.15)          # accelerates
})

test_that("damage settings validate through the schema", {
  set.seed(1)
  schema <- read_config_yaml(file.path(repo_root, "config", "schema.yaml"))
  expect_setequal(RUST_HAZARD_MODES, unlist(schema$params$rust_hazard_mode$choices))
  expect_identical(schema$params$rust_hazard_mode$default, "instant")
  expect_identical(schema$params$damage_recovery$required_when$rust_hazard_mode, "damage")
  for (k in c("max_effect", "form", "half_sat", "hill")) {
    expect_identical(schema$params$damage_gain_dose_response[[k]]$required_when$rust_hazard_mode, "damage")
    expect_identical(schema$params$damage_mortality_dose_response[[k]]$required_when$rust_hazard_mode, "damage")
  }
  # Selecting damage without its settings is refused
  root <- tempfile("mnsim_cfg_"); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  f <- file.path(root, "species", "m_nodosa.yaml"); txt <- readLines(f)
  writeLines(append(txt, '  rust_hazard_mode: {value: "damage", prior: true, choices: [instant, damage], reason: "test"}',
                    after = match("params:", txt)), f)
  expect_error(load_species("m_nodosa", root = root, check_sources = FALSE), "needed when rust_hazard_mode is 'damage'")
})
