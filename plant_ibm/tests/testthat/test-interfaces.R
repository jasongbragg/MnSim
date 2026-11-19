# tests/testthat/test-interfaces.R
#
# Swappable fire response, rust pressure and recruitment
# (docs/design_swappable_interfaces.md; notes T12): registries match the
# schema, optional fields keep existing species unchanged, and the default
# paths give exactly what the code gave before the split.

repo_root <- normalizePath(file.path(proj_root, ".."))
schema    <- read_config_yaml(file.path(repo_root, "config", "schema.yaml"))

# ---- Registries and schema ----------------------------------------------------
test_that("each registry matches its schema choices", {
  set.seed(1)
  expect_setequal(setdiff(names(FIRE_RESPONSES), "legacy"), unlist(schema$params$fire_extent_mode$choices))
  expect_setequal(names(RECRUITMENT_CURVES), unlist(schema$params$recruitment_curve$choices))
  expect_setequal(names(RUST_PRESSURE_MODES), unlist(schema$rust_anchor$mode$choices))
  expect_identical(schema$params$recruitment_curve$default, "bevholt")
  expect_identical(schema$rust_anchor$mode$default, "sdm")
})

# ---- Optional fields in the loader (section 2.1) --------------------------------
cfg_root <- function(edit) {
  root <- tempfile("mnsim_cfg_")
  dir.create(root)
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  f <- file.path(root, "species", "m_nodosa.yaml")
  writeLines(edit(readLines(f)), f)
  root
}
# Insert lines into the species file's params block, just after 'params:'.
add_params <- function(...) function(txt) {
  i <- match("params:", txt)
  append(txt, c(...), after = i)
}
curve_entry <- function(v)
  sprintf('  recruitment_curve: {value: "%s", prior: true, choices: [bevholt, ricker], reason: "test"}', v)

test_that("M. nodosa sets no optional field, so its list has none of them", {
  set.seed(1)
  p <- load_species("m_nodosa", root = repo_root)
  expect_false(any(c("recruitment_curve", "ricker_r", "ricker_K") %in% names(p)))
  expect_null(load_rust_anchor("m_nodosa", root = repo_root)$mode)
})

test_that("setting an optional field to its default gives the same list and run ID", {
  set.seed(1)
  root <- cfg_root(add_params(curve_entry("bevholt")))
  on.exit(unlink(root, recursive = TRUE))
  p0 <- load_species("m_nodosa", root = repo_root)
  p1 <- load_species("m_nodosa", root = root)
  expect_identical(p1, p0)
  expect_identical(make_run_id(p1), make_run_id(p0))
})

test_that("the rust anchor's mode set to its default leaves the anchor unchanged", {
  set.seed(1)
  root <- cfg_root(function(txt) append(txt, '  mode: {value: "sdm", prior: true, choices: [sdm], reason: "test"}',
                                       after = match("rust_anchor:", txt)))
  on.exit(unlink(root, recursive = TRUE))
  expect_identical(load_rust_anchor("m_nodosa", root = root), load_rust_anchor("m_nodosa", root = repo_root))
})

test_that("an unknown option is rejected when the configuration is loaded", {
  set.seed(1)
  root <- cfg_root(add_params(sprintf('  recruitment_curve: {value: "logistic", prior: true, choices: [logistic], reason: "test"}')))
  on.exit(unlink(root, recursive = TRUE))
  expect_error(load_species("m_nodosa", root = root), "'logistic' is not one of: bevholt, ricker", fixed = TRUE)
})

test_that("an option's own parameters are required only when it is selected", {
  set.seed(1)
  r_entry <- '  ricker_r: {value: 1.5, prior: true, min: 0, max: 3, reason: "test"}'
  k_entry <- '  ricker_K: {value: 50, prior: true, min: 1, max: 100, reason: "test"}'
  root1 <- cfg_root(add_params(curve_entry("ricker"), r_entry))
  root2 <- cfg_root(add_params(r_entry))
  root3 <- cfg_root(add_params(curve_entry("ricker"), r_entry, k_entry))
  on.exit(unlink(c(root1, root2, root3), recursive = TRUE))
  expect_error(load_species("m_nodosa", root = root1),
               "ricker_K: missing: needed when recruitment_curve is 'ricker'", fixed = TRUE)
  expect_error(load_species("m_nodosa", root = root2),
               "ricker_r: only used when recruitment_curve is 'ricker'; remove it", fixed = TRUE)
  p <- load_species("m_nodosa", root = root3)
  expect_identical(p[c("recruitment_curve", "ricker_r", "ricker_K")],
                   list(recruitment_curve = "ricker", ricker_r = 1.5, ricker_K = 50))
})

# ---- Fire responses (section 3) ---------------------------------------------------
# apply_fire() as it was before the split (commit 1c2e028), kept here to show
# each response makes the same draws and gives the same result.
apply_fire_before_split <- function(pop, params, t) {
  alive_idx <- which(pop$alive); n_alive <- length(alive_idx)
  if (n_alive == 0) return(pop)
  mode <- if (is.null(params$fire_extent_mode)) "plant" else params$fire_extent_mode
  if (params$fire_p_fimp > 0 && mode == "stand") {
    f <- draw_fire_size(params$fire_p_fimp, params$fire_size_sd)
    if (stats::runif(1) >= f) { attr(pop, "fire_reached") <- FALSE; return(pop) }
    cover <- if (is.null(params$fire_stand_cover)) 1 else params$fire_stand_cover
    imp_idx <- alive_idx[rbinom(n_alive, 1, cover) == 1]
    if (length(imp_idx) == 0) return(pop)
    ceiling <- params$fire_kill_scalar * f
  } else if (params$fire_p_fimp > 0) {
    imp_idx <- alive_idx[rbinom(n_alive, 1, params$fire_p_fimp) == 1]
    if (length(imp_idx) == 0) return(pop)
    ceiling <- params$fire_kill_scalar * params$fire_p_fimp
  }
  if (params$fire_p_fimp > 0) {
    p_kill <- ceiling + (1 - ceiling) *
      hill_weight(pop$age[imp_idx], params$fire_kill_half_sat, params$fire_kill_hill)
    killed <- rbinom(length(imp_idx), 1, p_kill) == 1
    dead_idx <- imp_idx[killed]; surv_idx <- imp_idx[!killed]
  } else {
    killed <- rbinom(n_alive, 1, params$fire_kill_prob) == 1
    dead_idx <- alive_idx[killed]; surv_idx <- alive_idx[!killed]
  }
  pop$alive[dead_idx] <- FALSE; pop$death_year[dead_idx] <- t
  if (!is.null(pop$death_cause)) pop$death_cause[dead_idx] <- "fire"
  n_surv <- length(surv_idx)
  if (n_surv > 0) {
    rust_delay <- rust_resprout_delay(pop, surv_idx, params, t)
    base  <- if (identical(params$resprout_delay_dist, "fixed"))
      rep(as.integer(params$resprout_yrs_base), n_surv) else rpois(n_surv, params$resprout_yrs_base)
    delay <- pmax(base + rpois(n_surv, rust_delay), 1L)
    pop$resprout[surv_idx] <- TRUE; pop$flowering[surv_idx] <- FALSE
    pop$resprout_yrs_remain[surv_idx] <- delay; pop$resprout_total_yrs[surv_idx] <- delay
  }
  pop
}

test_that("each fire response gives exactly the pre-split result", {
  set.seed(1)
  base <- make_toy_params(fire_kill_prob = 0.3, fire_p_fimp = 0.4, fire_size_sd = 0.2,
                          fire_stand_cover = 0.9, fire_kill_scalar = 0.3)
  pop  <- make_toy_pop(300, base)
  pop$age <- pop$age + 5L
  cases <- list(legacy = modifyList(base, list(fire_p_fimp = 0)),
                plant  = modifyList(base, list(fire_extent_mode = "plant")),
                stand  = modifyList(base, list(fire_extent_mode = "stand")),
                absent = base[setdiff(names(base), "fire_extent_mode")])
  for (nm in names(cases)) for (seed in 1:5) {
    set.seed(seed); new <- apply_fire(pop, cases[[nm]], 10L);              r1 <- runif(1)
    set.seed(seed); old <- apply_fire_before_split(pop, cases[[nm]], 10L); r2 <- runif(1)
    expect_identical(new, old, info = paste(nm, seed))
    expect_identical(r1, r2, info = paste(nm, seed, "(same number of draws)"))
  }
  expect_identical(fire_response_name(cases$legacy), "legacy")
  expect_identical(fire_response_name(cases$absent), "plant")
  expect_error(fire_response_name(modifyList(base, list(fire_extent_mode = "patchy"))),
               "Unknown fire_extent_mode 'patchy'")
})

# ---- Recruitment (section 5) --------------------------------------------------------
test_that("expected_recruits() applies the selected curve to the flowering plants", {
  set.seed(1)
  p   <- make_toy_params()
  pop <- make_toy_pop(200, p)
  pop$flowering <- pop$alive & pop$age >= pop$age_first_flower
  N   <- sum(pop$flowering)
  expect_gt(N, 0)
  shade <- 1 - shade_suppression(canopy_density(pop), p)
  expect_equal(expected_recruits(pop, 1L, p), p$R_max * N / (p$K_half + N) * shade)
  pr <- modifyList(p, list(recruitment_curve = "ricker", ricker_r = 1.5, ricker_K = 50))
  expect_equal(expected_recruits(pop, 1L, pr), N * exp(1.5 * (1 - N / 50)) * shade)
  expect_error(expected_recruits(pop, 1L, modifyList(p, list(recruitment_curve = "logistic"))),
               "Unknown recruitment_curve 'logistic'")
})

test_that("a short run with Ricker recruitment runs and recruits", {
  set.seed(1)
  p <- make_toy_params(recruitment_curve = "ricker", ricker_r = 1, ricker_K = 100, n_years = 15L)
  res <- run_simulation(p)
  expect_equal(nrow(res$census), 15L)
  expect_gt(sum(res$census$births), 0)
})

# ---- Rust pressure (section 4) -------------------------------------------------------
test_that("cell_rust_pressure() defaults to the SDM mode", {
  set.seed(1)
  sites  <- data.frame(sdm_suitability = c(0.2, 0.5, 0.9))
  anchor <- list(site = "x", lon = 0, lat = 0, mu = 0.7, conc = 20, cap = 0.95)
  expect_identical(cell_rust_pressure(sites, 2L, anchor), rust_pressure_from_sdm(sites$sdm_suitability, 2L, anchor))
  expect_identical(cell_rust_pressure(sites, 2L, c(anchor, mode = "sdm")), cell_rust_pressure(sites, 2L, anchor))
  expect_error(cell_rust_pressure(sites, 2L, c(anchor, mode = "flat")), "Unknown rust-pressure mode 'flat'")
})

test_that("draw_annual_env() is zero before rust and one Beta draw after", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 2010, annual_env_alpha = 14, annual_env_beta = 6)
  set.seed(3); r0 <- draw_annual_env(p, 2009L); a <- runif(1)
  set.seed(3); b <- runif(1)
  expect_identical(r0, 0)
  expect_identical(a, b)                         # no draw before rust
  set.seed(3); r1 <- draw_annual_env(p, 2010L)
  set.seed(3); expect_identical(r1, rbeta(1L, 14, 6))
})
