# tests/testthat/test-rust-modes.R
#
# Rust-pressure modes (docs/design_rust_pressure_modes.md; notes T18): the
# anchor as a property of place, the constant and sdm_gradient modes, and rust
# arrival year and ramp per species.

repo_root <- normalizePath(file.path(proj_root, ".."))
.A <- load_rust_anchor("m_nodosa", root = repo_root)

cfg_copy <- function() {
  root <- tempfile("mnsim_cfg_")
  dir.create(root)
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  root
}
edit_file <- function(root, file, edit) {
  f <- file.path(root, file)
  writeLines(edit(readLines(f)), f)
}
add_species_lines <- function(...) function(txt) append(txt, c(...), after = match("params:", txt) - 1L)
add_params <- function(...) function(txt) append(txt, c(...), after = match("params:", txt))

# ---- The anchor: a property of place --------------------------------------------
test_that("the landscape anchor carries the M. nodosa CSV values", {
  set.seed(1)
  csv  <- utils::read.csv(file.path(repo_root, "data", "m_nodosa_params.csv"), colClasses = "character")
  land <- read_config_yaml(file.path(repo_root, "config", "landscape.yaml"))$rust_anchor
  for (k in names(land)) {
    row <- csv[csv$parameter == paste0("rust_anchor__", k), ]
    expect_identical(nrow(row), 1L)
    expect_true(config_values_match(row$value, land[[k]]$value), info = k)
  }
  expect_identical(names(.A), c("site", "lon", "lat", "mu", "conc", "cap", "cell_lon", "cell_lat"))
  # The stored M. nodosa table's anchor cell is the landscape anchor cell
  s <- load_site_params(drop_incomplete = FALSE)
  a <- s[s$rust_anchor_cell %in% TRUE, ]
  expect_equal(c(a$lon_centre, a$lat_centre), c(.A$cell_lon, .A$cell_lat))
})

test_that("a species file may choose the mode but not the anchor level", {
  set.seed(1)
  root <- cfg_copy(); on.exit(unlink(root, recursive = TRUE))
  edit_file(root, "species/m_nodosa.yaml",
            add_species_lines("rust_anchor:", '  mu: {value: 0.5, cite: m_nodosa_params.csv}'))
  expect_error(load_rust_anchor("m_nodosa", root = root), "property of place")
  root2 <- cfg_copy(); on.exit(unlink(root2, recursive = TRUE), add = TRUE)
  edit_file(root2, "species/m_nodosa.yaml", add_species_lines(
    "rust_anchor:", '  mode: {value: "constant", prior: true, choices: [sdm, constant, sdm_gradient], reason: "test"}'))
  a <- load_rust_anchor("m_nodosa", root = root2)
  expect_identical(a$mode, "constant")
  expect_identical(a[names(.A)], .A)
  # The species' params do not change with the mode
  expect_identical(load_species("m_nodosa", root = root2), load_species("m_nodosa", root = repo_root))
})

test_that("a species that does not occupy the anchor cell normalises to the same cell", {
  skip_if_not_installed("terra")
  set.seed(1)
  r <- terra::rast(xmin = 152.5, xmax = 153.5, ymin = -29.5, ymax = -28.5, res = 0.025, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  terra::values(r) <- 0.5 + 0.4 * (xy[, 2] + 29.5)
  s <- data.frame(cell_id = 1:3, lon_centre = c(153.25, 152.95, 153.05), lat_centre = c(-29.15, -28.75, -29.35))
  with_anchor <- suppressMessages(apply_sdm_to_sites(s, r, .A))
  without     <- suppressMessages(apply_sdm_to_sites(s[2:3, ], r, .A))
  expect_equal(without$annual_env_mean, with_anchor$annual_env_mean[2:3])
  expect_false(any(without$rust_anchor_cell))
  expect_error(suppressMessages(apply_sdm_to_sites(without, NULL, .A)), "does not occupy the rust anchor cell")
})

# ---- Modes ------------------------------------------------------------------------
test_that("constant mode gives mu in every cell, with or without an SDM value", {
  set.seed(1)
  sites <- data.frame(lon_centre = c(153.25, 150, 151), lat_centre = c(-29.15, -33, -31),
                      sdm_suitability = c(0.9, NA, 0.2))
  out <- suppressMessages(apply_sdm_to_sites(sites, NULL, c(.A, mode = "constant")))
  expect_equal(out$annual_env_mean, rep(.A$mu, 3))
  expect_equal(out$annual_env_alpha + out$annual_env_beta, rep(.A$conc, 3))
  expect_equal(out$sdm_rel[3], 0.2 / 0.9)
  # No SDM column at all, anchor not occupied: still mu
  rp <- cell_rust_pressure(data.frame(x = 1:2), NA_real_, c(.A, mode = "constant"))
  expect_equal(rp$annual_env_mean, rep(.A$mu, 2))
})

test_that("sdm_gradient's site table equals sdm's; the gradient is added per run", {
  set.seed(1)
  sites <- data.frame(sdm_suitability = c(0.2, 0.5, 0.9))
  expect_identical(cell_rust_pressure(sites, 0.5, c(.A, mode = "sdm_gradient")),
                   cell_rust_pressure(sites, 0.5, .A))
})

test_that("apply_rust_gradient: g = 0 is a no-op, the anchor latitude keeps its mean, the sign is as stated", {
  set.seed(1)
  p <- make_toy_params(); p$annual_env_alpha <- 14; p$annual_env_beta <- 6
  expect_identical(apply_rust_gradient(p, 0, -35, .A$cell_lat), p)
  expect_identical(apply_rust_gradient(p, NULL, -35, .A$cell_lat), p)
  for (g in c(-1, -0.3, 0.8)) {
    q <- apply_rust_gradient(p, g, .A$cell_lat, .A$cell_lat)
    expect_equal(q$annual_env_alpha, 14); expect_equal(q$annual_env_beta, 6)
    expect_identical(q$rust_gradient, g)
  }
  # 10 degrees south with g = -1: logit pressure down by exactly 1, concentration kept
  q <- apply_rust_gradient(p, -1, .A$cell_lat - 10, .A$cell_lat)
  expect_equal(stats::qlogis(q$annual_env_alpha / 20), stats::qlogis(0.7) - 1)
  expect_equal(q$annual_env_alpha + q$annual_env_beta, 20)
  # ... and up by 1 ten degrees north
  q <- apply_rust_gradient(p, -1, .A$cell_lat + 10, .A$cell_lat)
  expect_equal(stats::qlogis(q$annual_env_alpha / 20), stats::qlogis(0.7) + 1)
})

test_that("range plans: the gradient is hashed into run IDs and plans pair up", {
  set.seed(1)
  sites <- load_range_sites(SITE_PARAMS_PATH)
  base  <- get_default_params()
  ga    <- c(.A, mode = "sdm_gradient")
  draws <- list(grad = list(path = "rust_gradient", min = -1, max = 0, label = "g", note = "test"))
  a0 <- plan_range_run(5L, sites, base, param_draws = draws, anchor = ga, rust_gradient = 0)
  ad <- plan_range_run(5L, sites, base, param_draws = draws, anchor = ga)
  af <- plan_range_run(5L, sites, base, param_draws = draws, anchor = ga, rust_gradient = -0.5)
  keep <- c("cell_row", "cell_id", "sim_seed", "site_fire_prob", "grad")
  expect_identical(a0$row[keep], ad$row[keep]); expect_identical(a0$row[keep], af$row[keep])
  expect_identical(ad$row$rust_gradient, ad$row$grad)
  expect_identical(af$row$rust_gradient, -0.5)
  expect_identical(af$params$rust_gradient, -0.5)
  expect_null(a0$params$rust_gradient)
  expect_length(unique(c(a0$row$run_id, ad$row$run_id, af$row$run_id)), 3L)
  # g = 0 is the sdm run exactly
  s0 <- plan_range_run(5L, sites, base, param_draws = list(), anchor = .A)
  z0 <- plan_range_run(5L, sites, base, param_draws = list(), anchor = ga, rust_gradient = 0)
  expect_identical(s0$row$run_id, z0$row$run_id)
  expect_identical(s0$params, z0$params)
  expect_false("rust_gradient" %in% names(s0$row))
  # A gradient needs the sdm_gradient mode
  expect_error(plan_range_run(5L, sites, base, param_draws = list(), anchor = .A, rust_gradient = -0.5),
               "needs the species' rust_anchor mode 'sdm_gradient'")
})

test_that("ensemble_draws may draw rust_gradient, a drawable run setting, but not other run settings", {
  set.seed(1)
  root <- cfg_copy(); on.exit(unlink(root, recursive = TRUE))
  draw <- function(param) function(txt) c(txt, sprintf("  g_test:\n    param: %s\n    min: -1\n    max: 0\n    label: \"g\"\n    note: \"test\"\n    cite: m_nodosa_params.csv", param))
  edit_file(root, "species/m_nodosa.yaml", draw("rust_gradient"))
  d <- load_param_draws("m_nodosa", root = root, check_sources = FALSE)
  expect_identical(d$g_test$path, "rust_gradient")
  root2 <- cfg_copy(); on.exit(unlink(root2, recursive = TRUE), add = TRUE)
  edit_file(root2, "species/m_nodosa.yaml", draw("n_years"))
  expect_error(load_param_draws("m_nodosa", root = root2, check_sources = FALSE), "drawable run setting")
})

# ---- Timing: arrival year and ramp ------------------------------------------------
test_that("rust_arrival_year delays the start and never brings it forward", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 2010)
  expect_identical(apply_rust_timing(p), p)
  expect_identical(apply_rust_timing(modifyList(p, list(rust_arrival_year = 2018)))$rust_start_year, 2018)
  expect_identical(apply_rust_timing(modifyList(p, list(rust_arrival_year = 2000)))$rust_start_year, 2010)
  expect_identical(apply_rust_timing(modifyList(p, list(rust_start_year = Inf, rust_arrival_year = 2018)))$rust_start_year, Inf)
})

test_that("the ramp: zero before arrival, a linear rise, then the full draw, with no extra draws", {
  set.seed(1)
  p <- make_toy_params(rust_start_year = 10, annual_env_alpha = 14, annual_env_beta = 6)
  q <- modifyList(p, list(rust_ramp_years = 4))
  expect_identical(vapply(c(5, 10, 11, 12, 13, 14, 30), rust_ramp_factor, 1, params = q),
                   c(1, 0.25, 0.5, 0.75, 1, 1, 1))
  for (t in c(5, 10, 12, 20)) {
    set.seed(t); a <- draw_annual_env(p, t)
    set.seed(t); b <- draw_annual_env(q, t)
    expect_identical(b, a * rust_ramp_factor(q, t))
    expect_identical(runif(1), { set.seed(t); draw_annual_env(p, t); runif(1) })
  }
  expect_identical(draw_annual_env(q, 9), 0)
})

test_that("arrival and ramp run through the simulation; absent, runs are unchanged", {
  set.seed(1)
  p <- make_toy_params(n_years = 30L, rust_start_year = 5, annual_env_alpha = 14, annual_env_beta = 6)
  p$rust_pressure <- 1
  a <- run_simulation(p)
  b <- run_simulation(modifyList(p, list(rust_ramp_years = 0)))
  expect_identical(a, b)
  late <- run_simulation(modifyList(p, list(rust_arrival_year = 20)))
  expect_false(any(late$census$rust_active[late$census$year < 20]))
  expect_true(all(late$census$rust_active[late$census$year >= 20]))
  expect_identical(late$census[late$census$year < 20, "N_IUCN"],
                   run_simulation(modifyList(p, list(rust_start_year = 20)))$census[late$census$year < 20, "N_IUCN"])
})

test_that("the timing fields are optional species settings", {
  set.seed(1)
  schema <- read_config_yaml(file.path(repo_root, "config", "schema.yaml"))
  expect_identical(schema$params$rust_arrival_year$layer, "species")
  expect_identical(schema$params$rust_arrival_year$default, -Inf)
  expect_identical(schema$params$rust_ramp_years$default, 0)
  expect_true(isTRUE(schema$params$rust_gradient$drawable))
  m <- load_species("m_nodosa", root = repo_root)
  expect_false(any(c("rust_arrival_year", "rust_ramp_years", "rust_gradient") %in% names(m)))
  root <- cfg_copy(); on.exit(unlink(root, recursive = TRUE))
  edit_file(root, "species/m_nodosa.yaml", add_params(
    '  rust_ramp_years: {value: 0, prior: true, min: 0, max: 5, reason: "test"}'))
  expect_identical(load_species("m_nodosa", root = root, check_sources = FALSE), m)
})
