# tests/testthat/test-fire-regime.R
#
# Fire pipeline (docs/methodological_notes.md D1-D5): burnt fraction per
# cell-year, the empirical-Bayes fire-rate estimator, the per-cell site
# table, and how runs pick a fire regime.

# ---- Estimator ------------------------------------------------------------------

.sim_counts <- function(n, alpha, beta, T, seed) {
  set.seed(seed)
  T <- rep_len(T, n)
  lam <- stats::rgamma(n, alpha, beta)
  list(k = stats::rpois(n, lam * T), T = T)
}

test_that("the NB maximum-likelihood prior recovers the Gamma that generated the counts", {
  d <- .sim_counts(3000, alpha = 0.7, beta = 17, T = c(38, 11), seed = 1)
  p <- fit_prior_nb(d$k, d$T)
  expect_equal(unname(p["alpha"] / p["beta"]), 0.7 / 17, tolerance = 0.08)
  expect_equal(unname(p["alpha"]), 0.7, tolerance = 0.25)
})

test_that("the NB prior uses cells with no fire years; the old method of moments ignores them", {
  d <- .sim_counts(400, alpha = 1, beta = 20, T = 38, seed = 2)
  more0 <- list(k = c(d$k, rep(0L, 200)), T = c(d$T, rep(38, 200)))
  nb1 <- fit_prior_nb(d$k, d$T); nb2 <- fit_prior_nb(more0$k, more0$T)
  mm1 <- fit_prior_mm02(d$k, d$T); mm2 <- fit_prior_mm02(more0$k, more0$T)
  expect_lt(nb2["alpha"] / nb2["beta"], nb1["alpha"] / nb1["beta"])   # more zeros: lower mean rate
  expect_equal(mm1[c("alpha", "beta")], mm2[c("alpha", "beta")])     # unchanged
})

test_that("mm02 is the pre-Nov-2026 02 formula (method of moments on k/T, k > 0)", {
  k <- c(0, 0, 1, 2, 3, 1, 5, 2, 1, 4, 2, 3, 0); T <- 55
  lam <- (k / T)[k > 0]
  p <- fit_prior_mm02(k, rep(T, length(k)))
  expect_equal(unname(p["alpha"]), mean(lam)^2 / var(lam))
  expect_equal(unname(p["beta"]),  mean(lam) / var(lam))
})

.toy_cells <- function() data.frame(cell_id = 1:6, state = rep(c("NSW", "Qld"), each = 3))
.toy_cy <- function() data.frame(
  cell_id     = c(1, 1, 1, 2, 4, 4, 5),
  year        = c(2000, 2005, 2010, 2001, 2015, 2018, 2016),
  frac_burned = c(0.5, 0.2, 0.05, 0.3, 0.12, 0.4, 0.08))

test_that("estimate_fire_regime: fire years, burnt fractions, per-cell record length", {
  cells <- .toy_cells()
  T_obs <- setNames(c(38, 38, 38, 11, 11, 11), cells$cell_id)
  r <- estimate_fire_regime(.toy_cy(), cells, T_obs, prior = "fixed")
  expect_equal(r$n_wildfire_years_all, c(3, 1, 0, 2, 1, 0))
  expect_equal(r$n_wildfire_years_sig, c(2, 1, 0, 2, 0, 0))       # 0.05 and 0.08 below 10%
  expect_equal(r$mean_frac_burned[1], mean(c(0.5, 0.2)))
  expect_equal(r$T_obs, unname(T_obs))
  expect_equal(r$fire_prob_annual, pmax(FIRE_MIN, (PRIOR_ALPHA + r$n_wildfire_years_sig) / (PRIOR_BETA + r$T_obs)))
  # cells without a fire year get the median burnt fraction of cells with one
  expect_equal(r$fire_p_fimp[3], median(r$mean_frac_burned, na.rm = TRUE))
})

test_that("fire_min = 0 removes the rate floor (planned burns), prior_group fits one prior per group", {
  cells <- .toy_cells()
  T_obs <- setNames(rep(11, 6), cells$cell_id)
  r0 <- estimate_fire_regime(.toy_cy(), cells, T_obs, prior = "fixed", fire_min = 0)
  expect_true(all(r0$fire_prob_annual >= 0))
  d <- .sim_counts(600, alpha = 0.7, beta = 17, T = 38, seed = 3)
  big <- data.frame(cell_id = seq_along(d$k), state = rep(c("NSW", "Qld"), length.out = length(d$k)))
  cy <- do.call(rbind, lapply(which(d$k > 0), function(i)
    data.frame(cell_id = i, year = 1987 + seq_len(min(d$k[i], 38)) - 1L, frac_burned = 0.3)))
  r <- estimate_fire_regime(cy, big, setNames(d$T, big$cell_id), prior_group = big$state)
  pr <- attr(r, "priors")
  expect_setequal(pr$group, c("NSW", "Qld"))
  expect_true(all(pr$type == "nb_ml"))
})

test_that("prior_lr_test separates different regimes and accepts a shared one", {
  a <- .sim_counts(500, 0.7, 17, 38, seed = 4); b <- .sim_counts(500, 0.7, 17, 11, seed = 5)
  same <- prior_lr_test(c(a$k, b$k), c(a$T, b$T), rep(c("NSW", "Qld"), each = 500))
  expect_gt(same[["p"]], 0.001)
  c2 <- .sim_counts(500, 0.7, 3, 11, seed = 6)                      # ~6x the fire rate
  diff <- prior_lr_test(c(a$k, c2$k), c(a$T, c2$T), rep(c("NSW", "Qld"), each = 500))
  expect_lt(diff[["p"]], 1e-6)
})

# ---- Site table -------------------------------------------------------------------

.write_toy_tables <- function(dir) {
  cy <- .toy_cy()
  for (nm in c("nsw_wild", "nsw_planned", "nsw_all", "qld_wild", "qld_planned", "qld_all")) {
    d <- if (grepl("planned", nm)) cy[cy$cell_id %in% c(1, 4), ] else cy
    if (grepl("all", nm)) d$frac_burned <- pmin(1, d$frac_burned + 0.1)
    utils::write.csv(d, file.path(dir, paste0("cell_years_", nm, ".csv")), row.names = FALSE)
  }
}

test_that("site_fire_regime writes wildfire under the standard names, all/planned with suffixes", {
  dir <- tempfile("fire"); dir.create(dir); .write_toy_tables(dir)
  regs <- lapply(FIRE_REGIMES, function(sc) { sc$prior <- "fixed"; sc })
  out <- site_fire_regime(.toy_cells(), regimes = regs, dir = dir)
  expect_true(all(setdiff(SITE_FIRE_COLUMNS(), "n_wildfire_events") %in% names(out)))
  expect_equal(out$fire_T_obs, rep(c(length(NSW_YEARS), length(QLD_YEARS)), each = 3))
  expect_true(all(grepl("^NSW", out$fire_record[1:3])) && all(grepl("^Qld", out$fire_record[4:6])))
  # Qld years outside QLD_YEARS are not counted (cell 4's 2015/2018 are inside; 2016 is below 10%)
  expect_equal(out$n_wildfire_years_sig, c(2, 1, 0, 2, 0, 0))
  # planned: no rate floor; cells never planned-burnt get the prior mean only
  expect_true(all(out$fire_prob_annual_planned[c(2, 3, 5, 6)] < out$fire_prob_annual[c(2, 3, 5, 6)] + 1e-12))
  expect_true(all(out$n_fire_years_all >= out$n_wildfire_years_sig))
})

# ---- Burnt fraction per cell-year (geometry) -----------------------------------------

test_that("cell_year_fractions: union within a year, one row per year, land or whole-cell base", {
  skip_if_not_installed("sf")
  sq <- function(x0, y0, w, h) sf::st_polygon(list(matrix(c(x0, y0, x0 + w, y0, x0 + w, y0 + h,
                                                            x0, y0 + h, x0, y0), ncol = 2, byrow = TRUE)))
  sf::sf_use_s2(FALSE)
  cell <- sf::st_sf(cell_id = 1L, geometry = sf::st_sfc(sq(150, -30, 0.1, 0.1), crs = 4326))
  land <- sf::st_sf(state = "NSW", geometry = sf::st_sfc(sq(149, -31, 1.05, 2), crs = 4326))  # west half of the cell
  a_cell <- as.numeric(sf::st_area(sf::st_transform(cell, AREA_CRS)))
  a_land <- as.numeric(sf::st_area(suppressWarnings(sf::st_intersection(
    sf::st_geometry(sf::st_transform(cell, AREA_CRS)), sf::st_geometry(sf::st_transform(land, AREA_CRS))))))
  cell$area_cell_m2 <- a_cell; cell$area_land_m2 <- a_land
  fires <- sf::st_sf(year = c(2000L, 2000L, 2003L),
                     geometry = sf::st_sfc(sq(150, -30, 0.05, 0.1),      # west half
                                           sq(150, -30, 0.05, 0.1),      # same area again, same year
                                           sq(150, -30, 0.1, 0.05),      # south half, later year
                                           crs = 4326))
  r_cell <- suppressMessages(cell_year_fractions(fires, cell, 1999:2005, "cell", verbose = FALSE))
  expect_equal(r_cell$year, c(2000, 2003))
  expect_equal(r_cell$frac_burned, c(0.5, 0.5), tolerance = 0.01)      # overlap not double counted
  r_land <- suppressMessages(cell_year_fractions(fires, cell, 1999:2005, "land", state_poly = land,
                                                 verbose = FALSE))
  expect_equal(r_land$frac_burned[1], 1, tolerance = 0.01)            # the burnt half is all the land
  expect_equal(r_land$frac_burned[2], 0.5, tolerance = 0.01)
  expect_equal(nrow(suppressMessages(cell_year_fractions(fires, cell, 2001:2002, "cell", verbose = FALSE))), 0)
})

# ---- Runs: choosing the fire regime ---------------------------------------------------

.site_row <- function() data.frame(
  cell_id = 1, lon_centre = 153, lat_centre = -27, fire_state = "Qld",
  fire_prob_annual = 0.04, fire_p_fimp = 0.2, sd_frac_burned = 0.1,
  fire_prob_annual_all = 0.06, fire_p_fimp_all = 0.25, sd_frac_burned_all = NA,
  annual_env_alpha = 14, annual_env_beta = 6)

test_that("build_site_params uses wildfire by default, all fire on request, and scales the rate", {
  b <- make_toy_params(); b$fire_size_sd <- 0.125
  w <- build_site_params(b, .site_row())
  expect_equal(c(w$fire_prob_annual, w$fire_p_fimp, w$fire_size_sd), c(0.04, 0.2, 0.1))
  a <- build_site_params(b, .site_row(), fire_regime = "all")
  expect_equal(c(a$fire_prob_annual, a$fire_p_fimp), c(0.06, 0.25))
  expect_equal(a$fire_size_sd, 0.125)                                   # NA: keep the base value
  h <- build_site_params(b, .site_row(), fire_scale = 0.5)
  expect_equal(h$fire_prob_annual, 0.02)
  expect_equal(build_site_params(b, .site_row(), fire_scale = 0)$fire_prob_annual, 0)
  old <- .site_row()[, c("fire_prob_annual", "fire_p_fimp", "sd_frac_burned", "annual_env_alpha", "annual_env_beta")]
  expect_error(build_site_params(b, old, fire_regime = "all"), "rebuild")
  expect_error(build_site_params(b, .site_row(), fire_regime = "planned"))
})

test_that("fire scenarios pair up with the default run by run (same cell, seed and draws)", {
  s <- load_range_sites(SITE_PARAMS_PATH); b <- make_toy_params()
  for (i in c(2, 9)) {
    d <- plan_range_run(i, s, b)
    h <- plan_range_run(i, s, b, fire_scale = 0.5)
    expect_identical(d$row[c("cell_id", "sim_seed", "resprout_rust_mortality")],
                     h$row[c("cell_id", "sim_seed", "resprout_rust_mortality")])
    expect_equal(h$params$fire_prob_annual, 0.5 * d$params$fire_prob_annual)
    expect_equal(d$row$fire_regime, "wild"); expect_equal(h$row$fire_scale, 0.5)
    expect_false(identical(d$row$run_id, h$row$run_id))
  }
})
