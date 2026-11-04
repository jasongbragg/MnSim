# tests/testthat/test-site-params.R
#
# One rust-pressure anchor shared by params.R, the cell table and the Pegg
# validation (docs/methodological_notes.md, "Rust pressure").

.anchor_env <- new.env()
sys.source(file.path(proj_root, "params.R"), .anchor_env)
.A <- .anchor_env$RUST_ANCHOR

test_that("rust_pressure_from_sdm gives the anchor cell exactly mu, whatever the SDM scale", {
  sdm <- c(0.2, 0.5, 0.8, 0.4)
  for (scale in c(1, 0.5, 0.01)) {
    rp <- rust_pressure_from_sdm(sdm * scale, anchor_index = 3, anchor = .A)
    expect_equal(rp$sdm_rel[3], 1)
    expect_equal(rp$annual_env_mean[3], .A$mu)
    expect_equal(rp$annual_env_alpha[3] + rp$annual_env_beta[3], .A$conc)
    expect_equal(rp$annual_env_mean[1], .A$mu * 0.2 / 0.8)
  }
})

test_that("rust_pressure_from_sdm applies the cap and the numerical floor", {
  rp <- rust_pressure_from_sdm(c(0, 0.1, 2), anchor_index = 2, anchor = .A)
  expect_equal(rp$annual_env_mean[3], .A$cap)
  expect_equal(rp$annual_env_mean[1], 1 / (.A$conc * 10))
  expect_error(rust_pressure_from_sdm(c(0, 0.5), anchor_index = 1, anchor = .A), "anchor")
})

test_that("params.R default rust pressure is the anchor's Beta", {
  p <- .anchor_env$get_default_params()
  expect_equal(p$annual_env_alpha, .A$mu * .A$conc)
  expect_equal(p$annual_env_beta, (1 - .A$mu) * .A$conc)
})

test_that("the cell table's anchor cell carries the params.R default rust pressure", {
  s <- load_site_params(drop_incomplete = FALSE)
  a <- anchor_site_row(s, .A)
  expect_true(isTRUE(a$rust_anchor_cell))
  expect_equal(a$sdm_rel, 1)
  p <- .anchor_env$get_default_params()
  expect_equal(a$annual_env_alpha, p$annual_env_alpha, tolerance = 1e-3)
  expect_equal(a$annual_env_beta,  p$annual_env_beta,  tolerance = 1e-3)
  expect_equal(sum(s$rust_anchor_cell %in% TRUE), 1L)
})

test_that("Pegg runs pick up the anchor cell, and site = NULL leaves params alone", {
  base <- make_toy_params(); base$annual_env_alpha <- 40; base$annual_env_beta <- 8
  p <- pegg_site_params(base, "anchor")
  a <- anchor_site_row(load_site_params(drop_incomplete = FALSE), .A)
  expect_equal(p$annual_env_alpha, a$annual_env_alpha)
  expect_equal(p$pegg_site_cell, a$cell_id)
  expect_identical(pegg_site_params(base, NULL), base)
})

test_that("extract_sdm_cells averages over the cell's land, filling cells whose centre is at sea", {
  skip_if_not_installed("terra")
  r <- terra::rast(xmin = 150, xmax = 151, ymin = -31, ymax = -30, res = 0.025, crs = "EPSG:4326")
  terra::values(r) <- 800                                     # biomod2 0-1000 scale
  xy <- terra::crds(r, na.rm = FALSE)
  r[xy[, 1] > 150.52 & xy[, 2] < -30.5] <- NA                 # "sea" over cell 2's centre
  r[xy[, 1] < 150.5 & xy[, 2] > -30.5 & xy[, 1] < 150.25] <- 400   # west half of cell 1 lower
  lon <- c(150.25, 150.55); lat <- c(-30.25, -30.55)
  ctr <- extract_sdm_cells(r, lon, lat, "centre")
  avg <- extract_sdm_cells(r, lon, lat, "cell_mean")
  expect_true(is.na(ctr[2]))
  expect_equal(avg[2], 0.8)                                   # land part only
  expect_equal(avg[1], 0.6, tolerance = 1e-6)                 # half 0.4, half 0.8
})

test_that("apply_sdm_to_sites keeps the anchor cell at exactly mu with cell means", {
  skip_if_not_installed("terra")
  r <- terra::rast(xmin = 152.5, xmax = 153.5, ymin = -29.5, ymax = -28.5, res = 0.025, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  terra::values(r) <- 0.5 + 0.4 * (xy[, 2] + 29.5)             # rises northwards
  s <- data.frame(cell_id = 1:3, lon_centre = c(153.25, 152.95, 153.05),
                  lat_centre = c(-29.15, -28.75, -29.35))
  out <- suppressMessages(apply_sdm_to_sites(s, r, .A))
  expect_true(out$rust_anchor_cell[1])
  expect_equal(out$annual_env_mean[1], .A$mu)
  expect_gt(out$annual_env_mean[2], out$annual_env_mean[3])
})
