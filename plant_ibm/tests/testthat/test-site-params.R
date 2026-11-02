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
