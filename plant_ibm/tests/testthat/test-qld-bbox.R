# tests/testthat/test-qld-bbox.R
#
# The Queensland fire download area, worked out from the species grids
# (spatial/fire/R/bbox.R; methodological notes T10).

source(file.path(proj_root, "spatial/fire/R/bbox.R"), local = TRUE)
grids_root <- file.path(proj_root, "spatial/data")

test_that("qld_cell_extent() covers the M. nodosa cells near or in Queensland", {
  set.seed(1)
  ext <- qld_cell_extent(species_grid_files(grids_root))
  expect_equal(ext, c(xmin = 148.10, ymin = -29.20, xmax = 153.60, ymax = -24.10))
  # Inside the hand-set box the existing downloads were made with, so they still cover the cells.
  expect_true(bbox_contains(c(xmin = 147.9, ymin = -29.3, xmax = 154.2, ymax = -23.9), ext))
})

test_that("qld_bbox() adds the configured margin", {
  set.seed(1)
  ext <- c(xmin = 148.10, ymin = -29.20, xmax = 153.60, ymax = -24.10)
  expect_identical(LANDSCAPE$qld_bbox_margin, 0.25)
  expect_equal(qld_bbox(ext), c(xmin = 147.85, ymin = -29.45, xmax = 153.85, ymax = -23.85))
  expect_error(qld_bbox(NULL), "nothing to download")
})

test_that("a second, more northern species widens the extent", {
  set.seed(1)
  d <- tempfile("grids_"); dir.create(file.path(d, "sp_north"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))
  file.copy(file.path(grids_root, "m_nodosa"), d, recursive = TRUE)
  utils::write.csv(data.frame(cell_id = 1:2, lon_centre = c(145.75, 146.05), lat_centre = c(-17.25, -18.05)),
                   file.path(d, "sp_north", "grid.csv"), row.names = FALSE)
  ext <- qld_cell_extent(species_grid_files(d))
  expect_equal(ext, c(xmin = 145.70, ymin = -29.20, xmax = 153.60, ymax = -17.20))
})

test_that("check_qld_coverage() stops on a short download and warns when the area is unknown", {
  set.seed(1)
  d <- tempfile("dl_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  f <- file.path(d, c("ok.gpkg", "short.gpkg", "old.gpkg"))
  file.create(f)
  need <- c(xmin = 148.10, ymin = -29.20, xmax = 153.60, ymax = -24.10)
  record_download_bbox(f[1], qld_bbox(need))
  record_download_bbox(f[2], c(xmin = 148.5, ymin = -29.3, xmax = 154.2, ymax = -23.9))
  expect_identical(read_download_bbox(f[1]), qld_bbox(need))
  expect_silent(check_qld_coverage(f[1], need))
  expect_error(check_qld_coverage(f[1:2], need), "do not cover all cells.*short.gpkg")
  expect_warning(check_qld_coverage(f[c(1, 3)], need), "no recorded area.*old.gpkg")
  expect_silent(check_qld_coverage(file.path(d, "missing.gpkg"), need))
})
