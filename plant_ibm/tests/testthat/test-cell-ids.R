# tests/testthat/test-cell-ids.R
#
# The master grid and its cell ids (spatial/cell_ids.R; notes T13), the shared
# fire tables (spatial/data/fire/), and the NSW clipping fix (notes T14).

source(file.path(proj_root, "spatial", "cell_ids.R"), local = TRUE)
repo_root <- normalizePath(file.path(proj_root, ".."))
data_root <- file.path(proj_root, "spatial", "data")
grid      <- master_grid(load_landscape_spatial(root = repo_root))

test_that("the master grid covers Australia in whole 0.1-degree cells", {
  set.seed(1)
  expect_identical(grid[c("xmin", "xmax", "ymin", "ymax", "res")],
                   list(xmin = 110, xmax = 160, ymin = -45, ymax = -9, res = 0.1))
  expect_identical(c(grid$ncol, grid$nrow), c(500L, 360L))
})

test_that("cell_from_xy() numbers cells as terra does, and cell_centre() inverts it", {
  set.seed(1)
  skip_if_not_installed("terra")
  r <- master_raster(grid)
  x <- runif(2e4, 110, 160); y <- runif(2e4, -45, -9)
  expect_identical(cell_from_xy(x, y, grid), as.integer(terra::cellFromXY(r, cbind(x, y))))
  ids <- sample(grid$ncol * grid$nrow, 500)
  cc  <- cell_centre(ids, grid)
  expect_equal(unname(as.matrix(cc[, c("lon_centre", "lat_centre")])), unname(terra::xyFromCell(r, ids)))
  expect_identical(cell_from_xy(cc$lon_centre, cc$lat_centre, grid), as.integer(ids))
  expect_true(is.na(cell_from_xy(109.9, -30, grid)) && is.na(cell_from_xy(150, -8.9, grid)))
})

test_that("M. nodosa's cells carry their master-grid ids; renumbering kept the centres", {
  set.seed(1)
  g <- read.csv(file.path(data_root, "m_nodosa", "grid.csv"))
  expect_identical(g$cell_id, cell_from_xy(g$lon_centre, g$lat_centre, grid))
  s <- read.csv(file.path(data_root, "m_nodosa", "site_params.csv"))
  expect_identical(s$cell_id, cell_from_xy(s$lon_centre, s$lat_centre, grid))
  # The pre-T13 ids (range-box raster) still map one to one onto these cells.
  old <- old_box_cell_id(g$lon_centre, g$lat_centre,
                         list(lon_min = 148, lon_max = 155, lat_min = -38, lat_max = -24), 0.1)
  expect_false(anyNA(old))
  expect_identical(anyDuplicated(old), 0L)
})

test_that("the shared fire tables cover every species' cells, and only those", {
  set.seed(1)
  fc  <- read.csv(file.path(data_root, "fire", "cells.csv"))
  all <- all_species_cells(data_root)
  expect_setequal(fc$cell_id, all$cell_id)
  expect_identical(fc$cell_id, cell_from_xy(fc$lon_centre, fc$lat_centre, grid))
  for (f in Sys.glob(file.path(data_root, "fire", "cell_years_*.csv"))) {
    cy <- read.csv(f)
    expect_true(all(cy$cell_id %in% fc$cell_id), info = basename(f))
    st <- if (grepl("nsw", basename(f))) "NSW" else "Qld"
    expect_true(all(fc$state[match(cy$cell_id, fc$cell_id)] == st), info = basename(f))
  }
})

test_that("load_fire_cells() gives the species' cells in grid order, and stops if any are missing", {
  set.seed(1)
  skip_if_not_installed("sf")
  grid_csv <- file.path(data_root, "m_nodosa", "grid.csv")
  cells_csv <- file.path(data_root, "fire", "cells.csv")
  cells <- load_fire_cells(grid_csv, cells_csv)
  expect_identical(cells$cell_id, read.csv(grid_csv)$cell_id)
  short <- tempfile(fileext = ".csv")
  on.exit(unlink(short))
  write.csv(read.csv(cells_csv)[-1, ], short, row.names = FALSE)
  expect_error(load_fire_cells(grid_csv, short), "1 of the species' 375 cells are not in")
})

test_that("the loader rejects a range box outside the master grid, or a grid of part cells", {
  set.seed(1)
  root <- tempfile("mnsim_cfg_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  for (d in c("config", "species", "data")) file.copy(file.path(repo_root, d), root, recursive = TRUE)
  f <- file.path(root, "config", "landscape.yaml")
  txt <- readLines(f)
  i <- grep("^    lon_min:", txt)
  txt[i + 1L] <- "      value: 149"
  writeLines(txt, f)
  expect_error(load_species_spatial("m_nodosa", root = root), "must lie inside the master grid", fixed = TRUE)
  txt[i + 1L] <- "      value: 110.05"
  writeLines(txt, f)
  expect_error(load_landscape_spatial(root = root), "whole number of grid_deg cells", fixed = TRUE)
})

test_that("cell_year_fractions keeps each cell's own geometry when a cell has no land (notes T14)", {
  set.seed(1)
  skip_if_not_installed("sf")
  sf::sf_use_s2(FALSE)
  sq <- function(x0, y0, w, h) sf::st_polygon(list(matrix(c(x0, y0, x0 + w, y0, x0 + w, y0 + h,
                                                            x0, y0 + h, x0, y0), ncol = 2, byrow = TRUE)))
  # Three cells in a row; the middle one has no land in the state.
  cells <- sf::st_sf(cell_id = c(11L, 12L, 13L),
                     geometry = sf::st_sfc(sq(150, -30, 0.1, 0.1), sq(150.1, -30, 0.1, 0.1),
                                           sq(150.2, -30, 0.1, 0.1), crs = 4326))
  land  <- sf::st_sf(state = "NSW", geometry = sf::st_sfc(sq(149, -31, 1.1, 2), sq(150.2, -31, 1, 2), crs = 4326))
  ca <- sf::st_transform(cells, AREA_CRS)
  la <- sf::st_union(sf::st_transform(land, AREA_CRS))
  cells$area_cell_m2 <- as.numeric(sf::st_area(ca))
  cells$area_land_m2 <- c(cells$area_cell_m2[1], 0, cells$area_cell_m2[3])
  # One fire burning half of the third cell only.
  fires <- sf::st_sf(year = 2005L, geometry = sf::st_sfc(sq(150.2, -30, 0.05, 0.1), crs = 4326))
  r <- suppressMessages(cell_year_fractions(fires, cells, 2000:2010, "land", state_poly = land, verbose = FALSE))
  expect_identical(r$cell_id, 13L)            # before the fix: labelled 12, with cell 12's area
  expect_equal(r$frac_burned, 0.5, tolerance = 0.01)
})
