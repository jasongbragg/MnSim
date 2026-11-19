# spatial/cell_ids.R -- the master grid and its cell ids (notes T13).
#
# Every species uses one master grid: a 0.1-degree raster over all of
# Australia (config/landscape.yaml, spatial$grid_extent and grid_deg). A
# cell's id is its number in that raster (terra's numbering: row by row from
# the north-west corner, starting at 1), so it is the same whatever species
# occupies the cell. That lets the fire cell-year tables in
# spatial/data/fire/ be shared: they cover the union of all species' cells
# (all_species_cells()), and each species takes its own.
#
# Before T13 a cell's id was its number in a raster spanning the species'
# range box, so ids differed between species. The cell centres are the
# same either way (both grids sit on 0.1-degree steps); only the numbers
# changed. old_box_cell_id() keeps the old numbering, for the one-off
# renumbering and its test.
#
# Needs R/species_config.R (load_landscape_spatial()).

master_grid <- function(land = load_landscape_spatial()) {
  e   <- land$grid_extent
  res <- land$grid_deg
  list(xmin = e$lon_min, xmax = e$lon_max, ymin = e$lat_min, ymax = e$lat_max, res = res,
       ncol = as.integer(round((e$lon_max - e$lon_min) / res)),
       nrow = as.integer(round((e$lat_max - e$lat_min) / res)))
}

# The master grid as a terra raster (spatial/01: assigning records to cells).
master_raster <- function(grid = master_grid()) {
  terra::rast(xmin = grid$xmin, xmax = grid$xmax, ymin = grid$ymin, ymax = grid$ymax,
              ncols = grid$ncol, nrows = grid$nrow, crs = "EPSG:4326")
}

# Cell id of each point (NA outside the grid). Same numbering as
# terra::cellFromXY(master_raster(), ...), without needing terra.
cell_from_xy <- function(lon, lat, grid = master_grid()) {
  col <- floor((lon - grid$xmin) / grid$res + 1e-9)
  row <- floor((grid$ymax - lat) / grid$res + 1e-9)
  ok  <- !is.na(col) & !is.na(row) & col >= 0 & col < grid$ncol & row >= 0 & row < grid$nrow
  out <- rep(NA_integer_, length(lon))
  out[ok] <- as.integer(row[ok] * grid$ncol + col[ok] + 1)
  out
}

# Centre of each cell id.
cell_centre <- function(cell_id, grid = master_grid()) {
  i <- cell_id - 1L
  data.frame(cell_id    = cell_id,
             lon_centre = round(grid$xmin + (i %% grid$ncol + 0.5) * grid$res, 6),
             lat_centre = round(grid$ymax - (i %/% grid$ncol + 0.5) * grid$res, 6))
}

# Every species' occupied cells (spatial/data/<sp>/grid.csv), once each: the
# cells the shared fire tables must cover.
all_species_cells <- function(data_root = "spatial/data") {
  files <- Sys.glob(file.path(data_root, "*", "grid.csv"))
  if (!length(files)) stop("no species grids (", data_root, "/<sp>/grid.csv): run spatial/01_get_melaleuca_occurrences.R")
  g <- do.call(rbind, lapply(files, function(f) utils::read.csv(f)[, c("cell_id", "lon_centre", "lat_centre")]))
  g <- g[!duplicated(g$cell_id), ]
  g[order(g$cell_id), , drop = FALSE]
}

# The pre-T13 numbering: cell number in a raster spanning a range box.
old_box_cell_id <- function(lon, lat, box, res) {
  g <- list(xmin = box$lon_min, ymax = box$lat_max, res = res,
            ncol = as.integer(round((box$lon_max - box$lon_min) / res)),
            nrow = as.integer(round((box$lat_max - box$lat_min) / res)))
  cell_from_xy(lon, lat, g)
}
