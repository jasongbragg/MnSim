# spatial/fire/R/bbox.R -- the area of the Queensland fire downloads
# (docs/methodological_notes.md T10).
#
# The Qld fire data are a shared landscape layer, so one download serves every
# species. Its area is worked out from the grids of all species
# (spatial/data/<sp>/grid.csv), not typed in:
#   qld_cell_extent()  extent (cell edges) of the cells that may need Qld fire
#                      data: those north of QLD_BORDER_LAT - BORDER_BAND. Cells
#                      are assigned to a state by the state polygons, not by
#                      latitude, so the band makes this a safe superset.
#   qld_bbox()         that extent plus a margin (config/landscape.yaml,
#                      spatial$qld_bbox_margin): the area requested.
# Each finished download records the area it covers in '<file>.bbox'.
# check_qld_coverage() stops if a download does not cover qld_cell_extent(),
# e.g. after a more northern species was added; downloads from before this
# check have no record (they used xmin 147.9, ymin -29.3, xmax 154.2,
# ymax -23.9) and give a warning instead; record_download_bbox() can add one.
#
# Needs spatial/fire/fire_config.R (QLD_BORDER_LAT, BORDER_BAND, GRID_DEG, LANDSCAPE).

species_grid_files <- function(data_root = "spatial/data")
  Sys.glob(file.path(data_root, "*", "grid.csv"))

qld_cell_extent <- function(grids = species_grid_files(), border_lat = QLD_BORDER_LAT,
                            band = BORDER_BAND, deg = GRID_DEG) {
  if (!length(grids)) stop("no species grids (spatial/data/<sp>/grid.csv): run spatial/01_get_melaleuca_occurrences.R")
  g <- do.call(rbind, lapply(grids, function(f) utils::read.csv(f)[, c("lon_centre", "lat_centre")]))
  q <- g[g$lat_centre > border_lat - band, ]
  if (!nrow(q)) return(NULL)   # no species reaches Queensland
  h <- deg / 2
  round(c(xmin = min(q$lon_centre) - h, ymin = min(q$lat_centre) - h,
          xmax = max(q$lon_centre) + h, ymax = max(q$lat_centre) + h), 6)
}

qld_bbox <- function(extent = qld_cell_extent(), margin = LANDSCAPE$qld_bbox_margin) {
  if (is.null(extent)) stop("no species cells near or in Queensland: nothing to download")
  round(extent + c(-1, -1, 1, 1) * margin, 6)
}

bbox_file <- function(out) paste0(out, ".bbox")

record_download_bbox <- function(out, bbox) {
  utils::write.csv(as.data.frame(as.list(bbox[c("xmin", "ymin", "xmax", "ymax")])),
                   bbox_file(out), row.names = FALSE)
  invisible(bbox)
}

read_download_bbox <- function(out) {
  f <- bbox_file(out)
  if (!file.exists(f)) return(NULL)
  unlist(utils::read.csv(f)[1, c("xmin", "ymin", "xmax", "ymax")])
}

# Does 'outer' contain 'inner'? (both c(xmin, ymin, xmax, ymax))
bbox_contains <- function(outer, inner)
  outer[["xmin"]] <= inner[["xmin"]] && outer[["ymin"]] <= inner[["ymin"]] &&
  outer[["xmax"]] >= inner[["xmax"]] && outer[["ymax"]] >= inner[["ymax"]]

# Stop if any existing download misses part of 'need'; warn for downloads
# without a recorded area. Returns the files checked, invisibly.
check_qld_coverage <- function(files, need = qld_cell_extent()) {
  files <- files[file.exists(files)]
  if (is.null(need) || !length(files)) return(invisible(files))
  rec  <- lapply(files, read_download_bbox)
  none <- files[vapply(rec, is.null, logical(1))]
  if (length(none))
    warning(length(none), " Qld download(s) have no recorded area, so their coverage of ",
            paste(sprintf("%.2f", need), collapse = ", "), " (xmin, ymin, xmax, ymax) is not checked: ",
            paste(basename(none), collapse = ", "), call. = FALSE)
  short <- files[vapply(rec, function(b) !is.null(b) && !bbox_contains(b, need), logical(1))]
  if (length(short))
    stop("These Qld downloads do not cover all cells that need Qld fire data (",
         paste(sprintf("%.2f", need), collapse = ", "), "): ", paste(basename(short), collapse = ", "),
         ". Delete them and rerun spatial/fire/01_download_qld.R.", call. = FALSE)
  invisible(files)
}
