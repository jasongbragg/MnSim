# spatial/fire/R/cell_years.R -- burnt fraction of each grid cell in each year, the same
# way for both states.
#
# For a cell and year: union all of that year's burnt polygons (so overlaps
# are not double counted), optionally remove planned-burn areas, intersect
# with the cell (clipped to the cell's state), and divide by
#   area_base = "cell": the whole 0.1-degree cell (the pre-Nov-2026 02 convention), or
#   area_base = "land": the land part of the cell within its state.
# Areas in Australian Albers (EPSG:3577).

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

.polys_only <- function(g) {
  g <- g[!st_is_empty(g)]
  if (!length(g)) return(g)
  types <- as.character(st_geometry_type(g))
  if (any(types %in% c("GEOMETRYCOLLECTION")))
    g <- suppressWarnings(st_collection_extract(g, "POLYGON"))
  g
}

# fires: standard polygons (R/sources.R) for ONE state; cells: load_cells()
# rows for that state; remove: optional polygons (with year) whose area is
# subtracted from each year's fires (e.g. QPWS planned burns).
# Returns data.frame(cell_id, year, frac_burned) for every cell-year with
# any burnt area, plus attribute "years" = years examined.
cell_year_fractions <- function(fires, cells, years, area_base = c("land", "cell"),
                                remove = NULL, state_poly = NULL, verbose = TRUE) {
  area_base <- match.arg(area_base)
  fires <- st_transform(fires[fires$year %in% years, ], AREA_CRS)
  cells_a <- st_transform(cells, AREA_CRS)
  if (!is.null(state_poly)) {             # clip cells to their state's land
    sp <- st_union(st_transform(state_poly, AREA_CRS))
    clipped <- suppressWarnings(st_intersection(st_geometry(cells_a), sp))
    # keep cells with no land intersection as empty (rare)
    cells_geom <- clipped
  } else cells_geom <- st_geometry(cells_a)
  if (!is.null(remove)) remove <- st_transform(remove, AREA_CRS)
  base <- if (area_base == "cell") cells$area_cell_m2 else cells$area_land_m2

  out <- list()
  for (y in sort(unique(fires$year))) {
    fy <- fires[fires$year == y, ]
    g  <- .polys_only(st_geometry(fy))
    if (!length(g)) next
    u  <- st_union(st_make_valid(g))
    if (!is.null(remove)) {
      ry <- remove[remove$year == y, ]
      if (nrow(ry)) u <- suppressWarnings(st_difference(u, st_union(st_make_valid(st_geometry(ry)))))
    }
    if (length(u) == 0 || all(st_is_empty(u))) next
    hit <- which(lengths(st_intersects(cells_geom, u)) > 0)
    if (!length(hit)) next
    a <- vapply(hit, function(i) {
      inter <- suppressWarnings(st_intersection(cells_geom[i], u))
      if (!length(inter)) 0 else sum(as.numeric(st_area(.polys_only(inter))))
    }, numeric(1))
    keep <- a > 0
    if (any(keep))
      out[[length(out) + 1L]] <- data.frame(cell_id = cells$cell_id[hit[keep]], year = y,
                                            frac_burned = pmin(1, a[keep] / base[hit[keep]]))
    if (verbose) message(sprintf("    %d: %d cells burnt", y, sum(keep)))
  }
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(cell_id = integer(0), year = integer(0), frac_burned = numeric(0))
  attr(res, "years") <- years
  res
}

# ---- Queensland, one year at a time (memory-safe) ------------------------------
# Recent Qld years hold >100,000 scar polygons in the download area, so years
# are processed one at a time: read the year, keep only polygons touching an
# occupied cell, cut them into cell pieces, and union per cell. Returns the
# burnt fraction of each cell's Qld land for: all fire (scars), wild (scars
# minus that year's QPWS planned burns) and planned (QPWS planned burns alone).
qld_year_fractions <- function(y, cells_land, pb, scar_dir = file.path(FIRE_RAW_DIR, "qld_scars")) {
  empty <- data.frame(cell_id = integer(0), year = integer(0), frac_all = numeric(0),
                      frac_wild = numeric(0), frac_planned = numeric(0))
  f <- file.path(scar_dir, sprintf("qld_scars_%d.gpkg", y))
  if (!file.exists(f)) { warning("missing Qld scars for ", y); return(empty) }
  x <- st_read(f, quiet = TRUE)
  if (!nrow(x)) return(empty)
  x <- st_transform(st_geometry(x), AREA_CRS)
  footprint <- st_union(st_geometry(cells_land))
  x <- x[lengths(st_intersects(x, footprint)) > 0]
  if (!length(x)) return(empty)
  x <- st_make_valid(x)

  pieces <- function(g) {                   # g: sfc in AREA_CRS -> sf(cell_id, geometry)
    if (!length(g)) return(NULL)
    p <- suppressWarnings(st_intersection(st_sf(geometry = g), cells_land[, "cell_id"]))
    if (!nrow(p)) return(NULL)
    p[!st_is_empty(p), ]
  }
  sp <- pieces(x); rm(x)
  pp <- NULL
  if (!is.null(pb)) {
    pby <- st_transform(st_geometry(pb[pb$year == y, ]), AREA_CRS)
    if (length(pby)) pp <- pieces(st_make_valid(pby))
  }
  ids <- unique(c(if (!is.null(sp)) sp$cell_id, if (!is.null(pp)) pp$cell_id))
  if (!length(ids)) return(empty)
  base <- setNames(cells_land$area_land_m2, cells_land$cell_id)
  area_of <- function(g) if (!length(g) || all(st_is_empty(g))) 0 else sum(as.numeric(st_area(.polys_only(g))))
  rows <- lapply(ids, function(id) {
    u_s <- if (!is.null(sp) && any(sp$cell_id == id)) st_union(st_geometry(sp[sp$cell_id == id, ])) else NULL
    u_p <- if (!is.null(pp) && any(pp$cell_id == id)) st_union(st_geometry(pp[pp$cell_id == id, ])) else NULL
    a_all  <- if (is.null(u_s)) 0 else area_of(u_s)
    a_plan <- if (is.null(u_p)) 0 else area_of(u_p)
    a_wild <- if (is.null(u_s)) 0 else if (is.null(u_p)) a_all else
      area_of(suppressWarnings(st_difference(u_s, u_p)))
    b <- base[[as.character(id)]]
    data.frame(cell_id = id, year = y, frac_all = min(1, a_all / b),
               frac_wild = min(1, a_wild / b), frac_planned = min(1, a_plan / b))
  })
  do.call(rbind, rows)
}

# Cells clipped to their state's land, in AREA_CRS (for qld_year_fractions).
cells_on_land <- function(cells, state_poly) {
  ca <- st_transform(cells, AREA_CRS)
  sp <- st_union(st_transform(state_poly, AREA_CRS))
  g  <- suppressWarnings(st_intersection(st_geometry(ca), sp))
  ca <- ca[, c("cell_id", "area_land_m2")]
  st_geometry(ca) <- g
  ca[!st_is_empty(ca), ]
}
