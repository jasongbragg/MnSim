# spatial/fire/R/grid.R -- the M. nodosa grid cells as polygons, with state and land area.

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

# State polygons (NSW, Qld) in lon/lat. ozmaps if installed (as 02 uses);
# otherwise a latitude split of the Australian land polygon (maps package).
state_polygons <- function() {
  if (requireNamespace("ozmaps", quietly = TRUE)) {
    st <- ozmaps::ozmap_states
    st <- st[st$NAME %in% c("New South Wales", "Queensland"), ]
    st$state <- ifelse(st$NAME == "Queensland", "Qld", "NSW")
    return(st_transform(st[, "state"], 4326))
  }
  message("ozmaps not installed: approximating the border at ", QLD_BORDER_LAT, " deg S")
  aus <- st_as_sf(maps::map("world", "Australia", plot = FALSE, fill = TRUE))
  aus <- st_make_valid(st_transform(aus, 4326))
  aus <- st_union(aus)
  north <- st_as_sfc(st_bbox(c(xmin = 137, ymin = QLD_BORDER_LAT, xmax = 155, ymax = -10), crs = 4326))
  south <- st_as_sfc(st_bbox(c(xmin = 140, ymin = -38, xmax = 155, ymax = QLD_BORDER_LAT), crs = 4326))
  st_sf(state = c("Qld", "NSW"),
        geometry = c(st_intersection(aus, north), st_intersection(aus, south)))
}

# Cell polygons with: state (by cell centre), whole-cell area, and the land
# area within that state (the base for burnt fractions when area_base = "land").
load_cells <- function(grid_csv = MN_GRID_CSV) {
  g <- read.csv(grid_csv)
  h <- GRID_DEG / 2
  polys <- lapply(seq_len(nrow(g)), function(i) {
    x <- g$lon_centre[i]; y <- g$lat_centre[i]
    st_polygon(list(matrix(c(x - h, y - h, x + h, y - h, x + h, y + h, x - h, y + h, x - h, y - h),
                           ncol = 2, byrow = TRUE)))
  })
  cells <- st_sf(g, geometry = st_sfc(polys, crs = 4326))
  states <- state_polygons()

  ctr <- st_as_sf(g, coords = c("lon_centre", "lat_centre"), crs = 4326)
  hit <- st_intersects(ctr, states)
  cells$state <- vapply(seq_along(hit), function(i)
    if (length(hit[[i]])) states$state[hit[[i]][1]] else NA_character_, "")
  # Coastal cells whose centre is in the sea: nearest state
  miss <- is.na(cells$state)
  if (any(miss)) cells$state[miss] <- states$state[st_nearest_feature(ctr[miss, ], states)]

  cells_a <- st_transform(cells, AREA_CRS)
  states_a <- st_transform(states, AREA_CRS)
  cells$area_cell_m2 <- as.numeric(st_area(cells_a))
  cells$area_land_m2 <- vapply(seq_len(nrow(cells)), function(i) {
    s <- states_a[states_a$state == cells$state[i], ]
    a <- suppressWarnings(st_intersection(st_geometry(cells_a[i, ]), st_geometry(s)))
    if (length(a) == 0) 0 else sum(as.numeric(st_area(a)))
  }, numeric(1))
  cells$land_frac <- cells$area_land_m2 / cells$area_cell_m2
  cells
}

# Cells within BORDER_BAND degrees of latitude of the border, both sides.
border_band <- function(cells, border_lat = QLD_BORDER_LAT, band = BORDER_BAND) {
  abs(cells$lat_centre - border_lat) <= band
}
