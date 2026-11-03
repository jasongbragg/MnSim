# spatial/fire/R/sources.R -- read each state's fire record into one standard form:
#   sf polygons with columns  year (calendar), type ("wildfire" | "planned" |
#   "unknown"), source; lon/lat (EPSG:4326).
# Everything downstream (R/cell_years.R) treats both states identically.

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

# NSW NPWS Fire History: one polygon per fire, Label "YYYY-YY Wildfire|Prescribed Burn".
# Year = first year of the label (a 2019-20 fire counts as 2019).
read_nsw_npws <- function(shp = NPWS_SHP, bbox = NULL) {
  x <- st_read(shp, quiet = TRUE) |> st_make_valid() |> st_transform(4326)
  x$year <- suppressWarnings(as.integer(substr(x$Label, 1L, 4L)))
  x$type <- dplyr::case_when(
    grepl("Wildfire",   x$Label, ignore.case = TRUE) ~ "wildfire",
    grepl("Prescribed", x$Label, ignore.case = TRUE) ~ "planned",
    TRUE ~ "unknown")
  x <- x[!is.na(x$year), ]
  if (!is.null(bbox)) {
    box <- st_as_sfc(st_bbox(bbox, crs = 4326))
    x <- x[suppressWarnings(st_intersects(x, box, sparse = FALSE)[, 1]), ]
  }
  x$source <- "NSW NPWS"
  x[, c("year", "type", "source")]
}

# Queensland fire scars (downloaded per year). Type unknown: wildfire,
# planned and agricultural burning are not distinguished.
read_qld_scars <- function(years, dir = file.path(FIRE_RAW_DIR, "qld_scars")) {
  pieces <- lapply(years, function(y) {
    f <- file.path(dir, sprintf("qld_scars_%d.gpkg", y))
    if (!file.exists(f)) { warning("missing Qld scars for ", y); return(NULL) }
    x <- st_read(f, quiet = TRUE)
    if (nrow(x) == 0) return(NULL)
    x <- st_make_valid(st_transform(x, 4326))
    st_geometry(x) <- "geometry"
    x$year <- as.integer(y)                 # the layer's year (annual composite)
    x$type <- "unknown"; x$source <- "Qld scars"
    x[, c("year", "type", "source")]
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) stop("No Queensland scar files found in ", dir)
  do.call(rbind, pieces)
}

# QPWS fire history: planned burns (PB) and wildfires (WF) on QPWS / state land.
# Year = calendar year of ignition (falling back to out / completed year),
# to match the calendar-year scar composites.
read_qpws <- function(f = file.path(FIRE_RAW_DIR, "qpws_fire_history.gpkg")) {
  x <- st_read(f, quiet = TRUE) |> st_make_valid() |> st_transform(4326)
  st_geometry(x) <- "geometry"
  ig <- x$ignitiondate
  if (is.numeric(ig)) ig <- as.POSIXct(ig / 1000, origin = "1970-01-01", tz = "UTC")  # epoch ms
  yr <- suppressWarnings(as.integer(format(as.POSIXct(ig), "%Y")))
  if (!is.null(x$outyear))       yr <- dplyr::coalesce(yr, as.integer(x$outyear))
  if (!is.null(x$completedyear)) yr <- dplyr::coalesce(yr, as.integer(x$completedyear))
  x$year <- yr
  x$type <- dplyr::case_when(toupper(x$type) == "PB" ~ "planned",
                             toupper(x$type) == "WF" ~ "wildfire",
                             TRUE ~ "unknown")
  x$source <- "QPWS"
  x <- x[!is.na(x$year), ]
  x[, c("year", "type", "source")]
}
