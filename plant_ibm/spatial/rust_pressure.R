# spatial/rust_pressure.R
#
# Maps a myrtle rust SDM onto the pathogen vertex of the disease triangle
# (annual_env Beta parameters) for every occupied grid cell, anchored at one
# site (RUST_ANCHOR in params.R, from config/landscape.yaml).
#
# ASSUMPTION (docs/methodological_notes.md, "Rust pressure"):
#   mean annual pressure in cell i = mu * s_i / s_anchor, capped at 'cap'
# where s is the SDM value extracted the same way for every cell: the
# area-weighted mean of the raster over the grid cell (its land, since
# the SDM is NA at sea; notes B4, C2), and the anchor cell is a fixed
# master-grid cell (rust_anchor cell_lon, cell_lat; notes T18): rust pressure is
# a property of place, so every species normalises to the SDM in the same
# cell, occupied by that species or not. (Before T18 it was the occupied cell
# nearest the anchor site, which for M. nodosa is the same cell.) Until Nov 2026 s was read at the cell centre,
# which left 16 coastal cells (centre at sea) without a value; that reading
# remains available as extract = "centre". Dividing by the anchor cell's own SDM value means:
#   - the anchor cell always gets exactly mu (0.70), for any SDM, so the
#     field calibration at Yarringully and the Pegg et al. validation survive
#     SDM updates;
#   - a new SDM changes only the relative pattern across the range.
# The normaliser is computed from the data each time; it is never typed in.
#
# Functions only; sourced by 02_cross_reference_layers.R,
# 02c_update_rust_pressure.R and the tests. Grid size: config/landscape.yaml.

source("R/species_config.R")
GRID_DEG <- load_landscape_spatial()$grid_deg   # also set by spatial/fire/fire_config.R

# Row of the anchor cell in a site table (centre within 1e-6 degrees of the
# anchor's cell_lon, cell_lat), or NA if the species does not occupy it. An
# anchor without a cell (older lists) falls back to the occupied cell nearest
# the anchor site.
anchor_cell_index <- function(lon, lat, anchor = RUST_ANCHOR) {
  if (is.null(anchor$cell_lon)) {
    d2 <- (lon - anchor$lon)^2 + ((lat - anchor$lat))^2
    return(which.min(d2))
  }
  k <- which(abs(lon - anchor$cell_lon) < 1e-6 & abs(lat - anchor$cell_lat) < 1e-6)
  if (length(k)) k[1] else NA_integer_
}

# Core mapping, on SDM values already extracted per cell. s_anchor: the SDM in
# the anchor cell (by default the cell at anchor_index).
rust_pressure_from_sdm <- function(sdm, anchor_index, anchor = RUST_ANCHOR,
                                   s_anchor = sdm[anchor_index]) {
  if (length(s_anchor) != 1L || is.na(s_anchor) || s_anchor <= 0)
    stop("SDM value at the anchor cell is missing or <= 0; cannot normalise.")
  floor_env <- 1 / (anchor$conc * 10)    # numerical floor: alpha >= 0.1
  rel  <- sdm / s_anchor
  mean <- pmin(anchor$cap, pmax(floor_env, anchor$mu * rel))
  rust_pressure_columns(rel, mean, anchor)
}

rust_pressure_columns <- function(rel, mean, anchor) data.frame(
  sdm_rel          = rel,
  annual_env_mean  = mean,
  annual_env_alpha = round(mean * anchor$conc, 4L),
  annual_env_beta  = round((1 - mean) * anchor$conc, 4L)
)

# Rust-pressure modes: how each cell's mean annual pressure is set
# (design_swappable_interfaces.md section 4.1, design_rust_pressure_modes.md;
# notes T12, T18). Each takes the site table, the SDM in the anchor cell and
# the anchor, and returns the columns of rust_pressure_from_sdm(). The species'
# rust_anchor 'mode' picks one (absent = "sdm"); the names are the schema's
# rust_anchor$mode choices (tested).
#   sdm           mu x SDM / SDM_anchor, capped
#   constant      mu (capped) in every cell, SDM or not; sdm_rel kept for reference
#   sdm_gradient  the sdm columns; the run adds its latitude shift
#                 (apply_rust_gradient() in runner/site_params.R)
RUST_PRESSURE_MODES <- list(
  sdm = function(sites, s_anchor, anchor)
    rust_pressure_from_sdm(sites$sdm_suitability, NA_integer_, anchor, s_anchor = s_anchor),
  constant = function(sites, s_anchor, anchor) {
    rel <- if (is.null(sites$sdm_suitability) || is.na(s_anchor)) rep(NA_real_, nrow(sites))
           else sites$sdm_suitability / s_anchor
    rust_pressure_columns(rel, rep(min(anchor$cap, anchor$mu), nrow(sites)), anchor)
  },
  sdm_gradient = function(sites, s_anchor, anchor)
    rust_pressure_from_sdm(sites$sdm_suitability, NA_integer_, anchor, s_anchor = s_anchor)
)

rust_mode <- function(anchor = RUST_ANCHOR) if (is.null(anchor$mode)) "sdm" else anchor$mode

cell_rust_pressure <- function(sites, s_anchor, anchor = RUST_ANCHOR) {
  mode <- rust_mode(anchor)
  f <- RUST_PRESSURE_MODES[[mode]]
  if (is.null(f)) stop("Unknown rust-pressure mode '", mode, "'. Must be one of: ",
                       paste(names(RUST_PRESSURE_MODES), collapse = ", "), ".")
  f(sites, s_anchor, anchor)
}

# Grid cell polygons (terra) from cell centres; deg = grid size (config/landscape.yaml).
cell_polygons <- function(lon, lat, deg = GRID_DEG) {
  h <- deg / 2
  terra::vect(lapply(seq_along(lon), function(i)
    rbind(c(lon[i] - h, lat[i] - h), c(lon[i] + h, lat[i] - h), c(lon[i] + h, lat[i] + h),
          c(lon[i] - h, lat[i] + h), c(lon[i] - h, lat[i] - h))),
    type = "polygons", crs = "EPSG:4326")
}

# SDM value per cell, on a 0-1 scale (biomod2 writes 0-1000).
#   "cell_mean": area-weighted mean over the raster cells inside the
#                grid cell, ignoring NA (sea) -- the default
#   "centre":    the raster value at the cell centre (pre-Nov-2026)
extract_sdm_cells <- function(sdm, lon, lat, extract = c("cell_mean", "centre"),
                              deg = GRID_DEG) {
  extract <- match.arg(extract)
  sdm <- terra::classify(sdm[[1]], cbind(-Inf, -1e30, NA))   # nodata written as a huge negative
  vals <- if (extract == "centre") {
    terra::extract(sdm, terra::vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"),
                                    crs = "EPSG:4326"))[, 2L]
  } else {
    terra::extract(sdm, cell_polygons(lon, lat, deg), fun = mean, na.rm = TRUE, exact = TRUE)[, 2L]
  }
  vals[is.nan(vals)] <- NA_real_
  if (any(!is.na(vals)) && max(vals, na.rm = TRUE) > 1.5) vals <- vals / 1000
  vals
}

# Add / refresh the rust columns of a site table (one row per occupied cell,
# columns lon_centre, lat_centre). If 'sdm' (a terra SpatRaster) is given the
# SDM is re-extracted per cell (extract_sdm_cells); otherwise the existing
# sdm_suitability column is used. Prints the anchor diagnostics.
apply_sdm_to_sites <- function(sites, sdm = NULL, anchor = RUST_ANCHOR,
                               neighbourhood_deg = 0.15, warn_rel_diff = 0.10,
                               extract = "cell_mean") {
  if (!is.null(sdm)) {
    old <- sites$sdm_suitability
    sites$sdm_suitability <- extract_sdm_cells(sdm, sites$lon_centre, sites$lat_centre, extract)
    message(sprintf("  SDM per cell (%s): %d of %d cells without a value",
                    extract, sum(is.na(sites$sdm_suitability)), nrow(sites)))
    if (!is.null(old) && any(is.na(old) & !is.na(sites$sdm_suitability)))
      message(sprintf("  %d cells newly filled (no value before)",
                      sum(is.na(old) & !is.na(sites$sdm_suitability))))
  }
  k <- anchor_cell_index(sites$lon_centre, sites$lat_centre, anchor)
  s_anchor <- if (!is.na(k)) sites$sdm_suitability[k]
    else if (!is.null(sdm)) extract_sdm_cells(sdm, anchor$cell_lon, anchor$cell_lat, extract)
    else if (identical(rust_mode(anchor), "constant")) NA_real_
    else stop("The species does not occupy the rust anchor cell (", anchor$cell_lon, ", ", anchor$cell_lat,
              "), so its SDM value must be read from the raster: supply 'sdm'.")
  rp <- cell_rust_pressure(sites, s_anchor, anchor)
  sites[names(rp)] <- rp
  sites$rust_anchor_cell <- seq_len(nrow(sites)) %in% k

  message(sprintf(
    "  Rust anchor: %s -> cell %s (%.2f, %.2f)%s, SDM %.4f -> pressure %.2f; mode %s",
    anchor$site,
    if (is.na(k)) "(not occupied)" else if ("cell_id" %in% names(sites)) sites$cell_id[k] else k,
    if (is.na(k)) anchor$cell_lon else sites$lon_centre[k],
    if (is.na(k)) anchor$cell_lat else sites$lat_centre[k],
    if (is.na(k)) " read from the raster" else "", s_anchor, anchor$mu, rust_mode(anchor)))
  if (is.na(k)) k <- NULL

  # Is the anchor cell typical of its surroundings? (single-cell noise guard)
  if (!is.null(sdm) && !is.null(k)) {
    box <- terra::ext(sites$lon_centre[k] - neighbourhood_deg, sites$lon_centre[k] + neighbourhood_deg,
                      sites$lat_centre[k] - neighbourhood_deg, sites$lat_centre[k] + neighbourhood_deg)
    nb <- terra::values(terra::crop(sdm, box), mat = FALSE)
    nb <- nb[!is.na(nb) & nb > -1e30]
    if (length(nb) && max(nb) > 1.5) nb <- nb / 1000
    if (length(nb)) {
      rel_diff <- abs(sites$sdm_suitability[k] - mean(nb)) / mean(nb)
      message(sprintf("  Anchor neighbourhood (+/-%.2f deg): mean SDM %.4f (anchor differs by %.1f%%)",
                      neighbourhood_deg, mean(nb), 100 * rel_diff))
      if (rel_diff > warn_rel_diff)
        warning(sprintf("Anchor cell SDM differs from its neighbourhood by %.0f%%: the anchor may be unrepresentative.",
                        100 * rel_diff))
    }
  }
  n_cap <- sum(sites$annual_env_mean >= anchor$cap - 1e-12, na.rm = TRUE)
  message(sprintf("  Cells at the pressure cap (%.2f): %d of %d; cells with SDM above the anchor: %d",
                  anchor$cap, n_cap, sum(!is.na(sites$sdm_suitability)),
                  sum(sites$sdm_rel > 1, na.rm = TRUE)))
  sites
}
