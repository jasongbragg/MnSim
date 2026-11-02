# spatial/rust_pressure.R
#
# Maps a myrtle rust SDM onto the pathogen vertex of the disease triangle
# (annual_env Beta parameters) for every occupied grid cell, anchored at one
# site (RUST_ANCHOR in params.R).
#
# ASSUMPTION (docs/methodological_notes.md, "Rust pressure"):
#   mean annual pressure in cell i = mu * s_i / s_anchor, capped at 'cap'
# where s is the SDM value extracted the same way for every cell (at the
# 0.1-degree cell centre), and the anchor cell is the occupied cell nearest
# the anchor site. Dividing by the anchor cell's own SDM value means:
#   - the anchor cell always gets exactly mu (0.70), for any SDM, so the
#     field calibration at Yarringully and the Pegg et al. validation survive
#     SDM updates;
#   - a new SDM changes only the relative pattern across the range.
# The normaliser is computed from the data each time; it is never typed in.
#
# Functions only; sourced by 02_cross_reference_layers.R,
# 02c_update_rust_pressure.R and the tests.

# Index (row) of the occupied cell nearest the anchor site.
anchor_cell_index <- function(lon, lat, anchor = RUST_ANCHOR) {
  d2 <- (lon - anchor$lon)^2 + ((lat - anchor$lat))^2
  which.min(d2)
}

# Core mapping, on SDM values already extracted per cell.
rust_pressure_from_sdm <- function(sdm, anchor_index, anchor = RUST_ANCHOR) {
  s_anchor <- sdm[anchor_index]
  if (is.na(s_anchor) || s_anchor <= 0)
    stop("SDM value at the anchor cell is missing or <= 0; cannot normalise.")
  floor_env <- 1 / (anchor$conc * 10)    # numerical floor: alpha >= 0.1
  rel  <- sdm / s_anchor
  mean <- pmin(anchor$cap, pmax(floor_env, anchor$mu * rel))
  data.frame(
    sdm_rel          = rel,
    annual_env_mean  = mean,
    annual_env_alpha = round(mean * anchor$conc, 4L),
    annual_env_beta  = round((1 - mean) * anchor$conc, 4L)
  )
}

# Add / refresh the rust columns of a site table (one row per occupied cell,
# columns lon_centre, lat_centre). If 'sdm' (a terra SpatRaster) is given the
# SDM is re-extracted at cell centres; otherwise the existing
# sdm_suitability column is used. Prints the anchor diagnostics.
apply_sdm_to_sites <- function(sites, sdm = NULL, anchor = RUST_ANCHOR,
                               neighbourhood_deg = 0.15, warn_rel_diff = 0.10) {
  if (!is.null(sdm)) {
    pts  <- terra::vect(as.data.frame(sites[, c("lon_centre", "lat_centre")]),
                        geom = c("lon_centre", "lat_centre"), crs = "EPSG:4326")
    vals <- terra::extract(sdm, pts)[, 2L]
    vals[!is.na(vals) & vals < -1e30] <- NA_real_
    if (max(vals, na.rm = TRUE) > 1.5) vals <- vals / 1000   # biomod2 0-1000 scale
    sites$sdm_suitability <- vals
  }
  k <- anchor_cell_index(sites$lon_centre, sites$lat_centre, anchor)
  rp <- rust_pressure_from_sdm(sites$sdm_suitability, k, anchor)
  sites[names(rp)] <- rp
  sites$rust_anchor_cell <- seq_len(nrow(sites)) == k

  message(sprintf(
    "  Rust anchor: %s -> cell %s (%.2f, %.2f), SDM %.4f -> pressure %.2f",
    anchor$site, if ("cell_id" %in% names(sites)) sites$cell_id[k] else k,
    sites$lon_centre[k], sites$lat_centre[k], sites$sdm_suitability[k], anchor$mu))

  # Is the anchor cell typical of its surroundings? (single-cell noise guard)
  if (!is.null(sdm)) {
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
