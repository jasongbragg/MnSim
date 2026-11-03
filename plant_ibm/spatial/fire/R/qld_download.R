# spatial/fire/R/qld_download.R -- fetch Queensland fire polygons for the M. nodosa area.
#
# Both services are ArcGIS REST MapServer layers. Large years hold >100,000
# scar polygons in our bounding box, and deep offset paging fails with HTTP
# 500, so downloads work by object id instead:
#   1. one request returns every matching object id (returnIdsOnly; not
#      limited by maxRecordCount);
#   2. features are fetched in batches of ids (POST), each batch saved as it
#      arrives, so an interrupted year resumes where it stopped;
#   3. failed requests are retried with back-off; a batch that keeps failing
#      is split in half until it goes through;
#   4. outlines are simplified by the server to ~5 m (maxAllowableOffset;
#      half a Sentinel-2 pixel), which shrinks downloads a lot and changes
#      areas negligibly.
# Output: data/raw/qld_scars/qld_scars_<year>.gpkg (one per year) and
#         data/raw/qpws_fire_history.gpkg. Finished files are skipped on rerun.

suppressPackageStartupMessages({ library(sf); library(httr); library(jsonlite) })

ID_BATCH        <- 500L        # features per request
MAX_OFFSET_DEG  <- 0.00005     # ~5 m generalisation (degrees, outSR 4326)
RETRIES         <- 5L

.bbox_params <- function(bbox) list(
  geometry = paste(bbox[c("xmin", "ymin", "xmax", "ymax")], collapse = ","),
  geometryType = "esriGeometryEnvelope", inSR = 4326,
  spatialRel = "esriSpatialRelIntersects")

# POST with retries and exponential back-off; returns the response text.
.post_retry <- function(url, body, retries = RETRIES) {
  for (i in seq_len(retries)) {
    resp <- tryCatch(httr::POST(url, body = body, encode = "form", httr::timeout(600)),
                     error = function(e) e)
    if (!inherits(resp, "error") && httr::status_code(resp) == 200) {
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      if (!grepl('^\\s*\\{\\s*"error"', txt)) return(txt)
      msg <- tryCatch(jsonlite::fromJSON(txt)$error$message, error = function(e) "service error")
    } else {
      msg <- if (inherits(resp, "error")) conditionMessage(resp) else paste("HTTP", httr::status_code(resp))
    }
    wait <- 2^i
    message(sprintf("      attempt %d failed (%s); retrying in %d s", i, msg, wait))
    Sys.sleep(wait)
  }
  stop("request failed after ", retries, " attempts")
}

arcgis_object_ids <- function(layer_url, bbox, where = "1=1") {
  txt <- .post_retry(paste0(layer_url, "/query"),
                     c(list(where = where, returnIdsOnly = "true", f = "json"), .bbox_params(bbox)))
  ids <- jsonlite::fromJSON(txt)$objectIds
  if (is.null(ids)) integer(0) else sort(as.integer(ids))
}

# Fetch features for a set of ids; on persistent failure split the set.
arcgis_fetch_ids <- function(layer_url, ids, out_fields, max_offset = MAX_OFFSET_DEG) {
  body <- list(objectIds = paste(ids, collapse = ","), outFields = out_fields,
               returnGeometry = "true", outSR = 4326, f = "geojson",
               maxAllowableOffset = max_offset, geometryPrecision = 6)
  res <- tryCatch(.post_retry(paste0(layer_url, "/query"), body), error = function(e) NULL)
  if (is.null(res)) {
    if (length(ids) <= 10) stop("could not fetch ids ", paste(range(ids), collapse = "-"))
    h <- length(ids) %/% 2
    message(sprintf("      splitting a batch of %d ids", length(ids)))
    return(rbind(arcgis_fetch_ids(layer_url, ids[seq_len(h)], out_fields, max_offset),
                 arcgis_fetch_ids(layer_url, ids[-seq_len(h)], out_fields, max_offset)))
  }
  tf <- tempfile(fileext = ".geojson"); on.exit(unlink(tf))
  writeLines(res, tf)
  x <- sf::read_sf(tf)
  if (nrow(x) == 0) return(NULL)
  x
}

# Download one layer within the bbox to 'out', resumable via a partial folder.
download_layer <- function(layer_url, bbox, out, out_fields, batch = ID_BATCH, label = "") {
  if (file.exists(out)) { message(label, ": already downloaded"); return(invisible(TRUE)) }
  part_dir <- paste0(out, ".partial")
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  ids_file <- file.path(part_dir, "ids.rds")
  ids <- if (file.exists(ids_file)) readRDS(ids_file) else {
    i <- arcgis_object_ids(layer_url, bbox); saveRDS(i, ids_file); i }
  n_b <- ceiling(length(ids) / batch)
  message(sprintf("%s: %d features in %d batches", label, length(ids), n_b))
  for (b in seq_len(n_b)) {
    f <- file.path(part_dir, sprintf("batch_%05d.rds", b))
    if (file.exists(f)) next
    sel <- ids[((b - 1) * batch + 1):min(b * batch, length(ids))]
    saveRDS(arcgis_fetch_ids(layer_url, sel, out_fields), f)
    if (b %% 20 == 0 || b == n_b) message(sprintf("    %s: batch %d / %d", label, b, n_b))
  }
  parts <- lapply(list.files(part_dir, "^batch_.*\\.rds$", full.names = TRUE), readRDS)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  x <- if (length(parts)) do.call(rbind, parts) else NULL
  if (is.null(x)) x <- st_sf(year = integer(0), geometry = st_sfc(crs = 4326))
  st_write(x, out, quiet = TRUE)
  unlink(part_dir, recursive = TRUE)
  message(sprintf("%s: %d polygons saved", label, nrow(x)))
  invisible(TRUE)
}

# Fire scars, one file per year (layer id = year - 1986).
download_qld_scars <- function(years = QLD_SCAR_YEARS, bbox = QLD_BBOX,
                               out_dir = file.path(FIRE_RAW_DIR, "qld_scars")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (y in years)
    download_layer(paste0(QLD_SCARS_URL, "/", y - 1986L), bbox,
                   file.path(out_dir, sprintf("qld_scars_%d.gpkg", y)),
                   out_fields = "objectid,year,month,area_ha", label = as.character(y))
  invisible(TRUE)
}

# QPWS fire history (planned burns PB and wildfires WF on QPWS / state land).
download_qpws <- function(bbox = QLD_BBOX, out = file.path(FIRE_RAW_DIR, "qpws_fire_history.gpkg")) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  download_layer(QLD_QPWS_URL, bbox, out,
                 out_fields = "objectid,label,type,burnstatus,ignitiondate,outdate,outyear,completedyear,month,area_ha",
                 label = "QPWS")
}
