# spatial/fire/01_download_qld.R -- download Queensland fire polygons.
# Run from plant_ibm/. Needs internet; resumable (finished years are skipped,
# an interrupted year continues from its last batch). Writes to FIRE_RAW_DIR
# (spatial/data/fire/raw/, not kept in git):
#   qld_scars/qld_scars_<year>.gpkg   QRSC fire scars 1987-2025
#   qpws_fire_history.gpkg            QPWS planned burns and wildfires
# each with '<file>.bbox', the area downloaded: qld_bbox(), around the cells of
# every species that may need Qld fire data (spatial/fire/R/bbox.R; notes T10).
# A download that does not cover the cells (e.g. after adding a more northern
# species) is fetched again.
source("spatial/fire/fire_config.R")
source("spatial/fire/R/bbox.R")
source("spatial/fire/R/qld_download.R")
message("Qld download area: ", paste(sprintf("%.2f", qld_bbox()), collapse = ", "), " (xmin, ymin, xmax, ymax)")
download_qld_scars()
download_qpws()
