# spatial/fire/01_download_qld.R -- download Queensland fire polygons.
# Run from plant_ibm/. Needs internet; resumable (finished years are skipped,
# an interrupted year continues from its last batch). Writes to FIRE_RAW_DIR
# (spatial/data/fire/raw/, not kept in git):
#   qld_scars/qld_scars_<year>.gpkg   QRSC fire scars 1987-2025, M. nodosa area
#   qpws_fire_history.gpkg            QPWS planned burns and wildfires
source("spatial/fire/fire_config.R")
source("spatial/fire/R/qld_download.R")
download_qld_scars()
download_qpws()
