# spatial/02c_update_rust_pressure.R
#
# Re-run ONLY the rust step of 02_cross_reference_layers.R: refresh the SDM
# and rust-pressure columns of spatial/data/mn_site_params.csv, leaving the
# fire columns and weights untouched. Use after an SDM update.
#
#   cd plant_ibm
#   MN_SDM_TIFF=/path/to/AP_suitability_current.tif \
#     /opt/R/4.3.1/bin/Rscript spatial/02c_update_rust_pressure.R
#
# Without MN_SDM_TIFF (or if the file is missing) the existing
# sdm_suitability column is re-normalised, which only refreshes the derived
# columns (useful after changing RUST_ANCHOR in params.R).

source("params.R")
source("spatial/rust_pressure.R")

SITE_CSV <- "spatial/data/mn_site_params.csv"
SDM_TIFF <- Sys.getenv("MN_SDM_TIFF",
  "/home/jgb/projects/mrsdm/output/projections/AP_suitability_current.tif")

sites <- read.csv(SITE_CSV)
sdm <- if (file.exists(SDM_TIFF)) {
  message("Re-extracting SDM from ", SDM_TIFF)
  terra::rast(SDM_TIFF)
} else {
  message("No SDM raster at ", SDM_TIFF, "; re-normalising the stored sdm_suitability")
  NULL
}
sites <- apply_sdm_to_sites(sites, sdm, RUST_ANCHOR)
write.csv(sites, SITE_CSV, row.names = FALSE)
message("Updated ", SITE_CSV)
