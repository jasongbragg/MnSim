# spatial/01_get_melaleuca_occurrences.R
#
# Download, clean, and spatially thin occurrence records of the default
# species (config/run_defaults.yaml, 'species') from the Atlas
# of Living Australia (ALA) via the galah package. The species name, range
# box and grid size come from the configuration (notes T9):
#   species/<sp>.yaml         'species' (the ALA query) and spatial$range_box
#   config/landscape.yaml     spatial$grid_deg
# cell_id is the cell number in the master grid (config/landscape.yaml,
# spatial$grid_extent; spatial/cell_ids.R; notes T13): the same for every
# species, so the shared fire tables are keyed by it. The range box only
# limits the records.
#
# OUTPUTS (<sp> = the species file name, e.g. m_nodosa)
# -------
# spatial/data/<sp>/raw.csv            raw ALA download
# spatial/data/<sp>/clean.csv          after coordinate cleaning
# spatial/data/<sp>/thinned.csv        one record per grid cell
# spatial/data/<sp>/grid.csv           occupied cells: lon_centre, lat_centre,
#                                        n_raw, cell_id -- ready for cross-
#                                        referencing with fire/rust rasters
# spatial/figures/<sp>/records.png     3-panel diagnostic map
# spatial/figures/<sp>/year_histogram.png
#
# DEPENDENCIES
# ------------
# install.packages(c("galah","CoordinateCleaner","terra","sf",
#                    "ggplot2","ozmaps","dplyr","readr"))
#
# USAGE
# -----
# Run from the project root (MnSim/plant_ibm/) or set .root below.
# Set your ALA-registered email in galah_config() before running.
# The download (~few thousand records) takes 10-30 seconds.
#
# NOTES ON GRID RESOLUTION
# -------------------------
# grid_deg = 0.1 (config/landscape.yaml) (~10km at 34°S) is the recommended starting point.
# At this resolution the fire history polygons and rust SDM have meaningful
# variation; most occupied cells will have at least 1 record post-thinning.
# Increase to 0.25 if occupied cell count is too low for a useful ensemble;
# decrease to 0.05 for finer spatial parameterisation once FRI/rust layers
# are at that resolution.

# ---- Configuration --------------------------------------------------------
source("R/species_config.R")
source("spatial/cell_ids.R")                     # master grid (notes T13)
SP          <- load_species_spatial()            # species/<sp>.yaml
ALA_EMAIL   <- "your.email@institution.edu.au"  # EDIT: your ALA-registered email
GRID_DEG    <- load_landscape_spatial()$grid_deg # grid cell size in degrees (see note above)
YEAR_MIN    <- 1970   # exclude very old records (coordinate quality declines)
MAX_COORD_UNCERTAINTY <- 5000  # metres; NA records are retained (see below)
SPECIES     <- SP$name
BOX         <- SP$range_box                      # lon_min, lon_max, lat_min, lat_max
DATA_DIR    <- SP$data_dir
FIG_DIR     <- SP$fig_dir

# ---- Setup ----------------------------------------------------------------
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR))  dir.create(FIG_DIR,  recursive = TRUE)

pkgs <- c("galah","CoordinateCleaner","terra","sf",
          "ggplot2","ozmaps","dplyr","readr")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nRun: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ","), "))")

library(galah)
library(CoordinateCleaner)
library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================================
# STEP 1: Download from ALA via galah
# ============================================================================
message("Step 1: Downloading from ALA...")

galah_config(email = "jasongbragg@gmail.com", verbose = FALSE)

occ_raw <- tryCatch(
  galah_call() |>
    galah_identify(SPECIES) |>
    galah_filter(
      basisOfRecord == "HumanObservation" |
        basisOfRecord == "PreservedSpecimen" |
        basisOfRecord == "MachineObservation",
      year >= YEAR_MIN
    ) |>
    galah_select(
      decimalLongitude, decimalLatitude,
      year, month,
      coordinateUncertaintyInMeters,
      basisOfRecord,
      datasetName,
      institutionCode,
      stateProvince
    ) |>
    atlas_occurrences(),
  error = function(e) {
    stop("galah download failed: ", conditionMessage(e),
         "\nCheck your ALA_EMAIL and internet connection.")
  }
)

message(sprintf("  Downloaded %d records.", nrow(occ_raw)))

# Standardise column names (galah versions differ slightly)
names(occ_raw) <- gsub("^decimal", "", names(occ_raw))
names(occ_raw)[names(occ_raw) == "Longitude"] <- "longitude"
names(occ_raw)[names(occ_raw) == "Latitude"]  <- "latitude"
if (!"longitude" %in% names(occ_raw))
  stop("Could not find longitude column. Column names: ",
       paste(names(occ_raw), collapse = ", "))

write_csv(occ_raw, file.path(DATA_DIR, "raw.csv"))
message("  Saved: ", file.path(DATA_DIR, "raw.csv"))

# ============================================================================
# STEP 2: Basic filters before CoordinateCleaner
# ============================================================================
message("Step 2: Basic filters...")

occ <- occ_raw |>
  filter(!is.na(longitude), !is.na(latitude)) |>
  # Coordinate uncertainty: keep NA (many herbarium records lack this field)
  # and records within threshold. Very precise records (uncertainty = 0) are
  # suspicious but rare -- keep them.
  filter(is.na(coordinateUncertaintyInMeters) |
           coordinateUncertaintyInMeters <= MAX_COORD_UNCERTAINTY) |>
  # Restrict to the species' range box (species/<sp>.yaml, spatial$range_box)
  # (prevents erroneous records in WA or overseas from passing through)
  filter(longitude >= BOX$lon_min, longitude <= BOX$lon_max,
         latitude  >= BOX$lat_min, latitude  <= BOX$lat_max) |>
  # Remove duplicate lon/lat/year combinations (re-uploaded records, etc.)
  distinct(longitude, latitude, year, .keep_all = TRUE)

message(sprintf("  After basic filters: %d records (was %d).", nrow(occ), nrow(occ_raw)))

# ============================================================================
# STEP 3: Automated coordinate cleaning with CoordinateCleaner
# ============================================================================
message("Step 3: CoordinateCleaner...")

# CoordinateCleaner needs standard column names
occ_cc <- occ |>
  rename(decimallongitude = longitude, decimallatitude = latitude)

cc_flags <- clean_coordinates(
  occ_cc |> mutate(species = SPECIES),
  lon     = "decimallongitude",
  lat     = "decimallatitude",
  species = "species",
  tests  = c("capitals",      # within 10km of a national capital (use for
                               # country-level centroid contamination)
             "centroids",     # country/province centroid coordinates
             "equal",         # lat == lon (recording error)
             "gbif",          # GBIF headquarters coordinates
             "institutions",  # within 100m of a biodiversity institution
             "seas",          # in the ocean
             "zeros"),        # exactly (0,0)
  value  = "spatialvalid",    # return the full data frame with flag columns
  verbose = FALSE
)

# .summary column is TRUE when all tests pass
occ_clean <- occ_cc |>
  filter(cc_flags$.summary) |>
  rename(longitude = decimallongitude, latitude = decimallatitude)

n_flagged <- nrow(occ) - nrow(occ_clean)
message(sprintf("  CoordinateCleaner removed %d records (%.1f%%).",
                n_flagged, 100 * n_flagged / nrow(occ)))
message(sprintf("  Clean records: %d", nrow(occ_clean)))

write_csv(occ_clean, file.path(DATA_DIR, "clean.csv"))
message("  Saved: ", file.path(DATA_DIR, "clean.csv"))

# ============================================================================
# STEP 4: Grid thinning (one record per cell)
# ============================================================================
message(sprintf("Step 4: Grid thinning at %.2f degrees (~%.0f km)...",
                GRID_DEG, GRID_DEG * 111))

# The master grid: the same cells and ids for every species (notes T13)
r <- master_raster()

# Assign each record to a cell
coords  <- cbind(occ_clean$longitude, occ_clean$latitude)
cell_id <- cellFromXY(r, coords)
occ_clean$cell_id <- cell_id

# Count raw records per cell (before thinning)
cell_counts <- occ_clean |>
  count(cell_id, name = "n_raw")

# Thin: keep the record with smallest coordinate uncertainty per cell
# (or random if uncertainty is NA throughout)
occ_thinned <- occ_clean |>
  group_by(cell_id) |>
  arrange(coordinateUncertaintyInMeters, .na_last = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

message(sprintf("  Occupied cells: %d (of %d thinned to 1 record each).",
                nrow(occ_thinned), nrow(occ_thinned)))

write_csv(occ_thinned, file.path(DATA_DIR, "thinned.csv"))
message("  Saved: ", file.path(DATA_DIR, "thinned.csv"))

# ============================================================================
# STEP 5: Build occupied-cell grid table
# ============================================================================
message("Step 5: Building occupied-cell grid table...")

# Get cell centres for each occupied cell
cell_xy <- xyFromCell(r, occ_thinned$cell_id) |>
  as.data.frame() |>
  setNames(c("lon_centre", "lat_centre"))

occ_grid <- bind_cols(
  data.frame(cell_id = occ_thinned$cell_id),
  cell_xy
) |>
  left_join(cell_counts, by = "cell_id") |>
  arrange(lat_centre, lon_centre)

message(sprintf("  Grid table: %d occupied cells, total %d raw records.",
                nrow(occ_grid), sum(occ_grid$n_raw)))

# Summary by state
state_summary <- occ_clean |>
  count(stateProvince, name = "n_records") |>
  arrange(desc(n_records))
message("  Records by state:")
print(state_summary)

write_csv(occ_grid, file.path(DATA_DIR, "grid.csv"))
message("  Saved: ", file.path(DATA_DIR, "grid.csv"))

# ============================================================================
# STEP 6: Diagnostic maps
# ============================================================================
message("Step 6: Making diagnostic maps...")

# Base layers
aus_states <- ozmaps::ozmap_states
aus_coast  <- ozmaps::ozmap_country

# Bounding box for plot
xlim <- c(BOX$lon_min + 0.5, BOX$lon_max - 0.5)
ylim <- c(BOX$lat_min + 0.5, BOX$lat_max - 0.5)

base_map <- function() {
  list(
    geom_sf(data = aus_states, inherit.aes = FALSE,
            fill = "#f5f3ee", colour = "#cccccc", linewidth = 0.3),
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE),
    theme_minimal(base_size = 10),
    theme(panel.grid = element_line(colour = "#e8e8e8", linewidth = 0.2),
          axis.text  = element_text(size = 7),
          legend.position = "bottom")
  )
}

# Panel A: raw records coloured by decade
occ_raw_plot <- occ_raw |>
  filter(!is.na(longitude), !is.na(latitude),
         longitude >= BOX$lon_min, longitude <= BOX$lon_max,
         latitude  >= BOX$lat_min, latitude  <= BOX$lat_max) |>
  mutate(decade = paste0(floor(year / 10) * 10, "s"))

pA <- ggplot(occ_raw_plot, aes(x = longitude, y = latitude, colour = decade)) +
  base_map() +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_viridis_d(name = "Decade", option = "C") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = sprintf("Raw ALA records (n = %d)", nrow(occ_raw_plot)),
       x = NULL, y = NULL)

# Panel B: clean records
pB <- ggplot(occ_clean, aes(x = longitude, y = latitude)) +
  base_map() +
  geom_point(colour = "#2a78d6", size = 0.7, alpha = 0.6) +
  labs(title = sprintf("After coordinate cleaning (n = %d)", nrow(occ_clean)),
       x = NULL, y = NULL)

# Panel C: thinned grid cells, sized by raw record count
pC <- ggplot(occ_grid, aes(x = lon_centre, y = lat_centre, size = n_raw)) +
  base_map() +
  geom_point(colour = "#e34948", alpha = 0.8) +
  scale_size_continuous(name = "Raw records\nper cell",
                        range = c(0.8, 4),
                        breaks = c(1, 5, 20, 50)) +
  labs(title = sprintf("Grid thinned (%.2f°): %d occupied cells",
                       GRID_DEG, nrow(occ_grid)),
       x = NULL, y = NULL)

# Combine and save
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(pA, pB, pC, ncol = 3) +
    patchwork::plot_annotation(
      title    = paste(SPECIES, "— ALA occurrence records"),
      subtitle = sprintf("Downloaded %s | Year >= %d | Grid = %.2f° | Uncertainty <= %d m",
                          Sys.Date(), YEAR_MIN, GRID_DEG, MAX_COORD_UNCERTAINTY),
      theme = theme(plot.title = element_text(size = 12, face = "bold"),
                    plot.subtitle = element_text(size = 8)))
} else {
  fig <- pA  # fallback if patchwork not installed
}

ggsave(file.path(FIG_DIR, "records.png"), fig,
       width = 14, height = 5.5, dpi = 150, bg = "white")
message("  Saved: ", file.path(FIG_DIR, "records.png"))

# ---- Supplementary: year histogram ----------------------------------------
year_hist <- occ_clean |>
  filter(!is.na(year)) |>
  ggplot(aes(x = year)) +
  geom_histogram(binwidth = 5, fill = "#2a78d6", colour = "white") +
  geom_vline(xintercept = YEAR_MIN, linetype = "dashed", colour = "#888780") +
  labs(x = "Year", y = "Records",
       title = "Temporal coverage of cleaned records") +
  theme_minimal(base_size = 10)

ggsave(file.path(FIG_DIR, "year_histogram.png"), year_hist,
       width = 7, height = 3.5, dpi = 150, bg = "white")

# ============================================================================
# SUMMARY
# ============================================================================
cat(sprintf("\n=== %s occurrence summary ===\n", SPECIES))
cat(sprintf("  Raw download:        %5d records\n", nrow(occ_raw)))
cat(sprintf("  After basic filter:  %5d records\n", nrow(occ)))
cat(sprintf("  After CC cleaning:   %5d records\n", nrow(occ_clean)))
cat(sprintf("  Grid-thinned cells:  %5d occupied cells (%.2f° grid)\n",
            nrow(occ_grid), GRID_DEG))
cat(sprintf("  Total raw in cells:  %5d records\n", sum(occ_grid$n_raw)))
cat(sprintf("  Records per cell:    mean=%.1f  median=%.0f  max=%.0f\n",
            mean(occ_grid$n_raw), median(occ_grid$n_raw), max(occ_grid$n_raw)))
cat(sprintf("\nOutputs in %s/ and %s/\n", DATA_DIR, FIG_DIR))
cat("Next step: spatial/02_cross_reference_layers.R\n")
