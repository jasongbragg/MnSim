# spatial/fire/fire_config.R -- settings for the fire pipeline (NSW + Qld).
# Sourced by spatial/fire/*.R scripts, spatial/02_cross_reference_layers.R
# and the tests. Paths are relative to plant_ibm/ (run everything from there);
# the two external sources can be overridden with environment variables.
# Methods and the reasons for each choice: docs/methodological_notes.md, D1-D5.

# ---- Paths ---------------------------------------------------------------------
NPWS_SHP      <- Sys.getenv("NPWS_SHP", "/home/jgb/projects/melsims/utilities/fire/NPWSFireHistory.shp")
MN_GRID_CSV   <- "spatial/data/mn_grid.csv"          # occupied 0.1-degree cells (01)
FIRE_DATA_DIR <- "spatial/data/fire"                 # cell-year tables (kept in git)
FIRE_RAW_DIR  <- Sys.getenv("FIRE_RAW_DIR", file.path(FIRE_DATA_DIR, "raw"))  # Qld downloads (not in git)
FIRE_FIG_DIR  <- "spatial/figures/fire"
GRID_DEG      <- 0.1

# ---- Queensland services (ArcGIS REST) -----------------------------------------
QLD_SCARS_URL  <- "https://spatial-gis.information.qld.gov.au/arcgis/rest/services/Environment/FireScarMapping/MapServer"
QLD_QPWS_URL   <- "https://spatial-gis.information.qld.gov.au/arcgis/rest/services/Environment/ParksFireHistory/MapServer/0"
QLD_SCAR_YEARS <- 1987:2025      # layer id = year - 1986 (all downloaded; see QLD_YEARS for what is used)
QLD_BBOX       <- c(xmin = 147.9, ymin = -29.3, xmax = 154.2, ymax = -23.9)  # Qld cells + margin

AREA_CRS <- 3577                 # Australian Albers: all areas computed in this CRS

# ---- Observation periods (D1, D4) ------------------------------------------------
NSW_YEARS <- 1987:2024           # NPWS record (1970-2024 kept as a sensitivity check)
QLD_YEARS <- 2014:2024           # edited scar mapping only: 1987-2013 rejected (D4)
NSW_YEARS_ALL <- 1970:2024       # years processed into the NSW cell-year tables
QLD_YEARS_ALL <- 1987:2024       # years processed into the Qld cell-year tables

# ---- Fire-year definition and estimator (D2) ---------------------------------------
AREA_THRESH <- 0.10     # a cell-year is a fire year when >= 10% of the cell's land burns
PRIOR_ALPHA <- 2        # fixed Gamma prior, used only if the empirical-Bayes fit fails
PRIOR_BETA  <- 100
FIRE_MIN    <- 0.005    # floor on the annual rate (wildfire, all fire); 0 for planned burns
FIRE_MAX    <- 0.500
PFIMP_MIN   <- 0.05

# ---- Fire regimes ------------------------------------------------------------------
# 'fire': which fires count (D5)
#   "wild"    NSW wildfires; Qld scars minus QPWS planned burns   <- the default
#   "all"     every recorded fire (NSW wildfire + prescribed + unlabelled; Qld all scars)
#   "planned" NSW prescribed burns; QPWS planned burns
# 'prior': "nb_ml" (negative-binomial maximum likelihood, D2) or "mm02"
#   (method of moments, the pre-Nov-2026 02); 'prior_by': "pooled" or "state".
fire_scenario <- function(nsw = NSW_YEARS, qld = QLD_YEARS, fire = "wild", area_base = "land",
                          prior = "nb_ml", prior_by = "pooled", fire_min = FIRE_MIN)
  list(nsw = nsw, qld = qld, fire = fire, area_base = area_base, prior = prior,
       prior_by = prior_by, fire_min = fire_min)

# Written to mn_site_params.csv by 02: 'wild' under the standard column names,
# the others with a suffix (_all, _planned).
FIRE_REGIMES <- list(
  wild    = fire_scenario(),
  all     = fire_scenario(fire = "all"),
  planned = fire_scenario(fire = "planned", prior_by = "state", fire_min = 0)
)

# Sensitivity checks (spatial/fire/03_fire_checks.R); all wildfire unless stated.
FIRE_SENSITIVITY <- list(
  prior_state   = fire_scenario(prior_by = "state"),
  common_period = fire_scenario(nsw = QLD_YEARS),
  nsw_long      = fire_scenario(nsw = 1970:2024),
  qld_from2013  = fire_scenario(qld = 2013:2024),
  qld_automated = fire_scenario(qld = 1987:2024, prior_by = "state")  # rejected record
)

# Cell-year table for a state, fire type and area base (written by 02_fire_cell_years.R).
cy_table <- function(state, fire, area_base = "land") {
  if (area_base != "land") stop("only the land base is tabulated")
  if (state == "NSW") c(all = "nsw_all", wild = "nsw_wild", planned = "nsw_planned")[[fire]]
  else                c(all = "qld_all", wild = "qld_wild", planned = "qld_planned")[[fire]]
}

# ---- Border check (D4) -------------------------------------------------------------
QLD_BORDER_LAT <- -28.16         # used only if ozmaps is unavailable
BORDER_BAND    <- 1.0            # degrees of latitude either side of the border
