# spatial/02_cross_reference_layers.R
#
# Cross-reference 375 M. nodosa occupied cells with:
#   (A) NSW NPWS Fire History → fire_prob_annual and fire_p_fimp per cell
#   (B) Myrtle rust MaxEnt SDM → annual_env Beta parameters per cell
#
# FIRE METHODOLOGY
# ----------------
# Current implementation: WILDFIRE ONLY, area-weighted.
#
# Two parameters are estimated separately from the fire history:
#
#   fire_prob_annual = fraction of years 1970-2024 with any wildfire in cell
#                    = n_wildfire_years / OBS_PERIOD
#     Answers: HOW OFTEN does a fire event occur?
#
#   fire_p_fimp      = mean fractional area of cell burned per wildfire year
#                    = mean(area_intersection / area_cell) over fire years only
#     Answers: WHEN A FIRE OCCURS, what fraction of the population is affected?
#
# This separation is biologically meaningful: a cell may burn frequently but
# with small, patchy fires (low fire_p_fimp), or rarely but with landscape-
# scale fires that cover the whole population (fire_p_fimp → 1). Both
# parameters feed directly into the IBM.
#
# The variance of area fraction across fire years is also recorded as
# var_frac_burned. This could later be used to parameterise stochastic
# variation in fire_p_fimp across fire events (e.g., as a Beta distribution),
# but is currently used only as a diagnostic metric.
#
# PRESCRIBED BURNS: currently excluded. To reactivate, see the clearly
# marked block in Step 2 below. Note that including prescribed burns would
# require either (a) separate fire_p_fimp values for wildfire vs prescribed
# (biologically justified: prescribed burns are lower intensity and patchier),
# or (b) accepting that the area-weighted approach naturally discounts small
# burns (most prescribed burns contribute small fractions, reducing mean
# fire_p_fimp). The IBM currently has no mechanism to run two fire types per
# year, so option (b) is simpler when prescribed burns are reactivated.
#
# SDM CALIBRATION
# ---------------
# Yarringully NR anchor: lon=153.21967, lat=-29.09357, SDM=0.5736
# annual_env_mean(site) = MU_YARR × (SDM_site / SDM_YARR)
# Beta concentration (alpha+beta = CONC) held constant across sites.
# At Yarringully: Beta(14,6), mean=0.70, SD=0.10 — consistent with params.R
#
# OUTPUT
# ------
# spatial/data/mn_site_params.csv   one row per occupied cell, all IBM params
# spatial/figures/mn_site_params.png two-panel diagnostic map
#
# USAGE
# -----
# Edit the four paths in Configuration, then:
#   source("spatial/02_cross_reference_layers.R")

# ============================================================================
# Configuration — edit paths before running
# ============================================================================
FIRE_SHP    <- "/home/jgb/projects/melsims/utilities/fire/NPWSFireHistory.shp"  # edit to local path after download
#SDM_TIFF    <- "/home/jgb/projects/melsims/utilities/rust/maxent_MR.tiff"       # edit to local path
SDM_TIFF    <- "/home/jgb/projects/mrsdm/output/projections/AP_suitability_current.tif"
MN_GRID   <- "spatial/data/mn_grid.csv"
GRID_DEG  <- 0.1          # must match 01_get_melaleuca_occurrences.R

OBS_START <- 1970L
OBS_END   <- 2024L
OBS_PERIOD <- OBS_END - OBS_START + 1L    # 55 years

# Fire probability floor: 0.020 = 1 fire per 50 years.
# Literature basis: natural wildfire FRI in NSW/QLD coastal heathland and
# wallum scrub (the vegetation types where M. nodosa occurs) is typically
# 7-30 years (Hotspots Fire Project; Keith & Bradstock 1994). Long-unburnt
# heathland at 20-40 years is ecologically notable (Stanton et al. 2025
# Austral Ecology; Freestone et al. 2015 Aust J Bot). At 37+ years without
# fire, heathland begins transitioning to woodland. A floor of 50 years
# (fire_prob_annual = 0.020) represents the upper end of the plausible
# natural fire regime, not fire-absence. A 100-year floor would imply
# effective fire absence over human-relevant timescales, inconsistent with
# the ecology of this vegetation type.
FIRE_MAX      <- 0.500    # hard ceiling fire_prob_annual (~1 fire per 2 yr)
FIRE_MIN      <- 0.005    # hard floor (biological minimum; Bayesian mean handles the rest)
FIRE_UNBURNED <- 0.020    # retained as a comment reference; no longer used as a floor

# Bayesian Gamma-Poisson prior for fire rate (lambda = fire_prob_annual).
# Fires arrive as a Poisson process; conjugate prior is Gamma(alpha, beta).
# Posterior given k fires in T years: Gamma(alpha+k, beta+T)
# Posterior mean = (alpha + k) / (beta + T)
#
# Interpretation: PRIOR_ALPHA pseudo-fires in PRIOR_BETA pseudo-years of
# prior observation — equivalent to having "seen" 2 fires in 100 years.
# Prior mean = PRIOR_ALPHA / PRIOR_BETA = 0.020 (FRI = 50 yr).
#
# This regularises 0-fire cells toward FRI~78 yr rather than infinity,
# and multi-fire cells toward the prior from above. No arbitrary floor needed.
#
# These values are fixed priors; see Step 3b for the optional Empirical Bayes
# update that sets them from the observed data distribution across all cells.
PRIOR_ALPHA   <- 2        # pseudo-fires (prior strength)
PRIOR_BETA    <- 100      # pseudo-years (prior strength)
PFIMP_MIN     <- 0.05     # floor fire_p_fimp (minimum plausible impact fraction)
# PFIMP_DEFAULT (for cells with no significant wildfire) is set post-hoc
# to the range-wide median mean_frac_burned from cells that DID have fires.
# This avoids the inversion where unburned cells get higher p_fimp than
# cells with many tiny fires.

# Minimum fire size threshold: only fire years where the total burned fraction
# in the cell exceeds AREA_THRESH are counted toward FRI and fire_p_fimp.
# This filters out small fires that affect a negligible fraction of the
# M. nodosa population in the cell.
# AREA_THRESH = 0.05  → 5% of cell ≈ 500 ha  (moderate, recommended)
# AREA_THRESH = 0.10  → 10% of cell ≈ 1000 ha (conservative)
# AREA_THRESH = 0.02  → 2% of cell ≈ 200 ha   (lenient)
# All wildfire areas (below and above threshold) are still recorded as
# n_wildfire_years_all for sensitivity checking and method transparency.
AREA_THRESH <- 0.1

YNR_LON   <- 153.21967
YNR_LAT   <- -29.09357
# SDM_YARR  <- 0.5736       # MaxEnt at Yarringully (do not change)
SDM_YARR <- 0.993
MU_YARR   <- 0.70         # calibrated annual_env mean at Yarringully
CONC      <- 20L          # Beta concentration (alpha + beta), fixed

# ============================================================================
# Setup
# ============================================================================
if (!dir.exists("spatial/data"))    dir.create("spatial/data",    recursive = TRUE)
if (!dir.exists("spatial/figures")) dir.create("spatial/figures", recursive = TRUE)

pkgs <- c("sf","terra","dplyr","readr","ggplot2","ozmaps")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing: ", paste(missing_pkgs, collapse = ", "))

library(sf); library(terra); library(dplyr); library(readr); library(ggplot2)

# ============================================================================
# STEP 1: Load M. nodosa grid, build cell polygons
# ============================================================================
message("Step 1: Loading M. nodosa grid...")
mn_grid <- read_csv(MN_GRID, show_col_types = FALSE)
message(sprintf("  %d occupied cells at %.2f deg", nrow(mn_grid), GRID_DEG))

half <- GRID_DEG / 2
cell_polys <- lapply(seq_len(nrow(mn_grid)), function(i) {
  x <- mn_grid$lon_centre[i]; y <- mn_grid$lat_centre[i]
  st_polygon(list(matrix(c(x-half,y-half, x+half,y-half,
                             x+half,y+half, x-half,y+half,
                             x-half,y-half), ncol=2L, byrow=TRUE)))
})
mn_sf <- st_sf(mn_grid, geometry = st_sfc(cell_polys, crs = 4326L))

# ============================================================================
# STEP 2: Load and prepare fire history
# ============================================================================
message("Step 2: Loading fire history shapefile...")

# Disable S2: government shapefiles commonly have minor geometric invalidity
# (duplicate vertices, etc.) that fail S2 but work fine under GEOS
sf_use_s2(FALSE)
options(sf.use_s2 = FALSE)   # suppress planar coordinate warnings

fire_raw <- st_read(FIRE_SHP, quiet = TRUE)
message(sprintf("  %d polygons, CRS: EPSG:%s", nrow(fire_raw),
                st_crs(fire_raw)$epsg))

# Repair invalid geometries before any spatial operations
fire_raw <- st_make_valid(fire_raw)

# Reproject GDA94 -> WGS84 (EPSG:4283 -> 4326)
fire <- st_transform(fire_raw, 4326L)

# Extract fire year from Label: "YYYY-YY Wildfire/Prescribed Burn"
fire$fire_year <- as.integer(substr(fire$Label, 1L, 4L))
fire$fire_type <- dplyr::case_when(
  grepl("Wildfire",   fire$Label, ignore.case = TRUE) ~ "Wildfire",
  grepl("Prescribed", fire$Label, ignore.case = TRUE) ~ "Prescribed",
  TRUE ~ "Unknown"
)

# Filter to observation period
fire_obs <- fire[!is.na(fire$fire_year) &
                   fire$fire_year >= OBS_START &
                   fire$fire_year <= OBS_END, ]

# ---- WILDFIRE FILTER -------------------------------------------------------
# Currently using WILDFIRES ONLY.
# Rationale: the IBM fire event represents a landscape-scale disturbance
# calibrated against the 2020 Black Summer wildfire (Pegg et al. validation).
# Prescribed burns are typically patchy and low-intensity; including them
# would inflate fire frequency and conflate two ecologically different regimes.
#
# TO REACTIVATE PRESCRIBED BURNS:
#   1. Comment out the next line and uncomment the alternative below
#   2. Decide whether to run a single combined fire_p_fimp (the area-weighted
#      mean naturally discounts small prescribed burns) or add a separate
#      IBM mechanism for low-intensity fire (requires new model code)
fire_active <- fire_obs[fire_obs$fire_type == "Wildfire", ]
# fire_active <- fire_obs  # <-- reactivate this to include prescribed burns
# ---- END WILDFIRE FILTER ---------------------------------------------------

message(sprintf("  Wildfire polygons in period: %d (prescribed excluded: %d)",
                nrow(fire_active),
                sum(fire_obs$fire_type == "Prescribed")))

# Clip to M. nodosa bounding box before intersection (major speed-up)
bbox_mn <- st_as_sfc(st_bbox(mn_sf))
fire_clip <- fire_active[st_intersects(fire_active, bbox_mn, sparse=FALSE)[,1L], ]
message(sprintf("  After M. nodosa bbox clip: %d polygons", nrow(fire_clip)))

# ============================================================================
# STEP 3: Area-weighted intersection — fire polygons x grid cells
# ============================================================================
message("Step 3: Computing area-weighted fire statistics per cell...")
message("  (st_intersection on ~375 cells; expect 2-5 minutes)")

# Cell areas in m² (accounts for latitude-dependent cell size)
cell_areas_m2 <- as.numeric(st_area(mn_sf))

# Fast initial intersection test to get candidate pairs
idx_list <- st_intersects(mn_sf, fire_clip)

fire_stats <- lapply(seq_len(nrow(mn_sf)), function(i) {

  if (i %% 75L == 0L)
    message(sprintf("    Cell %d / %d  (%.0f%%)", i, nrow(mn_sf),
                    100 * i / nrow(mn_sf)))

  idx <- idx_list[[i]]
  no_data <- list(n_wildfire_years_all = 0L, n_wildfire_years_sig = 0L,
                  n_wildfire_events    = 0L,
                  # Bayesian posterior for 0 fires: Gamma(alpha, beta+T) mean
                  fire_prob_annual     = (PRIOR_ALPHA) / (PRIOR_BETA + OBS_PERIOD),
                  mean_frac_burned     = NA_real_,
                  sd_frac_burned       = NA_real_,
                  var_frac_burned      = NA_real_,
                  fire_p_fimp          = NA_real_)

  if (length(idx) == 0L) return(no_data)

  # Subset to candidate fire polygons, keep only needed columns
  fire_sub <- fire_clip[idx, c("fire_year", "fire_type", "geometry")]
  cell_i   <- mn_sf[i, "geometry"]
  cell_a   <- cell_areas_m2[i]

  # Intersection geometry: gives the actual burned area within this cell
  inter <- tryCatch(
    suppressWarnings(st_intersection(cell_i, fire_sub)),
    error = function(e) NULL
  )
  if (is.null(inter) || nrow(inter) == 0L) return(no_data)

  # Extract polygon parts only — but only if geometry is a collection type;
  # calling st_collection_extract on a plain POLYGON raises a warning
  geom_types <- unique(as.character(st_geometry_type(inter)))
  if (any(grepl("COLLECTION|GEOMETRYCOLLECTION|MULTI", geom_types))) {
    inter <- tryCatch(
      st_collection_extract(inter, "POLYGON"),
      error = function(e) inter
    )
  }
  if (nrow(inter) == 0L) return(no_data)

  # Fraction of cell area burned by each intersecting fire polygon
  inter_areas <- as.numeric(st_area(inter))
  frac        <- pmin(1, inter_areas / cell_a)

  # Aggregate by year: sum all fire fractions within a year, cap at 1
  by_year      <- tapply(frac, inter$fire_year,
                          function(x) min(1, sum(x, na.rm = TRUE)))
  annual_frac_all <- as.numeric(by_year)    # all wildfire years
  n_yrs_all       <- length(annual_frac_all)

  # Apply area threshold: only count years where a SUBSTANTIAL fraction
  # of the cell burned. Small fires affecting <AREA_THRESH of the cell
  # are unlikely to affect a meaningful fraction of the M. nodosa population
  # and should not inflate fire frequency or drive down FRI.
  # n_wildfire_years_all is retained for transparency and sensitivity checks.
  annual_frac_sig <- annual_frac_all[annual_frac_all >= AREA_THRESH]
  n_yrs_sig       <- length(annual_frac_sig)

  # Use SIGNIFICANT fire years for FRI and p_fimp
  var_f <- if (n_yrs_sig > 1L) var(annual_frac_sig) else NA_real_
  sd_f  <- if (n_yrs_sig > 1L) sd(annual_frac_sig)  else NA_real_

  list(
    n_wildfire_years_all = n_yrs_all,           # all wildfire years (for sensitivity)
    n_wildfire_years_sig = n_yrs_sig,           # significant-fire years (for FRI)
    n_wildfire_events    = length(frac),
    # Bayesian posterior mean: Gamma(PRIOR_ALPHA + k_sig, PRIOR_BETA + OBS_PERIOD)
    # Regularises 0-fire cells toward FRI ~ PRIOR_BETA/PRIOR_ALPHA years.
    # No arbitrary floor needed: the prior handles the sparse-data case.
    fire_prob_annual  = pmin(FIRE_MAX, pmax(FIRE_MIN,
                             (PRIOR_ALPHA + n_yrs_sig) / (PRIOR_BETA + OBS_PERIOD))),
    mean_frac_burned  = if (n_yrs_sig > 0L) mean(annual_frac_sig) else NA_real_,
    sd_frac_burned    = sd_f,
    var_frac_burned   = var_f,
    fire_p_fimp       = NA_real_    # filled post-hoc from significant fires
  )
})

fire_df <- do.call(rbind, lapply(fire_stats, as.data.frame))
mn_grid <- bind_cols(mn_grid, fire_df)

# ---- Optional Empirical Bayes: update prior from observed data ---------------
# Rather than fixing PRIOR_ALPHA=2, PRIOR_BETA=100 from theory, fit the prior
# from the distribution of observed k_i (significant fire years) across cells.
# This sets alpha and beta to match the empirical mean and variance of k/T,
# using the method-of-moments for the Gamma distribution:
#   mean_obs = alpha/beta  →  alpha = mean^2/var,  beta = mean/var
# The posterior mean then automatically reflects the range-wide fire regime
# rather than an externally assumed FRI=50 baseline.
#
# NOTE: with only 375 cells and many zeros, this estimate may be noisy.
# A better approach (future work) would compute k_i across ALL NSW heathland
# cells from the fire shapefile, giving thousands of data points from which
# to estimate alpha/beta. That spatial Empirical Bayes approach borrows
# strength from the full landscape rather than just the M. nodosa footprint.

k_obs <- mn_grid$n_wildfire_years_sig
lambda_obs <- k_obs[k_obs > 0] / OBS_PERIOD   # observed rates (non-zero cells only)

if (length(lambda_obs) >= 10L) {
  mu_obs  <- mean(lambda_obs)
  var_obs <- var(lambda_obs)
  if (var_obs > 0 && mu_obs^2 / var_obs > 0.5) {
    alpha_eb <- mu_obs^2 / var_obs
    beta_eb  <- mu_obs  / var_obs
    message(sprintf(
      "  Empirical Bayes prior: Gamma(%.2f, %.1f)  mean=%.3f (FRI=%.0f yr)",
      alpha_eb, beta_eb, alpha_eb/beta_eb, beta_eb/alpha_eb))
    message(sprintf(
      "  Fixed prior was:        Gamma(%.2f, %.1f)  mean=%.3f (FRI=%.0f yr)",
      PRIOR_ALPHA, PRIOR_BETA, PRIOR_ALPHA/PRIOR_BETA, PRIOR_BETA/PRIOR_ALPHA))
    message("  Recomputing fire_prob_annual with Empirical Bayes prior...")
    # Re-compute posteriors with fitted prior
    mn_grid <- mn_grid |>
      mutate(fire_prob_annual = pmin(FIRE_MAX, pmax(FIRE_MIN,
               (alpha_eb + n_wildfire_years_sig) / (beta_eb + OBS_PERIOD))))
  } else {
    message("  Empirical Bayes: insufficient variance in lambda_obs; keeping fixed prior.")
  }
} else {
  message("  Empirical Bayes: too few non-zero cells for stable fit; keeping fixed prior.")
}
# ---- End Empirical Bayes -----------------------------------------------------

# fire_p_fimp: for cells WITH significant fires, use mean_frac_burned.
# For cells WITHOUT significant fires (NA), use the range-wide median of
# cells that do have fire data -- this avoids the inversion where unburned
# cells get a higher p_fimp than cells with many tiny fires.
# (PFIMP_MIN floor still applies to all cells.)
pfimp_default <- median(mn_grid$mean_frac_burned, na.rm = TRUE)
message(sprintf("  PFIMP_DEFAULT from data median: %.3f", pfimp_default))

mn_grid <- mn_grid |>
  mutate(
    fire_p_fimp = pmin(1.0, pmax(PFIMP_MIN,
                       dplyr::coalesce(mean_frac_burned, pfimp_default))),
    FRI_years = round(1 / fire_prob_annual, 1)
  )

message(sprintf("\n  Wildfire summary (AREA_THRESH=%.2f, ~%.0f ha):",
                AREA_THRESH, AREA_THRESH * 10000))
message(sprintf("    Cells with ≥1 significant wildfire: %d (%.0f%%)  [all wildfires: %d]",
                sum(mn_grid$n_wildfire_years_sig > 0),
                100 * mean(mn_grid$n_wildfire_years_sig > 0),
                sum(mn_grid$n_wildfire_years_all > 0)))
message(sprintf("    fire_prob_annual:  median=%.3f  IQR=[%.3f, %.3f]",
                median(mn_grid$fire_prob_annual),
                quantile(mn_grid$fire_prob_annual, 0.25),
                quantile(mn_grid$fire_prob_annual, 0.75)))
message(sprintf("    FRI_years:         median=%.0f    IQR=[%.0f, %.0f]",
                median(mn_grid$FRI_years),
                quantile(mn_grid$FRI_years, 0.25),
                quantile(mn_grid$FRI_years, 0.75)))
message(sprintf("    fire_p_fimp:       median=%.3f  IQR=[%.3f, %.3f]",
                median(mn_grid$fire_p_fimp, na.rm=TRUE),
                quantile(mn_grid$fire_p_fimp, 0.25, na.rm=TRUE),
                quantile(mn_grid$fire_p_fimp, 0.75, na.rm=TRUE)))
message(sprintf("    mean_frac_burned:  median=%.3f  range=[%.3f, %.3f]",
                median(mn_grid$mean_frac_burned, na.rm=TRUE),
                min(mn_grid$mean_frac_burned, na.rm=TRUE),
                max(mn_grid$mean_frac_burned, na.rm=TRUE)))

# ============================================================================
# STEP 4: SDM suitability → disease triangle parameters
# ============================================================================
message("\nStep 4: Extracting MaxEnt SDM values...")

sdm      <- rast(SDM_TIFF)
mn_pts   <- vect(mn_grid, geom = c("lon_centre","lat_centre"), crs = "EPSG:4326")
sdm_vals <- terra::extract(sdm, mn_pts)[, 2L]
sdm_vals[sdm_vals < -1e30] <- NA_real_
mn_grid$sdm_suitability <- sdm_vals

# Verify Yarringully value
ynr_val <- terra::extract(sdm,
  vect(data.frame(lon=YNR_LON, lat=YNR_LAT),
       geom=c("lon","lat"), crs="EPSG:4326"))[,2L]
message(sprintf("  YNR check: SDM=%.4f (expected %.4f, diff=%.4f)",
                ynr_val, SDM_YARR, abs(ynr_val - SDM_YARR)))

# Map SDM → annual_env Beta parameters, anchored at Yarringully.
# annual_env_mean scales linearly with SDM suitability all the way down.
# The ONLY floor needed is numerical: Beta(alpha, beta) requires alpha > 0.
# With CONC=20, floor = 1/(CONC*10) = 0.005 gives alpha_min = 0.1, which is
# well within R's rbeta() numerical range. Sites with SDM near zero get
# annual_env_mean ≈ 0.005-0.02 — effectively rust-absent pressure without
# needing a separate rust_start_year=Inf for those cells.
# Raising this floor (e.g. to 0.10 as originally) would artificially impose
# substantial rust pressure at sites where the SDM says rust can barely persist.
FLOOR_ENV <- 1 / (CONC * 10)   # 0.005 with CONC=20; keeps alpha >= 0.1

mn_grid <- mn_grid |>
  mutate(
    annual_env_mean  = pmin(0.95, pmax(FLOOR_ENV,
                            MU_YARR * (sdm_suitability / SDM_YARR))),
    annual_env_alpha = round(annual_env_mean * CONC, 2L),
    annual_env_beta  = round((1 - annual_env_mean) * CONC, 2L)
  )

message(sprintf("  annual_env_mean: median=%.3f  IQR=[%.3f, %.3f]",
                median(mn_grid$annual_env_mean, na.rm=TRUE),
                quantile(mn_grid$annual_env_mean, 0.25, na.rm=TRUE),
                quantile(mn_grid$annual_env_mean, 0.75, na.rm=TRUE)))

# ============================================================================
# STEP 5: Save
# ============================================================================
message("\nStep 5: Saving site parameter table...")
write_csv(mn_grid, "spatial/data/mn_site_params.csv")
message("  Saved: spatial/data/mn_site_params.csv")

# ============================================================================
# STEP 5b: Ensemble weights — how much should each cell contribute?
# ============================================================================
# The 375 occupied cells are not equally representative of the species'
# abundance. Western/inland cells may have genuine but sparse M. nodosa
# populations; central coastal cells likely hold the bulk of the species.
# Running the ensemble with uniform sampling over cells would give equal
# weight to a single roadside record and a cell with 50 records.
#
# We estimate relative cell importance two ways and store both:
#
# weight_n_raw:
#   Proportional to n_raw (ALA raw record count per cell, from mn_grid.csv).
#   Simple and transparent; biased by survey effort but directionally correct
#   — heavily botanised areas with genuinely dense populations accumulate more
#   records.
#
# weight_kde:
#   Kernel density of the 5000+ cleaned occurrence records (before grid
#   thinning), evaluated at each cell centre. Smooths over sampling noise and
#   edge artifacts from isolated western records. Still effort-biased but
#   spatially more stable than raw counts.
#   Uses MASS::kde2d; bandwidth selected by 2D normal-reference rule.
#
# LIMITATION: neither metric is true abundance. Both are proxies. The KDE
# approach is the recommended default for the ensemble. If independent
# density data become available (plot surveys, remote sensing cover) these
# should replace these proxies.
#
# USAGE IN ENSEMBLE (03_build_ensemble_params.R):
#   # Weighted random sample of cells for a 5000-run ensemble
#   cell_draw <- sample(nrow(site_params), size = 5000, replace = TRUE,
#                       prob = site_params$weight_kde)
#   # OR: run a fixed number of seeds per cell and weight results post-hoc
#   #   weighted_mean_decline = sum(w * pct_decline) / sum(w)
#   #   where w = site_params$weight_kde

message("Step 5b: Computing ensemble weights...")

# Weight 1: n_raw (normalised)
mn_grid$weight_n_raw <- mn_grid$n_raw / sum(mn_grid$n_raw)

# Weight 2: kernel density of cleaned occurrence records (pre-thinning)
if (requireNamespace("MASS", quietly = TRUE)) {
  mn_clean_path <- "spatial/data/mn_clean.csv"
  if (file.exists(mn_clean_path)) {
    mn_clean <- read_csv(mn_clean_path, show_col_types = FALSE)

    # 2D kernel density over occurrence locations
    # Bandwidth: default normal-reference rule (MASS::kde2d h=bandwidth.nrd)
    bw_x <- MASS::bandwidth.nrd(mn_clean$longitude)
    bw_y <- MASS::bandwidth.nrd(mn_clean$latitude)
    kde  <- MASS::kde2d(
      mn_clean$longitude, mn_clean$latitude,
      h = c(bw_x, bw_y),
      n = 200L,
      lims = c(range(mn_grid$lon_centre) + c(-0.5, 0.5),
               range(mn_grid$lat_centre) + c(-0.5, 0.5))
    )

    # Bilinear interpolation at each cell centre
    # kde$z is a matrix (x rows, y cols); use approx to evaluate at points
    interp_kde <- function(px, py, kde_obj) {
      ix <- findInterval(px, kde_obj$x)
      iy <- findInterval(py, kde_obj$y)
      ix <- pmax(1L, pmin(ix, length(kde_obj$x) - 1L))
      iy <- pmax(1L, pmin(iy, length(kde_obj$y) - 1L))
      # Bilinear interpolation
      x1 <- kde_obj$x[ix];     x2 <- kde_obj$x[ix + 1L]
      y1 <- kde_obj$y[iy];     y2 <- kde_obj$y[iy + 1L]
      wx2 <- (px - x1) / (x2 - x1); wx1 <- 1 - wx2
      wy2 <- (py - y1) / (y2 - y1); wy1 <- 1 - wy2
      z11 <- kde_obj$z[cbind(ix,   iy  )]
      z12 <- kde_obj$z[cbind(ix,   iy+1L)]
      z21 <- kde_obj$z[cbind(ix+1L,iy  )]
      z22 <- kde_obj$z[cbind(ix+1L,iy+1L)]
      wx1*(wy1*z11 + wy2*z12) + wx2*(wy1*z21 + wy2*z22)
    }

    kde_vals <- interp_kde(mn_grid$lon_centre, mn_grid$lat_centre, kde)
    kde_vals <- pmax(kde_vals, 0)  # guard against tiny negatives at edges
    mn_grid$weight_kde <- kde_vals / sum(kde_vals)

    message(sprintf(
      "  KDE bandwidth: lon=%.4f, lat=%.4f (normal-reference rule)",
      bw_x, bw_y))
  } else {
    message("  mn_clean.csv not found — using weight_n_raw as KDE fallback")
    mn_grid$weight_kde <- mn_grid$weight_n_raw
  }
} else {
  message("  Package 'MASS' not available — using weight_n_raw as KDE fallback")
  mn_grid$weight_kde <- mn_grid$weight_n_raw
}

# Sanity check: effective sample size (ESS) = 1/Σ(w²)
# ESS close to n_cells → weights nearly uniform (western cells still prominent)
# ESS much lower → weights concentrated on a few high-density cells
ess_raw <- 1 / sum(mn_grid$weight_n_raw^2)
ess_kde <- 1 / sum(mn_grid$weight_kde^2)
message(sprintf(
  "  Effective sample size: weight_n_raw=%.0f / %.0f cells  |  weight_kde=%.0f / %.0f cells",
  ess_raw, nrow(mn_grid), ess_kde, nrow(mn_grid)))

# Re-save with weights included
write_csv(mn_grid, "spatial/data/mn_site_params.csv")
message("  Re-saved: spatial/data/mn_site_params.csv  (with ensemble weights)")

cat("\n=== Site parameter summary ===\n")
cat(sprintf("  Cells:                    %d\n", nrow(mn_grid)))
cat(sprintf("  Cells with significant wildfire:  %d (%.0f%%)  [all: %d]\n",
            sum(mn_grid$n_wildfire_years_sig > 0),
            100 * mean(mn_grid$n_wildfire_years_sig > 0),
            sum(mn_grid$n_wildfire_years_all > 0)))
cat(sprintf("  fire_prob_annual:  median=%.3f  range=[%.3f, %.3f]\n",
            median(mn_grid$fire_prob_annual),
            min(mn_grid$fire_prob_annual), max(mn_grid$fire_prob_annual)))
cat(sprintf("  FRI_years:         median=%.0f    range=[%.0f, %.0f]\n",
            median(mn_grid$FRI_years),
            min(mn_grid$FRI_years), max(mn_grid$FRI_years)))
cat(sprintf("  fire_p_fimp:       median=%.3f  range=[%.3f, %.3f]\n",
            median(mn_grid$fire_p_fimp, na.rm=TRUE),
            min(mn_grid$fire_p_fimp,    na.rm=TRUE),
            max(mn_grid$fire_p_fimp,    na.rm=TRUE)))
cat(sprintf("  annual_env_mean:   median=%.3f  range=[%.3f, %.3f]\n",
            median(mn_grid$annual_env_mean, na.rm=TRUE),
            min(mn_grid$annual_env_mean,    na.rm=TRUE),
            max(mn_grid$annual_env_mean,    na.rm=TRUE)))
cat("  Fixed across sites: microclim_alpha=20, microclim_beta=5\n")
cat("\nKey IBM connections:\n")
cat("  fire_prob_annual → IBM fire_prob_annual (year-to-year fire occurrence)\n")
cat("  fire_p_fimp      → IBM fire_p_fimp      (fraction of population impacted)\n")
cat("  annual_env_alpha/beta → IBM annual_env_alpha/beta (rust pressure)\n")
cat("\nNext: spatial/03_build_ensemble_params.R\n")

# ============================================================================
# STEP 6: Diagnostic maps
# ============================================================================
message("\nStep 6: Making diagnostic maps...")

aus_states <- ozmaps::ozmap_states
xlim <- c(148.5, 154.5); ylim <- c(-37.5, -24.5)

base_map <- list(
  geom_sf(data = aus_states, inherit.aes = FALSE,
          fill = "#f5f3ee", colour = "#cccccc", linewidth = 0.3),
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE),
  theme_minimal(base_size = 10),
  theme(panel.grid = element_line(colour="#e8e8e8", linewidth=0.2),
        legend.position = "bottom", axis.text = element_text(size=7))
)

# Map A: fire_p_fimp (fraction of population impacted per fire)
pA <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre, colour=fire_p_fimp)) +
  base_map +
  geom_point(size=2.5, alpha=0.9) +
  scale_colour_gradient(low="#ffffcc", high="#e34948", name="fire_p_fimp",
                         na.value="#cccccc") +
  labs(title = "fire_p_fimp: fraction of population impacted per fire event",
       subtitle = "Wildfire only, area-weighted. Grey = no wildfire recorded.",
       x=NULL, y=NULL)

# Map B: FRI (fire return interval)
mn_grid_fire <- mn_grid[mn_grid$n_wildfire_years_sig > 0, ]
pB <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre)) +
  base_map +
  geom_point(colour="#cccccc", size=2, alpha=0.7) +  # unburned cells
  geom_point(data=mn_grid_fire,
             aes(colour=FRI_years), size=2.5, alpha=0.9) +
  scale_colour_gradient(low="#e34948", high="#2a78d6", name="FRI (years)") +
  labs(title = "Fire return interval from NPWS wildfire history (1970-2024)",
       subtitle = sprintf("Median FRI = %.0f yr | %.0f%% cells with ≥1 significant wildfire (>%.0f%% of cell)",
                           median(mn_grid$FRI_years),
                           100*mean(mn_grid$n_wildfire_years_sig > 0),
                           AREA_THRESH*100),
       x=NULL, y=NULL)

# Map C: SDM rust suitability
pC <- ggplot(mn_grid[!is.na(mn_grid$sdm_suitability),],
             aes(x=lon_centre, y=lat_centre, colour=sdm_suitability)) +
  base_map +
  geom_point(size=2.5, alpha=0.9) +
  annotate("point", x=YNR_LON, y=YNR_LAT,
           shape=4, size=4, colour="#e34948", stroke=1.5) +
  annotate("text", x=YNR_LON-0.3, y=YNR_LAT, label="YNR",
           hjust=1, size=3, colour="#e34948") +
  scale_colour_gradient(low="#ffffcc", high="#006837", name="SDM\nsuitability") +
  labs(title="Myrtle rust MaxEnt suitability",
       subtitle=sprintf("YNR anchor=%.4f | median=%.3f",
                         SDM_YARR, median(mn_grid$sdm_suitability, na.rm=TRUE)),
       x=NULL, y=NULL)

# Map D: ensemble weight (KDE)
pD <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre, colour=weight_kde*1000)) +
  base_map +
  geom_point(size=2.5, alpha=0.9) +
  scale_colour_gradient(low="#f0f0f0", high="#3d0099",
                         name="KDE weight\n(×1000)") +
  labs(title="Ensemble weight (kernel density of ALA records)",
       subtitle=sprintf("ESS=%.0f / %d cells  |  coastal cells upweighted",
                         ess_kde, nrow(mn_grid)),
       x=NULL, y=NULL)

if (requireNamespace("patchwork", quietly=TRUE)) {
  fig <- patchwork::wrap_plots(pA, pB, pC, pD, ncol=2L) +
    patchwork::plot_annotation(
      title="M. nodosa site parameterisation — fire, rust and ensemble weight",
      theme=theme(plot.title=element_text(size=12, face="bold")))
  ggsave("spatial/figures/mn_site_params.png", fig,
         width=14, height=10, dpi=150, bg="white")
} else {
  ggsave("spatial/figures/mn_fire_pfimp.png", pA, width=7, height=6, dpi=150, bg="white")
  ggsave("spatial/figures/mn_fire_fri.png",   pB, width=7, height=6, dpi=150, bg="white")
  ggsave("spatial/figures/mn_sdm_rust.png",   pC, width=7, height=6, dpi=150, bg="white")
  ggsave("spatial/figures/mn_weights.png",    pD, width=7, height=6, dpi=150, bg="white")
}
message("  Saved: spatial/figures/mn_site_params.png")
message("Done.")
