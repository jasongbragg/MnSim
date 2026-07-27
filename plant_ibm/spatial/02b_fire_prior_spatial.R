# spatial/02b_fire_prior_spatial.R
#
# Option B: Spatially-varying Bayesian prior for fire return interval.
#
# Extends the Empirical Bayes approach in 02_cross_reference_layers.R by
# computing the prior mean from a BROADER GRID across all coastal NSW,
# smoothed spatially via a Gaussian focal kernel. This borrows information
# from neighbouring cells outside the M. nodosa footprint, giving a more
# stable prior where M. nodosa records are sparse.
#
# COMPARISON WITH OPTION A
# -------------------------
# Option A (in 02_cross_reference_layers.R):
#   Prior hyperparameters fit from the distribution of k_i across the
#   375 M. nodosa cells. Every cell shares the same prior mean.
#
# Option B (this script):
#   Prior mean is spatially VARYING: each M. nodosa cell gets a prior
#   derived from its local neighbourhood's fire frequency across the
#   broader coastal NSW landscape.
#
# If Options A and B give similar fire_prob_annual estimates, the data
# from the 375 cells are adequate. If they diverge substantially, the
# local data are too sparse to dominate, and the spatial prior is doing
# significant regularisation work.
#
# FIRE SIZE APPROXIMATION
# -----------------------
# For efficiency, the broader grid uses a fire size pre-filter based on
# total polygon area (AreaHa >= MIN_FIRE_HA) rather than the exact
# cell-specific coverage fraction used in 02_cross_reference_layers.R.
# This is an approximation: a 600ha fire at a cell corner covers much
# less than 5% of the cell. The smoothing step partially averages over
# this noise. For the M. nodosa cells themselves, the exact area-weighted
# estimate from 02_cross_reference_layers.R remains the primary estimate;
# Option B is used only to inform the spatial prior.
#
# USAGE
# -----
# Run AFTER 02_cross_reference_layers.R (reads mn_site_params.csv).
# Edit paths and SMOOTH_SIGMA before running.
#   source("spatial/02b_fire_prior_spatial.R")

# ============================================================================
# Configuration
# ============================================================================
FIRE_SHP    <- "/home/jgb/projects/melsims/utilities/fire/NPWSFireHistory.shp"  # edit to local path after download
SITE_PARAMS    <- "spatial/data/mn_site_params.csv"
GRID_DEG       <- 0.1
OBS_START      <- 1970L
OBS_END        <- 2024L
OBS_PERIOD     <- OBS_END - OBS_START + 1L

# Broader grid extent (coastal NSW + QLD fringe)
BROAD_XMIN     <- 148.0;  BROAD_XMAX <- 155.0
BROAD_YMIN     <- -38.0;  BROAD_YMAX <- -24.0

# Fire size pre-filter: only count fires larger than this (ha).
# AREA_THRESH from Option A = 0.05 of cell area ≈ 500 ha for a 10km cell.
MIN_FIRE_HA    <- 500

# Gaussian kernel bandwidth for spatial smoothing (degrees).
# 0.5° ≈ 50km: borrows from nearby cells, preserves regional gradients
# 1.0° ≈ 100km: more smoothing, less sensitivity to local data
# 0.3° ≈ 30km: mainly local, closer to raw data
SMOOTH_SIGMA   <- 0.5

# Bayesian concentration parameter: equivalent prior sample size
# Higher C → prior dominates over local data; lower C → local data dominates
# C = 55 (= OBS_PERIOD): prior and local data weighted equally
# C = 110: prior is twice as influential as the data
C_CONCENTRATION <- 55L

FIRE_MAX       <- 0.500
FIRE_MIN       <- 0.005

# ============================================================================
# Setup
# ============================================================================
if (!dir.exists("spatial/data"))    dir.create("spatial/data",    recursive = TRUE)
if (!dir.exists("spatial/figures")) dir.create("spatial/figures", recursive = TRUE)

pkgs <- c("sf","terra","dplyr","readr","ggplot2","ozmaps","patchwork")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing: ", paste(missing_pkgs, collapse=", "))

library(sf); library(terra); library(dplyr); library(readr)
library(ggplot2)

# ============================================================================
# STEP 1: Load fire history and pre-filter to large fires
# ============================================================================
message("Step 1: Loading and filtering fire history...")

sf_use_s2(FALSE)
fire_raw <- st_read(FIRE_SHP, quiet=TRUE) |> st_make_valid()
fire     <- st_transform(fire_raw, 4326L)

fire$fire_year <- as.integer(substr(fire$Label, 1L, 4L))
fire$fire_type <- dplyr::case_when(
  grepl("Wildfire",   fire$Label, ignore.case=TRUE) ~ "Wildfire",
  grepl("Prescribed", fire$Label, ignore.case=TRUE) ~ "Prescribed",
  TRUE ~ "Unknown"
)

# Filter: wildfires only, observation period, size threshold
fire_large <- fire[
  !is.na(fire$fire_year) &
  fire$fire_year >= OBS_START &
  fire$fire_year <= OBS_END &
  fire$fire_type == "Wildfire" &
  !is.na(fire$AreaHa) &
  fire$AreaHa >= MIN_FIRE_HA, ]

# Clip to broader bounding box
bbox_broad <- st_as_sfc(st_bbox(c(
  xmin=BROAD_XMIN, xmax=BROAD_XMAX,
  ymin=BROAD_YMIN, ymax=BROAD_YMAX), crs=4326L))
fire_large <- fire_large[
  suppressWarnings(st_intersects(fire_large, bbox_broad, sparse=FALSE)[,1L]), ]

message(sprintf("  Large wildfires (>=%d ha) in region: %d polygons, %d unique years",
                MIN_FIRE_HA, nrow(fire_large), length(unique(fire_large$fire_year))))

# ============================================================================
# STEP 2: Build broader reference raster and rasterise fire counts
# ============================================================================
message("Step 2: Rasterising fire history across broader grid...")

ref_rast <- rast(
  xmin=BROAD_XMIN, xmax=BROAD_XMAX,
  ymin=BROAD_YMIN, ymax=BROAD_YMAX,
  resolution=GRID_DEG, crs="EPSG:4326")

# For each year: rasterise presence of large wildfires, add to cumulative count
fire_count <- rast(ref_rast); values(fire_count) <- 0L
years_present <- sort(unique(fire_large$fire_year))

for (yr in years_present) {
  fy <- fire_large[fire_large$fire_year == yr, ]
  if (nrow(fy) == 0L) next
  # Rasterise: cell = 1 if any large wildfire touched it this year, 0 otherwise
  fy_v    <- vect(fy)
  r_yr    <- suppressWarnings(
    rasterize(fy_v, ref_rast, field=1L, fun="max", background=0L))
  fire_count <- fire_count + r_yr
}

# Convert to annual fire rate (= MLE of lambda, before smoothing/prior)
lambda_raw <- fire_count / OBS_PERIOD

message(sprintf("  Raw lambda across broader grid: mean=%.4f (FRI=%.0f yr)",
                mean(values(lambda_raw), na.rm=TRUE),
                1/mean(values(lambda_raw)[values(lambda_raw)>0], na.rm=TRUE)))

# ============================================================================
# STEP 3: Gaussian spatial smoothing
# ============================================================================
message(sprintf("Step 3: Gaussian smoothing (sigma=%.1f deg = ~%.0f km)...",
                SMOOTH_SIGMA, SMOOTH_SIGMA * 111))

# focalMat() takes d in the raster's native units (degrees for geographic CRS),
# NOT in cell counts. SMOOTH_SIGMA is already in degrees so pass it directly.
# At GRID_DEG=0.1, sigma=0.5 deg → 5-cell SD → ~30-cell weight matrix.
w <- focalMat(lambda_raw, d=SMOOTH_SIGMA, type="Gauss")
message(sprintf("  Weight matrix: %d x %d cells", nrow(w), ncol(w)))
lambda_smooth <- focal(lambda_raw, w=w, fun="sum", na.rm=TRUE)

message(sprintf("  Smoothed lambda: mean=%.4f (FRI=%.0f yr)  SD=%.4f",
                mean(values(lambda_smooth), na.rm=TRUE),
                1/mean(values(lambda_smooth)[values(lambda_smooth)>0], na.rm=TRUE),
                sd(values(lambda_smooth), na.rm=TRUE)))

# ============================================================================
# STEP 4: Extract smoothed prior at M. nodosa cells
# ============================================================================
message("Step 4: Extracting spatial prior at M. nodosa cells...")

mn_grid <- read_csv(SITE_PARAMS, show_col_types=FALSE)
mn_pts  <- vect(mn_grid, geom=c("lon_centre","lat_centre"), crs="EPSG:4326")

lambda_spatial <- extract(lambda_smooth, mn_pts)[, 2L]

# Guard: small positive floor (can't have exactly-zero prior mean)
lambda_spatial <- pmax(lambda_spatial, FIRE_MIN)

mn_grid$lambda_prior_B <- lambda_spatial

# ============================================================================
# STEP 5: Bayesian update with spatially-varying prior
# ============================================================================
message("Step 5: Computing Bayesian posteriors (Options A and B)...")

# Option B posterior: spatially-varying prior
# Prior for cell i: Gamma(lambda_prior_B_i * C, C)
# Posterior: Gamma(lambda_prior_B_i * C + k_i, C + T)
# Posterior mean = (lambda_prior_B_i * C + k_i) / (C + T)

k <- mn_grid$n_wildfire_years_sig   # significant fire years (from Option A run)

mn_grid <- mn_grid |>
  mutate(
    # Option B: spatially-varying prior
    fire_prob_B = pmin(FIRE_MAX, pmax(FIRE_MIN,
      (lambda_prior_B * C_CONCENTRATION + k) /
      (C_CONCENTRATION + OBS_PERIOD))),
    FRI_B       = round(1 / fire_prob_B, 1),

    # Option A for reference (fixed prior from 02_cross_reference_layers.R)
    # Extract from existing column if present, else note it needs recomputing
    fire_prob_A = if ("fire_prob_annual" %in% names(mn_grid))
                    fire_prob_annual else NA_real_,
    FRI_A       = round(1 / fire_prob_A, 1)
  )

# Residuals
mn_grid$delta_AB <- mn_grid$fire_prob_B - mn_grid$fire_prob_A

message(sprintf("\n  Option A — fire_prob_annual: median=%.4f (FRI=%.0f yr)",
                median(mn_grid$fire_prob_A, na.rm=TRUE),
                median(mn_grid$FRI_A, na.rm=TRUE)))
message(sprintf("  Option B — fire_prob_annual: median=%.4f (FRI=%.0f yr)",
                median(mn_grid$fire_prob_B, na.rm=TRUE),
                median(mn_grid$FRI_B, na.rm=TRUE)))
message(sprintf("  Mean |A-B|: %.4f  (max: %.4f)",
                mean(abs(mn_grid$delta_AB), na.rm=TRUE),
                max(abs(mn_grid$delta_AB), na.rm=TRUE)))

# ============================================================================
# STEP 6: Save updated site params
# ============================================================================
message("Step 6: Saving...")

write_csv(mn_grid, "spatial/data/mn_site_params.csv")
message("  Updated: spatial/data/mn_site_params.csv (added fire_prob_B, FRI_B, lambda_prior_B)")

# Also save the smoothed fire frequency raster for inspection
writeRaster(lambda_smooth, "spatial/data/fire_lambda_smooth.tif", overwrite=TRUE)
writeRaster(lambda_raw,    "spatial/data/fire_lambda_raw.tif",    overwrite=TRUE)
message("  Saved: spatial/data/fire_lambda_smooth.tif")

# ============================================================================
# STEP 7: Comparison figures
# ============================================================================
message("Step 7: Making comparison figures...")

aus_states <- ozmaps::ozmap_states
xlim <- c(148.5, 154.5); ylim <- c(-37.5, -24.5)
base_map <- list(
  geom_sf(data=aus_states, inherit.aes=FALSE,
          fill="#f5f3ee", colour="#cccccc", linewidth=0.3),
  coord_sf(xlim=xlim, ylim=ylim, expand=FALSE),
  theme_minimal(base_size=10),
  theme(panel.grid=element_line(colour="#e8e8e8", linewidth=0.2),
        legend.position="bottom", axis.text=element_text(size=7))
)

# Map A: Option A fire probability
pA <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre, colour=fire_prob_A)) +
  base_map +
  geom_point(size=2.5) +
  scale_colour_gradient2(low="#2a78d6", mid="#fdae6b", high="#e34948",
                          midpoint=0.05, name="fire_prob\n(Option A)",
                          na.value="#cccccc") +
  labs(title="Option A: Empirical Bayes (M. nodosa cells only)",
       subtitle=sprintf("median FRI = %.0f yr", median(mn_grid$FRI_A, na.rm=TRUE)),
       x=NULL, y=NULL)

# Map B: Option B fire probability
pB <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre, colour=fire_prob_B)) +
  base_map +
  geom_point(size=2.5) +
  scale_colour_gradient2(low="#2a78d6", mid="#fdae6b", high="#e34948",
                          midpoint=0.05, name="fire_prob\n(Option B)",
                          na.value="#cccccc") +
  labs(title="Option B: Spatial prior (smoothed landscape)",
       subtitle=sprintf("median FRI = %.0f yr  |  sigma=%.1f deg",
                         median(mn_grid$FRI_B, na.rm=TRUE), SMOOTH_SIGMA),
       x=NULL, y=NULL)

# Map C: Difference B - A
pC <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre, colour=delta_AB)) +
  base_map +
  geom_point(size=2.5) +
  scale_colour_gradient2(low="#2a78d6", mid="white", high="#e34948",
                          midpoint=0, name="B - A\n(fire_prob)") +
  labs(title="Difference: Option B minus Option A",
       subtitle="Red = B gives higher fire prob than A; blue = lower",
       x=NULL, y=NULL)

# Scatter: A vs B
pD <- ggplot(mn_grid, aes(x=fire_prob_A, y=fire_prob_B)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="#888780") +
  geom_point(colour="#2a78d6", alpha=0.6, size=1.5) +
  geom_smooth(method="lm", colour="#e34948", se=FALSE, linewidth=0.8) +
  scale_x_continuous(limits=c(0, FIRE_MAX)) +
  scale_y_continuous(limits=c(0, FIRE_MAX)) +
  labs(x="fire_prob_annual  Option A",
       y="fire_prob_annual  Option B",
       title="A vs B comparison",
       subtitle=sprintf("1:1 line shown dashed | C=%d | sigma=%.1f°",
                         C_CONCENTRATION, SMOOTH_SIGMA)) +
  theme_minimal(base_size=10)

if (requireNamespace("patchwork", quietly=TRUE)) {
  fig <- patchwork::wrap_plots(pA, pB, pC, pD, ncol=2L) +
    patchwork::plot_annotation(
      title="Fire return interval: Option A (Empirical Bayes) vs Option B (Spatial prior)",
      theme=theme(plot.title=element_text(size=11, face="bold")))
  ggsave("spatial/figures/fire_AB_comparison.png", fig,
         width=14, height=11, dpi=150, bg="white")
} else {
  ggsave("spatial/figures/fire_optA.png", pA, width=7, height=6, dpi=150, bg="white")
  ggsave("spatial/figures/fire_optB.png", pB, width=7, height=6, dpi=150, bg="white")
}
message("  Saved: spatial/figures/fire_AB_comparison.png")

# Also plot the raw vs smoothed landscape fire frequency
lambda_df <- as.data.frame(lambda_smooth, xy=TRUE)
names(lambda_df)[3] <- "lambda"
lambda_df <- lambda_df[!is.na(lambda_df$lambda) & lambda_df$lambda > 0, ]

p_landscape <- ggplot(lambda_df, aes(x=x, y=y, fill=lambda)) +
  geom_raster() +
  geom_point(data=mn_grid, aes(x=lon_centre, y=lat_centre),
             inherit.aes=FALSE, colour="white", size=0.8, alpha=0.7) +
  geom_sf(data=aus_states, inherit.aes=FALSE,
          fill=NA, colour="#333333", linewidth=0.3) +
  coord_sf(xlim=xlim, ylim=ylim, expand=FALSE) +
  scale_fill_gradient(low="#ffffcc", high="#e34948", name="lambda\n(fire/yr)") +
  labs(title="Smoothed wildfire frequency — broader coastal NSW landscape",
       subtitle=sprintf("White dots = M. nodosa cells  |  sigma=%.1f° = ~%.0f km",
                         SMOOTH_SIGMA, SMOOTH_SIGMA*111),
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom")

ggsave("spatial/figures/fire_landscape.png", p_landscape,
       width=7, height=8, dpi=150, bg="white")
message("  Saved: spatial/figures/fire_landscape.png")
message("\nDone. Key outputs: mn_site_params.csv (fire_prob_B added),")
message("  fire_AB_comparison.png, fire_landscape.png")
