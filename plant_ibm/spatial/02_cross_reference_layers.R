# spatial/02_cross_reference_layers.R
#
# Cross-reference 375 M. nodosa occupied cells with:
#   (A) NSW NPWS Fire History → fire_prob_annual and fire_p_fimp per cell
#   (B) Myrtle rust MaxEnt SDM → annual_env Beta parameters per cell
#
# FIRE (docs/methodological_notes.md D1-D5; code in spatial/fire/)
# ----
# Built beforehand by spatial/fire/02_fire_cell_years.R: the burnt fraction of
# each cell's land in each year, NSW (NPWS Fire History) and Queensland
# (QRSC fire scars, QPWS planned burns) on the same footing. Here
# site_fire_regime() turns those tables into each cell's fire regime:
#
#   fire_prob_annual = annual probability of a fire year (>= 10% of the cell
#                      burnt): Gamma-Poisson posterior mean, prior fitted by
#                      maximum likelihood over all cells (empirical Bayes)
#   fire_p_fimp      = mean burnt fraction in fire years (median of burnt
#                      cells where a cell has none)
#   sd_frac_burned   = between-fire sd of the burnt fraction
#
# Default regime = WILDFIRE ONLY: NSW wildfires 1987-2024; Qld fire scars
# minus QPWS planned burns, 2014-2024 (the automated 1987-2013 mapping is
# rejected, D4). All fire and planned burns are written alongside with
# suffixes _all and _planned (comparison and management scenarios, D5).
# Checks behind these choices: spatial/fire/03_fire_checks.R.
#
# SDM CALIBRATION
# ---------------
# Anchor (RUST_ANCHOR in params.R): Yarringully NR, mean pressure mu = 0.70.
# annual_env_mean(cell) = mu * SDM_cell / SDM_anchor_cell, capped at 0.95,
# with the anchor cell's SDM read from the raster each run (never typed in),
# so the anchor cell is always exactly mu: Beta(14, 6), as in params.R.
# Beta concentration (alpha + beta = conc) held constant across cells.
# See spatial/rust_pressure.R; to refresh only this step after an SDM
# update, run spatial/02c_update_rust_pressure.R.
#
# OUTPUT
# ------
# spatial/data/mn_site_params.csv   one row per occupied cell, all IBM params
#   (fire columns: see spatial/fire/R/regime.R, site_fire_regime())
# spatial/figures/mn_site_params.png two-panel diagnostic map
#
# USAGE
# -----
# Fire tables first (once, or when fire data change):
#   source("spatial/fire/01_download_qld.R")      # Qld downloads
#   source("spatial/fire/02_fire_cell_years.R")   # burnt fraction per cell-year
# Then, with SDM_TIFF set below:
#   source("spatial/02_cross_reference_layers.R")

# ============================================================================
# Configuration — edit paths before running
# ============================================================================
#SDM_TIFF    <- "/home/jgb/projects/melsims/utilities/rust/maxent_MR.tiff"       # edit to local path
SDM_TIFF    <- "/home/jgb/projects/mrsdm/output/projections/AP_suitability_current.tif"
MN_GRID   <- "spatial/data/mn_grid.csv"
GRID_DEG  <- 0.1          # must match 01_get_melaleuca_occurrences.R

# Fire settings (record years, fire-year threshold, prior, floors): see
# spatial/fire/fire_config.R. NPWS shapefile path: NPWS_SHP there (or the
# NPWS_SHP environment variable); needed only by spatial/fire/02_fire_cell_years.R.

# Rust anchor (site, mu, conc, cap) lives in params.R: RUST_ANCHOR.
source("params.R")
source("spatial/rust_pressure.R")
source("spatial/fire/fire_config.R")
for (f in c("grid", "estimate", "regime")) source(file.path("spatial/fire/R", paste0(f, ".R")))
YNR_LON   <- RUST_ANCHOR$lon
YNR_LAT   <- RUST_ANCHOR$lat

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
# STEP 1: Load M. nodosa grid
# ============================================================================
message("Step 1: Loading M. nodosa grid...")
mn_grid <- read_csv(MN_GRID, show_col_types = FALSE)
message(sprintf("  %d occupied cells at %.2f deg", nrow(mn_grid), GRID_DEG))
sf_use_s2(FALSE)

# ============================================================================
# STEP 2-3: Fire regime per cell (from spatial/data/fire/cell_years_*.csv)
# ============================================================================
message("Step 2: Fire regime per cell (wildfire; all fire and planned burns alongside)...")
need <- file.path(FIRE_DATA_DIR, paste0("cell_years_", c("nsw_wild", "nsw_planned", "nsw_all",
                                                         "qld_wild", "qld_planned", "qld_all"), ".csv"))
if (!all(file.exists(need)))
  stop("Fire cell-year tables missing: run spatial/fire/02_fire_cell_years.R first\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
fire_cells <- load_fire_cells()
if (!setequal(fire_cells$cell_id, mn_grid$cell_id))
  stop("spatial/data/fire/cells.rds is from a different grid: delete it and rerun spatial/fire/02_fire_cell_years.R")
fire_reg <- site_fire_regime(fire_cells)
pri <- attr(fire_reg, "priors")
for (i in seq_len(nrow(pri)))
  message(sprintf("  prior (%s, %s): Gamma(%.2f, %.1f), mean %.0f yr between fire years",
                  pri$regime[i], pri$group[i], pri$alpha[i], pri$beta[i], 1 / pri$prior_mean[i]))
mn_grid <- mn_grid[, setdiff(names(mn_grid), SITE_FIRE_COLUMNS())]
mn_grid <- left_join(mn_grid, fire_reg, by = "cell_id")

message(sprintf("\n  Wildfire summary (fire year = >= %.0f%% of the cell's land burnt):", 100 * AREA_THRESH))
for (st in c("NSW", "Qld")) {
  d <- mn_grid[mn_grid$fire_state == st, ]
  message(sprintf("    %s (%d cells, %s): >= 1 fire year in %.0f%% of cells; median FRI %.0f yr; median burnt fraction %.3f",
                  st, nrow(d), d$fire_record[1], 100 * mean(d$n_wildfire_years_sig > 0),
                  median(d$FRI_years), median(d$mean_frac_burned, na.rm = TRUE)))
}
message(sprintf("    fire_p_fimp: median=%.3f  IQR=[%.3f, %.3f]",
                median(mn_grid$fire_p_fimp), quantile(mn_grid$fire_p_fimp, 0.25),
                quantile(mn_grid$fire_p_fimp, 0.75)))

# ============================================================================
# STEP 4: SDM suitability → disease triangle parameters
# ============================================================================
message("\nStep 4: SDM suitability -> rust pressure (anchored)...")

sdm     <- rast(SDM_TIFF)
mn_grid <- apply_sdm_to_sites(mn_grid, sdm, RUST_ANCHOR)

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
# USAGE IN ENSEMBLE (now runner/range_ensemble.R, which uses weight_n_raw):
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
cat("\nNext: runscripts/range_ensemble.run.code (runner/range_ensemble.R)\n")

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
       subtitle = "Wildfire only (NSW 1987-2024, Qld 2014-2024). Cells without a fire year: range median.",
       x=NULL, y=NULL)

# Map B: FRI (fire return interval)
mn_grid_fire <- mn_grid[mn_grid$n_wildfire_years_sig > 0, ]
pB <- ggplot(mn_grid, aes(x=lon_centre, y=lat_centre)) +
  base_map +
  geom_point(colour="#cccccc", size=2, alpha=0.7) +  # unburned cells
  geom_point(data=mn_grid_fire,
             aes(colour=FRI_years), size=2.5, alpha=0.9) +
  scale_colour_gradient(low="#e34948", high="#2a78d6", name="FRI (years)") +
  labs(title = "Fire return interval (NSW NPWS 1987-2024; Qld scars 2014-2024)",
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
  labs(title="Myrtle rust SDM suitability",
       subtitle=sprintf("anchor cell=%.4f | median=%.3f",
                         mn_grid$sdm_suitability[mn_grid$rust_anchor_cell],
                         median(mn_grid$sdm_suitability, na.rm=TRUE)),
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
