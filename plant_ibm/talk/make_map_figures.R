# talk/make_map_figures.R
#
# Seminar maps: (1) wildfire history and the fire regime given to the model,
# (2) the myrtle rust SDM and the rust pressure given to the model.
#
#   cd plant_ibm
#   /opt/R/4.3.1/bin/Rscript talk/make_map_figures.R
#
# Source data (same files as spatial/02_cross_reference_layers.R; override
# with environment variables if they live elsewhere on this machine):
#   MN_FIRE_SHP  NSW NPWS Fire History shapefile
#   MN_SDM_TIFF  rust SDM projection (mrsdm/output/projections/...)
# Without them the script still runs, using the cell-level values already in
# spatial/data/ (the right-hand panels) and the 0.1-degree fire raster from
# 02b for the left fire panel; the full-resolution panels need the sources.
#
# OUTPUTS (talk/figures/)
#   fig6_fire_map.png   left: years with wildfire 1970-2024 (NPWS, ~1 km)
#                       right: mean fire interval used per M. nodosa cell
#   fig7_rust_sdm.png   left: SDM suitability with M. nodosa records
#                       right: mean annual rust pressure used per cell

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(terra); library(sf)
})

FIRE_SHP <- Sys.getenv("MN_FIRE_SHP",
  "/home/jgb/projects/melsims/utilities/fire/NPWSFireHistory.shp")
SDM_TIFF <- Sys.getenv("MN_SDM_TIFF",
  "/home/jgb/projects/mrsdm/output/projections/AP_suitability_current.tif")
FIG_DIR   <- Sys.getenv("MN_FIG_DIR", "talk/figures/")
CACHE_DIR <- Sys.getenv("MN_MAP_CACHE", "results/talk_maps/")
SITE_CSV  <- "spatial/data/mn_site_params.csv"
OCC_CSV   <- "spatial/data/mn_clean.csv"
LAMBDA_TIF <- "spatial/data/fire_lambda_raw.tif"

# Map window: NSW part of the range (set XLIM/YLIM to c(148, 154.2), c(-37.6, -24)
# to include south-east Queensland)
XLIM <- c(148.0, 154.0)
YLIM <- c(-37.6, -28.0)
FIRE_RES  <- 0.01          # degrees (~1 km) for the wildfire-count raster
OBS_START <- 1970L; OBS_END <- 2024L
YNR <- c(lon = 153.21967, lat = -29.09357)   # Yarringully NR (Pegg et al. site)

dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
sf_use_s2(FALSE)

INK <- "#2b2b2b"
theme_map <- function(base = 14) {
  theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold", colour = INK),
          plot.subtitle = element_text(colour = "#555555", size = base * 0.85),
          plot.caption = element_text(colour = "#777777", size = base * 0.65, hjust = 0),
          panel.grid = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "#eaf1f6", colour = NA),
          legend.title = element_text(size = base * 0.8))
}
wrap <- function(x, width = 120) paste(strwrap(x, width), collapse = "\n")

land <- ggplot2::map_data("world", region = "Australia")
base_layers <- function() list(
  geom_polygon(data = land, aes(long, lat, group = group), fill = "#f4f2ec",
               colour = "#b8b1a2", linewidth = 0.3, inherit.aes = FALSE))
top_layers <- function() list(
  annotate("point", x = YNR[["lon"]], y = YNR[["lat"]], shape = 24, size = 3.4,
           fill = "white", colour = INK, stroke = 0.9),
  annotate("label", x = YNR[["lon"]] - 0.15, y = YNR[["lat"]], hjust = 1,
           label = "Yarringully NR\n(Pegg et al.)", size = 3.2, label.size = 0,
           fill = scales::alpha("white", 0.8), lineheight = 0.9),
  coord_quickmap(xlim = XLIM, ylim = YLIM, expand = FALSE))

sites <- read.csv(SITE_CSV)
sites <- sites[sites$lon_centre >= XLIM[1] - 0.1 & sites$lon_centre <= XLIM[2] + 0.1 &
               sites$lat_centre >= YLIM[1] - 0.1 & sites$lat_centre <= YLIM[2] + 0.1, ]
occ <- read.csv(OCC_CSV)
occ <- occ[occ$longitude >= XLIM[1] & occ$longitude <= XLIM[2] &
           occ$latitude  >= YLIM[1] & occ$latitude  <= YLIM[2], ]
cell_size <- 0.1   # model grid (deg)

# ============================================================================
# Fire
# ============================================================================
wildfire_count_raster <- function() {
  cache <- file.path(CACHE_DIR, sprintf("wildfire_years_%g.tif", FIRE_RES))
  if (file.exists(cache)) return(rast(cache))
  message("Reading ", FIRE_SHP, " ...")
  f <- st_read(FIRE_SHP, quiet = TRUE) |> st_make_valid() |> st_transform(4326L)
  f$fire_year <- suppressWarnings(as.integer(substr(f$Label, 1L, 4L)))
  f <- f[!is.na(f$fire_year) & f$fire_year >= OBS_START & f$fire_year <= OBS_END &
         grepl("Wildfire", f$Label, ignore.case = TRUE), ]
  win <- st_as_sfc(st_bbox(c(xmin = XLIM[1], xmax = XLIM[2], ymin = YLIM[1], ymax = YLIM[2]),
                           crs = 4326L))
  f <- f[suppressWarnings(st_intersects(f, win, sparse = FALSE)[, 1L]), ]
  tmpl <- rast(ext(XLIM[1], XLIM[2], YLIM[1], YLIM[2]), resolution = FIRE_RES,
               crs = "EPSG:4326")
  cnt <- tmpl; values(cnt) <- 0
  for (y in sort(unique(f$fire_year))) {
    v <- vect(f[f$fire_year == y, ])
    cnt <- cnt + rasterize(v, tmpl, field = 1, background = 0, touches = FALSE)
  }
  names(cnt) <- "years_burnt"
  writeRaster(cnt, cache, overwrite = TRUE)
  cnt
}

fig_fire <- function() {
  have_shp <- file.exists(FIRE_SHP)
  if (have_shp) {
    r  <- wildfire_count_raster()
    df <- as.data.frame(r, xy = TRUE)
    df$x <- round(df$x, 6); df$y <- round(df$y, 6)
    df <- df[df$years_burnt > 0, ]
    df$cls <- cut(df$years_burnt, c(0, 1, 2, 3, 4, Inf),
                  labels = c("1", "2", "3", "4", "5+"))
    pA <- ggplot() + base_layers() +
      geom_raster(data = df, aes(x, y, fill = cls)) +
      scale_fill_manual(values = c("#fdd49e", "#fc8d59", "#e34a33", "#b30000", "#5c0000"),
                        name = "Years with\nwildfire") +
      top_layers() +
      labs(title = "Mapped wildfire, 1970-2024",
           subtitle = "NSW NPWS Fire History, ~1 km;\nwildfires only (prescribed burns excluded)")
  } else {
    r  <- crop(rast(LAMBDA_TIF), ext(XLIM[1], XLIM[2], YLIM[1], YLIM[2]))
    df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "rate"
    df <- df[df$rate > 0, ]
    pA <- ggplot() + base_layers() +
      geom_tile(data = df, aes(x, y, fill = rate), width = 0.1, height = 0.1) +
      scale_fill_gradientn(colours = c("#fdd49e", "#fc8d59", "#e34a33", "#7a0000"),
                           name = "Large wildfires\nper year") +
      top_layers() +
      labs(title = "Mapped wildfire, 1970-2024",
           subtitle = "NSW NPWS Fire History: large wildfires\n(>= 500 ha) per year, 0.1-degree cells")
  }

  s <- sites[!is.na(sites$fire_prob_annual), ]
  s$fri <- 1 / s$fire_prob_annual
  prior <- s$n_wildfire_years_sig == 0
  pB <- ggplot() + base_layers() +
    geom_tile(data = s[prior, ], aes(lon_centre, lat_centre),
              width = cell_size, height = cell_size, fill = "#bdbdbd",
              colour = "white", linewidth = 0.15) +
    geom_tile(data = s[!prior, ], aes(lon_centre, lat_centre, fill = fri),
              width = cell_size, height = cell_size, colour = "white", linewidth = 0.15) +
    scale_fill_gradientn(colours = c("#7a0000", "#e34a33", "#fc8d59", "#fdd49e"),
                         trans = "log10", breaks = c(5, 10, 20, 40),
                         name = "Mean years\nbetween fires") +
    top_layers() +
    labs(title = "Fire regime given to the model",
         subtitle = sprintf("%d M. nodosa cells (0.1 degree); grey = no fire year\nrecorded (estimate then: one fire per %.0f-%.0f yr)",
                            nrow(s), 1 / max(s$fire_prob_annual[prior]), 1 / min(s$fire_prob_annual[prior])))

  (pA | pB) + plot_annotation(
    caption = wrap(paste0(
      "Left: NSW only. Right: Bayesian estimate per cell from wildfire records, NSW NPWS 1987-2024 and Qld fire scars 2014-2024 (Gamma-Poisson, prior fitted across all cells; a cell with no fire year gets ",
      sprintf("one fire per %.0f-%.0f yr, longer where the record is longer).",
              1 / max(s$fire_prob_annual[prior]), 1 / min(s$fire_prob_annual[prior])),
      " Each model fire also burns a cell-specific share of plants, from the mean burnt fraction.",
      if (!have_shp) " Left panel at 0.1 degree; run with the NPWS shapefile for the 1-km version." else ""), 170)) &
    theme_map()
}

# ============================================================================
# Rust SDM
# ============================================================================
fig_sdm <- function() {
  have_sdm <- file.exists(SDM_TIFF)
  s <- sites[!is.na(sites$annual_env_mean), ]
  pB <- ggplot() + base_layers() +
    geom_tile(data = s, aes(lon_centre, lat_centre, fill = annual_env_mean),
              width = cell_size, height = cell_size, colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 0.75),
                         name = "Mean annual\nrust pressure") +
    top_layers() +
    labs(title = "Rust pressure given to the model",
         subtitle = "Mean annual pathogen draw per cell;\nanchored at Yarringully (0.70)")
  if (have_sdm) {
    r <- rast(SDM_TIFF)
    r <- crop(r, ext(XLIM[1] - 0.2, XLIM[2] + 0.2, YLIM[1] - 0.2, YLIM[2] + 0.2))
    r[r < -1e30] <- NA
    if (global(r, "max", na.rm = TRUE)[1, 1] > 1.5) r <- r / 1000   # biomod2 0-1000 scale
    if (ncell(r) > 2e6) r <- aggregate(r, fact = ceiling(sqrt(ncell(r) / 2e6)), fun = "mean")
    land_v <- vect(sf::st_as_sf(maps::map("world", "Australia", plot = FALSE, fill = TRUE)))
    r  <- mask(r, project(land_v, crs(r)))
    df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "suit"
    df$x <- round(df$x, 6); df$y <- round(df$y, 6)
    pA <- ggplot() + base_layers() +
      geom_raster(data = df, aes(x, y, fill = suit)) +
      scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1),
                           name = "Climatic\nsuitability") +
      geom_point(data = occ, aes(longitude, latitude), size = 0.35, colour = "#1b6fb3",
                 alpha = 0.6) +
      top_layers() +
      labs(title = "Myrtle rust species distribution model",
           subtitle = "A. psidii, current climate (biomod2 ensemble);\nblue: M. nodosa records")
    out <- (pA | pB)
  } else {
    s2 <- sites[!is.na(sites$sdm_suitability), ]
    pA <- ggplot() + base_layers() +
      geom_tile(data = s2, aes(lon_centre, lat_centre, fill = sdm_suitability),
                width = cell_size, height = cell_size, colour = "white", linewidth = 0.15) +
      scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1),
                           name = "Climatic\nsuitability") +
      top_layers() +
      labs(title = "Myrtle rust SDM at M. nodosa cells",
           subtitle = "Austropuccinia psidii suitability (biomod2 ensemble);\nrun with the SDM raster for the full surface")
    out <- pA
  }
  out + plot_annotation(
    caption = wrap(paste0(
      "Rust pressure scales linearly with SDM suitability, anchored at Yarringully NR (SDM 0.99, pressure 0.70 from field incidence); ",
      "it is the pathogen vertex of the disease triangle in each simulated year. ",
      "Suitability is near-saturated across most coastal cells, so the SDM mainly separates the coastal range from drier inland outliers."),
      if (file.exists(SDM_TIFF)) 170 else 85)) &
    theme_map()
}

ggsave(file.path(FIG_DIR, "fig6_fire_map.png"), fig_fire(), width = 13, height = 8.5,
       dpi = 200, bg = "white")
message("  wrote ", file.path(FIG_DIR, "fig6_fire_map.png"),
        if (!file.exists(FIRE_SHP)) "  [fallback: no NPWS shapefile]" else "")
ggsave(file.path(FIG_DIR, "fig7_rust_sdm.png"), fig_sdm(),
       width = if (file.exists(SDM_TIFF)) 13 else 7, height = 8.5,
       dpi = 200, bg = "white")
message("  wrote ", file.path(FIG_DIR, "fig7_rust_sdm.png"),
        if (!file.exists(SDM_TIFF)) "  [fallback: no SDM raster]" else "")
