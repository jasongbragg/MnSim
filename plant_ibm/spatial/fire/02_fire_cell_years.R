# spatial/fire/02_fire_cell_years.R -- burnt fraction of each occupied cell in
# each year, NSW and Qld on the same footing (docs/methodological_notes.md D1).
# Run from plant_ibm/ after 01_get_melaleuca_occurrences.R (the grid) and
# spatial/fire/01_download_qld.R.
#
# The tables are shared by all species (notes T13): they cover every cell in
# any species' grid (spatial/data/<sp>/grid.csv, numbered on the master grid,
# spatial/cell_ids.R). Each run adds only the cells not yet covered, so adding
# a species computes just its new cells. A cell's values depend on that cell
# alone, so adding cells never changes the rows already there.
#
# Writes to spatial/data/fire/ (kept in git):
#   cells.csv          the cells covered: cell_id, centre, state, land area
#                      (FIRE_CELLS_CSV; appended last, so a cell is listed only
#                      once all six tables hold its rows)
#   cell_years_<t>.csv one row per cell-year with any burnt area:
#     nsw_wild     NPWS wildfires                      } fraction of the cell's
#     nsw_planned  NPWS prescribed burns                } land within its state
#     nsw_all      all NPWS fire (incl. unlabelled)     }
#     qld_wild     Qld scars minus QPWS planned burns   }
#     qld_planned  QPWS planned burns                   }
#     qld_all      Qld scars, all fire                  }
#   qld_years/<batch>/qld_<y>.rds  per-year Qld cache for one batch of new
#                      cells (not in git)
# To recompute everything, delete cells.csv and the six tables. Slow: NSW a
# few minutes per table, Qld a few minutes per year (recent years are large).
source("spatial/fire/fire_config.R")
for (f in c("grid", "bbox", "sources", "cell_years", "regime")) source(file.path("spatial/fire/R", paste0(f, ".R")))
dir.create(FIRE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
cy_file <- function(name) file.path(FIRE_DATA_DIR, paste0("cell_years_", name, ".csv"))

# Add rows for 'ids' to a table, replacing any rows those cells already have
# there (a re-run after an interrupted batch does not duplicate them).
append_cy <- function(x, name, ids) {
  f   <- cy_file(name)
  old <- if (file.exists(f)) utils::read.csv(f) else NULL
  if (!is.null(old)) old <- old[!old$cell_id %in% ids, ]
  utils::write.csv(rbind(old, x[, c("cell_id", "year", "frac_burned")]), f, row.names = FALSE)
}

# ---- Which cells are new ------------------------------------------------------------
need <- all_species_cells()
done <- read_fire_cells()
new  <- need[!need$cell_id %in% done$cell_id, ]
message(sprintf("Cells: %d across all species grids; %d already in the tables; %d to add",
                nrow(need), if (is.null(done)) 0L else nrow(done), nrow(new)))
if (!nrow(new)) {
  message("Nothing to do: ", FIRE_DATA_DIR, "  (next: spatial/02_cross_reference_layers.R)")
} else {

cells  <- load_cells(new)
states <- state_polygons()
message(sprintf("New cells: %d NSW, %d Qld", sum(cells$state == "NSW"), sum(cells$state == "Qld")))
nsw_cells <- cells[cells$state == "NSW", ]
qld_cells <- cells[cells$state == "Qld", ]

# ---- NSW: wildfire, prescribed, all ------------------------------------------------
nsw_jobs <- list(nsw_wild    = "wildfire",
                 nsw_planned = "planned",
                 nsw_all     = c("wildfire", "planned", "unknown"))
if (nrow(nsw_cells)) {
  message("NSW: reading NPWS fire history")
  nsw <- read_nsw_npws(bbox = sf::st_bbox(nsw_cells))
  print(table(type = nsw$type, in_period = nsw$year %in% NSW_YEARS_ALL))
  for (nm in names(nsw_jobs)) {
    message(sprintf("NSW %s -> cells: %s", paste(nsw_jobs[[nm]], collapse = " + "), nm))
    append_cy(cell_year_fractions(nsw[nsw$type %in% nsw_jobs[[nm]], ], nsw_cells, NSW_YEARS_ALL,
                                  "land", state_poly = states[states$state == "NSW", ]),
              nm, nsw_cells$cell_id)
    invisible(gc())
  }
  rm(nsw); invisible(gc())
}

# ---- Queensland, one year at a time ------------------------------------------------
if (nrow(qld_cells)) {
  # The downloads must cover every cell that needs Qld fire data (bbox.R; notes T10).
  check_qld_coverage(c(file.path(FIRE_RAW_DIR, "qld_scars", sprintf("qld_scars_%d.gpkg", QLD_YEARS_ALL)),
                       file.path(FIRE_RAW_DIR, "qpws_fire_history.gpkg")))
  batch  <- substr(digest::digest(sort(qld_cells$cell_id)), 1, 10)
  yr_dir <- file.path(FIRE_DATA_DIR, "qld_years", batch); dir.create(yr_dir, recursive = TRUE, showWarnings = FALSE)
  qld_land_cells <- cells_on_land(qld_cells, states[states$state == "Qld", ])
  qpws <- read_qpws()
  pb   <- qpws[qpws$type == "planned", ]
  message(sprintf("QPWS: %d planned-burn polygons", nrow(pb)))
  for (y in QLD_YEARS_ALL) {
    out <- file.path(yr_dir, sprintf("qld_%d.rds", y))
    if (file.exists(out)) next
    t0 <- Sys.time()
    r <- qld_year_fractions(y, qld_land_cells, pb)
    saveRDS(r, out)
    message(sprintf("    %d: %d cells burnt (%.0f s)", y, sum(r$frac_all > 0),
                    as.numeric(Sys.time() - t0, units = "secs")))
    rm(r); invisible(gc())
  }
  qy <- do.call(rbind, lapply(QLD_YEARS_ALL, function(y) readRDS(file.path(yr_dir, sprintf("qld_%d.rds", y)))))
  pick <- function(col) { d <- qy[qy[[col]] > 0, c("cell_id", "year", col)]; names(d)[3] <- "frac_burned"; d }
  append_cy(pick("frac_wild"), "qld_wild", qld_cells$cell_id)
  append_cy(pick("frac_all"), "qld_all", qld_cells$cell_id)
  append_cy(pick("frac_planned"), "qld_planned", qld_cells$cell_id)
}

# ---- Record the cells as covered (last: only once their rows are written) -----------
add <- sf::st_drop_geometry(cells)[, FIRE_CELL_COLUMNS]
all <- rbind(done, add)
utils::write.csv(all[order(all$cell_id), ], FIRE_CELLS_CSV, row.names = FALSE)
message(sprintf("Done: %d cells in %s  (next: spatial/02_cross_reference_layers.R)", nrow(all), FIRE_CELLS_CSV))
}
