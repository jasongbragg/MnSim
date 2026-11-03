# spatial/fire/02_fire_cell_years.R -- burnt fraction of each occupied cell in
# each year, NSW and Qld on the same footing (docs/methodological_notes.md D1).
# Run from plant_ibm/ after 01_get_melaleuca_occurrences.R (the grid) and
# spatial/fire/01_download_qld.R. Writes to spatial/data/fire/:
#   cells.rds          cells with state and land area (cache)
#   cell_years_<t>.csv one row per cell-year with any burnt area:
#     nsw_wild     NPWS wildfires                      } fraction of the cell's
#     nsw_planned  NPWS prescribed burns                } land within its state
#     nsw_all      all NPWS fire (incl. unlabelled)     }
#     qld_wild     Qld scars minus QPWS planned burns   }
#     qld_planned  QPWS planned burns                   }
#     qld_all      Qld scars, all fire                  }
#   qld_years/qld_<y>.rds  per-year Qld cache (not in git)
# Existing tables are skipped; delete one to recompute it. Slow: NSW a few
# minutes per table, Qld a few minutes per year (recent years are large).
source("spatial/fire/fire_config.R")
for (f in c("grid", "sources", "cell_years", "regime")) source(file.path("spatial/fire/R", paste0(f, ".R")))
dir.create(FIRE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
cy_file  <- function(name) file.path(FIRE_DATA_DIR, paste0("cell_years_", name, ".csv"))
write_cy <- function(x, name) utils::write.csv(x, cy_file(name), row.names = FALSE)

cells  <- load_fire_cells()
states <- state_polygons()
message(sprintf("Cells: %d NSW, %d Qld", sum(cells$state == "NSW"), sum(cells$state == "Qld")))
nsw_cells <- cells[cells$state == "NSW", ]
qld_cells <- cells[cells$state == "Qld", ]

# ---- NSW: wildfire, prescribed, all ------------------------------------------------
nsw_jobs <- list(nsw_wild    = "wildfire",
                 nsw_planned = "planned",
                 nsw_all     = c("wildfire", "planned", "unknown"))
todo <- names(nsw_jobs)[!file.exists(vapply(names(nsw_jobs), cy_file, ""))]
if (!length(todo)) {
  message("NSW: cell-year tables exist; skipping")
} else {
  message("NSW: reading NPWS fire history")
  nsw <- read_nsw_npws(bbox = sf::st_bbox(nsw_cells))
  print(table(type = nsw$type, in_period = nsw$year %in% NSW_YEARS_ALL))
  for (nm in todo) {
    message(sprintf("NSW %s -> cells: %s", paste(nsw_jobs[[nm]], collapse = " + "), nm))
    write_cy(cell_year_fractions(nsw[nsw$type %in% nsw_jobs[[nm]], ], nsw_cells, NSW_YEARS_ALL,
                                 "land", state_poly = states[states$state == "NSW", ]), nm)
    invisible(gc())
  }
  rm(nsw); invisible(gc())
}

# ---- Queensland, one year at a time ------------------------------------------------
if (all(file.exists(cy_file(c("qld_wild", "qld_planned", "qld_all"))))) {
  message("Qld: cell-year tables exist; skipping")
} else {
  yr_dir <- file.path(FIRE_DATA_DIR, "qld_years"); dir.create(yr_dir, showWarnings = FALSE)
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
  write_cy(pick("frac_wild"), "qld_wild")
  write_cy(pick("frac_all"), "qld_all")
  write_cy(pick("frac_planned"), "qld_planned")
}
message("Done: ", FIRE_DATA_DIR, "  (next: spatial/02_cross_reference_layers.R)")
