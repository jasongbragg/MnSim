# spatial/fire/03_fire_checks.R -- the evidence behind the fire choices
# (docs/methodological_notes.md D1-D5). Run from plant_ibm/ after
# 02_fire_cell_years.R; not needed to build mn_site_params.csv.
#
#   1. every regime and sensitivity scenario: fire years, fitted prior, and a
#      test of one shared prior against one per state
#   2. border check: do NSW and Qld agree in cells within +/- 1 degree of
#      the border? (figures/fire/border_<scenario>.png)
#   3. Qld record year by year: the step after 2013, and the "chronic" cells
#      marked burnt in (nearly) every automated year (figures/fire/qld_by_year.png)
#   4. wildfire vs planned burns, per state
# Writes spatial/data/fire/fire_scenarios_summary.csv, fire_priors.csv,
# border_summary.csv and qld_chronic_cells.csv.
source("spatial/fire/fire_config.R")
for (f in c("grid", "estimate", "regime", "border"))
  source(file.path("spatial/fire/R", paste0(f, ".R")))
suppressPackageStartupMessages(library(ggplot2))
dir.create(FIRE_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
cells <- load_fire_cells()

# ---- 1. Scenarios -----------------------------------------------------------------
scen <- c(FIRE_REGIMES, FIRE_SENSITIVITY)
regs <- lapply(scen, run_fire_scenario, cells = cells)
summ <- do.call(rbind, Map(function(name, reg) {
  lr <- prior_lr_test(reg$n_wildfire_years_sig, reg$T_obs, reg$state)
  s <- do.call(rbind, lapply(split(reg, reg$state), function(d) data.frame(
    scenario = name, state = d$state[1], cells = nrow(d), T = d$T_obs[1],
    fire_years_per_100_cell_years = round(100 * sum(d$n_wildfire_years_sig) / sum(d$T_obs), 2),
    share_cells_with_fire = round(mean(d$n_wildfire_years_sig > 0), 2),
    median_FRI = stats::median(d$FRI_years),
    median_frac = round(stats::median(d$mean_frac_burned, na.rm = TRUE), 3))))
  s$same_prior_p <- if (is.null(lr)) NA else signif(lr[["p"]], 2)
  s
}, names(regs), regs))
pri <- do.call(rbind, Map(function(name, reg) { p <- attr(reg, "priors"); p$scenario <- name; p },
                          names(regs), regs))
utils::write.csv(summ, file.path(FIRE_DATA_DIR, "fire_scenarios_summary.csv"), row.names = FALSE)
utils::write.csv(pri, file.path(FIRE_DATA_DIR, "fire_priors.csv"), row.names = FALSE)
cat("Fitted priors (Gamma(alpha, beta) on the annual fire-year rate):\n")
print(transform(pri[, c("scenario", "group", "alpha", "beta", "type")],
                alpha = round(alpha, 2), beta = round(beta, 1),
                mean_FRI = round(pri$beta / pri$alpha)), row.names = FALSE)
cat("\nScenarios (same_prior_p: test of one shared prior vs one per state):\n")
print(summ, row.names = FALSE)

# ---- 2. Border check --------------------------------------------------------------
bs <- list()
for (name in names(regs)) {
  sc <- scen[[name]]
  bc <- border_compare(regs[[name]], cells, attr(regs[[name]], "cy"), list(NSW = sc$nsw, Qld = sc$qld))
  plot_border(bc, name, file.path(FIRE_FIG_DIR, sprintf("border_%s.png", name)))
  s <- as.data.frame(bc$summary); s$scenario <- name
  if (!is.null(bc$rate_ratio)) {
    s$qld_nsw_rate_ratio <- unname(bc$rate_ratio$estimate)
    s$rr_lo <- bc$rate_ratio$conf.int[1]; s$rr_hi <- bc$rate_ratio$conf.int[2]
  }
  bs[[name]] <- s
}
bs <- do.call(rbind, bs)
utils::write.csv(bs, file.path(FIRE_DATA_DIR, "border_summary.csv"), row.names = FALSE)
cat(sprintf("\nBorder band (+/- %.1f deg): Qld/NSW fire-year rate\n", BORDER_BAND))
rr <- unique(bs[, c("scenario", "qld_nsw_rate_ratio", "rr_lo", "rr_hi")])
print(transform(rr, qld_nsw_rate_ratio = round(qld_nsw_rate_ratio, 2),
                rr_lo = round(rr_lo, 2), rr_hi = round(rr_hi, 2)), row.names = FALSE)

# ---- 3. Queensland record year by year ---------------------------------------------
sig <- function(d) d[d$frac_burned >= AREA_THRESH, ]
nsw <- sig(read_cy("nsw_wild")); qld <- sig(read_cy("qld_wild"))
n_nsw <- sum(cells$state == "NSW"); n_qld <- sum(cells$state == "Qld")
share <- function(d, yrs, n) vapply(yrs, function(y) length(unique(d$cell_id[d$year == y])) / n, 0)
yrs <- QLD_YEARS_ALL
by_year <- rbind(data.frame(state = "NSW", year = yrs, share = share(nsw, yrs, n_nsw)),
                 data.frame(state = "Qld", year = yrs, share = share(qld, yrs, n_qld)))
p <- ggplot(by_year, aes(year, share, colour = state)) + geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = min(QLD_YEARS) - 0.5, linetype = "dotted") +
  scale_colour_manual(values = c(NSW = "#1d5fae", Qld = "#c0392b")) +
  labs(x = NULL, y = "Share of cells with a fire year (>= 10% burnt)",
       title = "Wildfire, all occupied cells, year by year",
       subtitle = sprintf("Qld used from %d (dotted line); earlier years are automated mapping",
                          min(QLD_YEARS))) +
  theme_minimal(base_size = 12)
ggsave(file.path(FIRE_FIG_DIR, "qld_by_year.png"), p, width = 9, height = 4.5, dpi = 160, bg = "white")

cd <- sf::st_drop_geometry(cells)
q <- cd[cd$state == "Qld", c("cell_id", "lon_centre", "lat_centre", "land_frac")]
cnt <- function(d, yy) as.integer(table(factor(d$cell_id[d$year %in% yy], levels = q$cell_id)))
auto <- QLD_YEARS_ALL[QLD_YEARS_ALL < min(QLD_YEARS) - 1L]
q$k_automated <- cnt(qld, auto)
q$fire_last_automated_year <- cnt(qld, min(QLD_YEARS) - 1L)
q$k_used <- cnt(qld, QLD_YEARS)
q$chronic <- q$k_automated >= round(0.75 * length(auto))
utils::write.csv(q, file.path(FIRE_DATA_DIR, "qld_chronic_cells.csv"), row.names = FALSE)
cat(sprintf("\nQld: %d of %d cells marked burnt in >= 75%% of the automated years %d-%d.\n",
            sum(q$chronic), nrow(q), min(auto), max(auto)))
cat(sprintf("  chronic cells: fire in %d %.0f%%, fire years %d-%d per cell %.2f\n",
            min(QLD_YEARS) - 1L, 100 * mean(q$fire_last_automated_year[q$chronic] > 0),
            min(QLD_YEARS), max(QLD_YEARS), mean(q$k_used[q$chronic])))
cat(sprintf("  other cells:   fire in %d %.0f%%, fire years %d-%d per cell %.2f\n",
            min(QLD_YEARS) - 1L, 100 * mean(q$fire_last_automated_year[!q$chronic] > 0),
            min(QLD_YEARS), max(QLD_YEARS), mean(q$k_used[!q$chronic])))

# ---- 4. Wildfire vs planned burns ---------------------------------------------------
cm <- data.frame(cell_id = regs$wild$cell_id, state = regs$wild$state,
                 wild = regs$wild$n_wildfire_years_sig,
                 planned = regs$planned$n_wildfire_years_sig[match(regs$wild$cell_id, regs$planned$cell_id)],
                 all = regs$all$n_wildfire_years_sig[match(regs$wild$cell_id, regs$all$cell_id)],
                 frac_wild = regs$wild$mean_frac_burned,
                 frac_planned = regs$planned$mean_frac_burned[match(regs$wild$cell_id, regs$planned$cell_id)])
cat("\nFire years by type (main period). Wild + planned can exceed all: both clear\n",
    "the 10% threshold in some cell-years.\n", sep = "")
print(do.call(rbind, lapply(split(cm, cm$state), function(d) data.frame(
  state = d$state[1], all = sum(d$all), wild = sum(d$wild), planned = sum(d$planned),
  planned_share = round(sum(d$planned) / sum(d$all), 2),
  cells_with_planned = sum(d$planned > 0),
  median_frac_wild = round(stats::median(d$frac_wild, na.rm = TRUE), 3),
  median_frac_planned = round(stats::median(d$frac_planned, na.rm = TRUE), 3)))), row.names = FALSE)
message("\nFigures: ", FIRE_FIG_DIR)
