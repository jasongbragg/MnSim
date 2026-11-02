# runner/range_ensemble.R
#
# Range-wide IUCN A4 ensemble for M. nodosa (replaces spatial/03_range_wide_iucn.R).
#
# DESIGN (docs/methodological_notes.md B1-B3)
#   run i has its own random-number stream (master_seed + i), so adding runs
#   never changes earlier ones. Each run draws
#     - a grid cell, with probability proportional to abundance (weight_n_raw);
#       the cell's fire regime and rust pressure (runner/site_params.R)
#     - a simulation seed (fire timing, annual rust years, demography)
#     - values for the uncertain parameters in 'param_draws'
#   and gives one A4 decline. A plan table (range_plan.csv) records every run's
#   cell, seed, draws and run_id; each run's result is cached as <run_id>.rds.
#
# USAGE (from plant_ibm/, after sourcing params.R, R/*.R, runner/params_utils.R,
#        runner/iucn_a4.R):
#   source("runner/range_ensemble.R")
#   run_range_ensemble(n_runs = 300, results_dir = "results/range_v1/", n_cores = 7)
#   R <- collect_range_ensemble("results/range_v1/")
#   range_category_shares(R$runs)            # overall and by region
#   plot_range_drivers(R, "results/range_v1/drivers.png")
# Figures for talks: talk/make_talk_figures.R (uses these functions).

source("runner/site_params.R")

# ---- Simulation phasing (calendar years) -----------------------------------
YEAR_SPINUP        <- 1600L
YEAR_GENTIME_START <- 1810L
YEAR_GENTIME_END   <- 1910L
YEAR_RUST          <- 2010L
YEAR_ASSESS        <- 2026L   # the present
YEAR_END           <- 2210L   # >= 2009 + 3*T_g; comparison capped at YEAR_ASSESS + 100
AVG_YEARS          <- 11L     # reference and comparison averaging windows

ASSESS <- iucn_assess_spec(
  year_spinup = YEAR_SPINUP, year_gentime_start = YEAR_GENTIME_START,
  year_gentime_end = YEAR_GENTIME_END, year_rust = YEAR_RUST,
  year_assess = YEAR_ASSESS, year_end = YEAR_END, avg_years = AVG_YEARS)

QLD_LAT <- -28.16   # approx. NSW/Qld border on the coast

# ---- Parameter uncertainty --------------------------------------------------
# Drawn independently per run, uniform on [min, max]; 'path' locates the value
# in the params list. list() switches parameter uncertainty off.
DEFAULT_PARAM_DRAWS <- list(
  resprout_rust_mortality = list(
    path  = c("rust_dose_response", "resprout", "max_effect"),
    min   = 0.32, max = 0.78,
    label = "Extra rust mortality of resprouting plants (max_effect)",
    note  = "values whose Pegg-scenario median lies in the 95% CI of 24/50 dead (0.34-0.63), at the Yarringully cell")
)

# Apply one run's parameter draws (uses the current RNG stream).
draw_params <- function(p, param_draws = DEFAULT_PARAM_DRAWS) {
  vals <- list()
  for (nm in names(param_draws)) {
    d <- param_draws[[nm]]
    v <- stats::runif(1, d$min, d$max)
    p[[d$path]] <- v
    vals[[nm]]  <- v
  }
  list(params = p, values = vals)
}

draws_caption <- function(param_draws = DEFAULT_PARAM_DRAWS) {
  if (length(param_draws) == 0L) return("No parameter draws (site and seed variation only).")
  paste(vapply(param_draws, function(d)
    sprintf("%s ~ U(%.2f, %.2f) [%s]", d$label, d$min, d$max, d$note), ""),
    collapse = "; ")
}

# Cells usable for runs (complete rust and fire values), with region.
load_range_sites <- function(path = SITE_PARAMS_PATH) {
  s <- load_site_params(path, drop_incomplete = TRUE)
  s$region <- ifelse(s$lat_centre > QLD_LAT, "Qld", "NSW")
  s
}

# ---- Planning ---------------------------------------------------------------
plan_range_run <- function(i, sites, base, param_draws = DEFAULT_PARAM_DRAWS,
                           master_seed = 20261022L) {
  set.seed(master_seed + i)
  k  <- sample.int(nrow(sites), 1L, prob = sites$weight_n_raw)
  sd <- sample.int(99999L, 1L)
  dr <- draw_params(base, param_draws)
  p  <- build_site_params(dr$params, sites[k, ])
  p$seed <- sd
  row <- data.frame(run = i, cell_row = k, cell_id = sites$cell_id[k],
                    lon = sites$lon_centre[k], lat = sites$lat_centre[k],
                    region = sites$region[k], sim_seed = sd,
                    site_fire_prob = sites$fire_prob_annual[k],
                    site_fire_p_fimp = sites$fire_p_fimp[k],
                    site_rust = sites$annual_env_alpha[k] /
                                (sites$annual_env_alpha[k] + sites$annual_env_beta[k]),
                    run_id = iucn_run_id(p, ASSESS))
  for (nm in names(dr$values)) row[[nm]] <- dr$values[[nm]]
  list(params = p, row = row)
}

# ---- Running ----------------------------------------------------------------
run_range_ensemble <- function(n_runs, results_dir,
                               n_cores     = max(1L, parallel::detectCores() - 1L),
                               param_draws = DEFAULT_PARAM_DRAWS,
                               master_seed = 20261022L,
                               base        = get_default_params(),
                               sites       = load_range_sites()) {
  dir.create(file.path(results_dir, "range"), recursive = TRUE, showWarnings = FALSE)
  base$fire_years <- integer(0)
  message(sprintf(
    "Range-wide: %d runs on %d cores (%d cells; %d cells without rust/fire values dropped, %.1f%% of abundance weight)",
    n_runs, n_cores, nrow(sites), attr(sites, "n_dropped"), 100 * attr(sites, "weight_dropped")))
  plans <- lapply(seq_len(n_runs), plan_range_run, sites = sites, base = base,
                  param_draws = param_draws, master_seed = master_seed)
  plan  <- do.call(rbind, lapply(plans, `[[`, "row"))
  write.csv(plan, file.path(results_dir, "range_plan.csv"), row.names = FALSE)

  rdir <- file.path(results_dir, "range")
  parallel::mclapply(seq_len(n_runs), function(i) {
    out <- file.path(rdir, paste0(plans[[i]]$row$run_id, ".rds"))
    if (file.exists(out)) return(invisible(NULL))
    tryCatch(
      suppressWarnings(run_iucn_a4(
        plans[[i]]$params,
        year_spinup = YEAR_SPINUP, year_gentime_start = YEAR_GENTIME_START,
        year_gentime_end = YEAR_GENTIME_END, year_rust = YEAR_RUST,
        year_assess = YEAR_ASSESS, year_end = YEAR_END, avg_years = AVG_YEARS,
        save_rds = out, plot = FALSE, verbose = FALSE)),
      error = function(e) message("  range run ", i, ": ", conditionMessage(e)))
    invisible(NULL)
  }, mc.cores = n_cores)
  invisible(plan)
}

# ---- Collecting ---------------------------------------------------------------
# Summaries (one row per finished run, joined to the plan) and N_IUCN / N_ref
# trajectories over 'traj_years'.
collect_range_ensemble <- function(results_dir, n_runs = NULL, traj_years = 1950:2130) {
  plan <- read.csv(file.path(results_dir, "range_plan.csv"))
  if (!is.null(n_runs)) plan <- plan[plan$run <= n_runs, ]
  rdir <- file.path(results_dir, "range")
  rows <- list(); traj <- list()
  for (j in seq_len(nrow(plan))) {
    f <- file.path(rdir, paste0(plan$run_id[j], ".rds"))
    if (!file.exists(f)) next
    r <- readRDS(f)
    s <- summarise_iucn_a4(r)
    rows[[length(rows) + 1L]] <- cbind(plan[j, setdiff(names(plan), "run_id")], s)
    cen <- r$census
    N <- cen$N_IUCN[match(traj_years, cen$year)]
    if (isTRUE(cen$extinct[nrow(cen)])) N[traj_years > max(cen$year)] <- 0
    traj[[length(traj) + 1L]] <- data.frame(run = plan$run[j], year = traj_years,
                                            rel = N / s$N_ref)
  }
  if (length(rows) == 0L) stop("No range-wide results found in ", rdir)
  df <- do.call(rbind, rows)
  df$category <- factor(vapply(df$pct_decline, iucn_category, ""),
                        levels = c("LC/NT", "VU", "EN", "CR"))
  message(sprintf("Collected %d of %d planned runs.", nrow(df), nrow(plan)))
  list(runs = df, traj = do.call(rbind, traj))
}

# Share of runs per IUCN category, overall and by region (rows: groups).
range_category_shares <- function(df) {
  one <- function(x) {
    t <- table(factor(vapply(x, iucn_category, ""), levels = c("LC/NT", "VU", "EN", "CR")))
    c(n = length(x), round(100 * as.numeric(t) / max(1, length(x)), 1),
      median = round(stats::median(x, na.rm = TRUE), 1))
  }
  out <- rbind(all = one(df$pct_decline),
               do.call(rbind, lapply(split(df$pct_decline, df$region), one)))
  colnames(out) <- c("n", "LC/NT", "VU", "EN", "CR", "median_decline")
  out
}

# ---- Diagnostic plot: what drives the decline? -------------------------------
# A: decline vs cell fire frequency, coloured by rust pressure.
# B: fire x rust surface, coloured by decline (carried over from 03).
# C: decline vs the drawn parameter(s).
plot_range_drivers <- function(R, save_path = NULL, param_draws = DEFAULT_PARAM_DRAWS) {
  df <- R$runs
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  th <- ggplot2::theme_minimal(base_size = 12)
  hl <- ggplot2::geom_hline(yintercept = c(30, 50, 80), colour = "#9a9a9a", linetype = "dashed")
  pA <- gg(df, aes(1 / site_fire_prob, pct_decline, colour = site_rust)) +
    hl + ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_viridis_c(option = "rocket", direction = -1, name = "Rust\npressure") +
    ggplot2::labs(x = "Mean years between fires (cell)", y = "Decline (%)",
                  title = "Decline vs fire frequency") + th
  pB <- gg(df, aes(1 / site_fire_prob, site_rust, colour = pct_decline)) +
    ggplot2::geom_jitter(alpha = 0.7, size = 1.8, width = 0.02, height = 0.005) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_gradient2(low = "#2a78d6", mid = "#fdae6b", high = "#d6452f",
                                    midpoint = 50, name = "Decline (%)") +
    ggplot2::labs(x = "Mean years between fires (cell)", y = "Mean rust pressure (cell)",
                  title = "Fire x rust") + th
  panels <- list(pA, pB)
  for (nm in intersect(names(param_draws), names(df))) {
    panels[[length(panels) + 1L]] <- gg(df, aes(.data[[nm]], pct_decline)) +
      hl + ggplot2::geom_point(alpha = 0.6, size = 1.6, colour = "#1d5fae") +
      { if (nrow(df) >= 30)
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                               colour = "#2b2b2b", linewidth = 0.8) } +
      ggplot2::labs(x = param_draws[[nm]]$label, y = "Decline (%)",
                    title = "Decline vs drawn parameter") + th
  }
  p <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title = sprintf("What drives the decline? (%d runs)", nrow(df)),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  if (!is.null(save_path))
    ggplot2::ggsave(save_path, p, width = 5.2 * length(panels), height = 5, dpi = 170, bg = "white")
  invisible(p)
}
