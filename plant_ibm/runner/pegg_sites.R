# runner/pegg_sites.R
#
# Phase 2 fit: post-fire mortality of resprouting M. nodosa at several
# monitored plots (Pegg et al. 2025), by observation interval, fitting both
# the level and the spread between plots (docs/methodological_notes.md V6).
#
# WHY
#   One plot (Yarringully, 48% dead after two years) fixes the level only.
#   Seven plots with nearly the same modelled rust pressure range from 0 to
#   56% dead. In the current model a plot's outcome varies only through the
#   two post-fire annual rust draws, which barely vary at concentration 20,
#   so it never predicts a plot with no deaths (Nov 2026 check). Fitting the
#   spread asks how much between-site variation is needed, and whether it is
#   year-to-year (Beta concentration, notes C3) or persistent (rust_site_sd).
#
# OBSERVATION MODEL
#   Data (runner/data/pegg2025_mnodosa_plots.csv): for each plot, intervals
#   [t1, t2] in months since the fire, with n_at_risk trees alive at t1 and
#   the deaths by t2. Tagged trees are adult resprouters, so the model cohort
#   is plants of reproductive age (age >= their age at first flowering) alive
#   before the fire and not killed by it ('cohort = "mature"'; "all" = every
#   resprouter, as in the single-site validation V2).
#   Model: the Pegg scenario (spin-up 1910, rust 2010, one fire in model
#   year 2020 burning every plant) at the plot's grid cell; the cohort's
#   survival after each post-fire year S(12), S(24), S(36) months; within a
#   year the death rate is taken as constant (log-linear interpolation), so
#   any interval maps to a mortality p = 1 - S(t2)/S(t1).
#   Likelihood of a plot = mean over seeds (simulated sites) of the product
#   over its intervals of Binomial(deaths | n_at_risk, p_seed). The mixture
#   over seeds is what makes between-plot spread count: a parameter set that
#   never produces a near-zero site gives the zero-death plots almost no
#   probability. Total log-likelihood = sum over plots.
#
# COST
#   The spin-up to 2019 does not depend on resprout rust mortality
#   (max_effect acts only on resprouting plants, and the only fire is in
#   2020), so each (cell, concentration, site sd, seed) is spun up once and
#   branched for every max_effect value (3 post-fire years each). Results are
#   cached per task in results_dir.
#
# USAGE: runscripts/pegg.sites.fit.run.code

source("runner/site_params.R")

PEGG_PLOTS_PATH <- "runner/data/pegg2025_mnodosa_plots.csv"
PEGG_SITES_VERSION <- "v1-interval-mixture"

load_pegg_plots <- function(path = PEGG_PLOTS_PATH) {
  d <- utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE)
  d$plot <- trimws(d$plot)
  d
}

# Nearest occupied cell (row of 'sites') for each plot.
pegg_plot_cells <- function(plots, sites) {
  u <- unique(plots[, c("plot", "lon", "lat")])
  u$cell_row <- vapply(seq_len(nrow(u)), function(i)
    which.min((sites$lon_centre - u$lon[i])^2 + (sites$lat_centre - u$lat[i])^2), 1L)
  u$cell_id <- sites$cell_id[u$cell_row]
  u
}

# Parameters for one simulated site: the cell's values, the Pegg fire
# scenario, a Beta concentration, and a persistent site effect drawn from
# the seed (no draw when site_sd = 0).
pegg_site_scenario <- function(base, site_row, seed, conc = NULL, site_sd = 0,
                               year_spinup = 1910L, year_rust = 2010L, year_fire = 2020L) {
  p <- build_site_params(base, site_row)
  if (!is.null(conc)) p <- set_rust_concentration(p, conc)
  p$rust_site_sd <- site_sd
  z <- if (site_sd > 0) { set.seed(1000003L + seed); stats::rnorm(1, 0, site_sd) } else 0
  p <- apply_rust_site_effect(p, z)
  p$seed             <- seed
  p$rust_start_year  <- year_rust
  p$fire_years       <- as.integer(year_fire)
  p$fire_prob_annual <- 0
  p$fire_p_fimp      <- 1
  p$fire_size_sd     <- 0
  p$fire_stand_cover <- 1
  p
}

# Cohort survival S(12k months), k = 1..n_post, for each max_effect value.
pegg_cohort_survival <- function(p, max_effects, n_post = 3L, cohort = c("mature", "all"),
                                 year_spinup = 1910L, year_fire = 2020L) {
  cohort <- match.arg(cohort)
  p_spin <- p
  p_spin$n_years <- as.integer(year_fire - year_spinup)
  st <- run_simulation(p_spin, year0 = year_spinup, verbose = FALSE)
  pre <- st$individuals$alive
  if (cohort == "mature")
    pre <- pre & st$individuals$age >= st$individuals$age_first_flower
  ids <- st$individuals$id[pre]
  S <- t(vapply(max_effects, function(me) {
    q <- p
    q$rust_dose_response$resprout$max_effect <- me
    q$n_years <- as.integer(n_post)
    r <- run_simulation(q, init_state = st, year0 = year_fire, verbose = FALSE)
    ind <- r$individuals[match(ids, r$individuals$id), ]
    resp <- !(ind$death_cause %in% "fire" & !ind$alive)
    n <- sum(resp)
    vapply(seq_len(n_post), function(k) {
      gone <- !ind$alive[resp] & !is.na(ind$death_year[resp]) &
              ind$death_year[resp] <= year_fire + k - 1L
      1 - sum(gone) / n
    }, 0)
  }, numeric(n_post)))
  if (length(max_effects) == 1L) S <- matrix(S, nrow = 1L)
  dimnames(S) <- list(paste0("me_", max_effects), paste0("S", 12 * seq_len(n_post)))
  attr(S, "n_cohort") <- length(ids)
  S
}

# Mortality over [t1, t2] months from annual survivals (constant rate
# within each post-fire year).
interval_mortality <- function(S_annual, t1, t2) {
  logS <- log(pmax(c(1, S_annual), 1e-12))
  m <- 12 * (seq_along(logS) - 1L)
  if (t2 > max(m)) stop("interval ends after the simulated years; increase n_post")
  f <- stats::approxfun(m, logS)
  1 - exp(f(t2) - f(t1))
}

# ---- Running a grid -----------------------------------------------------------
# grid: data.frame(conc, site_sd); one task per grid row x cell x seed.
run_pegg_plot_grid <- function(grid, max_effects, seeds, results_dir,
                               base = get_default_params(), plots = load_pegg_plots(),
                               sites = load_site_params(drop_incomplete = TRUE),
                               cohort = "mature",
                               n_cores = max(1L, parallel::detectCores() - 1L)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  pc <- pegg_plot_cells(plots, sites)
  cells <- unique(pc$cell_row)
  tasks <- expand.grid(g = seq_len(nrow(grid)), cell = cells, seed = seeds)
  message(sprintf("Pegg plots: %d grid points x %d cells x %d seeds = %d spin-ups (x %d max_effect branches)",
                  nrow(grid), length(cells), length(seeds), nrow(tasks), length(max_effects)))
  parallel::mclapply(seq_len(nrow(tasks)), function(j) {
    tk <- tasks[j, ]
    g <- grid[tk$g, ]
    p <- pegg_site_scenario(base, sites[tk$cell, ], tk$seed, conc = g$conc, site_sd = g$site_sd)
    key <- make_run_id(list(p = p, me = max_effects, cohort = cohort, v = PEGG_SITES_VERSION))
    f <- file.path(results_dir, paste0(key, ".rds"))
    if (file.exists(f)) return(invisible(NULL))
    S <- tryCatch(pegg_cohort_survival(p, max_effects, cohort = cohort),
                  error = function(e) { message("task ", j, ": ", conditionMessage(e)); NULL })
    if (!is.null(S))
      saveRDS(list(conc = g$conc, site_sd = g$site_sd, cell_id = sites$cell_id[tk$cell],
                   seed = tk$seed, z = p$rust_site_z, max_effects = max_effects,
                   S = S, cohort = cohort, version = PEGG_SITES_VERSION), f)
    invisible(NULL)
  }, mc.cores = n_cores)
  invisible(pc)
}

# Read cached tasks into one long table: conc, site_sd, cell_id, seed, z,
# max_effect, S12, S24, S36.
collect_pegg_plot_grid <- function(results_dir) {
  fs <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(fs, function(f) {
    x <- readRDS(f)
    if (!identical(x$version, PEGG_SITES_VERSION)) return(NULL)
    data.frame(conc = x$conc, site_sd = x$site_sd, cell_id = x$cell_id, seed = x$seed,
               z = x$z, max_effect = x$max_effects, x$S, row.names = NULL)
  })
  do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
}

# Log-likelihood of the plot data for each (conc, site_sd, max_effect), with
# per-plot predictions (mean deaths, 5-95% range, P(no deaths)).
pegg_plot_loglik <- function(sims, plots = load_pegg_plots(),
                             sites = load_site_params(drop_incomplete = TRUE)) {
  pc <- pegg_plot_cells(plots, sites)
  plots$cell_id <- pc$cell_id[match(plots$plot, pc$plot)]
  Scols <- grep("^S[0-9]+$", names(sims), value = TRUE)
  keys <- unique(sims[, c("conc", "site_sd", "max_effect")])
  per_plot <- list(); tot <- numeric(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, ]
    sk <- sims[sims$conc == k$conc & sims$site_sd == k$site_sd & sims$max_effect == k$max_effect, ]
    ll <- 0
    for (pl in unique(plots$plot)) {
      obs <- plots[plots$plot == pl, ]
      ss <- sk[sk$cell_id == obs$cell_id[1], ]
      if (!nrow(ss)) next
      # likelihood of this plot under each simulated site (seed)
      Ls <- vapply(seq_len(nrow(ss)), function(r) {
        S <- unlist(ss[r, Scols])
        prod(vapply(seq_len(nrow(obs)), function(j)
          stats::dbinom(obs$deaths[j], obs$n_at_risk[j],
                        interval_mortality(S, obs$t1_months[j], obs$t2_months[j])), 0))
      }, 0)
      ll <- ll + log(mean(Ls) + 1e-300)
      # predicted deaths over the plot's whole record, per simulated site
      pr <- vapply(seq_len(nrow(ss)), function(r) {
        S <- unlist(ss[r, Scols])
        1 - prod(1 - vapply(seq_len(nrow(obs)), function(j)
          interval_mortality(S, obs$t1_months[j], obs$t2_months[j]), 0))
      }, 0)
      n0 <- obs$n_at_risk[1]
      draws <- stats::rbinom(4000, n0, sample(pr, 4000, replace = TRUE))
      per_plot[[length(per_plot) + 1L]] <- data.frame(
        conc = k$conc, site_sd = k$site_sd, max_effect = k$max_effect, plot = pl,
        observed = sum(obs$deaths) / n0, pred_mean = mean(draws) / n0,
        pred_lo = stats::quantile(draws, .05) / n0, pred_hi = stats::quantile(draws, .95) / n0,
        p_zero = mean(draws == 0), n_sims = nrow(ss), row.names = NULL)
    }
    tot[i] <- ll
  }
  keys$loglik <- tot
  keys <- keys[order(-keys$loglik), ]
  list(loglik = keys, per_plot = do.call(rbind, per_plot))
}

# Log-likelihood surface (one panel per variance structure) and, for the
# best parameter set, predicted vs observed mortality per plot.
plot_pegg_plot_fit <- function(fit, save_path = NULL) {
  ll <- fit$loglik
  ll$structure <- ifelse(ll$site_sd > 0, sprintf("persistent site effect, sd %.2f", ll$site_sd),
                         sprintf("year-to-year only, conc %g", ll$conc))
  best <- ll[1, ]
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  pA <- gg(ll, aes(max_effect, loglik, colour = structure)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
    ggplot2::labs(x = "Extra rust mortality of resprouting plants (max_effect)",
                  y = "Log-likelihood (7 plots)", colour = NULL,
                  title = "Fit to post-fire mortality at the Pegg et al. plots") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  pp <- fit$per_plot[fit$per_plot$conc == best$conc & fit$per_plot$site_sd == best$site_sd &
                     fit$per_plot$max_effect == best$max_effect, ]
  pB <- gg(pp, aes(stats::reorder(plot, observed))) +
    ggplot2::geom_linerange(aes(ymin = pred_lo, ymax = pred_hi), colour = "#1d5fae", linewidth = 3, alpha = 0.4) +
    ggplot2::geom_point(aes(y = pred_mean), colour = "#1d5fae", size = 2) +
    ggplot2::geom_point(aes(y = observed), colour = "#c0392b", size = 3, shape = 18) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Share of tagged trees dead over the plot's record",
                  title = sprintf("Best: conc %g, site sd %.2f, max_effect %.2f",
                                  best$conc, best$site_sd, best$max_effect),
                  subtitle = "Blue: model, mean and 5-95% for a 50-tree plot; red: observed") +
    ggplot2::theme_minimal(base_size = 11)
  p <- patchwork::wrap_plots(pA, pB, nrow = 1, widths = c(1.2, 1))
  if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 13, height = 5.5, dpi = 150, bg = "white")
  p
}
