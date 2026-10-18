# runner/pegg_validation.R
#
# Validation of the IBM against Pegg et al. (2025), who followed resprouting
# Melaleuca nodosa after the Black Summer fires (Dec 2019 - Jan 2020) at
# sites in northern New South Wales. At the Jan 2022 follow-up (~2 yr
# post-fire), 48% of the trees that had resprouted were dead.
#
# TIMING (model year t = one annual step, Jan..Dec of t)
# ------------------------------------------------------
# Within a model year the order is: age -> resprout recovery -> FIRE ->
# flowering -> MORTALITY -> recruitment -> census. So:
#   year_fire = 2020   the fire step of model year 2020 is the Jan 2020 fire;
#                      the mortality step of 2020 is the first post-fire year.
#   year_survey = 2022 the Jan 2022 survey = state at the END of model year
#                      2021 (= year_survey - 1). Post-fire deaths therefore
#                      come from two mortality steps (2020 and 2021), matching
#                      the ~2 yr between fire and survey.
# The run stops at the end of year_survey - 1, so each plant's final state
# IS its state at the survey.
#
# SIMULATION DESIGN
# -----------------
# Spin-up (no fire): year_spinup -> year_rust - 1  (default 1910-2009)
# Rust period:       year_rust  ->                 (default from 2010)
# Fire:              year_fire, single deterministic event, fire_p_fimp = 1
#                    (every tree at the survey sites burned)
# End of run:        year_survey - 1               (default 2021)
#
# PRIMARY STATISTIC
# -----------------
# Fraction of resprouting fire survivors dead at the survey:
#
#   fraction_resprout_dead = postfire_dead / N_resprouted
#
#   N_prefire      plants alive at the end of year_fire - 1
#   fire_dead      pre-fire plants with death_cause == "fire" (killed in
#                  the fire itself)
#   N_resprouted   N_prefire - fire_dead (fire_p_fimp = 1, so every
#                  surviving plant was impacted and resprouted)
#   postfire_dead  pre-fire plants with death_cause == "mortality" and
#                  death_year in year_fire .. year_survey - 1
#
# Using death_cause matters: deaths from rust/background mortality in the
# fire year itself are post-fire deaths of resprouters, not fire deaths.
# (Earlier versions classified on death_year == year_fire, which moved those
# deaths out of both numerator and denominator and biased the statistic low.)
#
# OBSERVATION (Pegg et al. 2025, Yarringully NR)
# ----------------------------------------------
# 50 fire-affected trees, all resprouting when monitoring began (Mar 2020),
# assessed to Jul 2021 with a follow-up in Jan 2022: 24 of 50 (48%) dead.
# The 95% binomial interval for 24/50 is 0.34-0.63; plot_pegg_ensemble()
# shades it by default (observed_n = 50). The target is to land roughly in
# that range after two years, not to reproduce the monthly curve.
#
# FATES (mutually exclusive, sum to N_prefire)
# --------------------------------------------
#   (1) killed by the fire
#   (2) resprouted, then died before the survey
#   (3) alive and still resprouting at the survey
#   (4) alive and recovered (resprout = FALSE) at the survey
#
# USAGE
# -----
#   source all model files, then:
#   source("runner/params_utils.R")
#   source("runner/pegg_validation.R")
#
#   result <- run_pegg_validation(get_default_params())          # single run
#   df <- run_pegg_seeds(get_default_params(), seeds = 1:100,
#                        results_dir = "results/pegg_seeds_v2/")
#   plot_pegg_ensemble(df, observed = 0.48)
#
# Seed results are cached by pegg_run_id(), which hashes the params, the
# scenario years and PEGG_VALIDATION_VERSION, so results from earlier
# definitions are never silently reused.

PEGG_VALIDATION_VERSION <- "v2-death-cause-survey-timing"

pegg_run_id <- function(params, scenario) {
  make_run_id(list(params = params, scenario = scenario,
                   version = PEGG_VALIDATION_VERSION))
}

# ============================================================================
# compute_pegg_stats()
# ============================================================================
# Primary statistic and fate breakdown from a run that ENDED at the survey
# (end of year_survey - 1). year_spinup is needed to interpret founder
# birth_year encoding.
compute_pegg_stats <- function(pop, census, year_fire, year_survey, year_spinup) {

  if (is.null(pop$death_cause))
    stop("pop has no death_cause column; rerun with the current model code.")
  last_year <- year_survey - 1L

  # Calendar birth year (founders have birth_year = -(age at spinup))
  cal_birth <- ifelse(pop$birth_year < year_spinup,
                       year_spinup + pop$birth_year,
                       pop$birth_year)

  # Pre-fire cohort: alive at the end of year_fire - 1
  prefire_year  <- year_fire - 1L
  existed       <- cal_birth <= prefire_year
  alive_prefire <- existed & (pop$alive |
                               (!is.na(pop$death_year) & pop$death_year >= year_fire))
  N_prefire <- sum(alive_prefire)

  if (N_prefire == 0L) {
    warning("No pre-fire individuals found. Check year_fire and year_spinup.")
    return(list(N_prefire = 0L, fraction_dead = NA_real_,
                fraction_resprout_dead = NA_real_))
  }

  pf <- pop[alive_prefire, ]
  if (any(!is.na(pf$death_year) & pf$death_year > last_year))
    stop("Run continues past the survey; end it at year_survey - 1.")

  fire_dead       <- !pf$alive & pf$death_cause %in% "fire"
  post_fire_dead  <- !pf$alive & pf$death_cause %in% "mortality"
  resprouting_now <- pf$alive & pf$resprout
  recovered       <- pf$alive & !pf$resprout

  n_fire_dead     <- sum(fire_dead)
  n_postfire_dead <- sum(post_fire_dead)
  n_resprouting   <- sum(resprouting_now)
  n_recovered     <- sum(recovered)
  if (n_fire_dead + n_postfire_dead + n_resprouting + n_recovered != N_prefire)
    stop("Pegg fate categories do not partition the pre-fire cohort.")

  N_resprouted           <- N_prefire - n_fire_dead
  fraction_resprout_dead <- if (N_resprouted > 0L)
    n_postfire_dead / N_resprouted else NA_real_


  N_iucn_prefire <- census$N_IUCN[census$year == prefire_year]
  N_iucn_survey  <- census$N_IUCN[census$year == last_year]
  if (length(N_iucn_prefire) == 0L) N_iucn_prefire <- NA_real_
  if (length(N_iucn_survey)  == 0L) N_iucn_survey  <- NA_real_

  list(
    N_prefire              = N_prefire,
    N_resprouted           = N_resprouted,
    n_fire_dead            = n_fire_dead,
    n_postfire_dead        = n_postfire_dead,
    n_resprouting          = n_resprouting,
    n_recovered            = n_recovered,
    fraction_dead          = (n_fire_dead + n_postfire_dead) / N_prefire,
    fraction_fire          = n_fire_dead     / N_prefire,
    fraction_postfire      = n_postfire_dead / N_prefire,
    fraction_resprout      = n_resprouting   / N_prefire,
    fraction_recover       = n_recovered     / N_prefire,
    fraction_resprout_dead = fraction_resprout_dead,   # PRIMARY: of resprouters
    N_iucn_prefire         = N_iucn_prefire,
    N_iucn_assess          = N_iucn_survey
  )
}

# ============================================================================
# run_pegg_validation()
# ============================================================================
run_pegg_validation <- function(
    params,
    year_spinup = 1910L,
    year_rust   = 2010L,
    year_fire   = 2020L,
    year_survey = 2022L,   # Jan 2022 survey = end of model year 2021
    observed    = 0.48,    # Pegg et al.: 48% of resprouting trees dead by Jan 2022
    save_rds    = NULL,
    plot        = TRUE,
    verbose     = TRUE
) {
  scenario <- list(year_spinup = as.integer(year_spinup),
                   year_rust   = as.integer(year_rust),
                   year_fire   = as.integer(year_fire),
                   year_survey = as.integer(year_survey))

  p                  <- params
  p$rust_start_year  <- year_rust
  p$fire_years       <- as.integer(year_fire)
  p$fire_prob_annual <- 0
  # Run ends at the end of year_survey - 1 (last t = year0 + n_years - 1).
  p$n_years          <- as.integer(year_survey - year_spinup)
  # All trees at the survey sites burned: fire_p_fimp = 1 (not a free parameter).
  p$fire_p_fimp      <- 1.0

  if (verbose) message(sprintf(
    "Pegg validation: spin=%d  rust=%d  fire=%d  survey=Jan %d (end of %d)  p_fimp=1",
    year_spinup, year_rust, year_fire, year_survey, year_survey - 1L))

  res <- run_simulation(p, year0 = year_spinup, verbose = FALSE)

  stats <- compute_pegg_stats(res$individuals, res$census,
                               year_fire, year_survey, year_spinup)

  if (verbose) message(sprintf(
    "  fraction_resprout_dead=%.3f (observed=%.2f)  |  N_resprouted=%d  fire_dead=%.3f",
    stats$fraction_resprout_dead, observed,
    stats$N_resprouted, stats$fraction_fire))

  result <- c(
    list(
      year_spinup     = year_spinup,
      year_rust       = year_rust,
      year_fire       = year_fire,
      year_survey     = year_survey,
      year_assessment = year_survey - 1L,   # last model year (for plots)
      observed        = observed,
      census          = res$census,
      fire_p_fimp     = p$fire_p_fimp,
      version         = PEGG_VALIDATION_VERSION,
      run_id          = pegg_run_id(params, scenario),
      params          = params
    ),
    stats
  )

  if (!is.null(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, save_rds)
  }

  if (plot) plot_pegg_validation(result)
  invisible(result)
}

# ============================================================================
# summarise_pegg()
# ============================================================================
# Single-row data.frame for batch collection (same pattern as summarise_iucn_a4).
summarise_pegg <- function(result) {
  data.frame(
    run_id                 = result$run_id,
    fire_p_fimp            = result$fire_p_fimp,
    N_prefire              = result$N_prefire,
    N_resprouted           = result$N_resprouted,
    fraction_resprout_dead = result$fraction_resprout_dead,  # PRIMARY stat
    fraction_fire          = result$fraction_fire,
    fraction_postfire      = result$fraction_postfire,
    fraction_resprout      = result$fraction_resprout,
    fraction_recover       = result$fraction_recover,
    N_iucn_prefire         = result$N_iucn_prefire,
    N_iucn_assess          = result$N_iucn_assess,
    observed               = result$observed,
    residual               = result$fraction_resprout_dead - result$observed,
    stringsAsFactors       = FALSE
  )
}

# ============================================================================
# collect_pegg_results()
# ============================================================================
collect_pegg_results <- function(results_dirs) {
  if (is.character(results_dirs)) results_dirs <- as.list(results_dirs)
  files <- unlist(lapply(results_dirs, function(d)
    list.files(d, pattern = "\\.rds$", full.names = TRUE)))
  if (length(files) == 0L) { warning("No .rds files found."); return(data.frame()) }
  rows <- lapply(files, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || is.null(x$fraction_dead)) return(NULL)
    if (!identical(x$version, PEGG_VALIDATION_VERSION)) return(NULL)
    summarise_pegg(x)
  })
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0L) {
    warning("No results from the current Pegg validation version found.")
    return(data.frame())
  }
  out  <- do.call(rbind, rows)
  out  <- out[!duplicated(out$run_id), ]
  rownames(out) <- NULL
  message(sprintf("Collected %d Pegg validation runs (%s).",
                  nrow(out), PEGG_VALIDATION_VERSION))
  out
}

# ============================================================================
# run_pegg_seeds()
# ============================================================================
run_pegg_seeds <- function(
    params,
    seeds         = 1:20,
    results_dir   = "results/pegg_seeds/",
    n_cores       = max(1L, parallel::detectCores() - 1L),
    skip_existing = TRUE,
    verbose       = TRUE,
    ...           # passed to run_pegg_validation
) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  extra    <- list(...)
  defaults <- formals(run_pegg_validation)
  pick     <- function(nm) as.integer(if (!is.null(extra[[nm]])) extra[[nm]] else eval(defaults[[nm]]))
  scenario <- list(year_spinup = pick("year_spinup"), year_rust = pick("year_rust"),
                   year_fire   = pick("year_fire"),   year_survey = pick("year_survey"))
  n <- length(seeds)
  if (verbose) message(sprintf("Pegg ensemble: %d runs on %d core(s)", n, n_cores))

  parallel::mclapply(seeds, function(s) {
    p      <- params
    p$seed <- as.integer(s)
    id     <- pegg_run_id(p, scenario)
    out    <- file.path(results_dir, paste0(id, ".rds"))
    if (skip_existing && file.exists(out)) return(invisible(NULL))
    tryCatch(
      run_pegg_validation(p, ..., save_rds = out, plot = FALSE, verbose = FALSE),
      error = function(e) message(sprintf("  seed %d error: %s", s, conditionMessage(e)))
    )
    invisible(NULL)
  }, mc.cores = n_cores)

  df <- collect_pegg_results(results_dir)
  if (verbose && nrow(df) > 0)
    message(sprintf(
      "  fraction_resprout_dead: median=%.3f  [%.3f - %.3f]  (observed=%.2f)",
      median(df$fraction_resprout_dead, na.rm=TRUE),
      quantile(df$fraction_resprout_dead, 0.1, na.rm=TRUE),
      quantile(df$fraction_resprout_dead, 0.9, na.rm=TRUE),
      df$observed[1]))
  df
}

# ============================================================================
# plot_pegg_validation()  -- single run
# ============================================================================
plot_pegg_validation <- function(result, save_path = NULL,
                                  width_in = 12, height_in = 5, res = 150) {
  if (!is.null(save_path)) {
    grDevices::png(save_path, width=width_in, height=height_in,
                   units="in", res=res)
    on.exit(grDevices::dev.off(), add=TRUE)
  }
  if (requireNamespace("ggplot2", quietly=TRUE)) {
    .gg_pegg_single(result)
  } else {
    .base_pegg_single(result)
  }
}

.gg_pegg_single <- function(result) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  cen <- result$census
  yf  <- result$year_fire; ya <- result$year_assessment
  obs <- result$observed

  # Panel A: N_alive time series
  ymax <- max(cen$N_alive, na.rm=TRUE) * 1.1

  pA <- gg(cen, aes(x=year, y=N_alive)) +
    ggplot2::annotate("rect", xmin=result$year_spinup, xmax=result$year_rust,
                       ymin=0, ymax=ymax, fill="#dddddd", alpha=0.4) +
    ggplot2::annotate("text", x=(result$year_spinup+result$year_rust)/2,
                       y=ymax*0.97, label="spin-up\n(no fire)",
                       size=2.8, colour="#999", vjust=1) +
    ggplot2::geom_line(aes(y=N_IUCN), colour="#888780", linewidth=0.6,
                        linetype="dashed") +
    ggplot2::geom_line(colour="#333333", linewidth=0.8) +
    ggplot2::geom_vline(xintercept=result$year_rust,
                         colour="#e34948", linetype="dashed", linewidth=0.7) +
    ggplot2::geom_vline(xintercept=yf,
                         colour="#e08000", linewidth=1) +
    ggplot2::geom_vline(xintercept=ya,
                         colour="#2a78d6", linetype="dotted", linewidth=0.8) +
    ggplot2::annotate("text", x=result$year_rust+0.3, y=ymax*0.88,
                       label="rust", colour="#e34948", hjust=0, size=2.8) +
    ggplot2::annotate("text", x=yf+0.3, y=ymax*0.78,
                       label="fire\n(Black\nSummer)", colour="#e08000",
                       hjust=0, size=2.8) +
    ggplot2::annotate("text", x=ya+0.3, y=ymax*0.5,
                       label=sprintf("survey\n(Jan %d)", result$year_survey), colour="#2a78d6",
                       hjust=0, size=2.8) +
    ggplot2::scale_y_continuous(limits=c(0,ymax),
                                  expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="Year", y="Individuals (solid=N_alive, dashed=N_IUCN)",
                  title=sprintf(
                    "Pegg et al. validation  |  fire_p_fimp=%.2f  |  fraction_dead=%.3f (observed=%.2f)",
                    result$fire_p_fimp, result$fraction_dead, obs)) +
    ggplot2::theme_minimal(base_size=11)

  # Panel B: fate of resprouting survivors (denominator = N_resprouted)
  # Observed 48% is of this denominator, so this is the directly comparable panel.
  # Fire-direct deaths shown separately as an annotation.
  N_resp <- result$N_resprouted
  fates <- data.frame(
    fate = factor(c("Died\npost-fire","Still\nresprouting","Recovered\n& alive"),
                  levels=c("Died\npost-fire","Still\nresprouting","Recovered\n& alive")),
    frac = c(result$fraction_resprout_dead,
             result$n_resprouting / N_resp,
             result$n_recovered   / N_resp),
    fill = c("#fdae6b","#2a78d6","#a1d99b")
  )

  pB <- gg(fates, aes(x=fate, y=frac, fill=fate)) +
    ggplot2::geom_col(width=0.65, colour="white") +
    ggplot2::geom_hline(yintercept=obs, linetype="dashed",
                         colour="#e34948", linewidth=0.9) +
    ggplot2::annotate("text", x=0.55, y=obs+0.02,
                       label=sprintf("Pegg et al. %.0f%%\n(resprouters dead)", obs*100),
                       hjust=0, size=2.8, colour="#e34948") +
    ggplot2::annotate("text", x=2.5, y=0.97,
                       label=sprintf("fire-direct deaths\nexcluded from denom.\n(%.0f%% of pre-fire)",
                                      result$fraction_fire*100),
                       hjust=0.5, vjust=1, size=2.5, colour="#888780") +
    ggplot2::scale_fill_manual(
      values=setNames(fates$fill, levels(fates$fate)), guide="none") +
    ggplot2::scale_y_continuous(limits=c(0,1),
                                  expand=ggplot2::expansion(mult=c(0,0.05))) +
    ggplot2::labs(x=NULL, y="Fraction of resprouting survivors",
                  title=sprintf("Fate of resprouters at survey, Jan %d  (N=%d)",
                                result$year_survey, N_resp),
                  subtitle="Denominator = trees that survived fire and resprouted") +
    ggplot2::theme_minimal(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=9))

  if (!requireNamespace("scales", quietly=TRUE))
    pB <- pB + ggplot2::scale_y_continuous(limits=c(0,1),
                                              expand=ggplot2::expansion(mult=c(0,0.05)))

  if (requireNamespace("patchwork", quietly=TRUE)) {
    print(patchwork::wrap_plots(pA, pB, ncol=2, widths=c(1.8,1)))
  } else if (requireNamespace("gridExtra", quietly=TRUE)) {
    gridExtra::grid.arrange(pA, pB, ncol=2, widths=c(1.8,1))
  } else {
    print(pA); readline("Enter for fate panel..."); print(pB)
  }
}

.base_pegg_single <- function(result) {
  op <- par(mfrow=c(1,2), mar=c(4,4,3,1), mgp=c(2.3,0.7,0))
  on.exit(par(op), add=TRUE)
  cen <- result$census
  plot(cen$year, cen$N_alive, type="l", lwd=1.5,
       xlab="Year", ylab="N_alive",
       main=sprintf("Pegg validation (p_fimp=%.2f)", result$fire_p_fimp), las=1)
  lines(cen$year, cen$N_IUCN, lty=2, col="#888780")
  abline(v=result$year_rust, col="#e34948", lty=2)
  abline(v=result$year_fire, col="#e08000", lwd=2)
  abline(v=result$year_assessment, col="#2a78d6", lty=3)

  N_resp <- result$N_resprouted
  frac <- c(result$n_postfire_dead / N_resp,
            result$n_resprouting   / N_resp,
            result$n_recovered     / N_resp)
  cols <- c("#fdae6b","#2a78d6","#a1d99b")
  bp <- barplot(frac, col=cols, border="white", ylim=c(0,1),
                names.arg=c("Post-fire\ndead","Resprouting","Recovered"),
                ylab="Fraction of resprouting survivors",
                main=sprintf("Fate of resprouters, Jan %d (N=%d)",
                             result$year_survey, N_resp),
                las=1, cex.names=0.8)
  abline(h=result$observed, lty=2, lwd=1.5, col="#e34948")
  text(bp[1], result$observed+0.03,
       sprintf("Pegg %.0f%%", result$observed*100), adj=0, cex=0.8, col="#e34948")
}

# ============================================================================
# plot_pegg_ensemble()  -- across seeds
# ============================================================================
# Panel B shows the MEAN fate composition across runs (means of shares sum
# to 1; medians of each share need not). Panel A shades the 95% binomial
# interval of the observation (observed_n = 50 trees in Pegg et al.; set
# observed_n = NULL to omit it).
plot_pegg_ensemble <- function(df, observed = 0.48, observed_n = 50L,
                                save_path = NULL,
                                width_in = 10, height_in = 5, res = 150) {
  if (nrow(df) == 0L) stop("No results to plot.")
  n <- nrow(df)

  cat(sprintf("\n=== Pegg et al. ensemble (%d runs, observed=%.2f) ===\n", n, observed))
  cat(sprintf("  fraction_resprout_dead: median=%.3f  IQR=[%.3f, %.3f]  range=[%.3f, %.3f]\n",
              median(df$fraction_resprout_dead, na.rm=TRUE),
              quantile(df$fraction_resprout_dead, 0.25, na.rm=TRUE),
              quantile(df$fraction_resprout_dead, 0.75, na.rm=TRUE),
              min(df$fraction_resprout_dead, na.rm=TRUE),
              max(df$fraction_resprout_dead, na.rm=TRUE)))
  fm <- .pegg_fate_means(df)
  cat(sprintf("  mean fate of resprouters at survey: died=%.3f  still resprouting=%.3f  recovered=%.3f  (sum=%.3f)\n",
              fm[1], fm[2], fm[3], sum(fm)))
  ci <- .pegg_obs_ci(observed, observed_n)
  if (!is.null(ci))
    cat(sprintf("  observed 95%% binomial interval (n=%d): [%.3f, %.3f]; runs inside: %d%%\n",
                as.integer(observed_n), ci[1], ci[2],
                round(100*mean(df$fraction_resprout_dead >= ci[1] &
                               df$fraction_resprout_dead <= ci[2], na.rm=TRUE))))

  if (!is.null(save_path)) {
    grDevices::png(save_path, width=width_in, height=height_in,
                   units="in", res=res)
    on.exit(grDevices::dev.off(), add=TRUE)
  }

  if (requireNamespace("ggplot2", quietly=TRUE)) {
    .gg_pegg_ensemble(df, observed, ci)
  } else {
    .base_pegg_ensemble(df, observed)
  }
  invisible(NULL)
}

# Mean share of resprouters in each fate (died / still resprouting / recovered).
.pegg_fate_means <- function(df) {
  denom <- 1 - df$fraction_fire
  denom[denom <= 0] <- NA_real_
  c(mean(df$fraction_postfire / denom, na.rm=TRUE),
    mean(df$fraction_resprout / denom, na.rm=TRUE),
    mean(df$fraction_recover  / denom, na.rm=TRUE))
}

.pegg_obs_ci <- function(observed, n) {
  if (is.null(n)) return(NULL)
  as.numeric(stats::binom.test(round(observed * n), n)$conf.int)
}

.gg_pegg_ensemble <- function(df, observed, ci = NULL) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  n  <- nrow(df)

  # Panel A: histogram of fraction_resprout_dead
  pA <- gg(df, aes(x=fraction_resprout_dead)) +
    { if (!is.null(ci))
        ggplot2::annotate("rect", xmin=ci[1], xmax=ci[2], ymin=-Inf, ymax=Inf,
                          fill="#e34948", alpha=0.12) } +
    ggplot2::geom_histogram(binwidth=0.02, fill="#2a78d6",
                             colour="white", alpha=0.85, na.rm=TRUE) +
    ggplot2::geom_vline(xintercept=observed,
                         colour="#e34948", linewidth=1.2, linetype="dashed") +
    ggplot2::geom_vline(xintercept=median(df$fraction_resprout_dead, na.rm=TRUE),
                         colour="#333333", linewidth=1) +
    ggplot2::annotate("text", x=observed, y=Inf,
                       label=sprintf(" Pegg et al.\n %.0f%%\n(resprouters\ndead)", observed*100),
                       hjust=0, vjust=1.3, size=3, colour="#e34948") +
    ggplot2::annotate("text",
                       x=median(df$fraction_resprout_dead, na.rm=TRUE), y=Inf,
                       label=sprintf(" median\n %.0f%%",
                                      100*median(df$fraction_resprout_dead, na.rm=TRUE)),
                       hjust=0, vjust=1.3, size=3, colour="#333333") +
    ggplot2::scale_x_continuous(limits=c(0,1)) +
    ggplot2::labs(x="Fraction of resprouting survivors dead at survey",
                  y="Number of runs",
                  title=sprintf("Pegg et al. validation ensemble  (n=%d seeds)", n),
                  subtitle=paste0("Denominator = fire survivors that resprouted",
                                  if (!is.null(ci)) "  |  shaded = 95% CI of observation" else "")) +
    ggplot2::theme_minimal(base_size=11)

  # Panel B: fate breakdown of resprouting survivors (matched to observed denominator)
  # Guard against division by zero (runs where all trees killed by fire)
  fate_medians <- stats::setNames(.pegg_fate_means(df),
    c("Died\npost-fire", "Still\nresprouting", "Recovered\n& alive"))
  fate_df <- data.frame(
    fate = factor(names(fate_medians), levels=names(fate_medians)),
    frac = as.numeric(fate_medians),
    fill = c("#fdae6b","#2a78d6","#a1d99b")
  )

  pB <- gg(fate_df, aes(x=fate, y=frac, fill=fate)) +
    ggplot2::geom_col(width=0.65, colour="white") +
    ggplot2::geom_hline(yintercept=observed, linetype="dashed",
                         colour="#e34948", linewidth=0.8) +
    ggplot2::annotate("text", x=2, y=observed+0.03,
                       label=sprintf("Pegg %.0f%%", observed*100),
                       hjust=0.5, size=3, colour="#e34948") +
    ggplot2::scale_fill_manual(
      values=setNames(fate_df$fill, levels(fate_df$fate)), guide="none") +
    ggplot2::scale_y_continuous(limits=c(0,1),
                                  expand=ggplot2::expansion(mult=c(0,0.05))) +
    ggplot2::labs(x=NULL, y="Fraction of resprouting survivors",
                  title="Fate of resprouters at survey (mean over runs)",
                  subtitle="Shares sum to 1; fire-killed trees excluded") +
    ggplot2::theme_minimal(base_size=11)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    print(patchwork::wrap_plots(pA, pB, ncol=2))
  } else if (requireNamespace("gridExtra", quietly=TRUE)) {
    gridExtra::grid.arrange(pA, pB, ncol=2)
  } else {
    print(pA); readline("Enter for fate panel..."); print(pB)
  }
}

.base_pegg_ensemble <- function(df, observed) {
  op <- par(mfrow=c(1,2), mar=c(4,4,3,1), mgp=c(2.3,0.7,0))
  on.exit(par(op), add=TRUE)
  hist(df$fraction_resprout_dead, breaks=seq(0,1,by=0.05), col="#2a78d6",
       border="white", xlab="Fraction of resprouters dead", ylab="Runs",
       main=sprintf("Pegg ensemble (n=%d)", nrow(df)), las=1)
  abline(v=observed, col="#e34948", lwd=2, lty=2)
  abline(v=median(df$fraction_resprout_dead, na.rm=TRUE), col="#333333", lwd=2)
  legend("topright", bty="n", lwd=2, col=c("#e34948","#333333"),
         legend=c(sprintf("Pegg %.0f%%", observed*100), "median"))

  fate_m <- .pegg_fate_means(df)
  bp <- barplot(fate_m, col=c("#fdae6b","#2a78d6","#a1d99b"),
                border="white", ylim=c(0,1),
                names.arg=c("Post-fire\ndead","Resprouting","Recovered"),
                ylab="Fraction of resprouting survivors",
                main="Fate breakdown (resprouters)", las=1, cex.names=0.8)
  abline(h=observed, lty=2, col="#e34948", lwd=1.5)
}
