# runner/pegg_validation.R
#
# Validation of the IBM against Pegg et al. (2025), who measured mortality
# in resprouting Melaleuca nodosa following the Black Summer fires (Jan 2020)
# at sites in northern New South Wales. At a Jan 2022 follow-up (~2 yr
# post-fire), approximately 40% of pre-fire M. nodosa were dead.
#
# SIMULATION DESIGN
# -----------------
# Spin-up (no fire): year_spinup → year_rust (default 1910→2010, 100 yr)
# Rust period:       year_rust   → year_fire  (default 2010→2020, 10 yr)
# Fire:              year_fire                (single deterministic event)
# Post-fire:         year_fire+1 → year_end  (default 2021→2025)
#
# Fire intensity: fire_p_fimp controls what fraction of plants are impacted.
# For the Black Summer fires (very intense), values near 1.0 are appropriate.
# This is the primary sensitivity parameter to explore.
#
# PRIMARY STATISTIC
# -----------------
# Fraction of pre-fire individuals dead by the assessment date:
#
#   fraction_dead = (fire_dead + post_fire_dead) / N_prefire
#
# where N_prefire = individuals alive at END of (year_fire - 1), before fire.
# fire_dead      = N_prefire individuals with death_year == year_fire
# post_fire_dead = N_prefire individuals with death_year in
#                  (year_fire+1):year_assessment
#
# This directly matches the Pegg et al. measurement (dead pre-fire trees at
# the follow-up census), without confounding with new post-fire recruits.
#
# FATE BREAKDOWN
# --------------
# All N_prefire individuals are assigned one of four fates at year_assessment:
#   (1) Killed directly by fire (year_fire deaths)
#   (2) Resprouted then died: rust/background mortality during window
#   (3) Still resprouting at assessment
#   (4) Recovered: alive, resprout=FALSE
# These four are mutually exclusive and sum to N_prefire.
#
# USAGE
# -----
#   source all model files, then:
#   source("runner/params_utils.R")
#   source("runner/pegg_validation.R")
#
#   # Single run with figure
#   result <- run_pegg_validation(get_default_params())
#
#   # Ensemble across seeds
#   df <- run_pegg_seeds(get_default_params(), seeds = 1:50, n_cores = 8L)
#   plot_pegg_ensemble(df, observed = 0.40)
#
#   # Sensitivity to fire intensity
#   df_hi <- run_pegg_seeds(within(get_default_params(),
#                                   {fire_p_fimp <- 1.0}), seeds = 1:30)
#   df_lo <- run_pegg_seeds(within(get_default_params(),
#                                   {fire_p_fimp <- 0.7}), seeds = 1:30)

# ============================================================================
# compute_pegg_stats()
# ============================================================================
# Extract the primary validation statistic and fate breakdown from a
# completed run's pop data frame and census.
# year_spinup is needed to correctly interpret founder birth_year encoding.
compute_pegg_stats <- function(pop, census, year_fire, year_assessment,
                                year_spinup) {

  # Calendar birth year (founders have birth_year = -(age at spinup))
  cal_birth <- ifelse(pop$birth_year < year_spinup,
                       year_spinup + pop$birth_year,
                       pop$birth_year)

  # Pre-fire cohort: alive at end of (year_fire - 1)
  # A plant was alive at end of that year if it either:
  #   (a) is still alive at run end, OR
  #   (b) died in year_fire or later
  prefire_year <- year_fire - 1L
  existed      <- cal_birth <= prefire_year
  alive_prefire <- existed & (pop$alive |
                               (!is.na(pop$death_year) & pop$death_year >= year_fire))
  N_prefire <- sum(alive_prefire)

  if (N_prefire == 0L) {
    warning("No pre-fire individuals found. Check year_fire and year_spinup.")
    return(list(N_prefire = 0L, fraction_dead = NA_real_,
                fire_dead = 0L, post_fire_dead = 0L,
                resprouting_2022 = 0L, recovered_2022 = 0L))
  }

  pf_pop <- pop[alive_prefire, ]

  # Fate classification at year_assessment
  fire_dead       <- !is.na(pf_pop$death_year) & pf_pop$death_year == year_fire
  post_fire_dead  <- !is.na(pf_pop$death_year) &
                     pf_pop$death_year >  year_fire &
                     pf_pop$death_year <= year_assessment
  resprouting_now <- pf_pop$alive & pf_pop$resprout
  recovered       <- pf_pop$alive & !pf_pop$resprout

  n_fire_dead      <- sum(fire_dead)
  n_postfire_dead  <- sum(post_fire_dead)
  n_resprouting    <- sum(resprouting_now)
  n_recovered      <- sum(recovered)

  fraction_dead    <- (n_fire_dead + n_postfire_dead) / N_prefire
  fraction_fire    <- n_fire_dead     / N_prefire
  fraction_postfire<- n_postfire_dead / N_prefire
  fraction_resprout<- n_resprouting   / N_prefire
  fraction_recover <- n_recovered     / N_prefire

  # Primary validation statistic: fraction of RESPROUTING SURVIVORS that died
  # by year_assessment. Denominator = trees that survived fire and entered
  # resprout (= N_prefire - fire_dead). This matches the Pegg et al. metric:
  # of the trees that showed resprouting during the survey, 48% were dead
  # at Jan 2022. Trees killed outright by fire are excluded from both
  # numerator and denominator.
  N_resprouted           <- N_prefire - n_fire_dead
  fraction_resprout_dead <- if (N_resprouted > 0L)
    n_postfire_dead / N_resprouted else NA_real_

  # N_IUCN (reproductive adults) from census at prefire_year and assessment
  N_iucn_prefire  <- census$N_IUCN[census$year == prefire_year]
  N_iucn_assess   <- census$N_IUCN[census$year == year_assessment]
  if (length(N_iucn_prefire) == 0L) N_iucn_prefire <- NA_real_
  if (length(N_iucn_assess)  == 0L) N_iucn_assess  <- NA_real_

  list(
    N_prefire              = N_prefire,
    N_resprouted           = N_resprouted,
    n_fire_dead            = n_fire_dead,
    n_postfire_dead        = n_postfire_dead,
    n_resprouting          = n_resprouting,
    n_recovered            = n_recovered,
    fraction_dead          = fraction_dead,          # of all pre-fire trees
    fraction_fire          = fraction_fire,
    fraction_postfire      = fraction_postfire,
    fraction_resprout      = fraction_resprout,
    fraction_recover       = fraction_recover,
    fraction_resprout_dead = fraction_resprout_dead, # PRIMARY: of resprouters
    N_iucn_prefire         = N_iucn_prefire,
    N_iucn_assess          = N_iucn_assess
  )
}

# ============================================================================
# run_pegg_validation()
# ============================================================================
run_pegg_validation <- function(
    params,
    year_spinup     = 1910L,
    year_rust       = 2010L,
    year_fire       = 2020L,
    year_assessment = 2022L,
    year_end        = 2025L,
    observed        = 0.48,     # Pegg et al.: 48% of resprouting trees dead by Jan 2022
    save_rds        = NULL,
    plot            = TRUE,
    verbose         = TRUE
) {
  p                    <- params
  p$rust_start_year    <- year_rust
  p$fire_years         <- as.integer(year_fire)
  p$fire_prob_annual   <- 0
  p$n_years            <- year_end - year_spinup

  # All trees at the survey site were impacted by the Black Summer fires.
  # fire_p_fimp = 1 is fixed for this validation (not a free parameter).
  p$fire_p_fimp        <- 1.0

  if (verbose) message(sprintf(
    "Pegg validation: spin=%d  rust=%d  fire=%d  assess=%d  p_fimp=1 (fixed)",
    year_spinup, year_rust, year_fire, year_assessment))

  res <- run_simulation(p, year0 = year_spinup, verbose = FALSE)

  stats <- compute_pegg_stats(res$individuals, res$census,
                               year_fire, year_assessment, year_spinup)

  if (verbose) message(sprintf(
    "  fraction_resprout_dead=%.3f (observed=%.2f)  |  N_resprouted=%d  fire_dead=%.3f",
    stats$fraction_resprout_dead, observed,
    stats$N_resprouted, stats$fraction_fire))

  result <- c(
    list(
      year_spinup     = year_spinup,
      year_rust       = year_rust,
      year_fire       = year_fire,
      year_assessment = year_assessment,
      year_end        = year_end,
      observed        = observed,
      census          = res$census,
      fire_p_fimp     = p$fire_p_fimp,
      run_id          = make_run_id(p),
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
# Single-row data.frame for batch collection (same pattern as summarise_iucn_a3).
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
    summarise_pegg(x)
  })
  rows <- rows[!sapply(rows, is.null)]
  out  <- do.call(rbind, rows)
  out  <- out[!duplicated(out$run_id), ]
  rownames(out) <- NULL
  message(sprintf("Collected %d Pegg validation runs.", nrow(out)))
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
  n <- length(seeds)
  if (verbose) message(sprintf("Pegg ensemble: %d runs on %d core(s)", n, n_cores))

  parallel::mclapply(seeds, function(s) {
    p      <- params
    p$seed <- as.integer(s)
    id     <- make_run_id(p)
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
                       label="assessment\n(Jan 2022)", colour="#2a78d6",
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
                  title=sprintf("Fate of resprouters at %d  (N=%d)",
                                result$year_assessment, N_resp),
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
                main=sprintf("Fate of resprouters at %d (N=%d)",
                             result$year_assessment, N_resp),
                las=1, cex.names=0.8)
  abline(h=result$observed, lty=2, lwd=1.5, col="#e34948")
  text(bp[1], result$observed+0.03,
       sprintf("Pegg %.0f%%", result$observed*100), adj=0, cex=0.8, col="#e34948")
}

# ============================================================================
# plot_pegg_ensemble()  -- across seeds
# ============================================================================
plot_pegg_ensemble <- function(df, observed = 0.40,
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
  cat(sprintf("  Runs within 5%% of observed: %d%%\n",
              round(100*mean(abs(df$fraction_resprout_dead - observed) <= 0.05, na.rm=TRUE))))
  cat(sprintf("  median fate breakdown (of resprouting survivors): post-fire=%.3f  resprout=%.3f  recovered=%.3f\n",
              median(df$fraction_postfire / (1-df$fraction_fire), na.rm=TRUE),
              median(df$fraction_resprout / (1-df$fraction_fire), na.rm=TRUE),
              median(df$fraction_recover  / (1-df$fraction_fire), na.rm=TRUE)))

  if (!is.null(save_path)) {
    grDevices::png(save_path, width=width_in, height=height_in,
                   units="in", res=res)
    on.exit(grDevices::dev.off(), add=TRUE)
  }

  if (requireNamespace("ggplot2", quietly=TRUE)) {
    .gg_pegg_ensemble(df, observed)
  } else {
    .base_pegg_ensemble(df, observed)
  }
  invisible(NULL)
}

.gg_pegg_ensemble <- function(df, observed) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  n  <- nrow(df)

  # Panel A: histogram of fraction_resprout_dead
  pA <- gg(df, aes(x=fraction_resprout_dead)) +
    ggplot2::geom_histogram(binwidth=0.02, fill="#2a78d6",
                             colour="white", alpha=0.85) +
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
    ggplot2::labs(x="Fraction of resprouting survivors dead by assessment",
                  y="Number of runs",
                  title=sprintf("Pegg et al. validation ensemble  (n=%d seeds)", n),
                  subtitle="Denominator = trees that survived fire and resprouted  |  fire_p_fimp = 1 (fixed)") +
    ggplot2::theme_minimal(base_size=11)

  # Panel B: fate breakdown of resprouting survivors (matched to observed denominator)
  # Guard against division by zero (runs where all trees killed by fire)
  denom <- 1 - df$fraction_fire
  denom[denom <= 0] <- NA_real_
  fate_medians <- c(
    "Died\npost-fire"    = median(df$fraction_postfire / denom, na.rm=TRUE),
    "Still\nresprouting" = median(df$fraction_resprout / denom, na.rm=TRUE),
    "Recovered\n& alive" = median(df$fraction_recover  / denom, na.rm=TRUE)
  )
  fate_medians <- pmin(pmax(fate_medians, 0), 1)  # safety clamp
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
                  title="Median fate breakdown",
                  subtitle="Denominator = resprouting survivors (excl. fire-direct deaths)") +
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

  denom <- 1 - df$fraction_fire
  denom[denom <= 0] <- NA_real_
  fate_m <- c(median(df$fraction_postfire / denom, na.rm=TRUE),
              median(df$fraction_resprout / denom, na.rm=TRUE),
              median(df$fraction_recover  / denom, na.rm=TRUE))
  fate_m <- pmin(pmax(fate_m, 0), 1)
  bp <- barplot(fate_m, col=c("#fdae6b","#2a78d6","#a1d99b"),
                border="white", ylim=c(0,1),
                names.arg=c("Post-fire\ndead","Resprouting","Recovered"),
                ylab="Fraction of resprouting survivors",
                main="Fate breakdown (resprouters)", las=1, cex.names=0.8)
  abline(h=observed, lty=2, col="#e34948", lwd=1.5)
}
