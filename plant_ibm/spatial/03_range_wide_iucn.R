# spatial/03_range_wide_iucn.R
#
# Range-wide IUCN Criterion A4 assessment for Melaleuca nodosa.
#
# Samples N site-seed pairs (Option X: site drawn ∝ weight_n_raw, seed
# uniformly random), runs run_iucn_a4() for each, and compiles the
# distribution of A4 decline across the species' range. Each run gives one
# decline: mean N_IUCN over the 11 yr ending the year before rust arrival
# (2009) vs mean N_IUCN over the 11 yr ending 3*T_g later. Because sites are
# drawn in proportion to weight_n_raw, the plain median across runs is the
# abundance-weighted central tendency over space, fire and rust variation;
# that median is the value reported for the assessment.
#
# DESIGN
# ------
# Each run represents one realisation of the species at one location:
#   - Site drawn proportionally to weight_n_raw (coastal high-density
#     cells appear in proportion to their relative importance)
#   - Site-specific params override: fire_prob_annual, fire_p_fimp,
#     annual_env_alpha, annual_env_beta
#   - All other params: from params.R (working M. nodosa calibration)
#   - Seed: drawn uniformly; determines stochastic fire timing and
#     annual rust pressure draws within a run
#
# RESISTANCE NOTE
# ---------------
# resist_freq0 = c(0) in params.R (fully susceptible population).
# Site-specific resistance data are not yet available. When data arrive,
# add a resist_freq0 column to mn_site_params.csv and override per-site
# in the build_site_params() function below (search "RESISTANCE HOOK").
#
# RESTARTABILITY
# --------------
# skip_existing=TRUE means interrupted runs resume from where they stopped.
# Run the same command again to add more runs to an existing results_dir.
#
# USAGE
# -----
# Source all model and runner files first (runner/iucn_a4.R must be
# sourced before this file: ASSESS is built with iucn_assess_spec()), then:
#   source("spatial/03_range_wide_iucn.R")
#   df <- run_range_wide_iucn(n_runs=200, n_cores=8)
#   plot_range_wide_iucn(df)
#   plot_range_wide_iucn(df, save_path="outputs/range_wide_iucn.png")

# ============================================================================
# Configuration
# ============================================================================
SITE_PARAMS_PATH <- "spatial/data/mn_site_params.csv"

# IUCN simulation phasing (calendar years)
YEAR_SPINUP          <- 1600L   # 200-yr spin-up
YEAR_GENTIME_START   <- 1810L
YEAR_GENTIME_END     <- 1910L
YEAR_RUST            <- 2010L
YEAR_ASSESS          <- 2026L   # the present; A4 windows straddle this
YEAR_END             <- 2210L   # >= 2009 + 3*T_g (comparison window end), capped at YEAR_ASSESS + 100
AVG_YEARS            <- 11L     # length of the reference and comparison averaging windows

# Assessment settings in one object: used for the run AND its cache ID
ASSESS <- iucn_assess_spec(
  year_spinup = YEAR_SPINUP, year_gentime_start = YEAR_GENTIME_START,
  year_gentime_end = YEAR_GENTIME_END, year_rust = YEAR_RUST,
  year_assess = YEAR_ASSESS, year_end = YEAR_END, avg_years = AVG_YEARS)

# ============================================================================
# build_site_params(): override base params for one site
# ============================================================================
build_site_params <- function(base_params, site_row) {
  p <- base_params

  # Site-specific fire regime
  p$fire_prob_annual <- site_row$fire_prob_annual
  p$fire_p_fimp      <- site_row$fire_p_fimp
  p$fire_years       <- integer(0)   # must use fire_prob_annual for long runs

  # Site-specific rust pressure (from SDM, anchored at Yarringully)
  p$annual_env_alpha <- site_row$annual_env_alpha
  p$annual_env_beta  <- site_row$annual_env_beta
  # microclim_alpha/beta held constant (within-site individual variation,
  # assumed similar across all M. nodosa sites)

  # --- RESISTANCE HOOK -------------------------------------------------------
  # When site-level genetic resistance data are available, add a
  # resist_freq0 column to mn_site_params.csv and uncomment:
  # p$resist_freq0 <- c(site_row$resist_freq0)
  # For now: resist_freq0 = c(0) from params.R (fully susceptible everywhere)
  # ---------------------------------------------------------------------------

  p
}

# ============================================================================
# run_range_wide_iucn()
# ============================================================================
run_range_wide_iucn <- function(
    n_runs        = 200L,
    n_cores       = max(1L, parallel::detectCores() - 1L),
    results_dir   = "results/range_wide_iucn/",
    skip_existing = TRUE,
    seed_max      = 9999L,
    verbose       = TRUE
) {
  dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

  # Load site params
  site_params <- readr::read_csv(SITE_PARAMS_PATH, show_col_types=FALSE)
  message(sprintf("Loaded %d M. nodosa sites from %s", nrow(site_params), SITE_PARAMS_PATH))

  # Load base params (params.R must already be sourced)
  base_params <- get_default_params()
  base_params$rust_start_year <- YEAR_RUST    # overridden inside run_iucn_a4 anyway
  base_params$fire_years      <- integer(0)   # stochastic fire throughout

  # Sample N site-seed pairs: site ∝ weight_n_raw, seed uniform
  set.seed(42L)
  site_idx <- sample(nrow(site_params), size=n_runs, replace=TRUE,
                     prob=site_params$weight_n_raw)
  seeds    <- sample(seed_max, size=n_runs, replace=TRUE)

  message(sprintf(
    "Running %d site-seed pairs on %d core(s)  [sites drawn ∝ weight_n_raw]",
    n_runs, n_cores))

  parallel::mclapply(seq_len(n_runs), function(i) {
    site_row <- site_params[site_idx[i], ]
    p        <- build_site_params(base_params, site_row)
    p$seed   <- seeds[i]

    id  <- iucn_run_id(p, ASSESS)
    out <- file.path(results_dir, paste0(id, ".rds"))
    if (skip_existing && file.exists(out)) return(invisible(NULL))

    tryCatch(
      run_iucn_a4(
        p,
        year_spinup        = YEAR_SPINUP,
        year_gentime_start = YEAR_GENTIME_START,
        year_gentime_end   = YEAR_GENTIME_END,
        year_rust          = YEAR_RUST,
        year_assess        = YEAR_ASSESS,
        year_end           = YEAR_END,
        avg_years          = AVG_YEARS,
        save_rds           = out,
        plot               = FALSE,
        verbose            = FALSE
      ),
      error = function(e)
        message(sprintf("  run %d (site %d, seed %d): %s",
                        i, site_idx[i], seeds[i], conditionMessage(e)))
    )
    invisible(NULL)
  }, mc.cores=n_cores)

  # Collect and join site metadata
  df <- collect_iucn_results(results_dir)
  if (nrow(df) == 0L) { warning("No results collected."); return(df) }

  # Join site info back (run_id → site via seed+site params)
  # Reconstruct the mapping from the same random seed
  set.seed(42L)
  site_idx2 <- sample(nrow(site_params), size=n_runs, replace=TRUE,
                       prob=site_params$weight_n_raw)
  seeds2    <- sample(seed_max, size=n_runs, replace=TRUE)

  run_map <- do.call(rbind, lapply(seq_len(n_runs), function(i) {
    site_row <- site_params[site_idx2[i], ]
    p        <- build_site_params(base_params, site_row)
    p$seed   <- seeds2[i]
    data.frame(
      run_id            = iucn_run_id(p, ASSESS),
      lon_centre        = site_row$lon_centre,
      lat_centre        = site_row$lat_centre,
      site_fire_prob    = site_row$fire_prob_annual,
      site_fire_p_fimp  = site_row$fire_p_fimp,
      site_annual_env   = site_row$annual_env_alpha /
                           (site_row$annual_env_alpha + site_row$annual_env_beta),
      site_sdm          = site_row$sdm_suitability,
      weight_n_raw      = site_row$weight_n_raw,
      stringsAsFactors  = FALSE
    )
  }))
  run_map <- run_map[!duplicated(run_map$run_id), ]

  df <- dplyr::left_join(df, run_map, by="run_id")

  if (verbose) {
    cat(sprintf("\n=== Range-wide IUCN ensemble (%d runs) ===\n", nrow(df)))
    cat(sprintf("  pct_decline:   median=%.1f%%  IQR=[%.1f%%, %.1f%%]\n",
                median(df$pct_decline, na.rm=TRUE),
                quantile(df$pct_decline, 0.25, na.rm=TRUE),
                quantile(df$pct_decline, 0.75, na.rm=TRUE)))
    cat(sprintf("  T_g (yr):      median=%.1f\n", median(df$T_g, na.rm=TRUE)))
    shares <- function(x) round(100*prop.table(table(factor(
      vapply(x, iucn_category, ""), levels=c("LC/NT","VU","EN","CR")))))
    sh <- shares(df$pct_decline)
    cat(sprintf("  IUCN outcomes (mature individuals): LC/NT=%d%%  VU=%d%%  EN=%d%%  CR=%d%%\n",
                sh[[1]], sh[[2]], sh[[3]], sh[[4]]))
    if ("pct_decline_flowering" %in% names(df)) {
      sf <- shares(df$pct_decline_flowering)
      cat(sprintf("  IUCN outcomes (flowering only):     LC/NT=%d%%  VU=%d%%  EN=%d%%  CR=%d%%   median=%.1f%%\n",
                  sf[[1]], sf[[2]], sf[[3]], sf[[4]],
                  median(df$pct_decline_flowering, na.rm=TRUE)))
    }
  }

  df
}

# ============================================================================
# plot_range_wide_iucn()
# ============================================================================
plot_range_wide_iucn <- function(df, save_path=NULL,
                                  width_in=14, height_in=10, res=150) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 required for plotting")

  if (!is.null(save_path)) {
    grDevices::png(save_path, width=width_in, height=height_in,
                   units="in", res=res)
    on.exit(grDevices::dev.off(), add=TRUE)
  }

  gg <- ggplot2::ggplot; aes <- ggplot2::aes

  aus_states <- ozmaps::ozmap_states
  xlim <- c(148.5, 154.5); ylim <- c(-37.5, -24.5)
  base_map <- list(
    ggplot2::geom_sf(data=aus_states, inherit.aes=FALSE,
                     fill="#f5f3ee", colour="#cccccc", linewidth=0.3),
    ggplot2::coord_sf(xlim=xlim, ylim=ylim, expand=FALSE),
    ggplot2::theme_minimal(base_size=10),
    ggplot2::theme(panel.grid=ggplot2::element_line(colour="#e8e8e8",linewidth=0.2),
                   legend.position="bottom")
  )

  n <- nrow(df)

  # Panel A: histogram of pct_decline
  pA <- gg(df, aes(x=pct_decline)) +
    ggplot2::geom_histogram(binwidth=2, fill="#2a78d6", colour="white",
                             alpha=0.85, na.rm=TRUE) +
    ggplot2::geom_vline(xintercept=c(30,50,80),
                         colour=c("#eda100","#e34948","#8b0000"),
                         linetype="dashed", linewidth=0.8) +
    ggplot2::annotate("text", x=c(30,50,80), y=Inf,
                       label=c("VU","EN","CR"), hjust=-0.15, vjust=1.4,
                       size=3, colour=c("#eda100","#e34948","#8b0000")) +
    ggplot2::geom_vline(xintercept=median(df$pct_decline, na.rm=TRUE),
                         colour="#333333", linewidth=1) +
    ggplot2::scale_x_continuous(limits=c(min(0, 10*floor(min(df$pct_decline, na.rm=TRUE)/10)), 100)) +
    ggplot2::labs(x="Decline in mature individuals from pre-rust to 3 generations later (%)",
                  y="Number of runs",
                  title=sprintf("Range-wide IUCN A4 ensemble  (n=%d)", n),
                  subtitle=sprintf(
                    "LC/NT: %d%%   VU: %d%%   EN: %d%%   CR: %d%%   |   median = %.0f%%",
                    round(100*mean(df$pct_decline <  30, na.rm=TRUE)),
                    round(100*mean(df$iucn_vu==1 & df$iucn_en==0, na.rm=TRUE)),
                    round(100*mean(df$iucn_en==1 & df$iucn_cr==0, na.rm=TRUE)),
                    round(100*mean(df$iucn_cr==1, na.rm=TRUE)),
                    median(df$pct_decline, na.rm=TRUE))) +
    ggplot2::theme_minimal(base_size=11)

  # Panel B: map of median pct_decline per site location
  site_summary <- df |>
    dplyr::group_by(lon_centre, lat_centre) |>
    dplyr::summarise(median_decline = median(pct_decline, na.rm=TRUE),
                     n_runs = dplyr::n(), .groups="drop")

  pB <- gg(site_summary, aes(x=lon_centre, y=lat_centre,
                               colour=median_decline, size=n_runs)) +
    base_map +
    ggplot2::geom_point(alpha=0.85) +
    ggplot2::scale_colour_gradient2(low="#2a78d6", mid="#fdae6b", high="#e34948",
                                     midpoint=50, name="Median\ndecline (%)") +
    ggplot2::scale_size_continuous(name="Runs", range=c(1.5,5)) +
    ggplot2::labs(title="Median A4 decline by site",
                  x=NULL, y=NULL)

  # Panel C: fire × rust interaction
  pC <- gg(df, aes(x=site_fire_prob, y=site_annual_env, colour=pct_decline)) +
    ggplot2::geom_point(alpha=0.6, size=1.8, na.rm=TRUE) +
    ggplot2::scale_colour_gradient2(low="#2a78d6", mid="#fdae6b", high="#e34948",
                                     midpoint=50, name="Decline (%)") +
    ggplot2::scale_x_continuous(
      sec.axis=ggplot2::sec_axis(~1/., name="FRI (years)",
                                  breaks=1/c(0.5,0.1,0.05,0.02),
                                  labels=c("2","10","20","50"))) +
    ggplot2::labs(x="fire_prob_annual", y="annual_env_mean (rust pressure)",
                  title="Fire × rust interaction surface") +
    ggplot2::theme_minimal(base_size=11)

  # Panel D: IUCN category bars
  cat_df <- data.frame(
    category = factor(c("LC/NT","VU","EN","CR"),
                      levels=c("LC/NT","VU","EN","CR")),
    pct = c(mean(df$pct_decline <  30, na.rm=TRUE),
            mean(df$iucn_vu==1 & df$iucn_en==0, na.rm=TRUE),
            mean(df$iucn_en==1 & df$iucn_cr==0, na.rm=TRUE),
            mean(df$iucn_cr==1, na.rm=TRUE)) * 100,
    fill = c("#888780","#eda100","#e34948","#8b0000")
  )
  pD <- gg(cat_df, aes(x=category, y=pct, fill=category)) +
    ggplot2::geom_col(width=0.6, colour="white") +
    ggplot2::scale_fill_manual(values=setNames(cat_df$fill, cat_df$category),
                                guide="none") +
    ggplot2::scale_y_continuous(limits=c(0,100),
                                  expand=ggplot2::expansion(mult=c(0,0.05))) +
    ggplot2::labs(x=NULL, y="% of runs",
                  title="IUCN outcome distribution") +
    ggplot2::theme_minimal(base_size=11)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    print(patchwork::wrap_plots(pA, pB, pC, pD, ncol=2L) +
      patchwork::plot_annotation(
        title="M. nodosa range-wide IUCN Criterion A4 assessment",
        theme=ggplot2::theme(plot.title=ggplot2::element_text(size=12, face="bold"))))
  } else {
    for (p in list(pA,pB,pC,pD)) { print(p); readline("Enter...") }
  }

  invisible(df)
}
