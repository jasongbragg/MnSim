# runner/iucn_a4.R
#
# IUCN Criterion A4 assessment: the decline in mature individuals (N_IUCN)
# over a window of 10 years or three generations, whichever is longer, that
# includes both the past and the future (reaching at most 100 years ahead).
# Here the window is fixed: it starts at the arrival of rust.
#
# MATURE INDIVIDUALS
# ------------------
# N_IUCN = alive plants of reproductive age (age >= age_first_flower),
# including plants kept from flowering in a given year by fire (resprouting)
# or rust; they remain capable of reproduction (IUCN Guidelines v16 s4.3).
# The stricter flowering-only count (census$N_flowering) is assessed over the
# same window in every run and stored alongside as pct_decline_flowering.
#
# WHY A4
# ------
# Myrtle rust arrived in 2010 and the decline is ongoing, so the relevant
# reduction spans past and future. A3 (future-only, from the present) would
# miss the part of the decline that has already happened; A2 (past-only)
# would miss the projected part. A4 takes the most severe window that
# straddles the present.
#
# WINDOW DEFINITION (implemented in assess_a4())
# ----------------------------------------------
# One fixed window per run, anchored on the arrival of rust:
#   L          = max(10, 3 * T_g) years, rounded to whole years
#   t0         = year_rust - 1             (last pre-rust year)
#   reference  = mean N_IUCN over the avg_years (11) years ending at t0
#   comparison = mean N_IUCN over the avg_years years ending at t0 + L
#   decline    = 100 * (N_ref - N_compare) / N_ref
# The window starts in the past (2010) and ends in the future, as A4
# requires. If t0 + L would fall more than max_future (100) years after
# year_assess, the comparison end is capped there (with a warning).
# If the population goes extinct, years after the last census row count
# as N = 0.
#
# T_g, and hence L, is estimated within each run from that run's own
# pre-rust dynamics, then applied to the same run's post-rust trajectory.
#
# The assessment is one decline per run. Ensembles over seeds, sites and
# parameter draws are summarised as the frequency distribution of runs across
# the IUCN categories (LC/NT, VU, EN, CR), with the median and quantiles as
# the central tendency -- see run_iucn_seeds() / plot_iucn_ensemble() and
# spatial/03_range_wide_iucn.R.
#
# GENERATION TIME
# ---------------
# Mean maternal age of offspring born in a pre-rust window (IUCN
# definition: mean age of parents of the current cohort), from the
# natural (pre-disturbance) dynamics. See compute_gen_time().
#
# SIMULATION PHASING (calendar years, single continuous run)
# ----------------------------------------------------------
#   year_spinup        start of run (e.g. 1600)  ─┐ spin-up: discard
#   year_gentime_start                            ─┤ gen-time window
#   year_gentime_end                              ─┤ pre-rust equilibrium
#   year_rust          rust arrives (2010)        ─┤ reference: 11 yr ending 2009
#   year_assess        the present (2026)         ─┤ window must straddle this
#   year_end           end of run                 ─┘ must be >= 2009 + 3*T_g
#
# FIRE NOTE
# ---------
# Use fire_prob_annual > 0 and fire_years = integer(0) for these long runs.
#
# CACHING / RUN IDs
# -----------------
# iucn_run_id(params, assess) hashes the params AND the assessment settings
# (years, averaging window, criterion). Cached .rds files therefore never get
# reused when the assessment settings change, and old A3 results in a
# results directory are not mistaken for A4 results.
#
# SINGLE-RUN USAGE
# ----------------
#   source("runner/params_utils.R"); source("runner/iucn_a4.R")
#   # (source all model files and params first)
#   result <- run_iucn_a4(get_default_params())
#   plot_iucn_a4(result)
#   summarise_iucn_a4(result)

IUCN_CRITERION <- "A4"

# ============================================================================
# compute_gen_time()
# ============================================================================
# Mean maternal age at reproduction (IUCN generation length definition)
# over a specified window of births. Handles the birth_year encoding:
# founders store birth_year = -age (relative to year_spinup); recruits
# store birth_year = t (calendar year). Both are adjusted to calendar
# years using year_spinup before computing maternal age.
compute_gen_time <- function(pop, window_start, window_end, year_spinup) {

  # Adjust birth_year to calendar year:
  # Founders store birth_year = -age (relative to simulation start, so negative)
  # Recruits store birth_year = calendar year (>= year_spinup)
  cal_birth <- ifelse(pop$birth_year < year_spinup,
                       year_spinup + pop$birth_year,
                       pop$birth_year)

  # Offspring born in window with known mother
  in_window <- !is.na(pop$mother_id) &
               cal_birth >= window_start &
               cal_birth <= window_end
  offs <- pop[in_window, ]

  if (nrow(offs) == 0L) {
    warning("No offspring with known mothers in gen-time window [",
            window_start, ", ", window_end, "]. Returning NA.")
    return(list(T_g = NA_real_, mother_ages = numeric(0), n_offspring = 0L))
  }

  offs_cal_birth   <- cal_birth[in_window]
  mother_idx       <- match(offs$mother_id, pop$id)
  mother_cal_birth <- cal_birth[mother_idx]
  mother_alive     <- pop$alive[mother_idx]

  # FILTER 1: deceased mothers only.
  # Living mothers have incomplete reproductive histories -- they may still
  # produce offspring after window_end, so their mean age computed from births
  # seen so far is younger than their true lifetime average. Dead mothers have
  # complete histories, avoiding this downward bias on T_g.
  dead_mother <- !is.na(mother_alive) & !mother_alive
  n_excl_alive <- sum(!dead_mother, na.rm = TRUE)
  if (n_excl_alive > 0L)
    message(sprintf(
      "  compute_gen_time: excluded %d offspring of still-living mothers.",
      n_excl_alive))

  # FILTER 2: mothers born during the simulation (recruits, not founders).
  # Founders are initialised with arbitrary ages drawn from the equilibrium
  # age distribution; their birth_year is back-calculated as
  # year_spinup - age and does not represent an actual birth event in the
  # model. Their apparent maternal age therefore depends on the arbitrary
  # initial age assignment rather than real simulated dynamics.
  # Recruits (cal_birth >= year_spinup) have genuine birth events.
  born_in_sim <- !is.na(mother_cal_birth) & mother_cal_birth >= year_spinup
  n_excl_founder <- sum(dead_mother & !born_in_sim, na.rm = TRUE)
  if (n_excl_founder > 0L)
    message(sprintf(
      "  compute_gen_time: excluded %d offspring of founding mothers (arbitrary initial age).",
      n_excl_founder))

  eligible <- dead_mother & born_in_sim
  offs_cal_birth   <- offs_cal_birth[eligible]
  mother_cal_birth <- mother_cal_birth[eligible]

  if (length(offs_cal_birth) == 0L) {
    warning("No eligible offspring remain after filtering (dead + born-in-sim mothers). ",
            "Consider extending year_spinup further back or shifting the gen-time window later.")
    return(list(T_g = NA_real_, mother_ages = numeric(0), n_offspring = 0L))
  }

  mother_age <- offs_cal_birth - mother_cal_birth
  valid      <- !is.na(mother_age) & mother_age > 0L

  list(
    T_g         = mean(mother_age[valid]),
    mother_ages = mother_age[valid],
    n_offspring = sum(valid)
  )
}


# ============================================================================
# Assessment settings and run IDs
# ============================================================================
# The settings that define an assessment, in one list. Used both to run the
# assessment and to build the cache ID, so the two can never disagree.
iucn_assess_spec <- function(year_spinup        = 1600L,
                             year_gentime_start = 1810L,
                             year_gentime_end   = 1910L,
                             year_rust          = 2010L,
                             year_assess        = 2026L,
                             year_end           = 2210L,
                             avg_years          = 11L,
                             max_future         = 100L) {
  list(criterion          = IUCN_CRITERION,
       year_spinup        = as.integer(year_spinup),
       year_gentime_start = as.integer(year_gentime_start),
       year_gentime_end   = as.integer(year_gentime_end),
       year_rust          = as.integer(year_rust),
       year_assess        = as.integer(year_assess),
       year_end           = as.integer(year_end),
       avg_years          = as.integer(avg_years),
       max_future         = as.integer(max_future))
}

# Cache/file ID for one assessment run: hash of params + assessment settings.
iucn_run_id <- function(params, assess = iucn_assess_spec()) {
  make_run_id(list(params = params, assess = assess))
}

# IUCN category from a percent decline (A2/A3/A4 thresholds: 30/50/80%).
iucn_category <- function(pct_decline) {
  if (is.na(pct_decline))    return(NA_character_)
  if (pct_decline >= 80)     return("CR")
  if (pct_decline >= 50)     return("EN")
  if (pct_decline >= 30)     return("VU")
  "LC/NT"
}

# ============================================================================
# assess_a4()
# ============================================================================
# Pure function: census + generation time -> A4 decline over the fixed
# window anchored on rust arrival. No simulation.
#
# n_col picks the census count: "N_IUCN" (mature individuals: alive plants of
# reproductive age, including those suppressed from flowering by fire or
# rust; the headline metric) or "N_flowering" (plants flowering that year).
assess_a4 <- function(census, T_g, year_rust, year_assess,
                      avg_years = 11L, max_future = 100L, n_col = "N_IUCN") {

  na_out <- function(L, ...) {
    list(window_len = L, ref_start = NA_integer_, ref_end = NA_integer_,
         cmp_start = NA_integer_, cmp_end = NA_integer_,
         N_ref = NA_real_, N_compare = NA_real_,
         pct_decline = NA_real_, iucn_cat = NA_character_)
  }

  if (!n_col %in% names(census))
    stop(sprintf("census has no column '%s'.", n_col))

  if (is.na(T_g)) {
    warning("T_g is NA; using a 100-year window.")
    L <- 100L
  } else {
    L <- as.integer(round(max(10, 3 * T_g)))
  }
  a  <- as.integer(avg_years)
  t0 <- as.integer(year_rust) - 1L

  cmp_end <- t0 + L
  if (cmp_end > year_assess + max_future) {
    warning(sprintf(
      "3*T_g window would end in %d, more than %d yr after %d; capping at %d.",
      cmp_end, max_future, year_assess, year_assess + max_future))
    cmp_end <- as.integer(year_assess + max_future)
    L       <- cmp_end - t0
  }
  if (cmp_end <= year_assess)
    warning(sprintf("Window ends in %d, not after the present (%d): ",
                    cmp_end, year_assess),
            "this is not an A4 window (it does not include the future).")

  ref_start <- t0 - a + 1L
  cmp_start <- cmp_end - a + 1L

  yrs <- census$year
  N   <- census[[n_col]]

  # After extinction run_simulation() stops writing census rows. Treat the
  # missing years as N = 0.
  if (nrow(census) > 0 && isTRUE(census$extinct[nrow(census)]) &&
      max(yrs) < cmp_end) {
    extra <- seq(max(yrs) + 1L, cmp_end)
    yrs   <- c(yrs, extra)
    N     <- c(N, rep(0, length(extra)))
  }

  win_mean <- function(y0, y1) {
    sel <- yrs >= y0 & yrs <= y1
    if (sum(sel) < (y1 - y0 + 1L)) return(NA_real_)
    mean(N[sel])
  }
  N_ref <- win_mean(ref_start, t0)
  N_cmp <- win_mean(cmp_start, cmp_end)

  if (is.na(N_ref) || is.na(N_cmp)) {
    warning(sprintf("Census does not cover %d-%d and %d-%d. ",
                    ref_start, t0, cmp_start, cmp_end),
            "Extend year_spinup / year_end.")
    return(na_out(L))
  }

  pct <- if (N_ref > 0) 100 * (N_ref - N_cmp) / N_ref else NA_real_

  list(window_len  = L,
       ref_start   = ref_start, ref_end = t0,
       cmp_start   = cmp_start, cmp_end = cmp_end,
       N_ref       = N_ref,     N_compare = N_cmp,
       pct_decline = pct,
       iucn_cat    = iucn_category(pct))
}

# ============================================================================
# run_iucn_a4()
# ============================================================================
run_iucn_a4 <- function(
    params,
    year_spinup          = 1600L,
    year_gentime_start   = 1810L,
    year_gentime_end     = 1910L,
    year_rust            = 2010L,
    year_assess          = 2026L,
    year_end             = 2210L,
    avg_years            = 11L,   # years averaged at each end of the window
    max_future           = 100L,  # A4: window may reach at most this far ahead
    run_counterfactual   = FALSE, # also run with rust disabled (slower)
    save_rds             = NULL,  # path to save .rds; NULL = don't save
    plot                 = TRUE,
    verbose              = TRUE
) {
  assess <- iucn_assess_spec(year_spinup, year_gentime_start, year_gentime_end,
                             year_rust, year_assess, year_end,
                             avg_years, max_future)

  # --- Validate settings ---------------------------------------------------
  if (params$fire_prob_annual == 0 &&
      length(params$fire_years) > 0 &&
      max(params$fire_years) < (year_end - year_spinup) * 0.4) {
    warning(
      "fire_years appears to cover only part of the simulation. ",
      "Consider using fire_prob_annual > 0 and fire_years = integer(0) ",
      "for IUCN runs spanning ", year_end - year_spinup, " years."
    )
  }

  # --- Build params for rust run -----------------------------------------
  p_rust                 <- params
  p_rust$rust_start_year <- year_rust
  p_rust$n_years         <- year_end - year_spinup
  p_rust$fire_years      <- integer(0)   # prevent conflict with fire_prob_annual

  if (verbose) message(sprintf(
    "Running %d years (spin-up %d -> end %d, rust from %d, present %d)...",
    p_rust$n_years, year_spinup, year_end, year_rust, year_assess))

  res_rust <- run_simulation(p_rust, year0 = year_spinup, verbose = FALSE)

  # --- Counterfactual (optional) -----------------------------------------
  res_cfact <- NULL
  if (run_counterfactual) {
    if (verbose) message("Running counterfactual (no rust)...")
    p_cfact                 <- p_rust
    p_cfact$rust_start_year <- Inf
    res_cfact <- run_simulation(p_cfact, year0 = year_spinup, verbose = FALSE)
  }

  # --- Generation time ---------------------------------------------------
  if (verbose) message("Computing generation time...")
  gt  <- compute_gen_time(res_rust$individuals,
                          window_start = year_gentime_start,
                          window_end   = year_gentime_end,
                          year_spinup  = year_spinup)
  T_g <- gt$T_g

  # --- A4 assessment -----------------------------------------------------
  # The window length comes from THIS run's own T_g, estimated above from
  # its pre-rust maternal ages, and is then applied to its own trajectory.
  a4 <- assess_a4(res_rust$census, T_g, year_rust, year_assess,
                  avg_years = avg_years, max_future = max_future)
  # Same window, flowering plants only (kept for comparison).
  a4_flw <- if ("N_flowering" %in% names(res_rust$census))
    suppressWarnings(assess_a4(res_rust$census, T_g, year_rust, year_assess,
                               avg_years = avg_years, max_future = max_future,
                               n_col = "N_flowering"))
  else list(N_ref = NA_real_, N_compare = NA_real_,
            pct_decline = NA_real_, iucn_cat = NA_character_)
  a4_flw <- stats::setNames(a4_flw[c("N_ref", "N_compare", "pct_decline", "iucn_cat")],
                            c("N_ref_flowering", "N_compare_flowering",
                              "pct_decline_flowering", "iucn_cat_flowering"))

  if (verbose) message(sprintf(
    "T_g = %.1f yr | window = %d yr | N_ref (%d-%d) = %.0f -> N_compare (%d-%d) = %.0f | decline = %.1f%% [%s]",
    T_g, a4$window_len, a4$ref_start, a4$ref_end, a4$N_ref,
    a4$cmp_start, a4$cmp_end, a4$N_compare, a4$pct_decline, a4$iucn_cat))
  if (verbose) message(sprintf(
    "  (flowering plants only: %.0f -> %.0f | decline = %.1f%% [%s])",
    a4_flw$N_ref_flowering, a4_flw$N_compare_flowering,
    a4_flw$pct_decline_flowering, a4_flw$iucn_cat_flowering))

  # --- Assemble result ---------------------------------------------------
  result <- c(
    list(criterion = IUCN_CRITERION, T_g = T_g),
    a4,
    a4_flw,
    assess[setdiff(names(assess), "criterion")],
    list(
      mother_ages    = gt$mother_ages,
      n_offspring_gt = gt$n_offspring,
      census         = res_rust$census,
      census_cfact   = if (!is.null(res_cfact)) res_cfact$census else NULL,
      run_id         = iucn_run_id(params, assess),
      params         = params
    )
  )

  # --- Save and plot -------------------------------------------------------
  if (!is.null(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, save_rds)
    if (verbose) message("Result saved: ", save_rds)
  }

  if (plot) plot_iucn_a4(result)

  invisible(result)
}

# ============================================================================
# summarise_iucn_a4()
# ============================================================================
# One-row data.frame of the key scalars, for batch accumulation.
# Threshold flags (1 = met): iucn_vu >= 30%, iucn_en >= 50%, iucn_cr >= 80%.
.or_na <- function(x) if (is.null(x)) NA_real_ else x

summarise_iucn_a4 <- function(result) {
  data.frame(
    run_id         = result$run_id,
    criterion      = result$criterion,
    T_g            = result$T_g,
    window_len     = result$window_len,
    year_assess    = result$year_assess,
    ref_start      = result$ref_start,
    ref_end        = result$ref_end,
    cmp_start      = result$cmp_start,
    cmp_end        = result$cmp_end,
    N_ref          = result$N_ref,
    N_compare      = result$N_compare,
    pct_decline    = result$pct_decline,
    iucn_vu        = as.integer(result$pct_decline >= 30),
    iucn_en        = as.integer(result$pct_decline >= 50),
    iucn_cr        = as.integer(result$pct_decline >= 80),
    n_offspring_gt = result$n_offspring_gt,
    N_ref_flowering       = .or_na(result$N_ref_flowering),
    N_compare_flowering   = .or_na(result$N_compare_flowering),
    pct_decline_flowering = .or_na(result$pct_decline_flowering),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# collect_iucn_results()
# ============================================================================
# Merge all A4 result .rds files from one or more directories into a
# data.frame (one row per run, deduplicated by run_id). Files from other
# criteria (e.g. old A3 runs) are skipped with a message.
collect_iucn_results <- function(results_dirs) {
  if (is.character(results_dirs)) results_dirs <- as.list(results_dirs)

  files <- unlist(lapply(results_dirs, function(d)
    list.files(d, pattern = "\\.rds$", full.names = TRUE)))

  if (length(files) == 0L) {
    warning("No .rds files found.")
    return(data.frame())
  }

  n_other <- 0L
  rows <- lapply(files, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || !is.list(x) || is.null(x$pct_decline)) return(NULL)
    if (!identical(x$criterion, IUCN_CRITERION)) {
      n_other <<- n_other + 1L
      return(NULL)
    }
    summarise_iucn_a4(x)
  })

  if (n_other > 0L)
    message(sprintf("Skipped %d result file(s) not from criterion %s.",
                    n_other, IUCN_CRITERION))

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(data.frame())
  out  <- do.call(rbind, rows)
  out  <- out[!duplicated(out$run_id), ]
  rownames(out) <- NULL
  message(sprintf("Collected %d unique IUCN %s runs.", nrow(out), IUCN_CRITERION))
  out
}

# ============================================================================
# plot_iucn_a4()
# ============================================================================
# Two panels:
#   1. N_IUCN through time, with the reference window (green, 11 yr ending
#      the year before rust) and the comparison window (orange, 11 yr ending
#      3 generations later), their means as dashed lines, and the decline.
#   2. Maternal-age histogram used to estimate the generation time.
plot_iucn_a4 <- function(result, save_path = NULL,
                         width_in = 12, height_in = 8, res = 150) {

  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_iucn(result)
  } else {
    .base_iucn(result)
  }
}

.fmt_span <- function(a, b) if (is.na(a)) "NA" else sprintf("%d-%d", a, b)

# --- ggplot2 rendering ------------------------------------------------------
.gg_iucn <- function(result) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes

  cen     <- result$census
  yr_rust <- result$year_rust
  yr_now  <- result$year_assess
  T_g     <- result$T_g
  pct     <- result$pct_decline
  cat     <- result$iucn_cat
  ymax    <- max(cen$N_IUCN, na.rm = TRUE) * 1.15

  # --- Panel 1: N_IUCN with reference and comparison windows ---------------
  p1 <- gg(cen, aes(x = year, y = N_IUCN)) +
    ggplot2::annotate("rect",
      xmin = result$year_spinup, xmax = result$year_gentime_start,
      ymin = 0, ymax = ymax, fill = "#dddddd", alpha = 0.4) +
    ggplot2::annotate("text",
      x = (result$year_spinup + result$year_gentime_start) / 2,
      y = ymax * 0.97, label = "spin-up",
      size = 2.8, colour = "#999999", vjust = 1) +
    ggplot2::annotate("rect",
      xmin = result$year_gentime_start, xmax = result$year_gentime_end,
      ymin = 0, ymax = ymax, fill = "#c6dbef", alpha = 0.35) +
    ggplot2::annotate("text",
      x = (result$year_gentime_start + result$year_gentime_end) / 2,
      y = ymax * 0.97,
      label = paste0("gen-time\nwindow\nT_g=", round(T_g, 1), " yr"),
      size = 2.8, colour = "#2a78d6", vjust = 1)

  if (!is.na(result$ref_start)) {
    p1 <- p1 +
      ggplot2::annotate("rect",
        xmin = result$ref_start, xmax = result$ref_end,
        ymin = 0, ymax = ymax, fill = "#a1d99b", alpha = 0.5) +
      ggplot2::annotate("rect",
        xmin = result$cmp_start, xmax = result$cmp_end,
        ymin = 0, ymax = ymax, fill = "#fdae6b", alpha = 0.5)
  }

  p1 <- p1 +
    { if ("N_flowering" %in% names(cen))
        ggplot2::geom_line(aes(y = N_flowering), colour = "#7f7f7f",
                           linewidth = 0.3, alpha = 0.7) } +
    ggplot2::geom_line(colour = "#333333", linewidth = 0.7) +
    { fires <- cen[cen$fire_event, ]
      if (nrow(fires) > 0)
        ggplot2::geom_rug(data = fires, aes(x = year), inherit.aes = FALSE,
                          colour = "#e08000", linewidth = 0.5,
                          sides = "b", length = grid::unit(0.04, "npc"))
    } +
    { if (!is.null(result$census_cfact))
        ggplot2::geom_line(data = result$census_cfact,
                           aes(x = year, y = N_IUCN),
                           colour = "#2a78d6", linewidth = 0.6,
                           linetype = "dashed")
    } +
    ggplot2::geom_vline(xintercept = yr_rust, colour = "#e34948",
                        linewidth = 0.9, linetype = "dashed") +
    ggplot2::annotate("text", x = yr_rust + 1, y = ymax * 0.88,
                      label = paste0("Rust\n", yr_rust),
                      colour = "#e34948", hjust = 0, size = 3) +
    ggplot2::geom_vline(xintercept = yr_now, colour = "#4a3aa7",
                        linewidth = 0.9) +
    ggplot2::annotate("text", x = yr_now + 1, y = ymax * 0.74,
                      label = paste0("Present\n", yr_now),
                      colour = "#4a3aa7", hjust = 0, size = 3)

  if (!is.na(result$ref_start) && !is.na(pct)) {
    p1 <- p1 +
      ggplot2::geom_hline(yintercept = result$N_ref, colour = "#1a7a35",
                          linewidth = 0.6, linetype = "longdash") +
      ggplot2::geom_hline(yintercept = result$N_compare, colour = "#b85c00",
                          linewidth = 0.6, linetype = "longdash") +
      ggplot2::annotate("segment",
        x = result$cmp_end, xend = result$cmp_end,
        y = result$N_ref, yend = result$N_compare,
        colour = "#555555", linewidth = 0.6,
        arrow = grid::arrow(length = grid::unit(0.08, "inches"))) +
      ggplot2::annotate("text",
        x = result$cmp_end + 2, y = (result$N_ref + result$N_compare) / 2,
        label = sprintf("%.0f%% decline\n[%s]", pct, cat),
        hjust = 0, size = 3.5, colour = "#333333")
  }

  p1 <- p1 +
    ggplot2::scale_y_continuous(limits = c(0, ymax),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = "Year", y = "Mature individuals (N_IUCN)",
      title = sprintf(
        "IUCN A4  |  T_g = %.1f yr  |  3*T_g = %d yr  |  %.0f%% decline -> %s",
        T_g, result$window_len, pct, cat),
      subtitle = sprintf(
        "Reference %s (mean %.0f, green)   vs   comparison %s (mean %.0f, orange)   |   thin grey = flowering only   |   rug = fires",
        .fmt_span(result$ref_start, result$ref_end), result$N_ref,
        .fmt_span(result$cmp_start, result$cmp_end), result$N_compare)) +
    ggplot2::theme_minimal(base_size = 11)

  # --- Panel 2: maternal age histogram (gen time) --------------------------
  ma_df <- data.frame(age = result$mother_ages)
  p2 <- gg(ma_df, aes(x = age)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "#2a78d6", alpha = 0.7,
                            colour = "white") +
    ggplot2::geom_vline(xintercept = T_g, colour = "#e34948",
                        linewidth = 1, linetype = "dashed") +
    ggplot2::annotate("text", x = T_g + 0.5, y = Inf,
                      label = sprintf("T_g = %.1f yr", T_g),
                      hjust = 0, vjust = 1.4, size = 3.5, colour = "#e34948") +
    ggplot2::labs(
      x = "Maternal age at offspring birth (years)", y = "Count",
      title = "Generation time",
      subtitle = sprintf("n = %d offspring  |  window %d-%d",
                         result$n_offspring_gt,
                         result$year_gentime_start, result$year_gentime_end)) +
    ggplot2::theme_minimal(base_size = 11)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(p1, p2, nrow = 2, heights = c(2.5, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, nrow = 2, heights = c(2.5, 1))
  } else {
    print(p1); readline("Enter for gen-time panel..."); print(p2)
  }
}

# --- base R fallback --------------------------------------------------------
.base_iucn <- function(result) {
  op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  cen  <- result$census
  ymax <- max(cen$N_IUCN, na.rm = TRUE) * 1.15
  plot(cen$year, cen$N_IUCN, type = "l", lwd = 1.5, col = "#333333",
       ylim = c(0, ymax), xlab = "Year",
       ylab = "Mature individuals (N_IUCN)",
       main = sprintf("IUCN A4: %.0f%% decline [%s]  ref %s vs %s  T_g=%.1f yr",
                      result$pct_decline, result$iucn_cat,
                      .fmt_span(result$ref_start, result$ref_end),
                      .fmt_span(result$cmp_start, result$cmp_end),
                      result$T_g), las = 1)
  rect(result$year_spinup, 0, result$year_gentime_start, ymax,
       col = grDevices::adjustcolor("#dddddd", 0.4), border = NA)
  rect(result$year_gentime_start, 0, result$year_gentime_end, ymax,
       col = grDevices::adjustcolor("#c6dbef", 0.4), border = NA)
  if (!is.na(result$ref_start)) {
    rect(result$ref_start, 0, result$ref_end, ymax,
         col = grDevices::adjustcolor("#a1d99b", 0.5), border = NA)
    rect(result$cmp_start, 0, result$cmp_end, ymax,
         col = grDevices::adjustcolor("#fdae6b", 0.5), border = NA)
    abline(h = result$N_ref,     col = "#1a7a35", lty = 5)
    abline(h = result$N_compare, col = "#b85c00", lty = 5)
  }
  lines(cen$year, cen$N_IUCN, col = "#333333", lwd = 1.5)
  fires <- cen[cen$fire_event, ]
  if (nrow(fires) > 0) rug(fires$year, col = "#e08000", ticksize = 0.04, lwd = 1)
  abline(v = result$year_rust,   col = "#e34948", lty = 2, lwd = 1.5)
  abline(v = result$year_assess, col = "#4a3aa7", lwd = 1.5)

  hist(result$mother_ages, breaks = 30, col = "#2a78d6", border = "white",
       xlab = "Maternal age at offspring birth (years)",
       main = sprintf("Generation time  (T_g = %.1f yr,  n = %d offspring)",
                      result$T_g, result$n_offspring_gt))
  abline(v = result$T_g, col = "#e34948", lty = 2, lwd = 2)
}

# ============================================================================
# run_iucn_seeds()
# ============================================================================
# Run the A4 assessment for a fixed parameterisation across many seeds.
# Each run gives one decline (reference vs comparison window). Fire timing
# is stochastic, so report the central tendency (median) and spread across
# seeds. One .rds per seed is written to results_dir, named by
# iucn_run_id(), so the ensemble is restartable and extendable.
#
# USAGE
#   df <- run_iucn_seeds(get_default_params(), seeds = 1:50,
#                        results_dir = "results/iucn_seeds_a4/", n_cores = 8L,
#                        year_spinup = 1600L, year_assess = 2026L, year_end = 2210L)
#   plot_iucn_ensemble(df)
run_iucn_seeds <- function(
    params,
    seeds         = 1:20,
    results_dir   = "results/iucn_seeds_a4/",
    n_cores       = max(1L, parallel::detectCores() - 1L),
    skip_existing = TRUE,
    verbose       = TRUE,
    ...           # passed to run_iucn_a4 (year_spinup, year_assess, etc.)
) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  extra  <- list(...)
  assess <- do.call(iucn_assess_spec,
                    extra[intersect(names(extra), names(formals(iucn_assess_spec)))])

  n <- length(seeds)
  if (verbose) message(sprintf("IUCN %s seed ensemble: %d runs on %d core(s)",
                               IUCN_CRITERION, n, n_cores))

  parallel::mclapply(seeds, function(s) {
    p      <- params
    p$seed <- as.integer(s)
    out    <- file.path(results_dir, paste0(iucn_run_id(p, assess), ".rds"))
    if (skip_existing && file.exists(out)) return(invisible(NULL))
    tryCatch(
      run_iucn_a4(p, ..., save_rds = out, plot = FALSE, verbose = FALSE),
      error = function(e) message(sprintf("  seed %d error: %s", s, conditionMessage(e)))
    )
    invisible(NULL)
  }, mc.cores = n_cores)

  df <- collect_iucn_results(results_dir)
  if (verbose && nrow(df) > 0) {
    message(sprintf(
      "Ensemble summary (%d runs):  median decline = %.1f%%  [10-90%%: %.1f%% - %.1f%%]",
      nrow(df),
      median(df$pct_decline, na.rm = TRUE),
      quantile(df$pct_decline, 0.1, na.rm = TRUE),
      quantile(df$pct_decline, 0.9, na.rm = TRUE)))
  }
  df
}

# ============================================================================
# plot_iucn_ensemble()
# ============================================================================
# Panel A: histogram of A4 decline across runs, with VU/EN/CR thresholds and
#          the median (the central tendency used for the assessment).
# Panel B: N_ref and N_compare across runs.
plot_iucn_ensemble <- function(df, save_path = NULL,
                               width_in = 10, height_in = 5, res = 150) {
  if (nrow(df) == 0L) stop("No results to plot.")

  n   <- nrow(df)
  med <- median(df$pct_decline, na.rm = TRUE)
  cat(sprintf("\n=== IUCN %s ensemble: %d runs ===\n", IUCN_CRITERION, n))
  cat(sprintf("  pct_decline: median=%.1f%% [%s]  IQR=[%.1f%%, %.1f%%]  range=[%.1f%%, %.1f%%]\n",
              med, iucn_category(med),
              quantile(df$pct_decline, 0.25, na.rm = TRUE),
              quantile(df$pct_decline, 0.75, na.rm = TRUE),
              min(df$pct_decline, na.rm = TRUE),
              max(df$pct_decline, na.rm = TRUE)))
  cat(sprintf("  Outcome frequencies:  LC/NT=%d%%  VU=%d%%  EN=%d%%  CR=%d%%\n",
              round(100 * mean(df$pct_decline <  30, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 30 & df$pct_decline < 50, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 50 & df$pct_decline < 80, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 80, na.rm = TRUE))))
  cat(sprintf("  T_g: median=%.1f yr   3*T_g window: median=%.0f yr\n",
              median(df$T_g, na.rm = TRUE), median(df$window_len, na.rm = TRUE)))

  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_iucn_ensemble(df)
  } else {
    .base_iucn_ensemble(df)
  }

  cats <- factor(vapply(df$pct_decline, iucn_category, ""),
                 levels = c("LC/NT", "VU", "EN", "CR"))
  invisible(list(
    n              = n,
    category_share = as.list(prop.table(table(cats))),
    median_decline = med,
    median_cat     = iucn_category(med),
    q10            = unname(quantile(df$pct_decline, 0.1, na.rm = TRUE)),
    q90            = unname(quantile(df$pct_decline, 0.9, na.rm = TRUE)),
    pct_vu         = mean(df$iucn_vu, na.rm = TRUE),
    pct_en         = mean(df$iucn_en, na.rm = TRUE),
    pct_cr         = mean(df$iucn_cr, na.rm = TRUE)
  ))
}

.gg_iucn_ensemble <- function(df) {
  gg  <- ggplot2::ggplot; aes <- ggplot2::aes
  n   <- nrow(df)
  med <- median(df$pct_decline, na.rm = TRUE)

  cat_lv  <- c("LC/NT", "VU", "EN", "CR")
  cat_col <- c("LC/NT" = "#9ecae1", VU = "#eda100", EN = "#e34948", CR = "#8b0000")
  df$cat  <- factor(vapply(df$pct_decline, iucn_category, ""), levels = cat_lv)
  shares  <- as.integer(round(100 * table(df$cat) / max(1, sum(!is.na(df$cat)))))

  pA <- gg(df, aes(x = pct_decline, fill = cat)) +
    ggplot2::geom_histogram(binwidth = 2, boundary = 0, colour = "white",
                            alpha = 0.9, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = cat_col, drop = FALSE, name = NULL,
                               labels = sprintf("%s  %d%%", cat_lv, shares)) +
    ggplot2::geom_vline(xintercept = c(30, 50, 80),
                        colour = c("#eda100", "#e34948", "#8b0000"),
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::annotate("text", x = c(30, 50, 80), y = Inf,
                      label = c("VU", "EN", "CR"),
                      hjust = -0.1, vjust = 1.4, size = 3,
                      colour = c("#eda100", "#e34948", "#8b0000")) +
    ggplot2::geom_vline(xintercept = med, colour = "#333333", linewidth = 1) +
    ggplot2::annotate("text", x = med, y = Inf,
                      label = sprintf(" median\n %.1f%%", med),
                      hjust = 0, vjust = 1.3, size = 3, colour = "#333333") +
    # Declines can be negative (N grew); keep those runs on the axis.
    ggplot2::scale_x_continuous(
      limits = c(min(0, 10 * floor(min(df$pct_decline, na.rm = TRUE) / 10)), 100)) +
    ggplot2::labs(x = "Decline in mature individuals, reference -> comparison (%)",
                  y = "Number of runs",
                  title = sprintf(
                    "IUCN A4 ensemble (n = %d)  |  LC/NT %d%%  VU %d%%  EN %d%%  CR %d%%",
                    n, shares[1], shares[2], shares[3], shares[4]),
                  subtitle = sprintf("Share of runs by category; median decline %.0f%% (%s)",
                                     med, iucn_category(med))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  lv <- c("Reference (pre-rust)", "Comparison (3 gen on)")
  pop_long <- rbind(
    data.frame(window = lv[1], N = df$N_ref,     stringsAsFactors = FALSE),
    data.frame(window = lv[2], N = df$N_compare, stringsAsFactors = FALSE))
  pop_long$window <- factor(pop_long$window, levels = lv)

  pB <- gg(pop_long, aes(x = window, y = N, fill = window)) +
    ggplot2::geom_violin(alpha = 0.6, colour = NA, na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.size = 1,
                          colour = "#333333", na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("#a1d99b", "#fdae6b"), guide = "none") +
    ggplot2::labs(x = NULL, y = "Mature individuals (N_IUCN, 11-yr mean)",
                  title = "Reference vs comparison",
                  subtitle = "One value per run") +
    ggplot2::theme_minimal(base_size = 11)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1.6, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, ncol = 2, widths = c(1.6, 1))
  } else {
    print(pA); readline("Enter for Panel B..."); print(pB)
  }
}

.base_iucn_ensemble <- function(df) {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  med <- median(df$pct_decline, na.rm = TRUE)
  lo  <- min(0, 5 * floor(min(df$pct_decline, na.rm = TRUE) / 5))
  hist(df$pct_decline, breaks = seq(lo, 100, by = 5),
       col = "#2a78d6", border = "white",
       xlab = "Decline, reference to comparison (%)", ylab = "Runs",
       main = sprintf("IUCN A4 ensemble (n=%d), median %.0f%%", nrow(df), med), las = 1)
  abline(v = c(30, 50, 80), col = c("#eda100", "#e34948", "#8b0000"), lty = 2, lwd = 1.5)
  abline(v = med, col = "#333333", lwd = 2)

  boxplot(list("Reference\n(pre-rust)" = df$N_ref, "Comparison\n(3 gen on)" = df$N_compare),
          col = c("#a1d99b", "#fdae6b"), border = "#333333",
          ylab = "Mature individuals", main = "Reference vs comparison", las = 1)
}
