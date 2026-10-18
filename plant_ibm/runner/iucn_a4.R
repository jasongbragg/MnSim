# runner/iucn_a4.R
#
# IUCN Criterion A4 assessment: the largest decline in reproductive adults
# (N_IUCN) over ANY window of 10 years or three generations, whichever is
# longer, where the window must include both the past and the future and
# may extend at most 100 years into the future.
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
#   L           = max(10, 3 * T_g) years, rounded to whole years
#   year_assess = the present (assessment year)
#   Candidate windows [s, s + L] with
#       s     <  year_assess          (window includes the past)
#       s + L >  year_assess          (window includes the future)
#       s + L <= year_assess + 100    (at most 100 years into the future)
#   N(y)      = centred running mean of N_IUCN over [y - h, y + h]
#               (h = smooth_half; dampens single-year fire dips)
#   decline(s) = 100 * (N(s) - N(s + L)) / N(s)
#   A4 decline = max over candidate s.
# If the population goes extinct, years after the last census row count
# as N = 0, so windows ending after extinction score 100% decline.
#
# NOTE ON MAXIMISING OVER WINDOWS
# -------------------------------
# Taking the maximum over many overlapping windows of a noisy trajectory
# selects windows that start on a high point and end on a low point, which
# biases a single-run decline upward. The running mean (smooth_half) limits
# this. The full decline-by-start-year curve is returned (decline_curve)
# so the sensitivity can be inspected, and ensembles should be summarised
# across seeds rather than read from one run.
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
#   year_rust          rust arrives (2010)        ─┤
#   year_assess        the present (2026)         ─┤ A4 windows straddle this
#   year_end           end of run                 ─┘ must be >= year_assess + 100 + h
#
# FIRE NOTE
# ---------
# Use fire_prob_annual > 0 and fire_years = integer(0) for these long runs.
#
# CACHING / RUN IDs
# -----------------
# iucn_run_id(params, assess) hashes the params AND the assessment settings
# (years, smoothing, criterion). Cached .rds files therefore never get
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
                             smooth_half        = 5L,
                             max_future         = 100L) {
  list(criterion          = IUCN_CRITERION,
       year_spinup        = as.integer(year_spinup),
       year_gentime_start = as.integer(year_gentime_start),
       year_gentime_end   = as.integer(year_gentime_end),
       year_rust          = as.integer(year_rust),
       year_assess        = as.integer(year_assess),
       year_end           = as.integer(year_end),
       smooth_half        = as.integer(smooth_half),
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
# Pure function: census + generation time -> A4 decline. No simulation.
# Returns list(window_len, win_start, win_end, N_start, N_end, pct_decline,
#              iucn_cat, decline_curve).
assess_a4 <- function(census, T_g, year_assess,
                      smooth_half = 5L, max_future = 100L) {

  if (is.na(T_g)) {
    warning("T_g is NA; using a 100-year window.")
    L <- 100L
  } else {
    L <- as.integer(round(max(10, 3 * T_g)))
  }
  h <- as.integer(smooth_half)

  yrs <- census$year
  N   <- census$N_IUCN

  # After extinction run_simulation() stops writing census rows. Treat the
  # missing years as N = 0 so post-extinction windows score 100% decline.
  if (nrow(census) > 0 && isTRUE(census$extinct[nrow(census)])) {
    pad_to <- year_assess + max_future + h
    if (max(yrs) < pad_to) {
      extra <- seq(max(yrs) + 1L, pad_to)
      yrs   <- c(yrs, extra)
      N     <- c(N, rep(0, length(extra)))
    }
  }

  # Centred running mean; NA unless the full (2h + 1)-year window exists.
  n_bar <- function(y) {
    sel <- yrs >= y - h & yrs <= y + h
    if (sum(sel) < 2L * h + 1L) return(NA_real_)
    mean(N[sel])
  }

  s_min <- year_assess - L + 1L
  s_max <- min(year_assess - 1L, year_assess + max_future - L)
  s     <- seq(s_min, s_max)

  N_s <- vapply(s,     n_bar, numeric(1))
  N_e <- vapply(s + L, n_bar, numeric(1))
  pct <- ifelse(!is.na(N_s) & N_s > 0, 100 * (N_s - N_e) / N_s, NA_real_)

  curve <- data.frame(win_start = s, win_end = s + L,
                      N_start = N_s, N_end = N_e, pct_decline = pct)

  if (all(is.na(pct))) {
    warning("No A4 window could be evaluated: census does not cover ",
            s_min - h, "-", s_max + L + h, ". Extend year_spinup / year_end.")
    return(list(window_len = L, win_start = NA_integer_, win_end = NA_integer_,
                N_start = NA_real_, N_end = NA_real_, pct_decline = NA_real_,
                iucn_cat = NA_character_, decline_curve = curve))
  }
  if (any(is.na(pct)))
    warning(sum(is.na(pct)), " of ", length(pct),
            " candidate A4 windows could not be evaluated (census too short).")

  k <- which.max(pct)
  list(window_len    = L,
       win_start     = curve$win_start[k],
       win_end       = curve$win_end[k],
       N_start       = curve$N_start[k],
       N_end         = curve$N_end[k],
       pct_decline   = curve$pct_decline[k],
       iucn_cat      = iucn_category(curve$pct_decline[k]),
       decline_curve = curve)
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
    smooth_half          = 5L,    # running-mean half width for N(y)
    max_future           = 100L,  # A4: window may reach at most this far ahead
    run_counterfactual   = FALSE, # also run with rust disabled (slower)
    save_rds             = NULL,  # path to save .rds; NULL = don't save
    plot                 = TRUE,
    verbose              = TRUE
) {
  assess <- iucn_assess_spec(year_spinup, year_gentime_start, year_gentime_end,
                             year_rust, year_assess, year_end,
                             smooth_half, max_future)

  # --- Validate settings ---------------------------------------------------
  if (year_end < year_assess + max_future + smooth_half)
    warning("year_end (", year_end, ") < year_assess + max_future + smooth_half (",
            year_assess + max_future + smooth_half,
            "): some A4 windows cannot be evaluated.")
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
  a4 <- assess_a4(res_rust$census, T_g, year_assess,
                  smooth_half = smooth_half, max_future = max_future)

  if (verbose) message(sprintf(
    "T_g = %.1f yr | window = %d yr | worst window %d-%d | N %.0f -> %.0f | decline = %.1f%% [%s]",
    T_g, a4$window_len, a4$win_start, a4$win_end,
    a4$N_start, a4$N_end, a4$pct_decline, a4$iucn_cat))

  # --- Assemble result ---------------------------------------------------
  result <- c(
    list(criterion = IUCN_CRITERION, T_g = T_g),
    a4,
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
summarise_iucn_a4 <- function(result) {
  data.frame(
    run_id         = result$run_id,
    criterion      = result$criterion,
    T_g            = result$T_g,
    window_len     = result$window_len,
    year_assess    = result$year_assess,
    win_start      = result$win_start,
    win_end        = result$win_end,
    N_start        = result$N_start,
    N_end          = result$N_end,
    pct_decline    = result$pct_decline,
    iucn_vu        = as.integer(result$pct_decline >= 30),
    iucn_en        = as.integer(result$pct_decline >= 50),
    iucn_cr        = as.integer(result$pct_decline >= 80),
    n_offspring_gt = result$n_offspring_gt,
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
plot_iucn_a4 <- function(result, save_path = NULL,
                         width_in = 12, height_in = 10, res = 150) {

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

# --- ggplot2 rendering ------------------------------------------------------
.gg_iucn <- function(result) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes

  cen     <- result$census
  h       <- result$smooth_half
  yr_rust <- result$year_rust
  yr_now  <- result$year_assess
  T_g     <- result$T_g
  pct     <- result$pct_decline
  cat     <- result$iucn_cat
  ymax    <- max(cen$N_IUCN, na.rm = TRUE) * 1.15

  # --- Panel 1: N_IUCN time series with the worst A4 window ---------------
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
      label = paste0("gen-time\nwindow\nTᴳ=", round(T_g, 1), " yr"),
      size = 2.8, colour = "#2a78d6", vjust = 1)

  if (!is.na(result$win_start)) {
    p1 <- p1 +
      # Worst A4 window
      ggplot2::annotate("rect",
        xmin = result$win_start, xmax = result$win_end,
        ymin = 0, ymax = ymax, fill = "#fee6ce", alpha = 0.5) +
      # Averaging bands at each end of the window
      ggplot2::annotate("rect",
        xmin = result$win_start - h, xmax = result$win_start + h,
        ymin = 0, ymax = ymax, fill = "#a1d99b", alpha = 0.5) +
      ggplot2::annotate("rect",
        xmin = result$win_end - h, xmax = result$win_end + h,
        ymin = 0, ymax = ymax, fill = "#fdae6b", alpha = 0.5)
  }

  p1 <- p1 +
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

  if (!is.na(result$win_start)) {
    p1 <- p1 +
      ggplot2::annotate("segment",
        x = result$win_start, xend = result$win_end,
        y = result$N_start, yend = result$N_end,
        colour = "#888780", linewidth = 0.6, linetype = "dotted") +
      ggplot2::annotate("point", x = c(result$win_start, result$win_end),
                        y = c(result$N_start, result$N_end),
                        colour = c("#1a7a35", "#b85c00"), size = 2.5) +
      ggplot2::annotate("text",
        x = result$win_end + 2, y = (result$N_start + result$N_end) / 2,
        label = sprintf("%.0f%% decline\n[%s]", pct, cat),
        hjust = 0, size = 3.5, colour = "#333333")
  }

  p1 <- p1 +
    ggplot2::scale_y_continuous(limits = c(0, ymax),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = "Year", y = "Reproductive adults (N_IUCN)",
      title = sprintf(
        "IUCN A4  |  Tᴳ = %.1f yr  |  window = %d yr  |  worst window %s  |  %.0f%% decline → %s",
        T_g, result$window_len,
        if (is.na(result$win_start)) "NA"
        else sprintf("%d–%d", result$win_start, result$win_end),
        pct, cat),
      subtitle = sprintf(
        "Window start/end values are %d-yr running means (green / orange bands)   |   rug = fire events",
        2L * h + 1L)) +
    ggplot2::theme_minimal(base_size = 11)

  # --- Panel 2: decline for every candidate window -------------------------
  dc <- result$decline_curve
  p2 <- gg(dc, aes(x = win_start, y = pct_decline)) +
    ggplot2::geom_hline(yintercept = c(30, 50, 80),
                        colour = c("#eda100", "#e34948", "#8b0000"),
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_line(colour = "#333333", linewidth = 0.7, na.rm = TRUE) +
    { if (!is.na(result$win_start))
        ggplot2::annotate("point", x = result$win_start, y = pct,
                          colour = "#e34948", size = 2.5) } +
    ggplot2::geom_vline(xintercept = yr_rust, colour = "#e34948",
                        linetype = "dashed", linewidth = 0.6) +
    # Windows starting in a post-fire trough can show large negative
    # "declines"; clip the display so the threshold region stays readable.
    ggplot2::coord_cartesian(ylim = c(-100, 100)) +
    ggplot2::labs(x = "Window start year",
                  y = "Decline over window (%)",
                  title = "Decline for each A4 window (all straddle the present)",
                  subtitle = "Red point = maximum, used for the assessment") +
    ggplot2::theme_minimal(base_size = 11)

  # --- Panel 3: maternal age histogram (gen time) --------------------------
  ma_df <- data.frame(age = result$mother_ages)
  p3 <- gg(ma_df, aes(x = age)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "#2a78d6", alpha = 0.7,
                            colour = "white") +
    ggplot2::geom_vline(xintercept = T_g, colour = "#e34948",
                        linewidth = 1, linetype = "dashed") +
    ggplot2::annotate("text", x = T_g + 0.5, y = Inf,
                      label = sprintf("Tᴳ = %.1f yr", T_g),
                      hjust = 0, vjust = 1.4, size = 3.5, colour = "#e34948") +
    ggplot2::labs(
      x = "Maternal age at offspring birth (years)", y = "Count",
      title = "Generation time",
      subtitle = sprintf("n = %d offspring  |  window %d–%d",
                         result$n_offspring_gt,
                         result$year_gentime_start, result$year_gentime_end)) +
    ggplot2::theme_minimal(base_size = 11)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(p1, p2, p3, nrow = 3, heights = c(2.5, 1, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, p3, nrow = 3, heights = c(2.5, 1, 1))
  } else {
    print(p1); readline("Enter for next panel..."); print(p2)
    readline("Enter for next panel..."); print(p3)
  }
}

# --- base R fallback --------------------------------------------------------
.base_iucn <- function(result) {
  op <- par(mfrow = c(3, 1), mar = c(4, 4, 3, 2), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  cen <- result$census
  h   <- result$smooth_half
  ymax <- max(cen$N_IUCN, na.rm = TRUE) * 1.15
  plot(cen$year, cen$N_IUCN, type = "l", lwd = 1.5, col = "#333333",
       ylim = c(0, ymax), xlab = "Year",
       ylab = "Reproductive adults (N_IUCN)",
       main = sprintf("IUCN A4: %.0f%% decline [%s]  window %s  T_g=%.1f yr",
                      result$pct_decline, result$iucn_cat,
                      if (is.na(result$win_start)) "NA"
                      else paste0(result$win_start, "-", result$win_end),
                      result$T_g), las = 1)
  rect(result$year_spinup, 0, result$year_gentime_start, ymax,
       col = grDevices::adjustcolor("#dddddd", 0.4), border = NA)
  rect(result$year_gentime_start, 0, result$year_gentime_end, ymax,
       col = grDevices::adjustcolor("#c6dbef", 0.4), border = NA)
  if (!is.na(result$win_start)) {
    rect(result$win_start, 0, result$win_end, ymax,
         col = grDevices::adjustcolor("#fee6ce", 0.5), border = NA)
    rect(result$win_start - h, 0, result$win_start + h, ymax,
         col = grDevices::adjustcolor("#a1d99b", 0.5), border = NA)
    rect(result$win_end - h, 0, result$win_end + h, ymax,
         col = grDevices::adjustcolor("#fdae6b", 0.5), border = NA)
  }
  lines(cen$year, cen$N_IUCN, col = "#333333", lwd = 1.5)
  fires <- cen[cen$fire_event, ]
  if (nrow(fires) > 0) rug(fires$year, col = "#e08000", ticksize = 0.04, lwd = 1)
  abline(v = result$year_rust,   col = "#e34948", lty = 2, lwd = 1.5)
  abline(v = result$year_assess, col = "#4a3aa7", lwd = 1.5)

  dc <- result$decline_curve
  plot(dc$win_start, dc$pct_decline, type = "l", lwd = 1.5, ylim = c(-100, 100),
       xlab = "Window start year", ylab = "Decline (%)",
       main = "Decline for each A4 window", las = 1)
  abline(h = c(30, 50, 80), col = c("#eda100", "#e34948", "#8b0000"), lty = 2)
  if (!is.na(result$win_start))
    points(result$win_start, result$pct_decline, col = "#e34948", pch = 19)

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
# Fire timing is stochastic, so a single run's decline depends partly on
# where fires fall relative to the window ends; report the distribution
# across seeds. One .rds per seed is written to results_dir, named by
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
      "Ensemble summary (%d runs):  median decline = %.1f%%  [%.1f%% - %.1f%%]",
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
# Panel A: histogram of A4 decline with VU/EN/CR thresholds.
# Panel B: N at the start and end of each run's worst window.
plot_iucn_ensemble <- function(df, save_path = NULL,
                               width_in = 10, height_in = 5, res = 150) {
  if (nrow(df) == 0L) stop("No results to plot.")

  n <- nrow(df)
  cat(sprintf("\n=== IUCN %s ensemble: %d runs ===\n", IUCN_CRITERION, n))
  cat(sprintf("  pct_decline: median=%.1f%%  IQR=[%.1f%%, %.1f%%]  range=[%.1f%%, %.1f%%]\n",
              median(df$pct_decline, na.rm = TRUE),
              quantile(df$pct_decline, 0.25, na.rm = TRUE),
              quantile(df$pct_decline, 0.75, na.rm = TRUE),
              min(df$pct_decline, na.rm = TRUE),
              max(df$pct_decline, na.rm = TRUE)))
  cat(sprintf("  Outcome frequencies:  LC/NT=%d%%  VU=%d%%  EN=%d%%  CR=%d%%\n",
              round(100 * mean(df$pct_decline <  30, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 30 & df$pct_decline < 50, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 50 & df$pct_decline < 80, na.rm = TRUE)),
              round(100 * mean(df$pct_decline >= 80, na.rm = TRUE))))
  cat(sprintf("  T_g: median=%.1f yr   worst-window start: median=%.0f\n",
              median(df$T_g, na.rm = TRUE), median(df$win_start, na.rm = TRUE)))

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

  invisible(list(
    n              = n,
    pct_vu         = mean(df$iucn_vu, na.rm = TRUE),
    pct_en         = mean(df$iucn_en, na.rm = TRUE),
    pct_cr         = mean(df$iucn_cr, na.rm = TRUE),
    median_decline = median(df$pct_decline, na.rm = TRUE)
  ))
}

.gg_iucn_ensemble <- function(df) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  n  <- nrow(df)

  pA <- gg(df, aes(x = pct_decline)) +
    ggplot2::geom_histogram(binwidth = 2, fill = "#2a78d6", colour = "white",
                            alpha = 0.85, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(30, 50, 80),
                        colour = c("#eda100", "#e34948", "#8b0000"),
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::annotate("text", x = c(30, 50, 80), y = Inf,
                      label = c("VU", "EN", "CR"),
                      hjust = -0.1, vjust = 1.4, size = 3,
                      colour = c("#eda100", "#e34948", "#8b0000")) +
    ggplot2::geom_vline(xintercept = median(df$pct_decline, na.rm = TRUE),
                        colour = "#333333", linewidth = 1) +
    ggplot2::annotate("text",
                      x = median(df$pct_decline, na.rm = TRUE), y = Inf,
                      label = sprintf(" median\n %.1f%%",
                                      median(df$pct_decline, na.rm = TRUE)),
                      hjust = 0, vjust = 1.3, size = 3, colour = "#333333") +
    ggplot2::scale_x_continuous(limits = c(0, 100)) +
    ggplot2::labs(x = "Decline in reproductive adults over worst A4 window (%)",
                  y = "Number of runs",
                  title = sprintf("IUCN A4 ensemble  (n = %d seeds)", n),
                  subtitle = sprintf(
                    "LC/NT: %d%%   VU: %d%%   EN: %d%%   CR: %d%%",
                    round(100 * mean(df$pct_decline <  30, na.rm = TRUE)),
                    round(100 * mean(df$iucn_vu & !df$iucn_en, na.rm = TRUE)),
                    round(100 * mean(df$iucn_en & !df$iucn_cr, na.rm = TRUE)),
                    round(100 * mean(df$iucn_cr, na.rm = TRUE)))) +
    ggplot2::theme_minimal(base_size = 11)

  pop_long <- rbind(
    data.frame(window = "Window start", N = df$N_start, stringsAsFactors = FALSE),
    data.frame(window = "Window end",   N = df$N_end,   stringsAsFactors = FALSE))
  pop_long$window <- factor(pop_long$window, levels = c("Window start", "Window end"))

  pB <- gg(pop_long, aes(x = window, y = N, fill = window)) +
    ggplot2::geom_violin(alpha = 0.6, colour = NA, na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.size = 1,
                          colour = "#333333", na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("#a1d99b", "#fdae6b"), guide = "none") +
    ggplot2::labs(x = NULL, y = "Reproductive adults (N_IUCN)",
                  title = "N at ends of the worst window",
                  subtitle = "Each point = one seed run") +
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

  hist(df$pct_decline, breaks = seq(0, 100, by = 5),
       col = "#2a78d6", border = "white",
       xlab = "Decline over worst A4 window (%)", ylab = "Runs",
       main = sprintf("IUCN A4 ensemble (n=%d)", nrow(df)), las = 1)
  abline(v = c(30, 50, 80), col = c("#eda100", "#e34948", "#8b0000"), lty = 2, lwd = 1.5)
  abline(v = median(df$pct_decline, na.rm = TRUE), col = "#333333", lwd = 2)

  boxplot(list("Window\nstart" = df$N_start, "Window\nend" = df$N_end),
          col = c("#a1d99b", "#fdae6b"), border = "#333333",
          ylab = "Reproductive adults", main = "N at window ends", las = 1)
}
