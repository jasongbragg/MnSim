# runner/iucn_a3.R
#
# IUCN Criterion A3 assessment: projected decline in reproductive adults
# (N_IUCN) over the assessment period min(100 yr, 3 × generation time),
# comparing a pre-rust reference to a post-rust comparison point.
#
# GENERATION TIME
# ---------------
# Computed as the mean age of mothers of offspring born during a specified
# pre-rust window (standard IUCN definition: mean age of parents of the
# current cohort). Requires mother_id and birth_year to be tracked in the
# pop data frame -- both are recorded by make_recruit_rows() and
# create_population() in the current model.
#
# SIMULATION PHASING (calendar years, single continuous run)
# ----------------------------------------------------------
#   year_spinup      start of run (e.g. 1700)  ─┐ spin-up: discard
#   year_gentime_start                          ─┤ gen-time window: T_g computed here
#   year_gentime_end                            ─┤ pre-rust equilibrium
#   year_rust        rust arrives (e.g. 2010)  ─┤ reference window ends here
#   year_end         end of run (e.g. 2210)    ─┘ comparison window within here
#
# FIRE NOTE
# ---------
# For runs spanning 400+ years, params$fire_years (fixed simulation years)
# will only cover the early part of the run. Use fire_prob_annual > 0
# and fire_years = integer(0) for IUCN runs. run_iucn_a3() warns if
# fire_years appears to be set instead of fire_prob_annual.
#
# BATCH INTERFACE
# ---------------
# run_iucn_a3(..., save_rds = file.path(dir, paste0(run_id, ".rds")))
#   saves the full result list, identified by the run_id hash.
# summarise_iucn_a3(result) extracts a named numeric vector for
#   accumulation via collect_iucn_results() or write.csv().
# collect_iucn_results(dirs) merges all .rds files from one or more
#   result directories into a single data.frame.
#
# SINGLE-RUN USAGE
# ----------------
#   source("runner/iucn_a3.R")
#   # (source all model files and params first)
#   result <- run_iucn_a3(get_default_params())
#   plot_iucn_a3(result)
#   summarise_iucn_a3(result)

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
# run_iucn_a3()
# ============================================================================
run_iucn_a3 <- function(
    params,
    year_spinup          = 1700L,
    year_gentime_start   = 1810L,
    year_gentime_end     = 1910L,
    year_rust            = 2010L,
    year_end             = 2210L,
    ref_half_window      = 5L,    # reference: ±5 yr around (year_rust - 1)
    compare_half_window  = 5L,    # comparison: ±5 yr around comparison point
    run_counterfactual   = FALSE, # also run with rust disabled (slower)
    save_rds             = NULL,  # path to save .rds; NULL = don't save
    plot                 = TRUE,
    verbose              = TRUE
) {

  # --- Validate fire parameterisation -------------------------------------
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
  p_rust               <- params
  p_rust$rust_start_year <- year_rust
  p_rust$n_years         <- year_end - year_spinup
  p_rust$fire_years      <- integer(0)   # prevent conflict with fire_prob_annual

  if (verbose) message(sprintf(
    "Running %d years (spin-up %d → end %d, rust from %d)...",
    p_rust$n_years, year_spinup, year_end, year_rust))

  res_rust <- run_simulation(p_rust, year0 = year_spinup, verbose = FALSE)

  # --- Counterfactual (optional) -----------------------------------------
  res_cfact <- NULL
  if (run_counterfactual) {
    if (verbose) message("Running counterfactual (no rust)...")
    p_cfact               <- p_rust
    p_cfact$rust_start_year <- Inf
    res_cfact <- run_simulation(p_cfact, year0 = year_spinup, verbose = FALSE)
  }

  # --- Generation time ---------------------------------------------------
  if (verbose) message("Computing generation time...")
  gt <- compute_gen_time(res_rust$individuals,
                          window_start  = year_gentime_start,
                          window_end    = year_gentime_end,
                          year_spinup   = year_spinup)
  T_g <- gt$T_g

  # --- Assessment period -------------------------------------------------
  t_inf <- if (!is.na(T_g)) min(100, 3 * T_g) else 100

  # --- Reference N_IUCN (pre-rust window) --------------------------------
  ref_mid   <- year_rust - 1L
  ref_start <- ref_mid - ref_half_window
  ref_end   <- ref_mid + ref_half_window
  cen_ref   <- res_rust$census[res_rust$census$year >= ref_start &
                                 res_rust$census$year <= ref_end, ]
  N_ref     <- mean(cen_ref$N_IUCN, na.rm = TRUE)

  # --- Comparison N_IUCN (post-rust window) ------------------------------
  compare_mid   <- year_rust + round(t_inf)
  compare_start <- compare_mid - compare_half_window
  compare_end   <- compare_mid + compare_half_window
  cen_cmp       <- res_rust$census[res_rust$census$year >= compare_start &
                                     res_rust$census$year <= compare_end, ]
  N_compare     <- mean(cen_cmp$N_IUCN, na.rm = TRUE)

  # --- Decline and IUCN category -----------------------------------------
  pct_decline  <- 100 * (N_ref - N_compare) / N_ref
  iucn_cat <- if (pct_decline >= 80) { "CR"
              } else if (pct_decline >= 50) { "EN"
              } else if (pct_decline >= 30) { "VU"
              } else { "LC/NT" }

  if (verbose) message(sprintf(
    "T_g = %.1f yr  |  assessment period = %.0f yr  |  ",
    T_g, t_inf), appendLF = FALSE)
  if (verbose) message(sprintf(
    "N_ref = %.0f  ->  N_compare = %.0f  |  decline = %.1f%%  [%s]",
    N_ref, N_compare, pct_decline, iucn_cat))

  # --- Assemble result ---------------------------------------------------
  result <- list(
    # Assessment numbers
    T_g             = T_g,
    t_inf           = t_inf,
    N_ref           = N_ref,
    N_compare       = N_compare,
    pct_decline     = pct_decline,
    iucn_cat        = iucn_cat,
    # Window positions
    year_spinup     = year_spinup,
    year_gentime_start = year_gentime_start,
    year_gentime_end   = year_gentime_end,
    year_rust       = year_rust,
    year_end        = year_end,
    ref_start       = ref_start,
    ref_end         = ref_end,
    compare_start   = compare_start,
    compare_end     = compare_end,
    compare_mid     = compare_mid,
    # Gen time data
    mother_ages     = gt$mother_ages,
    n_offspring_gt  = gt$n_offspring,
    # Full outputs
    census          = res_rust$census,
    census_cfact    = if (!is.null(res_cfact)) res_cfact$census else NULL,
    # Identification
    run_id          = make_run_id(params),
    params          = params
  )

  # --- Save and plot -------------------------------------------------------
  if (!is.null(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, save_rds)
    if (verbose) message("Result saved: ", save_rds)
  }

  if (plot) plot_iucn_a3(result)

  invisible(result)
}

# ============================================================================
# summarise_iucn_a3()
# ============================================================================
# Extracts the key scalars as a named numeric vector for batch accumulation.
# Compatible with collect_iucn_results() and CSV output.
# IUCN thresholds flags (1 = threshold met):
#   iucn_vu: pct_decline >= 30%  (Vulnerable)
#   iucn_en: pct_decline >= 50%  (Endangered)
#   iucn_cr: pct_decline >= 80%  (Critically Endangered)
summarise_iucn_a3 <- function(result) {
  # Returns a single-row data.frame rather than a named vector.
  # A named vector would coerce all fields to character because run_id
  # is a string; data.frame preserves types correctly.
  data.frame(
    run_id         = result$run_id,
    T_g            = result$T_g,
    t_inf          = result$t_inf,
    N_ref          = result$N_ref,
    N_compare      = result$N_compare,
    pct_decline    = result$pct_decline,
    iucn_vu        = as.integer(result$pct_decline >= 30),  # Vulnerable
    iucn_en        = as.integer(result$pct_decline >= 50),  # Endangered
    iucn_cr        = as.integer(result$pct_decline >= 80),  # Critically Endangered
    n_offspring_gt = result$n_offspring_gt,
    year_compare   = result$compare_mid,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# collect_iucn_results()
# ============================================================================
# Merge all IUCN .rds files from one or more results directories.
# Each .rds must be a result list from run_iucn_a3().
# Returns a data.frame with one row per run, deduplicated by run_id.
collect_iucn_results <- function(results_dirs) {
  if (is.character(results_dirs)) results_dirs <- as.list(results_dirs)

  files <- unlist(lapply(results_dirs, function(d)
    list.files(d, pattern = "\\.rds$", full.names = TRUE)))

  if (length(files) == 0L) {
    warning("No .rds files found.")
    return(data.frame())
  }

  rows <- lapply(files, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || !is.list(x) || is.null(x$pct_decline)) return(NULL)
    summarise_iucn_a3(x)
  })

  rows <- rows[!sapply(rows, is.null)]
  out  <- do.call(rbind, rows)
  out  <- out[!duplicated(out$run_id), ]
  rownames(out) <- NULL
  message(sprintf("Collected %d unique IUCN runs.", nrow(out)))
  out
}

# ============================================================================
# plot_iucn_a3()
# ============================================================================
plot_iucn_a3 <- function(result, save_path = NULL,
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

# --- ggplot2 rendering ------------------------------------------------------
.gg_iucn <- function(result) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes

  cen       <- result$census
  yr_rust   <- result$year_rust
  yr_spin   <- result$year_spinup
  T_g       <- result$T_g
  t_inf     <- result$t_inf
  N_ref     <- result$N_ref
  N_compare <- result$N_compare
  pct       <- result$pct_decline
  cat       <- result$iucn_cat

  # --- Panel 1: N_IUCN time series -----------------------------------------
  cen$N_IUCN_smooth <- cen$N_IUCN  # could smooth if desired

  # Background annotation bands
  ymax <- max(cen$N_IUCN, na.rm = TRUE) * 1.15

  p1 <- gg(cen, aes(x = year, y = N_IUCN)) +
    # Spin-up shade
    ggplot2::annotate("rect",
      xmin = yr_spin, xmax = result$year_gentime_start,
      ymin = 0, ymax = ymax,
      fill = "#dddddd", alpha = 0.4) +
    ggplot2::annotate("text",
      x = (yr_spin + result$year_gentime_start) / 2,
      y = ymax * 0.97, label = "spin-up",
      size = 2.8, colour = "#999999", vjust = 1) +
    # Gen-time window
    ggplot2::annotate("rect",
      xmin = result$year_gentime_start, xmax = result$year_gentime_end,
      ymin = 0, ymax = ymax,
      fill = "#c6dbef", alpha = 0.35) +
    ggplot2::annotate("text",
      x = (result$year_gentime_start + result$year_gentime_end) / 2,
      y = ymax * 0.97,
      label = paste0("gen-time\nwindow\nT\u1d33=", round(T_g, 1), " yr"),
      size = 2.8, colour = "#2a78d6", vjust = 1) +
    # Reference window
    ggplot2::annotate("rect",
      xmin = result$ref_start, xmax = result$ref_end,
      ymin = 0, ymax = ymax,
      fill = "#a1d99b", alpha = 0.5) +
    # Comparison window
    ggplot2::annotate("rect",
      xmin = result$compare_start, xmax = result$compare_end,
      ymin = 0, ymax = ymax,
      fill = "#fdae6b", alpha = 0.5) +
    # Main trajectory
    ggplot2::geom_line(colour = "#333333", linewidth = 0.7) +
    # Fire events: rug marks along the bottom of the panel
    { fires <- cen[cen$fire_event, ]
      if (nrow(fires) > 0)
        ggplot2::geom_rug(data = fires, aes(x = year),
                           inherit.aes = FALSE,
                           colour = "#e08000", linewidth = 0.5,
                           sides = "b", length = grid::unit(0.04, "npc"))
    } +
    # Counterfactual if available
    { if (!is.null(result$census_cfact))
        ggplot2::geom_line(data = result$census_cfact,
                            aes(x = year, y = N_IUCN),
                            colour = "#2a78d6", linewidth = 0.6,
                            linetype = "dashed")
    } +
    # Rust arrival
    ggplot2::geom_vline(xintercept = yr_rust,
                         colour = "#e34948", linewidth = 0.9, linetype = "dashed") +
    ggplot2::annotate("text", x = yr_rust + 1, y = ymax * 0.88,
                       label = paste0("Rust\n", yr_rust),
                       colour = "#e34948", hjust = 0, size = 3) +
    # Reference and comparison levels
    ggplot2::geom_hline(yintercept = N_ref,
                         colour = "#1a7a35", linewidth = 0.8, linetype = "longdash") +
    ggplot2::geom_hline(yintercept = N_compare,
                         colour = "#b85c00", linewidth = 0.8, linetype = "longdash") +
    # Decline annotation
    ggplot2::annotate("segment",
      x = result$compare_mid, xend = result$compare_mid,
      y = N_compare, yend = N_ref,
      colour = "#888780", linewidth = 0.6,
      arrow = grid::arrow(ends = "both", length = grid::unit(0.08, "inches"))) +
    ggplot2::annotate("text",
      x = result$compare_mid + (result$year_end - yr_rust) * 0.02,
      y = (N_ref + N_compare) / 2,
      label = sprintf("%.0f%% decline\n[%s]", pct, cat),
      hjust = 0, size = 3.5, colour = "#333333") +
    # Labels
    ggplot2::annotate("text",
      x = result$ref_end + 1, y = N_ref * 1.03,
      label = sprintf("N_ref = %.0f", N_ref),
      hjust = 0, size = 2.8, colour = "#1a7a35") +
    ggplot2::annotate("text",
      x = result$compare_end + 1, y = N_compare * 0.97,
      label = sprintf("N_compare = %.0f", N_compare),
      hjust = 0, size = 2.8, colour = "#b85c00", vjust = 1) +
    ggplot2::scale_y_continuous(limits = c(0, ymax),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = "Year",
      y = "Reproductive adults (N_IUCN)",
      title = sprintf(
        "IUCN A3 assessment  |  T\u1d33 = %.1f yr  |  assessment period = %.0f yr  |  %s",
        T_g, t_inf, sprintf("%.0f%% projected decline  \u2192  %s", pct, cat)),
      subtitle = sprintf(
        "Reference: %d\u2013%d (mean = %.0f)   |   Comparison: %d\u2013%d (mean = %.0f)   |   \u25a1 fire events (rug)",
        result$ref_start, result$ref_end, N_ref,
        result$compare_start, result$compare_end, N_compare)
    ) +
    ggplot2::theme_minimal(base_size = 11)

  # --- Panel 2: maternal age histogram (gen time) --------------------------
  ma_df <- data.frame(age = result$mother_ages)

  p2 <- gg(ma_df, aes(x = age)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "#2a78d6", alpha = 0.7,
                             colour = "white") +
    ggplot2::geom_vline(xintercept = T_g,
                         colour = "#e34948", linewidth = 1, linetype = "dashed") +
    ggplot2::annotate("text", x = T_g + 0.5, y = Inf,
                       label = sprintf("T\u1d33 = %.1f yr", T_g),
                       hjust = 0, vjust = 1.4, size = 3.5, colour = "#e34948") +
    ggplot2::labs(
      x = "Maternal age at offspring birth (years)",
      y = "Count",
      title = "Generation time",
      subtitle = sprintf("n = %d offspring  |  window %d\u2013%d",
                          result$n_offspring_gt,
                          result$year_gentime_start, result$year_gentime_end)
    ) +
    ggplot2::theme_minimal(base_size = 11)

  # --- Assemble -----------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(p1, p2, nrow = 2, heights = c(2.5, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, nrow = 2, heights = c(2.5, 1))
  } else {
    print(p1); readline("Enter for gen time panel..."); print(p2)
  }
}

# --- base R fallback --------------------------------------------------------
.base_iucn <- function(result) {
  op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  cen <- result$census
  yr_rust <- result$year_rust

  # Panel 1
  ymax <- max(cen$N_IUCN, na.rm = TRUE) * 1.15
  plot(cen$year, cen$N_IUCN, type = "l", lwd = 1.5, col = "#333333",
       ylim = c(0, ymax), xlab = "Year",
       ylab = "Reproductive adults (N_IUCN)",
       main = sprintf("IUCN A3: %.0f%% projected decline [%s]  T_g=%.1f yr  period=%.0f yr",
                       result$pct_decline, result$iucn_cat,
                       result$T_g, result$t_inf), las = 1)

  rect(result$year_spinup, 0, result$year_gentime_start, ymax,
       col = grDevices::adjustcolor("#dddddd", 0.4), border = NA)
  rect(result$year_gentime_start, 0, result$year_gentime_end, ymax,
       col = grDevices::adjustcolor("#c6dbef", 0.4), border = NA)
  rect(result$ref_start, 0, result$ref_end, ymax,
       col = grDevices::adjustcolor("#a1d99b", 0.5), border = NA)
  rect(result$compare_start, 0, result$compare_end, ymax,
       col = grDevices::adjustcolor("#fdae6b", 0.5), border = NA)
  lines(cen$year, cen$N_IUCN, col = "#333333", lwd = 1.5)
  fires <- cen[cen$fire_event, ]
  if (nrow(fires) > 0)
    rug(fires$year, col = "#e08000", ticksize = 0.04, lwd = 1)
  abline(v = yr_rust, col = "#e34948", lty = 2, lwd = 1.5)
  abline(h = result$N_ref,     col = "#1a7a35", lty = 5, lwd = 1.2)
  abline(h = result$N_compare, col = "#b85c00", lty = 5, lwd = 1.2)

  # Panel 2
  hist(result$mother_ages, breaks = 30, col = "#2a78d6", border = "white",
       xlab = "Maternal age at offspring birth (years)",
       main = sprintf("Generation time  (T_g = %.1f yr,  n = %d offspring)",
                       result$T_g, result$n_offspring_gt))
  abline(v = result$T_g, col = "#e34948", lty = 2, lwd = 2)
}
