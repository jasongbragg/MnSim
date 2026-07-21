# runner/morris.R
#
# Morris elementary-effects sensitivity screening.
#
# WHAT THIS DOES
# --------------
# 1. Builds a Morris OAT design matrix from param_ranges.
# 2. Maps each design row → nested params list via vec_to_params().
# 3. Runs run_simulation() for each (in parallel), computes summary stats.
# 4. Calls sensitivity::tell() to compute elementary effects.
# 5. Writes TWO CSVs to results_dir:
#      morris_ee.csv      Elementary effects per (output, factor) -- the
#                         primary deliverable. Columns: output, factor,
#                         mu, mu_star, sigma. mu_star is the key screening
#                         statistic (mean |effect|); sigma diagnoses
#                         nonlinearity / interactions.
#      morris_runs.csv    All raw runs: one row per simulation, columns
#                         are the __ param values plus all summary stats.
#                         Useful for scatter plots and further diagnostics.
#
# DEPENDENCIES
# ------------
# Requires package 'sensitivity'. Run: install.packages("sensitivity")
# Uses parallel::mclapply (n_cores > 1 only works on Linux/macOS).
#
# USAGE
# -----
#   source("runner/params_utils.R")
#   source("runner/summarise.R")
#   source("runner/morris.R")
#   # load all model files...
#   source("parlib/ensemble/no_rust_no_fire/params_base.R")
#   base   <- get_default_params()
#   source("parlib/ensemble/no_rust_no_fire/param_ranges.R")
#   ranges <- get_param_ranges()
#
#   out <- run_morris(base, ranges, r = 20, n_cores = 8,
#                     results_dir = "results/morris_no_rust_no_fire/")

run_morris <- function(
    base_params,
    param_ranges,
    r               = 20L,
    levels          = 5L,
    n_summary_years = 50L,
    n_cores         = max(1L, parallel::detectCores() - 1L),
    results_dir     = "results/morris/",
    skip_existing   = TRUE
) {
  if (!requireNamespace("sensitivity", quietly = TRUE))
    stop("Package 'sensitivity' is required.  Run: install.packages('sensitivity')")

  keys <- param_ranges$keys
  mins <- param_ranges$mins
  maxs <- param_ranges$maxs
  k    <- length(keys)

  # --- Morris OAT design --------------------------------------------------
  m_design <- sensitivity::morris(
    model  = NULL,
    factors = keys,
    r      = r,
    design = list(type = "oat", levels = levels,
                   grid.jump = max(1L, ceiling(levels / 2L))),
    binf   = mins,
    bsup   = maxs
  )

  param_mat <- m_design$X          # n_runs × k matrix
  n_runs    <- nrow(param_mat)
  message(sprintf("Morris: %d runs, %d parameters (r=%d, levels=%d)",
                   n_runs, k, r, levels))

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Run one parameter set ----------------------------------------------
  # --- Run one parameter set ----------------------------------------------
  # stat_names_fallback: used to name the NA vector on failed runs.
  # Derived by running summarise_run on a dummy result from the base params
  # so it stays in sync with summarise_run() automatically.
  .dummy_names <- tryCatch({
    dummy <- run_simulation(
      within(base_params, { n_years <- 5L; N0 <- 50L }),
      verbose = FALSE
    )
    names(summarise_run(dummy, n_summary_years = 1L))
  }, error = function(e) {
    c("N_alive_mean","N_alive_cv","N_alive_trend",
      "N_IUCN_mean","adult_frac_mean","adult_frac_final",
      "juv_age_mean","adult_age_mean","recruit_rate",
      "births_mean","extinction")
  })
  n_stats_expected <- length(.dummy_names)

  run_one <- function(i) {
    vec <- setNames(as.numeric(param_mat[i, ]), keys)
    p   <- vec_to_params(base_params, vec)
    id  <- make_run_id(p)
    out <- file.path(results_dir, paste0(id, ".rds"))

    if (skip_existing && file.exists(out)) {
      return(readRDS(out)$stats)
    }

    res <- tryCatch(
      run_simulation(p, verbose = FALSE),
      error = function(e) {
        message(sprintf("  Run %d (id=%s) error: %s", i, id, conditionMessage(e)))
        NULL
      }
    )

    stats <- if (is.null(res)) {
      s <- rep(NA_real_, n_stats_expected)
      names(s) <- .dummy_names
      s
    } else {
      summarise_run(res, n_summary_years = n_summary_years)
    }

    saveRDS(list(run_id = id, params_vec = vec, stats = stats), out)
    stats
  }

  # --- Parallel execution -------------------------------------------------
  message(sprintf("Running on %d core(s)...", n_cores))
  stats_list <- parallel::mclapply(seq_len(n_runs), run_one, mc.cores = n_cores)
  stats_mat  <- do.call(rbind, stats_list)
  stat_names <- colnames(stats_mat)

  # --- Tell Morris --------------------------------------------------------
  sensitivity::tell(m_design, stats_mat)

  # --- Write morris_runs.csv ----------------------------------------------
  runs_df <- cbind(
    run_index = seq_len(n_runs),
    as.data.frame(param_mat,    stringsAsFactors = FALSE),
    as.data.frame(stats_mat,    stringsAsFactors = FALSE)
  )
  runs_csv <- file.path(results_dir, "morris_runs.csv")
  write.csv(runs_df, runs_csv, row.names = FALSE)
  message("Run-level results: ", runs_csv)

  # --- Write morris_ee.csv ------------------------------------------------
  # sensitivity::tell() stores m$ee as a 3D array [r, k, n_stats] when Y
  # is a matrix (multiple outputs), or a 2D matrix [r, k] for a single
  # output. Normalise to a list of matrices before extracting effects.
  ee_raw <- m_design$ee
  ee_list <- if (length(dim(ee_raw)) == 3L) {
    lapply(seq_len(dim(ee_raw)[3L]), function(j) ee_raw[, , j])
  } else {
    list(ee_raw)   # single output: wrap so the loop below is uniform
  }

  ee_rows <- do.call(rbind, lapply(seq_along(stat_names), function(j) {
    ee_j <- ee_list[[j]]    # r × k matrix of elementary effects
    data.frame(
      output  = stat_names[j],
      factor  = keys,
      mu      = colMeans(ee_j,      na.rm = TRUE),
      mu_star = colMeans(abs(ee_j), na.rm = TRUE),
      sigma   = apply(ee_j, 2L, sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  ee_csv <- file.path(results_dir, "morris_ee.csv")
  write.csv(ee_rows, ee_csv, row.names = FALSE)
  message("Elementary effects:  ", ee_csv)
  message("Done.")

  invisible(list(morris = m_design, ee = ee_rows, runs = runs_df,
                  ee_csv = ee_csv, runs_csv = runs_csv))
}

# --- Convenience plot -------------------------------------------------------
# mu_star vs sigma for one output statistic -- the standard Morris screening
# plot. Points in the top-right have large effects AND interactions/
# nonlinearity; points near the x-axis have large effects but are near-
# linear in that statistic.
plot_morris_ee <- function(ee_df, stat = "adult_frac_mean",
                            label_top = 5L, ...) {
  sub <- ee_df[ee_df$output == stat, ]
  if (nrow(sub) == 0) stop("stat '", stat, "' not found in ee_df")

  plot(sub$mu_star, sub$sigma,
       xlab = expression(mu*"*  (mean |elementary effect|)"),
       ylab = expression(sigma~"(SD of elementary effect)"),
       main = paste0("Morris screening: ", stat),
       pch = 16, col = "#2a78d6", ...)

  # Label the top-n most influential factors
  ord <- order(sub$mu_star, decreasing = TRUE)
  top <- head(ord, label_top)
  text(sub$mu_star[top], sub$sigma[top], labels = sub$factor[top],
       pos = 3L, cex = 0.8)

  abline(a = 0, b = 1, lty = 2, col = "#888780")   # mu* = sigma reference line
  invisible(sub)
}
