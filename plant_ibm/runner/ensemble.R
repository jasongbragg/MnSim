# runner/ensemble.R
#
# Generate parameter matrices and run large ensembles of simulations.
# Used for ABC, posterior predictive checks, and scenario sweeps.
#
# WORKFLOW
# --------
# 1. generate_params_matrix()   -- draw n parameter sets from the prior
# 2. ensemble_run()             -- run all sets, write one .rds per run
# 3. collect_results()          -- merge .rds files → single data.frame
# 4. write.csv(collect_results(...), "results.csv")
#
# SPLITTING ACROSS SERVERS
# ------------------------
# To split a 5000-run ensemble across three servers:
#   mat <- generate_params_matrix(ranges, n = 5000, seed = 42)
#   write.csv(mat, "ensemble_params.csv", row.names = FALSE)
#   # On each server, load the same CSV and run a slice:
#   mat  <- read.csv("ensemble_params.csv")
#   rows <- 1:1667          # or 1668:3333, or 3334:5000
#   ensemble_run(mat[rows, ], base, results_dir = "results/server1/")
#   # collect_results() merges all three results/ subdirectories

# --- generate_params_matrix --------------------------------------------------
# Draw n parameter sets from the prior defined in param_ranges.
# Returns a data.frame with one column per __ key.
#
# method = "lhs"       Latin hypercube (best coverage; requires package 'lhs')
# method = "random"    Independent uniform draws (simplest)
# method = "mvn_copula" Multivariate normal copula respecting param_ranges$Sigma
#                      (requires package 'MASS'; use when parameters are
#                      correlated, e.g. due to shared environmental drivers)
generate_params_matrix <- function(param_ranges, n, method = "lhs", seed = 42L) {
  set.seed(seed)
  keys <- param_ranges$keys
  mins <- param_ranges$mins
  maxs <- param_ranges$maxs
  k    <- length(keys)

  raw <- switch(method,

    lhs = {
      if (!requireNamespace("lhs", quietly = TRUE))
        stop("Package 'lhs' required for method='lhs'.  Run: install.packages('lhs')")
      lhs::randomLHS(n, k)                        # n × k, uniform [0,1]
    },

    random = {
      matrix(runif(n * k), nrow = n, ncol = k)    # n × k, uniform [0,1]
    },

    mvn_copula = {
      if (is.null(param_ranges$Sigma))
        stop("param_ranges$Sigma must be provided for method='mvn_copula'")
      if (!requireNamespace("MASS", quietly = TRUE))
        stop("Package 'MASS' required for method='mvn_copula'.")
      z <- MASS::mvrnorm(n, mu = rep(0, k), Sigma = param_ranges$Sigma)
      pnorm(z)                                     # Gaussian copula → uniform marginals
    },

    stop("Unknown method '", method, "'. Choose: 'lhs', 'random', 'mvn_copula'")
  )

  # Scale each column from [0,1] to [min, max]
  mat <- matrix(0, nrow = n, ncol = k, dimnames = list(NULL, keys))
  for (j in seq_len(k)) mat[, j] <- mins[j] + raw[, j] * (maxs[j] - mins[j])

  as.data.frame(mat, stringsAsFactors = FALSE)
}

# --- ensemble_run ------------------------------------------------------------
# Run all rows of params_matrix, writing one .rds per run to results_dir.
# skip_existing = TRUE: skip any run whose .rds already exists (restartable).
# Each .rds contains: list(run_id, params_vec, stats).
ensemble_run <- function(
    params_matrix,
    base_params,
    n_summary_years = 50L,
    n_cores         = max(1L, parallel::detectCores() - 1L),
    results_dir     = "results/ensemble/",
    skip_existing   = TRUE,
    verbose_every   = 100L   # message every N completed runs
) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  keys    <- names(params_matrix)
  n_runs  <- nrow(params_matrix)
  message(sprintf("Ensemble: %d runs on %d core(s)", n_runs, n_cores))

  run_one <- function(i) {
    vec <- setNames(as.numeric(params_matrix[i, ]), keys)
    p   <- vec_to_params(base_params, vec)
    id  <- make_run_id(p)
    out <- file.path(results_dir, paste0(id, ".rds"))

    if (skip_existing && file.exists(out)) return(invisible(NULL))

    res <- tryCatch(
      run_simulation(p, verbose = FALSE),
      error = function(e) {
        message(sprintf("  Run %d (id=%s) error: %s", i, id, conditionMessage(e)))
        NULL
      }
    )

    stats <- if (is.null(res)) {
      s <- rep(NA_real_, 11L)
      names(s) <- c("N_alive_mean","N_alive_cv","N_alive_trend",
                     "N_IUCN_mean","adult_frac_mean","adult_frac_final",
                     "juv_age_mean","adult_age_mean","recruit_rate",
                     "births_mean","extinction")
      s
    } else {
      summarise_run(res, n_summary_years = n_summary_years)
    }

    saveRDS(list(run_id = id, params_vec = vec, stats = stats), out)

    if (i %% verbose_every == 0L)
      message(sprintf("  Completed %d / %d", i, n_runs))

    invisible(NULL)
  }

  parallel::mclapply(seq_len(n_runs), run_one, mc.cores = n_cores)
  message("Ensemble runs complete.")
  invisible(results_dir)
}

# --- collect_results ---------------------------------------------------------
# Merge all .rds files in one or more results directories into a single
# data.frame. Handles multiple directories (e.g. after splitting across
# servers). Duplicate run_ids (same run on two servers) are deduplicated.
collect_results <- function(results_dirs) {
  if (is.character(results_dirs)) results_dirs <- as.list(results_dirs)

  files <- unlist(lapply(results_dirs, function(d)
    list.files(d, pattern = "\\.rds$", full.names = TRUE)))

  if (length(files) == 0L) {
    warning("No .rds files found in supplied directories.")
    return(data.frame())
  }

  rows <- lapply(files, function(f) {
    x <- readRDS(f)
    as.data.frame(
      c(list(run_id = x$run_id), as.list(x$params_vec), as.list(x$stats)),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$run_id), ]
  rownames(out) <- NULL
  message(sprintf("Collected %d unique runs from %d files.", nrow(out), length(files)))
  out
}
