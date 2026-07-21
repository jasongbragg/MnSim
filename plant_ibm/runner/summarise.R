# runner/summarise.R
#
# summarise_run(): reduce a run_simulation() result to a named numeric
# vector of scalar summary statistics.
#
# VERSIONING NOTE
# ---------------
# This file is intentionally separate from the model code. Changes to
# which statistics you compute, or how, should be recorded here (e.g. in
# git history). Different sets of summary statistics produce different
# ABC posteriors; keeping this file distinct makes it easy to see exactly
# what statistics any given ensemble was computed under.
#
# USAGE
# -----
#   res   <- run_simulation(params)
#   stats <- summarise_run(res, n_summary_years = 50)
#
# n_summary_years: how many years from the END of the simulation to use
# for computing time-averaged statistics. A burn-in of (n_years -
# n_summary_years) years is implicitly discarded. If the run has fewer
# years than n_summary_years, all available years are used with a warning.
#
# Returns a named numeric vector. NA values indicate the statistic could
# not be computed (e.g. no alive juveniles at end of run). An 'extinction'
# flag (0/1) is always included; downstream filtering can discard extinct
# runs before ABC comparison.

summarise_run <- function(res, n_summary_years = 50L) {

  census <- res$census
  n      <- nrow(census)

  if (n < n_summary_years) {
    warning(sprintf(
      "Run has %d census rows but n_summary_years=%d; using all rows.",
      n, n_summary_years))
    sub <- census
  } else {
    sub <- census[(n - n_summary_years + 1L):n, ]
  }

  # --- Final-year population snapshot -------------------------------------
  pop    <- res$individuals
  alive  <- pop[pop$alive, ]
  juv    <- alive[alive$age < alive$age_first_flower & !alive$resprout, ]
  adults <- alive[alive$age >= alive$age_first_flower & !alive$resprout, ]

  # --- Time-averaged statistics (over summary window) ---------------------
  N_alive_vals  <- sub$N_alive
  N_IUCN_vals   <- sub$N_IUCN
  births_vals   <- sub$births

  adult_frac <- N_IUCN_vals / pmax(N_alive_vals, 1L)
  recruit_rt <- births_vals / pmax(N_alive_vals, 1L)

  # Linear trend in N_alive (negative = declining)
  trend <- if (length(N_alive_vals) > 1L) {
    coef(lm(N_alive_vals ~ seq_along(N_alive_vals), na.action = na.omit))[2L]
  } else {
    0
  }

  # --- Assemble -----------------------------------------------------------
  c(
    N_alive_mean     = mean(N_alive_vals, na.rm = TRUE),
    N_alive_cv       = sd(N_alive_vals, na.rm = TRUE) /
                       max(mean(N_alive_vals, na.rm = TRUE), 1),
    N_alive_trend    = as.numeric(trend),
    N_IUCN_mean      = mean(N_IUCN_vals,  na.rm = TRUE),
    adult_frac_mean  = mean(adult_frac,   na.rm = TRUE),
    adult_frac_final = if (nrow(alive) > 0)
                         nrow(adults) / nrow(alive) else NA_real_,
    juv_age_mean     = if (nrow(juv)    > 0) mean(juv$age)    else NA_real_,
    adult_age_mean   = if (nrow(adults) > 0) mean(adults$age) else NA_real_,
    recruit_rate     = mean(recruit_rt, na.rm = TRUE),
    births_mean      = mean(births_vals, na.rm = TRUE),
    extinction       = as.numeric(tail(N_alive_vals, 1L) == 0L)
  )
}
