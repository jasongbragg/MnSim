# tests/regression/stats.R
#
# Second regression check (PLAN.md Phase 1): summary statistics over 20 seeds,
# compared within a tolerance. The exact-match tests (test-regression.R) fail
# on ANY change to the random-number stream, including harmless refactors that
# only reorder draws. These statistics tell such a refactor (an independent
# sample from the same model: means agree within sampling error) from a real
# change in behaviour (means move).
#
# Each statistic re-uses a scenario from scenarios.R with 20 different seeds.
# make_reference_stats.R freezes the per-seed values in
# ../tests/reference/summary_stats.rds, with the params each seed ran with.
#
# Tolerance: with the same code the means are identical (same seeds). After a
# draw-order change the new mean is an independent sample, so
# mean_new - mean_ref has sd sqrt(2) * sd_ref / sqrt(n). The test allows
# REGRESSION_STATS_Z (4) of those sds: a false alarm about 1 time in 16,000.
#
# What this check can see (tested 2026-11-19, 20 seeds of ~400 plants): a
# 20% rise in rust mortality fails rust_mortality, and the Pegg fraction moves
# by about its tolerance. Changes in adult mortality below ~5% of its rate (e.g.
# background lambda 60 -> 55) pass. It is a coarse guard for refactors that
# reorder random draws; the exact-match tests catch everything else.

REGRESSION_STATS_N <- 20L
REGRESSION_STATS_Z <- 4

REGRESSION_STATS <- list(

  # Mature plants only: the all-plant death rate is mostly juveniles, and
  # recruitment compensates for adult deaths, so it hardly moves when adult
  # mortality changes.
  mature_mortality = list(
    description = "no_fire scenario: mean annual death rate of mature plants (deaths / N_IUCN the year before), 1950-1999",
    scenario = "no_fire",
    seed0    = 2000L,
    stat = function(out) {
      ind <- out$individuals; c <- out$census; yrs <- 1950:1999
      d <- vapply(yrs, function(y) sum(ind$death_year %in% y & ind$age >= ind$age_first_flower), 0)
      mean(d / c$N_IUCN[match(yrs - 1L, c$year)])
    }
  ),

  pegg_resprout_dead = list(
    description = "pegg scenario: fraction of resprouting fire survivors dead at the Jan 2022 survey (fire + rust mortality)",
    scenario = "pegg",
    seed0    = 2100L,
    stat = function(out) out$fraction_resprout_dead
  ),

  # Fire off: fire timing dominates the between-seed noise of rust outcomes.
  # (A decline statistic was tried first and dropped: recruitment compensates,
  # so a 20% rise in rust mortality left the 2009-2060 decline unchanged.)
  rust_mortality = list(
    description = "rust_high scenario with fire off: mean annual death rate 2011-2060 under anchor rust pressure",
    scenario = "rust_high",
    modify   = function(p) { p$fire_prob_annual <- 0; p },
    seed0    = 2200L,
    stat = function(out) {
      c <- out$census
      i <- which(c$year >= 2011 & c$year <= 2060)
      mean(c$deaths[i] / c$N_alive[i - 1L])
    }
  )
)

# Params for each seed of a statistic, from a base params list.
regression_stats_params <- function(name, base = regression_base()) {
  s <- REGRESSION_STATS[[name]]
  p <- regression_params(s$scenario, base)
  if (!is.null(s$modify)) p <- s$modify(p)
  lapply(s$seed0 + seq_len(REGRESSION_STATS_N), function(sd) { p$seed <- sd; p })
}

# Per-seed values of a statistic, from a list of params (one per seed).
regression_stats_values <- function(name, params_list) {
  s <- REGRESSION_STATS[[name]]
  vapply(params_list, function(p) s$stat(regression_run(s$scenario, p)), numeric(1))
}

REGRESSION_STATS_FILE <- "summary_stats.rds"
