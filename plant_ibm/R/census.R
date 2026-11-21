# R/census.R
#
# One row per simulated year, with everything an IUCN Criterion A
# assessment and an evolutionary-rescue analysis both need:
#
# N_IUCN       mature individuals in the IUCN sense (Guidelines v16 s4.3):
#              alive plants that have reached reproductive age
#              (age >= age_first_flower), INCLUDING those whose flowering is
#              suppressed this year by fire (resprouting) or rust. A plant
#              temporarily kept from flowering is still capable of
#              reproduction; it is counted.
# N_flowering  the stricter count: alive plants actually flowering this year.
#              Kept for comparison (and as an alternative assessment metric:
#              assess_a4(n_metric = "N_flowering")).
# Neither is total abundance: juveniles are excluded from both.

#' Build the census row for year t.
census_year <- function(pop, resist_gt, t, n_rec, n_dead, params,
                         fire_this_year = FALSE, fire_in_cell = fire_this_year,
                         n_dead_fire = 0L) {
  alive   <- pop$alive
  N_alive <- sum(alive)

  freqs <- allele_freqs(resist_gt, alive)
  names(freqs) <- paste0("freq_locus", seq_along(freqs))

  row <- data.frame(
    year              = t,
    N_alive           = N_alive,
    N_IUCN            = sum(alive & pop$age >= pop$age_first_flower),
    N_flowering       = sum(alive & pop$flowering),
    N_juvenile        = sum(alive & pop$age < pop$age_first_flower & !pop$resprout),
    N_resprout        = sum(alive & pop$resprout),
    births            = n_rec,
    deaths            = n_dead,        # all deaths this year (fire + other causes)
    deaths_fire       = n_dead_fire,   # of which killed by fire
    mean_resist_score = if (N_alive > 0) mean(pop$resist_score[alive]) else NA_real_,
    fire_event        = fire_this_year,   # fire reached the simulated stand
    fire_in_cell      = fire_in_cell,     # fire occurred in the cell

    rust_active       = t >= params$rust_start_year,
    extinct           = N_alive == 0
  )

  row <- cbind(row, as.data.frame(t(freqs)))
  if (is_damage_mode(params)) row <- cbind(row, damage_census(pop))
  row
}

# Damage mode (notes T19): mean rust damage of living mature plants, and the
# share of them with damage >= DAMAGE_REPORT_THRESHOLD (a reporting setting,
# not biology), for matching crown transparency and dieback.
DAMAGE_REPORT_THRESHOLD <- 0.6
damage_census <- function(pop) {
  m <- pop$alive & pop$age >= pop$age_first_flower
  D <- pop$rust_damage[m]
  data.frame(mean_damage   = if (length(D)) mean(D) else NA_real_,
             share_damaged = if (length(D)) mean(D >= DAMAGE_REPORT_THRESHOLD) else NA_real_)
}

# Rust's reproductive hooks on (notes T21): the share of eligible plants rust
# kept from flowering this year (counted at flowering, before mortality), and
# the mean fruit success of flowering plants (NA when no plant flowered or the
# fruit hook is off). Reported only; the mature count is unchanged (notes A4).
repro_census <- function(flower_stats, fruit_success) {
  data.frame(share_flowering_suppressed =
               if (flower_stats[["eligible"]] > 0) flower_stats[["suppressed"]] / flower_stats[["eligible"]] else NA_real_,
             fruit_success = if (is.null(fruit_success)) NA_real_ else fruit_success)
}

#' Simple IUCN Criterion A style percent-decline calculator on the
#' mature-individual count (N_IUCN), between any two years present in a census
#' data frame (e.g. a fixed window, or 3-generation length).
pct_decline <- function(census, y0, y1) {
  n0 <- census$N_IUCN[census$year == y0]
  n1 <- census$N_IUCN[census$year == y1]
  if (length(n0) == 0 || length(n1) == 0 || n0 == 0) return(NA_real_)
  100 * (n0 - n1) / n0
}
