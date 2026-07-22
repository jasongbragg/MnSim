# R/census.R
#
# One row per simulated year, with everything an IUCN Criterion A
# assessment and an evolutionary-rescue analysis both need:
# N_IUCN is reproductive (flowering) adults ONLY, per the stated
# IUCN counting rule -- never total abundance.

#' Build the census row for year t.
census_year <- function(pop, resist_gt, t, n_rec, n_dead, params,
                         fire_this_year = FALSE) {
  alive   <- pop$alive
  N_alive <- sum(alive)

  freqs <- allele_freqs(resist_gt, alive)
  names(freqs) <- paste0("freq_locus", seq_along(freqs))

  row <- data.frame(
    year              = t,
    N_alive           = N_alive,
    N_IUCN            = sum(alive & pop$flowering),
    N_juvenile        = sum(alive & pop$age < pop$age_first_flower & !pop$resprout),
    N_resprout        = sum(alive & pop$resprout),
    births            = n_rec,
    deaths            = n_dead,
    mean_resist_score = if (N_alive > 0) mean(pop$resist_score[alive]) else NA_real_,
    fire_event        = fire_this_year,
    rust_active       = t >= params$rust_start_year,
    extinct           = N_alive == 0
  )

  cbind(row, as.data.frame(t(freqs)))
}

#' Simple IUCN Criterion A style percent-decline calculator on the
#' reproductive-adult count, between any two years present in a census
#' data frame (e.g. a fixed window, or 3-generation length).
pct_decline <- function(census, y0, y1) {
  n0 <- census$N_IUCN[census$year == y0]
  n1 <- census$N_IUCN[census$year == y1]
  if (length(n0) == 0 || length(n1) == 0 || n0 == 0) return(NA_real_)
  100 * (n0 - n1) / n0
}
