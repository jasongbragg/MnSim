# spatial/fire/R/estimate.R -- per-cell fire regime from cell-year burnt fractions.
#
# Same quantities as the pre-Nov-2026 spatial/02 (docs/methodological_notes.md
# D2), with each cell allowed its own observation length T (NSW and Qld
# records cover different years):
#   k   = number of years with frac_burned >= AREA_THRESH ("fire years")
#   fire_prob_annual = posterior mean (alpha + k) / (beta + T) under a
#         Gamma(alpha, beta) prior on the cell's fire-year rate, clamped to
#         [FIRE_MIN, FIRE_MAX]
#   mean_frac_burned / sd_frac_burned over fire years
#   fire_p_fimp = mean_frac_burned, or the median over cells with fires.
#
# How the prior is fitted (argument 'prior'):
#   "nb_ml"  (default) empirical Bayes by maximum likelihood: k ~ Poisson(lambda T),
#            lambda ~ Gamma(alpha, beta), so k is negative binomial; alpha and beta
#            maximise that likelihood over ALL cells, each with its own T. Uses
#            the cells with no fire years too, and separates Poisson noise from
#            real between-cell variation (important for Qld's short record).
#   "mm02"   the pre-Nov-2026 02's method of moments on k/T over cells with
#            k > 0 (kept for comparison; it reproduced the old table exactly).
#   "fixed"  PRIOR_ALPHA, PRIOR_BETA.
# 'fire_min': floor on the rate (FIRE_MIN for all fire or wildfire; 0 for
# planned burns, where many cells are simply never burnt by an agency).
# 'prior_group': optional vector (one per cell, e.g. state); the prior is then
# fitted separately in each group. NULL = one prior for all cells.
# Output columns match those 02 writes to mn_site_params.csv, plus T_obs and
# prior_group.

# Negative-binomial (Gamma-Poisson) log likelihood, exposure T per cell.
nb_loglik <- function(alpha, beta, k, T)
  sum(lgamma(alpha + k) - lgamma(alpha) - lgamma(k + 1) +
      alpha * log(beta / (beta + T)) + k * log(T / (beta + T)))

fit_prior_nb <- function(k, T) {
  if (sum(k) == 0 || length(k) < 5) return(NULL)
  m0 <- sum(k) / sum(T)
  nll <- function(p) -nb_loglik(exp(p[1]), exp(p[2]), k, T)
  start <- c(log(2), log(2 / m0))
  o <- tryCatch(stats::optim(start, nll, method = "L-BFGS-B",
                             lower = c(log(0.05), log(0.5)), upper = c(log(1e4), log(1e6))),
                error = function(e) NULL)
  if (is.null(o) || o$convergence != 0) return(NULL)
  c(alpha = exp(o$par[1]), beta = exp(o$par[2]), loglik = -o$value)
}

fit_prior_mm02 <- function(k, T) {
  lam <- (k / T)[k > 0]
  if (length(lam) < 10) return(NULL)
  m <- mean(lam); v <- var(lam)
  if (!(v > 0 && m^2 / v > 0.5)) return(NULL)
  c(alpha = m^2 / v, beta = m / v, loglik = nb_loglik(m^2 / v, m / v, k, T))
}

fit_prior <- function(k, T, prior = c("nb_ml", "mm02", "fixed")) {
  prior <- match.arg(prior)
  p <- switch(prior, nb_ml = fit_prior_nb(k, T), mm02 = fit_prior_mm02(k, T), fixed = NULL)
  if (is.null(p))
    return(list(alpha = PRIOR_ALPHA, beta = PRIOR_BETA, type = "fixed",
                loglik = nb_loglik(PRIOR_ALPHA, PRIOR_BETA, k, T)))
  list(alpha = p[["alpha"]], beta = p[["beta"]], type = prior, loglik = p[["loglik"]])
}

estimate_fire_regime <- function(cy, cells, T_obs, prior = c("nb_ml", "mm02", "fixed"),
                                 prior_group = NULL, fire_min = FIRE_MIN) {
  prior <- match.arg(prior)
  # cy: data.frame(cell_id, year, frac_burned); T_obs: named by cell_id (years observed)
  ids <- cells$cell_id
  per <- lapply(ids, function(id) {
    f   <- cy$frac_burned[cy$cell_id == id]
    sig <- f[f >= AREA_THRESH]
    data.frame(cell_id = id,
               n_wildfire_years_all = length(f),
               n_wildfire_years_sig = length(sig),
               mean_frac_burned = if (length(sig)) mean(sig) else NA_real_,
               sd_frac_burned   = if (length(sig) > 1) sd(sig) else NA_real_,
               var_frac_burned  = if (length(sig) > 1) var(sig) else NA_real_)
  })
  r <- do.call(rbind, per)
  r$T_obs <- as.numeric(T_obs[as.character(r$cell_id)])
  r$prior_group <- if (is.null(prior_group)) "all" else as.character(prior_group)

  priors <- do.call(rbind, lapply(split(seq_len(nrow(r)), r$prior_group), function(i) {
    p <- fit_prior(r$n_wildfire_years_sig[i], r$T_obs[i], prior)
    data.frame(group = r$prior_group[i[1]], alpha = p$alpha, beta = p$beta, type = p$type,
               loglik = p$loglik, cells = length(i), prior_mean = p$alpha / p$beta)
  }))
  g <- match(r$prior_group, priors$group)
  r$fire_prob_annual <- pmin(FIRE_MAX, pmax(fire_min,
    (priors$alpha[g] + r$n_wildfire_years_sig) / (priors$beta[g] + r$T_obs)))
  pf_default <- stats::median(r$mean_frac_burned, na.rm = TRUE)
  r$fire_p_fimp <- pmin(1, pmax(PFIMP_MIN, dplyr::coalesce(r$mean_frac_burned, pf_default)))
  r$FRI_years <- round(1 / r$fire_prob_annual, 1)
  rownames(priors) <- NULL
  attr(r, "priors") <- priors
  r
}

# Does one prior fit both states? Likelihood-ratio test of a shared NB prior
# against one per group (2 extra parameters).
prior_lr_test <- function(k, T, group) {
  p0 <- fit_prior_nb(k, T)
  p1 <- lapply(split(seq_along(k), group), function(i) fit_prior_nb(k[i], T[i]))
  if (is.null(p0) || any(vapply(p1, is.null, TRUE))) return(NULL)
  ll1 <- sum(vapply(p1, function(p) p[["loglik"]], 0))
  stat <- 2 * (ll1 - p0[["loglik"]])
  df <- 2 * (length(p1) - 1)
  c(lr = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE))
}
