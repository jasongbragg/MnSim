# R/recruitment.R

#' Beverton-Holt recruitment: saturates smoothly, never explodes.
bevholt <- function(N, params) {
  (params$R_max * N) / (params$K_half + N)
}

#' Ricker recruitment (alternative; needs params$ricker_r).
ricker <- function(N, params) {
  N * exp(params$ricker_r * (1 - N / params$K))
}

#' Sample n_rec mother/father ROW INDEX pairs from the flowering pool.
sample_parents <- function(pop, n_rec, params) {
  flower_idx <- which(pop$alive & pop$flowering)
  n_flower   <- length(flower_idx)
  if (n_flower == 0) return(NULL)

  mother_idx <- sample(flower_idx, n_rec, replace = TRUE)

  if (n_flower == 1) {
    return(list(mother_idx = mother_idx, father_idx = mother_idx))
  }

  p_self <- params$selfing_rate / (n_flower - 1 + params$selfing_rate)
  is_self <- runif(n_rec) < p_self

  father_idx <- integer(n_rec)
  father_idx[is_self] <- mother_idx[is_self]

  need <- which(!is_self)
  if (length(need) > 0) {
    candidate <- sample(flower_idx, length(need), replace = TRUE)
    clash <- candidate == mother_idx[need]
    while (any(clash)) {
      candidate[clash] <- sample(flower_idx, sum(clash), replace = TRUE)
      clash <- candidate == mother_idx[need]
    }
    father_idx[need] <- candidate
  }

  list(mother_idx = mother_idx, father_idx = father_idx)
}

#' Standing canopy density, weighted by recovery progress for resprouting
#' individuals. See detailed docstring history in TESTING.md.
#' Key design: uses age >= age_first_flower (not resprout status) to avoid
#' counting surviving juveniles as canopy after a fire.
canopy_density <- function(pop) {
  canopy_age <- pop$alive & pop$age >= pop$age_first_flower
  full_mask  <- canopy_age & !pop$resprout
  resp_mask  <- canopy_age &  pop$resprout

  total    <- pop$resprout_total_yrs[resp_mask]
  remain   <- pop$resprout_yrs_remain[resp_mask]
  progress <- ifelse(total > 0, 1 - remain / total, 1)
  progress <- pmin(pmax(progress, 0), 1)

  sum(full_mask) + sum(progress)
}

#' Fraction of bevholt()'s expected count lost to canopy shading.
shade_suppression <- function(N_canopy, params) {
  pmin(dose_response(N_canopy, params$shade_dose_response), 1)
}

#' Full recruitment step for year t.
#' annual_env_t: pathogen vertex of the disease triangle (passed from
#' simulate.R's annual draw). Default 0 → no rust effect on fecundity
#' (used by tests that call recruit() directly).
#'
#' Rust suppresses per-parent fecundity contribution via
#' rust_recruit_dose_response applied to each adult's effective
#' susceptibility. The "effective" flowering count fed to bevholt() is
#' the sum of each adult's (1 - rust_fecundity_suppression) weight,
#' so a fully resistant adult contributes 1.0 regardless of annual_env_t.
#'
#' Returns list(pop, resist_gt, n_rec).
recruit <- function(pop, resist_gt, t, params, annual_env_t = 0) {
  flower_mask <- pop$alive & pop$flowering
  n_flowering <- sum(flower_mask)
  N_canopy    <- canopy_density(pop)

  # Rust-suppressed fecundity: each flowering adult's contribution to the
  # recruit pool is weighted by (1 - suppression), where suppression is
  # driven by that adult's own effective_susceptibility. Off by default.
  if (params$rust_recruit_dose_response$max_effect > 0 &&
      t >= params$rust_start_year && annual_env_t > 0) {
    flower_idx <- which(flower_mask)
    eff_susc   <- (1 - pop$resist_score[flower_idx]) *
                   pop$env_susceptibility[flower_idx] * annual_env_t
    weights    <- 1 - vapply(eff_susc, dose_response,
                              numeric(1), cfg = params$rust_recruit_dose_response)
    n_flowering_eff <- sum(weights)
  } else {
    n_flowering_eff <- n_flowering
  }

  expected_rec <- bevholt(n_flowering_eff, params) *
                  (1 - shade_suppression(N_canopy, params))
  n_rec <- rpois(1, expected_rec)

  if (n_rec == 0 || n_flowering == 0) {
    return(list(pop = pop, resist_gt = resist_gt, n_rec = 0L))
  }

  parents <- sample_parents(pop, n_rec, params)
  if (is.null(parents)) {
    return(list(pop = pop, resist_gt = resist_gt, n_rec = 0L))
  }

  rec_gt    <- inherit_alleles(resist_gt, parents$mother_idx, parents$father_idx)
  rec_score <- resist_score_from_gt(rec_gt, params)

  new_rows              <- make_recruit_rows(n_rec, t, parents$mother_idx,
                                              parents$father_idx, pop, params)
  new_rows$resist_score <- rec_score

  pop       <- rbind(pop, new_rows)
  resist_gt <- rbind(resist_gt, rec_gt)

  list(pop = pop, resist_gt = resist_gt, n_rec = n_rec)
}
