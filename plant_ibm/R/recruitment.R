# R/recruitment.R

#' Append rows to a data frame column by column. Same result as
#' rbind(pop, new) for these all-atomic columns, without rbind.data.frame's
#' row-name bookkeeping (which dominated run time). Falls back to rbind if
#' the columns differ.
append_rows <- function(pop, new) {
  if (!setequal(names(pop), names(new))) return(rbind(pop, new))
  out <- lapply(names(pop), function(nm) c(pop[[nm]], new[[nm]]))
  names(out) <- names(pop)
  attr(out, "row.names") <- .set_row_names(nrow(pop) + nrow(new))
  class(out) <- "data.frame"
  out
}

#' Beverton-Holt recruitment: saturates smoothly, never explodes.
bevholt <- function(N, params) {
  (params$R_max * N) / (params$K_half + N)
}

#' Ricker recruitment (alternative): N * exp(r * (1 - N / K)), with r =
#' ricker_r and K = ricker_K. Peaks at N = K / r and falls at high N. Its scale
#' is not comparable with bevholt()'s R_max / K_half, so a species using it
#' needs its own calibrated values.
ricker <- function(N, params) {
  N * exp(params$ricker_r * (1 - N / params$ricker_K))
}

#' Recruitment curves: expected recruits from the effective number of
#' flowering plants (design_swappable_interfaces.md section 5; notes T12).
#' params$recruitment_curve picks one; absent means "bevholt". The names are
#' the schema's recruitment_curve choices (tested).
RECRUITMENT_CURVES <- list(bevholt = bevholt, ricker = ricker)

recruitment_curve <- function(params) {
  nm <- if (is.null(params$recruitment_curve)) "bevholt" else params$recruitment_curve
  f  <- RECRUITMENT_CURVES[[nm]]
  if (is.null(f)) stop("Unknown recruitment_curve '", nm, "'. Must be one of: ",
                       paste(names(RECRUITMENT_CURVES), collapse = ", "), ".")
  f
}

#' Sample n_rec mother/father ROW INDEX pairs from the flowering pool.
#' Indexes flower_idx explicitly (sample.int), so a pool of one plant is not
#' misread as sample(1:k); otherwise identical draws to sample(flower_idx, ...).
sample_parents <- function(pop, n_rec, params) {
  flower_idx <- which(pop$alive & pop$flowering)
  n_flower   <- length(flower_idx)
  if (n_flower == 0) return(NULL)

  mother_idx <- flower_idx[sample.int(n_flower, n_rec, replace = TRUE)]

  if (n_flower == 1) {
    return(list(mother_idx = mother_idx, father_idx = mother_idx))
  }

  p_self <- params$selfing_rate / (n_flower - 1 + params$selfing_rate)
  is_self <- runif(n_rec) < p_self

  father_idx <- integer(n_rec)
  father_idx[is_self] <- mother_idx[is_self]

  need <- which(!is_self)
  if (length(need) > 0) {
    candidate <- flower_idx[sample.int(n_flower, length(need), replace = TRUE)]
    clash <- candidate == mother_idx[need]
    while (any(clash)) {
      candidate[clash] <- flower_idx[sample.int(n_flower, sum(clash), replace = TRUE)]
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

#' Expected number of recruits in year t: the recruitment curve applied to the
#' effective number of flowering plants (rust-suppressed fecundity), times
#' (1 - shade suppression). Uses no random numbers. Later fecundity rules
#' (e.g. fruit infection) belong here.
expected_recruits <- function(pop, t, params, annual_env_t = 0) {
  flower_mask <- pop$alive & pop$flowering
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
    n_flowering_eff <- sum(flower_mask)
  }

  recruitment_curve(params)(n_flowering_eff, params) *
    (1 - shade_suppression(N_canopy, params))
}

#' Full recruitment step for year t.
#' annual_env_t: pathogen vertex of the disease triangle (passed from
#' simulate.R's annual draw). Default 0 → no rust effect on fecundity
#' (used by tests that call recruit() directly).
#'
#' The expected count comes from expected_recruits(): rust suppresses
#' per-parent fecundity via rust_recruit_dose_response applied to each
#' adult's effective susceptibility; the "effective" flowering count fed to
#' the recruitment curve is the sum of each adult's (1 - suppression) weight,
#' so a fully resistant adult contributes 1.0 regardless of annual_env_t.
#'
#' next_id: id for the first recruit (run_simulation keeps a running
#' counter because dead rows may have been archived); NULL = max(pop$id) + 1.
#' Returns list(pop, resist_gt, n_rec).
recruit <- function(pop, resist_gt, t, params, annual_env_t = 0, next_id = NULL) {
  n_flowering  <- sum(pop$alive & pop$flowering)
  expected_rec <- expected_recruits(pop, t, params, annual_env_t)
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
                                              parents$father_idx, pop, params,
                                              start_id = next_id)
  new_rows$resist_score <- rec_score
  if (!is.null(pop$genet_id)) {            # clonal mode: each seedling starts a genet
    new_rows$genet_id <- new_rows$id
    new_rows$origin   <- "seed"
  }

  pop       <- append_rows(pop, new_rows)
  resist_gt <- rbind(resist_gt, rec_gt)

  list(pop = pop, resist_gt = resist_gt, n_rec = n_rec)
}
