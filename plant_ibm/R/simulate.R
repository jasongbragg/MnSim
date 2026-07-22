# R/simulate.R
#
# The main annual loop. Every step operates on whole columns/vectors --
# there is no per-individual for-loop here. year0 labels runs with real
# calendar years (e.g. year0 = 2010), so that params$fire_years and
# params$rust_start_year can be specified as actual calendar years.
#
# init_state lets a run be CONTINUED from a prior run's end state
# (calibration → rust projection + counterfactual branching).
#
# ANNUAL PATHOGEN DRAW
# --------------------
# At the start of each year the pathogen vertex of the disease triangle
# is drawn once: annual_env_t ~ Beta(annual_env_alpha, annual_env_beta).
# This single draw is shared by ALL individuals in that year, representing
# site-wide inoculum pressure / weather favourability for rust. It is
# passed to nominate_deaths() and recruit() so that rust effects on
# mortality and fecundity experience the same environmental year.

run_simulation <- function(params, init_state = NULL, year0 = 1L, verbose = FALSE) {
  set.seed(params$seed)

  if (is.null(init_state)) {
    init      <- create_population(params)
    pop       <- init$pop
    resist_gt <- init$resist_gt
  } else {
    pop       <- init_state$individuals
    resist_gt <- init_state$resist_gt
  }

  n_years     <- params$n_years
  census_list <- vector("list", n_years)

  for (i in seq_len(n_years)) {
    t <- year0 + i - 1L

    # --- Annual pathogen vertex of the disease triangle ---------------------
    # Zero before rust arrives; a Beta draw once rust is active.
    # annual_env_alpha = annual_env_beta = 1 → Uniform(0,1): maximum
    # year-to-year uncertainty. Tune alpha/beta to represent site
    # pressure. Passes through to nominate_deaths() and recruit() so
    # both mortality and fecundity experience the same year's conditions.
    annual_env_t <- if (t >= params$rust_start_year) {
      rbeta(1L, params$annual_env_alpha, params$annual_env_beta)
    } else {
      0
    }

    # 1. age all living individuals
    pop$age[pop$alive] <- pop$age[pop$alive] + 1L

    # 2. resprout recovery (from prior fire events)
    rs_idx <- which(pop$alive & pop$resprout)
    if (length(rs_idx) > 0) {
      if (params$resprout_recovery == "geometric") {
        q <- params$resprout_prob_recovery
        recovered <- rs_idx[rbinom(length(rs_idx), 1L, q) == 1L]
      } else {
        pop$resprout_yrs_remain[rs_idx] <- pop$resprout_yrs_remain[rs_idx] - 1L
        recovered <- rs_idx[pop$resprout_yrs_remain[rs_idx] <= 0L]
      }
      if (length(recovered) > 0) {
        pop$resprout[recovered]            <- FALSE
        pop$resprout_yrs_remain[recovered] <- 0L
        pop$resprout_total_yrs[recovered]  <- 0L
      }
    }

    # 3. fire event this year?
    fire_this_year <- if (params$fire_prob_annual > 0) {
      rbinom(1L, 1L, params$fire_prob_annual) == 1L
    } else {
      t %in% params$fire_years
    }
    if (fire_this_year) pop <- apply_fire(pop, params, t)

    # 4. update flowering status, post-fire
    # Basic eligibility: alive, not resprouting, old enough
    eligible <- pop$alive & !pop$resprout & pop$age >= pop$age_first_flower

    # Rust-suppressed flowering: per-year Bernoulli on eligible individuals,
    # driven by each individual's effective_susceptibility and the year's
    # annual_env_t. Off by default (max_effect = 0).
    # NOTE: This is a per-year draw, so a plant can recover flowering in a
    # favourable year (low annual_env_t). Revisit if data suggest flowering
    # suppression should be more persistent once infection is established --
    # e.g. if heavily infected plants fail to set buds for multiple
    # consecutive seasons regardless of annual conditions.
    if (params$rust_flower_dose_response$max_effect > 0 &&
        t >= params$rust_start_year && annual_env_t > 0) {
      eff_susc    <- (1 - pop$resist_score) * pop$env_susceptibility * annual_env_t
      p_suppress  <- vapply(eff_susc, dose_response,
                             numeric(1), cfg = params$rust_flower_dose_response)
      not_suppressed <- rbinom(length(eff_susc), 1L, 1 - p_suppress) == 1L
      pop$flowering <- eligible & not_suppressed
    } else {
      pop$flowering <- eligible
    }

    # 5. mortality: union of all hazard components
    dead_idx <- nominate_deaths(pop, t, params, annual_env_t)
    pop$alive[dead_idx]      <- FALSE
    pop$death_year[dead_idx] <- t
    pop$flowering[dead_idx]  <- FALSE

    # 6. recruitment: Beverton-Holt, parent sampling, Mendelian inheritance
    rec       <- recruit(pop, resist_gt, t, params, annual_env_t)
    pop       <- rec$pop
    resist_gt <- rec$resist_gt

    # 7. census
    census_list[[i]] <- census_year(pop, resist_gt, t, rec$n_rec,
                                     length(dead_idx), params,
                                     fire_this_year = fire_this_year)

    if (verbose && i %% 10 == 0) {
      cat(sprintf("year %d: N=%d  N_IUCN=%d\n",
                  t, sum(pop$alive), sum(pop$alive & pop$flowering)))
    }

    if (sum(pop$alive) == 0) {
      census_list <- census_list[seq_len(i)]
      break
    }
  }

  census <- do.call(rbind, census_list)
  list(individuals = pop, resist_gt = resist_gt, census = census)
}
