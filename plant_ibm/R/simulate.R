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
#
# ARCHIVING DEAD PLANTS (speed)
# -----------------------------
# Dead plants are never needed again inside the loop, but kept in 'pop'
# they make every step scan and copy an ever-growing table (~400,000 rows
# against ~7,000 alive over a 610-year run). When dead rows outnumber the
# living (and exceed 1,000), they are moved, with their genotype rows, to an
# archive; at the end the archive and the living are recombined in id
# order, so the returned individuals and resist_gt are exactly what an
# unarchived run returns. Archiving uses no random numbers and preserves
# the order of living plants, so results are identical with it on or off
# (params$archive_dead = FALSE turns it off).

# Pathogen vertex of the disease triangle for year t: zero before rust
# arrives, then one Beta(annual_env_alpha, annual_env_beta) draw per year,
# independent between years. The place for year-to-year structure in rust
# pressure (design_swappable_interfaces.md section 4.2; notes T12).
# A species' rust_ramp_years (notes T18) scales the draw by
# rust_ramp_factor(), with no extra random numbers, so runs without it are
# unchanged draw for draw.
draw_annual_env <- function(params, t) {
  if (t >= params$rust_start_year) {
    x <- rbeta(1L, params$annual_env_alpha, params$annual_env_beta)
    if (!is.null(params$rust_ramp_years)) x <- x * rust_ramp_factor(params, t)
    x
  } else {
    0
  }
}

# Rust timing per species (docs/design_rust_pressure_modes.md section 3.3;
# notes T18). rust_arrival_year delays the start of rust for the species: the
# run's rust_start_year becomes the later of the two (so a counterfactual with
# rust_start_year = Inf stays rust-free). Applied once at the start of a run,
# so every rust term (hazards, flowering, recruitment, resprout delay, census)
# uses the same start.
apply_rust_timing <- function(params) {
  if (!is.null(params$rust_arrival_year))
    params$rust_start_year <- max(params$rust_start_year, params$rust_arrival_year)
  params
}

# Share of full rust pressure in year t: rising linearly over rust_ramp_years
# after the start (1/ramp in the first year, as fire_change_factor()), then 1.
rust_ramp_factor <- function(params, t) {
  ramp <- params$rust_ramp_years
  if (is.null(ramp) || ramp <= 0 || t < params$rust_start_year) return(1)
  min(1, (t - params$rust_start_year + 1) / ramp)
}

# Flowering status for the year, after fire: alive, not resprouting and old
# enough, unless rust suppresses it. Rust suppression is a per-year
# Bernoulli on eligible plants, driven by each plant's effective
# susceptibility and the year's annual_env_t; off by default (max_effect = 0).
# The place for later flowering rules (notes T12).
# NOTE: This is a per-year draw, so a plant can recover flowering in a
# favourable year (low annual_env_t). Revisit if data suggest flowering
# suppression should be more persistent once infection is established --
# e.g. if heavily infected plants fail to set buds for multiple
# consecutive seasons regardless of annual conditions.
update_flowering <- function(pop, t, params, annual_env_t) {
  eligible <- pop$alive & !pop$resprout & pop$age >= pop$age_first_flower
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
  pop
}

# Clonal mode (docs/design_clonal_stage.md; notes T17). Off unless the species
# sets clonal_mode: genet; R/clonal.R is then needed as well.
is_clonal <- function(params) identical(params$clonal_mode, "genet")

run_simulation <- function(params, init_state = NULL, year0 = 1L, verbose = FALSE) {
  set.seed(params$seed)
  params <- apply_rust_timing(params)
  clonal <- is_clonal(params)
  genets <- NULL

  if (is.null(init_state)) {
    init      <- create_population(params, year0 = year0)
    pop       <- init$pop
    resist_gt <- init$resist_gt
    if (clonal) {
      ci        <- clonal_init(pop, resist_gt, params)
      pop       <- ci$pop
      resist_gt <- ci$resist_gt
      genets    <- ci$genets
    }
  } else {
    pop       <- init_state$individuals
    resist_gt <- init_state$resist_gt
    genets    <- init_state$genets
    if (clonal && is.null(genets)) stop("clonal run continued from a state without genets")
  }
  # States saved before death_cause existed: add it (causes unknown).
  if (is.null(pop$death_cause)) pop$death_cause <- NA_character_

  archive_on  <- !isFALSE(params$archive_dead)
  arch_pop    <- list()
  arch_gt     <- list()
  next_id     <- max(pop$id) + 1L

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
    annual_env_t <- draw_annual_env(params, t)

    # 1. age all living individuals
    pop$age[pop$alive] <- pop$age[pop$alive] + 1L
    if (clonal) alive_start <- pop$alive

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
    # Probability: the cell's rate, times the species' exposure and any
    # fire-change scenario (fire_prob_at(); notes T15).
    fire_this_year <- if (params$fire_prob_annual > 0) {
      rbinom(1L, 1L, fire_prob_at(params, t)) == 1L
    } else {
      t %in% params$fire_years
    }
    # A fire in the cell may miss the stand ("stand" extent mode); the
    # census records whether it reached the simulated population.
    fire_reached <- FALSE
    n_dead_fire  <- 0L
    if (fire_this_year) {
      n_before <- sum(pop$alive)
      pop <- apply_fire(pop, fire_params_at(params, t), t)
      n_dead_fire <- n_before - sum(pop$alive)
      fire_reached <- !identical(attr(pop, "fire_reached"), FALSE)
      attr(pop, "fire_reached") <- NULL
      if (clonal && fire_reached) genets <- clonal_after_fire(pop, genets, params, t)
    }

    # 4. update flowering status, post-fire (rust may suppress it)
    pop <- update_flowering(pop, t, params, annual_env_t)

    # 5. mortality: union of all hazard components
    dead_idx <- nominate_deaths(pop, t, params, annual_env_t)
    pop$alive[dead_idx]      <- FALSE
    pop$death_year[dead_idx]  <- t
    pop$death_cause[dead_idx] <- "mortality"
    pop$flowering[dead_idx]   <- FALSE
    n_dead_other <- length(dead_idx)

    # 5b. clonal: rootstocks of genets with no mature ramet may die
    if (clonal) {
      rs <- clonal_rootstock_step(pop, genets, params, t)
      pop <- rs$pop; genets <- rs$genets
      n_dead_other <- n_dead_other + rs$n_dead
    }

    # 6. recruitment: Beverton-Holt, parent sampling, Mendelian inheritance
    rec       <- recruit(pop, resist_gt, t, params, annual_env_t, next_id = next_id)
    pop       <- rec$pop
    resist_gt <- rec$resist_gt
    next_id   <- next_id + rec$n_rec

    # 6b. clonal: seed recruits start genets; then root suckering
    n_suckers <- 0L
    if (clonal) {
      new_idx <- seq_len(rec$n_rec) + nrow(pop) - rec$n_rec
      genets  <- clonal_add_seed_genets(pop, resist_gt, genets, new_idx)
      sk <- clonal_sucker(pop, resist_gt, genets, t, params, next_id)
      pop <- sk$pop; resist_gt <- sk$resist_gt
      n_suckers <- sk$n_suckers
      next_id <- next_id + n_suckers
    }

    # 7. census
    census_list[[i]] <- census_year(pop, resist_gt, t, rec$n_rec,
                                     n_dead_other + n_dead_fire, params,
                                     fire_this_year = fire_reached,
                                     fire_in_cell   = fire_this_year,
                                     n_dead_fire    = n_dead_fire)
    if (clonal) {
      census_list[[i]] <- cbind(census_list[[i]], clonal_census(pop, genets, n_suckers))
      killed <- which(c(alive_start, logical(nrow(pop) - length(alive_start))) & !pop$alive)
      genets <- clonal_end_year(pop, genets, killed)
    }

    if (verbose && i %% 10 == 0) {
      cat(sprintf("year %d: N=%d  N_IUCN=%d\n",
                  t, sum(pop$alive), sum(pop$alive & pop$flowering)))
    }

    if (sum(pop$alive) == 0 && (!clonal || nrow(genets$df) == 0L)) {
      census_list <- census_list[seq_len(i)]
      break
    }

    # Archive dead rows once they outnumber the living.
    if (archive_on) {
      n_dead <- nrow(pop) - sum(pop$alive)
      if (n_dead > 1000L && n_dead > nrow(pop) - n_dead) {
        dead <- !pop$alive
        arch_pop[[length(arch_pop) + 1L]] <- pop[dead, , drop = FALSE]
        arch_gt[[length(arch_gt) + 1L]]   <- resist_gt[dead, , drop = FALSE]
        pop       <- pop[!dead, , drop = FALSE]
        resist_gt <- resist_gt[!dead, , drop = FALSE]
        rownames(pop) <- NULL
      }
    }
  }

  if (length(arch_pop)) {
    pop       <- do.call(rbind, c(arch_pop, list(pop)))
    resist_gt <- do.call(rbind, c(arch_gt, list(resist_gt)))
    ord       <- order(pop$id)
    pop       <- pop[ord, , drop = FALSE]
    resist_gt <- resist_gt[ord, , drop = FALSE]
  }
  rownames(pop) <- NULL

  census <- do.call(rbind, census_list)
  out <- list(individuals = pop, resist_gt = resist_gt, census = census)
  if (clonal) out$genets <- genets
  out
}
