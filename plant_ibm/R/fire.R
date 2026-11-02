# R/fire.R
#
# Fire runs in two stages over all currently alive individuals.
# Two modes, selected by params$fire_p_fimp:
#
# LEGACY MODE (fire_p_fimp = 0, default)
# ---------------------------------------
# Every alive individual faces flat P(die) = fire_kill_prob. Survivors
# always enter resprout. Used by existing scenarios, the Lefkovitch
# matrix-equivalence test, and the matrix builder (build_lefkovitch_fire).
# Preserves exact backward compatibility.
#
# TWO-STAGE AGE-DEPENDENT MODE (fire_p_fimp > 0)
# ------------------------------------------------
# Two variants, chosen by params$fire_extent_mode (docs/methodological_notes.md D3):
#
# "stand" (params.R default) -- a simulated population is a stand (~10^4
#   plants), much smaller than the 0.1-degree cell its fire data describe.
#   Each fire in the cell has its own size f, drawn from a Beta with mean
#   fire_p_fimp (the cell's mean burnt fraction) and sd fire_size_sd (the
#   cell's between-fire sd). The fire reaches the stand with probability f;
#   if it does, each plant burns with probability fire_stand_cover (a fixed
#   high value: unburnt patches), and intensity scales with the size of
#   this fire: old-plant kill ceiling = fire_kill_scalar * f. Big fires
#   are therefore both likelier to reach the stand and more lethal.
#   Over many fires the mean share of plants burnt matches the per-plant
#   variant, but delivered as occasional whole-stand fires.
#
# "plant" (default when fire_extent_mode is absent; the original model) --
#   every fire reaches the population; each plant is impacted with
#   probability fire_p_fimp, and the ceiling is fire_kill_scalar * fire_p_fimp.
#
# Stage 1 -- Impact ("plant" variant): each plant is impacted with P = fire_p_fimp.
#   Plants NOT impacted are completely unaffected (unburned patch).
# Stage 2 -- Kill vs resprout given impact:
#   p_kill(age) = c*p + (1 - c*p) * hill_weight(age, half_sat, hill)
#   where p = fire_p_fimp, c = fire_kill_scalar.
#   At age = 0: p_kill = 1 (seedlings always die if impacted).
#   At age → ∞: p_kill → c*fire_p_fimp (old-plant kill ceiling scales
#   with fire intensity, so hotter fires kill proportionally more large
#   trees, not just more total plants).
#   Non-killed impacted plants → resprout state.
#
# apply_fire() sets attr(pop, "fire_reached") to FALSE when a "stand" fire
# misses the stand; run_simulation() records that in the census.
#
# In both modes, resprout_total_yrs is recorded alongside
# resprout_yrs_remain for canopy_density()'s recovery-progress weighting.

# Size (burnt fraction) of one fire: Beta with the cell's mean and sd.
# sd <= 0 / NA, or mean at 0 or 1, gives the mean itself. The variance is
# capped below its Beta maximum m(1 - m) so the shape parameters stay > 0.
draw_fire_size <- function(m, s) {
  if (m >= 1) return(1)
  if (m <= 0) return(0)
  if (is.null(s) || is.na(s) || s <= 0) return(m)
  v <- min(s^2, 0.9 * m * (1 - m))
  k <- m * (1 - m) / v - 1
  stats::rbeta(1L, m * k, (1 - m) * k)
}

apply_fire <- function(pop, params, t) {
  alive_idx <- which(pop$alive)
  n_alive   <- length(alive_idx)
  if (n_alive == 0) return(pop)

  rust <- rust_modifiers(params, t)

  mode <- if (is.null(params$fire_extent_mode)) "plant" else params$fire_extent_mode

  if (params$fire_p_fimp > 0 && mode == "stand") {
    # --- Two-stage, stand variant ---
    f <- draw_fire_size(params$fire_p_fimp, params$fire_size_sd)
    if (stats::runif(1) >= f) {               # fire in the cell misses the stand
      attr(pop, "fire_reached") <- FALSE
      return(pop)
    }
    cover    <- if (is.null(params$fire_stand_cover)) 1 else params$fire_stand_cover
    impacted <- rbinom(n_alive, 1, cover) == 1
    imp_idx  <- alive_idx[impacted]
    if (length(imp_idx) == 0) return(pop)

    age_i   <- pop$age[imp_idx]
    ceiling <- params$fire_kill_scalar * f
    p_kill  <- ceiling + (1 - ceiling) *
               hill_weight(age_i, params$fire_kill_half_sat, params$fire_kill_hill)

    killed   <- rbinom(length(imp_idx), 1, p_kill) == 1
    dead_idx <- imp_idx[killed]
    surv_idx <- imp_idx[!killed]

  } else if (params$fire_p_fimp > 0) {
    # --- Two-stage age-dependent model ("plant" variant) ---
    impacted  <- rbinom(n_alive, 1, params$fire_p_fimp) == 1
    imp_idx   <- alive_idx[impacted]
    if (length(imp_idx) == 0) return(pop)

    age_i   <- pop$age[imp_idx]
    ceiling <- params$fire_kill_scalar * params$fire_p_fimp
    p_kill  <- ceiling + (1 - ceiling) *
               hill_weight(age_i, params$fire_kill_half_sat, params$fire_kill_hill)

    killed   <- rbinom(length(imp_idx), 1, p_kill) == 1
    dead_idx <- imp_idx[killed]
    surv_idx <- imp_idx[!killed]

  } else {
    # --- Legacy flat-probability model ---
    killed   <- rbinom(n_alive, 1, params$fire_kill_prob) == 1
    dead_idx <- alive_idx[killed]
    surv_idx <- alive_idx[!killed]
  }

  # Dead
  pop$alive[dead_idx]       <- FALSE
  pop$death_year[dead_idx]  <- t
  if (!is.null(pop$death_cause)) pop$death_cause[dead_idx] <- "fire"

  # Survivors → resprout
  n_surv <- length(surv_idx)
  if (n_surv > 0) {
    resist_surv <- pop$resist_score[surv_idx]
    rust_delay  <- (1 - resist_surv) * rust$delay_extra

    delay <- rpois(n_surv, params$resprout_yrs_base) + rpois(n_surv, rust_delay)
    delay <- pmax(delay, 1L)

    pop$resprout[surv_idx]            <- TRUE
    pop$flowering[surv_idx]           <- FALSE
    pop$resprout_yrs_remain[surv_idx] <- delay
    pop$resprout_total_yrs[surv_idx]  <- delay
  }

  pop
}
