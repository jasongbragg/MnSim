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
# Stage 1 -- Impact: each plant is impacted with P = fire_p_fimp.
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
# In both modes, resprout_total_yrs is recorded alongside
# resprout_yrs_remain for canopy_density()'s recovery-progress weighting.

apply_fire <- function(pop, params, t) {
  alive_idx <- which(pop$alive)
  n_alive   <- length(alive_idx)
  if (n_alive == 0) return(pop)

  rust <- rust_modifiers(params, t)

  if (params$fire_p_fimp > 0) {
    # --- Two-stage age-dependent model ---
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
  pop$alive[dead_idx]      <- FALSE
  pop$death_year[dead_idx] <- t

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
