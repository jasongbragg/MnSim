# R/matrix_model.R
#
# Deterministic matrix population models that correspond to specific IBM
# parameterisations, used by the matrix-equivalence test suite to validate
# that the IBM's per-year transition mechanics are correctly implemented.
#
# Three levels of equivalence (each corresponds to a test context):
#
#   Level 1 -- geometric: constant Weibull hazard (k=1), AFR=1, no fire,
#     no density dependence. Expected lambda = s*(1+F) in closed form.
#
#   Level 2 -- Leslie: age-structured Weibull (k>1) + delayed flowering
#     (AFR>1). Lambda is the dominant eigenvalue of the Leslie matrix built
#     analytically from the same params.
#
#   Level 3 -- Lefkovitch: 3-stage model (Juvenile, Adult, Resprout) with
#     annual fire and geometric resprout recovery. Because fire stochasticity
#     creates high run-to-run variance in an asymptotic lambda estimate, the
#     equivalence test uses FIXED fire_years and compares the year-by-year
#     IBM trajectory to the deterministic matrix projection step-by-step.
#     This tests the per-year transition mechanics directly, with no
#     dependence on convergence to a stable stage distribution.
#
# Timing convention (post-breeding census, matches IBM annual loop):
#   For an individual at age j in class X at census t:
#     - Ages to j+1 (step 1)
#     - Fire event if applicable (step 3)
#     - Flowers if alive & !resprout & age+1 >= AFR (step 4)
#     - Faces Weibull mortality at age j+1 (step 5)
#     - Produces F offspring per surviving adult (step 6)
#     - Offspring enter age-0 (J class) at census t+1
#
#   CRITICAL timing note for AFR=1 (used in Levels 1 and 3):
#     Age-0 juveniles (J class) age to 1 in the SAME year, which is >= AFR,
#     so they DO flower and contribute F*s offspring -- identical to adults.
#     Therefore M_nf[J, J] = F*s (not 0), and columns J and A are equal.
#     For AFR > 1 (Level 2), age-0 → age-1 < AFR, so J never contributes
#     offspring within its single-year class: M_leslie[1, 1] = 0 by timing.

# ----------------------------------------------------------------------
# 1. Leslie matrix (Levels 1 and 2)
# ----------------------------------------------------------------------

#' Build the post-breeding-census Leslie matrix implied by params.
#' Matrix dimensions: (max_age+1) x (max_age+1), indices corresponding
#' to age classes 0, 1, ..., max_age.
build_leslie <- function(params, max_age = 100) {
  AFR   <- as.integer(round(params$age_first_flower_mean))
  F_per <- params$R_max / params$K_half   # per-adult fecundity, low-N limit
  n     <- max_age + 1L
  L     <- matrix(0, nrow = n, ncol = n)

  for (j in 0:max_age) {
    # survival from age j: individual ages to j+1, faces hazard at j+1
    s_j <- max(1 - weibull_hazard(j + 1, params), 0)

    # fecundity: age-j individual ages to j+1. If j+1 >= AFR it flowers.
    if ((j + 1L) >= AFR) {
      L[1L, j + 1L] <- s_j * F_per
    }

    # survival subdiagonal (or terminal class accumulation)
    if (j < max_age) {
      L[j + 2L, j + 1L] <- s_j
    } else {
      L[n, n] <- s_j
    }
  }
  L
}

# ----------------------------------------------------------------------
# 2. Lefkovitch components (Level 3)
# ----------------------------------------------------------------------

#' No-fire-year 3-stage transition matrix: Juvenile (age-0), Adult (age>=1),
#' Resprout. Requires AFR=1 (one-year juvenile class) and weibull_k=1
#' (constant hazard, so s is the same for all ages and classes).
#'
#' With AFR=1, J individuals (age-0) age to 1 and flower in the same year,
#' so M_nf[J, J] = F*s (identical to M_nf[J, A]). Columns J and A are
#' therefore equal, and the system reduces to a 2-class (N, R) system with
#' dominant eigenvalue s*(1+F) when fire is absent.
build_lefkovitch_nofire <- function(params) {
  s   <- max(1 - weibull_hazard(2, params), 0)
  F   <- params$R_max / params$K_half
  q   <- params$resprout_prob_recovery
  AFR <- as.integer(round(params$age_first_flower_mean))

  # M[to_stage, from_stage], stages: 1=J, 2=A, 3=R
  f_J <- if (AFR <= 1L) F * s else 0  # J→offspring only if AFR=1
  matrix(c(
    f_J,  F * s,   F * s * q,
    s,    s,       s * q,
    0,    0,       s * (1 - q)
  ), nrow = 3L, ncol = 3L, byrow = TRUE)
}

#' Fire-year 3-stage transition matrix. In the IBM, fire (step 3) precedes
#' flowering (step 4), so nobody flowers in a fire year: no offspring.
#' All alive individuals either die (fire_kill_prob) or enter Resprout.
#' Survivors also face normal Weibull mortality (step 5), so the combined
#' per-individual survival rate in a fire year is (1 - fire_kill_prob) * s.
build_lefkovitch_fire <- function(params) {
  s   <- max(1 - weibull_hazard(2, params), 0)
  fkp <- params$fire_kill_prob
  matrix(c(
    0,           0,           0,
    0,           0,           0,
    (1-fkp)*s,   (1-fkp)*s,   (1-fkp)*s
  ), nrow = 3L, ncol = 3L, byrow = TRUE)
}

#' Expected-mixture Lefkovitch matrix, integrating fire and no-fire years
#' with annual fire probability params$fire_prob_annual. Used for computing
#' the asymptotic lambda and for consistency checks (e.g. fire reduces
#' lambda). NOT used for the primary equivalence test (which uses a fixed
#' fire sequence; see project_lefkovitch_seq).
build_lefkovitch <- function(params) {
  p_f  <- params$fire_prob_annual
  M_nf <- build_lefkovitch_nofire(params)
  M_f  <- build_lefkovitch_fire(params)
  (1 - p_f) * M_nf + p_f * M_f
}

#' Project a 3-stage Lefkovitch model year-by-year using the EXACT same
#' fire/no-fire sequence as the IBM, switching between M_nf and M_fire
#' according to fire_years. Returns a numeric vector of total N for each
#' year, starting with N(0) = sum(n0) and ending with N(n_years).
project_lefkovitch_seq <- function(M_nf, M_fire, fire_years, n0, n_years) {
  n   <- as.numeric(n0)
  out <- numeric(n_years + 1L)
  out[1L] <- sum(n)
  for (t in seq_len(n_years)) {
    M         <- if (t %in% fire_years) M_fire else M_nf
    n         <- as.numeric(M %*% n)
    out[t + 1L] <- sum(n)
  }
  out
}

# ----------------------------------------------------------------------
# 3. Generic projection and summary tools
# ----------------------------------------------------------------------

#' Dominant real eigenvalue of matrix M (= asymptotic population growth rate).
dominant_lambda <- function(M) {
  vals <- eigen(M, only.values = TRUE)$values
  max(Re(vals))
}

#' Project a matrix model forward n_years from initial stage vector n0.
#' Returns a data frame: year 0 = initial, year n_years = final.
project_matrix <- function(M, n0, n_years) {
  n_stages <- length(n0)
  out      <- matrix(0, nrow = n_years + 1L, ncol = n_stages + 1L)
  colnames(out) <- c("year", paste0("n", seq_len(n_stages)))
  out[1L, ] <- c(0, n0)
  n <- n0
  for (t in seq_len(n_years)) {
    n            <- M %*% n
    out[t + 1L, ] <- c(t, as.numeric(n))
  }
  as.data.frame(out)
}

#' Estimate the IBM's realised annual growth rate from its census, using
#' year-on-year N_alive ratios over years (burn+1) to n_years. The burn-in
#' period skips transient dynamics from the initial age distribution.
estimate_lambda_ibm <- function(census, burn = 50) {
  N    <- census$N_alive
  n_yr <- length(N)
  if (n_yr <= burn + 1L) return(NA_real_)
  N_use  <- N[(burn + 1L):n_yr]
  ratios <- N_use[-1L] / N_use[-length(N_use)]
  mean(ratios[is.finite(ratios) & ratios > 0], na.rm = TRUE)
}
