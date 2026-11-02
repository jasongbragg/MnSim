# R/individuals.R
#
# The `pop` data frame: one row per individual, dead or alive. Dead
# individuals are NEVER removed -- alive is flipped to FALSE -- so the
# full pedigree and selection history stays queryable. Genotypes live
# in the parallel `resist_gt` matrix (see genetics.R); row i of pop
# always corresponds to row i of resist_gt.
#
# Columns:
#   id, age, alive, birth_year, death_year, mother_id, father_id,
#   age_first_flower, resprout, resprout_yrs_remain, resprout_total_yrs,
#   flowering, resist_score, env_susceptibility, death_cause, founder
#   (birth_year: calendar year for everyone, with age = year - birth_year at
#    the end of each year. Founders get year0 - 1 - age (they age once before
#    their first census, like a recruit born the year before) and
#    founder = TRUE: their ages come from the initial age distribution, not
#    from simulated births.)
#   (Recruits enter at age 0 after the year's mortality step and first face
#    mortality at age 1, so R_max / K_half describe seedlings that survive
#    their first year; docs/methodological_notes.md E2.)
#   (death_cause: NA while alive; "fire" if killed by apply_fire(),
#    "mortality" if killed by nominate_deaths())

#' Create the founding population and its genotypes. year0 = first calendar
#' year of the run (run_simulation passes it); founders' birth_year = year0 - 1 - age.
#' Returns list(pop = data.frame, resist_gt = matrix).
create_population <- function(params, year0 = 0L) {
  n <- params$N0

  age_first_flower <- pmax(round(rnorm(n, params$age_first_flower_mean,
                                           params$age_first_flower_sd)), 1L)

  age <- rgeom(n, 1 / params$weibull_lambda)

  # Environment vertex of the disease triangle (see params.R header):
  # per-individual microsite susceptibility, drawn from Beta distribution,
  # permanently fixed (not heritable). Beta(1,1) = Uniform = no microsite
  # structure. Not zero-initialised -- it should vary even at founding so
  # selection acts on real variation from the start.
  env_susceptibility <- rbeta(n, params$microclim_alpha, params$microclim_beta)

  pop <- data.frame(
    id                   = seq_len(n),
    age                  = age,
    alive                = TRUE,
    birth_year           = as.integer(year0) - 1L - age,
    death_year           = NA_integer_,
    death_cause          = NA_character_,
    founder              = TRUE,
    mother_id            = NA_integer_,
    father_id            = NA_integer_,
    age_first_flower     = age_first_flower,
    resprout             = FALSE,
    resprout_yrs_remain  = 0L,
    resprout_total_yrs   = 0L,
    flowering            = FALSE,
    env_susceptibility   = env_susceptibility
  )

  pop$flowering <- pop$alive & !pop$resprout & pop$age >= pop$age_first_flower

  resist_gt <- init_resist_gt(n, params)
  pop$resist_score <- resist_score_from_gt(resist_gt, params)

  list(pop = pop, resist_gt = resist_gt)
}

#' Build the data-frame rows for n_rec new recruits.
make_recruit_rows <- function(n_rec, t, mother_idx, father_idx, pop, params) {
  if (n_rec == 0) return(NULL)

  start_id <- max(pop$id) + 1L
  age_first_flower <- pmax(round(rnorm(n_rec, params$age_first_flower_mean,
                                           params$age_first_flower_sd)), 1L)

  # Each recruit gets its own env_susceptibility from the same site
  # distribution as founders -- microclimate is determined by where the
  # seed lands, not by parentage.
  env_susceptibility <- rbeta(n_rec, params$microclim_alpha, params$microclim_beta)

  data.frame(
    id                  = start_id:(start_id + n_rec - 1L),
    age                 = 0L,
    alive               = TRUE,
    birth_year          = t,
    death_year          = NA_integer_,
    death_cause         = NA_character_,
    founder             = FALSE,
    mother_id           = pop$id[mother_idx],
    father_id           = pop$id[father_idx],
    age_first_flower    = age_first_flower,
    resprout            = FALSE,
    resprout_yrs_remain = 0L,
    resprout_total_yrs  = 0L,
    flowering           = FALSE,
    env_susceptibility  = env_susceptibility,
    resist_score        = NA_real_   # filled in by recruit() after genotypes are drawn
  )
}


#' Calendar birth years. Current populations store calendar years for
#' everyone (with a founder flag). Older saved populations (no founder
#' column) stored founders as -age relative to the start year; those are
#' converted using year_spinup.
calendar_birth_year <- function(pop, year_spinup) {
  if (!is.null(pop$founder) && !anyNA(pop$founder)) return(pop$birth_year)
  ifelse(pop$birth_year < year_spinup, year_spinup + pop$birth_year, pop$birth_year)
}

#' Founder flag, with the same legacy fallback (founders were the plants
#' whose stored birth year lies before the start year).
is_founder <- function(pop, year_spinup) {
  if (!is.null(pop$founder) && !anyNA(pop$founder)) return(pop$founder)
  calendar_birth_year(pop, year_spinup) < year_spinup
}
