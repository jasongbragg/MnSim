# runner/olney_validation.R
#
# Olney SF disease exclusion trial (Carnegie et al. 2016; Rhodamnia rubescens)
# as a validation of accumulating rust damage (docs/design_olney_validation.md;
# notes T19, T20). Needs a species in damage mode (rust_hazard_mode: damage).
#
# DESIGN
#   - Model years start in August: model year t runs Aug t to Aug t + 1, so the
#     end of model year t is compared with the trial's August t + 1 values.
#     Rust (detected Oct 2010) starts in model year 2010; the trial's first
#     assessment (Aug 2011, all arms alike) is the end of model year 2010.
#   - One run per arm (decision 1), all continued from the same stand at the
#     end of model year 2010, with the same yearly pressures (annual_env_series,
#     decision 2) and the same seed each year, so arms differ only by
#     treatment. Fungicide is run-level protection (rust_protection): the
#     year's pressure is multiplied by the unprotected share of the year.
#       untreated          none
#       treated            Aug 2011 - Oct 2014   (model 2011.0 - 2014.17)
#       partially_treated  Aug 2011 - Jun 2013   (model 2011.0 - 2012.83)
#   - Trial trees: 20 living trees of the stand at the start, in an age range
#     (the trial took trees 0.5-4 m; the age-height link is a species prior,
#     so the range is an argument), 10 untreated, 10 treated, of which 5
#     become partially treated.
#   - Crown transparency per tree T = b + (1 - b) D, dead = 1 (decision 4); the
#     arm's value is the mean over its trees. Compared once a year (decision 3).
#   - Stand mortality: share of the stand's trees alive at the start of rust
#     (untreated run) dead by the end of model year 2013 (~Oct 2014; target 53%).
# Fitting goes through the Phase 5 calibration code (CLAUDE.md rule 4).
#
# USAGE (from plant_ibm/, with R/*.R and the species' params loaded):
#   source("runner/olney_validation.R")
#   r <- run_olney_validation(base, n_reps = 20, tree_age = c(a, b), baseline = 0.32)
#   r$trees   # one row per replicate x arm x model year: mean transparency
#   r$stand   # one row per replicate: stand mortality

OLNEY_YEAR_RUST <- 2010        # model year in which rust starts (detected Oct 2010)
OLNEY_YEARS     <- 2011:2013   # trial model years (ends compared with Aug 2012, 2013, 2014)
OLNEY_ARMS <- list(
  untreated         = NULL,
  treated           = c(2011, 2014 + 2 / 12),
  partially_treated = c(2011, 2013 - 2 / 12)
)

# Crown transparency of trees from their rust damage; dead trees count as 1.
crown_transparency <- function(alive, damage, baseline) {
  ifelse(alive, baseline + (1 - baseline) * damage, 1)
}

# One replicate. base: the species' params (damage mode); spinup_years: years
# simulated before rust; tree_age: age range of trial trees at the start.
run_olney_replicate <- function(base, seed, tree_age, baseline, spinup_years = 100L,
                                n_trees = c(untreated = 10L, treated = 5L, partially_treated = 5L)) {
  if (!identical(base$rust_hazard_mode, "damage"))
    stop("The Olney validation needs a species in damage mode (rust_hazard_mode: damage)")
  set.seed(seed)
  years  <- OLNEY_YEAR_RUST:max(OLNEY_YEARS)
  series <- stats::setNames(stats::rbeta(length(years), base$annual_env_alpha, base$annual_env_beta), years)
  seeds  <- sample.int(99999L, spinup_years + length(years) + 1L)

  # Stand to the end of model year 2010 (one rust year)
  p <- base
  p$rust_start_year   <- OLNEY_YEAR_RUST
  p$annual_env_series <- series
  p$fire_years <- integer(0); p$fire_prob_annual <- 0
  p$seed    <- seeds[1]
  p$n_years <- spinup_years + 1L
  start <- run_simulation(p, year0 = OLNEY_YEAR_RUST - spinup_years)
  state <- start[c("individuals", "resist_gt")]
  ind   <- state$individuals

  # Trial trees: living, in the age range; randomly assigned to arms
  pool <- ind$id[ind$alive & ind$age >= tree_age[1] & ind$age <= tree_age[2]]
  if (length(pool) < sum(n_trees)) stop("Only ", length(pool), " living trees aged ",
                                        tree_age[1], "-", tree_age[2], " in the stand")
  picked <- sample(pool, sum(n_trees))
  arm_of <- rep(names(n_trees), n_trees)
  trees  <- list(untreated = picked[arm_of == "untreated"],
                 treated   = picked[arm_of %in% c("treated", "partially_treated")],
                 partially_treated = picked[arm_of == "partially_treated"])
  trees$treated <- setdiff(trees$treated, trees$partially_treated)   # 5 + 5 sprayed to Jun 2013
  # Before Jun 2013 the paper reports the 10 sprayed trees together; here each
  # arm is its own trees throughout.

  at_start <- ind$id[ind$alive]
  rows <- list()
  stand <- NA_real_
  for (arm in names(OLNEY_ARMS)) {
    s <- state
    q <- p
    q$rust_protection <- OLNEY_ARMS[[arm]]
    q$n_years <- 1L
    obs <- function(st, yr) {
      k <- match(trees[[arm]], st$individuals$id)
      tr <- crown_transparency(st$individuals$alive[k], st$individuals$rust_damage[k], baseline)
      data.frame(seed = seed, arm = arm, model_year = yr, obs = sprintf("Aug %d", yr + 1L),
                 mean_transparency = mean(tr), n_alive = sum(st$individuals$alive[k]),
                 n_trees = length(k))
    }
    rows[[length(rows) + 1L]] <- obs(s, OLNEY_YEAR_RUST)
    for (j in seq_along(OLNEY_YEARS)) {
      yr <- OLNEY_YEARS[j]
      q$seed <- seeds[1L + j]                              # same seed per year in every arm
      r <- run_simulation(q, init_state = s, year0 = yr)
      s <- r[c("individuals", "resist_gt")]
      rows[[length(rows) + 1L]] <- obs(s, yr)
    }
    if (arm == "untreated") {
      k <- match(at_start, s$individuals$id)
      stand <- mean(!s$individuals$alive[k])
    }
  }
  list(trees = do.call(rbind, rows),
       stand = data.frame(seed = seed, stand_mortality = stand, n_stand = length(at_start)))
}

run_olney_validation <- function(base, n_reps, tree_age, baseline, master_seed = 20261121L, ...) {
  reps <- lapply(seq_len(n_reps), function(i)
    run_olney_replicate(base, master_seed + i, tree_age, baseline, ...))
  list(trees = do.call(rbind, lapply(reps, `[[`, "trees")),
       stand = do.call(rbind, lapply(reps, `[[`, "stand")))
}

# The trial's August values for comparison (data/validation/; notes T20).
olney_observed <- function(path = file.path(mnsim_root(), "data", "validation",
                                            "carnegie2016_olney_crown_transparency.csv")) {
  d <- utils::read.csv(path, comment.char = "#")
  d <- d[substr(d$date, 6, 7) %in% c("07", "08"), ]      # Aug (the 2012 point is late Jul)
  d$model_year <- as.integer(substr(d$date, 1, 4)) - 1L
  d$mean_transparency <- d$mean_crown_transparency_pct / 100
  d[, c("treatment", "date", "model_year", "mean_transparency")]
}
