# R/clonal.R
#
# Clonal stage: ramets grouped into genets, root suckering and rootstock
# survival (docs/design_clonal_stage.md; notes T17). Used only when
# params$clonal_mode == "genet"; run_simulation() calls nothing here otherwise,
# so species without it make exactly the same draws as before.
#
# Each row of `pop` is a RAMET (stem), the IUCN mature individual, with a
# genet_id and an origin ("founder", "seed", "sucker"). Ramets of a genet share
# its genotype and microsite. The genet state is
#   genets = list(df = data.frame(genet_id, rootstock_alive, n_killed, env, established),
#                 gt = genotype matrix, one row per df row)
# n_killed: the genet's ramets killed last year (fire, mortality, rootstock);
# it drives next year's sucker pulse. Dead genets are dropped from the table.
#
# A genet is alive while its rootstock is. The rootstock dies
#   - with P = 1 - rootstock_fire_survival when a fire that reaches the stand
#     leaves the genet with no living ramet;
#   - with P = rootstock_hazard_no_canopy in each year the genet has no mature
#     ramet (reserves running down while suckers keep dying).
# Its remaining ramets (suckers attached to it) die with it.
# Only an ESTABLISHED genet (one that has had a mature ramet; founders count)
# keeps a rootstock once all its ramets are dead: a seedling that dies young
# leaves no root system to sucker from, so its genet ends with it.

CLONAL_MODES <- c("none", "genet")

clonal_param <- function(params, name, default) {
  v <- params[[name]]
  if (is.null(v)) default else v
}

# Group the founders into genets: every genet gets at least one ramet, the rest
# are spread at random. Ramets take their genet's genotype and microsite (those
# of its first ramet). Called after create_population()'s own draws, so the
# draws for the founders themselves are those of a non-clonal run.
clonal_init <- function(pop, resist_gt, params) {
  n <- nrow(pop)
  G <- as.integer(params$clonal_founder_genets)
  if (is.na(G) || G < 1L || G > n)
    stop("clonal_founder_genets must be between 1 and N0 (", n, "); got ", params$clonal_founder_genets)
  gid   <- c(seq_len(G), sample.int(G, n - G, replace = TRUE))
  first <- seq_len(G)                       # first ramet of genet g is row g
  pop$env_susceptibility <- pop$env_susceptibility[first[gid]]
  pop$resist_score       <- pop$resist_score[first[gid]]
  resist_gt              <- resist_gt[first[gid], , drop = FALSE]
  pop$genet_id <- gid
  pop$origin   <- "founder"
  genets <- list(
    df = data.frame(genet_id = seq_len(G), rootstock_alive = TRUE, n_killed = 0L,
                    env = pop$env_susceptibility[first], established = TRUE),
    gt = resist_gt[first, , drop = FALSE])
  list(pop = pop, resist_gt = resist_gt, genets = genets)
}

# Row of genets$df for each ramet (NA if its genet has died).
genet_row <- function(pop, genets) match(pop$genet_id, genets$df$genet_id)

# Living and mature living ramets per genet (vectors over genets$df rows).
genet_counts <- function(pop, genets) {
  gi <- genet_row(pop, genets)
  nb <- nrow(genets$df)
  live   <- pop$alive & !is.na(gi)
  mature <- live & pop$age >= pop$age_first_flower
  list(gi = gi,
       n_alive  = tabulate(gi[live],   nbins = nb),
       n_mature = tabulate(gi[mature], nbins = nb))
}

# Kill the living ramets of genets whose rootstock has just died.
kill_genet_ramets <- function(pop, dead_genet_ids, t) {
  idx <- which(pop$alive & pop$genet_id %in% dead_genet_ids)
  if (length(idx)) {
    pop$alive[idx]       <- FALSE
    pop$death_year[idx]  <- t
    pop$death_cause[idx] <- "rootstock"
    pop$flowering[idx]   <- FALSE
  }
  pop
}

# After a fire that reached the stand: genets left with no living ramet keep
# their rootstock with P = rootstock_fire_survival.
clonal_after_fire <- function(pop, genets, params, t) {
  cnt  <- genet_counts(pop, genets)
  bare <- which(genets$df$rootstock_alive & cnt$n_alive == 0L)
  if (length(bare)) {
    lost <- bare[stats::runif(length(bare)) >= params$rootstock_fire_survival]
    genets$df$rootstock_alive[lost] <- FALSE
  }
  genets
}

# Rootstocks of genets with no mature ramet die with P = rootstock_hazard_no_canopy;
# their remaining ramets die with them. Returns list(pop, genets, n_dead).
clonal_rootstock_step <- function(pop, genets, params, t) {
  cnt <- genet_counts(pop, genets)
  at_risk <- which(genets$df$rootstock_alive & cnt$n_mature == 0L)
  n_dead <- 0L
  h <- params$rootstock_hazard_no_canopy
  if (length(at_risk) && h > 0) {
    die <- at_risk[stats::runif(length(at_risk)) < h]
    if (length(die)) {
      genets$df$rootstock_alive[die] <- FALSE
      n0  <- sum(pop$alive)
      pop <- kill_genet_ramets(pop, genets$df$genet_id[die], t)
      n_dead <- n0 - sum(pop$alive)
    }
  }
  list(pop = pop, genets = genets, n_dead = n_dead)
}

# Seed recruits each start a genet (recruit() gives them genet_id = id).
clonal_add_seed_genets <- function(pop, resist_gt, genets, new_idx) {
  if (!length(new_idx)) return(genets)
  genets$df <- rbind(genets$df, data.frame(
    genet_id = pop$genet_id[new_idx], rootstock_alive = TRUE, n_killed = 0L,
    env = pop$env_susceptibility[new_idx], established = FALSE))
  genets$gt <- rbind(genets$gt, resist_gt[new_idx, , drop = FALSE])
  genets
}

# For each requested genet row in `want`, the id of one ramet drawn uniformly
# from that genet's members of `cand` (row indices into pop); NA if it has none.
# Vectorised: candidates sorted by genet, then one uniform draw per request.
pick_in_genet <- function(cand, gi, want, ids) {
  if (!length(want)) return(integer(0))
  nb   <- max(c(gi[cand], want, 0L), na.rm = TRUE)
  srt  <- cand[order(gi[cand])]
  n_g  <- tabulate(gi[cand], nbins = nb)
  off  <- c(0L, cumsum(n_g))[want]
  n_w  <- n_g[want]
  out  <- rep(NA_integer_, length(want))
  ok   <- n_w > 0L
  out[ok] <- ids[srt[off[ok] + floor(stats::runif(sum(ok)) * n_w[ok]) + 1L]]
  out
}

# Density dependence of suckering: the per-capita factor the species'
# recruitment curve applies at N canopy plants, relative to low density, so
# seed and suckers saturate against the same space (notes T17):
#   bevholt  K_half / (K_half + N)        ricker  exp(-ricker_r * N / ricker_K)
sucker_density_factor <- function(N, params) {
  curve <- if (is.null(params$recruitment_curve)) "bevholt" else params$recruitment_curve
  switch(curve,
    bevholt = params$K_half / (params$K_half + N),
    ricker  = exp(-params$ricker_r * N / params$ricker_K),
    stop("No sucker density factor for recruitment curve '", curve, "'"))
}

# Root suckering for year t. Each living genet produces Poisson suckers with mean
#   (sucker_rate x mature ramets + sucker_pulse x ramets killed last year
#    + sucker_rate_rootstock if it has no mature ramet)
#   x (1 - shade suppression) x sucker_density_factor(canopy).
# A sucker's parent (mother_id) is one of the genet's mature ramets, or any of
# its living ramets if none is mature, or NA for a bare rootstock.
# Returns list(pop, resist_gt, n_suckers).
clonal_sucker <- function(pop, resist_gt, genets, t, params, next_id) {
  g   <- genets$df
  cnt <- genet_counts(pop, genets)
  base <- params$sucker_rate * cnt$n_mature +
          params$sucker_pulse * g$n_killed +
          ifelse(cnt$n_mature == 0L, params$sucker_rate_rootstock, 0)
  N_canopy <- canopy_density(pop)
  lam <- base * (1 - shade_suppression(N_canopy, params)) * sucker_density_factor(N_canopy, params)
  lam[!g$rootstock_alive] <- 0
  n_s <- stats::rpois(nrow(g), lam)
  S <- sum(n_s)
  if (S == 0L) return(list(pop = pop, resist_gt = resist_gt, n_suckers = 0L))

  gr <- rep(seq_len(nrow(g)), n_s)          # genet row of each sucker
  live <- which(pop$alive & !is.na(cnt$gi))
  mature_live <- live[pop$age[live] >= pop$age_first_flower[live]]
  # Parent pool per genet: its mature ramets, else all its living ramets.
  use_mature <- cnt$n_mature[gr] > 0L
  mother <- rep(NA_integer_, S)
  mother[use_mature]  <- pick_in_genet(mature_live, cnt$gi, gr[use_mature],  pop$id)
  mother[!use_mature] <- pick_in_genet(live,        cnt$gi, gr[!use_mature], pop$id)

  microsite <- clonal_param(params, "sucker_microsite", "parent")
  env <- if (identical(microsite, "new"))
    stats::rbeta(S, params$microclim_alpha, params$microclim_beta) else g$env[gr]
  gt <- genets$gt[gr, , drop = FALSE]

  new <- data.frame(
    id                  = next_id:(next_id + S - 1L),
    age                 = 0L,
    alive               = TRUE,
    birth_year          = t,
    death_year          = NA_integer_,
    death_cause         = NA_character_,
    founder             = FALSE,
    mother_id           = mother,
    father_id           = NA_integer_,
    age_first_flower    = draw_age_first_flower(S, params),
    resprout            = FALSE,
    resprout_yrs_remain = 0L,
    resprout_total_yrs  = 0L,
    flowering           = FALSE,
    env_susceptibility  = env,
    resist_score        = resist_score_from_gt(gt, params),
    genet_id            = g$genet_id[gr],
    origin              = "sucker")
  if (!is.null(pop$rust_damage)) new$rust_damage <- 0        # damage mode (notes T19)
  list(pop = append_rows(pop, new), resist_gt = rbind(resist_gt, gt), n_suckers = S)
}

# End of year: record each genet's ramets killed this year (for next year's
# pulse), mark genets that now have a mature ramet as established, end
# unestablished genets with no living ramet, and drop dead genets.
clonal_end_year <- function(pop, genets, killed_idx) {
  gi <- match(pop$genet_id[killed_idx], genets$df$genet_id)
  genets$df$n_killed <- tabulate(gi[!is.na(gi)], nbins = nrow(genets$df))
  cnt <- genet_counts(pop, genets)
  genets$df$established <- genets$df$established | cnt$n_mature > 0L
  genets$df$rootstock_alive <- genets$df$rootstock_alive &
    (genets$df$established | cnt$n_alive > 0L)
  keep <- genets$df$rootstock_alive
  list(df = genets$df[keep, , drop = FALSE], gt = genets$gt[keep, , drop = FALSE])
}

# Census columns in clonal mode: living genets, genets with only a rootstock
# left, and suckers produced this year.
clonal_census <- function(pop, genets, n_suckers) {
  cnt <- genet_counts(pop, genets)
  live <- genets$df$rootstock_alive
  data.frame(N_genets = sum(live), N_rootstock_only = sum(live & cnt$n_alive == 0L),
             births_sucker = n_suckers)
}
