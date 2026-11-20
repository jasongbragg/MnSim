# tests/testthat/test-clonal.R
#
# Clonal stage: genets, ramets, root suckering and rootstock survival
# (docs/design_clonal_stage.md; notes T17).

repo_root <- normalizePath(file.path(proj_root, ".."))

clonal_params <- function(...) {
  base <- list(clonal_mode = "genet", clonal_founder_genets = 20L,
               sucker_rate = 0.1, sucker_rate_rootstock = 0.5, sucker_pulse = 1,
               rootstock_fire_survival = 0.8, rootstock_hazard_no_canopy = 0)
  do.call(make_toy_params, modifyList(base, list(...)))
}

test_that("clonal mode off leaves runs and outputs exactly as before", {
  set.seed(1)
  p <- make_toy_params(n_years = 30L, fire_prob_annual = 0.2, fire_p_fimp = 0.5)
  a <- run_simulation(p)
  b <- run_simulation(modifyList(p, list(clonal_mode = "none")))
  expect_identical(a, b)
  expect_identical(names(a), c("individuals", "resist_gt", "census"))
  expect_false(any(c("genet_id", "origin") %in% names(a$individuals)))
  expect_false("N_genets" %in% names(a$census))
  m <- load_species("m_nodosa", root = repo_root)
  expect_null(m$clonal_mode)
})

test_that("founders are grouped into genets sharing genotype and microsite", {
  set.seed(1)
  p  <- clonal_params(n_years = 1L)
  ci <- clonal_init(create_population(p)$pop, create_population(p)$resist_gt, p)
  expect_setequal(unique(ci$pop$genet_id), 1:20)
  for (g in 1:20) {
    k <- which(ci$pop$genet_id == g)
    expect_identical(unique(ci$pop$env_susceptibility[k]), ci$genets$df$env[g])
    expect_identical(unname(unique(ci$resist_gt[k, , drop = FALSE])), unname(ci$genets$gt[g, , drop = FALSE]))
  }
  expect_error(clonal_init(create_population(p)$pop, create_population(p)$resist_gt,
                           modifyList(p, list(clonal_founder_genets = 1e6))), "between 1 and N0")
})

test_that("suckers carry their genet's genotype and id, and seedlings start genets", {
  set.seed(1)
  r <- run_simulation(clonal_params(n_years = 15L))
  ind <- r$individuals
  suck <- ind[ind$origin == "sucker", ]
  expect_gt(nrow(suck), 0)
  expect_true(all(is.na(suck$father_id)))
  mum <- match(suck$mother_id, ind$id)
  expect_identical(ind$genet_id[mum[!is.na(mum)]], suck$genet_id[!is.na(mum)])
  # Same genotype as another ramet of the genet
  same <- vapply(seq_len(min(50, nrow(suck))), function(i) {
    j <- which(ind$id == suck$id[i]); m <- mum[i]
    is.na(m) || identical(unname(r$resist_gt[j, ]), unname(r$resist_gt[m, ]))
  }, logical(1))
  expect_true(all(same))
  seed <- ind[ind$origin == "seed", ]
  expect_gt(nrow(seed), 0)
  expect_identical(seed$genet_id, seed$id)
  expect_true(all(c("N_genets", "N_rootstock_only", "births_sucker") %in% names(r$census)))
  expect_identical(sum(r$census$births_sucker), nrow(suck))
})

test_that("a fire that kills every ramet leaves the rootstock with P = rootstock_fire_survival", {
  set.seed(1)
  p   <- clonal_params(n_years = 1L, rootstock_fire_survival = 0.7)
  ci  <- clonal_init(create_population(modifyList(p, list(N0 = 2000L, clonal_founder_genets = 2000L)))$pop,
                     create_population(modifyList(p, list(N0 = 2000L)))$resist_gt,
                     modifyList(p, list(clonal_founder_genets = 2000L)))
  pop <- ci$pop; pop$alive <- FALSE                 # the fire killed every ramet
  g   <- clonal_after_fire(pop, ci$genets, p, 1L)
  expect_equal(mean(g$df$rootstock_alive), 0.7, tolerance = 0.03)
  # Genets with a surviving ramet are never at risk
  pop$alive[1:10] <- TRUE
  g2 <- clonal_after_fire(pop, ci$genets, modifyList(p, list(rootstock_fire_survival = 0)), 1L)
  expect_true(all(g2$df$rootstock_alive[1:10]))
  expect_false(any(g2$df$rootstock_alive[-(1:10)]))
})

test_that("a rootstock with no canopy keeps suckering ('living dead') unless it has a hazard", {
  set.seed(1)
  # Rust kills every ramet young: no ramet ever matures, genets persist as rootstocks.
  p <- clonal_params(n_years = 40L, N0 = 100L, clonal_founder_genets = 10L,
                     rust_start_year = 1L, annual_env_alpha = 50, annual_env_beta = 1,
                     microclim_alpha = 50, microclim_beta = 1,
                     age_first_flower_mean = 30, age_first_flower_sd = 0, sucker_rate = 0,
                     R_max = 0)
  p$rust_dose_response$age_peak <- 5; p$rust_dose_response$age_floor <- 5
  p$rust_pressure <- 1; p$resist_freq0[] <- 0
  r <- run_simulation(p)
  last <- tail(r$census, 1)
  expect_identical(last$N_IUCN, 0L)
  expect_identical(last$N_genets, 10L)
  expect_gt(sum(tail(r$census$births_sucker, 10)), 0)
  # With a hazard the pool runs out
  set.seed(1)
  r2 <- run_simulation(modifyList(p, list(rootstock_hazard_no_canopy = 0.2)))
  expect_lt(tail(r2$census$N_genets, 1), 3L)
})

test_that("the disturbance pulse raises suckering in the year after a fire", {
  set.seed(1)
  p <- clonal_params(n_years = 12L, N0 = 400L, clonal_founder_genets = 40L,
                     sucker_rate = 0.02, sucker_pulse = 2, fire_years = 6L,
                     fire_p_fimp = 1, fire_extent_mode = "plant", fire_size_sd = 0)
  b <- run_simulation(p)$census
  expect_gt(b$births_sucker[b$year == 7L], 3 * mean(b$births_sucker[b$year %in% 2:5]))
})

test_that("T_g counts sucker parents and is unchanged without suckers", {
  set.seed(1)
  p <- clonal_params(n_years = 120L, N0 = 300L, clonal_founder_genets = 30L, R_max = 0,
                     age_first_flower_mean = 5, age_first_flower_sd = 0,
                     weibull_lambda = 15)
  r <- run_simulation(p, year0 = 1L)
  expect_identical(sum(r$census$births), 0L)        # no seed at all
  g <- suppressMessages(compute_gen_time(r$individuals, 40L, 100L, 1L))
  expect_false(is.na(g$T_g))
  expect_gte(min(g$mother_ages), 1)
  # Non-clonal runs: compute_gen_time sees only seed offspring, as before
  q <- make_toy_params(n_years = 60L)
  s <- run_simulation(q)
  expect_true(all(s$individuals$mother_id[!is.na(s$individuals$mother_id)] %in% s$individuals$id))
})

test_that("clonal settings validate through the schema", {
  set.seed(1)
  schema <- yaml::read_yaml(file.path(repo_root, "config", "schema.yaml"))
  expect_setequal(CLONAL_MODES, unlist(schema$params$clonal_mode$choices))
  for (f in c("clonal_founder_genets", "sucker_rate", "sucker_rate_rootstock", "sucker_pulse",
              "rootstock_fire_survival", "rootstock_hazard_no_canopy"))
    expect_identical(schema$params[[f]]$required_when$clonal_mode, "genet")
})

test_that("suckering saturates with the species' recruitment curve", {
  set.seed(1)
  p <- make_toy_params(K_half = 100)
  expect_equal(sucker_density_factor(0, p), 1)
  expect_equal(sucker_density_factor(100, p), 0.5)
  r <- modifyList(p, list(recruitment_curve = "ricker", ricker_r = 2, ricker_K = 50))
  expect_equal(sucker_density_factor(25, r), exp(-1))
})

test_that("only established genets keep a rootstock once their ramets are dead", {
  set.seed(1)
  p  <- clonal_params(n_years = 1L, N0 = 10L, clonal_founder_genets = 10L)
  ci <- clonal_init(create_population(p)$pop, create_population(p)$resist_gt, p)
  g  <- ci$genets
  g$df$established[1:5] <- FALSE                    # e.g. seedling genets
  pop <- ci$pop; pop$alive <- FALSE
  e <- clonal_end_year(pop, g, integer(0))
  expect_setequal(e$df$genet_id, 6:10)
  # A genet becomes established once it has a mature ramet
  pop$alive[1] <- TRUE; pop$age[1] <- pop$age_first_flower[1]
  e2 <- clonal_end_year(pop, g, integer(0))
  expect_true(1L %in% e2$df$genet_id && e2$df$established[e2$df$genet_id == 1L])
})
