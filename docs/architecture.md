# MnSim architecture

Written 2026-11-19 from the code on branch `multispecies` (commit 9d32066) for
PLAN.md Phase 1. It describes what the code does, not what the README says.
Where the two disagree, see "Documentation that no longer matches the code"
at the end. File:line references are to `plant_ibm/` unless stated otherwise.

---

## 1. Overview

MnSim is an individual-based model with an annual time step. Each plant is one
row of a data frame. Each year the plants age, recover from resprouting,
may burn, flower, die and reproduce. Myrtle rust acts through a "disease
triangle": host resistance × microsite × that year's pathogen pressure.
Above the model sit runners that:

- set up scenarios (calendar years, site inputs),
- assess IUCN Criterion A4 for each run (`runner/iucn_a4.R`),
- compare runs with field data (Pegg et al. 2025: `runner/pegg_validation.R`,
  `runner/pegg_sites.R`),
- run range-wide ensembles over grid cells, seeds and parameter draws
  (`runner/range_ensemble.R`).

There is no package. Scripts `source()` the files in `R/`, a parameter file
that defines `get_default_params()`, and whichever runners they need. The
working directory is always `plant_ibm/`.

```
plant_ibm/
  params.R            working M. nodosa parameters: get_default_params(), RUST_ANCHOR
  R/                  the model (sections 3-4)
    simulate.R          run_simulation(): the annual loop
    individuals.R       founders, recruit rows, birth-year bookkeeping
    genetics.R          resistance genotypes, Mendelian inheritance, resist_score
    mortality.R         hazards, dose_response(), nominate_deaths()
    fire.R              apply_fire(): occurrence handled in simulate.R, effects here
    recruitment.R       Beverton-Holt, parent sampling, canopy shading, recruit()
    census.R            census_year(), pct_decline()
    matrix_model.R      Leslie/Lefkovitch builders, used only by the tests
  runner/             scenarios and assessments (section 6)
  spatial/            per-cell inputs -> spatial/data/mn_site_params.csv (section 7)
  parlib/             alternative get_default_params() files (section 5.2)
  diagnostics/        plausibility plots; source params.R directly
  talk/               figures for talks
  tests/              testthat suite (section 9)
../runscripts/        driver scripts (*.run.code) for the workstation
../docs/methodological_notes.md   assumptions and decisions (IDs A1..V6 used below)
```

---

## 2. Data structures

### 2.1 `pop`: one row per plant, alive or dead

Built by `create_population()` (`R/individuals.R:27`) and extended by
`make_recruit_rows()` (`R/individuals.R:69`) through `append_rows()`
(`R/recruitment.R:7`). Rows are never deleted inside the model: dead plants
get `alive = FALSE`. (For speed they are moved to an archive during the run
and merged back, in id order, at the end; see 3.3.)

| Column | Meaning |
|---|---|
| `id` | unique integer; recruits get ids from a running counter |
| `age` | years; incremented at the start of each year (founders too) |
| `alive` | logical |
| `birth_year` | calendar year. For founders: `year0 - 1 - age` |
| `death_year`, `death_cause` | `NA` while alive; cause is `"fire"` (`apply_fire`) or `"mortality"` (`nominate_deaths`) |
| `founder` | `TRUE` for the initial population; these plants are excluded as mothers from the generation-time estimate |
| `mother_id`, `father_id` | pedigree (`NA` for founders) |
| `age_first_flower` | individual age at first flowering: `max(1, round(N(mean, sd)))`, drawn at birth |
| `resprout`, `resprout_yrs_remain`, `resprout_total_yrs` | post-fire resprouting state and countdown |
| `flowering` | recomputed every year (step 4) |
| `env_susceptibility` | microsite vertex of the disease triangle: Beta(`microclim_alpha`, `microclim_beta`), drawn at birth, not inherited |
| `resist_score` | host vertex: resistance in [0, 1] derived from the genotype |

Life stages are derived from these columns and are not stored:

- juvenile: `age < age_first_flower & !resprout`
- mature: `age >= age_first_flower`
- flowering: `flowering`
- resprouting: `resprout`

### 2.2 `resist_gt`: genotype matrix

An integer matrix with one row per plant and one column per locus, holding
allele doses 0, 1 or 2 (`R/genetics.R`). Row *i* always belongs to the same
plant as row *i* of `pop`, and the two tables are archived and reordered
together. `resist_score = min(1, Σ_j contribution_j)`, where a locus
contributes `effect_j` if the plant has two copies and
`effect_j × dominance_j` if it has one (`resist_score_from_gt`,
`R/genetics.R:73`). Loci are unlinked and inherited by Mendelian segregation
(`inherit_alleles`, `R/genetics.R:41`).

### 2.3 `census`: one row per simulated year

From `census_year()` (`R/census.R:18`):

| Column | Meaning |
|---|---|
| `year` | calendar year |
| `N_alive` | all living plants |
| `N_IUCN` | mature individuals: alive with `age >= age_first_flower`, **including** plants kept from flowering by fire or rust. This is the headline count (notes A4) |
| `N_flowering` | the stricter count, kept as an alternative |
| `N_juvenile`, `N_resprout` | stage counts |
| `births`, `deaths`, `deaths_fire` | `deaths` includes fire deaths |
| `mean_resist_score`, `freq_locus*` | selection trajectory |
| `fire_event` | fire reached the simulated stand |
| `fire_in_cell` | fire occurred in the grid cell |
| `rust_active`, `extinct` | flags |

`run_simulation()` returns `list(individuals = pop, resist_gt, census)`.
If the population dies out, the loop stops and the census ends at the
extinction year.

---

## 3. The annual loop: `run_simulation()` (`R/simulate.R:32`)

```
set.seed(params$seed)                                   # always; line 33
pop/resist_gt <- create_population(params, year0)       # or init_state
for t in year0 .. year0 + n_years - 1:
  0. annual_env_t ~ Beta(annual_env_alpha, annual_env_beta) if t >= rust_start_year, else 0
  1. age += 1 for living plants
  2. resprout recovery: countdown (default) or geometric
  3. fire? rbinom(fire_prob_annual) if > 0, else t %in% fire_years -> apply_fire()
  4. flowering <- alive & !resprout & age >= age_first_flower
       (minus rust suppression if rust_flower_dose_response$max_effect > 0)
  5. deaths <- nominate_deaths(pop, t, params, annual_env_t)
  6. recruits <- recruit(pop, resist_gt, t, params, annual_env_t)
  7. census_year(...)
  stop if extinct; archive dead rows if they outnumber the living (and > 1000)
merge archive, order by id
```

### 3.1 Consequences of the step order

- Fire comes **before** flowering and mortality. Plants that survive a fire
  are set to resprout and do not flower that year. They then face that
  year's mortality step, including the extra rust mortality of resprouting
  plants (Pegg timing, `runner/pegg_validation.R` header).
- Recruits are added **after** mortality, at age 0. They first face mortality
  at age 1, so `R_max`/`K_half` describe seedlings that survive their first
  year (notes E2).
- Flowering is decided before mortality, but plants that die are set to
  `flowering = FALSE` (`R/simulate.R:135`). `recruit()` samples parents from
  `alive & flowering`, so only plants that survive the year's mortality
  reproduce.
- One pathogen draw per year (`annual_env_t`) is shared by every plant and
  by mortality, flowering suppression and fecundity suppression.

### 3.2 Random-number streams (matters for exact-match regression tests)

Each run calls `set.seed(params$seed)` and draws in a fixed order:

1. Founders: `rnorm` (age at first flowering), `rgeom` (age), `rbeta`
   (microsite), then two `rbinom` per locus (genotypes).
2. Each year:
   - `rbeta` (pathogen, only once rust is on);
   - `rbinom` (geometric recovery only);
   - `rbinom` (fire occurrence, only if `fire_prob_annual > 0`);
   - fire draws: `rbeta` size, `runif` reach, `rbinom` cover and kill,
     `rpois` delays;
   - `rbinom` (flowering suppression, if on);
   - `rbinom(n_alive)` (deaths);
   - `rpois` (number of recruits);
   - `sample.int`/`runif` (parents and selfing);
   - `rbinom` per locus (inheritance; drawn even for homozygous parents,
     because `ifelse` evaluates both branches);
   - `rnorm`, `rbeta` (recruit attributes).

Consequences:

- The rust and no-rust runs use identical draws up to the year rust arrives
  (common random numbers; this is tested).
- Any refactor that adds, removes or reorders a draw changes every later
  output, even when the change is statistically neutral. Exact-match
  regression tests will catch such changes, and that is intended.
- Archiving uses no random numbers (tested to give identical output).
- Runners add draws **before** `run_simulation()` re-seeds: the range
  ensemble's parameter draws and the `rust_site_sd` site effect. These do not
  affect the model's own stream.

### 3.3 Archiving (`R/simulate.R:160`)

Dead rows (and their genotype rows) are moved out once there are more than
1000 of them and they outnumber the living. At the end they are restored and
everything is put back in id order. `params$archive_dead = FALSE` turns this
off. The output is identical either way.

---

## 4. Modules

### 4.1 Mortality (`R/mortality.R`)

`nominate_deaths()` (`:161`) gives every living plant one probability of
death and draws all deaths in a single `rbinom`:

```
p_death = 1 - (1 - p_age)(1 - p_sen)(1 - p_jd)(1 - p_rust)
```

| Term | Function | Applies to | Formula |
|---|---|---|---|
| `p_age` | `weibull_death_prob` (`:32`) | all | `1 - exp(-(H(a) - H(a-1)))`, `H(x) = (x/λ)^k`. With k = 1 this is constant: `1 - exp(-1/λ)` |
| `p_sen` | `senescence_hazard` (`:140`) | all | `dose_response(age, senescence_dose_response)`, a rising Hill curve |
| `p_jd` | `juv_decline_hazard` (`:108`) | juveniles (not resprouting) | `dose_response(N_canopy, juv_decline_dose_response) × hill_weight(age, juv_decline_age_half_sat, juv_decline_age_hill)` |
| `p_rust` | inline (`:198-215`) | all, once rust is on and `annual_env_t > 0` | `eff_susc × (rust_age_h + rust_resprout_h)` |

The rust term:

```
eff_susc        = (1 - resist_score) × env_susceptibility × annual_env_t
rust_age_h      = age_floor + (rust_pressure × age_peak - age_floor) × hill_weight(age, age_half_sat, age_hill)
rust_resprout_h = resprouting ? dose_response(rust_pressure, rust_dose_response$resprout) : 0
```

Density dependence acts only through canopy density, in `p_jd` and in
recruitment shading. There is no separate density-dependent mortality term:
`dd_hazard` and `dd_alpha` no longer exist.

`dose_response(pressure, cfg)` (`:52`) is the general response curve used
throughout the model. `cfg` has `max_effect`, `form`, `half_sat` and `hill`;
the forms are `linear`, `power`, `saturating`, `sigmoid` and `threshold`.
`hill_weight(x, half_sat, hill)` (`:43`) is the falling counterpart:
1 at 0, ½ at `half_sat`.

`rust_modifiers()` (`:91`) and `weibull_hazard()` (`:18`) are no longer used
by the model. They are kept for diagnostics and tests.

### 4.2 Fire (`R/fire.R:62`)

Whether a fire happens is decided in `simulate.R` (step 3). `apply_fire()`
then works in one of three modes:

| Mode | Selected by | Behaviour |
|---|---|---|
| legacy | `fire_p_fimp == 0` | every living plant dies with `fire_kill_prob`; survivors resprout. Used by the Lefkovitch matrix test and the unit-test parameter set |
| plant | `fire_p_fimp > 0`, `fire_extent_mode` absent or `"plant"` | each plant burns with probability `fire_p_fimp`; `p_kill(age) = c + (1-c)·hill_weight(age, fire_kill_half_sat, fire_kill_hill)`, `c = fire_kill_scalar × fire_p_fimp` |
| stand (params.R default) | `fire_p_fimp > 0`, `fire_extent_mode = "stand"` | fire size `f ~ Beta(mean fire_p_fimp, sd fire_size_sd)` (`draw_fire_size`, `:53`); the fire reaches the stand with probability `f`; each plant burns with probability `fire_stand_cover`; kill as above with `c = fire_kill_scalar × f` (notes D3) |

Burnt plants that survive go into resprout for
`base + Poisson(rust delay)` years, with a minimum of 1 (`:128`):

- `base` is fixed at `resprout_yrs_base` (`resprout_delay_dist = "fixed"`,
  notes D3), or is Poisson-distributed if that field is absent.
- The rust delay is `rust_resprout_delay()` (`R/mortality.R:84`):
  `dose_response(rust_pressure × (1 - resist) × microsite × mean pathogen,
  rust_dose_response$delay)`, and 0 before rust arrives.

Fire affects juveniles and adults alike. A juvenile that survives also
resprouts, but it does not count as canopy (see 4.3).

### 4.3 Recruitment (`R/recruitment.R:96`)

```
n_flowering_eff = Σ over flowering plants of (1 - rust fecundity suppression)   [suppression off by default]
expected        = bevholt(n_flowering_eff) × (1 - shade_suppression(N_canopy))
n_rec           ~ Poisson(expected)
```

- `bevholt(N) = R_max·N / (K_half + N)` (`:17`).
- `ricker()` (`:22`) is defined but not wired in: `recruit()` always calls
  `bevholt()`.
- `canopy_density()` (`:64`) is the number of living plants of flowering age
  that are not resprouting, plus resprouting plants of flowering age weighted
  by how far through recovery they are (`1 - remain/total`).
- `shade_suppression()` (`:78`) is `dose_response(N_canopy,
  shade_dose_response)`, capped at 1.
- Parents (`sample_parents`, `:29`) are drawn with replacement from flowering
  plants. The mother is chosen uniformly. The father is the mother itself
  with probability `selfing_rate / (n_flower - 1 + selfing_rate)`, otherwise
  a different flowering plant chosen uniformly.
- Recruits inherit genotypes from their parents but draw age at first
  flowering and microsite afresh.

### 4.4 Flowering (`R/simulate.R:107`)

A plant flowers if it is alive, not resprouting and has reached its own
`age_first_flower`. If `rust_flower_dose_response$max_effect > 0`, flowering
is also suppressed by a yearly Bernoulli draw with probability
`dose_response(eff_susc, ...)`. This is off by default.

### 4.5 Genetics (`R/genetics.R`)

Founding genotypes are drawn under Hardy–Weinberg at `resist_freq0`.
`n_loci = length(resist_locus_effect)`. Genotypes reach the rest of the model
only through `resist_score`. In the working parameters the population is
fully susceptible (`resist_freq0 = 0`; notes E3; real data in `private/`,
which must not be used).

### 4.6 Census and decline (`R/census.R`)

See 2.3. `pct_decline(census, y0, y1)` gives a simple two-year decline in
`N_IUCN`. The IUCN assessment itself is in `runner/iucn_a4.R`.

### 4.7 Founders (`R/individuals.R:27`)

- `N0` plants.
- Ages drawn `rgeom(1/weibull_lambda)`. This is a structural assumption,
  close to the stable age distribution only for a constant hazard.
- Microsite drawn from the site Beta; genotypes at `resist_freq0`.

Founders are only a starting point: runs spin up for centuries (section 6),
and generation time ignores founder mothers.

---

## 5. Parameters

### 5.1 Where parameter values come from

A run's parameter list is built in layers. Later layers override earlier ones.

1. **`get_default_params()`.** Whichever file was sourced last defines it:
   - `params.R`: working M. nodosa values, used by every runner and
     diagnostic;
   - `parlib/unittests/params.R`: used by the tests;
   - a `parlib/ensemble/...` file: used for Morris screening.
2. **Site overrides: `build_site_params()` (`runner/site_params.R:51`).**
   Values from one row of `spatial/data/mn_site_params.csv`:
   - `fire_prob_annual` (× `fire_scale`), `fire_p_fimp`, `fire_size_sd` for
     the chosen fire regime (`wild`/`all`);
   - `annual_env_alpha/beta` (rust pressure).

   `fire_years` is cleared. The resistance hook is commented out.
3. **Persistent site effect** (`apply_rust_site_effect`,
   `runner/site_params.R:98`). Shifts the logit of the mean pathogen
   pressure by z ~ N(0, `rust_site_sd`). Off when `rust_site_sd = 0`.
4. **Ensemble parameter draws** (`draw_params`, `runner/range_ensemble.R:57`).
   Currently only `rust_dose_response$resprout$max_effect` ~ U(0.32, 0.78).
5. **Scenario settings set by the runner:**
   - `rust_start_year`, `n_years` and `fire_years = integer(0)`
     (`run_iucn_a4`);
   - a single fire in 2020 with `fire_p_fimp = 1` and rust from 2010
     (Pegg scenarios);
   - `seed`.
6. **Flat vectors.** `vec_to_params()` (`runner/params_utils.R`) writes
   values into the nested list using `a__b__c` keys (Morris, `ensemble.R`).

`make_run_id()` hashes `deparse(params)` to name cache files. Any change to
the list (values, element order, or integer vs double type) gives a new id.

### 5.2 `parlib/`

| Path | Used by | State |
|---|---|---|
| `unittests/params.R` | `tests/testthat/helper-fixtures.R` (instead of `params.R`) | Turns everything off: legacy fire, Poisson resprout window of 3, Weibull k = 2, λ = 30, rust Beta(1,1), microsite Beta(1,1), `resist_freq0 = 0.05`, fixed fire years 40/80/120. Same schema as `params.R` since 2026-11-19 (checked by `test-params.R`; notes T5) |
| `diags/*/params.R` | — | Archived 2026-11-19 to `archive/parlib_diags_params.R` (five identical, unused copies of an older calibration) |
| `ensemble/no_rust_no_fire/params_base.R` + `param_ranges.R` | `runner/morris.R` via `runscripts/morris_sensitivity.run.code` | No-disturbance baseline plus prior ranges (`keys`/`mins`/`maxs`/`Sigma`) for 10 demographic parameters, each with a written rationale |

### 5.3 Parameter reference: where each parameter is used

Values: P = `params.R` (working M. nodosa); U = `parlib/unittests/params.R`.

| Parameter | P | U | Used in |
|---|---|---|---|
| `N0` | 10000 | 10000 | `R/individuals.R:28` |
| `n_years` | 200 | 200 | `R/simulate.R:51`; overridden by `runner/iucn_a4.R:317`, `runner/pegg_validation.R:208` |
| `seed` | 42 | 42 | `R/simulate.R:33`; set per run by the runners |
| `weibull_k` | 1 | 2 | `R/mortality.R:19,33` |
| `weibull_lambda` | 60 | 30 | `R/mortality.R:20,33`; founder ages `R/individuals.R:33` |
| `senescence_dose_response` | sigmoid, max 0.1, half_sat 40, hill 8 | off | `R/mortality.R:141` |
| `juv_decline_dose_response` | saturating, max 0.5, half_sat 5000 | off | `R/mortality.R:109` |
| `juv_decline_age_half_sat`, `_hill` | 2, 2 | 2, 1 | `R/mortality.R:110-111` |
| `R_max`, `K_half` | 1000, 4000 | 4000, 5000 | `R/recruitment.R:18`; `R/matrix_model.R:49,86` |
| `shade_dose_response` | off (max 0) | off | `R/recruitment.R:79` |
| `ricker_r`, `K` | absent | absent | `R/recruitment.R:23` (not wired in) |
| `age_first_flower_mean`, `_sd` | 5, 1.5 | 5, 1.5 | `R/individuals.R:30-31,74-75`; `R/matrix_model.R` |
| `selfing_rate` | 0.10 | 0.10 | `R/recruitment.R:40` |
| `fire_years` | none | 40, 80, 120 | `R/simulate.R:93` (used only when `fire_prob_annual == 0`) |
| `fire_prob_annual` | 1/50 | 0 | `R/simulate.R:90-91`; overridden per cell (`runner/site_params.R:57`) |
| `fire_p_fimp` | 0.6 | 0 | `R/fire.R:70-98` (selects the mode); per cell (`runner/site_params.R:58`); 1 in Pegg runs |
| `fire_extent_mode` | `"stand"` | absent → plant | `R/fire.R:68` |
| `fire_size_sd` | 0.146 | absent | `R/fire.R:72`; per cell (`runner/site_params.R:62`) |
| `fire_stand_cover` | 0.9 | absent → 1 | `R/fire.R:77` |
| `fire_kill_prob` | 0.35 | 0.35 | `R/fire.R:108` (legacy mode only) |
| `fire_kill_scalar`, `_half_sat`, `_hill` | 0.30, 2, 2 | 0.30, 5, 2 | `R/fire.R:83-100` |
| `resprout_yrs_base` | 2 | 3 | `R/fire.R:126` |
| `resprout_delay_dist` | `"fixed"` | absent → Poisson | `R/fire.R:125` |
| `resprout_recovery`, `resprout_prob_recovery` | countdown, 1/3 | countdown, 1/3 | `R/simulate.R:75-76` |
| `rust_start_year` | Inf (runners set 2010) | Inf | `R/simulate.R:63,120`; `R/mortality.R:85,92,199`; `R/recruitment.R:105`; `R/census.R:41` |
| `rust_pressure` | 1 | 1 | `R/mortality.R:88,206,212` (scales the age peak, the resprout term and the delay; not set per cell) |
| `microclim_alpha`, `_beta` | 20, 5 | 1, 1 | `R/individuals.R:40,80` |
| `annual_env_alpha`, `_beta` | 14, 6 (anchor μ 0.70 × conc 20) | 1, 1 | `R/simulate.R:64`; `R/mortality.R:86`; per cell (`runner/site_params.R:66-67`) |
| `rust_site_sd` | 0 | absent | `runner/site_params.R:91`; `runner/pegg_sites.R:70` |
| `rust_dose_response$age_*` | peak 0.3, floor 0.025, half_sat 2, hill 4 | 0.2, 0, 3, 2 | `R/mortality.R:204-208` |
| `rust_dose_response$resprout` | saturating, max 0.53, half_sat 0.2 | 0.30, 0.5 | `R/mortality.R:212`. Calibrated to Pegg (V4); the range ensemble draws U(0.32, 0.78) |
| `rust_dose_response$delay` | linear, max 1.79 | linear, 2 | `R/mortality.R:88` (via `rust_resprout_delay`) |
| `rust_flower_dose_response` | off | off | `R/simulate.R:119,123` |
| `rust_recruit_dose_response` | off | off | `R/recruitment.R:104,110` |
| `resist_locus_effect`, `_dominance`, `_freq0` | 1, 1, 0 | 1, 1, 0.05 | `R/genetics.R:22-27,79-81` |
| `archive_dead` | absent (on) | absent (on) | `R/simulate.R:46` |
| `RUST_ANCHOR` (global, not in the list) | Yarringully, μ 0.70, conc 20, cap 0.95 | — | `params.R:24`; `spatial/rust_pressure.R`; `runner/site_params.R:25` |

### 5.4 Constants in the code that are not parameters

These are a starting point for the Phase 2 inventory, not the full list:

- `R/individuals.R`
  - `:30`, `:74`: age at first flowering is at least 1.
  - `:33`: founder ages are geometric with rate `1/weibull_lambda`.
- `R/fire.R`
  - `:57`: the Beta variance of fire size is capped at 0.9·m(1−m).
  - `:128`: the resprout window is at least 1 year.
- `R/mortality.R:21`: age floor of 0.01 in `weibull_hazard`.
- `R/simulate.R:163`: archive threshold of 1000.
- Runner settings that are specific to M. nodosa or to the assessment:
  - calendar phasing 1600 / 1810–1910 / 2010 / 2026 / 2210 and 11-year
    averaging (`runner/iucn_a4.R:164`, `runner/range_ensemble.R:29-35`);
  - Pegg scenario 1910 / 2010 / 2020 / 2022, observed 0.48 of 50
    (`runner/pegg_validation.R:185`);
  - `SITE_PARAMS_PATH = "spatial/data/mn_site_params.csv"`;
  - `QLD_LAT`;
  - `DEFAULT_PARAM_DRAWS`;
  - the plot table `runner/data/pegg2025_mnodosa_plots.csv`.

---

## 6. Runners and workflows

### 6.1 IUCN A4: `run_iucn_a4()` (`runner/iucn_a4.R:284`)

1. Run one continuous simulation from `year_spinup` (1600) to `year_end`
   (2210), with rust from 2010. Optionally run a counterfactual with no rust
   from the same seed.
2. `compute_gen_time()` (`:90`) estimates T_g as the mean age of mothers at
   births in 1810–1910. Only mothers that are dead and were born during the
   simulation (not founders) are counted.
3. `assess_a4()` (`:206`), a pure function, computes the decline:
   - window L = max(10, round(3·T_g));
   - reference = mean `N_IUCN` over the 11 years ending in 2009;
   - comparison = mean over the 11 years ending in 2009 + L;
   - the comparison end is capped at the present (2026) + 100;
   - years after extinction count as N = 0.

   The same window is also assessed on `N_flowering`.
4. `iucn_category()` maps the decline to a category: ≥ 80 % CR, ≥ 50 % EN,
   ≥ 30 % VU, otherwise LC/NT.

Results are cached by `iucn_run_id(params, assess)`.
`run_iucn_seeds()` / `plot_iucn_ensemble()` run and plot many seeds.
There is no A3 code: `collect_iucn_results()` skips results that are not A4.

### 6.2 Pegg et al. 2025 validation

- `runner/pegg_validation.R`: one plot (Yarringully).
  - Spin-up from 1910; rust from 2010; one fire in model year 2020 burning
    every plant; the run ends in 2021.
  - `compute_pegg_stats()` (`:117`) computes the fraction of resprouters dead
    at the survey, using `death_cause`.
  - Site parameters come from the anchor cell (`pegg_site_params`).
- `runner/pegg_sites.R`: seven plots, mortality per observation interval.
  - Likelihood: a mixture over seeds, where each seed stands for one simulated
    site; this lets the fit reward spread between sites.
  - Each spin-up to 2019 is shared, then branched for every `max_effect`
    value.
  - Notes V6; this is ROADMAP Phase 2.

### 6.3 Range-wide ensemble: `runner/range_ensemble.R`

Run *i* uses its own RNG stream (`master_seed + i`). For each run it:

- draws a cell with probability proportional to abundance (`weight_n_raw`),
  then calls `build_site_params()`;
- draws a seed and parameter values;
- runs `run_iucn_a4()`.

A plan table and per-run `.rds` files make runs restartable. Fire scenarios
(`fire_regime`, `fire_scale`) reuse the same plan, so runs pair one-to-one.
Results are summarised as category shares overall and by region (NSW / Qld).

### 6.4 Other runners

- `runner/morris.R`: Morris screening on `summarise_run()` statistics.
- `runner/ensemble.R`: prior draws (LHS, random or copula) and batch runs,
  intended for ABC.
- `runner/summarise.R`: `summarise_run()`, which reduces a run to scalar
  summaries.
- `archive/main.R` (archived 2026-11-19): the old two-phase example (calibrate, then run rust and
  counterfactual from `init_state`).

`../runscripts/*.run.code` are hand-run driver scripts for the workstation.
They call `/opt/R/4.3.1/bin/Rscript`.

---

## 7. Spatial inputs (`spatial/`)

These scripts build `spatial/data/mn_site_params.csv`, one row per occupied
0.1° cell (375 cells). Each row holds the fire regime per source and regime,
the rust Beta parameters and the abundance weight.

- Occurrences: `01_get_melaleuca_occurrences.R`.
- Fire: `spatial/fire/`, from NSW NPWS and Qld fire-scar records. It gives
  the burnt fraction per cell-year and a negative-binomial prior fitted by
  maximum likelihood (notes D1–D5).
- Rust: `rust_pressure.R`. A cell's mean pressure is
  `μ · SDM_cell / SDM_anchor-cell`, capped at 0.95, so the anchor cell
  always gets μ = 0.70.

Everything here is specific to M. nodosa: the occurrence grid, the anchor,
and file names starting `mn_`.

---

## 8. What a multispecies version touches

Biology lives entirely in the parameter list, apart from the structural
constants in 5.4. To change species you need to change:

- `get_default_params()`;
- `RUST_ANCHOR` and the site table (`SITE_PARAMS_PATH`);
- the scenario years in the runners;
- the Pegg-specific validation code.

The model functions in `R/` read nothing from the global environment except
`params`. (`RUST_ANCHOR` is read only by `params.R` itself, the spatial code
and `anchor_site_row()`.)

Mechanisms that already exist but are switched off, and are relevant to
PLAN Phase 3:

- rust suppression of flowering (`rust_flower_dose_response`);
- rust suppression of fecundity (`rust_recruit_dose_response`);
- canopy shading of establishment (`shade_dose_response`);
- Ricker recruitment (defined, not wired in);
- three fire modes;
- resistance architecture with any number of loci.

Not present:

- a clonal or root-sucker stage;
- damage or defoliation history per plant (hazard is memoryless except for
  the resprout state);
- a constant rust-pressure mode separate from the SDM site table;
- an A3 window;
- committee-style projection.

---

## 9. Tests (`plant_ibm/tests/`)

Run from `plant_ibm/` with R 4.3.1:

```
cd plant_ibm && /opt/R/4.3.1/bin/Rscript tests/testthat.R
```

All 16 files pass (2026-11-19), in about 25 s. The system `Rscript` is
R 3.4.4 and fails because it cannot load `sf`. `renv` is not activated when
running from `plant_ibm/`, because `.Rprofile` is at the repo root.

`helper-fixtures.R`:

- sources `parlib/unittests/params.R`, not `params.R`, plus all of `R/` and
  most runners;
- defines the toy and matrix-equivalence fixtures.

Coverage, by file:

- `test-genetics`, `test-mortality`, `test-fire`, `test-recruitment`,
  `test-individuals`, `test-census`: unit and statistical checks of each
  function.
- `test-matrix-equivalence`: the IBM reduces to geometric, Leslie and
  Lefkovitch models with known answers.
- `test-simulate`: self-consistency. Same seed gives identical output; CRN
  between rust and counterfactual runs; `init_state` resume; archive on/off
  identical; extinction.
- `test-iucn-a4`: `compute_gen_time`, `assess_a4` windows and cap,
  categories, run ids, summaries.
- `test-pegg`, `test-pegg-sites`, `test-site-params`, `test-fire-regime`,
  `test-range-ensemble`: runner and spatial logic.
- `test-params`: the unit-test parameter set has everything off.

Only `test-site-params.R` and `test-pegg-sites.R:32` read the working
`params.R`, and they check values, not model output.

**No test compares output with stored results**, and `tests/reference/` (at
the repo root) is empty. The working M. nodosa configuration is not run end
to end by any test. That configuration is stand fire mode, a fixed resprout
window, senescence, juvenile decline, anchor rust Beta(14, 6) and microsite
Beta(20, 5).

---

## 10. Documentation that no longer matches the code

- `plant_ibm/README.md`:
  - says the IUCN N is `alive & flowering`; the code uses reproductive age
    (`N_IUCN`, notes A4);
  - describes `dd_hazard`, `dd_alpha` and `K` density dependence, which has
    been removed;
  - describes juvenile decline with `exp(-age/tau)`; the code uses
    `hill_weight`;
  - says legacy fire is the default; `params.R` uses stand mode;
  - says a resprouting adult counts fully as canopy; the code weights it by
    recovery progress;
  - lists `parlib` nowhere.
- `plant_ibm/TESTING.md`:
  - says "149 tests across 9 files, ~6.8 s"; there are 16 files and the run
    takes ~25 s;
  - lists `dd_hazard` tests that no longer exist.
- `R/fire.R` header: calls legacy mode the "default", which is true only for
  the unit-test parameters.
- `R/mortality.R`, comment above `rust_resprout_delay`: two docstrings run
  together (one for `rust_modifiers`, one for `rust_resprout_delay`).
- `parlib/unittests/params.R`: claims the same schema as `params.R` (see 5.2).
- `parlib/diags/*/params.R`: all five copies are headed
  `parlib/diags/demography/params.R`.
