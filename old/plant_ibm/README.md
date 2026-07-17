# Plant Population IBM: Myrtle Rust × Fire Disturbance

An individual-based model (IBM) in R projecting plant population dynamics
under myrtle rust (*Austropuccinia psidii*, arrived in Australia 2010) and
fire disturbance, built for an IUCN Red List Criterion A listing assessment.
The model tracks a full pedigree, genotype-dependent rust resistance, and
the interaction between fire and rust (a survivor's vulnerable post-fire
resprout window is lengthened by disease pressure, scaled down by their own
resistance genotype).

For test suite documentation, see **TESTING.md**.

## Contents

```
plant_ibm/
├── R/
│   ├── individuals.R     # create_population(), make_recruit_rows()
│   ├── genetics.R        # init_resist_gt(), inherit_alleles(),
│   │                     #   resist_score_from_gt(), allele_freqs()
│   ├── mortality.R       # weibull_hazard(), senescence_hazard(), hill_weight(),
│   │                     #   dd_hazard(), dd_age_weight(), dose_response(),
│   │                     #   rust_modifiers(), juv_decline_hazard(),
│   │                     #   nominate_deaths()
│   ├── fire.R            # apply_fire()
│   ├── recruitment.R     # bevholt(), ricker(), canopy_density(),
│   │                     #   shade_suppression(), sample_parents(), recruit()
│   ├── census.R          # census_year(), pct_decline()
│   ├── simulate.R        # run_simulation() -- the main annual loop
│   └── matrix_model.R    # Leslie/Lefkovitch matrix builders and
│                         #   projection tools, used by the test suite
│                         #   to validate the IBM against known-answer
│                         #   matrix models (see TESTING.md)
├── diagnostics/
│   ├── diag_rust_pressure.R          # pure functions: genotype/score
│   │                                 #   tables, mortality-vs-pressure,
│   │                                 #   dose-response shape comparison,
│   │                                 #   fire×rust compounding metric
│   ├── diag_age_structure.R          # pure functions: age-structure
│   │                                 #   table/summary, single- and
│   │                                 #   multi-scenario stacked plots
│   └── run_rust_pressure_diagnostics.R  # runnable script producing
│                                         #   the figures/numbers below
├── tests/
│   ├── testthat.R
│   └── testthat/          # see TESTING.md
├── params.R               # single source of truth for all parameters
└── main.R                 # example two-phase workflow (calibration →
                            #   rust projection + counterfactual)
```

## Core design decisions

### Data structure

- **`pop`**: a `data.frame`, one row per individual, alive or dead.
  Individuals are **never removed** — `alive` is flipped to `FALSE` — so
  the full pedigree and selection history stays queryable for as long as
  the run lasts.
- **`resist_gt`**: a parallel integer matrix, `n_individuals × n_loci`,
  allele doses `0`/`1`/`2`. **Row `i` of `pop` always corresponds to row
  `i` of `resist_gt`.** Both grow in lockstep on every `rbind()` of new
  recruits. Genotypes live outside the main table so that changing
  `n_loci` never changes the shape of the demography table.

### Individual columns

`id, age, alive, birth_year, death_year, mother_id, father_id,
age_first_flower, resprout, resprout_yrs_remain, flowering, resist_score`

### Life-history states (derived, never stored)

| State | Condition |
|---|---|
| Juvenile | alive, `age < age_first_flower`, `!resprout` |
| Reproductive adult | alive, `age >= age_first_flower`, `!resprout`, `flowering` |
| Resprout | alive, `resprout` |
| Dead | `!alive` |

**IUCN N** (the Criterion A count) is `sum(pop$alive & pop$flowering)` —
reproductive adults only, never total abundance.

### Genetic architecture (generalised to arbitrary `n_loci`)

Three equal-length vectors in `params.R` define the architecture:
`resist_locus_effect`, `resist_dominance`, `resist_freq0`. `resist_score`
is the sum of each locus's contribution (`effect × dominance_modifier`),
clamped at 1. This single clamp is what lets the same code represent:

- a single fully dominant locus (`resist_locus_effect = c(1.0)`,
  `resist_dominance = c(1.0)`) — heterozygote and homozygote both score 1
- a single fully recessive locus (`resist_dominance = c(0.0)`) — only the
  homozygote scores 1
- an oligogenic, additive, partial-dominance architecture where no single
  locus is sufficient on its own (e.g. `resist_locus_effect = c(0.4, 0.35,
  0.25)`, `resist_dominance = c(0.5, 0.5, 0.5)`)

Mendelian segregation in `inherit_alleles()` is ported from
jasongbragg/PlantPopGenFit (epi branch); loci are independent
(Hardy–Weinberg within locus, linkage equilibrium across loci) unless
future data justifies otherwise.

### Myrtle rust: a single tunable knob

`rust_pressure` (0 = none, 1 = reference humid-coastal site, >1 = more
severe) is the *one* cross-stage parameter representing geographic
variation in disease pressure. **How** that knob translates into a
realised effect is configured separately per life-history stage in
`rust_dose_response` (juvenile mortality, resprout mortality, resprout
delay), each with its own functional form — `"linear"`, `"power"`,
`"saturating"`, `"sigmoid"`, or `"threshold"` (see the docstring above
`dose_response()` in `mortality.R` for the exact definitions, and
`diagnostics/diag_rust_pressure.R::compare_dose_response_shapes()` to
visualise candidates before committing to one).

Rust is "switched on" from `rust_start_year` — before that, all three
effects are exactly zero, which is what makes the resistance locus
selectively neutral (freely drifting) in any pre-rust calibration run.

The rust-driven excess hazard (and resprout delay) for an individual is
always scaled by `(1 - resist_score)`: a fully resistant individual
experiences zero rust-specific excess hazard or delay, but is otherwise
demographically identical to a susceptible one.

### Fire model

Two ways to schedule fire:

- **`fire_years`**: a fixed list of calendar/sim years — used for
  calibration runs and any scenario where the historical fire record
  matters.
- **`fire_prob_annual > 0`**: a Bernoulli draw each year — used for
  ensemble runs (and overrides `fire_years` when set).

`apply_fire()` runs two vectorised passes: a kill pass (`rbinom` on
`fire_kill_prob`), then a resprout pass over survivors. Resprout recovery
has two modes (`params$resprout_recovery`):

- **`"countdown"`** (default, biologically realistic): duration drawn as
  `rpois(resprout_yrs_base) + rpois((1 - resist_score) * delay_extra)` at
  fire time, decremented each year.
- **`"geometric"`**: a constant per-year recovery probability
  (memoryless) — required for the Lefkovitch matrix-equivalence test (see
  TESTING.md); not intended for ordinary scenario runs.

### Canopy shading suppression of recruitment

Mature canopy suppresses how many seedlings successfully establish —
modelled as a discount on `recruit()`'s expected recruit count, not as
extra juvenile mortality, on the reasoning that establishment failure
and early seedling death aren't meaningfully distinguishable in the
field: both are just "a seed that never became a standing individual."

Two distinct counts drive two distinct processes, deliberately kept
separate:

- **`n_flowering`** (alive & flowering) drives `bevholt()`'s seed
  *production* — unchanged from before this mechanism existed.
- **`canopy_density(pop)`** (alive & `age >= age_first_flower`,
  regardless of current flowering/resprout status) drives shading
  *suppression*, via `shade_suppression()` — which reuses `dose_response()`
  (the same generic function `rust_modifiers()` uses) fed a raw density
  count instead of a normalised pressure.

```r
expected_rec <- bevholt(n_flowering, params) * (1 - shade_suppression(canopy_density(pop), params))
n_rec        <- rpois(1, expected_rec)
```

`max_effect = 0` by default — the mechanism exists but is a no-op unless
`params$shade_dose_response` is explicitly parameterised, so it doesn't
change any scenario that doesn't ask for it.

**Why `canopy_density()` checks age, not resprout status**: a fire that
actually kills an adult opens the canopy; a fire that merely sends a
*juvenile* into `resprout` doesn't, because that juvenile had no canopy
to lose in the first place. But `apply_fire()` sends survivors of *every*
stage into `resprout` — juveniles included (see its docstring). An
earlier version of this mechanism counted `resprout | age >= AFR` as
canopy, which let surviving juveniles re-labelled as resprouts silently
offset the canopy genuinely lost to adult mortality — verified against
an actual fire event before being caught (see TESTING.md for the
numbers). `canopy_density()` counts `age >= age_first_flower` alone,
so a resprouting *adult* still counts (real residual structure) but a
resprouting *juvenile* never does.

**Calibration note**: `shade_dose_response$half_sat` is in the same raw
count units as `canopy_density()` — not a normalised 0–1 pressure like
`rust_pressure`. It needs to be set relative to the actual canopy density
your population reaches; the placeholder default (`5000`) won't produce
a meaningful contrast for a population that only ever reaches a few
thousand.

### Juvenile-decline mortality (bathtub-shaped, canopy-dependent height)

A single Weibull hazard is monotonic in age for any fixed `weibull_k` —
it cannot represent "high mortality right after establishment, declining
over a couple of years to a background floor," which real survivorship
data for this kind of species generally shows. `juv_decline_hazard()`
adds a second, independent additive hazard term for true juveniles
(age `< age_first_flower`, not resprouting) on top of the existing
`weibull_hazard()` + `dd_hazard()` background:

```r
height(N_canopy) = dose_response(N_canopy, params$juv_decline_dose_response)
hazard(age)       = height * exp(-age / params$juv_decline_tau)
```

This is deliberately a *different* question from `shade_dose_response`
above: that asks "does a seed become a standing seedling at all" (wired
into `recruit()`); this asks "given a standing seedling, how much extra
risk does it carry in its first couple of years" (wired into
`nominate_deaths()`). Both can be active at once — they aren't redundant.

**Only the peak height scales with canopy density, not the duration.**
`juv_decline_tau` is fixed regardless of `N_canopy` — the vulnerable
*window* doesn't change, just how much is lost within it. At `N_canopy =
0` (a fresh post-fire site), `dose_response()` returns exactly `0`
regardless of `max_effect`, so the term vanishes entirely and a
post-fire seedling flush faces only the ordinary background hazard.
There is **no mortality floor at zero shade** — an earlier version of
this discussion considered one (e.g. "50% even with no shade") and
explicitly decided against it; if you want one later, `height()` would
need to become `base_height + dose_response(...)` instead of
`dose_response(...)` alone.

**No additional senescence is layered on top by default.** If you want
old-age removal to come from fire recurrence rather than a deterministic
aging hazard, set `weibull_k` close to `1` (near-constant background,
exactly the Level-1 matrix-equivalence case in TESTING.md) rather than
the default `2`. If instead you want a deterministic senescence term
but find Weibull's `k`/`lambda` awkward to calibrate, see
`senescence_dose_response` below — a more interpretable, bounded
alternative.

### Senescence mortality (Hill function of age, replacing Weibull)

Weibull (`weibull_k > 1`) is unbounded and its shape is set jointly by
two parameters (`k`, `lambda`) with no single number that directly
answers "at what age does senescence become noticeable." `senescence_hazard()`
offers a more interpretable alternative: a *rising* Hill function of age,
reusing `dose_response()`'s `"sigmoid"` form directly with age as the
pressure axis:

```r
hazard(age) = max_effect * age^hill / (age^hill + half_sat^hill)
```

`hazard(0) = 0`. `hazard(half_sat)` is **exactly** `max_effect / 2` — so
`half_sat = 30` means "by age 30, senescence-driven mortality is halfway
to its eventual ceiling," a single clearly-named number instead of a
`k`/`lambda` pair. `hazard(age)` asymptotes to `max_effect` as age grows,
a hard ceiling unlike Weibull, which keeps climbing forever for any
`k > 1`. `hill` controls steepness of the transition — larger `hill`
keeps the hazard negligible for longer below `half_sat`, then rises more
sharply around it.

**This is meant to *replace* `weibull_k > 1`, not stack with it.** To
switch a scenario over: set `weibull_k = 1` *and*
`senescence_dose_response$max_effect > 0`. **Be aware of what `weibull_k
= 1` actually leaves behind, though**: `weibull_hazard(age) =
(k/lambda)*(age/lambda)^(k-1)` reduces to a *constant* `1/weibull_lambda`
at every age when `k = 1` — it does **not** go to zero. At the default
`weibull_lambda = 30`, that's a flat ~3.3% annual hazard for every
individual at every age, stacking with `senescence_hazard()` via the
union formula. This is sometimes exactly what you want (a deliberate
non-senescent background rate — predation, random misfortune — *plus*
age-driven senescence on top), but if you want `senescence_hazard()` to
be the *only* source of age-related mortality, also raise
`weibull_lambda` to something very large (e.g. `1e9`, the convention
already used throughout `helper-fixtures.R`/`test-mortality.R` for "this
term should contribute essentially nothing") rather than relying on
`k = 1` alone to zero it out.

`max_effect = 0` by default: off unless explicitly parameterised, so
the matrix-equivalence validation (which builds its survival schedule
from `weibull_hazard()` only, via `build_leslie()` in `matrix_model.R`)
is completely unaffected by this mechanism's existence.

**Calibration note, and an honest result**: with `weibull_k = 1` and
`juv_decline` calibrated so mature-canopy first-year mortality reaches
~99% (vs. near-zero at `N_canopy = 0`), juvenile fraction at a fire-free
equilibrium dropped from 67% to 55% with the *default* `dd_alpha = 0.4`
active — a real, moderate effect, but **total population and absolute
adult count both shrink substantially too** (population more than
halves) — this is not a free improvement. An earlier check of this same
calibration with `dd_alpha` accidentally left at `0` (no density
regulation at all) showed adults nearly doubling instead; that result
doesn't hold once density dependence is restored to its normal strength,
and is not representative of the mechanism's real effect — see
TESTING.md if extending this further.

### Two-phase simulation workflow (see `main.R`)

```r
# Phase 1: pre-rust calibration (fire only, rust off) -- run long enough
# to reach a quasi-equilibrium age structure
res_calib <- run_simulation(params_calib, year0 = 1710, verbose = TRUE)
init_state <- list(individuals = res_calib$individuals, resist_gt = res_calib$resist_gt)

# Phase 2a: rust-impacted projection, from the SAME starting point
res_rust  <- run_simulation(params_rust,  init_state = init_state, year0 = 2010)

# Phase 2b: counterfactual -- identical starting point and fire regime,
# rust never arrives. Isolates the rust-attributable decline.
res_cfact <- run_simulation(params_cfact, init_state = init_state, year0 = 2010)
```

`run_simulation()` always calls `set.seed(params$seed)` first, so 2a and
2b consume *identical* random draws right up until `rust_modifiers()`
starts returning a nonzero effect — a paired (common-random-numbers)
comparison, not two independent stochastic runs. This property is locked
in by a dedicated regression test (see TESTING.md).

## Parameters reference (`params.R`)

| Group | Key parameters |
|---|---|
| Run control | `N0`, `n_years`, `seed` |
| Age-dependent mortality | `weibull_k` (>1 senescence, =1 constant, <1 front-loaded juvenile risk), `weibull_lambda`; `senescence_dose_response` (Hill-function alternative to `weibull_k>1`, bounded, `half_sat`=age at half-max severity); `juv_decline_dose_response`/`juv_decline_age_half_sat`/`juv_decline_age_hill` (bathtub-shaped juvenile mortality, canopy-dependent height) |
| Density dependence | `K`, `dd_alpha` |
| Recruitment | `R_max`, `K_half` (Beverton–Holt, driven by flowering-adult count), `shade_dose_response` (canopy suppression of establishment, driven by `canopy_density()`) |
| Flowering | `age_first_flower_mean`, `age_first_flower_sd` |
| Mating system | `selfing_rate` (0 = obligate outcrossing, 1 = no selfing penalty) |
| Fire | `fire_years`, `fire_prob_annual`, `fire_kill_prob`, `resprout_yrs_base`, `resprout_recovery`, `resprout_prob_recovery` |
| Myrtle rust | `rust_start_year`, `rust_pressure`, `rust_dose_response` (per-stage: `juv`, `resprout`, `delay`) |
| Genetics | `resist_locus_effect`, `resist_dominance`, `resist_freq0` |

Nothing is hard-coded into the simulation functions — every scenario is
built by copying `get_default_params()` and overriding individual
elements.

## Running a scenario

```bash
Rscript main.R
```

produces the two-phase calibration → rust/counterfactual comparison
described above, printing the projected percent decline in reproductive
adults for both scenarios and (if `ggplot2` is available) plotting
abundance and resistance-allele-frequency trajectories. Falls back to
base graphics if `ggplot2` isn't installed.

## Diagnostics: biological plausibility, not validation

```bash
Rscript diagnostics/run_rust_pressure_diagnostics.R
```

produces, under `outputs/`:

- **`rust_pressure_mortality.png`** — annual mortality probability vs
  `rust_pressure`, by resistance class (susceptible / heterozygous /
  resistant / population mean), for juvenile and resprout stages.
- **`dose_response_shapes.png`** — all five functional forms normalised
  to `max_effect = 1`, for choosing per-stage shapes.
- **`cumulative_resprout_risk.png`** — P(dies before resuming flowering |
  survived the fire), the fire × rust compounding metric. This is
  explicitly an **approximation** (treats resprout duration as fixed at
  its expected value and the per-year hazard as constant across that
  window) — the simulation engine itself draws both stochastically every
  fire event; the diagnostic exists to make the interaction legible to a
  reader, not to substitute for the IBM.

These diagnostics answer "does the parameterisation imply something
biologically sensible?" — a different question from "is the code
correct?", which is what the test suite (TESTING.md) is for.

### Age structure at the end of a run

`diagnostics/diag_age_structure.R` provides standalone functions (no
dependency on the rest of the model — they only read `pop`'s own
columns) for inspecting the age structure of any population snapshot,
typically `res$individuals` at the end of a `run_simulation()` call:

```r
source("diagnostics/diag_age_structure.R")

summarize_age_structure(res$individuals)        # numeric summary by stage
age_structure_table(res$individuals)            # long data frame, for custom plotting
plot_age_structure(res$individuals)             # stacked bar, by life stage

# Compare two scenarios at the same final year side by side, e.g. a rust
# projection against its no-rust counterfactual:
plot_age_structure_compare(list(
  rust                   = res_rust$individuals,
  no_rust_counterfactual = res_cfact$individuals
))
```

`plot_age_structure(pop, age_breaks = seq(0, 70, by = 5))` bins ages
instead of using single-year bars — useful once a population spans
several decades and single-year bars become unreadable.

**`plot_age_structure_compare()` defaults to a shared y-axis across
panels (`shared_y = TRUE`).** Earlier versions independently rescaled
each panel — which meant the tallest bar in *every* panel always filled
its own axis, regardless of the actual counts, so a genuinely large
difference in juvenile dominance between two scenarios could be
invisible by eye (both panels just "looked like the same shape," each
against its own ruler). Caught this exact way: a real ~30 percentage
point drop in juvenile fraction between two calibrations looked
unchanged in the chart until the axes were forced to match (see
TESTING.md for the side-by-side). Set `shared_y = FALSE` only if
absolute scale genuinely isn't the point of the comparison and you want
each panel's internal shape judged purely on its own terms.

**Pull the number, don't eyeball the chart, when it matters.**
`summarize_age_structure(pop)$by_stage$prop` gives the juvenile/adult
fraction directly — immune to any plotting-scale issue, shared-axis or
not. Treat the chart as illustration, the summary as the source of truth.

## What the IBM can do that a matrix model cannot

The test suite validates the IBM against matrix population models in
regimes simple enough for both to apply (see TESTING.md). Genuinely
individual-based behaviour — and therefore *not* coverable by that
validation — includes:

- individual genetic resistance variation and its selective dynamics
- density dependence in mortality or recruitment
- individual variation in age-at-first-flower
- resistance-scaled resprout delay (the fire × rust interaction)

These are the domain of the diagnostics above, and of the full
`run_simulation()` integration tests.

## Extensibility notes

- **More loci / different architecture**: edit `resist_locus_effect`,
  `resist_dominance`, `resist_freq0` in `params.R` — nothing else needs
  to change, including the test suite's genetics tests, which are
  written generically over `n_loci`.
- **Alternative recruitment function**: `ricker()` is already provided
  in `recruitment.R` as a drop-in alternative to `bevholt()` (needs
  `params$ricker_r`); swapping it changes the recruitment count's
  response to crowding but not the rest of the pipeline.
- **Geographic ensemble** (planned next step): `rust_pressure` and
  `fire_prob_annual` drawn from geographic distributions across many
  thousands of replicate runs, using the stochastic `fire_prob_annual`
  fire mode already implemented.
