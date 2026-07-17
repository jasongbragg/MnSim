# Testing Guide

This documents the `tests/testthat/` suite for the plant population IBM.
For the model itself, see **README.md**.

## Running the suite

```bash
Rscript tests/testthat.R
```

`testthat` isn't on CRAN-via-pip in this environment; on a fresh Ubuntu
box it installs via `apt-get install r-cran-testthat`. The project isn't
a formal R package (no `DESCRIPTION`), so there's no `test_check()` to
hand off to — `tests/testthat.R` calls `testthat::test_dir()` directly,
and `tests/testthat/helper-fixtures.R` (auto-loaded first) sources every
`R/*.R` file and `params.R` itself, then defines the fixtures shared
across test files.

Current suite: **149 tests across 9 files**, full run **~6.8 seconds**.

```
census: ....................
fire: .................
genetics: .................
individuals: ..............
matrix-equivalence: ..............................................
mortality: ............................................................................
params: .......................
recruitment: .................................
simulate: ..............
```

## Test categories

Every test in the suite falls into one of five categories:

1. **Deterministic edge cases** — inputs chosen so the expected output is
   exact, not statistical (e.g. `SS × SS` always produces `SS`;
   `fire_kill_prob = 1` kills everyone; `dd_hazard(K, params) == dd_alpha`
   exactly).
2. **Statistical checks** — stochastic functions tested at large `N`
   against their theoretical expectation, with a tolerance sized to the
   sampling noise actually present (see the fixture pitfall below — this
   bit us once already).
3. **Structural invariants** — properties that must hold regardless of
   parameters or randomness: `pop` and `resist_gt` stay row-aligned,
   dead individuals always have a recorded `death_year`, alive
   individuals never have an `NA` `resist_score`.
4. **Regression / reproducibility** — identical seed and params give
   byte-identical results; the rust-vs-counterfactual common-random-
   numbers property (see below) is locked in by a dedicated test so a
   future refactor can't silently break it.
5. **Integration smoke tests** — `run_simulation()` runs end-to-end
   without error and produces internally consistent output across a
   multi-year run, including the extinction path.

| File | Tests | Covers |
|---|---|---|
| `test-genetics.R` | 13 | Mendelian segregation, `resist_score_from_gt` across architectures, allele frequencies |
| `test-mortality.R` | 52 | `weibull_hazard`, `hill_weight`, `dd_age_weight`, `senescence_hazard`, `dd_hazard`, `dose_response` (all 5 forms), `rust_modifiers`, `juv_decline_hazard`, disease-triangle integration, `nominate_deaths` |
| `test-fire.R` | 13 | `apply_fire` legacy flat-probability mode AND two-stage age-dependent mode, resist-scaled delay |
| `test-recruitment.R` | 24 | `bevholt`, `sample_parents` (selfing), `canopy_density`, `shade_suppression`, `recruit` end-to-end |
| `test-individuals.R` | 7 | `create_population`, `make_recruit_rows`, `env_susceptibility` column |
| `test-census.R` | 8 | `census_year` columns, `pct_decline` |
| `test-params.R` | 12 | `get_default_params()` schema sanity |
| `test-simulate.R` | 8 | `run_simulation` integration, reproducibility, extinction, common-random-numbers |
| `test-matrix-equivalence.R` | 12 | IBM-vs-matrix validation (see below) |

## Matrix-model equivalence: the validation centrepiece

The single most rigorous test available for this kind of model: contrive
an IBM parameterisation that **degenerates exactly to a known-answer
matrix population model**, run both, and check they agree. A mismatch
here is a genuine implementation error, not a calibration question —
there's no "tune the tolerance until it passes" available, because the
matrix answer is fixed by the chosen parameters, not fit to the IBM.

`R/matrix_model.R` provides three levels, each corresponding to
successively richer matrix structure:

### Level 1 — Geometric

Constant hazard (`weibull_k = 1`), `age_first_flower_mean = 1` (everyone
flowers after their first year), no density dependence, no fire. The
system degenerates to a single-class geometric growth model with a
**closed-form** answer: `λ = s·(1 + F)`, where `s` is annual survival and
`F` is per-adult fecundity. Tested two ways: IBM's realised asymptotic
growth rate against the closed form directly, and against the dominant
eigenvalue of `build_leslie()` (confirming the matrix builder's
degenerate case is itself correct before trusting it elsewhere).

### Level 2 — Leslie

Age-structured Weibull hazard (`weibull_k = 2`) with delayed flowering
(`age_first_flower_mean = 5`). No closed form here — the benchmark is the
dominant eigenvalue of the full age-structured Leslie matrix built by
`build_leslie()` directly from the same Weibull parameters. Tested at two
independent `(F, AFR)` combinations, including one with sub-replacement
fecundity, to rule out the agreement being an artefact of one specific
parameter choice.

### Level 3 — Lefkovitch

Adds the Resprout life-history stage, driven by fire. This is where the
validation got genuinely interesting — see "What went wrong" below.

### Comparison methodology: asymptotic λ vs. year-by-year trajectory

Levels 1–2 compare **asymptotic growth rate**: let the IBM run past a
burn-in period (so the initial age distribution's transient effects wash
out), estimate its realised year-on-year growth rate
(`estimate_lambda_ibm()`), and compare to the matrix's dominant
eigenvalue (`dominant_lambda()`).

Level 3 instead compares a **year-by-year trajectory** under a *fixed*
fire schedule (`build_lefkovitch_nofire()` / `build_lefkovitch_fire()` /
`project_lefkovitch_seq()`), starting both the IBM and the matrix from
the *exact same* initial stage vector. This wasn't the first design
tried — see below for why.

## What went wrong while building this (and why it's worth keeping)

### 1. Stochastic fire makes an asymptotic-λ comparison unreliable

The first attempt at Level 3 used `fire_prob_annual` (a Bernoulli draw
each year) and compared the IBM's asymptotic λ to the dominant eigenvalue
of the probability-weighted mean matrix `(1−p)·M_nofire + p·M_fire`.
Across seeds, the final-year IBM/matrix ratio ranged from **0.85 to
1.03** — far too noisy to be a meaningful pass/fail test.

This isn't a bug; it's a real result in stochastic demography. For a
matrix model with random year-to-year switching, the realised long-run
("stochastic") growth rate is **not in general equal to** the dominant
eigenvalue of the mean matrix — switching variance itself depresses
realised growth relative to that deterministic benchmark (Tuljapurkar's
result). An asymptotic-λ comparison under random fire was conflating "did
I implement the per-year mechanics correctly" with "did this particular
draw of fire timing land close to its theoretical long-run average" —
two different questions, and the test only needed to ask the first one.

**Fix**: switch to a *fixed* `fire_years` schedule and project the matrix
through the identical year-by-year sequence (`project_lefkovitch_seq()`).
Both trajectories then experience exactly the same fire events in the
same years, so any remaining mismatch is genuine demographic sampling
noise in the IBM, not switching-process noise. Final-year ratios across 7
seeds then settled to **[0.965, 1.033]** — comfortably inside the 0.08
tolerance the test actually uses.

### 2. A real timing bug, caught by the validation working as intended

`build_lefkovitch_nofire()` initially set the Juvenile column's fecundity
entry to `0`, reasoning that juveniles don't reproduce. But with
`age_first_flower_mean = 1`, an age-0 juvenile **ages to 1 within the
same simulation year** — which already meets the flowering threshold — so
it flowers and reproduces exactly like an existing adult. The IBM was
correct; the matrix builder had the bug. The Level-3 trajectory comparison
caught this immediately as a large, systematic mismatch (the matrix
column gave roughly half the expected offspring contribution). Fixed by
setting `M_nf[J, J] = F·s` whenever `AFR <= 1`, with the timing
subtlety documented directly in the function's docstring and covered by
its own dedicated test (`"build_lefkovitch_nofire: with AFR=1, Juvenile
and Adult columns are identical"`).

### 3. `K_half` needs to be far above any population size the run reaches

A second systematic (not random) bias appeared once population size grew
into the tens of thousands: the IBM trailed the matrix by ~5–15%, growing
with population size. Beverton–Holt recruitment is `R = F·N / (1 + N/K_half)`
— the matrix model's fecundity term `F` assumes the **linear, low-density
limit** `R ≈ F·N`, which only holds while `N << K_half`. With
`K_half = 1e6` and `N` reaching several tens of thousands, `N/K_half` was
large enough (a few percent) to measurably suppress realised recruitment
below the matrix's assumption — not a bug, but a **confound built into
the test's own parameter choice**. Fixed by setting `K_half = 1e9` in all
three fixture levels, two-to-three orders of magnitude above any `N` the
runs reach, collapsing the bias to <1%.

### 4. The fixture itself was quietly weakening every statistical test

`make_toy_pop(n, params)` took an `n` argument but never actually set
`params$N0 <- n` before calling `create_population()` — every call
silently created `params$N0` (200) individuals regardless of what was
requested. A test asking for `make_toy_pop(5000, p)` to get a tight
sampling-noise margin was actually running on 200, with 5× the intended
standard error. This had been quietly true for *every* statistical test
in the suite; it only surfaced when one assertion's tolerance happened to
be tight enough to fail on the resulting noise. Fixed in
`helper-fixtures.R`; **if extending this fixture, check that any new `n`
or `N0`-like argument actually flows through to `params$N0`** before
trusting a statistical test's tolerance.

### 5. A real bug in the canopy-shading mechanism, caught before it ever ran

When `recruit()`'s canopy-shading suppression (`README.md`'s "Canopy
shading" section) was first wired up, `N_canopy` was defined as
`resprout | age >= age_first_flower` — alive adults *or* resprouting
survivors, on the reasoning that a resprouting tree still has residual
trunk structure even though it isn't currently flowering. Reasonable in
isolation, but `apply_fire()` sends survivors of *every* stage into
`resprout` — juveniles included (see its own docstring: "applies
uniformly to juveniles and adults"). So a fire that killed half the adult
canopy also pushed every surviving juvenile into `resprout`, and those
juveniles got counted as canopy under the original formula — silently
offsetting the very canopy loss the mechanism was built to detect.

This wasn't caught by reasoning about the code; it was caught by
instrumenting an actual fire event and looking at the numbers
(`N_canopy` going *up* across a fire year, 1561 → 1789, instead of down)
before trusting the mechanism. Fixed by dropping `resprout` from the
condition entirely — `canopy_density(pop)` is just `age >=
age_first_flower`, so a resprouting *adult* still counts (real residual
structure) but a resprouting *juvenile* never does, regardless of how it
got into `resprout`. `test-recruitment.R`'s `"canopy_density counts
alive flowering-age individuals, excluding resprouting juveniles"` test
encodes exactly this scenario as a deterministic six-row fixture, so a
future change can't silently reintroduce the conflation.

**General lesson, same shape as pitfall #2**: any time a new "is this
individual currently in stage X" flag is reused as a proxy for something
else (here: "does this individual have a canopy"), check what *other*
pathways can set that same flag — `grep` for every place the flag gets
assigned, not just the place you're adding logic to.

### 6. An isolating test setting masqueraded as a real result

While calibrating `juv_decline_hazard` to a target first-year mortality
rate, `dd_alpha` was deliberately set to `0` to isolate the new term's
effect from density-dependent mortality during a quick sanity check.
That check looked great — juvenile fraction fell from 67% to 38%, adults
*nearly doubled*. The "doubling" was real for that run, but it wasn't a
property of the mechanism; it was a property of having no density
regulation at all, which lets removed-juvenile "room" translate directly
into uncapped adult growth (population growing toward whatever `bevholt()`
and weak hazards allow, unconstrained by `dd_hazard`).

Re-running the identical calibration with `dd_alpha` restored to its
default (`0.4`) and no fire confound gave a much more modest, and more
honest, result: juvenile fraction 67% → 55%, but **total population and
absolute adult count both shrink substantially** rather than the adult
population growing. Same mechanism, same parameters, opposite story on
adults — entirely because of one isolating setting left in from an
earlier, narrower check.

**General lesson**: a parameter set to an extreme value "just to isolate
one effect" is doing real work in the background — re-run the
final read with every other parameter back at its intended value before
trusting a result, especially one that sounds like a clean win.

### 7. A real plotting bug, caught because it made a real improvement invisible

`plot_age_structure_compare()` originally used
`facet_wrap(~ scenario, scales = "free_y")` (ggplot2) and independent
per-panel auto-scaling (base R fallback) — each panel's y-axis
rescaled to fit its own tallest bar. This is a genuine bug for the
function's main use case: comparing relative juvenile dominance between
two scenarios. With independent rescaling, age 0 always looks like it's
hitting the top of the chart, in *every* panel, regardless of the
actual counts — because each panel's tallest bar always fills its own
axis. Two scenarios with a real, substantial difference in juvenile
fraction can render as "the same shape," because each is being judged
against a different ruler.

This is exactly what happened in practice: a calibration that genuinely
dropped juvenile fraction from 67% to 46% (a real, verified number from
`summarize_age_structure()`) looked completely unchanged in the
side-by-side chart, because both panels independently rescaled to their
own tallest bar. The fix shares one y-axis across all panels by default
(`shared_y = TRUE` in `plot_age_structure_compare()`), computed as the
tallest *stacked* bar across every panel, not just one. `shared_y =
FALSE` is preserved as an explicit opt-out for the (rarer) case where
absolute scale genuinely isn't the point of the comparison.

**General lesson**: a chart that auto-scales per-panel can hide the
exact effect you're trying to detect, especially when "did this go up
or down" is itself the question. When a visual comparison is the
thing a user will actually look at and trust, treat its scaling choice
as part of the correctness of the function, not a cosmetic detail —
and always have a non-visual number (here, `summarize_age_structure()`'s
`prop` column) to check the chart against.

### 8. "Turn off Weibull" advice that didn't actually turn it off

When `senescence_hazard()` (a Hill-function alternative to Weibull
senescence) was introduced, the documented recommendation for switching
a scenario over was "set `weibull_k = 1`." This is true in the narrow
sense that `k = 1` removes the age-*dependent* shape -- but
`weibull_hazard(age) = (k/lambda)*(age/lambda)^(k-1)` reduces to a
*constant* `1/lambda` at every age when `k = 1`, not zero. At the
project's own default `weibull_lambda = 30`, that's a flat ~3.3% annual
hazard, stacking with the new senescence term via the union formula --
exactly what someone following the documented advice would NOT expect
if they read "turn off Weibull" as "make this term contribute nothing."

The actual test code for the new mechanism (written in the same session)
never had this bug -- every integration test isolated the new term with
`weibull_lambda = 1e9`, not `weibull_k = 1` alone, which is the correct
way to make a Weibull-based hazard genuinely negligible. The bug was
purely in the prose recommendation (README.md, params.R comments, and
the chat explanation), caught only because the user asked "does
`weibull_k = 1` actually turn Weibull off... not so sure" rather than
taking the documentation at face value.

**General lesson**: "setting parameter X removes effect Y" is a claim
that deserves the same scrutiny as a code change, especially for any
parameter whose value at one extreme (here, `k=1`) might *look* like
"off" without actually evaluating to zero. When writing this kind of
recommendation, compute the actual resulting value at the documented
defaults before asserting what it does -- and prefer phrasing that
states the literal resulting formula ("reduces to a constant `1/lambda`")
over a shorthand label ("flat background") that can be misread as
negligible.

## Validated parameter sets (for reference — see `helper-fixtures.R`)

| Level | `weibull_k` | `weibull_lambda` | AFR | `F` | `N0` | `n_years` | Tolerance | Result |
|---|---|---|---|---|---|---|---|---|
| 1 Geometric | 1.0 | 10 | 1 | 0.15 | 1,000 | 80 | 0.015 | diff ≈ 0.0009 |
| 2 Leslie | 2.0 | 25 | 5 | 0.10 | 1,000 | 80 | 0.020 | diff ≈ 0.0015 |
| 3 Lefkovitch | 1.0 | 10 | 1 | 0.18 | 20,000 | 40 | 0.08 (ratio) | ratios ∈ [0.965, 1.033] across 7 seeds |

All three use `K_half = 1e9` (density independence in practice) and a
fully-resistant population (`resist_freq0 = c(1.0)`) so genetic variation
doesn't introduce a mortality-rate confound unrelated to what's being
tested.

## What matrix equivalence deliberately does *not* cover

Matrix population models cannot represent individual-level heterogeneity.
The following are therefore tested elsewhere (unit tests on the relevant
function, or the diagnostics in README.md), never by equivalence:

- individual genetic resistance variation (`test-genetics.R`)
- density-dependent mortality or recruitment
- individual variation in age-at-first-flower
- resistance-scaled resprout delay — the fire × rust interaction
  (`test-fire.R`'s `"resist_score proportionally reduces..."` test)
- canopy-density-dependent recruitment suppression (`test-recruitment.R`'s
  `canopy_density`/`shade_suppression`/`recruit` tests)
- canopy-density-dependent juvenile mortality height (`test-mortality.R`'s
  `juv_decline_hazard`/`nominate_deaths` integration tests)
- Hill-function senescence (`senescence_hazard`) — off by default
  (`max_effect = 0`); `build_leslie()` builds its survival schedule from
  `weibull_hazard()` only, so this mechanism's existence cannot affect
  the Level-2 equivalence regardless of how it's parameterised elsewhere
- the disease-triangle rust formulation (`effective_susceptibility =
  host × environment × pathogen`) — gated by `rust_start_year = Inf` in
  every matrix fixture, same protection as the original rust mechanism
- the two-stage age-dependent fire model — gated by `fire_p_fimp = 0`
  (legacy mode) in every matrix fixture; `build_lefkovitch_fire()` only
  ever encodes the legacy flat-`fire_kill_prob` mechanics

## Adding new tests

- New statistical test? Check the sample size you're actually getting
  (see pitfall #4 above) before picking a tolerance. Also check the
  *magnitude* of what you're estimating — a large expected count (e.g. a
  Poisson mean in the hundreds) already has small relative noise per
  draw, and replicating far more than that noise warrants just slows the
  suite down for no statistical benefit (the original shading tests did
  exactly this: 2,000 replicates where 200 gave the same ~10σ safety
  margin at a tenth of the runtime).
- New matrix-equivalence level? If it involves any stochastic event
  schedule (fire, disturbance), prefer a fixed schedule + trajectory
  comparison over an asymptotic-λ comparison unless you've confirmed the
  random-switching growth rate genuinely converges to the mean-matrix
  eigenvalue for your specific setup (it often won't, exactly).
- New genetic architecture? `test-genetics.R`'s `resist_score` tests are
  already written generically over `n_loci` — extending
  `resist_locus_effect`/`resist_dominance` in a test's `make_toy_params()`
  call exercises the new architecture without new test code.
- New mechanism reusing an existing per-individual flag as a proxy for
  something else? See pitfall #5 — check every code path that can set
  that flag, not just the one you're adding logic around.
- New isolating parameter set ("just to check X in isolation")? See
  pitfall #6 — re-run with every other parameter restored to its
  intended value before trusting the result, especially a clean-looking
  one.
- New comparison plot? See pitfall #7 — if the chart's job is to show
  whether something changed, per-panel auto-scaling can hide exactly
  that. Default to a shared scale, and always have `summarize_*()`-style
  numbers to check the chart against.
