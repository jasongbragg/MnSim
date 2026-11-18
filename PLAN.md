# PLAN.md — Generalising MnSim to multiple rust-affected species

<!-- Place at the repo root. Tick boxes as work completes. Who does each item: [agent], [you], [both]. -->

Species order: *R. rubescens* first (best calibration data), then *R. psidioides*, then
*A. beckleri* and *D. humile*. *B. leptopetala* is optional, as a low-impact contrast species with
parameters set rather than fitted.

---

## Phase 0 — Setup  [you]
- [ ] Install Claude Code; check with `claude --version`
- [ ] Create the `multispecies` branch; make `data/ private/ species/ tests/reference/ docs/`
- [ ] Add `private/` to `.gitignore`; confirm no unpublished resistance data is tracked
- [ ] Copy `rust_species_params.csv` into `data/`; add `CLAUDE.md`, `PLAN.md`, `.claude/settings.json`
- [ ] Install R packages (testthat, yaml, abc, data.table, future.apply, renv); `renv::init()`
- [ ] Time one M. nodosa run and record wall time and peak RAM here: ________

## Phase 1 — Regression harness  [agent writes, you check]
Goal: a safety net, so any refactor that changes M. nodosa's behaviour is caught.
- [ ] [agent] Read `plant_ibm/` and write `docs/architecture.md`: the modules, the main loop, where each parameter is used
- [ ] [you] Check the architecture note is correct
- [ ] [agent] Choose 3–4 small fixed-seed scenarios (no fire, fire, high rust, low rust); save outputs to `tests/reference/`
- [ ] [you] Make the reference files read-only (`chmod -R a-w tests/reference`)
- [ ] [agent] Write testthat tests that re-run each scenario and compare (exact match for the same seed)
- [ ] [agent] Add a quick test of summary statistics (mean mortality over 20 seeds within tolerance) as a second check
- **Done when:** the test suite passes on unchanged code, runs in under 5 minutes, and is committed.

## Phase 2 — Refactor to species configs  [agent]
Goal: no M. nodosa constants left in the code; everything is in `species/m_nodosa.yaml`.
- [ ] [agent] List every hard-coded biological constant (with file:line) in `docs/constants_inventory.md`
- [ ] [you] Review the inventory: which are biology (go into YAML) and which are numerical settings (stay in code)
- [ ] [agent] Write a YAML schema and loader with validation (required fields, types, ranges)
- [ ] [agent] Move the constants into `species/m_nodosa.yaml`; the model reads them from config
- [ ] [agent] Make swappable interfaces for: fire response, rust pressure, reproduction/recruitment
- **Done when:** the Phase 1 tests pass with identical outputs, and grepping for the old constants finds nothing.

## Phase 3 — New model structure  [agent proposes, you decide]
For each item: design note in `docs/` → [you] approve → implement → add tests.
- [ ] **Clonal sucker stage** (*R. psidioides*): canopy death leaves a persisting root-sucker pool with no seed input; suckers get infected at about 30 cm
- [ ] **Rust-pressure modes:** (a) SDM-scaled (current); (b) saturating/constant. Evidence: no climate predictor of infection for *R. psidioides* (Fensham et al. 2020); strong NSW vs SE Qld contrast for *A. beckleri* (Manea et al. 2019 vs Pegg et al. 2017)
- [ ] **Accumulating-damage hazard:** death risk rises with consecutive defoliation events (Olney exclusion trial: ~6-month lag, death in <4 yr; Table 10 resurveys)
- [ ] **Fruit-infection suppression:** fecundity → ~0 under high rust (Tucki Tucki NR; no fruit seen for *R. psidioides* since 2010)
- [ ] Check that M. nodosa outputs are unchanged when each new option is switched off (Phase 1 tests still pass)

## Phase 4 — Species parameter files  [agent, you set priors]
- [ ] [agent] Write `scripts/csv_to_yaml.R`: turns `role == parameter` rows into `species/<sp>.yaml`, each value carrying its CSV row reference
- [ ] [agent] List every gap (value NA, or low confidence) per species in `docs/parameter_gaps.md`
- [ ] [you] Set prior ranges for the gaps, at minimum generation length, age at first reproduction and longevity, with a reason for each
- [ ] [agent] Build `species/r_rubescens.yaml`, then `species/r_psidioides.yaml`; the loader validates both
- [ ] [agent] Smoke test: each species runs 50 years without errors and gives plausible output

## Phase 5 — Calibration  [agent builds, you judge]
Targets are the `role == calibration_target` rows in the CSV.
- [ ] [agent] Write `scripts/targets.R`: turns each species' targets into summary statistics with tolerances
  - *R. rubescens:* 12% dead at ~3.5 yr after arrival (n=669, 43 sites) and the spread across sites; Table 10 paired resurveys 2014→2016; Olney stand 53% dead in <4 yr; mean crown transparency ~0.76
  - *R. psidioides:* 57% dead at ~3.5 yr (n=297, 18 sites) and the spread across sites; by 2018, 23% of populations extinct and 61% reduced to suckers
  - *A. beckleri / D. humile:* Ryans Road 2016→2017 (13→44%, 36→73%); NSW populations with little mortality (Manea et al. 2019)
- [ ] [agent] Rejection ABC with the `abc` package; draw from the priors; distance on the summary statistics
- [ ] [agent] Benchmark: runs per hour on the laptop → estimate the cost of a full ABC; if it is more than about an overnight run, write it as a script for rcn01
- [ ] [you] Check the posterior predictive plots: does the fitted model reproduce the between-site variance as well as the mean?
- [ ] Decide whether rejection ABC is enough or whether ABC-SMC is needed

## Phase 6 — IUCN comparison and sensitivity  [both]
- [ ] [agent] Implement "committee mode": a constant exponential decline from the observed 3.5-yr mortality over 3 generations (capped at 100 yr). Check it reproduces NSW TSSC 96.3–98.8% for *R. rubescens*
- [ ] [agent] For each species, the distribution of IUCN categories under: (1) committee mode, (2) IBM with an A3 window, (3) IBM with an A4 window (fixed window anchored at rust arrival in 2010)
- [ ] [agent] Global sensitivity analysis (e.g. Sobol or PRCC) on the posterior plus prior-only parameters
- [ ] [you] Interpret the results: which unmeasured parameters most change the category? These set the fieldwork priorities.

## Phase 7 — Independent check  [fresh agent session]
- [ ] Start a new Claude Code session with no history. Prompt: *"Check every value in species/*.yaml against data/rust_species_params.csv and the cited source. Report mismatches, unsourced values and priors without a reason. Do not edit anything."*
- [ ] [you] Resolve the findings; update `docs/methodological_notes.md`

---

## Decisions log
| Date | Decision | Why |
|---|---|---|
| | | |

## Known data issues (from the literature review)
- *R. psidioides*: NSW TSSC 99% projection can't be reproduced from its stated rate; the same report also says ~10% decline in <5 yr, which conflicts with 57% mortality
- *R. rubescens*: 3 generations at 40 yr = 120 yr, above the IUCN 100-yr cap
- *D. humile* 2016 mortality: 36% (Pegg et al. 2017 PLoS ONE) vs 53.85% (PBCRC report Table 12)
- *D. humile*: the NE Qld form is rated relatively tolerant; confirm which form is being modelled
- Manea et al. 2019: a single October survey in drought with no active infection; treat as weak evidence
