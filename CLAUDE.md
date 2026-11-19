# CLAUDE.md — MnSim multispecies branch

<!-- Place at the repo root. Claude Code reads this at the start of every session. -->

## What this project is
MnSim is an individual-based model (IBM), written in R, of a Myrtaceae population under myrtle rust
(*Austropuccinia psidii*), fire and demographic processes. It is used for IUCN Criterion A
assessments. It was built for *Melaleuca nodosa*. This branch generalises it to other
rust-affected species: *Rhodamnia rubescens*, *Rhodomyrtus psidioides*, *Archirhodomyrtus beckleri*,
*Decaspermum humile* and *Backhousia leptopetala*.

The work plan, with checklists, is in `PLAN.md`. Work on one phase at a time and tick items off as they are done.

## Repo layout
- `plant_ibm/` — the IBM code
- `species/` — one YAML parameter file per species (created in Phase 2/4): parameters, rust anchor,
  assessment windows, range box, range-wide draws; `<sp>_validation.yaml` holds its validation set-up
- `config/` — `schema.yaml` (every parameter and setting: layer, type, bounds), `landscape.yaml`
  (shared landscape defaults and spatial settings), `assessment.yaml` (shared IUCN settings),
  `run_defaults.yaml` (default species, run settings); loaded by `plant_ibm/R/species_config.R`
- `plant_ibm/spatial/data/<sp>/` — per-species occurrence and site tables, fire checks
- `plant_ibm/spatial/data/fire/` — fire cell-year tables shared by all species (`cells.csv` lists the
  cells covered); cells are numbered on one master grid (`plant_ibm/spatial/cell_ids.R`)
- `data/rust_species_params.csv` — literature parameters and calibration targets, one row per value,
  each with a source. Columns: species, ibm_module, parameter, role, value, units, evidence_type,
  confidence, source, notes. `role` is one of: parameter | calibration_target | committee_method |
  committee_projection.
- `data/m_nodosa_params.csv` — the same columns for *M. nodosa*, whose values come mostly from this
  project's own calibration and assumptions (sources cite `docs/methodological_notes.md` entries).
  Extra `evidence_type` labels: calibration | model_assumption | model_choice.
- `tests/` — testthat tests; `tests/reference/` holds frozen M. nodosa outputs
- `docs/methodological_notes.md` — running log of assumptions and decisions
- `private/` — unpublished data (git-ignored)

## Hard rules
1. **Never read, copy, summarise or commit anything in `private/`.** It holds unpublished resistance
   data. Where code needs resistance data, use placeholder values from the species YAML.
2. **Never edit files in `tests/reference/`.** If a regression test fails, the code is wrong, not
   the reference. Stop and report the difference. Don't "fix" it by changing expected values.
3. **Run the tests before every commit:** `cd plant_ibm && /opt/R/4.3.1/bin/Rscript tests/testthat.R`.
   Do not commit while tests are failing. (Run from `plant_ibm/`: the helper sources model code by
   relative path. Use R 4.3.1: the system `Rscript` is R 3.4.4 and cannot load `sf`.)
4. **Do not hand-tune parameters to hit calibration targets.** All fitting goes through the
   calibration code (Phase 5), so it is logged and reproducible.
5. **Every parameter value in a species YAML must cite its row** (species + parameter) in
   `data/rust_species_params.csv` or `data/m_nodosa_params.csv`, or be marked `prior: true` with the
   range and a reason. No unsourced numbers.
6. **Log every modelling assumption or structural decision** as a dated entry in
   `docs/methodological_notes.md`.
7. **Propose before building new model structure.** For new modules (clonal stage, hazard
   accumulation, rust-pressure modes), write a short design note in `docs/` and wait for approval
   before coding.

## Conventions
- R, base pipe `|>`; data.table where the existing code already uses it; match the existing style.
- Fixed seeds in all tests (`set.seed()` at the start of each test).
- Keep test runs small (population ≤ 500, ≤ 20 replicates) so the suite runs in a few minutes on a laptop.
- Large ensembles are not run on the laptop. Write them as scripts with a `--reps` argument,
  to be run on the workstation (rcn01).

## IUCN framing (do not change without asking)
- M. nodosa uses Criterion A4 with one fixed window per run, and window length taken from the
  generation time inferred in that run.
- Reduction windows for A3/A4 are capped at 100 years (IUCN guidelines).
- The main output is the frequency distribution of runs across IUCN categories, not a single central value.
- "Mature individuals" count living plants of reproductive age, including those whose flowering is
  suppressed. Keep the flowering-only counting code as an alternative.

## Rust pressure
Rust pressure for M. nodosa comes from the `mrsdm` SDM, normalised to the value at Yarringully NR.
Keep that anchor. New species may use a different mode (see PLAN.md Phase 3).
