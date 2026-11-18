# Review: existing tests and parlib against PLAN.md Phases 1–2

Written 2026-11-19 on branch `multispecies` (commit 9d32066). It goes with
`docs/architecture.md`; section numbers below refer to that note.

## What already exists

- **Test suite.** `plant_ibm/tests/testthat/` has 16 test files. They cover:
  - each model function;
  - the model checked against matrix models with known answers;
  - the IUCN A4 window and generation-time logic;
  - the Pegg runners, the fire pipeline and the range ensemble.

  All pass in about 25 s when run from `plant_ibm/` with
  `/opt/R/4.3.1/bin/Rscript tests/testthat.R`.
- **No frozen regression tests.** The reproducibility tests run the same seed
  twice in one session. Nothing is compared with stored outputs, and
  `tests/reference/` is empty.
- **The working M. nodosa parameters are never run end to end by a test.**
  - The tests load `parlib/unittests/params.R`, where every mechanism is off.
  - That file also predates the current fire and resprout settings. It lacks
    `resprout_delay_dist`, `fire_extent_mode`, `fire_size_sd`,
    `fire_stand_cover` and `rust_site_sd`, so the code fills in old defaults.
  - The real settings in `params.R` have no output test: stand fire mode,
    fixed 2-year resprout window, senescence, juvenile decline, anchor rust
    Beta(14, 6), microsite Beta(20, 5).
  - Phase 1 has to close exactly this gap.
- **`parlib/`:**
  - `unittests/` is test infrastructure, not species biology.
  - `diags/*` is five byte-identical copies of an older calibration. Nothing
    uses them: the diagnostics source `params.R` directly.
  - `ensemble/no_rust_no_fire/param_ranges.R` defines prior ranges, in
    `keys`/`mins`/`maxs`/`Sigma` form, with a written reason for each. That is
    a ready model for the Phase 4 priors. It is used by `runner/morris.R`.
- **Already in place for Phase 2:**
  - All biology is in one nested list returned by `get_default_params()`.
  - `runner/params_utils.R` already maps between that list and flat keys.
  - Response curves are already pluggable through `dose_response`.
  - Three fire modes exist (legacy, plant, stand).
  - Rust suppression of flowering (`rust_flower_dose_response`) and of seed
    set (`rust_recruit_dose_response`) exist but are off. Together they
    partly cover PLAN Phase 3's fruit-infection item.
  - `ricker()` is defined but not wired into `recruit()`.

## Overlaps with ROADMAP.md

- ROADMAP Phase 2 (fitting the Pegg plots, `runner/pegg_sites.R`) overlaps
  PLAN Phase 5 (ABC calibration).
- ROADMAP Phase 3 (Morris, Sobol) overlaps PLAN Phase 6 (sensitivity
  analysis).
- ROADMAP Phase 4 ("configurable paths in place of hard-coded ones") overlaps
  PLAN Phase 2.
- ROADMAP's item on rust effects on flowering and seed set overlaps PLAN
  Phase 3's fruit suppression.
- The resistance hook (ROADMAP 1c) has to stay on placeholder values under
  the `private/` rule.

## What PLAN.md (and CLAUDE.md) should change

1. **The test command in CLAUDE.md rule 3 doesn't work.** From the repo root,
   `Rscript -e 'testthat::test_dir("tests")'` runs the empty `tests/` folder.
   The system `Rscript` is also R 3.4.4, which fails to load `sf`. Pick one
   runner, for example a root script that runs `plant_ibm/tests` plus the new
   regression tests, and state the R version.
2. **The Phase 1 reference scenarios must use `params.R`, not the unit-test
   parameters.**
   - Include a shortened A4-style run and a Pegg-scenario run, since those
     are the outputs used for listing.
   - Set "high/low rust" through `annual_env_alpha/beta`, using the anchor
     cell versus a low-pressure cell. `rust_pressure` is not varied by site.
   - Exact matching depends on the R version and the RNG settings. Run the
     tests under the R version from item 1.
   - `renv` is not activated when running from `plant_ibm/`, because
     `.Rprofile` is at the repo root. `abc` and `future.apply` are also
     missing from `renv.lock`.
   - Any refactor that reorders random draws will break exact matching.
     Section 3.2 of the architecture note lists the draw order.
   - The "mean mortality over 20 seeds" check can use the existing Pegg
     statistic (`compute_pegg_stats`) or `summarise_run()`.
3. **Rule 2 needs a procedure for deliberate reference updates.** The
   ROADMAP's M. nodosa work (Pegg per-tree refit, range_v2) will change
   `params.R` after the references are frozen. Updating the references should
   be a step you do yourself, logged in `methodological_notes.md`.
4. **The Phase 2 inventory must cover `runner/` and `spatial/`, not just
   `R/`.**
   - Items to include: `RUST_ANCHOR`, the 1600/1810–1910/2010/2026/2210
     phasing, the Pegg values (0.48, 2020/2022), `mn_site_params.csv`,
     `QLD_LAT` and `DEFAULT_PARAM_DRAWS`.
   - Split the inventory into three layers: species biology (YAML), site
     inputs, and assessment settings.
   - Section 5.4 of the architecture note is a starting list.
5. **The YAML should mirror the nested list exactly**, so the loader returns
   a list that is `identical()` to the current `get_default_params()`.
   - YAML reads `10000` as an integer, not a double. That changes
     `make_run_id` hashes and invalidates the result caches.
   - Values computed in code (`1/50`, the anchor Beta parameters) need a
     defined form in the YAML.
6. **Rule 5 can't be met for `m_nodosa.yaml` as things stand.**
   `data/rust_species_params.csv` has no M. nodosa rows. Either add them, or
   allow citations to `methodological_notes.md` IDs (E1, V4, D3, …).
7. **parlib needs decisions:**
   - keep `unittests` as a test fixture and add the missing fields (a
     validating loader would otherwise reject it);
   - remove or archive `diags/*` (your call);
   - adopt the `param_ranges.R` format for priors.
8. **PLAN Phase 6 needs new A3 code.** `runner/iucn_a4.R` assesses A4 only,
   and `collect_iucn_results()` skips anything that is not A4.

## Out-of-date documentation (not edited)

`plant_ibm/README.md` and `plant_ibm/TESTING.md` no longer match the code:

- they describe density dependence (`dd_hazard`, `dd_alpha`, `K`) that has
  been removed;
- the README gives the old IUCN N definition (flowering plants only);
- the README says legacy fire is the default;
- TESTING.md says "149 tests across 9 files".

Section 10 of `docs/architecture.md` has the full list.
