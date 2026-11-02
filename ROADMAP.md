# MnSim roadmap: to a completed IUCN A4 assessment for *Melaleuca nodosa*

Status as of November 2026. Tick items off as they land; add the commit or
date after each. Assumptions and decisions: `docs/methodological_notes.md`.

## Where things stand

- Model: individual-based, annual step; fire (kill or resprout), myrtle rust
  via the disease triangle (host resistance x microsite x annual pathogen),
  Beverton-Holt recruitment with Mendelian inheritance of resistance.
- Assessment: IUCN criterion A4, one fixed window per run (11-yr mean before
  rust arrival vs 11-yr mean 3 generations later, T_g from the same run);
  mature individuals = alive plants of reproductive age (flowering-only
  count kept alongside).
- Validation: Pegg et al. 2025 (Yarringully NR), 24/50 resprouters dead by
  Jan 2022; resprout rust mortality (`rust_dose_response$resprout$max_effect`)
  0.53 at the Yarringully cell's rust pressure (Nov 2026).
- Range-wide ensemble (`plant_ibm/talk/make_talk_figures.R`): cells drawn by
  abundance x seed x parameter draws. Preliminary result (Oct 2026, 300 runs):
  99% of runs threatened, split VU 45% / EN 54%, median decline 51%.
- Not yet included: Queensland fire data (prior only), resistance (all plants
  fully susceptible), rust effects on seed set (option, off).

## Phase 0: housekeeping (before any fitting)

- [x] **Pegg runs use the Yarringully cell's rust pressure.** Single anchor
      (`RUST_ANCHOR` in `params.R`); SDM normalised at the anchor cell;
      `max_effect` recalibrated to 0.53 (range 0.32-0.78). Nov 2026.
      Notes: C2, V4.
- [ ] **Separate fire extent from fire intensity.** `fire_p_fimp` is the
      burnt fraction of a cell but also scales the old-plant kill ceiling
      (c x p), so partly burnt cells get gentler fires. Model change.
- [ ] **Retire `spatial/03_range_wide_iucn.R`** in favour of the run design in
      `talk/make_talk_figures.R` (per-run RNG streams; plan table instead of
      the seed/join reconstruction).
- [ ] **Smaller model checks:** yearly census deaths exclude fire deaths;
      Weibull hazard used as a probability; rust delay to resprouting ignores
      the disease triangle; no first-year seedling mortality; inconsistent
      `birth_year` encoding for founders.
- [ ] **Speed.** Dead plants are never removed from the population table,
      so runs slow as they go (~45 s per A4 run). Archive dead rows (keeping
      what `compute_gen_time()` and the Pegg statistics need). Needed before
      Phases 2-3, which take thousands of runs.

## Phase 1: data inputs

### 1a. Queensland fire
- [ ] Acquire Qld Landsat/Sentinel-2 fire-scar rasters; QPWS fire polygons
      if available (to separate wildfire from planned burns, which the
      rasters do not distinguish).
- [ ] Restructure `02` into one adapter per source, each producing burnt
      fraction per grid cell per year; the existing Gamma-Poisson estimator
      runs downstream of all of them.
- [ ] Calibrate the Qld product against NSW NPWS records in border cells so a
      "wildfire" means the same thing in both states.
- [ ] Re-derive `fire_prob_annual` and burnt fraction for the ~109 Qld cells
      (about 10% of abundance weight).

### 1b. Rust SDM (mrsdm)
- [ ] Run `evaluate_sdm.R`: Boyce index, spatial block CV, clamping check.
- [ ] Address saturation: suitability is ~1 in almost every coastal cell, so
      the SDM separates coast from inland but gives no gradient along the
      coast. Options: steeper suitability-to-pressure mapping; calibration
      against field incidence at several sites; humidity / leaf-wetness
      covariates.
- [ ] Fill the 16 cells with no SDM value (1.7% of abundance), currently
      dropped.
- [ ] Optional: carry the spread among SDM ensemble members into the runs.

### 1c. Resistance prevalence (once experimentally confirmed)
- [ ] Map measured resistance onto the model's genetics: number of loci,
      dominance, effect sizes (`resist_score`).
- [ ] Interpolate site frequencies to grid cells; enable `RESISTANCE HOOK`
      in `build_site_params()`.
- [ ] Draw resistance frequencies from their measurement uncertainty.

## Phase 2: fitting to data

- [ ] Assemble targets as coarse summaries, not detailed curves:
  - Pegg: ~48% of resprouters dead after 2 yr; ~15% of survivors
    flowering at Jan 2022 (4 of 26).
  - Rust incidence levels (Pegg Fig. 3, % trees infected), for rust pressure.
  - Resistance prevalence (Phase 1c).
  - Stand structure / age at first flowering, if available.
  - Share of trees killed outright by fire (literature, if reported).
- [ ] Decide which parameters are free; fix the rest from literature, with
      sources noted in `params.R`.
- [ ] Fit by simulation (ABC or history matching). Accepted parameter sets
      replace the hand-set ranges in `PARAM_DRAWS`.

## Phase 3: sensitivity analysis

Aim: which measurements would most improve confidence in the projection.

- [ ] Morris screening (`runner/morris.R`) on the A4 outputs: decline and
      P(threatened). **Run an early screen with the current model alongside
      Phase 1**, to help prioritise data effort.
- [ ] Variance decomposition (Sobol) on the top 5-8 parameters, likely via a
      Gaussian-process emulator trained on simulation runs. Partition
      outcome uncertainty into parameters, place (cell), and chance (fire
      timing, rust years).
- [ ] Ranked list of what to measure: for each uncertain quantity, how much
      an exact value would narrow the category distribution.

## Phase 4: final assessment and write-up

- [ ] Final ensemble: fitted parameters, all data inputs, resistance.
      Category distribution overall and by region.
- [ ] Robustness: flowering-only vs all mature individuals; with/without
      resistance; no-rust counterfactual (rust vs fire attribution);
      averaging-window length.
- [ ] Methods write-up; README update; configurable paths in place of
      hard-coded ones; tests for the spatial pipeline; remove `old/` and
      stale parameter files.
