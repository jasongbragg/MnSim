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
- Fire: NSW (NPWS 1987-2024) and Queensland (fire scars 2014-2024) on one
  footing, wildfire by default; all fire and planned burns kept per cell for
  comparison and management scenarios (Nov 2026; notes D1-D5).
- Not yet included: resistance (all plants fully susceptible), rust effects
  on seed set (option, off).

## Phase 0: housekeeping (before any fitting)

- [x] **Pegg runs use the Yarringully cell's rust pressure.** Single anchor
      (`RUST_ANCHOR` in `params.R`); SDM normalised at the anchor cell;
      `max_effect` recalibrated to 0.53 (range 0.32-0.78). Nov 2026.
      Notes: C2, V4.
- [x] **Separate fire extent from fire intensity.** Per-fire size drawn
      from the cell's burnt-fraction distribution; the fire reaches the
      stand with probability equal to its size, burns 90% of plants, and
      kills old plants up to 0.30 x size. Nov 2026. Notes: D3.
- [x] **Retire `spatial/03_range_wide_iucn.R`**: replaced by
      `runner/range_ensemble.R` (per-run RNG streams, plan table, category
      shares by region, drivers plot); runscript
      `runscripts/range_ensemble.run.code`. Nov 2026.
- [x] **Smaller model checks:** census counts fire deaths (`deaths_fire`);
      Weibull hazard converted to an annual probability; resprout delay
      through the disease triangle; recruits documented as first-year
      survivors (no code change); calendar birth years + `founder` flag.
      Pegg median 0.492 after the changes (max_effect unchanged). Nov 2026.
      Notes: C5, E1, E2, E2b.
- [x] **Speed.** Dead plants archived once they outnumber the living, and
      recruits appended column-wise instead of with rbind.data.frame:
      ~6.8x faster (610-yr run 53 s -> 8 s; full A4 run ~8-9 s), with
      identical output to the previous code (tested). Also fixed parent
      sampling when only one plant flowers. Nov 2026.

## Phase 1: data inputs

### 1a. Queensland fire
- [x] Acquire Qld fire scars (QRSC, 1987-2025) and QPWS fire history
      (planned burns): `spatial/fire/01_download_qld.R`. Nov 2026.
- [x] Restructure fire processing into one adapter per source, each giving
      the burnt fraction per cell per year (`spatial/fire/`); estimator
      downstream of all of them, prior now fitted by maximum likelihood.
      Nov 2026. Notes: D1, D2.
- [x] Border check: records agree using Qld 2014-2024; automated 1987-2013
      Qld mapping rejected. Nov 2026. Notes: D4.
- [x] Fire regime for all 375 cells (Qld from data, not the prior);
      wildfire default, all fire and planned burns alongside. Nov 2026.
      Notes: D5.
- [x] Refresh `fire_size_sd` in `params.R`: 0.125 -> 0.146 (range-wide
      median `sd_frac_burned`, new table). Nov 2026.
- [ ] Rerun the range ensemble with the new site table (`results/range_v2/`).
- [ ] Management scenarios: planned burning as a second fire process (own
      frequency, extent and possibly severity), e.g. no burning vs burning
      a target site every N years; suppression via `fire_scale` already
      available. Notes: D5.

### 1b. Rust SDM (mrsdm)
- [x] Run `evaluate_sdm.R`: Boyce 0.93, block CV TSS 0.47-0.60, clamping
      negligible, coastal climate gradient not passed on by the SDM. Nov 2026.
- [x] Saturation: current SDM kept (sharp decline with aridity, uniform
      along the coast). A random-background refit gives a coastal gradient
      that follows survey effort, not climate (`probe_sdm_bias.R`). Nov 2026.
      Notes: C4. Optional later: target-group background; incidence data;
      coastal gradient as an uncertain term in Phase 3.
- [x] Fill the 16 cells with no SDM value (1.7% of abundance): SDM averaged
      over each cell's land instead of read at the centre. Nov 2026. Notes: B4.
- [ ] Optional: carry the spread among SDM ensemble members into the runs.

### 1c. Resistance prevalence (once experimentally confirmed)
Status Nov 2026: range-wide seedling screening under way (unpublished);
resistance appears rare and clustered; single-locus model kept.
Parked until the data are in; expected to change results modestly. Notes: E3.
- [ ] Map measured resistance onto the model's genetics: number of loci,
      dominance, effect sizes (`resist_score`).
- [ ] Interpolate site frequencies to grid cells; enable `RESISTANCE HOOK`
      in `build_site_params()`.
- [ ] Draw resistance frequencies from their measurement uncertainty.

## Phase 2: fitting to data

Aim: fit the level and the spread between sites, using as many monitored
sites as possible, at coarse resolution (annual or per observation
interval; no monthly dynamics). Notes: V6.

- [x] Fire-only resprout window fixed at 2 years (was Poisson, mean 3).
      Nov 2026. Notes: D3.
- [x] Multi-plot fitting framework (`runner/pegg_sites.R`): seven Pegg et
      al. plots, interval mortality, mixture likelihood over simulated
      sites; persistent site effect on rust pressure (`rust_site_sd`) as an
      alternative to year-to-year variation. Interim run: persistent effect
      strongly preferred, best at the grid edge. Nov 2026.
- [ ] Full grid on jupiter (`runscripts/pegg.sites.fit.run.code`).
- [ ] Per-tree data from G. Pegg (requested Nov 2026): coordinates, fire
      dates, per-visit status, infection, flowering/capsules, seedlings,
      later revisits. Replace the interim plot table; refit.
- [ ] Decide where the site effect acts (rust pressure vs tolerance) once
      infection records are in; carry the fitted rust_site_sd and max_effect
      into the range ensemble (already wired; 0 until fitted).
- [ ] Flowering (after mortality): share of surviving resprouters
      flowering by post-fire year against a pre-rust baseline of ~0.7
      (Hewitt et al. 2014a); rust effect on flowering/seed set (now off);
      consequences for recruitment.
- [ ] Other targets as available: rust incidence levels; resistance
      prevalence (1c); stand structure; fire-killed fraction.

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
