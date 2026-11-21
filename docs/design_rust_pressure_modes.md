# Design note: rust-pressure modes

Status: **approved 2026-11-21 with the recommendations in section 6; built** (notes T18). Two
changes made while building: the anchor cell is fixed in the landscape file (`cell_lon`,
`cell_lat`), so species that do not occupy it normalise to the same cell; and the gradient's sign
is set so that g is the change in logit pressure per 10° **south**, shift = g × (lat_anchor − lat_i) / 10
(the formula in §3.2 below, with g < 0, would have raised pressure to the south).

## 1. Purpose

Rust pressure is the pathogen vertex of the disease triangle (C1): one Beta draw per year per cell,
shared by all plants there. It is now set for M. nodosa only. The new species raise three
questions it cannot answer.

| Evidence | Question for the model |
|---|---|
| *R. psidioides*: no climate predictor of infection across rainfall 941–2182 mm, Tmax 18–26 °C, elevation 2–926 m (Fensham et al. 2020) | Should pressure follow the SDM at all, or be flat across the range? |
| *A. beckleri*: clear past infection at Knorrit SF, little or none at Kerewong and Coopernook SF, all south of Port Macquarie (Manea et al. 2019), against heavy impact in SE Qld (Pegg et al. 2017) | Is there a north–south gradient the SDM misses (C4), and does it matter for the listing? |
| *U. australis*: rust first seen on the species in 2018–19, after a long drought, eight years after rust reached NSW | When does pressure start, and does it build up gradually? |

There is also a question of principle. The pressure level is set at Yarringully (0.70, C2) inside
the M. nodosa species file. For other species: is that a fact about the place, or about
M. nodosa?

**Constraint.** M. nodosa is unchanged, bit for bit. Its site table, params and run IDs stay as they
are, and all new settings are optional (T12 §2.1).

## 2. What exists now

- **Cell level** (`spatial/rust_pressure.R`). The mean annual pressure in cell i is
  mu × SDM_i / SDM_anchor, capped at 0.95. mu = 0.70 at the anchor cell (Yarringully), and the
  Beta concentration is 20. The single mode, `sdm`, is chosen by `rust_anchor$mode`
  (registry `RUST_PRESSURE_MODES`, T12). `02c_update_rust_pressure.R` writes `annual_env_alpha`
  and `annual_env_beta` to the species' site table.
- **The SDM saturates on the coast** (C4). 313 of 375 M. nodosa cells have suitability ≥ 0.9, so
  pressure is ≈ 0.70 along the whole coast and falls only towards the dry interior.
- **Run level.** `rust_site_sd` gives each run a persistent site effect on the logit of the cell
  mean (V6 favours a large one). The yearly draw is iid Beta (`draw_annual_env()`).
- **Timing.** One `rust_start_year` (2010) for every run. Pressure is zero before it and at full
  level from it.
- **Host side.** Species differences in how hard rust hits are in `rust_dose_response` and the
  `rust_pressure` multiplier (1 for M. nodosa), not in the environment term.

## 3. Proposal

### 3.1 Pressure is a property of place, shared by all species

- **The principle.** The anchor value (Yarringully, 0.70) and the SDM describe the pathogen and
  environment, not M. nodosa. Every species in a cell sees the same yearly pressure, and species
  differ only through the host side (`rust_dose_response`, `rust_pressure`, resistance), which
  Phase 5 fits.
- **In practice.**
  - The `rust_anchor` block moves from `species/m_nodosa.yaml` to `config/landscape.yaml`,
    keeping its values and citations (`m_nodosa_params.csv` rows).
  - A species file may still give a `rust_anchor$mode` (§3.2), but not its own anchor level.
  - M. nodosa's `RUST_ANCHOR` is the same list as now, so nothing downstream changes.
  - Like the fire tables (T13), every species' cells are on the master grid, and the same cell
    gets the same pressure whichever species it is computed for.
- **What it assumes.** The 0.70 was inferred from M. nodosa incidence, so it already mixes host and
  place. Calling it place fixes a convention: M. nodosa's host terms absorb the rest, and other
  species' host terms are fitted against that same scale. Logged as an assumption.

### 3.2 Spatial modes (the cell's mean pressure)

The registry gains two modes. Each returns the same columns, so everything downstream of
`annual_env_*` is unchanged.

| Mode | Cell mean | For |
|---|---|---|
| `sdm` (default, now) | mu × SDM_i / SDM_anchor, capped | M. nodosa: rust follows climate, mainly coast vs interior |
| `constant` (new) | mu in every occupied cell | species whose infection shows no climate pattern (*R. psidioides*), or whose range lies wholly in saturated SDM cells |
| `sdm_gradient` (new) | the `sdm` value with a shift on the logit scale: logit(mean_i) = logit(sdm mean_i) + g × (lat_i − lat_anchor) / 10 | testing whether a north–south gradient (C4 open item; *A. beckleri*) changes the category |

- **`constant`** has no settings. For the rainforest species it will often give almost the same
  pressures as `sdm`, because their cells are coastal and the SDM is saturated there. It states
  the assumption directly, and it does not depend on an SDM fitted to other hosts' records.
- **`sdm_gradient`**: g is the change in logit pressure per 10° of latitude. The anchor latitude
  keeps exactly mu, so V4 and the Pegg calibration are untouched. g is not known: no data give a
  gradient (C4 found the presence-background gradient follows survey effort). So g is a **run
  setting drawn per run** in the range ensemble, from a range that includes 0 (for example
  U(−1, 0), from no gradient to roughly halving the odds of pressure 10° south). The output is
  how much the category shares move with g. It is a sensitivity question, not an estimate. The
  gradient is applied when a run is built (like `rust_site_sd`), not stored in the site table.
- **`rust_site_sd`** still applies on top of any mode.

### 3.3 Timing: arrival year and build-up

- **`rust_arrival_year`** (species, optional; default: the assessment's `year_rust`, 2010). The
  first year with rust pressure for that species' runs. *U. australis* would cite 2018 (first
  detected 2018–19).
- **`rust_ramp_years`** (species, optional; default 0 = step, as now). Pressure rises linearly
  from 0 to its full level over this many years after arrival. It represents inoculum building up
  locally, and detection lags. It is applied inside `draw_annual_env()`, as a multiplier on the
  year's draw, so the Beta draw itself is unchanged.
- **The A4 window is not moved.** Under A4 it is anchored to the run's own reference period
  (`year_gentime_*`, `ref_*` from the assessment). A later arrival makes the decline start later
  within the same window. Whether a species should use its own `year_rust` in the assessment block
  is the assessment's choice (`species/<sp>.yaml` `assessment`), and this note does not change it.

### 3.4 Where it lives in code

- `spatial/rust_pressure.R`: add `constant` to `RUST_PRESSURE_MODES`. `sdm_gradient` uses the
  `sdm` columns in the site table, and adds the latitude shift when the run is built.
- `runner/site_params.R`: `apply_rust_gradient(p, g, lat, lat_anchor)`, next to
  `apply_rust_site_effect()`. In `plan_range_run()`, g is drawn from the species'
  `ensemble_draws`, or fixed with a runner argument `rust_gradient = g`, and hashed into the run ID.
- `R/simulate.R`: `draw_annual_env()` honours `rust_arrival_year` (falling back to
  `rust_start_year`) and `rust_ramp_years`. It draws no extra random numbers, so a run with the
  ramp absent is identical draw for draw.
- Schema:
  - `rust_anchor` moves to the landscape layer, with `mode` choices
    `[sdm, constant, sdm_gradient]`;
  - `rust_arrival_year` and `rust_ramp_years` are optional species fields;
  - `rust_gradient` is an optional run setting, default 0.

## 4. Not in this note

- **Year-to-year structure** (autocorrelated or regional wet years). `draw_annual_env()` is ready
  for it (T12 §4.2), but nothing in the data asks for it yet.
- **A new rust SDM.** C4 decided to keep the current one. A coastal gradient enters here only as
  the uncertain term g.
- **Host susceptibility** per species (`rust_dose_response`, `rust_pressure`). That is Phase 4
  priors and Phase 5 fitting.
- **Accumulating damage** (repeated defoliation raising death risk). That is a separate Phase 3
  note.

## 5. Verification

- The full suite passes with `tests/reference/` untouched. M. nodosa's site table is
  byte-identical after rerunning `02c` with the anchor moved to the landscape file.
- New tests:
  - `constant` gives mu in every cell;
  - `sdm_gradient` with g = 0 equals `sdm`; the anchor cell keeps mu for any g; the shift has the
    stated sign and size;
  - the ramp gives zero pressure before arrival, a linear rise, then the full draw; it is absent
    for M. nodosa and the run is identical draw for draw;
  - the gradient is hashed into the run ID, and paired plans (g vs 0) share cells, seeds and draws;
  - the anchor read from the landscape file equals the old `RUST_ANCHOR`.
- The decision is logged as a T entry once approved.

## 6. Decisions needed

1. **Pressure as a property of place** (§3.1): one anchor in the landscape file, shared by all
   species. *Recommended.* The alternative is a separate anchor and level in each species file.
   That would need a field-incidence level per species, which no species has; and the level would
   be confounded with host susceptibility anyway.
2. **The modes to add:** `constant` and `sdm_gradient` (*recommended*), or `constant` only.
3. **The gradient as a per-run sensitivity draw** from a range including 0 (*recommended*), or a
   fixed scenario value set by the runner.
4. **Species arrival year and an optional ramp** (§3.3). *Recommended*, for *U. australis*. The
   alternative is the arrival year only, with no ramp.
5. **The default ramp of 0** (a step, as now) for every species unless its file cites a reason.
   *Recommended.*
