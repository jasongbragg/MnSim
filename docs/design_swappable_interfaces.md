# Design note: swappable fire response, rust pressure and recruitment

Status: **approved 2026-11-20 with the recommendations in section 8; built** (notes T12).
PLAN.md Phase 2, last item.

## 1. Purpose and scope

Phase 3 will add new mechanisms: rust-pressure modes, a clonal sucker stage, fruit-infection
suppression, and possibly other recruitment curves. Each mechanism should plug into a fixed place
in the code, chosen per species by a setting in its YAML. It should not be a new `if` branch
inside the annual loop.

This note covers **only the plumbing**. In Phase 2:
- the existing behaviour is moved behind three seams, with no new biology;
- the only new options are those whose code already exists: the three fire variants, and
  `ricker()`.

The Phase 3 mechanisms each get their own design note, and plug into these seams.

**Hard constraint:** M. nodosa must stay bit-for-bit identical. That means:
- the same model output, so every exact-match regression test passes;
- the same random-number draws in the same order (architecture §3.2);
- the same `get_default_params()` list, so the frozen-params check passes and run IDs and cached
  results stay valid.

## 2. The common pattern

Each seam is the same three things.

1. **A selector field** with a fixed set of choices (`choices:` in `config/schema.yaml`), set in
   the species YAML.
2. **A registry.** This is a named list of functions in the module's R file, all with the same
   signature, for example:

   ```r
   RECRUITMENT_CURVES <- list(bevholt = bevholt, ricker = ricker)
   ```

   The model looks up `REGISTRY[[params$selector]]` and calls it.
3. **A test** that the schema's `choices` and `names(REGISTRY)` are the same set. This keeps the
   two from drifting apart.

**Why a registry, rather than S3 classes or a function stored in `params`:**
- `params` is a plain list that is hashed into run IDs and read from YAML, so it cannot hold
  functions.
- A named list of functions matches the existing style, and costs one lookup per call.

### 2.1 Optional fields (needed for the hard constraint)

Today every schema field is required. Adding a selector to the schema would therefore add it to
M. nodosa's params list. That would break the frozen-params check and change every run ID.

Proposal: a schema field may be declared `optional: true, default: <value>`.
- If a species file omits it, the loader leaves it **out** of the params list.
- If a species file sets it to its default, the loader also leaves it out. This keeps the list
  canonical: writing the default explicitly gives the same run ID as omitting it.
- Model code reads the field with the default as a fallback, `params$x %||% default`. The code
  already does this for `fire_extent_mode` and `resprout_delay_dist`.

M. nodosa's YAML sets none of the new fields, so its list is unchanged.

Parameters that only one option uses (for example `ricker_r` for Ricker) are also optional. The
loader requires them only when that option is selected. Like every other value, they must cite a
CSV row or be a prior (rule 5).

## 3. Fire response (`R/fire.R`)

**Today.** `apply_fire()` (`R/fire.R:62`) has three inline branches, chosen by two fields:
- `fire_p_fimp == 0` → legacy (flat `fire_kill_prob`);
- otherwise `fire_extent_mode`: `"stand"` or `"plant"`. An absent field means `"plant"`.

After the branch, all three share the same code for killing plants and for the resprout delays.

**Proposal.** Split `apply_fire()` into a selector and a shared tail.
- **Response functions.** Each one decides who is reached, who dies and who survives. Signature:

  ```r
  fire_response(pop, alive_idx, params, t) -> list(dead_idx, surv_idx, reached)
  ```

  ```r
  FIRE_RESPONSES <- list(legacy = fire_legacy, plant = fire_plant, stand = fire_stand)
  ```

  The three bodies are cut out of today's branches unchanged, so they make the same draws in the
  same order.
- **Selector.** The existing `fire_extent_mode` stays the selector, so no new field is needed.
  `fire_p_fimp == 0` still selects `legacy`, whatever the mode.
- **Shared tail.** It stays in `apply_fire()`: record the deaths, then give survivors their
  resprout delay (`rust_resprout_delay()`, fixed or Poisson base).

Later options plug in as new response functions, for example a seasonal or intensity-class
response, or a `sucker` response that sends burnt plants to a root-sucker pool (Phase 3 clonal
stage). Each would add a schema choice and a registry entry.

**Not in scope.** Fire *occurrence* (`fire_prob_annual` / `fire_years` in `R/simulate.R`) and the
cell fire regimes. Those are landscape inputs, not species responses.

## 4. Rust pressure

Rust pressure enters at two places.

### 4.1 Site level: a cell's mean annual pressure (`spatial/rust_pressure.R`)

**Today.** `rust_pressure_from_sdm()` (`spatial/rust_pressure.R:34`) sets the pressure of each grid
cell to

  mean = min(cap, max(floor, mu × SDM_cell / SDM_anchor))

and writes `annual_env_alpha` / `annual_env_beta` to the site table. `build_site_params()` copies
these columns into a run's params.

**Proposal.**
- Add an optional field `mode` to the `rust_anchor` block of the species YAML. It defaults to
  `sdm`, and M. nodosa omits it, so `RUST_ANCHOR` is unchanged.
- Add a registry of cell-pressure functions:

  ```r
  RUST_PRESSURE_MODES <- list(sdm = rust_pressure_from_sdm)
  ```

  Each function takes the site table and the anchor, and returns the same four columns.
- The site-table builder (`02c`) calls the selected function. Everything downstream of the
  `annual_env_*` columns stays the same.

Phase 3's "saturating/constant" mode (PLAN Phase 3: *R. psidioides*, *A. beckleri*) then becomes
one more function and one more choice. Its biology belongs to the Phase 3 design note, not here.

### 4.2 Annual draw (`R/simulate.R:63`)

**Today.** Each year after rust arrives, the loop makes one inline draw,
`rbeta(1, annual_env_alpha, annual_env_beta)`.

**Proposal.** Move that draw into a function, `draw_annual_env(params, t)`, behind an optional
selector `annual_env_mode` that defaults to `beta_iid`. It makes the same single `rbeta` call.

This gives a place for year-to-year structure later, such as autocorrelated or regional wet
years. No new mode is proposed now. If you would rather not add a selector without a second
option, I would extract the function and leave out the selector.

## 5. Reproduction and recruitment (`R/recruitment.R`, `R/simulate.R`)

**Today.**
- `recruit()` (`R/recruitment.R:96`) works out the expected number of recruits as

  bevholt(effective flowering plants) × (1 − shade suppression)

  The first term includes rust suppression of fecundity.
- It then draws the number of recruits from a Poisson, samples parents, and makes the recruits
  inherit their genotypes.
- `ricker()` (`R/recruitment.R:22`) is defined but never called.
- Flowering suppression by rust sits inline in the annual loop (`R/simulate.R:109`).

**Proposal.**
- **Recruitment curve.** Add an optional selector `recruitment_curve`, defaulting to `bevholt`:

  ```r
  RECRUITMENT_CURVES <- list(bevholt = bevholt, ricker = ricker)
  ```

  It is called with the same argument as today, the effective number of flowering plants. Shade
  and rust modifiers stay outside the curve, as multipliers, so they apply to every curve.
- **Ricker.** `ricker()` needs `ricker_r` and `K`, which are not in the schema. Both would be
  optional, and required only when `ricker` is selected (§2.1).

  Its current form, N·exp(r(1 − N/K)), returns recruits in proportion to flowering plants. Its
  scale is not comparable with `R_max` / `K_half`. A species that selects it needs its own
  calibration, which is Phase 5. **Recommendation:** wire Ricker in as a choice, since PLAN asks for
  it, but parameterise no species with it now.
- **Expected-recruits function.** Move the expected-recruits calculation out of `recruit()` into
  `expected_recruits(pop, t, params, annual_env_t)`. That gives one place where later fecundity
  rules plug in, for example fruit-infection suppression (Phase 3).
- **Flowering function.** Move the inline flowering-suppression block into
  `update_flowering(pop, t, params, annual_env_t)`. It makes the same draws, and gives the same
  kind of plug-in point for flowering.
- **Not here.** A second recruitment *pathway*, such as root suckers with no seed input, is a new
  state and part of the clonal-stage note (Phase 3), not a curve choice.

## 6. What does not change

- **Mortality.** `nominate_deaths()` already uses `dose_response()`, whose form can be swapped.
  Accumulating-damage hazard is a Phase 3 item with its own note.
- **The annual step order and the census.**
- **Resistance genetics.**

## 7. Verification

- The full suite must pass with **no change to `tests/reference/`**. That covers:
  - the exact-match regression tests;
  - the frozen-params check;
  - the 20-seed statistics check.
- New tests:
  - each registry matches its schema `choices`;
  - an unknown selector value is rejected when the configuration is loaded;
  - setting a default explicitly gives the same params list and run ID as omitting it;
  - an option's own parameters are required only when it is selected;
  - with `recruitment_curve: ricker` on a test fixture, `expected_recruits()` matches
    hand-computed values;
  - each fire response, called on a fixed population and seed, matches today's branch.
- `docs/architecture.md` §4 and §8 are updated. The decision is logged in the methodological notes
  as T12 once approved.

## 8. Decisions needed

1. **Optional schema fields (§2.1).**
   - *Recommended:* fields that may be omitted, with defaults left out of the params list, so
     M. nodosa's list and run IDs are unchanged.
   - *Alternative:* write every selector into `m_nodosa.yaml`. That changes the frozen params, so
     it needs a reference update (T3 procedure, which you install) and invalidates cached results.
2. **The fire selector.** Keep `fire_extent_mode` as the selector, with `fire_p_fimp == 0` →
   legacy as now (*recommended*). The alternative is a new explicit `fire_response` field that
   includes `legacy`.
3. **Annual-draw selector (§4.2).** Add `annual_env_mode` now, or only extract the function.
   *Recommended:* only extract the function, and add the selector when a second mode is designed.
4. **Ricker (§5).** Wire it in as an option with no species using it (*recommended*), or leave it
   out until a species needs it.
