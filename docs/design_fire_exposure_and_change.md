# Design note: species fire exposure and fire-change scenarios

Status: **approved 2026-11-20 with the recommendations in section 6; built** (notes T15).

## 1. Purpose

Two needs, raised for the rainforest species (*R. rubescens*, *R. psidioides*, *A. beckleri*,
*D. humile*):

1. **Exposure.** A cell's fire rate is all fire recorded in a roughly 10 km cell, whatever the
   vegetation. Rainforest and its margins burn far less than the surrounding cell, so the cell rate
   overstates these species' fire. They need a way to run with fire much lower, or off. That should
   be justified from the vegetation type, not fitted.
2. **Change.** Fire in these vegetation types may become more common, as in 2019–20. We want to
   ask whether more fire, together with rust, changes the IUCN category. Data cannot pin the size
   of the change, so this is a scenario: a baseline compared with one or more stated increases. It
   is not a prediction.

**Constraint.** M. nodosa is unchanged, bit for bit. All new settings are optional fields (T12,
§2.1). When they are absent the model makes the same draws as now, and the params list and run IDs
are unchanged.

## 2. What exists now

- **Fire per cell.** Each cell's `fire_prob_annual` and `fire_p_fimp` come from the cell's record
  (D1–D5). `build_site_params()` copies them into a run.
- **`fire_scale`.** A runner argument (`run_range_ensemble(..., fire_scale = s)`) that multiplies
  the cell's fire probability for **the whole run**. `fire_scale = 0` turns fire off.
- **Pairing.** Ensembles that differ only in `fire_scale` pair up run by run, because the cell,
  seed and parameter draws do not depend on it.
- **No change over time.** `R/simulate.R:123` draws a fire each year with the same probability.
  `run_iucn_a4()` clears `fire_years`, so assessment runs cannot be given a fire schedule.

## 3. Proposal

### 3.1 Species fire exposure (species biology)

- **New field.** An optional species parameter `fire_exposure`, a number from 0 to 1 with default
  1. It multiplies the cell's fire rate for the whole run, spin-up included:

  P(fire in year t) = min(1, fire_prob_annual × fire_exposure × ...)

- **Where it lives.** It is a property of the species' habitat, so it goes in the species YAML
  under rule 5. It must either cite a CSV row or be a prior, with a range and a reason. For
  example, "rainforest and wet margins; NSW fire history shows such patches rarely burn".
- **Range of runs.** Draws over its prior range in the range-wide ensemble, using
  `ensemble_draws`, give results across exposures from almost none upwards. Turning fire off is
  `fire_exposure: 0`, or a range starting at 0.
- **`fire_scale` stays as it is.** It remains a management or sensitivity scenario applied at run
  time. Exposure is biology and belongs in the species file. The two multiply together.
- **Default.** M. nodosa omits the field, so it defaults to 1 and nothing changes.

### 3.2 Fire-change scenario (a scenario setting)

**New settings.** Optional settings that change the fire regime from a calendar year:

| Setting | Meaning | Default |
|---|---|---|
| `fire_change_year` | first year with the changed regime | none (no change) |
| `fire_change_rate` | multiplier on the fire probability from that year | 1 |
| `fire_change_size` | multiplier on the mean burnt fraction `fire_p_fimp` (bigger fires; capped below 1) | 1 |
| `fire_change_ramp` | years over which the multipliers rise linearly to full size (0 = step) | 0 |

**How a run uses them.**
- Fire occurrence becomes `min(1, fire_prob_annual × fire_exposure × fire_scale × rate(t))`.
- The fire response is given `fire_p_fimp × size(t)`, so a "stand" fire is more likely to reach
  the stand and has a higher kill ceiling (D3). Here `rate(t)` and `size(t)` are 1 before the
  change year and then ramp up to their full values.
- In code, a small `fire_regime_at(params, t)` function, called each year in `run_simulation()`,
  gives the year's probability and the parameters passed to `apply_fire()`.
- The fire draw stays one `rbinom` per year, as now.

**Where the settings come from.**
- They are scenario settings, not biology: runner arguments, like `fire_scale`.
  `run_range_ensemble(..., fire_change = list(year = 2026, rate = 2, size = 1.5))` writes them into
  each run's params, so they are hashed into run IDs and cached results never mix scenarios.
- Ensembles that differ only in the scenario still pair up run by run.

**Timing.**
- The change year must come after the generation-time window ends (validated). Otherwise the
  window would see changed fire and T_g would no longer describe the pre-rust population (A3, T11).
- The default proposal is 2026, the assessment year. The cell records run to 2024, so they
  already include 2019–20, and a change from 2026 means "more than the recent record".
- The A4 window runs from 2009 to 2009 + 3 T_g, capped at 2126, so a change from 2026 falls inside
  it.

**Reporting.** Each scenario gives the usual spread of runs across categories. A comparison plot
shows baseline against change for the same cells and seeds, so the difference comes from fire
alone.

## 4. Not in this note

- **Fire response of the rainforest species.** Kill, resprouting and recovery parameters are
  species biology. They need priors in Phase 4, and may matter as much as fire frequency, since
  many rainforest Myrtaceae are thought to resprout weakly after hot fire.
- **Fire changing with climate in a graded or spatial way.** For example, rate tied to a projected
  fire-weather index per cell. That would need a projection data layer, and the scenario
  multipliers here are its simple stand-in.

## 5. Verification

- The full suite passes with `tests/reference/` untouched. All fields are absent for M. nodosa.
- New tests:
  - `fire_exposure` 0 gives no fire, and 0.5 halves the yearly probability;
  - before `fire_change_year` a run is identical to the baseline, draw for draw;
  - after it, the probability and `fire_p_fimp` are multiplied, ramped if asked;
  - a change year inside the generation-time window is rejected;
  - the scenario changes the run ID;
  - paired ensembles share cells, seeds and draws.
- The decision is logged as a T entry once approved.

## 6. Decisions needed

1. **Exposure as a species parameter** (§3.1), with `fire_scale` kept as a separate run-time
   scenario. *Recommended.* The alternative is to use `fire_scale` alone for the rainforest species.
2. **The shape of the change.** A step with an optional linear ramp (*recommended*), or a step only.
3. **Fire size as well as frequency** (`fire_change_size`). *Recommended,* since 2019–20 was
   unusual for its extent and severity as much as its frequency.
4. **Default change year.** 2026, the assessment year (*recommended*), or 2020, so that the record
   period is treated as already changed.
