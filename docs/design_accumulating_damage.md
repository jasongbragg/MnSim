# Design note: accumulating rust damage

Status: **approved 2026-11-21 with the recommendations in section 6; built** (notes T19, which
records each decision required and the choice made). The fire-reset setting is named
`damage_fire_reset` (`reset` | `keep`), since the schema has no logical type.

## 1. Purpose

The rust hazard has no memory. Each year a plant dies with a probability set by that year's
pressure, its susceptibility and its age; what rust did to it in earlier years does not matter.
The rainforest species die differently. They lose their crowns over successive infection events,
then stop flushing and die.

| Evidence | What it says about the hazard |
|---|---|
| *R. rubescens*, Olney SF exclusion trial (Carnegie et al. 2016): untreated trees near fully defoliated by Dec 2013, stopped flushing and dead by Aug 2014; 53% dead < 4 yr after rust established | Death follows defoliation built up over several years, not a constant yearly risk |
| Same trial, partially treated trees: no change for ~6 months after protection stopped, then incidence and severity rose sharply, then crown loss | There is a lag between infection and loss of canopy |
| *R. rubescens*, paired resurveys 2014 → 2016 (Pegg et al. 2017, Table 10): 0.10 → 0.50 dead at Bongil Bongil, 0.23 → 0.50 at Royal NP, 0.00 → 0.31 at Tallebudgera 2, 0.73 → 0.92 at Gold Creek, 0.25 → 0.30 at Tallebudgera 1 | Annual mortality *among survivors* rose at 4 of 5 sites (e.g. Bongil: ~3%/yr to 2014, ~25%/yr after) |
| *R. psidioides*: 57% dead ~3.5 yr after arrival (Carnegie et al. 2016), > 80% by 2016 (sub-sample, low confidence) | ~21%/yr at first, then ≥ 32%/yr among survivors |
| Crown transparency: *R. psidioides* mean 0.95, *R. rubescens* 0.76 (healthy understorey 0.30–0.35); branch dieback 0.91–1.00 in *A. beckleri* and *D. humile* (Pegg et al. 2017) | Most living plants carry heavy damage: there is a damaged-but-alive state to match |
| *U. australis*: branchlet dieback and death of small stems, then increasing loss of larger stems; ~10% of genets "on the brink" (Kooyman in litt. 2024) | The same progression in a clonal species |

The current model cannot produce rising mortality among survivors. With a steady mean pressure,
a memoryless hazard gives constant risk, and selection removes the most susceptible plants first,
so survivors' risk *falls*. The accelerating Table 10 sites are therefore outside what any
parameter values can give. That is the case for a new structure rather than refitting.

**Constraint.** M. nodosa is unchanged, bit for bit. The new mode is an optional species setting.
Off, it adds no column, no census field and no draws, and the params list and run IDs stay as they are.

## 2. What exists now

- **Rust hazard** (`nominate_deaths()`, `R/mortality.R`):
  p_rust = eff_susc × (age term + resprout term), where
  eff_susc = (1 − resist_score) × env_susceptibility × annual_env_t.
  - The age term falls from `rust_pressure × age_peak` in seedlings to `age_floor` in adults.
  - The resprout term is extra risk while recovering from fire.
- **Combination.** p_rust is combined with background, senescence and juvenile decline as
  independent risks, then one Bernoulli draw per plant.
- **No per-plant rust history.** The nearest things are the resprout state (fire) and the clonal
  `n_killed` pulse (T17).
- **Plug-in point.** `dose_response()` gives the bounded forms (linear, saturating, sigmoid, …)
  used throughout.

## 3. Proposal

### 3.1 A damage state per plant

- **New column.** `rust_damage` D in [0, 1], the share of the crown lost to rust. It starts at 0
  for every new plant: founder, seedling or sucker.
- **Annual update**, after the year's pressure is drawn and before mortality:

      D ← D × (1 − damage_recovery)                    recovery of the canopy
      D ← D + (1 − D) × gain                           new damage on the remaining crown
      gain = dose_response(rust_pressure × eff_susc, damage_gain_dose_response)

  - Damage in year t uses the same disease triangle as the hazard: host × microsite × that year's
    pressure.
  - The update is deterministic, so it uses **no random numbers**.
  - A run of bad years drives D towards 1. Light years let it recover at `damage_recovery`.
  - Before rust arrives, and in years with no pressure, gain is 0.
- **What the lag needs.** Crown loss takes several infection years to build up, so death lags
  infection. That is the multi-year form of the Olney lag. The ~6-month lag itself is below the
  model's annual step.

### 3.2 Death from damage

- **The hazard.** p_damage = dose_response(D, damage_mortality_dose_response). A steep sigmoid
  with half_sat near full crown loss gives "stopped flushing and died".
- **Combination.** It is added as one more independent risk in the union in `nominate_deaths()`.
  The existing p_rust term stays (decision 2): it carries the quick death of seedlings and young
  suckers, which the data also show. A species that dies mainly by damage sets a low `age_floor`.
- **Order.** The mortality draw is the same one Bernoulli per plant as now. Only its probability
  changes, so the random stream keeps its shape.

### 3.3 Fire and resprouting

- **Top-kill clears the damage.** When fire top-kills a plant and it resprouts, its canopy is
  new, so D is reset to 0 (decision 3).
- **What stays.** The existing resprout term still makes the resprouting plant more vulnerable,
  which is how "resprouting gives no escape" (*R. psidioides*) is already represented.
- **The alternative.** D carried through the fire, standing for depleted reserves.

### 3.4 Clonal species

- **Suckers.** They start with D = 0.
- **Rootstock.** It keeps its own `rootstock_hazard_no_canopy` (T17).
- **Not here.** Linking rootstock exhaustion to the damage of its ramets. A later refinement if
  *R. psidioides* sucker pools need it.

### 3.5 Census and outputs

In damage mode the census adds:
- `mean_damage`: mean D among living mature plants;
- `share_damaged`: the share of them with D ≥ `damage_report_threshold` (fixed at 0.6, a census
  setting, not biology).

These sit alongside the existing columns, like the clonal ones, and are matched to crown
transparency and dieback targets in Phase 5. Transparency includes a healthy baseline
(0.30–0.35), so it maps to D through that baseline.

### 3.6 Settings

All are optional, with schema `required_when: {rust_hazard_mode: damage}`.

| Setting | Layer | Meaning |
|---|---|---|
| `rust_hazard_mode` | species | `instant` (default, current) or `damage` |
| `damage_gain_dose_response` {max_effect, form, half_sat, hill} | species | share of the remaining crown lost per year, from rust_pressure × eff_susc |
| `damage_recovery` | species | share of damage recovered per year |
| `damage_mortality_dose_response` {max_effect, form, half_sat, hill} | species | annual death risk from D |
| `damage_reset_on_resprout` | species | `true` (default) or `false` (decision 3) |

They are biology, so rule 5 applies: each is a prior with a range and a reason, or cites a CSV
row, and is then fitted in Phase 5. Likely targets:
- the Table 10 trajectories, as conditional mortality among 2014 survivors;
- Olney's 53% in < 4 yr;
- *R. psidioides* 57% at 3.5 yr;
- crown transparency.

## 4. Not in this note

- **Effects of damage on flowering and fruit.** Crown loss very likely cuts fecundity. That is
  the fruit-infection suppression note. It can read D, which is why D is a column and not a
  hidden local.
- **Growth or size.** Stem age still stands in for size.
- **Sub-annual timing** (the ~6-month lag) and seasonal infection.
- **Damage-dependent infection risk** (damaged plants flushing more and so more exposed). This
  could make gain depend on D, but no data separate it.
- **M. nodosa.** It stays in `instant` mode. Its Pegg data are a two-year post-fire endpoint (V3),
  which do not constrain a damage state. Testing damage mode on M. nodosa would be a separate
  sensitivity run, not a change to its calibration.

## 5. Verification

- The full suite passes with `tests/reference/` untouched. Instant mode is the default, and
  M. nodosa sets nothing.
- New tests:
  - in instant mode a run is identical to the current code, draw for draw, with no new columns;
  - the damage update makes no random draws (the stream after it is unchanged);
  - D stays in [0, 1], rises under constant pressure towards 1, and decays at `damage_recovery`
    when pressure is 0;
  - D is 0 before rust, and in new plants;
  - top-kill with resprout resets D to 0 when the reset is on, and keeps it when it is off;
  - **the case for the structure**: with constant pressure, the death rate among survivors rises
    over the first years in damage mode, but in instant mode it does not rise;
  - `mean_damage` and `share_damaged` appear only in damage mode; the schema's `required_when`
    rules hold.
- A smoke run with illustrative values, set by hand to show the shape only, checks that a
  population declines with a lag after rust arrives.
- The decision is logged as a T entry once approved.

## 6. Decisions needed

1. **A continuous damage state D** with deterministic gain and recovery (*recommended*). The
   alternative is a counter of consecutive high-pressure years (PLAN's "consecutive defoliation
   events"). That is simpler, but it resets abruptly after one light year, and it needs a
   threshold for a "defoliation event" that no data give. D behaves like a smoothed version of
   the counter.
2. **Keep the instant hazard alongside** damage (*recommended*): seedlings and suckers are killed
   quickly, while adults die by damage. The alternative is damage only, which loses the rapid
   death of small plants that the *R. rubescens* and *R. psidioides* data report.
3. **Fire resets damage** for a resprouting plant (*recommended*: its canopy is new; the existing
   resprout term carries the extra risk). The alternative carries D through fire as depleted
   reserves.
4. **Include recovery** (`damage_recovery`, prior including 0) (*recommended*: *A. beckleri* at
   Knorrit had seedlings without dieback after the last infection event, and pressure varies
   from year to year). The alternative is no recovery, so damage only accumulates.
5. **Fecundity and flowering effects of D** go into the fruit-infection note, not here
   (*recommended*).
