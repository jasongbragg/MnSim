# Design note: Olney exclusion trial as a validation of accumulating damage

Status: **approved 2026-11-21 with the recommendations in section 6; hooks and runner built**
(notes T20, which records each decision required and the choice made). Run once the
*R. rubescens* species file exists (Phase 4). One change while building: protection scales the
year's pressure itself, so it covers every rust effect, not only damage gain and the instant hazard.

## 1. Purpose

The damage mode (T19) has a recovery rate, a damage gain and a death risk from damage. Field
mortality data alone cannot separate them: a high gain with fast recovery and a low gain with
slow recovery can give the same death toll. The Olney SF disease exclusion trial
(Carnegie et al. 2016) separates them, because fungicide works like a switch on rust.

- **Treated trees.** Rust is removed, so crown recovery is the only process at work. That gives
  `damage_recovery` directly. It is a rare record of canopy tissue turnover with rust removed.
- **Untreated trees.** They test the full chain: gain, partial recovery in seasons when a flush
  escapes infection, defoliation, then death.
- **Partially treated trees.** Spraying stopped in June 2013, so these trees show damage
  building up again from a nearly healthy crown.

**Species.** The trial is *Rhodamnia rubescens* only. The same paper surveys *R. rubescens* and
*Rhodomyrtus psidioides* across their ranges. For *R. psidioides* it gives crown transparency and
mortality per site, but there is no exclusion trial (section 3.5).

**Constraint.** This is a validation set-up (layer V, like the Pegg validation), not species
biology. M. nodosa is unchanged, and the hook it needs is optional and off by default.

## 2. The data (extracted 2026-11-21)

- **`data/validation/carnegie2016_olney_crown_transparency.csv`.** Mean crown transparency per
  treatment and date, Aug 2011 – Oct 2014, digitised from the plotted means in Fig. 3 (21 dates
  for untreated and treated, 12 for partially treated; ± ~1 percentage point).
  - Design: 20 trees 0.5–4 m tall. 10 were sprayed with triadimenol monthly from Aug 2011 to
    Oct 2014 and 10 were not sprayed. In Jun 2013, five of the sprayed trees stopped being sprayed.
  - Rust was first detected at Olney in Oct 2010.
  - Dead trees are scored 100%.
- **`data/validation/carnegie2016_site_mortality.csv`.** Percentage dead per site in 2014:
  *R. rubescens*, 43 sites (Table 2); *R. psidioides*, 18 sites (Table 4). For between-site
  spread.
- **New rows in `data/rust_species_params.csv`:**
  - healthy crown transparency (*R. rubescens* 0.30–0.35, *R. psidioides* 0.25–0.35);
  - the three trial trajectories in summary form;
  - pointers to both files.

  The existing Olney 53% row is corrected. That figure is a separate count of 100 trees in the
  stand (Table 2, "Olney SF 1"), not the trial trees.
- **Not available in the paper:**
  - the number of trial trees that died (the paper only says "some untreated trees" were dead by
    Aug 2014);
  - tree-level series;
  - the incidence, severity and flush series (Supplementary Fig. 1).

  Tree-level data from the authors (NSW DPI) would allow a proper likelihood. Worth asking.

**What the series already show** (taking a healthy baseline of 0.30):

| Arm | Crown transparency | Damage D = (T − 0.30) / 0.70 |
|---|---|---|
| Treated | 0.66 (Aug 2011) → 0.44 (Dec 2011) → 0.36 (Aug 2012) → ~0.31 (mid 2013) | 0.51 → 0.20 → 0.09 → ~0.02: annual recovery share ~0.8–0.95 |
| Untreated | 0.66 → ~0.83 (mid 2012) → ~0.88 (mid 2013) → 0.93 (Oct 2014), dipping to 0.79 in Dec 2013 | rising, with partial recovery when a flush escapes |
| Partially treated | 0.31–0.33 (Jun–Dec 2013) → 0.38 (Apr 2014) → 0.66 (Aug 2014) | ~0 → ~0.5 in ~8 months once infection took off |

Two caveats:
- The treated trees rose to 0.44 in Aug 2014. That may be a lapse in control: the paper notes
  that spraying less often than monthly failed.
- Recovery under fungicide includes protection of existing leaves as well as new flush.

## 3. Proposal

### 3.1 The treatment hook

An optional run setting, `rust_protection`: a list of (start, end) periods in calendar years,
with fractions allowed.
- **What it does.** In year t, the rust input to both damage gain and the instant hazard is
  multiplied by 1 − (the share of year t that is protected).
- **Where.** It is applied in `update_rust_damage()` and in `nominate_deaths()`'s rust term, as a
  scalar.
- **Random numbers.** It makes no draws.
- **The arms.**
  - Treated: 2011.6–2014.8.
  - Partially treated: 2011.6–2013.45.
  - Untreated: none.

  Each arm is a separate run (decision 1).

### 3.2 Shared yearly pressure across arms

All three arms had the same weather. An optional run setting, `annual_env_series` (one value per
year), replaces the yearly Beta draw. The validation runner draws the series once per replicate
and gives it to all three arms, so they differ only by treatment (decision 2). When the setting
is absent, `draw_annual_env()` is unchanged.

### 3.3 From damage to crown transparency

- Crown transparency per tree is T = b + (1 − b) × D for living trees, and 1 for dead trees.
- The model's trial value is the mean of T over the arm's trees.
- b is a prior from the new CSV rows: [0.30, 0.35] for *R. rubescens* (decision 4).

### 3.4 Runs and comparison

- **Trees.** Each arm is a cohort of 10 trees (5 for partially treated) of the trial's size,
  inside the Olney stand. The stand is simulated from before rust arrived, so competition and the
  instant hazard act as they would in the field. Replicates give the spread of a 5- or 10-tree
  mean.
- **Rust.** It starts in Oct 2010 (`rust_start_year` 2010; with an annual step, the 2010–11
  season is the first year of rust). The trial's starting value (0.66 in Aug 2011, the same in
  all arms) is itself a target: about one season of infection should bring D to ~0.5.
- **Resolution: annual** (decision 3). Observations are compared once a year, at the August
  values (the trial's start month) of 2011–2014. Oct 2014 is the end point. The ~6-month lag and
  the seasonal dips are below the model's step and are not targets.
- **Targets** for the Phase 5 calibration code (rule 4: nothing fitted by hand):
  - the three arms' August trajectories;
  - the stand's 53% dead by Oct 2014;
  - together with the Table 10 resurveys, the Table 2 site spread, and mean crown transparency
    0.76 with 79% of trees > 0.60.
- **Code.** `runner/olney_validation.R`, with its set-up in
  `species/r_rubescens_validation.yaml` (block `olney`).

### 3.5 *R. psidioides*

- **What the paper gives.** Crown transparency (mean 0.95, 82% of trees > 0.90) and mortality at
  18 sites, against a healthy baseline of 0.25–0.35. There is no trial.
- **Recovery.** Its recovery rate is not identified by its own data. Most of its trees sit near
  full damage, which says little about recovery.
- **Proposal** (decision 5). Its `damage_recovery` prior is the *R. rubescens* posterior range,
  logged as an assumption: both are evergreen rainforest-margin Myrtaceae and are similarly
  susceptible. Its gain and death risk are fitted to its own survey and site data.

## 4. Not in this note

- A sub-annual time step.
- Leaf-size effects. Smaller immature leaves on untreated trees (Table 1) are a further cost of
  damage that is not modelled.
- Flowering and fruit. None were produced in the trial. That belongs to the fruit-infection note.
- Tree-level fitting. That needs the authors' data.

## 5. Verification

- The full suite passes with `tests/reference/` untouched. Both hooks are absent for M. nodosa,
  and runs without them are identical draw for draw.
- New tests:
  - full protection gives no damage gain and no instant rust hazard; partial-year protection
    scales both;
  - `annual_env_series` replaces the draw and makes no rbeta draws;
  - the transparency mapping is correct, with dead trees = 1;
  - the three arms of one replicate share the same yearly pressures;
  - (once the species file exists) a smoke run gives the trial's qualitative pattern: the treated
    arm falls towards b, the untreated arm rises, and the partially treated arm turns upwards
    after protection ends.
- The decision is logged as a T entry once approved.

## 6. Decisions needed

1. **One run per arm**, with a run-level protection schedule (*recommended*: the trial trees do
   not interact in any way that matters over three years). The alternative is a per-plant
   protection flag within one run, which needs a new column.
2. **A shared yearly pressure series across arms** (*recommended*: the arms had the same weather,
   so comparisons between arms are not blurred by different draws). The alternative is
   independent draws per arm.
3. **Annual comparison at the August values** (*recommended*). The alternative is a sub-annual
   step, which would be a much larger change.
4. **Transparency = b + (1 − b) × D**, dead = 1, with b a prior in [0.30, 0.35]
   (*recommended*). The alternative is fixing b at 0.30.
5. ***R. psidioides* recovery prior** from the *R. rubescens* posterior, logged as an assumption
   (*recommended*). The alternative is a wide independent prior.
6. **Order of work.** Build the two hooks and the runner now, with toy tests, and run the
   validation once the *R. rubescens* species file exists in Phase 4 (*recommended*).
   Alternatively, wait and build everything in Phase 5.
