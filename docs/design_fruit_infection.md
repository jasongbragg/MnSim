# Design note: rust suppression of flowering and fruit

Status: **approved 2026-11-21 with the recommendations in section 8, except decision 5: the
alternative mature count `N_fruiting_capable` is not built** (the IUCN count stays as in notes A4;
a prospective "will never reproduce" test is a later item). Built: notes T21, which records each
decision and the choice made.

## 1. Purpose

Rust attacks the reproductive chain as well as the canopy. For several of the study species,
losing seed may matter as much as adult deaths. Recruitment drives recovery after rust, and fruit
failure takes susceptible mothers out of the gene pool, so it also matters for the IUCN
projection. The model has a basic version of this (section 3). This note proposes how to extend it
for the species that need it, so it is implemented properly when their species files are built
and calibrated.

Rust can act at three stages:
1. **Flowering (bud set).** The plant does not flower: bud or flush tissue is infected, or the
   crown is too damaged to produce reproductive growth.
2. **Fruit maturation.** Flowers form, but fruit is infected and aborts or does not mature
   (Tucki Tucki).
3. **Seedlings.** Infected seedlings die at the cotyledon stage. The age-dependent rust hazard
   already covers this. It is out of scope here.

## 2. Evidence (data/rust_species_params.csv)

| Species | Evidence | Stage | Note |
|---|---|---|---|
| *R. rubescens* | Tucki Tucki NR: all fruit on untreated branches infected, none matured (`fruit_survival_on_infected_branches` = 0) | fruit | Branch-level fungicide trial. Fruit fails even where the tree lives. |
| *R. rubescens* | Severe decline in fruit production since 2012 (`fruit_production_trend`) | both | Expert opinion, low confidence. |
| *R. rubescens* | Olney trial trees produced no flowers or fruit during the trial (Carnegie et al. 2016; notes T20) | flowering | Trees 0.5–4 m, so possibly below reproductive size. Weak evidence. |
| *R. psidioides* | No mature fruit seen since rust arrived; no flowering or seed since 2010 (`fruit_production_post_rust` = 0) | both | Calibration target, high confidence. |
| *U. australis* | No significant fruiting since 2018; uredinia on fruit (`fruit_production_post_rust`) | both | Seed recruitment was already ~0 before rust (no verified wild seedlings in >40 yr), so this changes little for its projection. |
| *A. beckleri* | No flowering or fruit in any stand, Oct 2018 (`flowering_fruiting_Oct2018`) | — | A single survey that may reflect phenology. Weak; not a target. |
| *M. nodosa* | 4 of 26 survivors flowering at Yarringully, Jan 2022 (notes A4, V5) | flowering | Secondary check only. Both hooks are off for *M. nodosa*. |

The *R. psidioides* and *R. rubescens* evidence separates the two stages. Tucki Tucki is about
fruit on trees that still flower. "No flowering since 2010" is about bud set. Rust in fruit is
direct infection of reproductive tissue. It is not the same process as canopy dieback, though
both are worse on badly affected trees.

## 3. What the model has now

- **Flowering:** `update_flowering()` (`R/simulate.R`). Each year, each eligible plant fails to
  flower with probability `dose_response(eff_susc, rust_flower_dose_response)`, where eff_susc =
  (1 − resist) × microsite × that year's pressure. It is a fresh Bernoulli draw each year, and it
  is off for *M. nodosa* (max_effect 0). Plants kept from flowering still count as mature
  (notes A4).
- **Fecundity:** `expected_recruits()` (`R/recruitment.R`). Each flowering plant contributes
  1 − `dose_response(eff_susc, rust_recruit_dose_response)` to the effective number of flowering
  plants fed to the recruitment curve. This is fruit maturation in all but name. It makes no
  draws and is off for *M. nodosa*.

**Gaps found when checking the code:**
1. **Site pressure is ignored.** Both hooks leave out the cell's `rust_pressure`, which the
   hazard and damage terms include. In `sdm` and `sdm_gradient` mode, a low-pressure site loses
   as much fruit as a high-pressure one. That is inconsistent with the rest of the model and with
   T18.
2. **Fruit failure does not affect who the parents are.** Fecundity weights shrink the number of
   recruits, but `sample_parents()` still draws mothers uniformly from all flowering plants. A
   susceptible plant whose fruit all fails is as likely to be a mother as a resistant one.
   Selection for resistance through fecundity, likely strong in these species, is therefore
   missing.
3. **Suppression is never persistent.** Flowering is redrawn every year, so a heavily affected
   plant flowers again in a light rust year. The code's own NOTE in `update_flowering()` flags
   this. Field reports ("no flowering since 2010") suggest it persists on damaged trees. The
   damage state D (notes T19, decision 5) was kept so this note could use it.

## 4. Proposal

All parts are optional and off by default. *M. nodosa* is unchanged draw for draw.

### 4.1 Site pressure in both hooks (fixes gap 1)
- Both hooks use `rust_pressure × eff_susc`, like the hazard and damage terms.
- **Applies to:** all species with the hooks switched on.
- **Not a change for *M. nodosa*:** its hooks are off, so nothing changes. The identity is tested.

### 4.2 Fruit failure weights the choice of mothers (fixes gap 2)
- Mothers are drawn from flowering plants with probability proportional to their fruit-success
  weight (1 − suppression), the same weight that sets the number of recruits.
- **Fathers:** drawn as now, from all flowering plants. Pollen comes from flowers, which exist,
  and the evidence is about fruit, not pollen.
- **Draws:** when the hook is off, all weights are 1 and the draw must stay identical
  (`sample.int` with uniform weights draws the same numbers only if no `prob` is passed, so pass
  `prob` only when the hook is on).
- **Applies to:** *R. rubescens*, *R. psidioides*.

### 4.3 Flowering driven by crown damage (fixes gap 3)
- A new optional setting, `rust_flower_driver`: `pressure` (the default, current behaviour) or
  `damage`.
  - **`damage`:** each year, the probability that a plant fails to flower is
    `dose_response(D, rust_flower_dose_response)`, where D is its crown damage. Damaged crowns
    stay damaged, so suppression persists until the crown recovers, and fires reset it through
    the resprout path.
  - Requires `rust_hazard_mode: damage` (the loader checks this).
- **Why D:** reproductive buds form on healthy growth, and a crown that is mostly dead does not
  set buds. The Olney treated arm shows crowns recover within one to two years once rust is
  removed, so persistence comes from the damage dynamics, not from a new parameter.
- **Fruit maturation stays on the year's pressure** (4.1). Fruit is infected directly in the
  season it develops (Tucki Tucki), whatever the state of the crown.
- **Applies to:** *R. psidioides* and *R. rubescens*, if they use damage mode.

### 4.4 Phenology
No sub-annual timing. Whether a species' fruiting season falls in the rust season (*R. rubescens*
fruit from Dec, warm and wet; *U. australis* Apr–Jul) is absorbed into each species' fruit
dose-response priors. The reason is recorded with the prior.

### 4.5 Calibration and reporting
- **Census:** two new columns, written only when either hook is on:
  - `share_flowering_suppressed`: the share of eligible plants kept from flowering by rust;
  - `fruit_success`: mean fruit-success weight over flowering plants.

  This keeps *M. nodosa* census tables unchanged.
- **Targets for Phase 5** (rule 4: nothing fitted by hand):
  - *R. psidioides*: fruit production ≈ 0 since rust arrived;
  - *R. rubescens*: fruit maturation ≈ 0 under high pressure (Tucki Tucki);
  - *U. australis*: little or no fruiting since 2018.

  How "≈ 0" becomes a tolerance belongs to `scripts/targets.R`.

## 5. IUCN question (not changed here)

IUCN Guidelines s4.3: mature individuals are those capable of reproduction, and individuals that
will never produce new recruits should not be counted. Under notes A4, a plant kept from flowering
by rust still counts as mature. Temporary suppression (the current per-year draw) fits that
reading. Persistent suppression (4.3), with no fruit for a decade, is closer to "will never
produce recruits".

**This note does not change the mature count** (IUCN framing: do not change without asking). It
proposes a third alternative count, `N_fruiting_capable`, alongside `N_mature` and `N_flowering`:
mature plants whose fruit-success weight is above a threshold. It is assessed over the same window
for comparison, as `N_flowering` is now. Whether any of these should become the main count is your
decision. It may matter a lot for *R. psidioides*.

## 6. Not in this note

- Pollinator loss or pollen limitation.
- Seed predation (*U. australis*): it belongs to the species' recruitment priors.
- A seed bank. All these species have short-lived, desiccation-intolerant seed.
- Effects of fruit failure on adult survival (reserve trade-offs).

## 7. Verification

- The full suite passes with `tests/reference/` untouched. With both hooks off and the driver
  absent, runs are identical draw for draw.
- New tests:
  - site pressure scales both hooks (and pressure 0 gives no suppression);
  - mothers follow fruit-success weights (a plant with weight 0 is never a mother but can be a
    father), and uniform weights give the same draws as now;
  - in `damage` mode, flowering suppression follows D and persists while D stays high;
  - `rust_flower_driver: damage` without damage mode is a loader error;
  - the census columns appear only when a hook is on.
- The decision is logged as a T entry, listing each decision required and the choice made.

## 8. Decisions needed

1. **Site pressure in both hooks** (*recommended*: matches the hazard and damage terms; no effect
   on *M. nodosa*). The alternative is leaving it out, so suppression would not vary with site.
2. **Mothers weighted by fruit success; fathers unweighted** (*recommended*). Alternatives: weight
   both, or neither (as now: no selection through fecundity).
3. **Optional flowering driver `damage`; fruit maturation stays on the year's pressure**
   (*recommended*). Alternatives: both stages on D; or a separate "reproductive suppression"
   state with its own recovery rate (one more parameter with no data to identify it).
4. **Phenology in the priors, not the time step** (*recommended*).
5. **Census columns and the alternative mature count `N_fruiting_capable`**, reported but not
   used for the main assessment (*recommended*). The main count stays as in notes A4 unless you
   decide otherwise.
6. **Order of work.** Build 4.1–4.3 and 4.5 now with toy tests, then set the per-species priors
   in Phase 4 (*recommended*). Alternatively, wait until the species files exist.
