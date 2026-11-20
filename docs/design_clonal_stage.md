# Design note: clonal stage (genets, ramets and root suckers)

Status: **approved 2026-11-21 with the recommendations in section 6; built** (notes T17, which also
records two refinements found in smoke runs: suckers saturate with the recruitment curve, and only
established genets keep a bare rootstock).

## 1. Purpose

Three species persist or spread by root suckers. The model cannot represent that now.

| Species | What the data say | What the model must do |
|---|---|---|
| *U. australis* | Reproduction is by suckers only: no wild seedling in >40 yr, and seed viability is ~20 weeks. There are ~300 genets and ~10,000 stems before rust. Stems live >1,500 yr. The 2019–20 fires killed stems by size class (100/80/52% at Nightcap), but at least 5 genets were lost. IUCN counts the ramet as the mature individual. | Make new ramets from existing genets. Kill stems by size in fire, while most genets survive. Count ramets. |
| *R. psidioides* | 61% of 66 populations were reduced to root suckers by 2018. Suckers are infected at ~30 cm, and 60% of 80+ suckers were in severe decline. | When rust kills the canopy, leave a rootstock that keeps suckering. Its suckers die young ("living dead"), until the rootstock runs out. |
| *R. rubescens* | All seedlings and suckers are killed by rust. 2 of 3 suspected clonal clumps persist. | The same machinery, at low sucker rates. |

**Constraint.** M. nodosa is unchanged, bit for bit. Clonal mode is an optional species setting
(T12 §2.1). When it is off, no new random numbers are drawn, no columns are added to the census,
and the params list and run IDs stay as they are.

## 2. What exists now

- **Individuals.** One row per plant in `pop`, and one row per plant in the genotype matrix
  `resist_gt`. Every plant comes from seed (`recruit()`): two parents from the flowering pool,
  Mendelian inheritance, and a fresh microsite draw (`env_susceptibility`).
- **Fire.** `apply_fire()` kills a plant or sets it resprouting. Kill risk falls with age
  (`fire_kill_half_sat`, `fire_kill_hill`) towards a ceiling set by fire size. Resprouting is the
  same plant recovering its canopy, not a new stem.
- **Rust.** The rust hazard decays with age from a peak (`rust_dose_response$age_peak`), so young
  plants are hit hardest. That already gives "suckers die young" once suckers exist.
- **Mature individuals.** `N_IUCN` counts living plants with age ≥ `age_first_flower`.
- **Generation length.** `compute_gen_time()` is the mean age of mothers of offspring born in
  the gen-time window. That matches the IUCN definition (average age of parents of the current
  cohort), but it counts seed only. With no seedlings, T_g is NA and the A4 window cannot be set.

## 3. Proposal

### 3.1 Structure: ramet rows grouped into genets

- **Ramets.** Each row of `pop` stays one individual, now read as a **ramet** (stem), with a new
  column `genet_id`.
  - Ramets of one genet share its genotype row, copied into `resist_gt`, so `genetics.R` is
    unchanged.
  - Census counts, mortality, fire and flowering all stay per ramet, which is the IUCN mature
    individual.
- **Genets.** A small `genets` table (`genet_id`, `rootstock_alive`, `years_without_canopy`)
  carries what belongs to the genet rather than to a stem: the rootstock.
  - It is returned with the run's state, so `init_state` continuation still works.
- **Seed recruits** start their own genet (`genet_id = id`), so sexual and clonal reproduction
  coexist.
- **Founders.** The `N0` founding ramets are grouped into `clonal_founder_genets` genets. Each
  genet takes one genotype and one microsite.
- **Why rows rather than a stem count per genet.** Fire kill by stem size (the U. australis
  quadrats) and maturity by stem size need per-stem age, and the existing per-plant code already
  does both. A count model would mean rewriting fire, rust and census.

### 3.2 Suckering (new step after recruitment)

Each year, a living genet with its rootstock alive produces suckers:

  E[suckers] = (sucker_rate × n_canopy_ramets + sucker_pulse × n_ramets_killed_last_year)
               × (1 − shade_suppression)

- `n_canopy_ramets` counts the genet's mature ramets.
- **Disturbance pulse.** Ramets killed by fire or rust in the previous year trigger suckering.
  This is the observed response to fire and slashing in U. australis, and to canopy death in
  R. psidioides.
- **Rootstock floor.** A genet with no canopy still produces suckers at `sucker_rate_rootstock`.
  This is the R. psidioides sucker pool.
- **Density.** Suckers are limited by the same canopy shading as seedlings (`shade_suppression`),
  so seed and suckers compete for the same space.
- **New ramets** enter at age 0 with the genet's genotype. They take the parent ramet's microsite
  (`env_susceptibility`), because a sucker grows beside its parent; a fresh draw is the
  alternative (decision 3). Their `mother_id` is a parent ramet sampled from the genet's ramets,
  weighted towards mature ones. A new `origin` column records "founder", "seed" or "sucker".
- **No new rust code.** Suckers face the existing hazards. Rust's age peak kills them young, as
  the R. psidioides and R. rubescens observations require.

### 3.3 Rootstock survival (genet death)

A genet dies when its rootstock dies. The rootstock is at risk in two ways:
- **When fire kills all its ramets:** it survives with P = `rootstock_fire_survival`. For
  U. australis this is a prior. At least 5 of the 34 burnt monitored genets were lost (Nightcap), so at most about 0.85.
- **Each year without a mature ramet:** an extra hazard `rootstock_hazard_no_canopy`, representing
  reserves running down while suckers keep dying.
  - This sets how long a "living dead" sucker pool lasts.
  - R. psidioides populations reduced to suckers are the calibration target.

Ramets otherwise die as now. The genet survives while any ramet or the rootstock is alive.

### 3.4 Census and assessment

- **Census.** `N_IUCN` is unchanged: mature ramets. Clonal mode adds `N_genets` (living genets)
  and `N_rootstock_only` (genets with no living ramet). Like flowering-only, these are reported
  alongside, not assessed. The columns appear only when clonal mode is on.
- **Generation length.** T_g becomes the mean age of the parents of *all* new individuals in the
  gen-time window, seed and sucker. For a sucker, the parent is the ramet chosen in §3.2. This
  keeps the IUCN definition, with the ramet as the individual. For species with seed only, it is
  the current calculation exactly. For U. australis, T_g then depends on how old a ramet is when
  it suckers, which is unmeasured. Hence decision 4.
- **The 100-yr cap** applies as now. Long T_g simply hits it.

### 3.5 Settings (all optional; schema `required_when: {clonal_mode: genet}`)

| Setting | Layer | Meaning |
|---|---|---|
| `clonal_mode` | species | `none` (default) or `genet` |
| `clonal_founder_genets` | species | number of founding genets (U. australis ~300) |
| `sucker_rate` | species | suckers per mature ramet per year |
| `sucker_rate_rootstock` | species | suckers per year from a genet with no mature ramet |
| `sucker_pulse` | species | suckers per ramet killed in the previous year |
| `rootstock_fire_survival` | species | P a genet's rootstock survives a fire that kills all its ramets |
| `rootstock_hazard_no_canopy` | species | annual rootstock death risk while the genet has no mature ramet |
| `sucker_microsite` | species | `parent` (default) or `new` |

All are biology, so they go in the species YAML under rule 5: each cites its CSV row or is a prior
with a range and reason. Most will be priors, then fitted in Phase 5.

## 4. Not in this note

- **Explicit space.** There are no coordinates. Genets compete only through shared canopy shading,
  so the spread of a clone across a gap is not modelled.
- **Size as a state.** Stem age stands in for size, as it does now for fire kill and maturity.
  If the U. australis quadrats need stem height or DBH explicitly, that is a later note.
- **Rust damage that accumulates in the rootstock** (the PLAN "accumulating-damage hazard"). That
  is a separate note. `rootstock_hazard_no_canopy` is its simple stand-in here.
- **Somatic mutation** and any genetic change within a genet.

## 5. Verification

- The full suite passes with `tests/reference/` untouched. Clonal mode is absent for M. nodosa.
- New tests:
  - with clonal mode off, a run is identical to the current code, draw for draw;
  - suckers carry their genet's genotype and genet id, and get `origin = "sucker"`;
  - a fire that kills every ramet ends the genet with P = 1 − `rootstock_fire_survival`;
  - a genet with no mature ramet and `rootstock_hazard_no_canopy = 0` keeps producing suckers
    indefinitely ("living dead"); with a hazard, the pool lasts about 1 / hazard years;
  - the disturbance pulse raises suckering in the year after a fire;
  - T_g counts sucker parents, and equals the current value when there are no suckers;
  - `N_genets` and `N_rootstock_only` appear only in clonal mode.
- A short smoke run per species shape (clonal-only, mixed seed and sucker) checks that the outputs
  are plausible.
- The decision is logged as a T entry once approved.

## 6. Decisions needed

1. **Ramet rows grouped into genets** (§3.1). *Recommended.* The alternative is genet rows with a
   stem count, which is cheaper but loses stem-size fire kill and per-stem maturity.
2. **Rootstock as a genet state** with fire survival and a no-canopy hazard (§3.3).
   *Recommended*, since both R. psidioides ("reduced to suckers") and U. australis (genets
   survive fire that kills their stems) need it. The alternative is no rootstock: a genet dies
   with its last ramet.
3. **Suckers keep the parent's microsite** (*recommended*: they grow beside it), or draw a new one.
   Keeping it means a genet in a bad microsite stays in it.
4. **Generation length for clonal species.** Either (a) mean parent age over seed and sucker
   offspring, emergent as now (*recommended*, the IUCN definition with the ramet as the
   individual), or (b) a fixed `generation_length` from the species file with a prior and reason.
   (b) is simpler but leaves the window to an unmeasured number. Under (a) that number is still
   unmeasured, but it follows from the sucker parameters that are fitted.
5. **The order of work.** Build this before the U. australis and R. psidioides YAML files in
   Phase 4 (*recommended*), since neither species can run without it.
