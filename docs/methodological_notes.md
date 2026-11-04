# MnSim: methodological notes

A running record of the assumptions, decisions and known limitations behind
the *Melaleuca nodosa* IUCN assessment. The seed of the eventual methods
section. Add to it whenever an assumption is made or changed; give each entry
an ID so code comments and discussions can point at it, and note the date.

Status tags: **[decided]** a deliberate choice; **[assumption]** taken as true,
not tested; **[calibrated]** set from data; **[open]** known issue or gap.
Work still to do is tracked in `ROADMAP.md`.

---

## A. Assessment framing

**A1. IUCN criterion A4** [decided, Oct 2026]
Myrtle rust arrived in 2010 and the decline is ongoing, so the reduction spans
past and future. A2 (past only) would miss the projected decline; A3 (future
only, from the present) would miss the decline already under way.

**A2. One fixed window per run** [decided, Oct 2026]
Reference = mean number of mature individuals over the 11 years ending 2009
(the year before rust arrival). Comparison = mean over the 11 years ending
2009 + L, where L = max(10, 3 x T_g) rounded to whole years. Decline =
100 x (reference - comparison) / reference.
Not the worst window: searching all windows for the maximum decline picks
whichever window happens to end just after a fire, biasing declines upward
(by 7-10 points in a 40-run comparison, Oct 2026).
Details: the comparison end is capped at 100 years after the present (2026);
years after extinction count as zero; 11-year means smooth fire-driven
fluctuations.

**A3. Generation time from each run's own dynamics** [decided]
T_g = mean age of mothers of offspring born 1810-1910 (before rust),
excluding offspring of founding plants (whose ages are arbitrary) and of
mothers still alive at the end of the window. Each run's window length L
uses its own T_g, so parameter changes that alter generation time change
the window consistently. With the current calibration T_g is about 23-24 yr,
so L is about 70 yr.

**A4. Mature individuals** [decided, Oct 2026]
Counted as plants alive and of reproductive age (age >= age at first
flowering), including plants not flowering in a given year because they are
resprouting after fire or suppressed by rust. Basis: IUCN Guidelines v16 s4.3,
mature individuals are those *capable of* reproduction; temporarily
suppressed plants remain capable.
The stricter flowering-only count (`N_flowering`) is kept and assessed over
the same window in every run (`pct_decline_flowering`), for comparison.
Field illustration: at Yarringully in Jan 2022, 4 of 26 surviving trees were
flowering (Pegg et al. 2025); a flowering count would read that as an 85%
loss of mature individuals among survivors.

**A5. Reporting the ensemble** [decided, Oct 2026]
The main result is the distribution of runs across IUCN categories (LC/NT,
VU, EN, CR), not a single value. The median and quantiles of decline are
reported as the central tendency and spread.

**A6. Rust arrival** [assumption]
Rust is present from 2010 in every grid cell (first Australian detection,
NSW Central Coast, 2010). No spatial spread lag is modelled; the
distribution model governs how hard rust hits, not when it arrived.

## B. Range-wide ensemble

**B1. What a run is** [decided]
One run = one grid cell (0.1 degree) x one simulation seed (fire timing,
annual rust years, demographic chance) x one set of values for the uncertain
parameters. Each run has its own random-number stream, so adding runs never
changes earlier ones (`runner/range_ensemble.R`, which replaced
`spatial/03_range_wide_iucn.R` in Nov 2026; same cells, seeds and draws as
the October talk ensemble).

**B2. Cells drawn in proportion to abundance** [assumption]
Cells are drawn with probability proportional to `weight_n_raw`, the number of
occurrence records in the cell. Records reflect survey effort as well as
abundance. The record-density alternative (`weight_kde`) was recommended in
`02` but is not yet used. [open] Choose and justify one; ideally replace both
with an independent abundance estimate.

**B3. Parameter uncertainty** [decided, Nov 2026]
Currently one parameter is drawn per run: the extra rust mortality of
resprouting plants (`max_effect`), uniform on 0.32-0.78, the range whose Pegg
scenario medians lie in the 95% interval of the observation (V4). To be
replaced by fitted parameter sets (ROADMAP Phase 2).

**B4. Cells without inputs** [decided, Nov 2026]
Until Nov 2026, 16 cells (1.7% of abundance weight; 14 NSW, 2 Qld, all
coastal, up to 19 records each) had no SDM value: the SDM was read at the
cell centre, which fell at sea. The SDM is now averaged over the raster
cells inside each 0.1-degree cell, ignoring sea (`extract_sdm_cells()`,
`spatial/rust_pressure.R`), which gives every coastal cell a value from its
land. Cells still lacking inputs would be dropped and reported, not lost.

## C. Myrtle rust

**C1. Disease triangle** [assumption]
Each plant's effective susceptibility in a year = (1 - resistance score) x
microsite susceptibility x the year's pathogen pressure. Resistance is
genetic and inherited; microsite susceptibility is fixed at birth
(Beta(20, 5), the same distribution in every cell); pathogen pressure is one
draw per year shared by all plants in the cell (Beta, see C2).

**C2. Rust pressure anchored at one site** [calibrated + assumption, Nov 2026]
The absolute level of rust pressure is set at Yarringully NR, where field
incidence gives a mean annual pressure of 0.70 (`RUST_ANCHOR` in `params.R`,
the only place this value is set).
Each cell's mean pressure = 0.70 x (SDM in the cell / SDM in the anchor
cell), capped at 0.95; the annual draw is Beta(mean x 20, (1 - mean) x 20).
Two assumptions: (i) pressure is proportional to climatic suitability
(straight line through zero); (ii) Yarringully's field incidence sets the
level for the whole range.
Implementation: the anchor cell is the occupied cell nearest Yarringully
(cell 3623, 153.25 E, 29.15 S; the reserve's point location falls in an
unoccupied cell). Its SDM value is read from the raster every time, with the
same extraction as every other cell (mean over the cell's land since Nov
2026, B4), so it always gets exactly 0.70 whatever
SDM is supplied. A new SDM therefore changes only the relative pattern;
the level, and the Pegg calibration (V4), carry over
(`spatial/rust_pressure.R`, `spatial/02c_update_rust_pressure.R`).
Guards: a warning if the anchor cell's SDM differs from its surrounding
+/-0.15 degree by more than 10%; a count of cells at the 0.95 cap. With the
current SDM no cell exceeds the anchor (anchor SDM 0.993).
[open] Confirm and cite the field-incidence source for 0.70.

**C3. Year-to-year variability of rust pressure** [assumption, open]
The Beta concentration (alpha + beta = 20) sets how much rust pressure varies
between years. No data constrain it yet; it affects the spread among runs
(e.g. Pegg-scenario spread is wider at 20 than at 48). Candidate for fitting
(Phase 2), e.g. against year-to-year incidence at Yarringully.

**C4. No rust gradient along the coast: current SDM kept** [decided, Nov 2026]
Suitability is about 1 in almost every coastal cell, so rust pressure is
about 0.70 throughout the core range and falls only towards the arid
interior. Numbers (Nov 2026): 313 of 375 cells have suitability >= 0.9;
none exceeds the anchor cell (0.993); median pressure 0.70 in every latitude
band from 24 to 38 S. In the range ensemble, fire frequency, not rust,
explains nearly all differences between coastal cells.
Decision: keep the current SDM as the rust layer. It captures what the data
support -- a sharp decline in rust with aridity -- and treating the coast
as uniform reflects that intensity along the coast is not resolved by the
data (year-to-year variation, C3, is large).
Why the SDM saturates (`mrsdm/code/evaluate_sdm.R`, Nov 2026): pseudo-
absences are drawn 250-500 km from records, so the model contrasts the humid
coast with the dry interior; GBM and RF responses are flat across coastal
climate (e.g. max summer temperature below ~31 C, rainfall 1000-3000 mm).
Coastal climate does vary (annual mean 17-21 C, Spearman rho with latitude
0.95) but the occurrence model has no reason to pass that on. Other
diagnostics: continuous Boyce index 0.93 (P/E > 1 only above suitability
~0.83); random-CV TSS 0.95-0.99 is inflated by the distant pseudo-absences;
spatial-block CV TSS 0.47-0.60 (AUC 0.74-0.98), from only 2 of 4 folds
because the others lacked presences or absences; clamping at 3 of 358
M. nodosa cells (current), 1 (future). GLM response curves show near-
separation (step and U shapes), an SDM-level caveat, not one that affects
MnSim given the saturation.
Bias probe (`mrsdm/code/probe_sdm_bias.R`, Nov 2026): refitting with
background drawn at random within the study area (as in KarinaGuo/MQ_SDM)
does give a coastal pattern (suitability 0.64-0.92, 10th-90th percentile),
but across the 313 coastal cells it follows survey effort, not climate:
Spearman rho with log distance to the nearest city -0.40, with local rust
record density 0.55, with annual mean temperature 0.16, with latitude 0.07.
In a joint regression on coastal cells, distance to city explains R2 0.16
and temperature 0.03 (coefficients -0.30 and 0.14, standardised). So the
gradient a presence-background model produces here is largely where people
look; it is not used.
[assumption] Rust pressure is uniform along the coast at the anchor level.
[open] If incidence data across several coastal sites become available,
test this directly; a target-group background (TG_FILE in the probe) would
be a stronger effort control. A coastal gradient could also enter Phase 3
as an uncertain term (strength drawn per run, including zero) to see whether
it would matter for the listing.

**C5. Rust effects in the model** [assumption]
Rust raises death risk, most for seedlings and young plants (declining with
age), with an extra component for resprouting plants (`max_effect`, C/V4);
it lengthens the time to recover from resprouting. Effects on flowering and
seed set are implemented but switched off (no data).
Resprout delay (Nov 2026): the extra years a burnt plant spends resprouting
(Poisson mean) = 1.79 x (1 - resistance) x microsite x the site's MEAN annual
pressure, i.e. through the disease triangle, using the site mean rather than
the fire year's draw because the delay spans several years. 1.79 = 1/(0.8 x
0.70) keeps the previous mean of one extra year for a fully susceptible plant
at the anchor site; the delay now also varies between plants (microsite) and
between cells (rust pressure). [assumption] Magnitude not constrained by data.

## D. Fire

**D1. Fire records and processing** [decided, Nov 2026]
Two state records, put on one footing (`spatial/fire/`, developed in a
sandpit and moved into MnSim in Nov 2026):
- *NSW*: NPWS Fire History (agency fire perimeters, labelled wildfire or
  prescribed burn; year = first year of the season label, so 2019-20 counts
  as 2019). Used 1987-2024 (1970-2024 kept as a sensitivity check).
- *Queensland*: QRSC fire scars (Landsat/Sentinel-2 burnt-area mapping,
  annual composites, fire type not recorded), downloaded from the Qld
  government ArcGIS service for the M. nodosa area; and the QPWS fire
  history (planned burns and wildfires on QPWS and state land), used to
  separate planned burns. Used 2014-2024 (D4).
Each cell's burnt fraction in a year = area of the union of that year's
burnt polygons (overlapping fires not double counted) / the area of the
cell's land within its own state, in Australian Albers (EPSG:3577). A cell
belongs to the state containing its centre (coastal cells: nearest state).
A cell-year is a *fire year* when >= 10% of the cell's land burns.
Changes from the pre-Nov-2026 `02`: union instead of summing a year's
polygons; land base instead of the whole 0.1-degree cell (sea included);
records start 1987 instead of 1970. With the old conventions the new code
reproduced `02`'s NSW table exactly (fire-year counts 255/256 cells
identical; burnt fractions within 0.03).
Outputs: `spatial/data/fire/cell_years_<state>_<type>.csv` (state nsw/qld;
type wild/planned/all), built by `spatial/fire/02_fire_cell_years.R`.

**D2. Fire frequency** [calibrated, Nov 2026]
Annual probability of a fire year per cell = posterior mean
(alpha + k) / (beta + T) for k fire years in T years observed, under a
Gamma(alpha, beta) prior on the cell's rate (Gamma-Poisson), clamped to
[0.005, 0.5]. T is the cell's own record length (38 yr NSW, 11 yr Qld).
The prior is fitted by maximum likelihood across all cells (empirical
Bayes: k is negative binomial given alpha, beta and T). This replaces
`02`'s method of moments on k/T over cells with k > 0, which ignored the
cells without fire and treated chance variation in short records as real
differences between cells.
One prior is shared by both states: fitted separately they are almost the
same (wildfire: NSW Gamma(0.68, 16.3), Qld Gamma(0.74, 19.8); likelihood-
ratio test p = 0.85). Shared prior Gamma(0.69, 16.7): mean 24 yr between
fire years, worth about 17 years of observation, so the 11-yr Qld records
lean on it (k = 0, 1, 2 in 2014-2024 -> about 40, 16, 10 yr between fires).
Burnt fraction per fire: mean and sd over the cell's fire years
(`fire_p_fimp`, `sd_frac_burned`; cells without a fire year get the median
mean burnt fraction of cells with one; sd from cells with >= 2 fire years,
else the range-wide value in `params.R`).
Fire occurrence in a run: a random draw each year with that probability.
Values (wildfire, Nov 2026): NSW 4.2 and Qld 3.7 fire years per 100
cell-years; median interval 32 yr (NSW), 40 yr (Qld); median burnt fraction
in a fire year 0.27 (NSW), 0.19 (Qld).
[assumption] Fire years are a Poisson process with a constant rate per cell
(no trend, no dependence on time since fire); 10% threshold.

**D3. Fire extent and severity** [decided + assumption, Nov 2026]
Scale: a simulated population is a stand or a few stands (~10,000 plants),
far smaller than the 0.1-degree cell (~10,000 ha) that the fire history
describes, and smaller than most fires counted (>= 10% of the cell).
Each fire in the cell has its own size f (burnt fraction of the cell),
drawn from a Beta with the cell's mean burnt fraction and between-fire sd
(`sd_frac_burned`; cells with fewer than two fire years use the range-wide
median, 0.146 in `params.R`: 105 cells with >= 2 wildfire years, NSW + Qld,
Nov 2026; previously 0.125 from the old NSW table). Then:
- *Extent*: the fire reaches the stand with probability f; if it does, each
  plant burns with probability 0.9 (`fire_stand_cover`, fixed: unburnt
  patches within the stand). Most fires miss; those that arrive burn most of
  the stand.
- *Intensity*: linked linearly to that fire's size. Burnt plants die with an
  age-dependent probability: seedlings almost always, older plants at a
  ceiling of 0.30 x f. Big fires are therefore both likelier to reach the
  stand and more lethal (extent and intensity correlated between fires,
  through one parameter).
- Survivors resprout and do not flower for a few years (Poisson, mean 3,
  plus a rust-dependent delay).
Over many fires the mean share of plants burnt equals the cell's mean burnt
fraction, as before; the change is that it comes as occasional whole-stand
fires rather than a scattered quarter of the stand every time, and that a
fire's lethality no longer depends on the cell's long-run average.
The Pegg scenario is the limiting case: f = 1, full cover, ceiling 0.30
(unchanged, so the V4 calibration stands; Pegg median 0.474 after the change).
The previous behaviour (every fire reaches the population; per-plant
probability and ceiling both set by the cell mean) remains available as
`fire_extent_mode = "plant"`.
Effect (60 paired range-wide runs, same cells/seeds/draws, Nov 2026): fires
reaching the stand fell from ~22 to ~5 per run (1600-2210); median decline
48.3% vs 52.7% (paired difference median -2 points), but the spread widened
(sd 12.5 vs 7.3 points; 95th percentile 78% vs 66%), giving VU 60% / EN 35% /
CR 5% vs VU 38% / EN 62% / CR 0%. Whether a large fire reaches the stand
during the rust period now drives much of the between-run variation. With
the median near the 50% EN threshold, category shares are sensitive to small
shifts; read them with that in mind.
[assumption] Linear size-intensity link; fixed within-stand cover 0.9; fire
sizes at a point follow the cell's distribution of burnt fractions.
[open] Severity data (e.g. NSW fire extent and severity mapping) could
replace the linear link; the ceiling 0.30 is not directly calibrated
(the share of trees killed outright at Yarringully is not reported).

**D4. Queensland fire record: edited years only** [decided, Nov 2026]
Only 2014-2024 of the Qld fire scars are used. Evidence (border check,
`spatial/fire/03_fire_checks.R`): over 1987-2024, Qld cells within 1 degree
of the border had 4.96 times the NSW fire-year rate (95% CI 3.43-7.41).
Year by year, 27-44% of Qld cells had a fire year in every year 1987-2013,
dropping to 0-8% from 2014 (NSW level). 31 of 119 Qld cells were marked burnt
in all 26 years 1987-2012 (and 94% of them in 2013) at burnt fractions just
above the threshold (typically 0.11-0.25), yet have fewer fire years after
2014 than other cells: repeated misclassification of the same patches in the
automated mapping (QRSC: automated before ~2013, ~80% detection and up to 30%
false scars; edited since; Sentinel-2 from 2017). 2013 still behaves like
the automated years (starting there raises the Qld rate by ~60%; test of a
shared prior p = 2.5e-5).
With 2014-2024 the records agree at the border: Qld/NSW rate 0.96
(0.52-1.72) over 2013-2024, 0.55 (0.25-1.12) in the main scenario (NSW
1987-2024). Over the common period 2014-2024 alone the ratio is 0.33
(0.14-0.76), driven by 2019 (43% of NSW border cells burnt vs 4% in Qld);
without 2019 the rates are similar (NSW ~0.03, Qld 0.021 per cell-year).
Read as one real fire season (2019-20) on the NSW side, and a reason to keep
NSW's longer record rather than cut both states to 2014-2024.
[assumption] Fire regimes are stable enough that 11 recent years describe
Qld cells (with the shared prior, D2). [open] A longer, consistent Qld record
(e.g. reprocessed Landsat burnt area) would help.

**D5. Wildfire as the default; planned burns kept separately** [decided, Nov 2026]
The default fire regime is wildfire only: NSW wildfires; Qld fire scars
minus QPWS planned burns. Planned burns matter (they are 56% of Qld fire
years at M. nodosa cells in 2014-2024, 14% in NSW over 1987-2024), but their
records are less consistent: separate agency programmes; QPWS planned burns
are mapped burn boundaries (which can include unburnt patches) and cover
QPWS and state land only, so Qld "wildfire" still contains planned and
agricultural burning elsewhere (the scars do not record fire type). The
planned share in Qld is therefore a lower bound.
All three regimes are written to `mn_site_params.csv` for every cell:
wildfire under the standard names (`fire_prob_annual`, `fire_p_fimp`,
`sd_frac_burned`, ...), all fire and planned burns with suffixes `_all` and
`_planned`. Runs choose with `fire_regime` ("wild" default, "all") and can
scale the fire probability with `fire_scale` (e.g. 0.5 = suppression
halving fire frequency) (`build_site_params()`, `run_range_ensemble()`).
Scenario ensembles reuse the default's cells, seeds and parameter draws, so
runs pair up.
Planned-burn estimates: one prior per state (agency practice differs: test
of a shared prior p < 0.001) and no rate floor (most cells are never burnt
by an agency; NSW Gamma(0.07, 10.6), Qld Gamma(0.14, 4.9)). Planned burns
burn smaller fractions than wildfire in NSW (median 0.19 vs 0.27); similar
in Qld (0.19). Wildfire + planned rates slightly exceed the all-fire rate
(both clear the threshold in some cell-years), so "all" is estimated
directly rather than summed.
NSW record types (NPWS polygons 1970-2024 within the bounding box of the NSW
cells): 18,471 wildfire, 12,084 prescribed burn, none unlabelled.
[open] Management scenarios with planned burning as its own process (e.g.
burning a target site every 8 or 15 years on top of wildfire) need a second
fire process in the model, with its own extent and possibly a lower
severity; not yet implemented (ROADMAP).

## E. Demography and genetics (current values, `params.R`)

**E1.** Background mortality is a constant annual hazard (Weibull k = 1,
scale 60 yr), plus a senescence term rising with age and a juvenile
mortality term that increases with canopy density. [assumption]
The Weibull hazard is converted to an annual probability by integrating it
over the plant's year of age: p = 1 - exp(-(H(a) - H(a - 1))),
H(x) = (x / scale)^k (Nov 2026; previously the hazard itself was used as the
probability, a small overestimate: 0.0167 vs 0.0165 per year at scale 60).
**E2.** Recruitment: Beverton-Holt seedling production from flowering plants
(R_max 1000, K_half 4000), parents sampled from flowering plants; selfing rate
0.10; mean age at first flowering 5 yr. [assumption]
Recruits enter at age 0 after the year's mortality step and first face
mortality at age 1, so there is no explicit first-year seedling mortality:
recruitment numbers are interpreted as seedlings that survive their first
year (Nov 2026, decided; adding age-0 mortality would require recalibrating
R_max and K_half, for which there are no data either).
**E2b.** Bookkeeping (Nov 2026): every plant's birth_year is a calendar year
(age = year - birth_year at year end); founding plants carry founder = TRUE
and are excluded as mothers from the generation-time estimate (A3). The
census counts all deaths, with fire deaths also given separately
(deaths_fire); previously fire deaths were missing from the yearly total.
Dead plants are archived during a run for speed and restored at the end;
this changes no result (verified identical, Nov 2026).
**E3.** Resistance: multilocus, Mendelian inheritance; currently starting
frequency 0 (fully susceptible) everywhere, pending the range-wide resistance
data (ROADMAP Phase 1c). [open]
**E4.** Simulation phasing: spin-up from 1600; generation time measured
1810-1910; rust from 2010; present 2026. [decided]

## V. Validation against Pegg et al. (2025)

**V1. The observation** Yarringully NR: 50 fire-affected trees, all
resprouting when monitoring began (March 2020); 24 (48%) dead at the
January 2022 follow-up; 95% binomial interval 0.34-0.63. Monthly values in
`talk/data/` are read off their Fig. 3a (approximate); only 24/50 is exact.

**V2. The model scenario** [decided, Oct-Nov 2026]
Spin-up 1910, rust from 2010, one fire in 2020 burning every tree, run at the
Yarringully anchor cell's rust pressure (C2). The cohort is plants alive at
the end of 2019 that the fire did not kill; deaths are classed by cause, so
rust or background deaths in the fire year count as post-fire deaths. The
January 2022 survey = the end of model year 2021 (two post-fire years).

**V3. What is fitted** [decided, Oct 2026]
Only the two-year endpoint. In the field, deaths accelerated in the second
year (about 10% in year 1, then about 42% of the remainder); the model's
yearly risk is roughly constant (about a third each year). Matching that
shape would mean fitting a ramp to a single 50-tree series, judged
overfitting; being in the right place after two years is the target.

**V4. Calibrated resprout rust mortality** [calibrated, Nov 2026]
`max_effect` = 0.53 reproduces a median of 48% dead (40 seeds per value,
Yarringully cell). Medians span the observation's 95% interval for values
0.32-0.78, which is the range drawn in the ensemble (B3).
History: 0.55 (Oct) -> 0.45 (tuned while the Pegg runs wrongly used rust
pressure 0.83) -> 0.53 (Nov, at the anchor cell's 0.70).

**V5. Secondary check** About 15% of surviving trees flowering in January
2022 (4 of 26) against about 13% of surviving resprouters recovered (able to
flower) in the model. Not used for fitting.

## Decision log

| Date | ID | Decision |
|---|---|---|
| Oct 2026 | A1 | Criterion A4 (switched from A3) |
| Oct 2026 | A2 | One fixed window per run; no worst-window search |
| Oct 2026 | A4 | Count all mature individuals; keep flowering-only alongside |
| Oct 2026 | A5 | Report category distribution across runs |
| Oct 2026 | V2 | Pegg: deaths by cause; survey = end of model year 2021 |
| Oct 2026 | V3 | Fit the two-year endpoint only |
| Nov 2026 | C2 | Single rust anchor; SDM normalised at the anchor cell |
| Nov 2026 | V4 | `max_effect` 0.53; ensemble range 0.32-0.78 |
| Nov 2026 | D3 | Fire: per-fire size; stand reached with P = size; cover 0.9; kill ceiling linear in size |
| Nov 2026 | C5, E1, E2 | Resprout delay via the disease triangle; Weibull hazard -> annual probability; recruits = first-year survivors |
| Nov 2026 | D1 | NSW + Qld fire records on one footing: union of a year's fires, land base, fire year >= 10% of land; NSW 1987-2024 |
| Nov 2026 | D2 | Fire-rate prior by maximum likelihood over all cells, shared by both states |
| Nov 2026 | D4 | Qld fire scars 2014-2024 only; automated 1987-2013 mapping rejected (border check) |
| Nov 2026 | D5 | Wildfire-only default; all fire and planned burns kept per cell for comparison and management scenarios |
| Nov 2026 | B4, C2 | SDM averaged over each cell's land, not read at the centre: fills the 16 coastal cells |
| Nov 2026 | C4 | Current SDM kept; rust uniform along the coast; random-background gradient shown to follow survey effort |
