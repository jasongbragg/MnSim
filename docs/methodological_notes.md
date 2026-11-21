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
unoccupied cell). (Since T18 this cell is fixed in `config/landscape.yaml`,
as `cell_lon`/`cell_lat`, for all species; its master-grid id is 100933.) Its SDM value is read from the raster every time, with the
same extraction as every other cell (mean over the cell's land since Nov
2026, B4), so it always gets exactly 0.70 whatever
SDM is supplied. A new SDM therefore changes only the relative pattern;
the level, and the Pegg calibration (V4), carry over
(`spatial/rust_pressure.R`, `spatial/02c_update_rust_pressure.R`).
Guards: a warning if the anchor cell's SDM differs from its surrounding
+/-0.15 degree by more than 10%; a count of cells at the 0.95 cap. With the
current SDM no cell exceeds the anchor (anchor SDM 0.993).
[open] Confirm and cite the field-incidence source for 0.70.

**C3. Year-to-year variability of rust pressure** [assumption, open; being fitted, V6]
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
type wild/planned/all), built by `spatial/fire/02_fire_cell_years.R`. Since T13 the tables are
shared by all species; see T14 for an error in the NSW tables.

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
- Survivors resprout: they do not flower, and are exposed to the extra rust
  mortality of resprouting plants, for a fixed 2 years from fire alone, plus
  a rust-dependent Poisson delay (decided Nov 2026; previously Poisson with
  mean 3, at least 1 year; `resprout_delay_dist = "poisson"` restores it).
  [assumption] Two years is a judgement, not a measured value; because it
  sets how long resprouters are exposed to rust mortality, it is fixed
  before the Phase 2 mortality fit (V6).
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
**E3.** Resistance: Mendelian inheritance; currently one locus
(`resist_locus_effect = 1`, `resist_dominance = 1`) at starting frequency 0
(fully susceptible) everywhere. [open, data pending]
Status (Nov 2026): range-wide seedling screening for resistance is under
way (unpublished); resistance appears rare and spatially clustered, and a
single-locus model is consistent with the work so far. Site samples are
small, so most site-level estimates will be imprecise (often zero).
Implications: the single-locus
assumption is reasonable; setting resistance to 0 is a close
approximation for most of the range, so including the data is expected to
change results modestly. Planned use (ROADMAP 1c): site prevalence ->
allele frequency (given dominance) -> per-cell `resist_freq0` through the
`RESISTANCE HOOK` in `build_site_params()`, drawn from its binomial
uncertainty. Until then, keep resistance confined to `params.R` and that
hook (no other hard-coding of resistance = 0).
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

**V6. Several plots: level and spread, by interval** [in progress, Nov 2026]
Why: Yarringully alone fixes the level of resprout rust mortality but not
how much sites differ. Pegg et al. (2025) followed seven M. nodosa plots of
50 tagged resprouting trees each, all at about the same modelled rust
pressure (0.70), with very different outcomes: Yarringully 48% dead (24/50,
to Jan 2022), Bundjalung 0% (to Jan 2022, despite heavy infection and
dieback), and five mid-north coast plots (the paper's "Central Coast";
tagged Sep 2020, final assessment May 2021) Kiwarrak 56%, Wallabi Point 20%,
Knappinghat 10%, Saltwater and McClymont Creek 0% (inferred from "deaths at
three of the five"). Fitting the spread asks how much between-site
variation the model needs and of what kind (C3).
Data: `runner/data/pegg2025_mnodosa_plots.csv`, one row per observation
interval (months since fire, trees at risk, deaths). INTERIM: coordinates
approximate except Yarringully; fire months approximate; the Yarringully
split between post-fire years read off Fig. 3a; two zero plots inferred.
Per-tree, per-visit records, coordinates, fire dates and infection records
requested from G. Pegg (Nov 2026).
Observation model (`runner/pegg_sites.R`): cohort = plants of reproductive
age alive before the fire and not killed by it (tagged trees are adults;
V2 used every resprouter); the Pegg scenario at each plot's cell gives the
cohort's survival after each post-fire year; a constant death rate within a
year maps that to any observed interval. Fitting is by interval, not two-
year totals, because the windows differ (the mid-north coast plots start 10
months after the fire), and whether deaths come in year 1 or year 2 matters.
Likelihood per plot = mean over simulated sites (seeds) of the binomial
probability of its deaths in each interval, so a parameter set that never
produces a near-zero site gives the zero-death plots almost no probability.
Two sources of between-site variation are compared: year-to-year (Beta
concentration of the annual rust draw, C3) and a persistent site effect
(`rust_site_sd`: logit-normal shift in a site's mean annual pressure, drawn
once per site; 0 by default, already carried through the range ensemble).
Interim result (Nov 2026; 20 seeds per cell; resprout window 2 yr, D3):
log-likelihood (7 plots), best max_effect per structure --
current model (conc 20, no site effect) -70.5 (max_effect 0.3);
year-to-year only: conc 5 -60.5, conc 2 -50.6 (0.6);
persistent site effect sd 1 -48.9 (0.8), sd 2 -34.0 (0.8, best, at the edge
of both ranges tried). The current model predicts every plot in a narrow
band (e.g. Bundjalung 26-50% dead for a 50-tree plot) and never zero deaths.
Even at the best set the chance of a 50-tree plot with no deaths is only
0-3%: adult resprouters die at ~4-5% a year in the model even without the
extra resprout rust mortality, near the most that zero deaths in ~20 months
allows. Single-site check: at the Yarringully cell (conc 20) the new window
and adult cohort give ~41% (max_effect 0.45) and ~50% (0.6) dead over its
record, so the V4 level is roughly unchanged.
Reading (provisional): the plots favour a large persistent site effect over
year-to-year variation, with a higher max_effect to keep the level. If
confirmed this matters for the assessment: a persistent effect implies
lasting low-rust sites (refugia) across the range. Caveats: interim data,
few seeds, optimum at the grid edge. Bundjalung's heavy infection without
deaths suggests the site effect may act on tolerance (mortality given
infection) rather than on pressure; the requested infection records can
separate the two. Next: the full grid (`runscripts/pegg.sites.fit.run.code`)
with wider ranges and 60 seeds; then the per-tree data.
[open] Background mortality of adult resprouters (~4-5%/yr) may be too high
to allow zero-death plots; revisit once the data are in.

## T. Testing (multispecies branch)

**T1. Regression references for M. nodosa** [decided, 2026-11-19]
Six fixed-seed scenarios (`plant_ibm/tests/regression/scenarios.R`):
- no fire;
- fixed fire years in stand mode;
- rust at the anchor pressure, Beta(14, 6);
- the same at a low-rust cell (7006, mean 0.148);
- the Pegg scenario;
- a shortened A4 run with counterfactual.

How the references work:
- The scenarios are built from `params.R`. The outputs are frozen in
  `tests/reference/` by `tests/regression/make_reference.R`, which never
  overwrites.
- Each reference stores the exact params list it ran with. Re-running the
  frozen params tests the code. Rebuilding the params from today's
  `params.R` tests whether the calibration has changed. A ROADMAP
  recalibration therefore shows up as a params change, not a code failure.

Test scale [assumption]: `params.R` settles at ~13,000 plants, but test runs
are capped at 500 (CLAUDE.md). The density-dependent parameters are
multiplied by 0.03, giving ~400 plants:
- `R_max` and `K_half` (counts of flowering plants);
- `half_sat` of juvenile decline and of shading (counts of canopy plants).

Per-capita dynamics stay close to the full-size model: the Pegg scenario
gives 52% of resprouters dead (observed 48%), and the A4 run gives a 52%
decline (EN), with T_g 22.7 yr. These runs check the code, not the
assessment.

Exact matching holds for R 4.3.1 with the default RNG. Any change to the
order of random draws breaks it (`docs/architecture.md` §3.2).

**T2. Second check: 20-seed summary statistics** [decided, 2026-11-19]
`plant_ibm/tests/regression/stats.R`; baseline frozen in
`tests/reference/summary_stats.rds` by `make_reference_stats.R`. Three
statistics:
- the annual death rate of mature plants without fire (1950-1999);
- the Pegg fraction of resprouters dead;
- the annual death rate under anchor rust with fire off (2011-2060).

Each takes the mean over 20 seeds. The test passes if the new mean lies
within 4 x sqrt(2) x sd/sqrt(20) of the baseline.

Purpose: a refactor that only reorders random draws fails the exact match
but should pass here.

Tested by deliberately changing a scratch copy of the code:

| Change | Exact-match tests | Statistics |
|---|---|---|
| Extra unused draw each year | fail | pass |
| Rust mortality x 1.2 | fail | `rust_mortality` fails |
| Senescence x 1.5 | fail | `mature_mortality` fails |
| Background lambda 60 -> 55 | — | pass |

The lambda 60 -> 55 change shifts mature mortality by ~2.5%, while the
tolerance is ~5%.

Two candidates were dropped because they hardly respond to real changes:
the all-plant death rate (mostly juveniles, offset by recruitment) and the
2009-2060 decline (a 20% rise in rust mortality left it unchanged). This is
a coarse guard; the exact-match tests are the main check.

**T3. Procedure for updating the references** [decided, 2026-11-19]
Full procedure: `docs/reference_updates.md`.
- Which test fails says what changed, and whether an update is allowed:
  - recalibration;
  - intended behaviour change;
  - draw-order refactor;
  - R upgrade.
  An unintended behaviour change is a bug, never an update.
- The 20-seed stats baseline is kept across draw-order refactors, so small shifts cannot add up
  unseen. It is regenerated only for deliberate behaviour changes or an R upgrade.
- The agent generates candidates into scratch and compares them with
  `tests/regression/compare_references.R`. The owner logs the reason here (T4, T5, ...), installs
  the whole set and commits it with the change, prefixed `REFERENCE UPDATE vN:`. The previous set is
  tagged `ref-vN`. The current set is v1 (commit `2cf1bb8`).
- Phase 4 species get subdirectories of `tests/reference/`.

Comparison tested on scratch sets, 2026-11-19:
- A calibration change (rust resprout `max_effect` 0.53 -> 0.54) was classed row 1.
- An extra unused draw per year was classed row 3: statistics within 0.2 tolerance.
- The same draw change moved the single-seed A4 run from EN to CR. Single-seed outputs after a
  random-stream change are a new draw, not evidence of a behaviour change.

**T4. Configuration layers for species** [decided, 2026-11-19]
From the constants inventory (`docs/constants_inventory.md` §7). Values move out of the code into
four layers:
- species biology (`species/<sp>.yaml`);
- landscape (one shared file);
- assessment settings (`iucn_assess_spec`);
- validation data (per species, kept separate from the biology).

What goes where:
- Landscape holds what is the same for every species in a cell: fire regime defaults, fire
  pipeline settings, grid size, the rust SDM.
- The rust anchor stays with the species, because it is calibrated to that species' field data.
- Rust arrival is 2010 for all species [assumption].
- Spin-up and generation-time windows are species settings. They must be long relative to T_g, which
  is emergent, so a run warns when they are short relative to its own T_g. M. nodosa
  (T_g ≈ 23 yr) uses a gen-time window of ≈ 4.4 T_g and spin-up before it of ≈ 9 T_g. The proposed
  warnings are below 3 T_g and 5 T_g.
- `main.R` and `02b_fire_prior_spatial.R` moved to `plant_ibm/archive/`.

**T5. `parlib/`** [decided, 2026-11-19]
- `parlib/unittests/params.R` stays as the unit-test fixture. It now has the five fields it lacked
  (`resprout_delay_dist`, `fire_extent_mode`, `fire_size_sd`, `fire_stand_cover`, `rust_site_sd`),
  each set to the value the code fell back to when the field was absent: "poisson", "plant", 0, 1,
  0. Test behaviour is therefore unchanged. `test-params.R` now checks that its fields match
  `params.R`.
- The five identical `parlib/diags/*/params.R` (an older, unused calibration) are archived as
  `plant_ibm/archive/parlib_diags_params.R`.
- Phase 4 priors use the `parlib/ensemble/no_rust_no_fire/param_ranges.R` format: `__` keys,
  min, max, an optional correlation matrix `Sigma`, and a written rationale for each range. This
  is the format `runner/ensemble.R` reads. In the species YAML each prior carries the same fields
  (`prior: true`, min, max, reason; CLAUDE.md rule 5) and is converted to keys/mins/maxs.
  `DEFAULT_PARAM_DRAWS` (`runner/range_ensemble.R`, inventory S7) uses a different layout; it
  should move to this format when the YAML loader is written.

**T6. Sources for M. nodosa parameter values** [decided, 2026-11-19]
M. nodosa values sit in their own file, `data/m_nodosa_params.csv`. It has the same columns as
`data/rust_species_params.csv`, and CLAUDE.md rule 5 now names both files.
- It has one row per `params.R` value (69, with `RUST_ANCHOR`), using `a__b__c` parameter names,
  plus 3 calibration targets (V1, V5, V6).
- Run settings are left out: `N0`, `n_years`, `seed`, `fire_years`, `rust_start_year`.
- Sources cite the entries in this file.
- New `evidence_type` labels: calibration (fitted, e.g. V4), model_assumption, model_choice (a rule
  or form).

19 active values could not be traced to any note, README or commit. They are marked as priors
[assumption]: `source` = "prior: ... T6", and notes start `PRIOR [min, max]` (numbers) or
`PRIOR {a, b}` (alternative forms or rules), then a reason. Where the range is not taken from
`parlib/ensemble/no_rust_no_fire/param_ranges.R` it was proposed on 2026-11-19 and is open to review.
The point values are unchanged: the priors document uncertainty and nothing was refitted. The 19:
- senescence (all four);
- juvenile decline (`max_effect`, form, `half_sat`, age `half_sat` and hill);
- `age_first_flower_sd`;
- fire kill age `half_sat` and hill;
- the resprout recovery rule;
- rust hazard by age (four);
- the resprout rust curve's form and `half_sat`.

`rust_dose_response__resprout__half_sat` is confounded with the V4-calibrated `max_effect` while
`rust_pressure` = 1.

18 rows stay marked `UNTRACED`:
- 16 are inactive (switched-off terms, the legacy fire mode, shape settings the chosen form ignores);
- 2 are the single-site fire defaults, which belong to the landscape layer (T4).

**T7. YAML schema and loader** [decided, 2026-11-19]
`plant_ibm/R/species_config.R` assembles the params list from four files and validates it:
- `config/schema.yaml`: every field's layer, type, hard bounds or choices, and the list order;
- `species/<sp>.yaml`: biology and rust anchor;
- `config/landscape.yaml`: shared fire defaults (T4);
- `config/run_defaults.yaml`: `N0`, `n_years`, `seed`, `fire_years`, `rust_start_year`.

For M. nodosa, `load_species("m_nodosa")` is `identical()` to `get_default_params()`, so
`make_run_id()` hashes and result caches are unchanged. `load_rust_anchor()` is identical to
`RUST_ANCHOR`. `params.R` is not yet switched over (next PLAN item).

Defined forms:
- Every YAML number is parsed by R's `as.numeric()`, then coerced to the schema type. A YAML
  `10000` therefore becomes a double, not an integer.
- Fractions may be written as `1/3` or `1/50`, and `.inf` reads as Inf.
- `annual_env_alpha` and `annual_env_beta` are written `{derived: rust_anchor}` and computed from
  the anchor exactly as `params.R` did. This matters: `(1 - 0.7) * 20` is not exactly 6.

Rule 5 at load time:
- Each species value carries `cite: <csv>`. The loader checks that the row (species + `__` path)
  exists once and that the value matches.
- Otherwise the value is a prior (`prior: true`, min/max or choices, reason). The loader checks
  the value lies within its range and that the range matches the CSV's PRIOR note.
- Landscape values need a `source` or must be a prior. Run settings are not cited.
- Fields the schema does not know, or in the wrong file, are errors. All problems are reported
  together.

`species/m_nodosa.yaml` was generated from `data/m_nodosa_params.csv`. A test checks that every
CSV parameter row is used.

The 18 rows still `UNTRACED` after T6 are now priors:
- The 16 inactive values are fixed at their value (range [v, v]), with the reason "inactive".
  This avoids implying uncertainty about settings that have no effect.
- Shade `max_effect` is the exception: [0, 0.95], with the upper bound from `param_ranges.R`.
- The single-site fire defaults are bounded by the fire pipeline's floor and cap:
  `fire_prob_annual` [0.005, 0.5], `fire_p_fimp` [0.05, 1]. Their canonical values are in
  `config/landscape.yaml`.

**T8. `params.R` reads the species configuration** [decided, 2026-11-19]
`plant_ibm/params.R` no longer holds values. When sourced, it loads `species_config("m_nodosa")`
once and keeps its old interface:
- `get_default_params()` returns the loaded list;
- `RUST_ANCHOR` is the loaded anchor.

Every runner, diagnostic and runscript that sources it is unchanged. It finds the repo from its own
location, so it can be sourced from `plant_ibm/`, `diagnostics/`, the repo root, with `chdir`, or
by `sys.source` into an environment; `MNSIM_ROOT` overrides the location. A configuration error
stops the source.

Checked before the switch: the new and old `params.R` give identical `get_default_params()` and
`RUST_ANCHOR` in all five ways of sourcing. All regression tests pass: frozen outputs, and
"params.R still builds the frozen params".

The comments that documented values in the old file are covered by the CSV notes and sections C,
D, E and V of this file. Loading takes about 0.3 s per source.

Still to do (PLAN Phase 2): the constants in `runner/` and `spatial/` (inventory sections 3–4).
Done in T9.

**T9. Runner and spatial constants moved to the configuration** [decided, 2026-11-19]
The constants listed in `docs/constants_inventory.md` sections 3–4 now come from configuration files,
following the four layers of T4. The loader (`R/species_config.R`) reads and checks them.

| What | Was | Now |
|---|---|---|
| Rust year, present, end of run, 11-yr averaging, 100-yr cap, criterion | `iucn_a4.R` (3 copies), `range_ensemble.R` | `config/assessment.yaml` (shared) |
| Spin-up 1600, gen-time window 1810–1910 | same | `assessment:` block of `species/m_nodosa.yaml` (species settings, decision 3) |
| Category thresholds 30/50/80 | `iucn_a4.R` (5 places), `range_ensemble.R` | `config/assessment.yaml`, `categories:`; read as `IUCN_CATEGORIES` |
| Pegg years, observed 0.48 of 50, site, fire overrides, plot file, `n_post` | `pegg_validation.R`, `pegg_sites.R` (twice) | `species/m_nodosa_validation.yaml` (layer V); read as `PEGG`. The observed value is read from its CSV row; the rust year is the assessment's |
| `DEFAULT_PARAM_DRAWS` (max_effect U(0.32, 0.78)) | `range_ensemble.R` | `ensemble_draws:` block of the species file; the range must match `RANGE [0.32, 0.78]` in the notes of the parameter's CSV row |
| ALA query, range box (3 copies) | `spatial/01_…` | the species file's `species:` name and `spatial: range_box:` |
| Grid size (5 copies), NSW/Qld border latitude (2 copies), rust SDM path (3 copies) | `spatial/`, `range_ensemble.R` | `config/landscape.yaml`, `spatial:` (`MN_SDM_TIFF` still overrides the SDM path) |
| `mn_*` data files and figures | `spatial/data/mn_*.csv`, `spatial/figures/mn_*.png` | `spatial/data/m_nodosa/` and `spatial/figures/m_nodosa/` (`git mv`; names without the prefix) |
| Age at first flowering draw (2 copies) | `R/individuals.R` | one helper, `draw_age_first_flower()`, same draws |

Further rules and changes:
- **Default species.** `config/run_defaults.yaml` names it (`species: m_nodosa`). The runners and
  `params.R` use it when no species is given, so the code names no species.
- **Settings need a source.** CLAUDE.md rule 5 covers parameter values (`params`, `rust_anchor`),
  which must cite a CSV row or be priors. Settings outside the params list (assessment windows,
  range box, validation set-up, grid size) are not parameter values. Each needs a `source:` string,
  as landscape values do, or must take its value from a CSV row (`cite:` + `row:`). The schema
  (`config/schema.yaml`, new sections) declares their types and bounds. The loader also checks that
  the assessment years are in order, that the thresholds increase, and that a draw range contains
  the parameter's value.
- **Repo root.** `mnsim_root()` now searches up from the working directory for
  `config/schema.yaml`, so tests and scripts find the configuration from any folder inside the repo.

Identity: every moved value is the same number and type as before (integer years, double
thresholds), so run IDs and stored results are unchanged. Tests check:
- `iucn_assess_spec()`, `DEFAULT_PARAM_DRAWS` and the Pegg settings are `identical()` to the old
  literals;
- every regression test passes (the a4 and pegg scenarios use the defaults that moved).

Not moved (inventory: keep in code), and open points:
- numerical settings and model structure (inventory section 2; V9–V10, S8–S10, P8–P11, P14);
- the fire-pipeline settings in `spatial/fire/fire_config.R` (record years, threshold, floors,
  regimes; D1–D5): already in one shared file.
- ~~Open: `QLD_BBOX`~~ (the Qld download extent, drawn by hand around the M. nodosa cells): now
  derived from the grids, T10.
- ~~**Open: cell ids and fire tables.**~~ Decided in T13 (master grid, shared tables). Was: `cell_id`
  is the cell number in a raster spanning the species' range box. The fire cell-year tables and `spatial/data/fire/cells.rds` are built for the
  M. nodosa grid. A second species needs its own (or a common, fixed) grid before the fire tables
  can be shared.
- `talk/` keeps its M. nodosa references (talk figures); only its paths and the renamed constants
  were updated.

**T10. Queensland fire download area derived from the species grids** [decided, 2026-11-19]
The hand-set `QLD_BBOX` (147.9, −29.3, 154.2, −23.9) is replaced by `qld_bbox()` in
`spatial/fire/R/bbox.R`. It is worked out as follows:
- take the cells, over **all** species grids (`spatial/data/*/grid.csv`), that lie north of the
  border latitude minus the border band (−28.16 − 1°). Cells are assigned to a state by the state
  polygons, so this is a safe superset of the Queensland cells;
- take their extent at the cell edges;
- add a margin, `qld_bbox_margin`: 0.25° in `config/landscape.yaml`.

The Qld fire data are a shared landscape layer (T4 decision 2), so one download serves every species.
For M. nodosa the cells span 148.10, −29.20, 153.60, −24.10, inside the old box. The downloads
already made still cover them, and the old box's eastern reach to 154.2° was not needed.

A guard stops a too-small download from being used without notice:
- each finished download records its area in `<file>.bbox`;
- `download_layer()` fetches a layer again if its recorded area does not cover the cells (a
  finished download used to be skipped whatever its area);
- `02_fire_cell_years.R` stops before computing the Qld tables if a download misses part of the
  cells;
- downloads made before this change have no record. They give a warning, not an error;
  `record_download_bbox(file, c(xmin = 147.9, ymin = -29.3, xmax = 154.2, ymax = -23.9))` records
  the old box.

No model output changes: the area only selects which polygons are downloaded, and the cell-year
tables in git are untouched. Still open (T9): cell ids depend on each species' range box (decided
in T13).

**T11. Warning when a window is short relative to the run's T_g** [decided, 2026-11-20]
Carries out T4 decision 3 (inventory §7). T_g comes out of each run, so it cannot set the spin-up
or the gen-time window in advance. Each A4 run is checked afterwards instead, by
`check_tg_windows()` in `runner/iucn_a4.R`:
- gen-time window (`year_gentime_end − year_gentime_start`) must be at least **3 T_g**;
- spin-up before it (`year_gentime_start − year_spinup`) must be at least **5 T_g**.

The thresholds are assessment settings shared by all species: `window_checks` in
`config/assessment.yaml`. They are set from M. nodosa's windows, which were judged adequate: with
T_g ≈ 23 yr, its gen-time window is ≈ 4.4 T_g and its spin-up before it ≈ 9 T_g. So M. nodosa
runs pass with room to spare. A species with T_g above ≈ 33 yr would fail the 3 T_g test with
M. nodosa's 100-yr window, and one above 42 yr would also fail the 5 T_g spin-up test.

How a short window shows up:
- `run_iucn_a4()` gives one warning per short window, of class `mnsim_tg_window_warning`, so a
  caller can muffle just these. The shortened regression A4 scenario does this: its spin-up is
  100 yr ≈ 4.4 T_g on purpose.
- `summarise_iucn_a4()` adds a column `tg_windows_ok`. It is computed from the result's own T_g and
  years, so results saved before this change are checked too. It is not stored in the result,
  which keeps the result identical to the frozen reference.
- `collect_iucn_results()` and `collect_range_ensemble()` print how many runs had a short window.
  `range_ensemble.R` suppresses run warnings, so this count is where such runs are reported.

The check is a warning, not an error. The assessment still uses the run, because the window
length L = 3 T_g is itself still valid; what is in doubt is whether T_g was estimated from an
equilibrium population. No model output changes.

**T12. Swappable fire response, rust pressure and recruitment** [decided, 2026-11-20]
Design note `docs/design_swappable_interfaces.md`, approved with its recommendations. Phase 3
mechanisms will plug into fixed places, chosen per species in its YAML, rather than adding branches
to the annual loop. Each place is a named list of functions (a registry) chosen by a setting whose
schema `choices` are the registry's names; a test checks they match.
- **Fire response** (`R/fire.R`): the three branches of `apply_fire()` are now `fire_legacy`,
  `fire_plant` and `fire_stand` in `FIRE_RESPONSES`. `fire_extent_mode` still chooses, and
  `fire_p_fimp == 0` still means legacy. Recording deaths and the resprout delays stay shared.
- **Rust pressure**: per cell, `RUST_PRESSURE_MODES` in `spatial/rust_pressure.R`, chosen by an
  optional `mode` in the species' rust anchor (default `sdm`, the only option). Per year, the Beta
  draw is now `draw_annual_env()`, with no selector until a second mode is designed.
- **Recruitment** (`R/recruitment.R`): optional `recruitment_curve`, `bevholt` (default) or `ricker`.
  Ricker uses `ricker_r` and `ricker_K` (the field `K` renamed; `ricker()` was never called). No
  species uses Ricker: its scale is not comparable with `R_max`/`K_half`, so a species choosing it
  needs its own calibrated values (Phase 5). The expected-recruit count is `expected_recruits()`,
  and flowering (with rust suppression) is `update_flowering()`: the places for later fecundity and
  flowering rules.
- **Optional schema fields** (`optional: true`, with a `default` or `required_when`). Omitted, or
  set to the default, they are left out of the params list, so existing species' lists and run IDs
  are unchanged; the code falls back to the default. `required_when: {field: value}` makes a field
  required when that option is chosen and refuses it otherwise. Their values follow rule 5 like any
  other.

No new biology, and no M. nodosa change: the split moves code without adding, removing or
reordering random draws. Tests: the full suite with `tests/reference/` untouched (exact-match
regression, frozen params, 20-seed statistics); each fire response against a copy of the pre-split
`apply_fire()` over several seeds, including the number of draws; registries against the schema;
the optional-field rules; Ricker's expected recruits computed by hand.

**T13. One master grid for all species; shared fire tables; per-species fire prior** [decided, 2026-11-20]
Options were proposed in conversation and approved with the recommendations (A + X + per-species
prior).
- **Master grid.** A cell's `cell_id` is its number in one fixed 0.1° raster over all of
  Australia: 110–160 E, 45–9 S, 500 × 360 cells (`config/landscape.yaml`, `spatial$grid_extent`;
  `spatial/cell_ids.R`). The extent includes Tasmania, the Torres Strait islands and Lord Howe
  Island. It was made Australia-wide, not just eastern Australia, so that any later species fits
  without renumbering ("master list"). The loader checks that each species' range box lies inside
  it. Before T13, ids were numbers in a raster spanning each species' range box. The cell centres
  are unchanged: both grids sit on 0.1° steps.
- **Renumbering (one-off).** M. nodosa's ids were mapped from the old box raster
  (148–155 E, 38–24 S) to the master grid by cell centre. Only the `cell_id` field was rewritten,
  in `spatial/data/m_nodosa/{grid,thinned,site_params}.csv`, the six cell-year tables and
  `qld_chronic_cells.csv`; every other byte is unchanged, and row order is kept (range plans pick
  rows by position).
- **Shared tables.** The cell-year tables in `spatial/data/fire/` cover the union of all species'
  cells. `cells.csv` (in git; replaces the `cells.rds` cache) lists the cells covered, with state
  and land area. `02_fire_cell_years.R` adds only cells not yet covered. Each cell's burnt fraction
  depends only on that cell, so adding cells leaves existing rows unchanged. A table's rows for a
  batch are replaced, not duplicated, if a run is interrupted, and `cells.csv` is updated last.
- **Per-species prior.** `load_fire_cells()` returns the species' own cells. The fire regime and
  its empirical-Bayes prior (D2) are fitted over those cells only, as before, so adding a species
  does not change another species' site table. The fire checks (`03`), which are also
  per species, now write to `spatial/data/<sp>/fire_checks/` and figures to
  `spatial/figures/<sp>/fire/`.

Check: `cells.csv`, built with `load_cells()`, agrees with the site table's states and with
`qld_chronic_cells.csv`'s land fractions. Rebuilding with `02` + `02c` then reproduces
`site_params.csv` byte for byte, and `03` reproduces all four check files byte for byte. M. nodosa's
inputs are unchanged. Pegg site runs store the cell id in their params (`pegg_site_cell`), so their
run IDs and cache change once; range-ensemble run IDs do not use the id. Tests:
`tests/testthat/test-cell-ids.R`.

**T14. NSW cell-year tables: fire records attached to the wrong cells** [found and corrected 2026-11-20]
Found while checking T13's incremental build. `cell_year_fractions()` clipped the cells to NSW
land with `st_intersection()`, which drops a cell with no land in the state. One M. nodosa NSW
cell has none: id 125413, position 23 of 256 in grid order. Every later cell was then matched to
the next cell's clipped geometry, and its burnt fraction was divided by its own land area. The
committed `cell_years_nsw_wild.csv` is reproduced exactly by this shift: all 3,369 rows, values
within 3e-15. **234 of the 256 NSW cells hold another cell's fire record.** Because `grid.csv` is
not sorted by location, the "next" cell is often far away. `nsw_planned` and `nsw_all` went
through the same code. The Qld tables are not affected: `cells_on_land()` would have stopped on a
length mismatch, and the Qld pieces carry their `cell_id`.

The code is fixed: the pieces are matched back by position, and a cell with no land stays empty;
`test-cell-ids.R` covers this.

**Rebuilt (approved 2026-11-20).** The three NSW tables were recomputed with the fixed code (about
30 minutes), then `02`, `02c` and `03` were rerun.
- *Tables.* Same rows (3,369 / 1,773 / 4,355): the same fires, now on the right cells. No rows
  fall on the zero-land cell.
- *Site table.* Every NSW cell's fire values changed, with each cell's annual fire probability
  moving by a factor of 0.2–3.9 (5–95%; median ×1.12). Qld cells' records are unchanged, but their
  estimates moved slightly because the prior is pooled over both states (D2). Non-fire columns are
  identical.
- *Across the range.* The abundance-weighted mean fire probability went from 0.033 to 0.036, and
  the mean burnt fraction `fire_p_fimp` from 0.287 to 0.243.
- *Rust anchor cell (Yarringully).* Fire probability 0.067 → 0.082, `fire_p_fimp` 0.28 → 0.44.
- *Prior (wildfire, pooled).* Gamma(0.69, 16.7) → Gamma(0.82, 20.6), mean 0.041 → 0.040 per year.
- *D2 and D4 still hold.* The shared-prior test gives p = 0.85 → 0.89. The border check still
  rejects the automated Qld record: Qld/NSW rate ratio 4.96 → 4.66 (95% CI 3.25–6.88), against
  0.96 → 0.90 for the edited 2014–2024 record.

What must be rerun, because their fire inputs changed (each takes minutes to under an hour on an
8-core machine):
- range-wide ensembles and talk figures;
- the Pegg-site fit (V6);
- the Pegg anchor-cell runs behind V4. The Pegg scenario sets the 2020 fire itself, but the
  spin-up fire regime comes from the anchor cell.

No parameter was refitted here (CLAUDE.md rule 4). If V4's `max_effect` range needs revisiting, that
goes through the calibration code. The single-site regression references use `params.R` defaults,
not the site table, so they are unaffected. D1's earlier cross-check against the previous `02`
reported 255/256 cells matching; why it did not catch this has not been traced.

**Reruns with the rebuilt tables (2026-11-20, commit 105e31f; results in `plant_ibm/results/`, not
tracked).** All three ran on one 8-core machine in about 70 minutes in total.
- *Range-wide ensemble* (`talk/make_talk_figures.R`, `results/talk_v3/`; 300 runs, 128 cells):
  median decline 47.9% (IQR 45.5–54.9%; 90% 37.5–78.5%). Categories: LC/NT 1%, VU 66%, EN 29%,
  CR 4%; NSW VU 67% / EN 28% / CR 5%, Qld (25 runs) VU 56% / EN 40% / CR 0%. Median T_g 23.7 yr,
  window 71 yr. Close to the earlier stand-mode result (60 runs: median 48.3%, VU 60% / EN 35% /
  CR 5%). The median still sits near the 50% EN threshold, so shares stay sensitive to small shifts.
- *Pegg scenario at the anchor cell* (100 runs, max_effect drawn from the B3 range): median 51.9%
  of resprouters dead by Jan 2022 (90% 36.3–65.7%), against 48% observed (24/50).
- *V4 sweep* (`max_effect` 0.2–0.9, 40 seeds each, Yarringully cell): the median reaches 48% at
  max_effect ≈ 0.49, against 0.53 before; medians inside the observation's 95% interval
  (0.35–0.62) for 0.35–0.70, against 0.32–0.78 before. The level moved little; the range narrowed
  slightly. Nothing is changed here (rule 4): `max_effect` = 0.53 and the B3 range
  U(0.32, 0.78) stay until the calibration code refits them.
- *V6 full grid* (`results/pegg_sites_t14/`; 12 grid points × 5 cells × 60 seeds): best
  log-likelihood −31.0 at site_sd 4, max_effect 1.0, then sd 4 / 1.2 (−31.6) and sd 3 / 1.0
  (−31.9). Year-to-year variation alone fits worse. The optimum is again at the edge of the
  site_sd range, so the V6 reading stands (a large persistent site effect), and the grid needs to
  go wider. At the best set the model still predicts every plot in a broad band (each plot's
  90% range 0–62% or 4–92% dead) with a 1–6% chance of a zero-death plot.

**T15. Species fire exposure and fire-change scenarios** [decided, 2026-11-20]
Design note `docs/design_fire_exposure_and_change.md`, approved with its recommendations. It is
aimed at the rainforest species, whose habitat burns far less than the ~10 km cell it sits in, and
at the question of whether more fire, together with rust, changes their category.
- **`fire_exposure`** (species, optional, 0–1, default 1) multiplies the cell's fire probability
  for the whole run, spin-up included. It is habitat biology, so it goes in the species YAML under
  rule 5: it either cites a CSV row or is a prior with a reason (e.g. how rarely rainforest patches
  burn in the NSW record). It can be drawn over its range with `ensemble_draws`. Fire off =
  exposure 0. The runner's `fire_scale` stays separate, as a management or sensitivity scenario;
  the two multiply.
- **Fire-change scenario** (`fire_change_year`, `_rate`, `_size`, `_ramp`; set by runners, e.g.
  `run_range_ensemble(..., fire_change = list(year = 2026, rate = 2, size = 1.5))`). From the change
  year, the fire probability is multiplied by `rate` and the mean burnt fraction `fire_p_fimp` by
  `size` (kept below 1), so stand fires reach the stand more often and kill more (D3). A `ramp`
  > 0 makes the rise linear over that many years; 0 is a step.
  - The settings are written into params, so they are hashed into run IDs.
  - Cells, seeds and parameter draws do not depend on the scenario, so scenario ensembles pair run
    by run with the baseline.
  - A run is identical to the baseline, draw for draw, before the change year.
  - `run_iucn_a4()` refuses a change year at or before the end of the gen-time window, so T_g
    describes the pre-rust population.
  - The suggested year is 2026. The cell records (to 2024) already include 2019–20, so the scenario
    means "more fire than the recent record".
- **Code.** `fire_prob_at()` and `fire_params_at()` in `R/fire.R`, called each year by
  `run_simulation()`. When the settings are absent they return exactly the old values, so
  M. nodosa is unchanged.

These are scenarios, not predictions: the size of a future fire increase in these vegetation types
is not estimated from data. They are reported as the change in the category distribution against
the paired baseline. The fire *response* of the rainforest species (kill, resprouting) is species
biology and needs Phase 4 priors. Tests: `tests/testthat/test-fire-change.R`.

**T16. *Uromyrtus australis* added** [decided, 2026-11-20]
Added at the user's request: 37 rows in `data/rust_species_params.csv`. Sources: NSW TSSC
conservation assessment (Saunders 2025) and final determination (8 May 2026, Critically
Endangered, B1ab(iii,v)); EPBC Conservation Advice (2021); Pegg et al. (2025) Austral Ecol;
Pegg et al. (2017) CRC2063 Table 13; Makinson (2018). Much of the recent monitoring is
R. Kooyman *in litt.* (2024–25), as cited by the committee. The primary reports (Kooyman 2005 MSc
thesis; Kooyman 2021 and 2024 SoS reports; DPE/DCCEEW SoS report cards 2021–23) were not found
online.
Why it fits this project:
- *Fire × rust.* It is the clearest case of the interaction. Its refugial rainforest
  rarely burns, but 2019–20 burnt ~30% of the range and 29% of monitored genets. Stem kill by
  size class is measured in three quadrats. Rust then infected 55–85% of reshoots, and no
  significant fruiting has occurred since 2018. This is the case `fire_exposure` and the
  fire-change scenario (T15) were built for.
- *Criterion A is open.* The committee found Criterion A Data Deficient because generation length
  is unknown, so a model-based A4 assessment would fill a stated gap.
What it needs before it can run (not built; rule 7):
- *Clonal stage.* Reproduction is by root suckers only: no wild seedling seen in >40 yr, and seed
  viability is ~20 weeks. The mature individual is the ramet (IUCN), but genets carry
  persistence. The current IBM recruits from seed, so the Phase 3 clonal-stage design note must
  cover genet → ramet structure.
- *Generation length.* Stems live >1,500 yr and age at first suckering is unknown, so T_g from the
  run's age structure needs a clonal definition. Any 3-generation window will hit the 100-yr cap.
- *Rust timing.* Rust was first seen on the species in 2018–19, eight years after it reached NSW.
  The rust start year is a species setting, not 2010.
- *Spatial.* EOO is 89 km²; the range (~13 × 11 km) covers a few master-grid cells, so range-wide
  draws reduce to one or two cells' fire and rust.

**T17. Clonal stage: genets, ramets and root suckers** [decided, 2026-11-21]
Design note `docs/design_clonal_stage.md`, approved with its recommendations. Needed by
*U. australis* (suckering only), *R. psidioides* (populations reduced to root suckers) and
*R. rubescens* (clonal clumps). Code in `R/clonal.R`; it is optional (`clonal_mode: genet`), and
absent it makes no draws, so M. nodosa is unchanged (regression passes).
- **Structure.** Each row is a ramet, the IUCN mature individual, with `genet_id` and `origin`.
  Ramets share their genet's genotype and microsite (decision 3: suckers stay beside the parent).
  Founders are grouped into `clonal_founder_genets` genets; each seedling starts a genet.
- **Suckering.** Each living genet produces Poisson suckers each year from its mature ramets
  (`sucker_rate`), from ramets killed the year before (`sucker_pulse`: response to fire, rust
  dieback or slashing), and, with no mature ramet, from the bare rootstock
  (`sucker_rate_rootstock`). Suckers face the existing hazards, so rust's young-plant peak kills
  them early.
- **Rootstock.** It survives a fire that leaves the genet with no ramet with
  P = `rootstock_fire_survival`, and dies with P = `rootstock_hazard_no_canopy` in each year
  without a mature ramet, taking its remaining ramets with it. That hazard sets how long a
  "living dead" sucker pool lasts.
- **Generation length** (decision 4a). T_g is the mean age of the parents of all new
  individuals, seed and sucker; a sucker's parent is a ramet of its genet (a mature one if any).
  Without suckers this is the old calculation exactly.
- **Census.** Adds `N_genets`, `N_rootstock_only` and `births_sucker`, reported alongside
  `N_IUCN` (ramets), not assessed.
Two refinements, found in smoke runs before any species used the code:
- *Suckers saturate with the recruitment curve.* The design limited suckers by canopy
  shading only. With M. nodosa's shading (never full suppression) a mixed seed + sucker run grew
  from 1,500 to ~260,000 stems in 250 years. Suckers are now also multiplied by the recruitment
  curve's per-capita saturation at the current canopy (`sucker_density_factor()`: Beverton-Holt
  K_half / (K_half + N); Ricker exp(−r N / K)). This reuses the species' own K, adds no
  parameter, and makes seed and suckers compete for the same space. The same run then settles
  near 21,000 stems.
- *Only established genets keep a bare rootstock.* A genet keeps a rootstock with no living ramet
  only once it has had a mature ramet (founders count). Otherwise every seedling that died young
  left a persisting rootstock (~12,000 in the same run).
[assumption] Pollination between ramets of one genet counts as outcrossing in `sample_parents()`
(within-genet selfing is not distinguished). [open] Stem age stands in for stem size (fire kill,
maturity); U. australis quadrat data are by height class. Smoke runs (illustrative values, not
priors): a clonal-only run (150 genets) holds ~1,000 mature stems for 250 years and declines
under rust; T_g ≈ 23 yr from sucker parents. Tests: `tests/testthat/test-clonal.R`.

**T18. Rust-pressure modes; the anchor as a property of place; rust timing per species** [decided, 2026-11-21]
Design note `docs/design_rust_pressure_modes.md`, approved with its recommendations (all five
decisions). Needed by *R. psidioides* (no climate predictor of infection), *A. beckleri* (a
possible north–south gradient, C4) and *U. australis* (rust first seen on it in 2018–19). All new
settings are optional; M. nodosa is unchanged (its site table is byte-identical after rerunning
the rust step; params, run IDs and the regression tests are unchanged).
- **The anchor is a property of place** (decision 1). The `rust_anchor` block (site, level 0.70,
  concentration 20, cap 0.95) moved from `species/m_nodosa.yaml` to `config/landscape.yaml`, with
  the same values; they still come from the `rust_anchor__*` rows of `data/m_nodosa_params.csv`
  (a test checks that they agree). A species file may set only the mode. Every species therefore
  sees the same yearly pressure in a cell; species differ through their host terms
  (`rust_dose_response`, `rust_pressure`, resistance), fitted in Phase 5.
  [assumption] The 0.70 was inferred from M. nodosa incidence, so it mixes host and place.
  Calling it place is a convention: other species' host terms are fitted on M. nodosa's scale.
- **A fixed anchor cell** (found while building, not in the note). The anchor cell had been "the
  occupied cell nearest the anchor site", which depends on the species. It is now a fixed
  master-grid cell, `cell_lon`/`cell_lat` = 153.25, −29.15 (cell 100933, the cell M. nodosa has
  always used; new CSV rows `rust_anchor__cell_lon`, `rust_anchor__cell_lat`). A species that
  does not occupy it normalises to the SDM read from the raster for that cell, with the same
  extraction, so a shared cell gets the same pressure whichever species it is computed for
  (tested). Without the raster, such a species' table cannot be re-normalised (error).
- **Modes** (decision 2; `RUST_PRESSURE_MODES`, chosen by `rust_anchor$mode`):
  `sdm` (default, unchanged); `constant` (μ in every cell, SDM or not; `sdm_rel` kept for
  reference); `sdm_gradient` (the `sdm` table, plus a per-run latitude shift).
- **The gradient** (decision 3). `apply_rust_gradient()` (`runner/site_params.R`) shifts the
  cell's logit mean pressure by g × (lat_anchor − lat) / 10, keeping the Beta concentration:
  g is the change in logit pressure per 10° **south** of the anchor cell (g < 0: less rust to the
  south; g = −1 roughly divides the odds by e ≈ 2.7 at 10° south). This corrects the sign in the
  design note's formula, which with its example range U(−1, 0) would have raised pressure to the
  south. The anchor cell's latitude keeps its mean exactly, so V4 and the Pegg calibration are
  untouched. g = `rust_gradient`, an optional run setting (default 0, left out of the list, so
  g = 0 gives the `sdm` run and its run ID exactly). It is drawn per run by listing it in a
  species' `ensemble_draws` (the schema now allows drawing run settings marked `drawable`), or
  fixed for every run with `run_range_ensemble(rust_gradient = g)`, which still makes the draw so
  plans pair up. A non-zero gradient with any other mode is refused. The plan table gets a
  `rust_gradient` column for `sdm_gradient` species only.
- **Timing** (decisions 4, 5). Optional species settings `rust_arrival_year` (default −Inf) and
  `rust_ramp_years` (default 0 = step). At the start of a run, `apply_rust_timing()` sets
  `rust_start_year` to the later of the run's rust year and the species' arrival, so every rust
  term (hazards, flowering, recruitment, resprout delay, the census's `rust_active`) starts
  together, and a counterfactual with no rust stays rust-free. The ramp multiplies the year's
  draw by min(1, (t − start + 1) / ramp), the same convention as `fire_change_ramp`, with no extra
  random numbers. It also scales the mean pressure used for the resprout delay. The A4 window is
  not moved.
Not done: year-to-year structure, a new SDM, host susceptibility per species (Phase 4/5).
Tests: `tests/testthat/test-rust-modes.R`.

**T19. Accumulating rust damage** [decided, 2026-11-21]
Design note `docs/design_accumulating_damage.md`, approved with its recommendations. Why: the
rust hazard had no memory, so under steady pressure survivors' death rate stays flat or falls
(selection removes the most susceptible first). The data show it rising: *R. rubescens* paired
resurveys 2014 → 2016 (Pegg et al. 2017, Table 10) rose at 4 of 5 sites (Bongil Bongil ~3%/yr to
2014, ~25%/yr after); *R. psidioides* ~21%/yr to 2014, ≥ 32%/yr among survivors to 2016; at Olney
trees lost their crowns, stopped flushing and died. No parameter values of the old hazard give
that. Optional (`rust_hazard_mode: damage`); in the default `instant` mode nothing new runs, so
M. nodosa is unchanged (regression passes).

Decisions required by the note, and what was decided:
1. *How damage is represented.* Options: a continuous damage state, or a counter of consecutive
   high-pressure years (PLAN's "consecutive defoliation events"). **Decided: a continuous state**
   D in [0, 1], the share of the crown lost. The counter resets abruptly after one light year
   and needs a "defoliation event" threshold that no data give; D is in effect a smoothed counter.
2. *Whether the existing (instant) rust hazard stays.* Options: keep it alongside damage, or
   damage only. **Decided: keep it.** It carries the quick death of seedlings and young suckers
   that the *R. rubescens* and *R. psidioides* data report; adults die mainly by damage. A species
   that dies mainly by damage sets a low `age_floor`.
3. *What fire does to damage.* Options: reset it for a top-killed plant that resprouts, or carry
   it through as depleted reserves. **Decided: reset** (`damage_fire_reset: reset`, the default;
   `keep` is available). The resprouting canopy is new; the existing resprout term carries its
   extra risk ("resprouting gives no escape", *R. psidioides*).
4. *Whether damage can recover.* Options: a recovery rate, or damage that only accumulates.
   **Decided: include recovery** (`damage_recovery`, whose prior range must include 0). Evidence:
   *A. beckleri* at Knorrit had seedlings without dieback after the last infection event, and
   pressure varies between years.
5. *Effects of damage on flowering and fruit.* Options: here, or in the fruit-infection note.
   **Decided: the fruit-infection note.** D is kept as a column so that note can read it.

Implementation (`R/mortality.R`, `R/fire.R`, `R/census.R`, `R/simulate.R`):
- Each year, after fire and flowering and before mortality (`update_rust_damage()`):
  D ← D (1 − `damage_recovery`); then D ← D + (1 − D) × gain, with
  gain = `dose_response(rust_pressure × eff_susc, damage_gain_dose_response)`, capped at 1, and
  eff_susc = (1 − resist) × microsite × that year's pressure. Deterministic: no random numbers.
- Death risk `dose_response(D, damage_mortality_dose_response)`, combined as one more independent
  risk with the existing union. The mortality draw is unchanged (one Bernoulli per plant).
- Founders, seedlings and suckers start at D = 0; a run continued from an instant-mode state
  starts every plant at 0.
- Census adds `mean_damage` and `share_damaged` (mature plants with D ≥ 0.6,
  `DAMAGE_REPORT_THRESHOLD`, a reporting setting) in damage mode only, for crown transparency and
  dieback targets (transparency includes a healthy baseline of 0.30–0.35).
[assumption] Damage gain does not depend on age or on D itself (damaged plants flushing more
and so more exposed is not modelled; no data separate it). [open] The ~6-month lag at Olney is
below the annual step; only the multi-year lag is represented. [open] Test against the Olney exclusion
trial (Carnegie et al. 2016): canopy recovery of fungicide-treated trees is the best available
evidence on `damage_recovery` (canopy tissue turnover with rust removed), and the untreated and
partially treated trees test gain, lag and the damage death risk (PLAN Phase 5). All settings are Phase 4 priors,
fitted in Phase 5 (targets: Table 10 conditional mortality, Olney 53% in < 4 yr, *R. psidioides*
57% at 3.5 yr, crown transparency). Check (test, illustrative values): under constant pressure
survivors' annual death rate stays ~11% in instant mode and climbs from 17% to 50–60% in damage
mode. Tests: `tests/testthat/test-damage.R`.

**T20. Olney exclusion trial as a validation of accumulating damage** [decided, 2026-11-21]
Design note `docs/design_olney_validation.md`, approved with its recommendations. The Olney SF
disease exclusion trial (Carnegie et al. 2016; *R. rubescens* only) separates the damage-mode
processes (T19) that field mortality alone cannot: fungicide switches rust off, so treated trees
show recovery alone, untreated trees the full chain, partially treated trees damage rebuilding
from a near-healthy crown.
Data (extracted 2026-11-21, `data/validation/`): crown transparency per arm and date, Aug 2011 –
Oct 2014, digitised from the plotted means in Fig. 3 (± ~1 point); percentage dead per site for
*R. rubescens* (43 sites, Table 2) and *R. psidioides* (18 sites, Table 4); summary and pointer
rows in `data/rust_species_params.csv`, with healthy crown transparency (0.30–0.35 and 0.25–0.35).
The existing Olney 53% row was corrected: it is a separate count of the stand (Table 2, "Olney
SF 1"), not the trial trees. Reading of the treated arm: damage (baseline 0.30) fell from ~0.51 to
~0.20 in 4 months and ~0.02 by mid 2013, an annual recovery share of ~0.8–0.95.

Decisions required by the note, and what was decided:
1. *How treatment arms are run.* Options: one run per arm with run-level protection, or one run
   with a per-plant protection flag. **Decided: one run per arm**, each continued from the same
   stand at the trial's start; the trial trees do not interact in any way that matters over three
   years, and no new column is needed.
2. *Yearly pressure across arms.* Options: the same series for all arms, or independent draws.
   **Decided: the same series** (`annual_env_series`), plus the same seed per year, so arms differ
   only by treatment (they shared the weather).
3. *Time resolution.* Options: annual comparison, or a sub-annual step. **Decided: annual**, at the
   trial's August values (model years start in August; the end of model year t is compared with
   August t + 1). The ~6-month lag and seasonal dips are not targets.
4. *Damage to crown transparency.* Options: T = b + (1 − b) D with b a prior, or b fixed at 0.30.
   **Decided: b a prior in [0.30, 0.35]**; dead trees count as 1; the arm's value is the mean over
   its trees.
5. *R. psidioides recovery.* Options: its prior taken from the *R. rubescens* posterior, or a wide
   independent prior. **Decided: from the R. rubescens posterior**, logged here as an
   [assumption]: both are evergreen rainforest-margin Myrtaceae and similarly susceptible, and
   *R. psidioides*' own data (nearly all trees near full damage) cannot identify recovery. Its gain
   and damage death risk are fitted to its own survey and site data.
6. *Order of work.* Options: build the hooks and runner now with toy tests, or wait for Phase 5.
   **Decided: build now**; run the validation once the *R. rubescens* species file exists (Phase 4).

Implementation:
- `rust_protection` (runner setting): protected periods (start, end pairs, in model years). The
  year's pressure is multiplied by the unprotected share of the year, in `draw_annual_env()`
  after the draw (no extra random numbers), and the same factor scales the mean pressure in the
  resprout delay. The note proposed applying it to damage gain and the instant hazard only;
  scaling the year's pressure is simpler and also covers rust effects on flowering and
  recruitment, as fungicide would.
- `annual_env_series` (runner setting): a pressure per year, named by year, used instead of the
  Beta draw (no rbeta). Both settings are set by runners only, not by config files, so they are
  not in the schema (like `rust_site_z`); absent, runs are unchanged.
- `runner/olney_validation.R`: stand simulated to the end of model year 2010 (one rust year), 20
  living trees in a given age range drawn as the trial trees (10 untreated, 5 treated to Oct 2014,
  5 treated to Jun 2013), each arm then continued year by year to 2013; outputs mean transparency
  per arm and year, and stand mortality (untreated run) by the end of model year 2013 (target
  53%). The trial trees' age range is an argument, because the age-height link is a species
  prior; it goes into `species/r_rubescens_validation.yaml` (block `olney`) with the species file.
[open] Ask the authors (NSW DPI) for tree-level data and the incidence, severity and flush
series (Supplementary Fig. 1), which would allow a tree-level likelihood.
Tests: `tests/testthat/test-olney.R`.

**T21. Rust suppression of flowering and fruit** [decided, 2026-11-21]
Design note `docs/design_fruit_infection.md`, approved with its recommendations, except that the
alternative mature count is not built (decision 5). Evidence: *R. rubescens* fruit on untreated
branches all infected, none matured (Tucki Tucki NR); *R. psidioides* no mature fruit, flowering
or seed since 2010; *U. australis* no significant fruiting since 2018 (seed recruitment ~0 even
before rust). Checking the existing hooks (T12) found three gaps: they ignored the cell's rust
pressure; fruit failure shrank the number of recruits but not who the mothers were, so there was
no selection for resistance through fecundity; and flowering suppression was redrawn every year,
so it never persisted.

Decisions required by the note, and what was decided:
1. *Site pressure in the hooks.* Options: include the cell's `rust_pressure`, or leave it out.
   **Decided: include it**, as in the hazard and damage terms (T18): both hooks now use
   `rust_pressure × eff_susc`. No effect on *M. nodosa*, whose hooks are off.
2. *Who fruit failure acts on.* Options: weight mothers by fruit success; weight mothers and
   fathers; weight neither. **Decided: mothers only.** Pollen comes from flowers, which exist; the
   evidence is about fruit.
3. *What drives flowering suppression.* Options: an optional `damage` driver (crown damage D) with
   fruit maturation on the year's pressure; both stages on D; or a separate reproductive state
   with its own recovery rate. **Decided: optional driver**, `rust_flower_driver: damage`
   (default `pressure`); fruit maturation stays on the year's pressure, because fruit is infected
   directly in the season it develops. A separate state would add a parameter no data identify.
   [assumption] Buds are set on last season's crown: flowering in year t reads D at the end of
   year t − 1 (flowering comes before the damage update in the year's order).
4. *Phenology.* Options: in each species' priors, or a sub-annual step. **Decided: in the priors**
   (whether fruiting falls in the rust season goes into the fruit dose-response priors, with the
   reason recorded).
5. *Reporting and the mature count.* Options: census columns plus an alternative count
   `N_fruiting_capable` (mature plants with fruit success above a threshold), or census columns
   only. **Decided: census columns only; the IUCN mature count is unchanged (A4).** The user
   considered this question earlier and judged it cleaner and more defensible to decide when a
   living plant will never reproduce, rather than counting by this year's fruit success. That is
   a later item: quantify, per plant, the probability of reaching a reproducing state, and treat a
   plant as not mature only when it falls below a threshold. Example: *R. psidioides* ramets in a
   phase of repeated suckering after canopy death, which may stay alive for years with little
   chance of ever fruiting again. [open] (PLAN.md, Phase 6.)
6. *Order of work.* Options: build now with toy tests, or wait for the species files. **Decided:
   build now**; per-species priors in Phase 4. Fruit ≈ 0 is a key Phase 5 target for
   *R. psidioides* and *R. rubescens*.

Implementation:
- `fruit_success()` (`R/recruitment.R`): each flowering plant's weight,
  1 − `dose_response(rust_pressure × eff_susc, rust_recruit_dose_response)`; NULL when off.
  `expected_recruits()` sums the weights; `recruit()` passes them to `sample_parents()` as mother
  weights. When the hook is off, no weights are passed and draws are unchanged.
- `update_flowering()` (`R/simulate.R`): driver `pressure` or `damage` (`rust_flower_driver()`,
  `RUST_FLOWER_DRIVERS`). The `damage` driver works in years of no pressure as well, so
  suppression persists while the crown stays damaged; it needs damage mode (a loader error
  otherwise). Probabilities are capped at 1, now that site pressure can exceed 1.
- Census (`repro_census()`, only when either hook is on): `share_flowering_suppressed` (share of
  eligible plants kept from flowering by rust, counted at flowering, before mortality) and
  `fruit_success` (mean weight over flowering plants).
- Change in behaviour when the hooks are on: earlier runs with the hooks on (diagnostics only; no
  species uses them) would differ, through site pressure and mother weights.
Tests: `tests/testthat/test-fruit.R`.

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
| Nov 2026 | D3 | Fire-only resprout window fixed at 2 years (was Poisson, mean 3), set before the Phase 2 fit |
| 2026-11-20 | T16 | *Uromyrtus australis* added to the species set; needs the clonal stage before it runs |
| 2026-11-21 | T17 | Clonal stage: ramets in genets, root suckering, rootstock survival; optional, off for M. nodosa |
| 2026-11-21 | T20 | Olney exclusion trial as damage-mode validation: one run per arm, shared yearly pressures, annual comparison; R. psidioides recovery prior from R. rubescens |
| 2026-11-21 | T19 | Accumulating rust damage: per-plant crown damage with recovery, death risk from damage, reset by fire; optional, off for M. nodosa |
| 2026-11-21 | T18 | Rust anchor is a property of place (landscape file, fixed anchor cell); modes constant and sdm_gradient; species arrival year and ramp |
| Nov 2026 | V6 | Phase 2: fit level and between-plot spread of post-fire mortality at seven Pegg plots, by interval |
| 2026-11-19 | T1 | Six fixed-seed regression scenarios from `params.R` at test scale (density parameters x 0.03); references store their params |
| 2026-11-19 | T2 | Second check: mature mortality, Pegg fraction dead and rust mortality over 20 seeds, tolerance 4 sd of the difference |
| 2026-11-19 | T3 | Reference-update procedure (docs/reference_updates.md); stats baseline kept across draw-order refactors; owner installs references |
| 2026-11-19 | T4 | Four config layers (species, landscape, assessment, validation); rust arrival 2010 for all species; windows are species settings with T_g-relative warnings |
| 2026-11-19 | T5 | parlib: unittests fixture completed (schema test added), diags archived, param_ranges format for Phase 4 priors |
| 2026-11-19 | T6 | M. nodosa values in data/m_nodosa_params.csv (rule 5 names both CSVs); 19 untraced active values marked as priors |
| 2026-11-19 | T7 | YAML schema + loader (config/, species/m_nodosa.yaml, R/species_config.R); identical to get_default_params(); rule 5 checked at load; remaining 18 untraced rows made priors |
| 2026-11-19 | T8 | params.R switched to the loader (same interface, identical output) |
| 2026-11-19 | T9 | Runner and spatial constants moved to config (assessment.yaml, species assessment/spatial/ensemble_draws blocks, m_nodosa_validation.yaml, landscape spatial); default species in run_defaults.yaml; data in spatial/data/m_nodosa/ |
| 2026-11-19 | T10 | Qld fire download area derived from all species grids (cells north of border − band, + 0.25° margin); downloads record their area; coverage checked |
| 2026-11-20 | T11 | A4 runs warn when the gen-time window < 3 T_g or the spin-up before it < 5 T_g (config/assessment.yaml); summaries flag them (`tg_windows_ok`) |
| 2026-11-20 | T12 | Swappable parts: registries for fire response, cell rust pressure and recruitment curve (Ricker wired in, unused); `draw_annual_env()`, `expected_recruits()`, `update_flowering()`; optional schema fields keep existing species unchanged |
| 2026-11-20 | T13 | One Australia-wide master grid (110–160 E, 45–9 S, 0.1°) for cell ids; fire tables shared over all species' cells and built incrementally; fire prior fitted per species |
| 2026-11-20 | T14 | Error found and corrected: NSW cell-year tables were shifted by one cell after a zero-land cell (234/256 NSW cells held another cell's record); tables, site table and fire checks rebuilt; D2/D4 still hold; range, Pegg-site and V4 runs to redo |
| 2026-11-20 | T15 | Species `fire_exposure` (optional, habitat share of cell fire rate); fire-change scenario (from a year: rate and fire size multiplied, optional ramp), set by runners, paired with baseline |
| 2026-11-21 | T21 | Rust suppression of flowering and fruit: site pressure in both hooks; mothers weighted by fruit success; optional flowering driver `damage` (crown damage D); census columns; mature count unchanged, a prospective "will never reproduce" test left for later |
