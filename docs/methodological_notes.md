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
changes earlier ones (`talk/make_talk_figures.R`).

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

**B4. Cells without inputs** [open]
16 cells (1.7% of abundance weight) lack SDM, and therefore rust, values. They
are dropped and reported, not silently lost.

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
same extraction as every other cell, so it always gets exactly 0.70 whatever
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

**C4. SDM saturation** [open]
Suitability is about 1 in almost every coastal cell (median 0.99), so the
SDM separates the coastal range from drier inland cells but gives no gradient
along the coast. Rust pressure is therefore about 0.70 almost everywhere in
the core range.

**C5. Rust effects in the model** [assumption]
Rust raises death risk, most for seedlings and young plants (declining with
age), with an extra component for resprouting plants (`max_effect`, C/V4);
it lengthens the time to recover from resprouting. Effects on flowering and
seed set are implemented but switched off (no data).
[open] The resprout delay uses rust pressure directly, not the disease
triangle.

## D. Fire

**D1. Fire history source** [decided]
NSW NPWS Fire History, 1970-2024, wildfires only. Prescribed burns are
excluded: they are typically patchy and low intensity, while the model's fire
event is a landscape-scale disturbance (validated against the 2019-20 Black
Summer wildfire). A cell-year counts as a fire year when wildfire burns at
least 10% of the cell.

**D2. Fire frequency** [calibrated]
Annual fire probability per cell = Bayesian (Gamma-Poisson) estimate from
the fire years, with a regional prior (empirical Bayes). A cell with no
recorded wildfire gets about one fire per 47 years. Fire occurrence in a run
is a random draw each year with that probability.

**D3. Fire extent and severity** [assumption, open]
Each model fire affects each plant with probability `fire_p_fimp`, the cell's
mean burnt fraction. Affected plants die with an age-dependent probability
(seedlings almost always; older plants at a ceiling of 0.30 x `fire_p_fimp`),
otherwise resprout and do not flower for a few years (Poisson, mean 3, plus
a rust-dependent delay).
[open] `fire_p_fimp` stands for both extent and intensity, so partly burnt
cells also get gentler fires. To be separated (ROADMAP Phase 0).

**D4. Queensland** [open]
No Queensland fire history yet: the about 109 Queensland cells (about 10% of
abundance weight) use the regional prior for frequency and extent. Results
are reported separately for NSW and Queensland cells until fixed (ROADMAP
Phase 1a).

## E. Demography and genetics (current values, `params.R`)

**E1.** Background mortality is a constant annual hazard (Weibull k = 1,
scale 60 yr), plus a senescence term rising with age and a juvenile
mortality term that increases with canopy density. [assumption]
**E2.** Recruitment: Beverton-Holt seedling production from flowering plants
(R_max 1000, K_half 4000), parents sampled from flowering plants; selfing rate
0.10; mean age at first flowering 5 yr. [assumption]
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
