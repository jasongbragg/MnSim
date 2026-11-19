# Constants inventory (PLAN.md Phase 2)

Status: **reviewed 2026-11-19**; decisions in §7 (methodological notes T4); **§3–4 carried out
2026-11-19** (§8; notes T9). Lists every hard-coded constant in `plant_ibm/` that a second
species could need to change, with `file:line`, and proposes a layer for each. Line numbers are for
commit `dbbde63`.

## How to read this

**Layers** (PLAN.md Phase 2), plus the categories needed for what does not fit them:

| Code | Layer | Where it would live |
|---|---|---|
| **B** | Species biology | `species/<sp>.yaml` |
| **S** | Site inputs | site table (`mn_site_params.csv` and the scripts that build it), rust anchor |
| **A** | Assessment settings | `iucn_assess_spec()` |
| **V** | Validation / calibration data | observations a species is checked or fitted against (for new species: CSV `calibration_target` rows). *Proposed fourth layer*: the Pegg values fit none of the three. |
| **M** | Model structure | code. A rule, not a number. Where a species might need a different rule, that is a Phase 2 interface or Phase 3 item, not YAML. |
| **N** | Numerical / run setting | code (floors, caps, seeds, cache sizes) |
| **I** | Infrastructure | paths, URLs, output folders, plot styling |

**Varies by species?**
- **yes**: must differ for a new species.
- **maybe**: a sensible shared default, but a species might need its own value.
- **no**: shared by all species (IUCN rules, numerics, regional fire data).

Not inventoried line by line:
- plot styling (colours, sizes);
- `diagnostics/` and `talk/` (both source `params.R`; check them if they are kept);
- `tests/`;
- `parlib/` (see `docs/phase1_2_review.md`).

---

## 1. The params list (`plant_ibm/params.R`)

Everything here is already a parameter and moves to the YAML as a block. Phase 2's point is that
**the list mixes layers**. Entries that are not species biology are marked, because they should not
end up in `m_nodosa.yaml` as if they were.

| Line | Entry | Value | Layer | Species? | Note |
|---|---|---|---|---|---|
| 24–31 | `RUST_ANCHOR` (site, lon, lat, mu, conc, cap) | Yarringully, 0.70, 20, 0.95 | **S** | yes | Anchor is per species (CLAUDE.md "Rust pressure"). `conc` is "not yet constrained by data". |
| 36 | `N0` | 10000 | N | maybe | Stand size: a run setting, but `R_max`, `K_half` and the two `half_sat`s are scaled to it (notes T1). |
| 37–38 | `n_years`, `seed` | 200, 42 | N | no | Run settings; every runner overrides them. |
| 41–42 | `weibull_k`, `weibull_lambda` | 1, 60 | B | yes | |
| 44–49 | `senescence_dose_response` | 0.1, sigmoid, 40, 8 | B | yes | |
| 51–58 | `juv_decline_dose_response`, `juv_decline_age_*` | 0.5, saturating, 5000, 2; 2, 2 | B | yes | `half_sat` is a count of canopy plants: scaled to stand size. |
| 61–62 | `R_max`, `K_half` | 1000, 4000 | B | yes | Counts of flowering plants: scaled to stand size. |
| 64–69 | `shade_dose_response` | off (0) | B | yes | |
| 72–73 | `age_first_flower_mean`, `_sd` | 5, 1.5 | B | yes | |
| 76 | `selfing_rate` | 0.10 | B | yes | |
| 81 | `fire_years` | `integer(0)` | scenario | no | Set by runners. |
| 82 | `fire_prob_annual` | `1/50` | **S** | yes | Replaced per cell by `build_site_params`. YAML needs a defined form for `1/50` (PLAN Phase 2). |
| 88–91 | `resprout_yrs_base`, `resprout_delay_dist`, `resprout_recovery`, `resprout_prob_recovery` | 2, "fixed", "countdown", 1/3 | B (+M choice) | yes | `_dist` and `_recovery` choose a rule. `1/3` is used only by "geometric". |
| 92 | `fire_kill_prob` | 0.35 | B | maybe | Legacy fire mode only (matrix tests). |
| 93 | `fire_p_fimp` | 0.6 | **S** | yes | Replaced per cell. |
| 98 | `fire_extent_mode` | "stand" | M choice | maybe | |
| 99 | `fire_size_sd` | 0.146 | **S** | maybe | Range-wide median, used where a cell has < 2 fire years. |
| 100 | `fire_stand_cover` | 0.9 | S/B | maybe | Unburnt patches in a stand: depends on habitat (rainforest margins vs heath). |
| 101–103 | `fire_kill_scalar`, `_half_sat`, `_hill` | 0.30, 2, 2 | B | yes | |
| 109 | `rust_start_year` | `Inf` | A / scenario | no | Set by runners from `year_rust`. |
| 110 | `rust_pressure` | 1 | B | maybe | Multiplier on the pathogen; overlaps `rust_dose_response$age_peak`. |
| 111–112 | `microclim_alpha`, `_beta` | 20, 5 | B | yes | Assumed the same at all sites (`runner/site_params.R:67`). |
| 115–116 | `annual_env_alpha`, `_beta` | computed from `RUST_ANCHOR` (14, 6) | **S** | yes | Replaced per cell. Computed in code: YAML needs a defined form (PLAN Phase 2). |
| 120 | `rust_site_sd` | 0 | S | maybe | Between-site spread (notes V6). |
| 122–145 | `rust_dose_response` (age_*, resprout, delay) | see file | B | yes | `delay$max_effect` 1.79 is derived from `0.8 * 0.70` (microsite mean × anchor mu), so it should be recomputed if either changes. |
| 147–159 | `rust_flower_`, `rust_recruit_dose_response` | off (0) | B | yes | |
| 165–167 | `resist_locus_effect`, `_dominance`, `_freq0` | 1, 1, 0 | B | yes | Resistance data are private (CLAUDE.md rule 1). Placeholders only. |

**Fields read by the code but absent from `params.R`**:
- `archive_dead` (`R/simulate.R:46`, default on);
- `ricker_r` and `K` (`R/recruitment.R:23`, Ricker, not wired in);
- `pegg_site_cell` and `rust_site_z`, written by runners.

The YAML schema should declare them or drop them.

---

## 2. Model code (`plant_ibm/R/`)

The core modules read almost everything from `params`. What is left is mostly structure and numerics.

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| R1 | `individuals.R:30`, `:74` | `pmax(round(rnorm(...)), 1)` | Age at first flowering: rounded Normal, at least 1 | M | maybe | Keep. Duplicated in two places: one helper. |
| R2 | `individuals.R:33` | `rgeom(n, 1/weibull_lambda)` | Founder ages: geometric, from the background hazard only | M | maybe | Keep; spin-up washes it out. A long-lived species with senescence needs a longer spin-up (A1). |
| R3 | `individuals.R:84` | age 0 | Recruits enter at age 0 after mortality (notes E2) | M | no | Keep |
| R4 | `mortality.R:21` | 0.01 | Age floor in `weibull_hazard` (diagnostic only) | N | no | Keep |
| R5 | `mortality.R:86` | `alpha/(alpha+beta)` | Resprout delay uses the site's *mean* pressure, not the year's draw | M | no | Keep |
| R6 | `mortality.R:170` | `age < AFR & !resprout` | Who counts as a juvenile | M | no | Keep |
| R7 | `mortality.R:218` | union formula | Hazards combine as independent causes | M | no | Keep |
| R8 | `fire.R:57` | 0.9 | Cap on Beta variance of fire size | N | no | Keep |
| R9 | `fire.R:68`, `:77`, `:125` | "plant", 1, Poisson | Fallbacks when `fire_extent_mode`, `fire_stand_cover`, `resprout_delay_dist` are absent | N | no | Remove once the YAML loader requires these fields. |
| R10 | `fire.R:84`, `:99` | `ceiling + (1-ceiling)*hill` | Kill curve: certain death at age 0, falling to a ceiling | M | maybe | Keep; candidate for the fire-response interface (PLAN Phase 2) |
| R11 | `fire.R:128` | 1 | Resprouting lasts at least 1 year | M | maybe | Keep; flag for clonal/resprouting species (Phase 3) |
| R12 | `recruitment.R:40` | `s/(n-1+s)` | Selfing probability per offspring | M | no | Keep |
| R13 | `recruitment.R:64–75` | `age >= AFR`, progress weight | Canopy = plants of flowering age; resprouters weighted by recovery | M | maybe | Keep |
| R14 | `recruitment.R:79` | 1 | Shade suppression capped at 1 | N | no | Keep |
| R15 | `recruitment.R:118` | Poisson | Recruit count is Poisson around Beverton–Holt | M | maybe | Keep; recruitment interface (PLAN Phase 2) |
| R16 | `simulate.R:63–67` | 0 | No pathogen before `rust_start_year` | M | no | Keep |
| R17 | `simulate.R:163` | 1000 | Archive dead rows above this count | N | no | Keep |
| R18 | `simulate.R:150` | 10 | Verbose print every 10 years | I | no | Keep |
| R19 | `genetics.R:26–28`, `:51–53`, `:84` | HWE, 0.5, clamp 1 | Diploid Mendelian inheritance; resistance capped at 1 | M | no | Keep |
| R20 | `census.R:30` | `age >= AFR` | Mature individuals for IUCN (CLAUDE.md "IUCN framing") | A | no | Keep (fixed by CLAUDE.md) |
| R21 | `matrix_model.R:47`, `:170` | 100, 50 | Leslie max age, burn-in for λ estimate | N | no | Keep (test support) |

---

## 3. Runners (`plant_ibm/runner/`)

### 3.1 Assessment settings

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| A1 | `iucn_a4.R:164–171` | 1600, 1810–1910, 2010, 2026, 2210, 11, 100 | `iucn_assess_spec()` defaults: spin-up, gen-time window, rust, present, end, averaging, A4 cap | A | **maybe** | Spin-up length (410 yr) and the gen-time window (100 yr) must be long relative to the species' generation time; a longer-lived species needs longer ones. Make them species-level defaults (YAML `assessment:` block) or derive from life history. |
| A2 | `iucn_a4.R:286–293` | same values | `run_iucn_a4()` defaults | A | — | **Duplicate of A1**: default to `iucn_assess_spec()` |
| A3 | `iucn_a4.R:207` | 11, 100 | `assess_a4()` defaults | A | — | **Duplicate of A1** |
| A4 | `range_ensemble.R:30–36` | same values | `YEAR_*`, `AVG_YEARS` | A | — | **Duplicate of A1** |
| A5 | `iucn_a4.R:167`, `range_ensemble.R:33` | 2010 | Rust arrival | A | maybe | Shared by species in NSW; may differ by region (later arrival further north). If so it becomes a site input. |
| A6 | `iucn_a4.R:81` | "A4" | Criterion | A | maybe | CLAUDE.md fixes A4 for M. nodosa; other species may use A3 (PLAN Phase 6). |
| A7 | `iucn_a4.R:191–193` | 30/50/80 | Category thresholds (A2–A4) | A | no | Keep; one constant |
| A8 | `iucn_a4.R:417–419`, `:786`, `:845`; `range_ensemble.R:197` | 30/50/80 | Same thresholds in summaries and plots | A | — | **Duplicates of A7** |
| A9 | `iucn_a4.R:223` | `max(10, 3*T_g)` | 10 yr or 3 generations (IUCN) | A | no | Keep |
| A10 | `iucn_a4.R:221` | 100 | Window when `T_g` is NA | A | no | Keep |
| A11 | `iucn_a4.R:226` | `year_rust - 1` | Reference window ends the year before rust | A | no | Keep |
| A12 | `iucn_a4.R:149` | maternal age, founders and live mothers excluded | Generation-length definition | A | no | Keep |
| A13 | `iucn_a4.R:306` | 0.4 | Warn if fixed fire years cover < 40 % of the run | N | no | Keep |

### 3.2 Pegg et al. 2025 validation (M. nodosa only)

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| V1 | `pegg_validation.R:187–190` | 1910, 2010, 2020, 2022 | Spin-up, rust, fire, survey years | V | **yes** | Per-species validation config (YAML `validation:` block or a CSV) |
| V2 | `pegg_validation.R:191`, `:493` | 0.48, n = 50 | Observed fraction dead, plot size | V | **yes** | Same; should cite the CSV/notes row. Duplicated twice. |
| V3 | `pegg_validation.R:206–213` | `fire_p_fimp` 1, `size_sd` 0, cover 1 | Every surveyed tree burnt | V | yes | Same config. Duplicated in `pegg_sites.R:76–79`. |
| V4 | `pegg_validation.R:193` | "anchor" | Validation runs at the rust anchor cell | V/S | yes | Same config |
| V5 | `pegg_validation.R:80`; `pegg_sites.R:46` | "v3-anchor-site", "v1-interval-mixture" | Cache version tags | I | no | Keep |
| V6 | `pegg_sites.R:45` | `runner/data/pegg2025_mnodosa_plots.csv` | Plot observations | V | yes | Per-species path in the validation config |
| V7 | `pegg_sites.R:67`, `:85` | 1910, 2010, 2020 | Pegg years again | V | — | **Duplicate of V1** |
| V8 | `pegg_sites.R:84` | `n_post = 3` | Post-fire years simulated | V | yes | Validation config (follows the survey length) |
| V9 | `pegg_sites.R:109`, `:118` | 12 | Months per year | N | no | Keep |
| V10 | `pegg_sites.R:71`, `:117`, `:192`, `:200` | 1000003, 1e-12, 1e-300, 4000 | Seed offset, log floors, predictive draws | N | no | Keep |

### 3.3 Range-wide runs and site overrides

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| S1 | `site_params.R:8` | `spatial/data/mn_site_params.csv` | Site table path | S | **yes** | Per species (`species/<sp>.yaml` → path, or `spatial/data/<sp>/`) |
| S2 | `site_params.R:28` | nearest cell | Anchor-cell fallback for old tables | S | no | Keep |
| S3 | `site_params.R:40` | "wild", "all" | Fire regimes | S | no | Keep |
| S4 | `site_params.R:57` | 1 | Fire probability capped at 1 | N | no | Keep |
| S5 | `site_params.R:67–75` | — | microclim constant across sites; resistance hook commented out | B / S | yes | Resistance per site needs private data: placeholder hook only |
| S6 | `range_ensemble.R:43` | −28.16 | `QLD_LAT`, labels cells NSW/Qld | S | no | **Duplicate** of `QLD_BORDER_LAT` (`spatial/fire/fire_config.R:74`): one constant |
| S7 | `range_ensemble.R:48–54` | `max_effect` U(0.32, 0.78) | `DEFAULT_PARAM_DRAWS`: prior range from the Pegg fit | B | **yes** | YAML prior block (rule 5: `prior: true` with range and reason). Use the `parlib/ensemble/.../param_ranges.R` format if adopted (PLAN Phase 2). |
| S8 | `range_ensemble.R:89`, `:114` | 20261022 | Master seed | N | no | Keep |
| S9 | `range_ensemble.R:91`, `:92` | `weight_n_raw`, 99999 | Cell weights; seed range | S / N | no | Keep |
| S10 | `range_ensemble.R:151` | 1950:2130 | Years kept for trajectory plots | I | no | Keep |

### 3.4 Other runners

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| O1 | `ensemble.R:33`, `:78–82`; `morris.R:43–46`; `summarise.R:29` | seed 42, `r` 20, `levels` 5, 50 summary years | Sensitivity-analysis settings | N | no | Keep |
| O2 | `ensemble.R:93`; `morris.R:83` | `n_years` 5, `N0` 50 | Dummy run to get column names | N | no | Keep |
| O3 | `params_utils.R:70` | 10000 | Hash fallback truncation | N | no | Keep |
| O4 | `main.R:35–75` | 1710, 300 yr, fire years 1750…1995, rust 2010 | Demo script scenario | I | — | Not used by any runner. Archive or leave as an example. |

---

## 4. Spatial pipeline (`plant_ibm/spatial/`)

These build the site table. Most constants are regional (NSW/Qld fire records) and would be shared by
every species in the same region. What is species-specific is the occurrence query, the range extent
and the file names.

| ID | Where | Value | What | Layer | Species? | Proposal |
|---|---|---|---|---|---|---|
| P1 | `01_get_melaleuca_occurrences.R:41`, `:138`, `:294` | "Melaleuca nodosa" | Species queried from ALA | B (identity) | **yes** | YAML `species:` name; read by the script |
| P2 | `01_…:121–122`, `:174`, `:260–261` | lon 148–155, lat −38 to −24 | Range bounding box | S | **yes** | Per species; one constant (now three copies) |
| P3 | `01_…:39`, `:40` | 1970, 5000 m | Occurrence filters (year, coordinate uncertainty) | S | maybe | Shared default |
| P4 | `01_…:37` | placeholder email | ALA login | I | no | Keep (placeholder, not a real address) |
| P5 | `01_…:38`; `02_…:58`; `fire/fire_config.R:13`; `rust_pressure.R:47`, `:59` | 0.1 | Grid cell size (degrees) | S | no | **Five copies**: one constant |
| P6 | `01_…:104`, `:164`, `:199`, `:229`; `02_…:57`, `:141`, `:188`, `:249`; `02c_…:18`; `fire/fire_config.R:9` | `spatial/data/mn_*.csv` | Species-named data files | S / I | **yes** | Per-species data folder |
| P7 | `02_…:56`; `02c_…:20` | `/home/jgb/projects/mrsdm/.../AP_suitability_current.tif` | Rust SDM (pathogen, not host) | S | no | Shared; absolute local path. Make configurable (env var, as `NPWS_SHP`). |
| P8 | `02_…:196–201` | `n = 200`, ±0.5° | KDE for cell weights | S | no | Keep |
| P9 | `rust_pressure.R:35` | `1/(conc*10)` | Floor on mean pressure (alpha ≥ 0.1) | N | no | Keep |
| P10 | `rust_pressure.R:61`, `:69` | −1e30, > 1.5 → /1000 | SDM nodata and unit handling | N | no | Keep |
| P11 | `rust_pressure.R:78` | 0.15°, 0.10 | Anchor neighbourhood check | N | no | Keep |
| P12 | `fire/fire_config.R:8`, `:10–12`, `:16–17` | paths, Qld service URLs | Fire data sources | I | no | Keep |
| P13 | `fire/fire_config.R:18–19`, `:24–27` | 1987–2025; Qld bbox; NSW 1987–2024, Qld 2014–2024 | Fire record years and Qld download extent | S | no (bbox: maybe) | Keep. The Qld bbox covers M. nodosa cells plus a margin: a more northern species needs a wider one. |
| P14 | `fire/fire_config.R:21` | 3577 | Area CRS (Albers) | N | no | Keep |
| P15 | `fire/fire_config.R:30–35` | 0.10, Gamma(2, 100), 0.005, 0.5, 0.05 | Fire-year threshold, fallback prior, floors and caps | S | no | Keep (notes D2) |
| P16 | `fire/fire_config.R:51–64` | regimes and sensitivity runs | Fire scenarios | S | no | Keep |
| P17 | `fire/fire_config.R:74–75`; `fire/R/grid.R:19–20` | −28.16, 1.0; state boxes | Border latitude, band, state bounding boxes | S | no | Border latitude duplicates S6 |
| P18 | `fire/R/qld_download.R:20–22` | 500, 0.00005, 5 | Download batching | I | no | Keep |
| P19 | `02b_fire_prior_spatial.R` (whole file) | — | Superseded fire prior (header says so) | — | — | Archive. It overwrites `mn_site_params.csv` if run. |

---

## 5. Summary

**Must become per-species (beyond the params list):**
- rust anchor (`RUST_ANCHOR`);
- site table path and the `mn_*` data files (S1, P6);
- occurrence query and range box (P1, P2);
- parameter-uncertainty prior (S7);
- the whole Pegg validation set-up (V1–V4, V6, V8).

The validation set-up is the clearest case for a fourth layer, **V**.

**Probably per-species:**
- spin-up and generation-time windows (A1), because they scale with generation time;
- rust arrival year (A5), if it differs by region;
- criterion (A6).

**Duplicates to remove in Phase 2:**
- phasing years: A1–A4, four copies;
- category thresholds: A7–A8, six copies;
- Pegg years and fire overrides: V1, V3, V7;
- grid size: P5, five copies;
- border latitude: S6, P17;
- range box: P2, three copies;
- age at first flowering draw: R1.

**Stays in code:** structure (M) and numerics (N). Structural rules that a new species might
challenge go to the Phase 2 interfaces or Phase 3:
- the kill curve (R10);
- the minimum resprout time (R11);
- Poisson recruitment (R15);
- founder ages (R2).

**Mixed layers in `params.R`:** `fire_prob_annual`, `fire_p_fimp`, `fire_size_sd`,
`annual_env_alpha/beta` and `RUST_ANCHOR` are site inputs carried in the species list as defaults for
single-site runs. Proposal: in `m_nodosa.yaml`, keep them under a separate `site_defaults:` block, so
the biology block holds only biology. `N0`, `n_years`, `seed`, `fire_years` and `rust_start_year` are
run settings. They need defaults, but not in the species file.

## 6. For the "done when" grep (PLAN Phase 2)

After the refactor, these should find nothing outside config, tests and docs:
```
grep -rnE '1600L|1810L|1910L|2010L|2026L|2210L|2020L|2022L|0\.48|nodosa|mn_site_params|mn_grid|Yarringully|-28\.16|GRID_DEG *<-' plant_ibm/R plant_ibm/runner plant_ibm/spatial
```

## 7. Decisions (2026-11-19)

1. **Validation is a fourth layer, V**, kept separate from the species biology: its own per-species
   file or block. This covers the Pegg years, the observed values, the scenario overrides and the plot
   data (§3.2).
2. **Landscape properties go in a separate file**, shared by all species: values that are the same
   whatever species occupies a cell. This covers:
   - the fire regime defaults and the fire pipeline settings (P13–P17);
   - the grid size (P5) and the rust SDM (P7).

   The rust anchor (`RUST_ANCHOR`) is calibrated to a species' own field data, so it stays with the
   species. The same goes for the occurrence query, the range box and the `mn_*` site table.
3. **Spin-up and generation-time windows (A1) are species settings.** Where they are set, a comment
   must say they need to be long relative to T_g. T_g is emergent (it comes out of each run), so it
   cannot set other values in advance. Instead, a run warns when a window is short relative to that
   run's T_g. The M. nodosa values can set the thresholds: they were judged adequate. With T_g ≈ 23 yr
   (test-scale A4 run, notes T1):
   - the gen-time window is 100 yr ≈ 4.4 T_g;
   - spin-up before that window is 210 yr ≈ 9 T_g;
   - spin-up before rust is 410 yr ≈ 18 T_g.

   Proposed thresholds, to be confirmed when built: warn if the gen-time window is < 3 T_g, or if
   spin-up before it is < 5 T_g.
4. **Rust arrival is 2010 for all species** (A5). It stays an assessment setting, defined once.
5. **`main.R` (O4) and `02b_fire_prior_spatial.R` (P19) archived** to `plant_ibm/archive/`.

## 8. Carried out (2026-11-19, notes T9)

- **Assessment** (A1–A5, A7–A8): one source. Shared settings and thresholds are in
  `config/assessment.yaml`; spin-up and gen-time window are in the species file's `assessment`
  block. `iucn_assess_spec()`, `run_iucn_a4()`, `assess_a4()`, `range_ensemble.R` and the plots all
  read `IUCN_ASSESSMENT` / `IUCN_CATEGORIES`.
- **Validation** (V1–V4, V6–V8): `species/m_nodosa_validation.yaml`, read as `PEGG`.
- **Range-wide** (S1, S6, S7): the site table path comes from the default species. The border
  latitude comes from `config/landscape.yaml`. The draws are in the species file's
  `ensemble_draws` block.
- **Spatial** (P1, P2, P5, P6, P7, P17):
  - the species name and range box come from the species file;
  - grid size, border latitude and SDM path come from `config/landscape.yaml`;
  - the data files are in `spatial/data/m_nodosa/`.
- **R1:** `draw_age_first_flower()`.
- **Not moved:** the items marked Keep. Two remain open (notes T9): `QLD_BBOX` (P13) and the
  grid/cell ids for a second species.
- **§6 grep:** no hits on code lines outside config, tests and docs. The remaining hits are:
  - comments that describe M. nodosa, mostly in the Pegg runners, which are M. nodosa validation code;
  - the `GRID_DEG <-` lines, which now read the configuration;
  - the `QLD_BBOX` comment.
