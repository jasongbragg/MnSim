# plant_ibm/archive

Files kept for reference only; nothing sources them. Archived 2026-11-19
(docs/constants_inventory.md, question 5).

- `main.R`: the old two-phase example (calibration, then rust projection and
  counterfactual from `init_state`). Not used by any runner. It sources model
  files by paths relative to `plant_ibm/`.
- `02b_fire_prior_spatial.R`: superseded spatial fire prior (Option B), replaced
  by `spatial/fire/` (methodological notes D2). Do not run it: it overwrites
  `spatial/data/mn_site_params.csv`.
- `parlib_diags_params.R`: the older calibration (senescence 0.2/30,
  `fire_p_fimp` 0.8, resprout `max_effect` 0.80, rust from year 1) that sat in
  five byte-identical copies under `parlib/diags/*/params.R` (md5
  79f2efd93e55095de071907cc07ba334). Nothing used them: the diagnostics source
  `params.R`. Archived 2026-11-19 (notes T5).
