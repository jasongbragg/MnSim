# plant_ibm/archive

Files kept for reference only; nothing sources them. Archived 2026-11-19
(docs/constants_inventory.md, question 5).

- `main.R`: the old two-phase example (calibration, then rust projection and
  counterfactual from `init_state`). Not used by any runner. It sources model
  files by paths relative to `plant_ibm/`.
- `02b_fire_prior_spatial.R`: superseded spatial fire prior (Option B), replaced
  by `spatial/fire/` (methodological notes D2). Do not run it: it overwrites
  `spatial/data/mn_site_params.csv`.
