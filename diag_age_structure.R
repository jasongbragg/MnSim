source("params.R"); source("R/individuals.R"); source("R/genetics.R")
source("R/mortality.R"); source("R/fire.R"); source("R/recruitment.R")
source("R/census.R"); source("R/simulate.R")
source("diagnostics/diag_age_structure.R")

# Build a scenario -- this is your no-resistance / rust calibration
p <- get_default_params()
p$resist_freq0    <- c(0)
p$rust_start_year <- 1
p$rust_dose_response$juv$max_effect      <- 0.275
p$rust_dose_response$resprout$max_effect <- 0.206
p$rust_dose_response$delay$max_effect    <- 0
p$n_years <- 30

# Run it
res_rust <- run_simulation(p, year0 = 1, verbose = FALSE)

# 1. Numeric summary, printed to console
summarize_age_structure(res_rust$individuals)

# 2. Single-population plot, saved to disk
dir.create("outputs", showWarnings = FALSE)
plot_age_structure(res_rust$individuals, save_path = "outputs/age_structure_rust.png")

# 3. Compare against the no-rust counterfactual at the same final year
p_cfact <- p
p_cfact$rust_start_year <- Inf
res_cfact <- run_simulation(p_cfact, year0 = 1, verbose = FALSE)

plot_age_structure_compare(
  list(rust = res_rust$individuals, no_rust_counterfactual = res_cfact$individuals),
  save_path = "outputs/age_structure_compare.png"
)
