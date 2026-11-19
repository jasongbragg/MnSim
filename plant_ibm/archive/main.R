# main.R
#
# Workflow:
#   Phase 1 -- pre-rust calibration: run with fire only, rust switched
#              off (rust_start_year = Inf), for long enough to reach a
#              quasi-equilibrium age structure. This is the run you'd
#              tune against real demographic data for the species.
#   Phase 2a -- rust projection: continue from the calibration end-state,
#              starting in 2010 (when myrtle rust arrived in Australia),
#              with rust switched on.
#   Phase 2b -- counterfactual: continue from the SAME calibration
#              end-state, same fire regime, but rust left off. This
#              isolates the rust-attributable decline from the
#              fire-attributable decline.
#
# fire_years and rust_start_year are real calendar years throughout,
# because year0 is threaded through run_simulation() as an offset.

source("R/individuals.R")
source("R/genetics.R")
source("R/mortality.R")
source("R/fire.R")
source("R/recruitment.R")
source("R/census.R")
source("R/simulate.R")
source("params.R")

have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (have_ggplot) library(ggplot2)

# ---------------------------------------------------------------------
# Phase 1: pre-rust calibration burn-in
# ---------------------------------------------------------------------
params_calib <- get_default_params()
params_calib$n_years         <- 300
params_calib$fire_years      <- c(1750, 1820, 1880, 1930, 1970, 1995)
params_calib$rust_start_year <- Inf

res_calib <- run_simulation(params_calib, year0 = 1710, verbose = TRUE)

cat(sprintf("\nCalibration complete (1710-2009). N_alive = %d, N_IUCN = %d\n",
            tail(res_calib$census$N_alive, 1), tail(res_calib$census$N_IUCN, 1)))

init_state <- list(individuals = res_calib$individuals, resist_gt = res_calib$resist_gt)

# ---------------------------------------------------------------------
# Phase 2a: rust-impacted projection, 2010 onward
# ---------------------------------------------------------------------
params_rust <- get_default_params()
params_rust$n_years         <- 50
params_rust$fire_years      <- c(2015, 2030, 2045)
params_rust$rust_start_year <- 2010

res_rust <- run_simulation(params_rust, init_state = init_state, year0 = 2010, verbose = TRUE)

# ---------------------------------------------------------------------
# Phase 2b: counterfactual -- identical starting point and fire regime,
# rust never arrives. Isolates rust-attributable decline.
# ---------------------------------------------------------------------
params_cfact <- params_rust
params_cfact$rust_start_year <- Inf

res_cfact <- run_simulation(params_cfact, init_state = init_state, year0 = 2010, verbose = TRUE)

# ---------------------------------------------------------------------
# Compare scenarios
# ---------------------------------------------------------------------
res_rust$census$scenario  <- "rust"
res_cfact$census$scenario <- "no_rust_counterfactual"
compare <- rbind(res_rust$census, res_cfact$census)

y_last <- tail(params_rust$n_years, 1) + 2010 - 1

decline_rust  <- pct_decline(res_rust$census,  2010, y_last)
decline_cfact <- pct_decline(res_cfact$census, 2010, y_last)

cat(sprintf(
  "\nProjected decline in reproductive adults, 2010-%d:\n  rust scenario:          %.1f%%\n  counterfactual (no rust): %.1f%%\n",
  y_last, decline_rust, decline_cfact))

# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------
if (have_ggplot) {

  p1 <- ggplot(compare, aes(x = year, y = N_IUCN, color = scenario)) +
    geom_line(linewidth = 1) +
    labs(title = "Reproductive adult abundance (IUCN Criterion A count)",
         x = "Year", y = "N flowering adults", color = "Scenario") +
    theme_minimal()
  print(p1)

  p2 <- ggplot(res_rust$census, aes(x = year, y = freq_locus1)) +
    geom_line(linewidth = 1, color = "darkred") +
    labs(title = "Resistance allele frequency, locus 1 (rust scenario)",
         x = "Year", y = "Frequency of R allele") +
    theme_minimal()
  print(p2)

} else {
  message("ggplot2 not available -- falling back to base graphics.")

  plot(res_rust$census$year, res_rust$census$N_IUCN, type = "l", col = "firebrick",
       xlab = "Year", ylab = "N flowering adults",
       main = "Reproductive adult abundance (IUCN Criterion A count)",
       ylim = range(c(res_rust$census$N_IUCN, res_cfact$census$N_IUCN)))
  lines(res_cfact$census$year, res_cfact$census$N_IUCN, col = "steelblue")
  legend("topright", legend = c("rust", "no_rust_counterfactual"),
         col = c("firebrick", "steelblue"), lty = 1)

  plot(res_rust$census$year, res_rust$census$freq_locus1, type = "l", col = "darkred",
       xlab = "Year", ylab = "Frequency of R allele",
       main = "Resistance allele frequency, locus 1 (rust scenario)")
}
