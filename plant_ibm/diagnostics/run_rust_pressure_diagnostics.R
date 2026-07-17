# diagnostics/run_rust_pressure_diagnostics.R
#
# Run this to produce the numbers and figure for communicating, to
# collaborators or in an article, exactly what the resistance genetics
# and rust_pressure parameterisation imply for mortality risk.
#
#   Rscript diagnostics/run_rust_pressure_diagnostics.R
#
# Run from the plant_ibm/ project root (paths below are relative to it).

source("R/individuals.R")
source("R/genetics.R")
source("R/mortality.R")
source("R/fire.R")
source("R/recruitment.R")
source("R/census.R")
source("R/simulate.R")
source("params.R")
source("diagnostics/diag_rust_pressure.R")

params <- get_default_params()

cat("================================================================\n")
cat(" Resistance genetics summary (Hardy-Weinberg expectation)\n")
cat("================================================================\n")
cat("Architecture: ", length(params$resist_locus_effect), "locus/loci\n")
cat("resist_locus_effect:", paste(params$resist_locus_effect, collapse = ", "), "\n")
cat("resist_dominance:   ", paste(params$resist_dominance, collapse = ", "), "\n")
cat("resist_freq0:        ", paste(params$resist_freq0, collapse = ", "), "\n\n")

summ <- summarize_resistance(params)
cat(sprintf("Fraction fully resistant (score = 1):        %.1f%%\n", 100 * summ$frac_fully_resistant))
cat(sprintf("Fraction with any protection (score > 0):     %.1f%%\n", 100 * summ$frac_any_protection))
cat(sprintf("Fraction above score = 0.5:                   %.1f%%\n", 100 * summ$frac_above_threshold))
cat(sprintf("Mean resistance score, population-wide:        %.3f\n", summ$mean_population_score))

cat("\nGenotype -> score -> expected frequency table:\n")
print(summ$table, row.names = FALSE)

cat("\nRepresentative resistance classes used in the figure below:\n")
print(round(representative_resist_scores(params), 3))

cat("\n================================================================\n")
cat(" Mortality vs rust pressure, by resistance class\n")
cat("================================================================\n")

out_path <- "outputs/rust_pressure_mortality.png"
dir.create("outputs", showWarnings = FALSE)
plot_mortality_vs_pressure(params, pressure_seq = seq(0, 2, by = 0.1), save_path = out_path)
cat(sprintf("\nFigure saved to: %s\n", out_path))

cat("\n================================================================\n")
cat(" Candidate dose-response shapes (for choosing per-stage forms)\n")
cat("================================================================\n")
shapes_path <- "outputs/dose_response_shapes.png"
compare_dose_response_shapes(half_sat = 1, hill = 2, save_path = shapes_path)
cat(sprintf("Figure saved to: %s\n", shapes_path))
cat("All curves normalised to max_effect = 1 -- this compares SHAPE only.\n")
cat("Edit params$rust_dose_response$<juv|resprout|delay>$form to change\n")
cat("which shape a given life-history stage actually uses.\n")

cat("\n================================================================\n")
cat(" Fire x rust compounding: cumulative resprout mortality risk\n")
cat("================================================================\n")
cum_path <- "outputs/cumulative_resprout_risk.png"
cum_df <- cumulative_resprout_mortality_risk(params)
plot_cumulative_resprout_risk(params, save_path = cum_path)
cat(sprintf("Figure saved to: %s\n", cum_path))

at_p1 <- cum_df[cum_df$rust_pressure == 1, ]
cat("\nAt rust_pressure = 1 (reference), P(dies before resuming flowering):\n")
for (i in seq_len(nrow(at_p1))) {
  cat(sprintf("  %-16s expected resprout duration = %.1f yr, cumulative risk = %.1f%%\n",
              at_p1$class[i], at_p1$expected_duration[i], 100 * at_p1$cumulative_risk[i]))
}
cat("\nThis is the number that makes the fire x rust INTERACTION legible:\n")
cat("it isn't just 'mortality is X% higher per year' -- the window an\n")
cat("individual is exposed to that elevated hazard is ALSO longer, and\n")
cat("the two effects compound multiplicatively across a single fire event.\n")

cat("\nReference point (rust_pressure = 1) corresponds to the max_effect\n")
cat("entries in params$rust_dose_response (params.R).\n")
cat("Scale rust_pressure down for drier/lower-disease sites, up for\n")
cat("wetter/more humid sites, when building the geographic ensemble.\n")
