# diagnostics/diag_rust_pressure.R
#
# Biological-intuition diagnostics for how myrtle rust pressure and
# genetic resistance jointly determine mortality risk. These are pure
# functions (no side effects, no execution) so they can be reused both
# from an interactive session and from future testthat tests of
# biological behaviour -- see run_rust_pressure_diagnostics.R for the
# runnable script that calls these and produces numbers/figures.
#
# Intended use: communicating to collaborators (or readers of an
# article) exactly what the genetic architecture and rust_pressure
# parameterisation imply, without anyone having to read the simulation
# code.

# ----------------------------------------------------------------------
# 1. Genotype -> resistance score table, and population-level summary
# ----------------------------------------------------------------------

#' Enumerate every possible genotype combination for the configured
#' architecture (3^n_loci rows, one per locus per individual: 0/1/2),
#' with its resist_score and its expected Hardy-Weinberg frequency
#' given params$resist_freq0. Loci are assumed independent (matching
#' the independence assumption in genetics.R::init_resist_gt()).
genotype_score_table <- function(params) {
  n_loci    <- length(params$resist_locus_effect)
  loci_grid <- expand.grid(rep(list(0:2), n_loci))
  colnames(loci_grid) <- paste0("locus", seq_len(n_loci))
  gt_matrix <- as.matrix(loci_grid)

  scores <- resist_score_from_gt(gt_matrix, params)

  freq <- rep(1, nrow(gt_matrix))
  for (j in seq_len(n_loci)) {
    p   <- params$resist_freq0[j]
    g   <- gt_matrix[, j]
    f_g <- ifelse(g == 2L, p^2, ifelse(g == 1L, 2 * p * (1 - p), (1 - p)^2))
    freq <- freq * f_g
  }

  data.frame(loci_grid, resist_score = scores, expected_freq = freq)
}

#' Population-level resistance summary under Hardy-Weinberg expectation
#' at the configured resist_freq0 -- makes "what fraction of individuals
#' are genetically resistant" a concrete, reportable number.
summarize_resistance <- function(params, threshold = 0.5) {
  tab <- genotype_score_table(params)
  list(
    frac_fully_resistant  = sum(tab$expected_freq[tab$resist_score >= 0.999]),
    frac_any_protection   = sum(tab$expected_freq[tab$resist_score > 0]),
    frac_above_threshold  = sum(tab$expected_freq[tab$resist_score >= threshold]),
    mean_population_score = sum(tab$expected_freq * tab$resist_score),
    table                  = tab
  )
}

# ----------------------------------------------------------------------
# 2. Representative resistance classes for plotting
# ----------------------------------------------------------------------

#' Reduce the (potentially large, 3^n_loci) genotype space to four
#' classes that stay meaningful regardless of n_loci or architecture:
#'   susceptible      - SS at every locus            (worst case)
#'   heterozygous      - Rs at every locus            (carrier, every locus)
#'   resistant         - RR at every locus            (best case)
#'   population_mean  - the HWE population-average individual at the
#'                       configured resist_freq0
#' For a single dominant locus, "heterozygous" and "resistant" coincide
#' -- seeing that on the figure IS the demonstration that one copy is
#' enough.
representative_resist_scores <- function(params) {
  n_loci <- length(params$resist_locus_effect)

  gt_s <- matrix(0L, nrow = 1, ncol = n_loci)
  gt_h <- matrix(1L, nrow = 1, ncol = n_loci)
  gt_r <- matrix(2L, nrow = 1, ncol = n_loci)

  pop_mean <- summarize_resistance(params)$mean_population_score

  c(
    susceptible      = resist_score_from_gt(gt_s, params),
    heterozygous     = resist_score_from_gt(gt_h, params),
    resistant        = resist_score_from_gt(gt_r, params),
    population_mean = pop_mean
  )
}

# ----------------------------------------------------------------------
# 3. Mortality vs rust pressure, by resistance class
# ----------------------------------------------------------------------

#' Build a long data frame of realised annual mortality probability
#' across a range of rust_pressure values, for each representative
#' resistance class, for one hazard component ("juvenile" or
#' "resprout"). Uses whatever functional form is configured for that
#' stage in params$rust_dose_response -- so this figure automatically
#' reflects a switch from "linear" to "saturating"/"sigmoid"/etc, with
#' no changes needed here. include_baseline = TRUE adds the age/density
#' baseline hazard so curves show absolute annual mortality risk, not
#' just the rust-attributable increment -- this is normally what you
#' want for a figure aimed at collaborators, since "how much extra
#' risk" is less immediately interpretable than "what's the actual
#' annual death probability".
mortality_vs_pressure <- function(params, pressure_seq = seq(0, 2, by = 0.1),
                                   hazard_type = c("juvenile", "resprout"),
                                   include_baseline = TRUE,
                                   baseline_age = NULL, baseline_N = NULL) {
  hazard_type <- match.arg(hazard_type)
  scores <- representative_resist_scores(params)
  stage  <- switch(hazard_type, juvenile = "juv", resprout = "resprout")

  df <- expand.grid(rust_pressure = pressure_seq,
                     class = factor(names(scores), levels = names(scores)))
  df$resist_score <- scores[as.character(df$class)]
  df$rust_extra   <- (1 - df$resist_score) *
                       dose_response(df$rust_pressure, params$rust_dose_response[[stage]])

  if (include_baseline) {
    age <- if (!is.null(baseline_age)) baseline_age else
      switch(hazard_type,
        juvenile = 1,
        resprout = round(params$age_first_flower_mean) + 2)
    N <- if (!is.null(baseline_N)) baseline_N else params$K
    df$baseline <- weibull_hazard(age, params) + dd_hazard(N, params)
  } else {
    df$baseline <- 0
  }

  df$p_death     <- pmin(df$baseline + df$rust_extra, 1)
  df$hazard_type <- hazard_type
  df
}

#' Plot mortality probability vs rust pressure, one line per resistance
#' class, juvenile and resprout hazard side by side. This is the figure
#' meant for sharing with collaborators or a manuscript methods
#' section. Uses ggplot2 if available, otherwise a base R fallback.
plot_mortality_vs_pressure <- function(params, pressure_seq = seq(0, 2, by = 0.1),
                                        save_path = NULL) {
  df_juv <- mortality_vs_pressure(params, pressure_seq, "juvenile")
  df_rsp <- mortality_vs_pressure(params, pressure_seq, "resprout")
  df <- rbind(df_juv, df_rsp)
  df$hazard_type <- factor(df$hazard_type, levels = c("juvenile", "resprout"),
                            labels = c("Juvenile mortality",
                                       "Post-fire resprout mortality"))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = rust_pressure, y = p_death, color = class)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~ hazard_type) +
      ggplot2::labs(x = "Rust pressure (1 = reference / humid coastal)",
                     y = "Annual mortality probability",
                     color = "Resistance class",
                     title = "Mortality risk vs rust pressure, by genetic resistance class") +
      ggplot2::theme_minimal()
    if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 9, height = 4.5, dpi = 150)
    return(p)
  }

  message("ggplot2 not available -- producing base R panel plot instead.")
  if (!is.null(save_path)) png(save_path, width = 1500, height = 750, res = 150)
  on.exit(if (!is.null(save_path)) dev.off(), add = TRUE)

  par(mfrow = c(1, 2))
  classes <- levels(df$class)
  cols    <- c("#1B7837", "#5AAE61", "#762A83", "#B35806")[seq_along(classes)]

  for (ht in levels(df$hazard_type)) {
    sub <- df[df$hazard_type == ht, ]
    plot(NULL, xlim = range(sub$rust_pressure), ylim = c(0, max(sub$p_death, 0.01)),
         xlab = "Rust pressure", ylab = "Annual mortality probability", main = ht)
    for (i in seq_along(classes)) {
      s <- sub[sub$class == classes[i], ]
      lines(s$rust_pressure, s$p_death, col = cols[i], lwd = 2)
    }
    legend("topleft", legend = classes, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8)
  }
  invisible(df)
}

# ----------------------------------------------------------------------
# 4. Comparing candidate functional forms, before committing one
# ----------------------------------------------------------------------

#' Overlay several candidate dose_response() shapes on one plot, all
#' NORMALISED to max_effect = 1, so only the SHAPE is being compared --
#' not the absolute magnitude (which is a separate, per-stage decision
#' already in params.R). Use this to decide, e.g., whether juvenile
#' mortality should be "linear" or "saturating" before editing
#' params$rust_dose_response$juv.
compare_dose_response_shapes <- function(pressure_seq = seq(0, 3, by = 0.05),
                                          half_sat = 1, hill = 2,
                                          save_path = NULL) {
  forms <- c("linear", "power", "saturating", "sigmoid", "threshold")
  df <- do.call(rbind, lapply(forms, function(f) {
    cfg <- list(max_effect = 1, form = f, half_sat = half_sat, hill = hill)
    data.frame(rust_pressure = pressure_seq,
               effect = dose_response(pressure_seq, cfg),
               form = f)
  }))
  df$form <- factor(df$form, levels = forms)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = rust_pressure, y = effect, color = form)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
      ggplot2::labs(x = "Rust pressure", y = "Realised effect (normalised, max_effect = 1)",
                     title = "Candidate dose-response shapes",
                     subtitle = sprintf("half_sat = %.2g, hill = %.2g", half_sat, hill),
                     color = "Functional form") +
      ggplot2::theme_minimal()
    if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 7, height = 4.5, dpi = 150)
    return(p)
  }

  message("ggplot2 not available -- producing base R plot instead.")
  if (!is.null(save_path)) png(save_path, width = 1100, height = 700, res = 150)
  on.exit(if (!is.null(save_path)) dev.off(), add = TRUE)
  cols <- c("#1B7837", "#5AAE61", "#762A83", "#B35806", "#2166AC")
  plot(NULL, xlim = range(pressure_seq), ylim = c(0, max(df$effect)),
       xlab = "Rust pressure", ylab = "Realised effect (normalised)",
       main = "Candidate dose-response shapes")
  abline(v = 1, lty = 2, col = "grey60")
  for (i in seq_along(forms)) {
    s <- df[df$form == forms[i], ]
    lines(s$rust_pressure, s$effect, col = cols[i], lwd = 2)
  }
  legend("topleft", legend = forms, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8)
  invisible(df)
}

# ----------------------------------------------------------------------
# 5. Fire x rust compounding: cumulative risk across a resprout episode
# ----------------------------------------------------------------------

#' The headline fire x rust INTERACTION number. A single fire-survival
#' event isn't fully described by "resprout mortality is X% higher" --
#' rust ALSO lengthens how long an individual spends in the vulnerable
#' resprout window (delay_extra), so the elevated per-year hazard gets
#' applied for longer. Those two effects compound multiplicatively,
#' not additively, which is exactly the kind of interaction that's easy
#' to under-communicate if you only ever look at the two rust effects
#' separately.
#'
#' Returns P(an individual that survived the fire's kill-pass goes on
#' to die during the resprout window, before resuming flowering), as a
#' function of rust_pressure, per representative resistance class.
#'
#' Approximation, clearly scoped to this diagnostic only: treats
#' resprout duration as fixed at its expected value
#' (resprout_yrs_base + delay_extra) and the per-year hazard as constant
#' across that window. The simulation engine itself does NOT use this
#' approximation -- it draws both duration and yearly mortality
#' stochastically every fire event (see fire.R, mortality.R).
cumulative_resprout_mortality_risk <- function(params, pressure_seq = seq(0, 2, by = 0.1),
                                                age_for_baseline = NULL) {
  scores <- representative_resist_scores(params)
  age <- if (!is.null(age_for_baseline)) age_for_baseline else
    round(params$age_first_flower_mean) + 2
  baseline <- weibull_hazard(age, params) + dd_hazard(params$K, params)

  df <- expand.grid(rust_pressure = pressure_seq,
                     class = factor(names(scores), levels = names(scores)))
  df$resist_score <- scores[as.character(df$class)]

  resprout_extra <- dose_response(df$rust_pressure, params$rust_dose_response$resprout)
  delay_extra    <- dose_response(df$rust_pressure, params$rust_dose_response$delay)

  df$p_death_per_year  <- pmin(baseline + (1 - df$resist_score) * resprout_extra, 1)
  df$expected_duration <- params$resprout_yrs_base + (1 - df$resist_score) * delay_extra
  df$cumulative_risk   <- 1 - (1 - df$p_death_per_year) ^ df$expected_duration
  df
}

#' Plot the fire x rust compounding metric vs rust pressure, one line
#' per resistance class. This is the figure that makes the INTERACTION
#' (not just the two separate rust effects) legible to a collaborator
#' or reader.
plot_cumulative_resprout_risk <- function(params, pressure_seq = seq(0, 2, by = 0.1),
                                           save_path = NULL) {
  df <- cumulative_resprout_mortality_risk(params, pressure_seq)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = rust_pressure, y = cumulative_risk, color = class)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Rust pressure",
                     y = "P(dies before resuming flowering | survived the fire)",
                     color = "Resistance class",
                     title = "Fire x rust compounding: cumulative mortality risk during resprout recovery") +
      ggplot2::theme_minimal()
    if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 7.5, height = 4.5, dpi = 150)
    return(p)
  }

  message("ggplot2 not available -- producing base R plot instead.")
  if (!is.null(save_path)) png(save_path, width = 1200, height = 700, res = 150)
  on.exit(if (!is.null(save_path)) dev.off(), add = TRUE)
  classes <- levels(df$class)
  cols <- c("#1B7837", "#5AAE61", "#762A83", "#B35806")[seq_along(classes)]
  plot(NULL, xlim = range(df$rust_pressure), ylim = c(0, 1),
       xlab = "Rust pressure", ylab = "P(dies before resuming flowering)",
       main = "Fire x rust compounding: cumulative resprout mortality risk")
  for (i in seq_along(classes)) {
    s <- df[df$class == classes[i], ]
    lines(s$rust_pressure, s$cumulative_risk, col = cols[i], lwd = 2)
  }
  legend("topleft", legend = classes, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8)
  invisible(df)
}
