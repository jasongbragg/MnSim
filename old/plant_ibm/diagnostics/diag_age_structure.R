# diagnostics/diag_age_structure.R
#
# Diagnostics for the age structure of a population snapshot -- typically
# the end of a run (res$individuals from run_simulation()), but works on
# any `pop` data frame with the standard alive/age/age_first_flower/
# resprout columns. Unlike diag_rust_pressure.R, these functions have NO
# dependency on the rest of the model (no genetics, no hazard functions)
# -- they only read pop's own columns -- so they can be sourced standalone,
# in any order, with no params.R needed.
#
# Life-history stage classification matches the convention documented in
# README.md exactly:
#   Juvenile            alive & age < age_first_flower & !resprout
#   Reproductive adult   alive & age >= age_first_flower & !resprout
#   Resprout             alive & resprout
# (These three are exhaustive and mutually exclusive for alive
# individuals -- there's no "alive, not resprouting, of flowering age,
# but not flowering" case in this model.)

# ----------------------------------------------------------------------
# 1. Pure data-prep function
# ----------------------------------------------------------------------

#' Long data frame of alive-individual counts by age and life-history
#' stage, for plotting or printing. If age_breaks is NULL (default),
#' one row per single year of age; otherwise ages are grouped via
#' cut(alive$age, breaks = age_breaks, right = FALSE) -- useful for a
#' long-lived population where single-year bars would be unreadable.
age_structure_table <- function(pop, age_breaks = NULL) {
  alive <- pop[pop$alive, , drop = FALSE]

  stage <- ifelse(alive$resprout, "Resprout",
             ifelse(alive$age < alive$age_first_flower, "Juvenile",
                    "Reproductive adult"))
  stage <- factor(stage, levels = c("Juvenile", "Reproductive adult", "Resprout"))

  age_grp <- if (!is.null(age_breaks)) {
    cut(alive$age, breaks = age_breaks, right = FALSE, include.lowest = TRUE)
  } else {
    factor(alive$age, levels = sort(unique(alive$age)))
  }

  tab <- as.data.frame(table(age = age_grp, stage = stage), stringsAsFactors = FALSE)
  names(tab)[names(tab) == "Freq"] <- "n"
  # table() always emits every age x stage combination, including zeros --
  # keep them (a zero-height bar segment is informative, not noise) but
  # restore both columns as properly-LEVELED factors: stringsAsFactors =
  # FALSE above silently demotes them to character, which would otherwise
  # make levels(tab$stage) NULL downstream (a real bug caught while
  # testing -- the base-R barplot path needs stages as factor levels).
  # Also preserves age_grp's level ORDER (numeric ages or chronological
  # bins), not table()'s default alphabetical-on-character ordering.
  tab$age   <- factor(tab$age,   levels = levels(age_grp))
  tab$stage <- factor(tab$stage, levels = levels(stage))
  tab[order(tab$age, tab$stage), ]
}

#' Numeric summary companion to age_structure_table() -- the figures
#' you'd quote in a methods section or report to a collaborator without
#' needing the figure itself.
summarize_age_structure <- function(pop) {
  alive <- pop[pop$alive, , drop = FALSE]
  stage <- ifelse(alive$resprout, "Resprout",
             ifelse(alive$age < alive$age_first_flower, "Juvenile",
                    "Reproductive adult"))
  stage <- factor(stage, levels = c("Juvenile", "Reproductive adult", "Resprout"))

  n_by_stage <- as.numeric(table(stage))
  by_stage <- data.frame(
    stage      = levels(stage),
    n          = n_by_stage,
    prop       = n_by_stage / nrow(alive),
    mean_age   = as.numeric(tapply(alive$age, stage, mean))[seq_along(levels(stage))],
    median_age = as.numeric(tapply(alive$age, stage, median))[seq_along(levels(stage))]
  )

  list(
    N_alive    = nrow(alive),
    mean_age   = mean(alive$age),
    median_age = median(alive$age),
    max_age    = max(alive$age),
    by_stage   = by_stage
  )
}

# ----------------------------------------------------------------------
# 2. Single-population plot
# ----------------------------------------------------------------------

#' Stacked bar chart of the age structure of one population snapshot,
#' coloured by life-history stage. Pass age_breaks for a long-lived
#' population where single-year bars would be too fine to read.
plot_age_structure <- function(pop, age_breaks = NULL, save_path = NULL,
                                title = "Age structure at end of run") {
  tab <- age_structure_table(pop, age_breaks)
  x_label <- if (is.null(age_breaks)) "Age (years)" else "Age bin"

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = age, y = n, fill = stage)) +
      ggplot2::geom_col(width = 0.9, color = "white", linewidth = 0.1) +
      ggplot2::labs(x = x_label, y = "Number of alive individuals",
                    fill = "Life stage", title = title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (is.null(age_breaks)) 0 else 45,
                                                          hjust = if (is.null(age_breaks)) 0.5 else 1))
    if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 7.5, height = 4.5, dpi = 150)
    return(p)
  }

  message("ggplot2 not available -- producing base R stacked bar plot instead.")
  if (!is.null(save_path)) png(save_path, width = 1200, height = 720, res = 150)
  on.exit(if (!is.null(save_path)) dev.off(), add = TRUE)
  .base_r_stage_barplot(tab, x_label = x_label, title = title)
  invisible(tab)
}

#' Shared base-R stacked-barplot helper (single panel). age_levels lets
#' a panel of plot_age_structure_compare() share one consistent x-axis
#' across scenarios with different observed age ranges. ylim lets
#' plot_age_structure_compare() force a shared y-axis across panels too
#' (see shared_y there for why this matters).
.base_r_stage_barplot <- function(tab, x_label, title, age_levels = NULL, ylim = NULL) {
  stages <- levels(tab$stage)
  ages   <- if (!is.null(age_levels)) age_levels else levels(tab$age)

  mat <- matrix(0, nrow = length(stages), ncol = length(ages),
                 dimnames = list(stages, ages))
  for (i in seq_len(nrow(tab))) {
    si <- match(as.character(tab$stage[i]), stages)
    ai <- match(as.character(tab$age[i]), ages)
    if (!is.na(ai)) mat[si, ai] <- tab$n[i]
  }

  cols <- c("#5AAE61", "#2166AC", "#B35806")[seq_along(stages)]
  bp_args <- list(height = mat, col = cols, border = NA,
                   las = if (is.null(age_levels)) 1 else 2,
                   cex.names = 0.7, xlab = x_label, ylab = "Number of alive individuals",
                   main = title)
  if (!is.null(ylim)) bp_args$ylim <- ylim
  do.call(barplot, bp_args)
  legend("topright", legend = stages, fill = cols, bty = "n", cex = 0.8)
}

# ----------------------------------------------------------------------
# 3. Multi-population comparison plot
# ----------------------------------------------------------------------

#' Compare age structure across several population snapshots side by
#' side -- e.g. the rust-projection end state vs. its no-rust
#' counterfactual at the same final year, or several time points within
#' one run. pop_list is a named list of `pop` data frames; names become
#' facet/panel labels.
#'
#' shared_y = TRUE (default): every panel uses the SAME y-axis scale.
#' This matters more than it sounds like it should -- with independently
#' rescaled axes (the old default here), age 0 looks like it's hitting
#' the top of the chart in EVERY panel regardless of the actual counts,
#' because each panel's tallest bar always fills its own axis. That
#' makes it impossible to tell apart "juvenile dominance didn't change"
#' from "juvenile dominance dropped a lot, but you can't see it because
#' both panels were independently rescaled to look the same height."
#' Set shared_y = FALSE only if you specifically want each panel's
#' internal shape compared on its own terms, with absolute scale
#' deliberately discarded -- e.g. comparing two populations of very
#' different sizes where total N truly isn't the point. For "did this
#' parameter change actually reduce the juvenile fraction" questions
#' (the typical use of this function), always use the default.
plot_age_structure_compare <- function(pop_list, age_breaks = NULL, save_path = NULL,
                                        title = "Age structure comparison",
                                        shared_y = TRUE) {
  if (is.null(names(pop_list)) || any(names(pop_list) == "")) {
    names(pop_list) <- paste0("scenario_", seq_along(pop_list))
  }

  tabs <- lapply(names(pop_list), function(nm) {
    t <- age_structure_table(pop_list[[nm]], age_breaks)
    t$scenario <- nm
    t
  })
  df <- do.call(rbind, tabs)
  df$scenario <- factor(df$scenario, levels = names(pop_list))
  x_label <- if (is.null(age_breaks)) "Age (years)" else "Age bin"

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = age, y = n, fill = stage)) +
      ggplot2::geom_col(width = 0.9, color = "white", linewidth = 0.1) +
      ggplot2::facet_wrap(~ scenario, scales = if (shared_y) "fixed" else "free_y") +
      ggplot2::labs(x = x_label, y = "Number of alive individuals",
                    fill = "Life stage", title = title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (is.null(age_breaks)) 0 else 45,
                                                          hjust = if (is.null(age_breaks)) 0.5 else 1))
    if (!is.null(save_path)) ggplot2::ggsave(save_path, p, width = 4 * length(pop_list) + 2,
                                              height = 4.5, dpi = 150, limitsize = FALSE)
    return(p)
  }

  message("ggplot2 not available -- producing base R panel plot instead.")
  if (!is.null(save_path)) png(save_path, width = 600 * length(pop_list), height = 720, res = 150)
  on.exit(if (!is.null(save_path)) dev.off(), add = TRUE)

  # Shared age axis across panels so bars are visually comparable
  all_ages <- levels(df$age)

  # Shared y-axis (default): find the tallest STACKED bar (sum across
  # stages at a given age) across ALL panels, use that as every panel's
  # ylim. Without this, barplot()'s automatic per-call scaling silently
  # reproduces the same misleading effect described above.
  panel_ylim <- NULL
  if (shared_y) {
    stacked_totals <- tapply(df$n, list(df$scenario, df$age), sum)
    panel_ylim <- c(0, max(stacked_totals, na.rm = TRUE) * 1.05)
  }

  par(mfrow = c(1, length(pop_list)))
  for (nm in names(pop_list)) {
    sub <- df[df$scenario == nm, ]
    .base_r_stage_barplot(sub, x_label = x_label, title = nm, age_levels = all_ages,
                           ylim = panel_ylim)
  }
  invisible(df)
}
