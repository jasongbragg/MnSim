# tests/regression/compare.R
#
# Functions to compare two sets of regression references (an old set and a
# candidate set), for the reference-update procedure in
# docs/reference_updates.md, step 1. Used by compare_references.R (the
# command-line report) and by test-regression.R (first_difference()).
#
# Only reads reference files; never writes them.

# Where two objects first differ, as a readable message, or NULL if identical.
# For data frames it names the first differing column and row (with the year,
# if there is a year column). In lists the census is checked first, so the
# message gives the first year that differs.
first_difference <- function(new, ref, path = "output") {
  if (identical(new, ref)) return(NULL)
  if (is.data.frame(ref) && is.data.frame(new)) {
    if (!identical(names(new), names(ref)))
      return(sprintf("%s: columns differ (new: %s)", path, paste(setdiff(names(new), names(ref)), collapse = ", ")))
    if (nrow(new) != nrow(ref))
      return(sprintf("%s: %d rows, reference %d", path, nrow(new), nrow(ref)))
    for (col in names(ref)) {
      a <- new[[col]]; b <- ref[[col]]
      if (!identical(a, b)) {
        i <- which(!(a == b | (is.na(a) & is.na(b))) | xor(is.na(a), is.na(b)))[1]
        yr <- if ("year" %in% names(ref)) sprintf(" (year %s)", ref$year[i]) else ""
        return(sprintf("%s$%s: first difference at row %s%s: new %s, reference %s",
                       path, col, i, yr, format(a[i]), format(b[i])))
      }
    }
    return(sprintf("%s: attributes differ", path))
  }
  if (is.list(ref) && is.list(new)) {
    if (!identical(names(new), names(ref)))
      return(sprintf("%s: elements differ (new only: %s; reference only: %s)", path,
                     paste(setdiff(names(new), names(ref)), collapse = ", "),
                     paste(setdiff(names(ref), names(new)), collapse = ", ")))
    for (k in order(names(ref) != "census")) {
      d <- first_difference(new[[k]], ref[[k]],
                            paste0(path, "$", if (is.null(names(ref))) k else names(ref)[k]))
      if (!is.null(d)) return(d)
    }
    return(sprintf("%s: attributes differ", path))
  }
  sprintf("%s: new %s, reference %s", path,
          paste(utils::head(format(new), 5), collapse = " "),
          paste(utils::head(format(ref), 5), collapse = " "))
}

# Reads every reference file in a directory: list(scenarios = list(name = rds
# contents), stats = summary_stats.rds contents or NULL, dir).
load_reference_set <- function(dir) {
  if (!dir.exists(dir)) stop("no such directory: ", dir)
  files <- setdiff(list.files(dir, pattern = "\\.rds$", full.names = TRUE),
                   file.path(dir, "summary_stats.rds"))
  sc <- lapply(files, readRDS)
  names(sc) <- vapply(sc, `[[`, "", "name")
  stats_path <- file.path(dir, "summary_stats.rds")
  list(dir = dir, scenarios = sc,
       stats = if (file.exists(stats_path)) readRDS(stats_path) else NULL)
}

# A nested params list as a flat named list of leaves ("a$b$c").
flatten_params <- function(x, prefix = "") {
  if (is.list(x) && !is.data.frame(x) && length(x) > 0) {
    nms <- names(x); if (is.null(nms)) nms <- rep("", length(x))
    out <- list()
    for (k in seq_along(x)) {
      nm <- if (nzchar(nms[k])) nms[k] else as.character(k)
      out <- c(out, flatten_params(x[[k]], if (nzchar(prefix)) paste0(prefix, "$", nm) else nm))
    }
    return(out)
  }
  stats::setNames(list(x), prefix)
}

format_value <- function(x, width = 50) {
  s <- if (is.null(x)) "(absent)"
       else if (is.function(x)) paste(deparse(x), collapse = " ")
       else if (length(x) == 0) "(empty)"
       else paste(format(x, digits = 6), collapse = ", ")
  if (nchar(s) > width) paste0(substr(s, 1, width - 3), "...") else s
}

# Params that differ between the two sets: one row per (parameter, old, new),
# with the scenarios where that difference occurs.
compare_params <- function(old, new) {
  rows <- list()
  for (nm in intersect(names(old$scenarios), names(new$scenarios))) {
    a <- flatten_params(old$scenarios[[nm]]$params)
    b <- flatten_params(new$scenarios[[nm]]$params)
    for (k in union(names(a), names(b))) {
      if (!identical(a[[k]], b[[k]]))
        rows[[length(rows) + 1]] <- data.frame(parameter = k, old = format_value(a[[k]]),
                                               new = format_value(b[[k]]), scenario = nm)
    }
  }
  if (length(rows) == 0)
    return(data.frame(parameter = character(), old = character(), new = character(), scenarios = character()))
  d <- do.call(rbind, rows)
  key <- paste(d$parameter, d$old, d$new, sep = "\r")
  agg <- d[!duplicated(key), c("parameter", "old", "new")]
  agg$scenarios <- vapply(key[!duplicated(key)], function(k) paste(d$scenario[key == k], collapse = ", "), "")
  rownames(agg) <- NULL
  agg
}

# Headline numbers of one scenario output (named list; numbers or strings).
scenario_metrics <- function(out) {
  m <- list()
  cen <- out$census
  if (!is.null(cen)) {
    ext <- cen$year[as.logical(cen$extinct)]
    m <- c(m, list(
      N_alive_max   = max(cen$N_alive),
      N_alive_final = utils::tail(cen$N_alive, 1),
      N_IUCN_final  = utils::tail(cen$N_IUCN, 1),
      births        = sum(cen$births),
      deaths        = sum(cen$deaths),
      deaths_fire   = sum(cen$deaths_fire),
      fires_reached = sum(as.logical(cen$fire_event)),
      resist_final  = round(utils::tail(cen$mean_resist_score, 1), 4),
      extinct_year  = if (length(ext)) min(ext) else NA))
  }
  pegg <- c("N_prefire", "fraction_dead", "fraction_fire", "fraction_postfire", "fraction_resprout_dead")
  a4   <- c("T_g", "window_len", "N_ref", "N_compare", "pct_decline", "iucn_cat",
            "pct_decline_flowering", "iucn_cat_flowering")
  for (k in c(pegg, a4)) if (!is.null(out[[k]])) m[[k]] <- out[[k]]
  if (!is.null(out$census_cfact))
    m$cfact_N_alive_final <- utils::tail(out$census_cfact$N_alive, 1)
  m
}

# Old vs new headline numbers for every scenario in both sets.
compare_outputs <- function(old, new) {
  rows <- list()
  for (nm in intersect(names(old$scenarios), names(new$scenarios))) {
    a <- scenario_metrics(old$scenarios[[nm]]$output)
    b <- scenario_metrics(new$scenarios[[nm]]$output)
    for (k in union(names(a), names(b))) {
      x <- a[[k]]; y <- b[[k]]
      num <- is.numeric(x) && is.numeric(y) && !is.na(x) && !is.na(y)
      ch  <- if (num && x != 0) sprintf("%+.1f%%", 100 * (y - x) / abs(x))
             else if (identical(x, y)) "" else "changed"
      if (num && x == y) ch <- ""
      rows[[length(rows) + 1]] <- data.frame(
        scenario = nm, metric = k,
        old = format_value(x, 20), new = format_value(y, 20), change = ch)
    }
  }
  do.call(rbind, rows)
}

# Old vs new 20-seed statistics, with the difference in units of the test
# tolerance (computed from the old baseline, as in test-regression.R).
compare_stats <- function(old, new) {
  if (is.null(old$stats) || is.null(new$stats)) return(NULL)
  rows <- list()
  for (nm in intersect(names(old$stats$stats), names(new$stats$stats))) {
    a <- old$stats$stats[[nm]]; b <- new$stats$stats[[nm]]
    z   <- if (!is.null(old$stats$meta$z)) old$stats$meta$z else 4
    tol <- z * sqrt(2) * a$sd / sqrt(length(a$values))
    rows[[length(rows) + 1]] <- data.frame(
      statistic = nm, old_mean = a$mean, old_sd = a$sd, new_mean = b$mean, new_sd = b$sd,
      tolerance = tol, diff_in_tol = (b$mean - a$mean) / tol,
      verdict = if (abs(b$mean - a$mean) <= tol) "within" else "OUTSIDE")
  }
  do.call(rbind, rows)
}

# The whole comparison as a list; report_comparison() prints it.
compare_reference_sets <- function(old, new) {
  only_old <- setdiff(names(old$scenarios), names(new$scenarios))
  only_new <- setdiff(names(new$scenarios), names(old$scenarios))
  common   <- intersect(names(old$scenarios), names(new$scenarios))
  out_diff <- vapply(common, function(nm)
    { d <- first_difference(new$scenarios[[nm]]$output, old$scenarios[[nm]]$output); if (is.null(d)) "" else d }, "")
  meta_of <- function(s) if (length(s$scenarios)) s$scenarios[[1]]$meta else list()
  list(old_dir = old$dir, new_dir = new$dir,
       old_meta = meta_of(old), new_meta = meta_of(new),
       old_stats_meta = old$stats$meta, new_stats_meta = new$stats$meta,
       only_old = only_old, only_new = only_new,
       output_differs = out_diff,
       params = compare_params(old, new),
       outputs = compare_outputs(old, new),
       stats = compare_stats(old, new),
       stats_present = c(old = !is.null(old$stats), new = !is.null(new$stats)))
}

# Which row of the table in docs/reference_updates.md the differences point to.
# Only an indication: whether the code changed behaviour is shown by the code
# tests run against the old set, not by comparing the two sets.
classify_comparison <- function(cmp) {
  notes <- character()
  same_R <- identical(cmp$old_meta$R_version, cmp$new_meta$R_version)
  if (!same_R)
    notes <- c(notes, sprintf("R version changed (%s -> %s): row 4, R upgrade.",
                              cmp$old_meta$R_version, cmp$new_meta$R_version))
  if (!identical(cmp$old_meta$scenarios_md5, cmp$new_meta$scenarios_md5))
    notes <- c(notes, "scenarios.R changed: the scenarios themselves differ, so outputs are not comparable like for like.")
  if (length(cmp$only_old) || length(cmp$only_new))
    notes <- c(notes, "The two sets do not contain the same scenarios.")
  n_out <- sum(nzchar(cmp$output_differs))
  stats_out <- !is.null(cmp$stats) && any(cmp$stats$verdict == "OUTSIDE")
  if (nrow(cmp$params) > 0) {
    notes <- c(notes, sprintf("%d parameter difference(s): calibration changed, row 1 (or row 2 if the code changed too).",
                              nrow(cmp$params)))
  } else if (n_out > 0) {
    if (is.null(cmp$stats))
      notes <- c(notes, "Params identical, outputs differ: a code change (row 2 or 3). Generate the candidate stats to tell which.")
    else if (stats_out)
      notes <- c(notes, "Params identical, outputs differ, a statistic moved beyond tolerance: behaviour change, row 2.")
    else
      notes <- c(notes, "Params identical, outputs differ, statistics within tolerance: consistent with a draw-order change, row 3 (keep the stats baseline).")
  } else if (same_R) {
    notes <- c(notes, "Params and outputs identical: nothing to update.")
  }
  notes
}

md_table <- function(d) {
  if (is.null(d) || nrow(d) == 0) return("(none)")
  esc <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")))
}

# The comparison as Markdown lines, for the update report.
report_comparison <- function(cmp, all_metrics = FALSE) {
  m <- function(x) if (length(x) == 0) "-" else paste(format(x), collapse = " ")
  meta_rows <- data.frame(
    field = c("created", "R_version", "git_commit", "git_dirty", "scenarios_md5", "scale"),
    old   = vapply(c("created", "R_version", "git_commit", "git_dirty", "scenarios_md5", "scale"),
                   function(k) m(cmp$old_meta[[k]]), ""),
    new   = vapply(c("created", "R_version", "git_commit", "git_dirty", "scenarios_md5", "scale"),
                   function(k) m(cmp$new_meta[[k]]), ""))
  outs <- cmp$outputs
  if (!all_metrics && !is.null(outs)) outs <- outs[outs$change != "" | outs$metric %in%
                                                     c("N_alive_final", "fraction_resprout_dead", "pct_decline", "iucn_cat", "T_g"), ]
  stats <- cmp$stats
  if (!is.null(stats)) {
    for (k in c("old_mean", "old_sd", "new_mean", "new_sd", "tolerance")) stats[[k]] <- signif(stats[[k]], 4)
    stats$diff_in_tol <- sprintf("%+.2f", stats$diff_in_tol)
  }
  od <- cmp$output_differs
  c("# Reference comparison",
    "",
    sprintf("- old: `%s`", cmp$old_dir),
    sprintf("- new: `%s`", cmp$new_dir),
    "",
    "## Indication (docs/reference_updates.md table)",
    "",
    paste0("- ", classify_comparison(cmp)),
    "",
    "## Provenance",
    "",
    md_table(meta_rows),
    "",
    "## Scenarios",
    "",
    if (length(cmp$only_old)) sprintf("- only in old: %s", paste(cmp$only_old, collapse = ", ")),
    if (length(cmp$only_new)) sprintf("- only in new: %s", paste(cmp$only_new, collapse = ", ")),
    paste0("- **", names(od), "**: ", ifelse(nzchar(od), od, "identical")),
    "",
    "## Parameters that differ",
    "",
    md_table(cmp$params),
    "",
    sprintf("## Outputs (%s)", if (all_metrics) "all metrics" else "changed metrics, plus key results; Pegg observed 0.48"),
    "",
    paste("One seed per scenario. If the random stream changed, these are a new random draw, not evidence of a",
          "behaviour change: a draw-order change alone moved the A4 run from EN to CR (tested 2026-11-19).",
          "Judge behaviour from the 20-seed statistics below."),
    "",
    md_table(outs),
    "",
    "## 20-seed statistics (difference in units of the test tolerance)",
    "",
    if (is.null(stats))
      sprintf("Not compared: summary_stats.rds %s.",
              if (!cmp$stats_present[["new"]]) "absent from the new set (baseline kept?)" else "absent from the old set")
    else md_table(stats))
}
