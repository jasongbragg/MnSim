# tests/regression/compare_references.R
#
# Compares a candidate set of regression references with the current set and
# prints a Markdown report: provenance, which scenarios differ and where,
# parameters that differ, headline outputs (population, Pegg fraction, A4
# category, ...) and the 20-seed statistics in units of the test tolerance.
# Step 1 of the reference-update procedure (docs/reference_updates.md).
#
#   cd plant_ibm
#   /opt/R/4.3.1/bin/Rscript tests/regression/compare_references.R NEW_DIR
#   /opt/R/4.3.1/bin/Rscript tests/regression/compare_references.R OLD_DIR NEW_DIR [--all] [--out FILE]
#
# OLD_DIR defaults to ../tests/reference. --all lists every metric, not only
# those that changed. --out also writes the report to FILE. Reads only; never
# writes to either directory.

args  <- commandArgs(trailingOnly = TRUE)
all_m <- "--all" %in% args
out_f <- NULL
if ("--out" %in% args) {
  i <- match("--out", args)
  if (i == length(args)) stop("--out needs a file name")
  out_f <- args[i + 1]
  args  <- args[-c(i, i + 1)]
}
dirs <- setdiff(args, "--all")
if (length(dirs) == 1) dirs <- c(file.path("..", "tests", "reference"), dirs)
if (length(dirs) != 2)
  stop("usage: compare_references.R [OLD_DIR] NEW_DIR [--all] [--out FILE]")

source("tests/regression/compare.R")

cmp <- compare_reference_sets(load_reference_set(dirs[1]), load_reference_set(dirs[2]))
rep <- report_comparison(cmp, all_metrics = all_m)
writeLines(rep)
if (!is.null(out_f)) {
  if (normalizePath(dirname(out_f)) %in% normalizePath(dirs))
    stop("--out must not write into a reference directory")
  writeLines(rep, out_f)
  message("Wrote ", out_f)
}
