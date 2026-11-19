# tests/regression/make_reference_stats.R
#
# Freezes the 20-seed summary statistics (tests/regression/stats.R) for the
# second regression check.
#
#   cd plant_ibm
#   /opt/R/4.3.1/bin/Rscript tests/regression/make_reference_stats.R [out_dir]
#
# out_dir defaults to ../tests/reference (repo root). Never overwrites: stops
# if summary_stats.rds already exists. Regenerate only deliberately, by the
# project owner (PLAN.md Phase 1; CLAUDE.md rule 2).

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else file.path("..", "tests", "reference")

for (f in c("individuals", "genetics", "mortality", "fire", "recruitment", "census", "simulate"))
  source(file.path("R", paste0(f, ".R")))
source("runner/params_utils.R")
source("runner/iucn_a4.R")
source("runner/pegg_validation.R")
source("tests/regression/scenarios.R")
source("tests/regression/stats.R")

path <- file.path(out_dir, REGRESSION_STATS_FILE)
if (file.exists(path))
  stop(path, " already exists. Refusing to overwrite. Regenerate only deliberately, into an empty directory.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

git <- function(...) tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE),
                              error = function(e) NA_character_)
meta <- list(
  created    = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  R_version  = R.version.string,
  RNGkind    = RNGkind(),
  git_commit = git("rev-parse", "HEAD"),
  git_dirty  = length(git("status", "--porcelain", "--", ".")) > 0,
  stats_md5  = unname(tools::md5sum("tests/regression/stats.R")),
  n          = REGRESSION_STATS_N,
  z          = REGRESSION_STATS_Z
)

base <- regression_base()
out  <- list()
for (nm in names(REGRESSION_STATS)) {
  ps <- regression_stats_params(nm, base)
  tm <- system.time(v <- regression_stats_values(nm, ps))[["elapsed"]]
  out[[nm]] <- list(description = REGRESSION_STATS[[nm]]$description,
                    params = ps, values = v, mean = mean(v), sd = stats::sd(v))
  message(sprintf("%-22s mean %8.4f  sd %8.4f  tolerance +/- %.4f  (%.1f s)", nm, mean(v), sd(v),
                  REGRESSION_STATS_Z * sqrt(2) * sd(v) / sqrt(length(v)), tm))
}
saveRDS(list(stats = out, meta = meta), path)
message("Wrote ", path, "  md5 ", unname(tools::md5sum(path)))
