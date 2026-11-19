# tests/regression/make_reference.R
#
# Generates the frozen M. nodosa reference outputs for the regression tests.
#
#   cd plant_ibm
#   /opt/R/4.3.1/bin/Rscript tests/regression/make_reference.R [out_dir]
#
# out_dir defaults to ../tests/reference (repo root). The script NEVER
# overwrites: it stops if any reference file already exists. Once created,
# the references are read-only (chmod -R a-w tests/reference) and are only
# regenerated deliberately by the project owner (PLAN.md Phase 1; CLAUDE.md
# rule 2), into an empty directory.

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else file.path("..", "tests", "reference")

for (f in c("individuals", "genetics", "mortality", "fire", "recruitment", "census", "simulate"))
  source(file.path("R", paste0(f, ".R")))
source("runner/params_utils.R")
source("runner/iucn_a4.R")
source("runner/pegg_validation.R")
source("tests/regression/scenarios.R")

names_sc <- names(REGRESSION_SCENARIOS)
paths    <- vapply(names_sc, regression_reference_path, "", ref_dir = out_dir)
if (any(file.exists(paths)))
  stop("Reference files already exist in ", out_dir, ": ",
       paste(basename(paths[file.exists(paths)]), collapse = ", "),
       ". Refusing to overwrite. Regenerate only deliberately, into an empty directory.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

git <- function(...) tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE),
                              error = function(e) NA_character_)
meta <- list(
  created      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  R_version    = R.version.string,
  RNGkind      = RNGkind(),
  git_commit   = git("rev-parse", "HEAD"),
  git_dirty    = length(git("status", "--porcelain", "--", ".")) > 0,
  scenarios_md5 = unname(tools::md5sum("tests/regression/scenarios.R")),
  scale        = REGRESSION_SCALE
)

base <- regression_base()
rows <- list()
for (nm in names_sc) {
  p  <- regression_params(nm, base)
  tm <- system.time(out <- regression_run(nm, p))[["elapsed"]]
  saveRDS(list(name = nm, description = REGRESSION_SCENARIOS[[nm]]$description,
               params = p, output = out, meta = meta),
          paths[[nm]])
  cen <- if (!is.null(out$census)) out$census else NULL
  rows[[nm]] <- data.frame(
    scenario  = nm,
    seconds   = round(tm, 1),
    N_alive_max  = if (!is.null(cen)) max(cen$N_alive) else NA,
    N_alive_last = if (!is.null(cen)) tail(cen$N_alive, 1) else NA,
    fires_reached = if (!is.null(cen)) sum(cen$fire_event) else NA,
    md5       = unname(tools::md5sum(paths[[nm]])))
  message(sprintf("%-10s %5.1f s", nm, tm))
}
manifest <- do.call(rbind, rows)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
print(manifest)
str(meta)
