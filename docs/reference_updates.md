# Procedure: updating the regression references

Status: **accepted 2026-11-19** (PLAN.md Phase 1; methodological notes T3). Decisions are at the end.

## Why a procedure is needed

The references in `tests/reference/` are frozen (CLAUDE.md rule 2; read-only on disk; `Edit` denied
in `.claude/settings.json`). But some changes are meant to alter model output, and after them the
regression tests fail by design. The first ones are already scheduled: the ROADMAP M. nodosa work
(Pegg per-tree refit, `range_v2` site table) will change `params.R`.

The procedure must keep two things true:
1. A failing regression test is never "fixed" by regenerating the references to match the new output.
   Each update is a deliberate, recorded decision by the project owner.
2. Each reference set can be traced to a commit, an R version and a reason.

## What a failure means

The tests separate the kinds of change. The table says which ones justify an update.

| What fails | What changed | Update the references? |
|---|---|---|
| params tests only (code tests pass) | `params.R`, e.g. a recalibration | Yes, if the calibration change was intended and logged |
| code tests and stats tests | Model behaviour | Only if the behaviour change was intended (an approved design note or a logged decision). Otherwise it is a bug: fix the code. |
| code tests only (stats pass) | Draw order or other RNG use; behaviour unchanged within sampling error | Yes for the exact-match files; **keep** the stats baseline (see below) |
| reference-file test (R version) | The R version | Yes, after the stats tests pass under the new R |
| stats tests only | Unusual: behaviour moved, but the fixed seeds still match | Investigate; do not update |

**Keep the stats baseline across refactors.** If `summary_stats.rds` were regenerated after every
refactor, small shifts could add up unseen, each one inside the tolerance. The stats baseline is
regenerated only when behaviour changes on purpose (rows 1, 2 and 4), never for a draw-order refactor
(row 3).

## Who does what

- **Agent:** makes the code or `params.R` change, runs the suite and writes a comparison report (step 1).
  The agent never writes to `tests/reference/`, never runs the generators with their default output
  directory and never changes file permissions there.
- **Project owner:** decides whether the change is intended, writes or approves the log entry,
  generates and installs the new references, and commits.

## Procedure

1. **Make the change and report** (agent).
   - Make the change; run the full suite. Only the expected regression tests may fail.
     Any other failure is fixed first.
   - Generate candidate references into a scratch directory, never `tests/reference/`:
     `Rscript tests/regression/make_reference.R <scratch>/candidate` (and `make_reference_stats.R`
     into the same directory, if the stats baseline is to change).
   - Check that the candidates reproduce themselves:
     `MNSIM_REFERENCE_DIR=<scratch>/candidate Rscript tests/testthat.R` passes.
   - Write a short comparison of old and new, per scenario:
     - final and maximum N_alive, and fires reached;
     - the Pegg fraction of resprouters dead (observed 0.48);
     - the A4 category and generation time;
     - the three statistics' means, old and new, in units of their tolerance.
     `tests/regression/compare_references.R <scratch>/candidate --out <scratch>/report.md` writes this
     report and says which row of the table the differences point to. Single-seed numbers are a new
     random draw whenever the random stream changes. A draw-order change alone moved the A4 run
     from EN to CR, so judge behaviour from the statistics.

2. **Decide and log** (owner).
   - Add a dated entry to `docs/methodological_notes.md`, section T, numbered T3, T4, …. It records:
     - the reason for the change;
     - which row of the table above applies;
     - the expected and observed direction of change, from the comparison;
     - the reference version number (below).
   - Say whether the stats baseline is kept or regenerated, and why.

3. **Generate and install** (owner, in a terminal, not through the agent).
   ```
   cd plant_ibm
   rm -rf /tmp/mnsim_ref && /opt/R/4.3.1/bin/Rscript tests/regression/make_reference.R /tmp/mnsim_ref
   /opt/R/4.3.1/bin/Rscript tests/regression/make_reference_stats.R /tmp/mnsim_ref   # only if regenerating stats
   MNSIM_REFERENCE_DIR=/tmp/mnsim_ref /opt/R/4.3.1/bin/Rscript tests/testthat.R      # must pass
   cd ..
   chmod -R u+w tests/reference
   cp /tmp/mnsim_ref/* tests/reference/        # summary_stats.rds is copied only if regenerated
   chmod -R a-w tests/reference
   cd plant_ibm && /opt/R/4.3.1/bin/Rscript tests/testthat.R                          # must pass
   ```
   Generate fresh into an empty directory and copy the whole set. Never regenerate a single
   scenario: all files in a set share one `meta` (commit, R version, `scenarios.R` md5).

4. **Commit together** (owner, or the agent when asked).
   - Commit the change, the new references and the log entry together, so the suite passes at every
     commit (CLAUDE.md rule 3).
   - Start the commit message with `REFERENCE UPDATE vN:` and cite the notes entry.
   - Tag the last commit with the previous set (`ref-v1`, `ref-v2`, …). Old references stay in git
     history and are easy to recover with `git checkout ref-v1 -- tests/reference`.

## Versioning and traceability

- Number the reference sets. The current set is **v1** (2026-11-19, commit `2cf1bb8`).
- Small changes to the generators make sets self-describing. They are to be made at the first update:
  - add `ref_version` to `meta` and to `manifest.csv`;
  - add md5 sums of `params.R` and `R/*.R` to `meta`. At generation time the change is usually not
    committed yet (`git_dirty = TRUE`), so `git_commit` alone does not identify the code.
- Batch changes where possible. For example, one update after the Pegg refit and `range_v2` are both
  done, not one after each, so fewer sets need tracking.

## Adding scenarios (Phase 4 species)

New species will need their own reference scenarios. Adding files is not editing them, but the
generators refuse to run if any reference file exists. Proposed:
- give each species its own subdirectory (`tests/reference/<species>/`);
- generate a species' set with its own run of the generator into an empty subdirectory;
- keep the M. nodosa set at the top level, unchanged.

## Decisions (2026-11-19)

1. The table is accepted, including keeping the stats baseline across draw-order refactors.
2. The agent may generate candidates into scratch and write the comparison report. The owner alone
   installs references.
3. Version numbers, `ref-vN` git tags and the `REFERENCE UPDATE vN:` commit prefix, as proposed.
4. `compare_references.R` built now (2026-11-19). The generator `meta` additions (`ref_version`, code
   and `params.R` md5 sums) wait until the first update.
5. Per-species subdirectories for Phase 4 references.
