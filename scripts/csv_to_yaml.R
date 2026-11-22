#!/usr/bin/env Rscript
# scripts/csv_to_yaml.R
#
# Builds species/<sp>.yaml from the species' rows in data/ (PLAN.md Phase 4),
# and lists what is still missing (docs/parameter_gaps.md).
#
# THE CONTRACT (CLAUDE.md rule 5; notes T22)
#   The CSV is the single source. A model field is filled by a row whose
#   `parameter` is the field's path joined with '__' (as in
#   data/m_nodosa_params.csv), role `parameter`:
#     - a cited value: the row's value; the YAML entry cites the file;
#     - a prior: notes start "PRIOR [min, max]: reason" (or "PRIOR {a, b}: reason");
#     - derived from the rust anchor: evidence_type `derived` (annual_env_alpha/beta).
#   Literature rows under other names (age_first_reproduction, longevity, ...)
#   are evidence, not model fields: data/param_evidence_map.csv links them to
#   the fields they inform, and the gaps report lists them. Turning evidence
#   into a field row is a judgement (a prior range and its reason), recorded in
#   the row, never done by this script.
#   Settings outside params use rows named rust_anchor__mode,
#   assessment__<field> and spatial__range_box__<field>.
#   Which optional fields are needed follows the species' mode rows
#   (clonal_mode, rust_hazard_mode, recruitment_curve, ...): a field the schema
#   marks required_when is needed only when its condition holds; an optional
#   field with a default is written only if it has a row.
#   Range-wide draws (ensemble_draws) are not generated: they are set by hand
#   when the species' range runs are planned.
#
# USAGE (from the repository root, with R 4.3.1)
#   Rscript scripts/csv_to_yaml.R r_rubescens            # write species/r_rubescens.yaml if complete
#   Rscript scripts/csv_to_yaml.R r_rubescens --out f    # write elsewhere
#   Rscript scripts/csv_to_yaml.R --gaps                 # refresh the generated part of docs/parameter_gaps.md
# A species with gaps is not written; its gaps are printed.

SPECIES <- c(
  m_nodosa      = "Melaleuca nodosa",
  r_rubescens   = "Rhodamnia rubescens",
  r_psidioides  = "Rhodomyrtus psidioides",
  a_beckleri    = "Archirhodomyrtus beckleri",
  d_humile      = "Decaspermum humile",
  b_leptopetala = "Backhousia leptopetala",
  u_australis   = "Uromyrtus australis"
)
csv_file_for <- function(slug) if (slug == "m_nodosa") "m_nodosa_params.csv" else "rust_species_params.csv"

script_root <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(file.path(dirname(sub("^--file=", "", a[1])), ".."))
  else normalizePath(".")
}

# ---- Reading -----------------------------------------------------------------

read_species_csv <- function(slug, root) {
  d <- utils::read.csv(file.path(root, "data", csv_file_for(slug)), colClasses = "character")
  d[d$species == SPECIES[[slug]], ]
}

# Field rows (role parameter, named by a schema path), keyed by parameter.
field_rows <- function(d) {
  d <- d[d$role == "parameter", ]
  if (anyDuplicated(d$parameter)) stop("duplicate parameter rows: ",
                                       paste(unique(d$parameter[duplicated(d$parameter)]), collapse = ", "))
  split(d, d$parameter) |> lapply(as.list)
}

# The species' schema fields and their state: one row per leaf of layer
# species in params (plus rust_anchor$mode, assessment, spatial).
# status: row (filled) | gap (needed, no row) | default (optional, no row) |
#         inactive (required_when does not hold).
field_status <- function(schema, rows) {
  pv <- function(key, spec_of) {   # current value of a mode field: row, else default
    r <- rows[[key]]
    if (!is.null(r)) return(r$value)
    if (!is.null(spec_of$default)) as.character(spec_of$default) else NA_character_
  }
  out <- list()
  add <- function(block, path, spec) {
    key <- paste(c(if (block != "params") block, path), collapse = "__")
    st <- if (!is.null(rows[[key]])) "row" else if (isTRUE(spec$optional)) "default" else "gap"
    rw <- spec$required_when
    cond <- NA_character_
    if (!is.null(rw)) {
      f <- names(rw)[1]
      cond <- sprintf("%s = %s", f, rw[[f]])
      if (!identical(pv(f, schema$params[[f]]), as.character(rw[[f]])))
        st <- if (st == "row") "row_unused" else "inactive"
      else if (st == "default") st <- "gap"
    }
    out[[length(out) + 1L]] <<- data.frame(block = block, key = key, path = paste(path, collapse = "$"),
                                           type = spec$type, status = st, condition = cond)
  }
  for (lf in schema_leaves(schema$rust_anchor))
    if (identical(lf$spec$layer, "species")) add("rust_anchor", lf$path, lf$spec)
  for (lf in schema_leaves(schema$params))
    if (identical(lf$spec$layer, "species")) add("params", lf$path, lf$spec)
  for (b in c("assessment", "spatial"))
    for (lf in schema_leaves(schema[[b]]))
      if (identical(lf$spec$layer, "species")) add(b, lf$path, lf$spec)
  do.call(rbind, out)
}

# ---- Writing -----------------------------------------------------------------

typed_value <- function(raw, type) {
  if (type %in% c("double", "integer")) {
    n <- suppressWarnings(as.numeric(raw))
    return(if (is.na(n)) raw else n)            # keeps forms like "1/3"
  }
  if (type == "double_vector") return(as.numeric(trimws(strsplit(raw, ";")[[1]])))
  raw
}

param_entry <- function(row, spec, file) {
  if (identical(row$evidence_type, "derived") && !is.null(spec$derive))
    return(list(derived = "rust_anchor", cite = file))
  e  <- list(value = typed_value(row$value, spec$type))
  pr <- csv_prior_range(row$notes)
  if (!is.null(pr)) {
    e$prior <- TRUE
    if (!is.null(pr$choices)) e$choices <- as.list(pr$choices)
    else {                                       # raw text, so forms like 1/3 stay exact
      items <- trimws(strsplit(sub("^PRIOR \\[([^]]*)\\].*$", "\\1", row$notes), ",")[[1]])
      e$min <- typed_value(items[1], "double"); e$max <- typed_value(items[2], "double")
    }
    e$reason <- sub("^PRIOR [\\[{][^]}]*[]}]:?\\s*", "", row$notes)
  }
  e$cite <- file
  e
}

set_path <- function(x, path, v) {
  if (length(path) == 1L) { x[[path]] <- v; return(x) }
  x[[path[1]]] <- set_path(if (is.null(x[[path[1]]])) list() else x[[path[1]]], path[-1], v)
  x
}

# blocks: the blocks to build (tests rebuild params alone).
build_species_yaml <- function(slug, root, blocks = c("rust_anchor", "params", "assessment", "spatial")) {
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  rows   <- field_rows(read_species_csv(slug, root))
  st     <- field_status(schema, rows)
  st     <- st[st$block %in% blocks, ]
  gaps   <- st[st$status == "gap", ]
  unused <- st[st$status == "row_unused", ]
  if (nrow(unused)) warning("rows not used (their mode is not selected): ", paste(unused$key, collapse = ", "))
  if (nrow(gaps)) return(list(ok = FALSE, gaps = gaps, status = st))

  file <- csv_file_for(slug)
  specs <- list(rust_anchor = schema$rust_anchor, params = schema$params,
                assessment = schema$assessment, spatial = schema$spatial)
  doc <- list(species = SPECIES[[slug]])
  for (i in which(st$status == "row")) {
    b <- st$block[i]; path <- strsplit(st$path[i], "$", fixed = TRUE)[[1]]
    spec <- cfg_get(specs[[b]], path)
    row  <- rows[[st$key[i]]]
    e <- if (b %in% c("assessment", "spatial")) list(cite = file, row = st$key[i])
         else param_entry(row, spec, file)
    doc[[b]] <- set_path(if (is.null(doc[[b]])) list() else doc[[b]], path, e)
  }
  doc <- doc[intersect(c("species", "rust_anchor", "params", "assessment", "spatial"), names(doc))]
  header <- c(
    sprintf("# species/%s.yaml", slug), "#",
    sprintf("# %s: generated on %s by scripts/csv_to_yaml.R from data/%s.", SPECIES[[slug]], Sys.Date(), file),
    "# Do not edit by hand: change the CSV rows and regenerate (notes T22).",
    "# Each value cites its CSV row (species + parameter = the path joined with '__')",
    "# or is a prior whose range and reason come from that row (CLAUDE.md rule 5).", "")
  list(ok = TRUE, text = c(header, yaml::as.yaml(doc, indent.mapping.sequence = TRUE)), status = st)
}

# ---- Gaps report ---------------------------------------------------------------

# Field group for the report: top-level param name (dose responses as one line).
field_group <- function(st) ifelse(st$block == "params", sub("\\$.*", "", st$path),
                                   paste0(st$block, "$", sub("\\$.*", "", st$path)))

evidence_for <- function(group, slug, root, emap, d) {
  par <- emap$csv_parameter[emap$field == group]
  if (!length(par)) return(d[0, ])
  base <- sub("\\[.*$", "", d$parameter)
  d[base %in% par, ]
}

gaps_markdown <- function(root) {
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  emap   <- utils::read.csv(file.path(root, "data", "param_evidence_map.csv"), colClasses = "character",
                            comment.char = "#")
  slugs  <- setdiff(names(SPECIES), "m_nodosa")
  per <- lapply(slugs, function(s) {
    d <- read_species_csv(s, root)
    st <- field_status(schema, field_rows(d))
    st$group <- field_group(st)
    list(d = d, st = st)
  })
  names(per) <- slugs
  groups <- unique(per[[1]]$st$group)
  groups <- groups[groups != "assessment$criterion"]
  ab <- c(r_rubescens = "R. rub", r_psidioides = "R. psi", a_beckleri = "A. bec", d_humile = "D. hum",
          b_leptopetala = "B. lep", u_australis = "U. aus")

  cell <- function(s, g) {
    st <- per[[s]]$st[per[[s]]$st$group == g, ]
    ev <- evidence_for(g, s, root, emap, per[[s]]$d)
    state <- if (all(st$status %in% c("row", "row_unused"))) "filled"
             else if (all(st$status %in% c("inactive"))) "n/a"
             else if (all(st$status %in% c("default", "inactive"))) "default"
             else "**gap**"
    ev <- ev[!ev$value %in% c("NA", ""), ]       # NA rows are listed below, not counted
    np <- sum(ev$role == "parameter"); nt <- sum(ev$role == "calibration_target")
    paste0(state, if (np + nt) sprintf(" (%s)", paste(c(if (np) paste0(np, "e"), if (nt) paste0(nt, "t")), collapse = " ")))
  }
  cond <- function(g) {
    c0 <- unique(stats::na.omit(per[[1]]$st$condition[per[[1]]$st$group == g]))
    if (length(c0)) paste0(" (if ", c0[1], ")") else ""
  }
  out <- c("<!-- BEGIN GENERATED by scripts/csv_to_yaml.R --gaps; do not edit between the markers -->",
           sprintf("_Generated %s._", Sys.Date()), "",
           "### Fields by species", "",
           "**filled**: has its field row; **gap**: needed and has none; **default**: optional, schema default applies;",
           "**n/a**: only needed under a mode the species does not (yet) use. In brackets: literature rows that",
           "bear on it with a value (e = parameter rows, t = calibration targets), listed per species below.", "",
           paste0("| Field | ", paste(ab[slugs], collapse = " | "), " |"),
           paste0("|---|", paste(rep("---", length(slugs)), collapse = "|"), "|"))
  for (g in groups)
    out <- c(out, paste0("| `", g, "`", cond(g), " | ", paste(vapply(slugs, cell, "", g = g), collapse = " | "), " |"))
  for (s in slugs) {
    out <- c(out, "", sprintf("### %s: evidence rows by field", SPECIES[[s]]), "")
    any_ev <- FALSE
    for (g in groups) {
      ev <- evidence_for(g, s, root, emap, per[[s]]$d)
      if (!nrow(ev)) next
      any_ev <- TRUE
      items <- sprintf("`%s` = %s (%s, %s)", ev$parameter, substr(gsub("\\|", "/", ev$value), 1, 60),
                       ev$role, ev$confidence)
      out <- c(out, sprintf("- `%s`: %s", g, paste(items, collapse = "; ")))
    }
    if (!any_ev) out <- c(out, "- none")
    n_gap <- sum(per[[s]]$st$status == "gap")
    out <- c(out, "", sprintf("Needed fields with no row: %d of %d.", n_gap,
                              sum(per[[s]]$st$status != "inactive" & per[[s]]$st$status != "default")))
  }
  c(out, "<!-- END GENERATED -->")
}

write_gaps <- function(root) {
  f <- file.path(root, "docs", "parameter_gaps.md")
  txt <- readLines(f)
  a <- grep("^<!-- BEGIN GENERATED", txt); b <- grep("^<!-- END GENERATED", txt)
  if (length(a) != 1L || length(b) != 1L || b < a) stop("docs/parameter_gaps.md needs one pair of GENERATED markers")
  writeLines(c(txt[seq_len(a - 1L)], gaps_markdown(root), txt[-seq_len(b)]), f)
  message("updated ", f)
}

# ---- Main ----------------------------------------------------------------------

main <- function(args) {
  root <- script_root()
  source(file.path(root, "plant_ibm", "R", "species_config.R"))
  assign("SPECIES_ROOT", root, envir = globalenv())
  if ("--gaps" %in% args) return(invisible(write_gaps(root)))
  slug <- args[!startsWith(args, "--")][1]
  if (is.na(slug) || !slug %in% names(SPECIES))
    stop("usage: csv_to_yaml.R <species> [--out file] | --gaps; species: ", paste(names(SPECIES), collapse = ", "))
  out <- if ("--out" %in% args) args[match("--out", args) + 1L] else file.path(root, "species", paste0(slug, ".yaml"))
  r <- build_species_yaml(slug, root)
  if (!r$ok) {
    message(sprintf("%s: %d needed fields have no row; nothing written:", slug, nrow(r$gaps)))
    message(paste0("  ", r$gaps$key, collapse = "\n"))
    quit(status = 1)
  }
  writeLines(r$text, out)
  message("wrote ", out)
}

if (sys.nframe() == 0L) main(commandArgs(TRUE))
