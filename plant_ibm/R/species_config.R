# R/species_config.R
#
# Species configuration: assembles the model's params list from YAML files
# and validates it (PLAN.md Phase 2; docs/methodological_notes.md T4, T7).
#
#   config/schema.yaml        every field: layer, type, bounds, order
#   species/<sp>.yaml         species biology, rust-pressure mode (layer "species")
#   config/landscape.yaml     shared landscape defaults and the rust anchor (layer "landscape")
#   config/run_defaults.yaml  run settings (layer "run"), and the default species
#   config/assessment.yaml    IUCN assessment settings (layer "assessment")
#   species/<sp>_validation.yaml  validation set-up (layer V; validation_config())
#
# The assembled list follows the schema's order, so for M. nodosa it is
# identical() to the old get_default_params() in params.R (and gives the same
# make_run_id() hashes).
#
# Numbers: every YAML number is parsed by R's as.numeric() (yaml handlers), then
# coerced to the schema type, so 10000 is a double where the schema says double.
# A value may also be written as a fraction ("1/3"). '.inf' is Inf.
#
# Paths are relative to the repo root: MNSIM_ROOT, or the nearest directory at
# or above the working directory that holds config/schema.yaml (scripts run
# from plant_ibm/, so usually its parent).
#
# Usage (species = NULL: the default species in config/run_defaults.yaml):
#   params <- load_species()
#   anchor <- load_rust_anchor()       # was RUST_ANCHOR
#   priors <- species_priors()         # data frame, one row per prior
#   assess <- load_assessment()        # iucn_assess_spec() defaults (notes T9)
#   cats   <- load_categories()        # IUCN thresholds c(VU, EN, CR)
#   chk    <- load_window_checks()     # window lengths, in T_g, below which a run warns
#   draws  <- load_param_draws()       # range-wide parameter draws (was DEFAULT_PARAM_DRAWS)
#   sp     <- load_species_spatial()   # name, range box, data and figure folders
#   land   <- load_landscape_spatial() # grid size, border latitude, rust SDM path
#   val    <- validation_config()      # e.g. val$pegg (layer V)

mnsim_root <- function() {
  env <- Sys.getenv("MNSIM_ROOT")
  if (nzchar(env)) return(env)
  d <- normalizePath(".", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "config", "schema.yaml"))) return(d)
    up <- dirname(d)
    if (identical(up, d)) break
    d <- up
  }
  normalizePath("..", mustWork = FALSE)   # not found: the old default (fails with a clear message)
}

# The species used when none is named (config/run_defaults.yaml, 'species').
default_species <- function(root = mnsim_root()) {
  s <- read_config_yaml(file.path(root, "config", "run_defaults.yaml"))$species
  if (!is.character(s) || length(s) != 1L || !nzchar(s))
    stop("config/run_defaults.yaml needs 'species: <name of a file in species/>'")
  s
}

read_config_yaml <- function(path) {
  if (!file.exists(path)) stop("config file not found: ", path)
  num <- function(x) {
    if (grepl("^[-+]?\\.(inf|Inf|INF)$", x)) return(if (startsWith(x, "-")) -Inf else Inf)
    if (grepl("^\\.(nan|NaN|NAN)$", x)) return(NaN)
    as.numeric(x)
  }
  yaml::read_yaml(path, handlers = list(float = num, int = num))
}

# A number written in a config or CSV: numeric, or a string "a/b" or "1.5".
config_number <- function(x) {
  if (is.numeric(x)) return(as.double(x))
  if (is.character(x) && length(x) == 1L) {
    m <- regmatches(x, regexec("^\\s*([-+0-9.eE]+)\\s*/\\s*([-+0-9.eE]+)\\s*$", x))[[1]]
    if (length(m) == 3L) return(as.numeric(m[2]) / as.numeric(m[3]))
    v <- suppressWarnings(as.numeric(x))
    if (!is.na(v)) return(v)
  }
  NA_real_
}

# Schema leaves (mappings with a 'type'), in order, with their paths.
schema_leaves <- function(s, path = character()) {
  if (is.list(s) && !is.null(s$type)) return(list(list(path = path, spec = s)))
  out <- list()
  for (k in names(s)) out <- c(out, schema_leaves(s[[k]], c(path, k)))
  out
}

cfg_get <- function(x, path) {
  for (p in path) { if (!is.list(x) || is.null(x[[p]])) return(NULL); x <- x[[p]] }
  x
}

cfg_set <- function(x, path, value) {
  if (length(path) == 1L) { x[path] <- list(value); return(x) }
  sub <- if (is.null(x[[path[1]]])) list() else x[[path[1]]]
  x[[path[1]]] <- cfg_set(sub, path[-1], value)
  x
}

# Paths of the leaves present in a config file, to find fields the schema lacks.
config_paths <- function(x, is_leaf, path = character()) {
  if (is_leaf(x) || !is.list(x) || is.null(names(x))) return(list(path))
  out <- list()
  for (k in names(x)) out <- c(out, config_paths(x[[k]], is_leaf, c(path, k)))
  out
}
is_entry <- function(x) is.list(x) && !is.null(names(x)) && any(c("value", "derived") %in% names(x))

# Coerce a raw value to a schema type; returns list(value, error).
coerce_config_value <- function(raw, spec) {
  bad <- function(msg) list(value = NULL, error = msg)
  switch(spec$type,
    double = {
      v <- config_number(raw)
      if (length(v) != 1L || is.na(v)) return(bad(sprintf("expected a number, got '%s'", format(raw))))
      list(value = v, error = NULL)
    },
    integer = {
      v <- config_number(raw)
      if (length(v) != 1L || is.na(v) || v != round(v)) return(bad("expected a whole number"))
      list(value = as.integer(v), error = NULL)
    },
    character = {
      if (!is.character(raw) || length(raw) != 1L) return(bad("expected a single string"))
      list(value = raw, error = NULL)
    },
    double_vector = {
      v <- vapply(as.list(unlist(raw)), config_number, numeric(1))
      if (length(v) < 1L || anyNA(v)) return(bad("expected one or more numbers"))
      list(value = as.double(v), error = NULL)
    },
    integer_vector = {
      v <- if (length(raw) == 0L) numeric(0) else vapply(as.list(unlist(raw)), config_number, numeric(1))
      if (anyNA(v) || any(v != round(v))) return(bad("expected whole numbers"))
      list(value = as.integer(v), error = NULL)
    },
    bad(sprintf("unknown schema type '%s'", spec$type)))
}

# Bounds and choices from the schema; returns an error string or NULL.
check_config_bounds <- function(v, spec) {
  if (!is.null(spec$choices) && !all(v %in% unlist(spec$choices)))
    return(sprintf("'%s' is not one of: %s", paste(v, collapse = ", "), paste(unlist(spec$choices), collapse = ", ")))
  if (is.numeric(v) && length(v)) {
    if (!is.null(spec$min) && any(v < config_number(spec$min)))
      return(sprintf("%s is below the minimum %s", paste(v, collapse = ", "), spec$min))
    if (!is.null(spec$max) && any(v > config_number(spec$max)))
      return(sprintf("%s is above the maximum %s", paste(v, collapse = ", "), spec$max))
  }
  NULL
}

# Does a CSV value equal a config value? Numbers to 1e-9 (relative), strings exactly.
config_values_match <- function(csv_value, v) {
  if (is.character(v)) return(identical(csv_value, v))
  n <- config_number(csv_value)
  !is.na(n) && length(v) == 1L && abs(n - v) <= 1e-9 * max(1, abs(v))
}

# Range or choices in a CSV prior note ("PRIOR [a, b]: ..." / "PRIOR {a, b}: ..."), or NULL.
csv_prior_range <- function(note) {
  m <- regmatches(note, regexec("^PRIOR ([\\[{])([^]}]*)[]}]", note))[[1]]
  if (length(m) != 3L) return(NULL)
  items <- trimws(strsplit(m[3], ",")[[1]])
  if (m[2] == "[") list(min = config_number(items[1]), max = config_number(items[2]))
  else list(choices = items)
}

# "RANGE [a, b]" anywhere in a CSV note (range of a range-wide draw), or NULL.
csv_draw_range <- function(note) {
  m <- regmatches(note, regexec("RANGE \\[([^]]*)\\]", note))[[1]]
  if (length(m) != 2L) return(NULL)
  items <- trimws(strsplit(m[2], ",")[[1]])
  c(config_number(items[1]), config_number(items[2]))
}

# Cached reader for the CSV files in data/.
csv_reader <- function(root) {
  csvs <- list()
  function(f) {
    if (is.null(csvs[[f]])) {
      p <- file.path(root, "data", f)
      csvs[[f]] <<- if (file.exists(p)) utils::read.csv(p, colClasses = "character") else NA
    }
    csvs[[f]]
  }
}

# Settings outside the params list (config/schema.yaml, notes T9). 'files' maps
# each layer to list(data = <the block in that file>, label = "<file>: <block>").
# An entry is {value, source} or {cite, row} (value from that CSV row of
# 'species_name'). Returns the values in schema order; problems go to err().
is_setting <- function(x) is.list(x) && !is.null(names(x)) && any(c("value", "cite") %in% names(x))
read_settings <- function(schema_block, files, err, csv_of = NULL, species_name = NULL) {
  out    <- list()
  leaves <- schema_leaves(schema_block)
  for (lf in leaves) {
    f     <- files[[lf$spec$layer]]
    where <- paste0(f$label, "$", paste(lf$path, collapse = "$"))
    e     <- cfg_get(f$data, lf$path)
    if (is.null(e)) { err(where, "missing"); next }
    if (!is_setting(e)) { err(where, "expected a mapping with 'value' and 'source' (or 'cite' and 'row')"); next }
    if (!is.null(e$cite)) {
      tab <- if (is.null(csv_of)) NA else csv_of(e$cite)
      if (is.null(e$row)) { err(where, "'cite' needs 'row: <parameter>'"); next }
      if (!is.data.frame(tab)) { err(where, sprintf("cited file data/%s not found", e$cite)); next }
      row <- tab[tab$species %in% species_name & tab$parameter == e$row, ]
      if (nrow(row) != 1L) {
        err(where, sprintf("%d rows in data/%s for %s / %s (need exactly 1)", nrow(row), e$cite,
                           format(species_name), e$row)); next
      }
      raw <- row$value
    } else {
      if (is.null(e$source) || !nzchar(e$source)) err(where, "needs a 'source'")
      raw <- e$value
    }
    cv <- coerce_config_value(raw, lf$spec)
    if (!is.null(cv$error)) { err(where, cv$error); next }
    b <- check_config_bounds(cv$value, lf$spec); if (!is.null(b)) err(where, b)
    out <- cfg_set(out, lf$path, cv$value)
  }
  for (layer in names(files)) {
    known <- vapply(Filter(function(lf) identical(lf$spec$layer, layer), leaves),
                    function(lf) paste(lf$path, collapse = "$"), "")
    for (p in config_paths(files[[layer]]$data, is_setting))
      if (length(p) && !paste(p, collapse = "$") %in% known)
        err(paste0(files[[layer]]$label, "$", paste(p, collapse = "$")),
            "not a setting of this file in config/schema.yaml")
  }
  out
}

config_error <- function(what, errs)
  stop(sprintf("Invalid configuration for %s (%d problem%s):\n  - %s", what,
               length(errs), if (length(errs) > 1L) "s" else "", paste(errs, collapse = "\n  - ")),
       call. = FALSE)

# Read one species' configuration; check_sources also checks citations and
# priors against the CSV files (CLAUDE.md rule 5).
species_config <- function(species = NULL, root = mnsim_root(), check_sources = TRUE) {
  if (is.null(species)) species <- default_species(root)
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  sp     <- read_config_yaml(file.path(root, "species", paste0(species, ".yaml")))
  land   <- read_config_yaml(file.path(root, "config", "landscape.yaml"))
  files  <- list(species   = sp$params,
                 landscape = land$params,
                 run       = read_config_yaml(file.path(root, "config", "run_defaults.yaml"))$params)
  errs   <- character()
  err    <- function(where, msg) errs <<- c(errs, sprintf("%s: %s", where, msg))
  csv_of <- csv_reader(root)
  priors <- list()

  # One entry (species or landscape layer): value, citation or prior. Returns the value.
  read_entry <- function(entry, spec, where, key, layer, anchor = NULL) {
    if (!is_entry(entry)) { err(where, "expected a mapping with 'value' (or 'derived')"); return(NULL) }
    if (!is.null(entry$derived)) {
      if (!identical(entry$derived, "rust_anchor") || is.null(spec$derive)) {
        err(where, "'derived' is only allowed as 'derived: rust_anchor' on fields the schema marks derivable"); return(NULL)
      }
      if (is.null(anchor)) return(NULL)
      v <- if (spec$derive == "anchor_alpha") anchor$mu * anchor$conc else (1 - anchor$mu) * anchor$conc
    } else {
      cv <- coerce_config_value(entry$value, spec)
      if (!is.null(cv$error)) { err(where, cv$error); return(NULL) }
      v <- cv$value
    }
    b <- check_config_bounds(v, spec); if (!is.null(b)) err(where, b)

    is_prior <- isTRUE(entry$prior)
    if (is_prior) {
      if (is.null(entry$reason) || !nzchar(entry$reason)) err(where, "a prior needs a reason")
      if (!is.null(entry$choices)) {
        if (!all(v %in% unlist(entry$choices))) err(where, "value is not among the prior's choices")
        priors[[length(priors) + 1L]] <<- data.frame(key = key, layer = layer, min = NA_real_, max = NA_real_,
          choices = paste(unlist(entry$choices), collapse = ", "), reason = entry$reason)
      } else {
        lo <- config_number(entry$min); hi <- config_number(entry$max)
        if (is.na(lo) || is.na(hi) || lo > hi) err(where, "a prior needs numeric min <= max (or choices)")
        else {
          if (is.numeric(v) && any(v < lo | v > hi)) err(where, sprintf("value outside its prior [%s, %s]", lo, hi))
          b <- check_config_bounds(c(lo, hi), spec); if (!is.null(b)) err(where, paste("prior range:", b))
          priors[[length(priors) + 1L]] <<- data.frame(key = key, layer = layer, min = lo, max = hi,
            choices = NA_character_, reason = entry$reason)
        }
      }
    }
    if (layer == "landscape" && !is_prior && is.null(entry$source))
      err(where, "a landscape value needs a 'source' or must be a prior")
    if (layer == "species" && !is_prior && is.null(entry$cite))
      err(where, "no citation: give 'cite: <csv file>' or mark it as a prior (CLAUDE.md rule 5)")

    if (check_sources && !is.null(entry$cite)) {
      tab <- csv_of(entry$cite)
      if (!is.data.frame(tab)) err(where, sprintf("cited file data/%s not found", entry$cite))
      else {
        row <- tab[tab$species == sp$species & tab$parameter == key & tab$role == "parameter", ]
        if (nrow(row) != 1L)
          err(where, sprintf("%d rows in data/%s for %s / %s (need exactly 1)", nrow(row), entry$cite, sp$species, key))
        else {
          if (length(v) == 1L && !config_values_match(row$value, v))
            err(where, sprintf("value %s differs from data/%s (%s)", format(v), entry$cite, row$value))
          pr <- csv_prior_range(row$notes)
          if (is_prior && !is.null(pr)) {
            same <- if (!is.null(pr$choices)) setequal(pr$choices, unlist(entry$choices))
                    else isTRUE(all.equal(c(pr$min, pr$max), c(config_number(entry$min), config_number(entry$max))))
            if (!same) err(where, sprintf("prior range differs from data/%s notes", entry$cite))
          }
          if (is_prior && is.null(pr)) err(where, sprintf("marked prior, but data/%s does not mark it PRIOR", entry$cite))
          if (!is_prior && !is.null(pr)) err(where, sprintf("data/%s marks it PRIOR; mark it prior here too", entry$cite))
        }
      }
    }
    v
  }

  # Optional fields (schema 'optional: true'; design_swappable_interfaces.md
  # section 2.1, notes T12). An optional field that is absent, or set to its
  # schema 'default', is left out of the list, so adding one to the schema
  # changes no existing species' list or run IDs. Model code reads it with the
  # default as fallback. 'required_when: {field: value}' makes it required
  # when another field (or that field's default) has that value.
  opt_default <- function(spec) {
    if (is.null(spec$default)) return(NULL)
    coerce_config_value(spec$default, spec)$value
  }
  is_default <- function(v, spec) !is.null(spec$default) && identical(v, opt_default(spec))
  for (lf in c(schema_leaves(schema$rust_anchor), schema_leaves(schema$params)))
    if (isTRUE(lf$spec$optional) && is.null(lf$spec$default) && is.null(lf$spec$required_when))
      err(paste(lf$path, collapse = "$"), "schema: an optional field needs a 'default' or 'required_when'")

  # --- Rust anchor -----------------------------------------------------------
  # A property of place (notes T18): site, cell and level from the landscape
  # file; the species file may only choose the mode.
  anchor <- list()
  anchor_src <- list(species   = list(data = sp$rust_anchor,   file = paste0("species/", species, ".yaml")),
                     landscape = list(data = land$rust_anchor, file = "config/landscape.yaml"))
  for (lf in schema_leaves(schema$rust_anchor)) {
    src <- anchor_src[[lf$spec$layer]]
    if (is.null(src)) { err(paste(c("rust_anchor", lf$path), collapse = "$"), "schema: bad layer"); next }
    where <- paste0(src$file, ": rust_anchor$", paste(lf$path, collapse = "$"))
    e <- cfg_get(src$data, lf$path)
    if (is.null(e)) { if (!isTRUE(lf$spec$optional)) err(where, "missing"); next }
    v <- read_entry(e, lf$spec, where, paste(c("rust_anchor", lf$path), collapse = "__"), lf$spec$layer)
    if (!is.null(v) && !is_default(v, lf$spec)) anchor <- cfg_set(anchor, lf$path, v)
  }
  anchor_ok <- all(c("mu", "conc") %in% names(anchor))

  # --- Params, in schema order ------------------------------------------------
  leaves <- schema_leaves(schema$params)
  params <- list()
  fname  <- c(species = paste0("species/", species, ".yaml"),
              landscape = "config/landscape.yaml", run = "config/run_defaults.yaml")
  for (lf in leaves) {
    layer <- lf$spec$layer
    if (is.null(layer) || !layer %in% names(files)) { err(paste(lf$path, collapse = "$"), "schema: bad layer"); next }
    where <- paste0(fname[[layer]], ": ", paste(lf$path, collapse = "$"))
    e <- cfg_get(files[[layer]], lf$path)
    if (is.null(e)) { if (!isTRUE(lf$spec$optional)) err(where, "missing"); next }   # an empty list ([]) is present, not NULL
    v <- if (layer == "run") {
      cv <- coerce_config_value(e, lf$spec)
      if (!is.null(cv$error)) { err(where, cv$error); NULL } else {
        b <- check_config_bounds(cv$value, lf$spec); if (!is.null(b)) err(where, b); cv$value }
    } else read_entry(e, lf$spec, where, paste(lf$path, collapse = "__"), layer,
                      anchor = if (anchor_ok) anchor else NULL)
    if (!is.null(v) && isTRUE(lf$spec$optional) && is_default(v, lf$spec)) next
    if (!is.null(v) || layer == "run") params <- cfg_set(params, lf$path, v)
  }
  # Optional fields that one option of another field needs (e.g. ricker_r when
  # recruitment_curve is ricker): required then, and refused otherwise.
  for (lf in leaves) {
    rw <- lf$spec$required_when
    if (is.null(rw)) next
    f    <- names(rw)[1]
    cur  <- cfg_get(params, f)
    if (is.null(cur)) cur <- opt_default(schema$params[[f]])
    need <- identical(cur, rw[[f]])
    have <- !is.null(cfg_get(files[[lf$spec$layer]], lf$path))
    where <- paste0(fname[[lf$spec$layer]], ": ", paste(lf$path, collapse = "$"))
    if (need && !have) err(where, sprintf("missing: needed when %s is '%s'", f, rw[[f]]))
    if (!need && have) err(where, sprintf("only used when %s is '%s'; remove it", f, rw[[f]]))
  }

  # --- Fields the schema does not know, or in the wrong file -----------------
  known <- vapply(leaves, function(lf) paste(c(lf$spec$layer, lf$path), collapse = "$"), "")
  for (layer in names(files)) {
    leaf_fn <- if (layer == "run") function(x) !is.list(x) || length(x) == 0L else is_entry
    for (p in config_paths(files[[layer]], leaf_fn)) {
      if (!length(p)) next
      if (!paste(c(layer, p), collapse = "$") %in% known)
        err(paste0(fname[[layer]], ": ", paste(p, collapse = "$")),
            "not a field of this layer in config/schema.yaml")
    }
  }
  for (layer in names(anchor_src)) {
    ok <- vapply(Filter(function(lf) identical(lf$spec$layer, layer), schema_leaves(schema$rust_anchor)),
                 function(lf) paste(lf$path, collapse = "$"), "")
    for (p in config_paths(anchor_src[[layer]]$data, is_entry))
      if (length(p) && !paste(p, collapse = "$") %in% ok)
        err(paste0(anchor_src[[layer]]$file, ": rust_anchor$", paste(p, collapse = "$")),
            if (layer == "species") "not a species setting: the rust anchor's site, cell and level are a property of place, set in config/landscape.yaml (notes T18)"
            else "not a landscape field of rust_anchor in config/schema.yaml")
  }
  extra <- setdiff(names(land), c("rust_anchor", "params", "spatial"))
  if (length(extra)) err("config/landscape.yaml", paste("unknown top-level keys:", paste(extra, collapse = ", ")))
  extra <- setdiff(names(sp), c("species", "rust_anchor", "params", "assessment", "spatial", "ensemble_draws"))
  if (length(extra)) err(fname[["species"]], paste("unknown top-level keys:", paste(extra, collapse = ", ")))
  if (is.null(sp$species) || !is.character(sp$species)) err(fname[["species"]], "needs 'species: <name>'")

  # --- Cross-field checks -----------------------------------------------------
  nl <- vapply(params[c("resist_locus_effect", "resist_dominance", "resist_freq0")], length, 1L)
  if (length(unique(nl)) > 1L) err("resistance", "resist_locus_effect, resist_dominance and resist_freq0 need the same length (one per locus)")
  if (identical(params$rust_flower_driver, "damage") && !identical(params$rust_hazard_mode, "damage"))
    err(paste0(fname[["species"]], ": rust_flower_driver"), "'damage' needs rust_hazard_mode: damage (notes T21)")

  # --- Assessment, spatial range and parameter draws (notes T9) ---------------
  af   <- read_config_yaml(file.path(root, "config", "assessment.yaml"))
  spf  <- function(block) list(data = sp[[block]], label = paste0(fname[["species"]], ": ", block))
  assessment <- read_settings(schema$assessment,
    list(assessment = list(data = af$assessment, label = "config/assessment.yaml: assessment"),
         species    = spf("assessment")), err)
  extra <- setdiff(names(af), c("assessment", "categories", "window_checks"))
  if (length(extra)) err("config/assessment.yaml", paste("unknown top-level keys:", paste(extra, collapse = ", ")))
  if (length(assessment) == length(schema_leaves(schema$assessment))) {
    a <- assessment
    if (!(a$year_spinup < a$year_gentime_start && a$year_gentime_start <= a$year_gentime_end &&
          a$year_gentime_end < a$year_rust && a$year_rust <= a$year_assess && a$year_assess < a$year_end))
      err("assessment", "years must run spin-up < gen-time start <= gen-time end < rust <= present < end")
  }
  cats <- read_settings(schema$categories,
    list(assessment = list(data = af$categories, label = "config/assessment.yaml: categories")), err)
  categories <- unlist(cats)
  if (length(categories) == 3L && is.unsorted(categories, strictly = TRUE))
    err("config/assessment.yaml: categories", "thresholds must increase: VU < EN < CR")
  window_checks <- read_settings(schema$window_checks,
    list(assessment = list(data = af$window_checks, label = "config/assessment.yaml: window_checks")), err)
  spatial <- read_settings(schema$spatial, list(species = spf("spatial")), err)
  if (!is.null(spatial$range_box$lon_min) && !is.null(spatial$range_box$lat_min) &&
      (spatial$range_box$lon_min >= spatial$range_box$lon_max || spatial$range_box$lat_min >= spatial$range_box$lat_max))
    err(paste0(fname[["species"]], ": spatial$range_box"), "needs lon_min < lon_max and lat_min < lat_max")

  draws <- list()
  draw_fields <- c("param", "min", "max", "label", "note", "cite")
  for (nm in names(sp$ensemble_draws)) {
    d     <- sp$ensemble_draws[[nm]]
    where <- paste0(fname[["species"]], ": ensemble_draws$", nm)
    if (!is.list(d) || !all(draw_fields %in% names(d))) {
      err(where, paste("needs", paste(draw_fields, collapse = ", "))); next
    }
    if (length(setdiff(names(d), draw_fields))) err(where, paste("unknown fields:", paste(setdiff(names(d), draw_fields), collapse = ", ")))
    path <- strsplit(d$param, "__", fixed = TRUE)[[1]]
    spec <- cfg_get(schema$params, path)
    if (is.null(spec$type) || spec$type != "double" ||
        !(identical(spec$layer, "species") || isTRUE(spec$drawable))) {
      err(where, sprintf("'%s' is not a numeric species parameter or a drawable run setting", d$param)); next
    }
    lo <- config_number(d$min); hi <- config_number(d$max)
    if (is.na(lo) || is.na(hi) || lo > hi) { err(where, "needs numeric min <= max"); next }
    b <- check_config_bounds(c(lo, hi), spec); if (!is.null(b)) err(where, b)
    v <- cfg_get(params, path)
    if (is.numeric(v) && (v < lo || v > hi)) err(where, sprintf("the parameter's value %s is outside [%s, %s]", v, lo, hi))
    if (check_sources) {
      tab <- csv_of(d$cite)
      row <- if (is.data.frame(tab)) tab[tab$species == sp$species & tab$parameter == d$param & tab$role == "parameter", ] else NULL
      if (is.null(row) || nrow(row) != 1L) err(where, sprintf("no single row for %s in data/%s", d$param, d$cite))
      else if (!isTRUE(all.equal(csv_draw_range(row$notes), c(lo, hi))))
        err(where, sprintf("range differs from 'RANGE [min, max]' in the notes of data/%s", d$cite))
    }
    draws[[nm]] <- list(path = path, min = lo, max = hi, label = d$label, note = d$note)
  }

  # The range box must lie inside the master grid (notes T13).
  ext <- tryCatch(load_landscape_spatial(root)$grid_extent, error = function(e) NULL)
  b   <- spatial$range_box
  if (!is.null(ext) && length(b) == 4L &&
      (b$lon_min < ext$lon_min || b$lon_max > ext$lon_max || b$lat_min < ext$lat_min || b$lat_max > ext$lat_max))
    err(paste0(fname[["species"]], ": spatial$range_box"),
        "must lie inside the master grid (config/landscape.yaml, spatial$grid_extent)")

  if (length(errs)) config_error(sprintf("species '%s'", species), errs)

  list(species = sp$species, params = params, rust_anchor = anchor,
       priors = if (length(priors)) do.call(rbind, priors) else NULL,
       assessment = assessment, categories = categories, window_checks = window_checks,
       spatial = c(list(name = sp$species, id = species), spatial,
                   list(data_dir = file.path("spatial", "data", species),
                        fig_dir  = file.path("spatial", "figures", species))),
       draws = draws)
}

load_species          <- function(species = NULL, ...) species_config(species, ...)$params
load_rust_anchor      <- function(species = NULL, ...) species_config(species, ...)$rust_anchor
species_priors        <- function(species = NULL, ...) species_config(species, ...)$priors
load_assessment       <- function(species = NULL, ...) species_config(species, ...)$assessment
load_categories       <- function(species = NULL, ...) species_config(species, ...)$categories
load_window_checks    <- function(species = NULL, ...) species_config(species, ...)$window_checks
load_param_draws      <- function(species = NULL, ...) species_config(species, ...)$draws
load_species_spatial  <- function(species = NULL, ...) species_config(species, ...)$spatial

# Shared landscape settings for the spatial pipeline (config/landscape.yaml,
# 'spatial'). The MN_SDM_TIFF environment variable overrides rust_sdm.
load_landscape_spatial <- function(root = mnsim_root()) {
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  lf     <- read_config_yaml(file.path(root, "config", "landscape.yaml"))
  errs   <- character()
  err    <- function(where, msg) errs <<- c(errs, sprintf("%s: %s", where, msg))
  out    <- read_settings(schema$landscape_spatial,
              list(landscape = list(data = lf$spatial, label = "config/landscape.yaml: spatial")), err)
  e <- out$grid_extent
  if (!is.null(e$lon_min) && !is.null(e$lat_min) && !is.null(out$grid_deg)) {
    if (e$lon_min >= e$lon_max || e$lat_min >= e$lat_max)
      err("config/landscape.yaml: spatial$grid_extent", "needs lon_min < lon_max and lat_min < lat_max")
    n <- c(e$lon_max - e$lon_min, e$lat_max - e$lat_min) / out$grid_deg
    if (any(abs(n - round(n)) > 1e-6))
      err("config/landscape.yaml: spatial$grid_extent", "must be a whole number of grid_deg cells wide and high")
  }
  if (length(errs)) config_error("config/landscape.yaml", errs)
  out$rust_sdm <- Sys.getenv("MN_SDM_TIFF", out$rust_sdm)
  out
}

# A species' validation set-up (layer V): species/<sp>_validation.yaml, one
# block per validation in config/schema.yaml ('validation'). NULL if the
# species has no validation file. The Pegg block also gets the assessment's
# year_rust (one rust year for everything; notes T4).
validation_config <- function(species = NULL, root = mnsim_root()) {
  if (is.null(species)) species <- default_species(root)
  path <- file.path(root, "species", paste0(species, "_validation.yaml"))
  if (!file.exists(path)) return(NULL)
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  v      <- read_config_yaml(path)
  name   <- read_config_yaml(file.path(root, "species", paste0(species, ".yaml")))$species
  label  <- paste0("species/", species, "_validation.yaml")
  errs   <- character()
  err    <- function(where, msg) errs <<- c(errs, sprintf("%s: %s", where, msg))
  if (!identical(v$species, name)) err(label, sprintf("'species' must be \"%s\", as in species/%s.yaml", name, species))
  extra <- setdiff(names(v), c("species", names(schema$validation)))
  if (length(extra)) err(label, paste("unknown top-level keys:", paste(extra, collapse = ", ")))
  csv_of <- csv_reader(root)
  out <- list()
  for (blk in intersect(names(schema$validation), names(v)))
    out[[blk]] <- read_settings(schema$validation[[blk]],
                   list(validation = list(data = v[[blk]], label = paste0(label, ": ", blk))),
                   err, csv_of, name)
  if (length(errs)) config_error(label, errs)
  if (!is.null(out$pegg)) out$pegg$year_rust <- load_assessment(species, root = root)$year_rust
  out
}
