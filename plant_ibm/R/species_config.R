# R/species_config.R
#
# Species configuration: assembles the model's params list from YAML files
# and validates it (PLAN.md Phase 2; docs/methodological_notes.md T4, T7).
#
#   config/schema.yaml        every field: layer, type, bounds, order
#   species/<sp>.yaml         species biology + rust anchor (layer "species")
#   config/landscape.yaml     shared landscape defaults (layer "landscape")
#   config/run_defaults.yaml  run settings (layer "run")
#
# The assembled list follows the schema's order, so for M. nodosa it is
# identical() to the old get_default_params() in params.R (and gives the same
# make_run_id() hashes).
#
# Numbers: every YAML number is parsed by R's as.numeric() (yaml handlers), then
# coerced to the schema type, so 10000 is a double where the schema says double.
# A value may also be written as a fraction ("1/3"). '.inf' is Inf.
#
# Paths are relative to the repo root: MNSIM_ROOT, or the parent of the
# working directory (scripts run from plant_ibm/).
#
# Usage:
#   params <- load_species("m_nodosa")
#   anchor <- load_rust_anchor("m_nodosa")      # was RUST_ANCHOR
#   priors <- species_priors("m_nodosa")        # data frame, one row per prior

mnsim_root <- function() Sys.getenv("MNSIM_ROOT", normalizePath("..", mustWork = FALSE))

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

# Read one species' configuration; check_sources also checks citations and
# priors against the CSV files (CLAUDE.md rule 5).
species_config <- function(species = "m_nodosa", root = mnsim_root(), check_sources = TRUE) {
  schema <- read_config_yaml(file.path(root, "config", "schema.yaml"))
  sp     <- read_config_yaml(file.path(root, "species", paste0(species, ".yaml")))
  files  <- list(species   = sp$params,
                 landscape = read_config_yaml(file.path(root, "config", "landscape.yaml"))$params,
                 run       = read_config_yaml(file.path(root, "config", "run_defaults.yaml"))$params)
  errs   <- character()
  err    <- function(where, msg) errs <<- c(errs, sprintf("%s: %s", where, msg))
  csvs   <- list()
  csv_of <- function(f) {
    if (is.null(csvs[[f]])) {
      p <- file.path(root, "data", f)
      csvs[[f]] <<- if (file.exists(p)) utils::read.csv(p, colClasses = "character") else NA
    }
    csvs[[f]]
  }
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

  # --- Rust anchor -----------------------------------------------------------
  anchor <- list()
  for (lf in schema_leaves(schema$rust_anchor)) {
    where <- paste0("species/", species, ".yaml: rust_anchor$", paste(lf$path, collapse = "$"))
    e <- cfg_get(sp$rust_anchor, lf$path)
    if (is.null(e)) { err(where, "missing"); next }
    v <- read_entry(e, lf$spec, where, paste(c("rust_anchor", lf$path), collapse = "__"), "species")
    if (!is.null(v)) anchor <- cfg_set(anchor, lf$path, v)
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
    if (is.null(e)) { err(where, "missing"); next }   # an empty list ([]) is present, not NULL
    v <- if (layer == "run") {
      cv <- coerce_config_value(e, lf$spec)
      if (!is.null(cv$error)) { err(where, cv$error); NULL } else {
        b <- check_config_bounds(cv$value, lf$spec); if (!is.null(b)) err(where, b); cv$value }
    } else read_entry(e, lf$spec, where, paste(lf$path, collapse = "__"), layer,
                      anchor = if (anchor_ok) anchor else NULL)
    if (!is.null(v) || layer == "run") params <- cfg_set(params, lf$path, v)
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
  for (p in config_paths(sp$rust_anchor, is_entry))
    if (length(p) && !paste(p, collapse = "$") %in% vapply(schema_leaves(schema$rust_anchor),
                                                          function(lf) paste(lf$path, collapse = "$"), ""))
      err(paste0(fname[["species"]], ": rust_anchor$", paste(p, collapse = "$")), "not in config/schema.yaml")
  extra <- setdiff(names(sp), c("species", "rust_anchor", "params"))
  if (length(extra)) err(fname[["species"]], paste("unknown top-level keys:", paste(extra, collapse = ", ")))
  if (is.null(sp$species) || !is.character(sp$species)) err(fname[["species"]], "needs 'species: <name>'")

  # --- Cross-field checks -----------------------------------------------------
  nl <- vapply(params[c("resist_locus_effect", "resist_dominance", "resist_freq0")], length, 1L)
  if (length(unique(nl)) > 1L) err("resistance", "resist_locus_effect, resist_dominance and resist_freq0 need the same length (one per locus)")

  if (length(errs))
    stop(sprintf("Invalid configuration for species '%s' (%d problem%s):\n  - %s", species,
                 length(errs), if (length(errs) > 1L) "s" else "", paste(errs, collapse = "\n  - ")),
         call. = FALSE)

  list(species = sp$species, params = params, rust_anchor = anchor,
       priors = if (length(priors)) do.call(rbind, priors) else NULL)
}

load_species     <- function(species = "m_nodosa", ...) species_config(species, ...)$params
load_rust_anchor <- function(species = "m_nodosa", ...) species_config(species, ...)$rust_anchor
species_priors   <- function(species = "m_nodosa", ...) species_config(species, ...)$priors
