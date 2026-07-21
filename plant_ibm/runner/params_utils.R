# runner/params_utils.R
#
# Utilities for translating between the nested params list that
# run_simulation() expects and the flat named numeric vector that
# Morris, ensemble runners, and CSV files use.
#
# NAMING CONVENTION
# -----------------
# Nested list paths are encoded with double-underscore (__) separators.
# Splitting on __ and navigating with [[ handles arbitrary depth:
#   "weibull_lambda"                        -> params$weibull_lambda
#   "senescence_dose_response__max_effect"  -> params$senescence_dose_response$max_effect
#   "rust_dose_response__resprout__max_effect"
#                                           -> params$rust_dose_response$resprout$max_effect
# Top-level scalars have no __ at all.
#
# RUN IDENTIFIER
# --------------
# make_run_id() hashes the full params list to a 12-character hex string.
# Two calls with identical params (including params$seed) always produce
# the same ID; any change to any field produces a different one. The ID
# is used as the filename for cached result .rds files, making ensemble
# runs restartable and cross-machine result merging collision-free.

# --- vec_to_params -----------------------------------------------------------
# Overwrite specific fields in a base params list from a flat named numeric
# vector. The names of vec are __ keys; the values overwrite the
# corresponding nested fields. Fields not in vec are left unchanged.
vec_to_params <- function(base, vec) {
  p <- base
  for (k in names(vec)) {
    parts <- strsplit(k, "__", fixed = TRUE)[[1]]
    p     <- .set_nested(p, parts, as.numeric(vec[[k]]))
  }
  p
}

# Recursive helper: set value at arbitrary depth in a nested list.
.set_nested <- function(lst, parts, value) {
  if (length(parts) == 1L) {
    lst[[parts]] <- value
    return(lst)
  }
  lst[[parts[1]]] <- .set_nested(lst[[parts[1]]], parts[-1L], value)
  lst
}

# --- params_to_vec -----------------------------------------------------------
# Extract a named numeric vector of specified fields from a nested params
# list. keys is a character vector of __ names.
params_to_vec <- function(params, keys) {
  vapply(keys, function(k) {
    parts <- strsplit(k, "__", fixed = TRUE)[[1]]
    val   <- params
    for (p in parts) val <- val[[p]]
    as.numeric(val)
  }, numeric(1))
}

# --- make_run_id -------------------------------------------------------------
# Deterministic 12-character hex identifier for a params list.
# Uses digest::xxhash64 if available; falls back to a simple character
# sum otherwise (less collision-resistant but always available).
make_run_id <- function(params) {
  txt <- paste(deparse(params), collapse = "\n")
  if (requireNamespace("digest", quietly = TRUE)) {
    return(substr(digest::digest(txt, algo = "xxhash64"), 1L, 12L))
  }
  # Fallback: modular sum of UTF-8 code points, formatted as 12-hex
  chars <- utf8ToInt(substr(txt, 1L, 10000L))
  sprintf("%012x", abs(sum(as.numeric(chars))) %% (16L^12L))
}
