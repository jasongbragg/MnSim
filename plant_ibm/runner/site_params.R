# runner/site_params.R
#
# Shared helpers for site-specific runs: the per-cell table written by
# spatial/02_cross_reference_layers.R and how a cell's values override the
# base parameters. Used by the range-wide ensemble, the talk figures and the
# Pegg et al. validation (which runs at the Yarringully anchor cell).

SITE_PARAMS_PATH <- "spatial/data/mn_site_params.csv"

# Read the cell table; optionally drop cells lacking rust or fire values.
load_site_params <- function(path = SITE_PARAMS_PATH, drop_incomplete = TRUE) {
  s <- read.csv(path)
  if (drop_incomplete) {
    ok <- !is.na(s$annual_env_alpha) & !is.na(s$annual_env_beta) &
          !is.na(s$fire_prob_annual) & !is.na(s$fire_p_fimp)
    attr_n <- sum(!ok); attr_w <- sum(s$weight_n_raw[!ok]) / sum(s$weight_n_raw)
    s <- s[ok, ]
    attr(s, "n_dropped") <- attr_n; attr(s, "weight_dropped") <- attr_w
  }
  s
}

# The row of the rust anchor cell (flagged by 02 / 02c; nearest cell to the
# anchor coordinates as a fallback for older tables).
anchor_site_row <- function(sites, anchor = RUST_ANCHOR) {
  if ("rust_anchor_cell" %in% names(sites) && any(sites$rust_anchor_cell %in% TRUE))
    return(sites[which(sites$rust_anchor_cell %in% TRUE)[1], ])
  d2 <- (sites$lon_centre - anchor$lon)^2 + (sites$lat_centre - anchor$lat)^2
  sites[which.min(d2), ]
}

# Fire regimes a run can use (docs/methodological_notes.md D5):
#   "wild" (default)  wildfire only: the standard columns of the site table
#   "all"             all recorded fire (wildfire + planned burns): *_all columns
# 'fire_scale' multiplies the annual fire probability (e.g. 0.5 = wildfire
# suppression halving fire frequency; 1 = as recorded).
FIRE_REGIME_CHOICES <- c("wild", "all")

site_fire_values <- function(site_row, fire_regime = "wild") {
  fire_regime <- match.arg(fire_regime, FIRE_REGIME_CHOICES)
  sfx <- if (fire_regime == "wild") "" else paste0("_", fire_regime)
  cols <- paste0(c("fire_prob_annual", "fire_p_fimp", "sd_frac_burned"), sfx)
  miss <- setdiff(cols, names(site_row))
  if (length(miss))
    stop("site table has no ", paste(miss, collapse = ", "),
         ": rebuild it with spatial/02_cross_reference_layers.R (fire regime '", fire_regime, "')")
  list(prob = site_row[[cols[1]]], p_fimp = site_row[[cols[2]]], size_sd = site_row[[cols[3]]])
}

# Override base params with one cell's values.
build_site_params <- function(base_params, site_row, fire_regime = "wild", fire_scale = 1) {
  p <- base_params

  # Site-specific fire regime
  fv <- site_fire_values(site_row, fire_regime)
  stopifnot(is.numeric(fire_scale), length(fire_scale) == 1L, fire_scale >= 0)
  p$fire_prob_annual <- min(1, fv$prob * fire_scale)
  p$fire_p_fimp      <- fv$p_fimp
  # Between-fire sd of burnt fraction: the cell's own where it had >= 2
  # fire years, else the base value (range-wide median).
  if (!is.null(fv$size_sd) && !is.na(fv$size_sd))
    p$fire_size_sd   <- fv$size_sd
  p$fire_years       <- integer(0)   # must use fire_prob_annual for long runs

  # Site-specific rust pressure (SDM, normalised at the rust anchor cell)
  p$annual_env_alpha <- site_row$annual_env_alpha
  p$annual_env_beta  <- site_row$annual_env_beta
  # microclim_alpha/beta held constant (within-site individual variation,
  # assumed similar across all M. nodosa sites)

  # --- RESISTANCE HOOK -------------------------------------------------------
  # When site-level resistance data are available, add a resist_freq0 column
  # to mn_site_params.csv and uncomment:
  # p$resist_freq0 <- c(site_row$resist_freq0)
  # For now: resist_freq0 from params.R (fully susceptible everywhere)
  # ---------------------------------------------------------------------------

  p
}
