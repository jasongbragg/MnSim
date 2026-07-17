# diagnostics/diag_fire_flowering.R
#
# Diagnostic: fire fates and return-to-flowering after fire, under the
# current model parameterisation.
#
# PANEL A -- Fire fate proportions by age at fire
#   For each age 0-max_age, what fraction of alive plants are killed,
#   enter resprout, or are unaffected? Shows the two-stage fire model
#   (fire_p_fimp > 0) vs the legacy flat model (fire_p_fimp = 0) and
#   makes the age-selectivity of the Hill-function immediately visible.
#   Calls hill_weight() from R/mortality.R directly.
#
# PANEL B -- Expected years to next flowering after fire, by age at fire
#   For a survivor of fire: how many years until it can flower again?
#   Two components:
#     (1) the resprout window (resprout_yrs_base, Poisson mean)
#     (2) the remaining juvenile period if burned before AFR
#         (max(0, AFR - fire_age - 1), since age increments to fire_age+1
#          in the fire year before mortality/resprout assignment)
#   A vertical dashed line marks age_first_flower_mean. Lines show no
#   rust vs susceptible plant with rust (the delay term is driven by
#   rust_pressure via dose_response(), NOT by annual_env_t -- the resprout
#   delay is not part of the disease triangle, only mortality and
#   fecundity suppression are). Calls dose_response() from R/mortality.R.
#
# PANEL C -- Proportion of resprouting survivors flowering, year by year
#   Among plants that survived a fire and entered the resprout state,
#   what fraction are flowering by year k post-fire? Uses the Poisson CDF
#   (ppois from base R) of the window-length distribution, conditioned on
#   the age-at-first-flower constraint. One panel per entry in fire_ages
#   (default: one juvenile, one adult). Calls dose_response() for the
#   rust delay. Optional overlay of observed data points.
#
# USAGE (from project root):
#   source("diagnostics/diag_fire_flowering.R")
#   diag_fire_flowering()
#   diag_fire_flowering(save_path = "outputs/fire_flowering.png")
#
#   # Supply observed proportion-flowering data for Panel C
#   diag_fire_flowering(
#     fire_ages = c(5, 25),
#     obs_data  = list(
#       data.frame(year_since_fire = 0:5, prop_flower = c(0,0,0.1,0.4,0.7,0.9)),
#       data.frame(year_since_fire = 0:5, prop_flower = c(0,0.2,0.6,0.85,0.95,0.98))
#     )
#   )

if (file.exists("params.R")) {
  .root <- "."
} else if (file.exists("../params.R")) {
  .root <- ".."
} else {
  stop("Source from project root or diagnostics/ directory.")
}
source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))   # hill_weight(), dose_response()

# --- Panel A: fire fate proportions ----------------------------------------
.fire_fates <- function(ages, params) {
  do.call(rbind, lapply(ages, function(age) {
    if (params$fire_p_fimp > 0) {
      # Two-stage model: calls hill_weight() from R/mortality.R
      ceiling <- params$fire_kill_scalar * params$fire_p_fimp
      p_kill_imp <- ceiling + (1 - ceiling) *
                    hill_weight(age, params$fire_kill_half_sat,
                                params$fire_kill_hill)
      p_killed    <- params$fire_p_fimp * p_kill_imp
      p_resprout  <- params$fire_p_fimp * (1 - p_kill_imp)
      p_unaffected <- 1 - params$fire_p_fimp
    } else {
      # Legacy flat-probability model
      p_killed    <- params$fire_kill_prob
      p_resprout  <- 1 - params$fire_kill_prob
      p_unaffected <- 0
    }
    data.frame(age = age, Killed = p_killed, Resprout = p_resprout,
               Unaffected = p_unaffected, stringsAsFactors = FALSE)
  }))
}

# --- Panel B: expected years to next flowering after fire ------------------
# Rust delay: driven by rust_pressure via dose_response(), NOT annual_env_t.
# For a plant with resist_score r: extra_delay = (1-r) * dose_response(...)
# Note apply_fire() enforces pmax(total_delay, 1L), so minimum window = 1 yr.
.years_to_flower <- function(ages, params,
                              resist_scores = c(0, 1)) {
  AFR <- params$age_first_flower_mean

  # Rust delay for a fully susceptible plant (resist=0); scaled by resist below
  delay_max <- if (!is.infinite(params$rust_start_year)) {
    dose_response(params$rust_pressure,    # R/mortality.R
                  params$rust_dose_response$delay)
  } else { 0 }

  do.call(rbind, lapply(resist_scores, function(r) {
    extra_delay <- (1 - r) * delay_max
    window_expected <- params$resprout_yrs_base + extra_delay
    lbl <- if (r == 0 && extra_delay > 0) "Susceptible (resist=0)"
           else if (r == 1)               "Resistant (resist=1)"
           else paste0("resist=", r)

    data.frame(
      fire_age = ages,
      years    = pmax(window_expected,
                      pmax(0, AFR - ages - 1)),
      label    = lbl,
      stringsAsFactors = FALSE
    )
  }))
}

# --- Panel C: proportion of resprouting survivors flowering, by year -------
# P(flowering by year k | entered resprout) =
#     ppois(k, window_lambda) * I(fire_age + 1 + k >= AFR)
# where window_lambda = resprout_yrs_base + extra_rust_delay.
# The pmax(delay,1) floor in apply_fire() means no recovery in year 0
# is possible -- ppois(0, lambda) ~ exp(-lambda) is already tiny for
# lambda >= 3, but we set year-0 to 0 explicitly for clarity.
.prop_flowering <- function(fire_ages, window_extra, params,
                             resist_scores = c(0, 1)) {
  AFR      <- params$age_first_flower_mean
  wlen     <- params$resprout_yrs_base

  delay_max <- if (!is.infinite(params$rust_start_year)) {
    dose_response(params$rust_pressure, params$rust_dose_response$delay)
  } else { 0 }

  do.call(rbind, lapply(fire_ages, function(fa) {
    years_seq    <- seq(-1L, wlen + window_extra)
    afr_wait     <- max(0L, AFR - fa - 1L)  # years post-fire to reach AFR

    do.call(rbind, lapply(resist_scores, function(r) {
      extra        <- (1 - r) * delay_max
      window_lam   <- wlen + extra
      lbl <- if (r == 0 && extra > 0) "Susceptible (resist=0)"
             else if (r == 1)          "Resistant (resist=1)"
             else paste0("resist=", r)

      p_vals <- sapply(years_seq, function(yr) {
        if (yr < 0L)  {
          # pre-fire: adults flower, juveniles do not
          if (fa >= AFR) 1.0 else 0.0
        } else if (yr == 0L) {
          0.0  # fire year: all survivors just entered resprout
        } else {
          age_ok <- yr >= afr_wait
          if (!age_ok) 0.0 else ppois(yr, window_lam)
        }
      })

      data.frame(fire_age = fa, year_since_fire = years_seq,
                 p_flower = p_vals, resist = r, label = lbl,
                 stringsAsFactors = FALSE)
    }))
  }))
}

# =============================================================================
# diag_fire_flowering()
# =============================================================================
diag_fire_flowering <- function(
    params       = get_default_params(),
    age_range    = 0:60,          # ages for Panels A and B
    fire_ages    = c(4, 25),      # reference ages for Panel C
    resist_scores = c(0, 1),      # resist_scores for Panels B and C
    window_extra = 4,             # years beyond expected window in Panel C
    obs_data     = NULL,          # list of data.frames, one per fire_ages entry
                                  # columns: year_since_fire, prop_flower
    save_path    = NULL,
    width_in     = 15,
    height_in    = 5,
    res          = 150
) {
  A_df <- .fire_fates(age_range, params)
  B_df <- .years_to_flower(age_range, params, resist_scores)
  C_df <- .prop_flowering(fire_ages, window_extra, params, resist_scores)

  fire_mode  <- if (params$fire_p_fimp > 0) "two-stage" else "legacy (flat)"
  rust_on    <- !is.infinite(params$rust_start_year)
  delay_note <- if (!rust_on) "rust_start_year=Inf: delay=0 regardless of resist"
                else paste0("rust_pressure=", params$rust_pressure)

  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_fire_flowering(A_df, B_df, C_df, params, fire_ages, fire_mode,
                        delay_note, rust_on, obs_data)
  } else {
    .base_fire_flowering(A_df, B_df, C_df, params, fire_ages, fire_mode,
                          delay_note, rust_on, obs_data)
  }

  invisible(list(fates = A_df, years_to_flower = B_df, flowering = C_df))
}

# =============================================================================
# ggplot2 rendering
# =============================================================================
.gg_fire_flowering <- function(A_df, B_df, C_df, params, fire_ages, fire_mode,
                                delay_note, rust_on, obs_data) {

  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes
  AFR <- params$age_first_flower_mean
  wlen <- params$resprout_yrs_base

  # -- Palette ---------------------------------------------------------------
  fate_pal <- c(Killed = "#e34948", Resprout = "#2a78d6", Unaffected = "#888780")
  n_rs     <- length(unique(B_df$label))
  rs_pal   <- grDevices::colorRampPalette(c("#e34948", "#2a78d6"))(n_rs)
  names(rs_pal) <- unique(B_df$label)

  # -- Panel A: fate proportions ---------------------------------------------
  A_long <- utils::stack(A_df[, c("Killed","Resprout","Unaffected")])
  A_long$age  <- rep(A_df$age, 3)
  A_long$fate <- as.character(A_long$ind)
  A_long$fate <- factor(A_long$fate, levels = c("Unaffected","Resprout","Killed"))

  pA <- gg(A_long, aes(x = age, y = values, fill = fate)) +
    ggplot2::geom_area(position = "stack", alpha = 0.85) +
    ggplot2::scale_fill_manual(values = fate_pal, name = "Fate") +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = c(0, .02))) +
    ggplot2::geom_vline(xintercept = AFR, linetype = "dashed",
                         colour = "white", linewidth = 0.7) +
    ggplot2::annotate("text", x = AFR + 0.5, y = 0.97,
                       label = paste0("AFR=", AFR), hjust = 0, size = 2.8,
                       colour = "white") +
    ggplot2::labs(x = "Age at fire (years)", y = "Proportion of alive plants",
                  title = paste0("Fire fates by age  [", fire_mode, "]")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # -- Panel B: years to next flowering --------------------------------------
  B_df$label <- factor(B_df$label, levels = unique(B_df$label))

  pB <- gg(B_df, aes(x = fire_age, y = years, colour = label)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = AFR, linetype = "dashed",
                         colour = "#888780", linewidth = 0.7) +
    ggplot2::annotate("text", x = AFR + 0.5,
                       y = max(B_df$years) * 0.97,
                       label = paste0("AFR=", AFR),
                       hjust = 0, size = 2.8, colour = "#888780") +
    ggplot2::scale_colour_manual(values = rs_pal, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .06))) +
    ggplot2::labs(x = "Age at fire (years)",
                  y = "Expected years until next flowering",
                  title = paste0("Years to next flowering  [", delay_note, "]")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # -- Panel C: proportion flowering over time -------------------------------
  C_df$panel <- paste0("Burned at age ", C_df$fire_age)
  C_df$label <- factor(C_df$label, levels = unique(C_df$label))

  strip_df <- data.frame(
    panel = paste0("Burned at age ", fire_ages),
    xmin = 0, xmax = wlen
  )

  pC <- gg(C_df, aes(x = year_since_fire, y = p_flower, colour = label)) +
    ggplot2::geom_rect(data = strip_df,
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                        fill = "#d0e4f4", alpha = 0.4, colour = NA,
                        inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 0,    linetype = "dotted", colour = "#888") +
    ggplot2::geom_vline(xintercept = wlen, linetype = "dashed",
                         colour = "#333", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ panel, ncol = length(fire_ages)) +
    ggplot2::scale_colour_manual(values = rs_pal, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = c(0, .03)),
                                  labels = scales::percent_format(accuracy = 1)) +
    ggplot2::annotate("text", x = wlen + 0.15, y = 0.15,
                       label = paste0("mean\nwindow\n(yr ", wlen, ")"),
                       hjust = 0, size = 2.5, colour = "#333") +
    ggplot2::labs(x = "Year since fire  (shaded = expected resprout window)",
                  y = "Proportion of resprout survivors flowering",
                  title = "Return to flowering post-fire") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"))

  # optional: add scales:: check
  if (!requireNamespace("scales", quietly = TRUE)) {
    pC <- pC + ggplot2::scale_y_continuous(limits = c(0, 1),
                                              expand = ggplot2::expansion(mult = c(0, .03)))
  }

  # optional observed data overlay
  if (!is.null(obs_data)) {
    for (i in seq_along(fire_ages)) {
      if (i > length(obs_data) || is.null(obs_data[[i]])) next
      od       <- obs_data[[i]]
      od$panel <- paste0("Burned at age ", fire_ages[i])
      pC <- pC +
        ggplot2::geom_point(data = od, inherit.aes = FALSE,
                             aes(x = year_since_fire, y = prop_flower),
                             colour = "black", size = 2.5, shape = 16)
      if (all(c("prop_lo","prop_hi") %in% names(od))) {
        pC <- pC +
          ggplot2::geom_errorbar(data = od, inherit.aes = FALSE,
                                  aes(x = year_since_fire,
                                      ymin = prop_lo, ymax = prop_hi),
                                  width = 0.25, colour = "black")
      }
    }
  }

  # -- assemble --------------------------------------------------------------
  n_C <- length(fire_ages)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, pC,
                                 ncol = 2 + n_C,
                                 widths = c(1, 1, rep(1, n_C))))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, pC, ncol = 2 + n_C)
  } else {
    print(pA)
    readline("Enter for Panel B ... ")
    print(pB)
    readline("Enter for Panel C ... ")
    print(pC)
  }
}

# =============================================================================
# base R fallback
# =============================================================================
.base_fire_flowering <- function(A_df, B_df, C_df, params, fire_ages, fire_mode,
                                  delay_note, rust_on, obs_data) {

  n_panels <- 2L + length(fire_ages)
  op <- par(mfrow = c(1L, n_panels), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  AFR  <- params$age_first_flower_mean
  wlen <- params$resprout_yrs_base
  rs_labels <- unique(B_df$label)
  n_rs  <- length(rs_labels)
  rs_col <- grDevices::colorRampPalette(c("#e34948","#2a78d6"))(n_rs)

  # -- Panel A ---------------------------------------------------------------
  # Stacked proportions as filled polygons
  ages <- A_df$age
  plot(NA, xlim = range(ages), ylim = c(0, 1),
       xlab = "Age at fire (years)", ylab = "Proportion of alive plants",
       main = paste0("Fire fates  [", fire_mode, "]"), las = 1)
  # stack from bottom: killed, then resprout, then unaffected
  y0 <- rep(0, nrow(A_df))
  y1 <- A_df$Killed
  polygon(c(ages, rev(ages)), c(y0, rev(y1)),
          col = grDevices::adjustcolor("#e34948", 0.85), border = NA)
  y2 <- y1 + A_df$Resprout
  polygon(c(ages, rev(ages)), c(y1, rev(y2)),
          col = grDevices::adjustcolor("#2a78d6", 0.85), border = NA)
  polygon(c(ages, rev(ages)), c(y2, rev(rep(1, nrow(A_df)))),
          col = grDevices::adjustcolor("#888780", 0.5), border = NA)
  abline(v = AFR, lty = 2, col = "white")
  legend("right", bty = "n", cex = 0.75,
         legend = c("Killed","Resprout","Unaffected"),
         fill = c("#e34948","#2a78d6","#888780"))

  # -- Panel B ---------------------------------------------------------------
  ymax_B <- max(B_df$years, na.rm = TRUE) * 1.06
  plot(NA, xlim = range(B_df$fire_age), ylim = c(0, ymax_B),
       xlab = "Age at fire (years)", ylab = "Expected years to next flowering",
       main = paste0("Years to flowering\n[", delay_note, "]"), las = 1)
  abline(v = AFR, lty = 2, col = "#888780")
  for (j in seq_along(rs_labels)) {
    sub <- B_df[B_df$label == rs_labels[j], ]
    lines(sub$fire_age, sub$years, col = rs_col[j], lwd = 2)
  }
  legend("topright", bty = "n", cex = 0.75, legend = rs_labels,
         col = rs_col, lwd = 2)

  # -- Panel C (one per fire_age) -------------------------------------------
  yC <- c(0, 1)
  for (k in seq_along(fire_ages)) {
    fa   <- fire_ages[k]
    sub  <- C_df[C_df$fire_age == fa, ]
    xr   <- range(sub$year_since_fire)

    plot(NA, xlim = xr, ylim = yC,
         xlab = "Year since fire", ylab = "Proportion flowering",
         main = paste0("Return to flowering\nBurned at age ", fa), las = 1)
    rect(0, 0, wlen, 1, col = grDevices::adjustcolor("#2a78d6", 0.08), border = NA)
    abline(v = 0,    lty = 3, col = "#888")
    abline(v = wlen, lty = 2, col = "#333")
    mtext(paste0("mean window yr ", wlen), at = wlen, side = 3, cex = 0.6)

    for (j in seq_along(rs_labels)) {
      et_sub <- sub[sub$label == rs_labels[j], ]
      et_sub <- et_sub[order(et_sub$year_since_fire), ]
      lines(et_sub$year_since_fire, et_sub$p_flower, col = rs_col[j], lwd = 2)
    }

    if (!is.null(obs_data) && k <= length(obs_data) && !is.null(obs_data[[k]])) {
      od <- obs_data[[k]]
      points(od$year_since_fire, od$prop_flower, pch = 16, cex = 1.3)
      if (all(c("prop_lo","prop_hi") %in% names(od)))
        arrows(od$year_since_fire, od$prop_lo, od$year_since_fire, od$prop_hi,
               angle = 90, code = 3, length = 0.06)
    }
    legend("topleft", bty = "n", cex = 0.75, legend = rs_labels,
           col = rs_col, lwd = 2)
  }
}
