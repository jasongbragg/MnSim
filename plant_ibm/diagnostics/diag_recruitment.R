# diagnostics/diag_recruitment.R
#
# Diagnostic: recruitment function shape, canopy suppression, and rust
# fecundity effect -- calling bevholt(), shade_suppression(), and
# dose_response() from the model code directly.
#
# PANEL A -- Beverton-Holt curve in isolation
#   bevholt(n_flowering, params) vs n_flowering. Asymptote at R_max,
#   half-saturation at K_half, Poisson ±1 SD band. A vertical line marks
#   the equilibrium n_flowering (from a short fire-free run, or supplied
#   by the user). Answers: is the population operating in the linear regime
#   (below K_half, every adult matters) or saturated (above K_half, losses
#   of adults barely affect recruitment)?
#
# PANEL B -- BH curve × canopy shade suppression
#   Expected recruits = bevholt × (1 - shade_suppression(N_canopy)) for
#   N_canopy = {0, low, medium, high}. N_canopy = 0 recovers the raw BH
#   curve (Panel A). The equilibrium operating point (n_flowering_ref,
#   N_canopy_ref) is marked as a point. Answers: does the canopy effect
#   dominate recruitment regulation, or does BH saturation?
#
# PANEL C -- Shade suppression calibration
#   shade_suppression(N_canopy) from 0 to a maximum canopy density.
#   Dashed lines at half_sat (50% suppression point) and at the
#   equilibrium N_canopy_ref. Answers the calibration question directly:
#   "at my typical canopy density, how much is establishment suppressed?"
#
# PANEL D -- Rust fecundity multiplier
#   Per-adult fecundity weight = 1 - dose_response(eff_susc, cfg) vs
#   effective_susceptibility. Flat line at 1.0 when max_effect = 0 (off),
#   labelled accordingly. Calls dose_response() from R/mortality.R.
#   Shows the parameter exists and what turning it on would do.
#
# EQUILIBRATION RUN
#   By default (run_equilibration = TRUE) a short fire-free simulation is
#   run internally to find where the population settles, setting
#   n_flowering_ref and N_canopy_ref. Fire is turned off for this run (to
#   avoid moving the canopy during the window); all mortality and rust
#   terms use the supplied params as-is. Override by passing
#   n_flowering_ref and/or N_canopy_ref directly.
#
# USAGE (from project root, after loading params):
#   source("diagnostics/diag_recruitment.R")
#   diag_recruitment()
#   diag_recruitment(save_path = "outputs/recruitment.png")
#
#   # Skip equilibration; supply reference point from field knowledge
#   diag_recruitment(run_equilibration = FALSE,
#                    n_flowering_ref = 800, N_canopy_ref = 1200)
#
#   # Show rust fecundity effect by turning it on temporarily
#   p <- get_default_params()
#   p$rust_recruit_dose_response$max_effect <- 0.6
#   diag_recruitment(params = p)

if (file.exists("params.R")) {
  .root <- "."
} else if (file.exists("../params.R")) {
  .root <- ".."
} else {
  stop("Source from project root or diagnostics/ directory.")
}
source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))     # dose_response(), hill_weight()
source(file.path(.root, "R/genetics.R"))
source(file.path(.root, "R/individuals.R"))
source(file.path(.root, "R/recruitment.R"))   # bevholt(), shade_suppression(),
                                               # canopy_density()
source(file.path(.root, "R/fire.R"))
source(file.path(.root, "R/census.R"))
source(file.path(.root, "R/simulate.R"))      # run_simulation() for equilibration

# =============================================================================
# diag_recruitment()
# =============================================================================
diag_recruitment <- function(
    params            = get_default_params(),
    n_flowering_max   = NULL,    # NULL: 3 * K_half
    N_canopy_values   = NULL,    # NULL: {0, ref*0.3, ref*0.7, ref} from equil
    n_flowering_ref   = NULL,    # NULL: from equilibration run
    N_canopy_ref      = NULL,    # NULL: from equilibration run
    run_equilibration = TRUE,    # run a short fire-free sim to find ref point
    equil_years       = 100L,
    save_path         = NULL,
    width_in          = 12,
    height_in         = 10,
    res               = 150
) {

  # --- equilibration run (fire-free) to find operating point ---------------
  if (run_equilibration && (is.null(n_flowering_ref) || is.null(N_canopy_ref))) {
    message("Running ", equil_years, "-year fire-free equilibration...")
    p_eq                  <- params
    p_eq$fire_years       <- integer(0)
    p_eq$fire_prob_annual <- 0
    p_eq$n_years          <- equil_years
    res_eq <- run_simulation(p_eq, year0 = 1L, verbose = FALSE)
    pop_eq <- res_eq$individuals
    if (is.null(n_flowering_ref))
      n_flowering_ref <- sum(pop_eq$alive & pop_eq$flowering)
    if (is.null(N_canopy_ref))
      N_canopy_ref <- canopy_density(pop_eq)
    message(sprintf("  Equilibrium: n_flowering = %d,  N_canopy = %.0f",
                    n_flowering_ref, N_canopy_ref))
  }

  ref_known <- !is.null(n_flowering_ref) && !is.null(N_canopy_ref)

  # --- derived axis limits -------------------------------------------------
  n_max  <- if (!is.null(n_flowering_max)) n_flowering_max
            else 3 * params$K_half
  nc_max <- if (!is.null(N_canopy_ref))    N_canopy_ref * 2
            else params$shade_dose_response$half_sat * 3

  # --- Panel A data: BH curve ----------------------------------------------
  n_seq  <- seq(0, n_max, length.out = 300)
  bh_exp <- sapply(n_seq, bevholt, params = params)  # R/recruitment.R
  bh_lo  <- pmax(0, bh_exp - sqrt(bh_exp))
  bh_hi  <- bh_exp + sqrt(bh_exp)

  A_df <- data.frame(n = n_seq, expected = bh_exp, lo = bh_lo, hi = bh_hi)

  # --- Panel B data: BH × shade at several N_canopy values ----------------
  if (is.null(N_canopy_values)) {
    ref <- if (!is.null(N_canopy_ref)) N_canopy_ref
           else params$shade_dose_response$half_sat
    N_canopy_values <- round(c(0, ref * 0.3, ref * 0.7, ref))
  }
  N_canopy_labels <- ifelse(N_canopy_values == 0, "N_canopy = 0 (open)",
                             paste0("N_canopy = ", format(N_canopy_values,
                                                           big.mark = ",")))

  B_df <- do.call(rbind, lapply(seq_along(N_canopy_values), function(i) {
    nc  <- N_canopy_values[i]
    sup <- shade_suppression(nc, params)    # R/recruitment.R
    data.frame(
      n        = n_seq,
      expected = bh_exp * (1 - sup),
      N_canopy = N_canopy_labels[i],
      nc_val   = nc,
      stringsAsFactors = FALSE
    )
  }))
  B_df$N_canopy <- factor(B_df$N_canopy, levels = N_canopy_labels)

  # --- Panel C data: shade suppression vs N_canopy -------------------------
  nc_seq   <- seq(0, nc_max, length.out = 300)
  shade_df <- data.frame(
    N_canopy   = nc_seq,
    suppressed = sapply(nc_seq, shade_suppression, params = params)  # R/recruitment.R
  )

  # --- Panel D data: rust fecundity multiplier -----------------------------
  cfg <- params$rust_recruit_dose_response
  es_seq  <- seq(0, 1, length.out = 200)
  rust_wt <- 1 - sapply(es_seq, dose_response, cfg = cfg)  # R/mortality.R
  rust_off <- cfg$max_effect == 0
  D_df <- data.frame(eff_susc = es_seq, weight = rust_wt)

  # --- flags ---------------------------------------------------------------
  shade_on <- params$shade_dose_response$max_effect > 0
  hs       <- params$shade_dose_response$half_sat

  # --- render --------------------------------------------------------------
  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_recruitment(A_df, B_df, shade_df, D_df, params,
                    n_flowering_ref, N_canopy_ref, ref_known,
                    shade_on, hs, rust_off, nc_max)
  } else {
    .base_recruitment(A_df, B_df, shade_df, D_df, params,
                      n_flowering_ref, N_canopy_ref, ref_known,
                      shade_on, hs, rust_off, nc_max)
  }

  invisible(list(bh_curve    = A_df,
                 bh_shaded   = B_df,
                 shade_curve = shade_df,
                 rust_weight = D_df,
                 n_flowering_ref = n_flowering_ref,
                 N_canopy_ref    = N_canopy_ref))
}

# =============================================================================
# ggplot2 rendering
# =============================================================================
.gg_recruitment <- function(A_df, B_df, shade_df, D_df, params,
                             n_flowering_ref, N_canopy_ref, ref_known,
                             shade_on, hs, rust_off, nc_max) {
  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes

  R_max   <- params$R_max
  K_half  <- params$K_half
  nc_pal  <- grDevices::colorRampPalette(c("#d9f0d3","#1b7837"))(nlevels(B_df$N_canopy))

  # -- Panel A: BH curve ----------------------------------------------------
  pA <- gg(A_df, aes(x = n, y = expected)) +
    ggplot2::geom_ribbon(aes(ymin = lo, ymax = hi),
                          fill = "#2a78d6", alpha = 0.15) +
    ggplot2::geom_line(colour = "#2a78d6", linewidth = 1) +
    ggplot2::geom_hline(yintercept = R_max, linetype = "dashed",
                         colour = "#888780", linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = K_half, linetype = "dotted",
                         colour = "#888780", linewidth = 0.7) +
    ggplot2::annotate("text", x = K_half, y = R_max * 0.05,
                       label = paste0("K_half=", format(K_half, big.mark=",")),
                       hjust = -0.1, size = 3, colour = "#888780") +
    ggplot2::annotate("text", x = max(A_df$n) * 0.02, y = R_max * 1.02,
                       label = paste0("R_max=", R_max),
                       hjust = 0, vjust = 0, size = 3, colour = "#888780") +
    { if (ref_known)
        ggplot2::geom_vline(xintercept = n_flowering_ref,
                             colour = "#e34948", linewidth = 0.8)
    } +
    { if (ref_known)
        ggplot2::annotate("text", x = n_flowering_ref, y = R_max * 0.5,
                           label = paste0(" equil.\n n_fl=", n_flowering_ref),
                           hjust = 0, size = 2.8, colour = "#e34948")
    } +
    ggplot2::scale_y_continuous(limits = c(0, R_max * 1.08),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "Flowering adults (n_flowering)",
                  y = "Expected recruits / year",
                  title = "Beverton-Holt recruitment curve",
                  subtitle = "Shaded band: ±1 Poisson SD") +
    ggplot2::theme_minimal(base_size = 11)

  # -- Panel B: BH × shade --------------------------------------------------
  pB <- gg(B_df, aes(x = n, y = expected, colour = N_canopy)) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::scale_colour_manual(values = nc_pal, name = NULL) +
    { if (ref_known)
        ggplot2::annotate("point", x = n_flowering_ref,
                           y = bevholt(n_flowering_ref, params) *
                               (1 - shade_suppression(N_canopy_ref, params)),
                           colour = "#e34948", size = 3, shape = 16)
    } +
    { if (ref_known)
        ggplot2::annotate("text", x = n_flowering_ref,
                           y = bevholt(n_flowering_ref, params) *
                               (1 - shade_suppression(N_canopy_ref, params)),
                           label = " operating\n point",
                           hjust = 0, vjust = 1, size = 2.8, colour = "#e34948")
    } +
    ggplot2::scale_y_continuous(limits = c(0, R_max * 1.08),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "Flowering adults (n_flowering)",
                  y = "Expected recruits / year",
                  title = "Recruitment with canopy shade suppression") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 9))

  # -- Panel C: shade suppression vs N_canopy ------------------------------
  shade_note <- if (!shade_on) " [shade off: max_effect=0]" else ""
  pC <- gg(shade_df, aes(x = N_canopy, y = suppressed)) +
    ggplot2::geom_line(colour = "#008300", linewidth = 1) +
    ggplot2::geom_hline(yintercept = params$shade_dose_response$max_effect * 0.5,
                         linetype = "dotted", colour = "#888780") +
    ggplot2::geom_vline(xintercept = hs, linetype = "dashed",
                         colour = "#888780", linewidth = 0.7) +
    ggplot2::annotate("text", x = hs, y = 0.02,
                       label = paste0("half_sat=", format(hs, big.mark=",")),
                       hjust = -0.05, size = 3, colour = "#888780") +
    { if (ref_known)
        ggplot2::geom_vline(xintercept = N_canopy_ref,
                             colour = "#e34948", linewidth = 0.8)
    } +
    { if (ref_known)
        ggplot2::annotate("text", x = N_canopy_ref,
                           y = params$shade_dose_response$max_effect * 0.7,
                           label = paste0(" equil.\n N_canopy=",
                                           round(N_canopy_ref)),
                           hjust = 0, size = 2.8, colour = "#e34948")
    } +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = c(0, .03))) +
    ggplot2::labs(x = "Canopy density (N_canopy)",
                  y = "Fraction of recruitment suppressed",
                  title = paste0("Shade suppression", shade_note),
                  subtitle = paste0("max_effect=", params$shade_dose_response$max_effect,
                                     "  half_sat=", format(hs, big.mark=","))) +
    ggplot2::theme_minimal(base_size = 11)

  # -- Panel D: rust fecundity multiplier ----------------------------------
  rust_note <- if (rust_off) " [off: max_effect=0]" else
    paste0("  max_effect=", params$rust_recruit_dose_response$max_effect)
  pD <- gg(D_df, aes(x = eff_susc, y = weight)) +
    ggplot2::geom_line(colour = "#e34948", linewidth = 1,
                        linetype = if (rust_off) "dashed" else "solid") +
    ggplot2::scale_y_continuous(limits = c(0, 1.05),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "Effective susceptibility",
                  y = "Per-adult fecundity weight",
                  title = paste0("Rust fecundity suppression", rust_note),
                  subtitle = "1 - dose_response(eff_susc, rust_recruit_dose_response)") +
    ggplot2::theme_minimal(base_size = 11)

  # -- assemble -------------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, pC, pD, ncol = 2))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, pC, pD, ncol = 2)
  } else {
    for (p in list(pA, pB, pC, pD)) {
      print(p); readline("Enter for next panel ...")
    }
  }
}

# =============================================================================
# base R fallback
# =============================================================================
.base_recruitment <- function(A_df, B_df, shade_df, D_df, params,
                               n_flowering_ref, N_canopy_ref, ref_known,
                               shade_on, hs, rust_off, nc_max) {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  R_max  <- params$R_max
  K_half <- params$K_half
  nc_pal <- grDevices::colorRampPalette(c("#d9f0d3","#1b7837"))(
              length(unique(B_df$N_canopy)))

  # -- Panel A --------------------------------------------------------------
  plot(A_df$n, A_df$expected, type = "n", ylim = c(0, R_max * 1.08),
       xlab = "Flowering adults", ylab = "Expected recruits / year",
       main = "Beverton-Holt curve", las = 1)
  polygon(c(A_df$n, rev(A_df$n)), c(A_df$lo, rev(A_df$hi)),
          col = grDevices::adjustcolor("#2a78d6", 0.15), border = NA)
  lines(A_df$n, A_df$expected, col = "#2a78d6", lwd = 2)
  abline(h = R_max, lty = 2, col = "#888780")
  abline(v = K_half, lty = 3, col = "#888780")
  mtext(paste0("R_max=", R_max), side = 3, line = -1.2, adj = 0.98,
        cex = 0.8, col = "#888780")
  if (ref_known) abline(v = n_flowering_ref, col = "#e34948", lwd = 1.5)

  # -- Panel B --------------------------------------------------------------
  nc_labels <- levels(B_df$N_canopy)
  plot(NA, xlim = range(B_df$n), ylim = c(0, R_max * 1.08),
       xlab = "Flowering adults", ylab = "Expected recruits / year",
       main = "Recruitment × shade suppression", las = 1)
  for (i in seq_along(nc_labels)) {
    sub <- B_df[B_df$N_canopy == nc_labels[i], ]
    lines(sub$n, sub$expected, col = nc_pal[i], lwd = 2)
  }
  if (ref_known) {
    op_y <- bevholt(n_flowering_ref, params) *
            (1 - shade_suppression(N_canopy_ref, params))
    points(n_flowering_ref, op_y, pch = 16, col = "#e34948", cex = 1.5)
  }
  legend("topleft", bty = "n", cex = 0.75, legend = nc_labels,
         col = nc_pal, lwd = 2)

  # -- Panel C --------------------------------------------------------------
  shade_note <- if (!shade_on) " [off]" else ""
  plot(shade_df$N_canopy, shade_df$suppressed, type = "l", col = "#008300",
       lwd = 2, ylim = c(0, 1),
       xlab = "Canopy density (N_canopy)", ylab = "Fraction suppressed",
       main = paste0("Shade suppression", shade_note), las = 1)
  abline(v = hs, lty = 2, col = "#888780")
  mtext(paste0("half_sat=", format(hs, big.mark=",")),
        side = 3, line = -1.2, adj = 0.98, cex = 0.8, col = "#888780")
  if (ref_known) abline(v = N_canopy_ref, col = "#e34948", lwd = 1.5)

  # -- Panel D --------------------------------------------------------------
  rust_note <- if (rust_off) " [off]" else ""
  plot(D_df$eff_susc, D_df$weight, type = "l",
       col = "#e34948", lwd = 2,
       lty = if (rust_off) 2 else 1,
       ylim = c(0, 1.05), xlim = c(0, 1),
       xlab = "Effective susceptibility",
       ylab = "Per-adult fecundity weight",
       main = paste0("Rust fecundity suppression", rust_note), las = 1)
}
