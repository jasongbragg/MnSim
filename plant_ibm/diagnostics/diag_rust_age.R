# diagnostics/diag_rust_age.R
#
# Diagnostic: rust hazard as a function of age, and its interaction with
# fire via the resprout bonus. Calls hill_weight(), dose_response(),
# weibull_hazard(), and senescence_hazard() from the model code directly.
# No simulation is required -- all curves are analytical.
#
# PANEL A -- Base rust age-hazard curve
#   rust_age_h(age) for several rust_pressure values. Shows the Hill-
#   function age-decay shape in isolation, before the disease triangle
#   multiplier is applied. Dashed lines mark age_floor (residual adult
#   hazard) and age_half_sat (where the curve is halfway between peak
#   and floor). Calls hill_weight() from R/mortality.R directly.
#
# PANEL B -- Rust probability by age × effective susceptibility
#   p_rust = eff_susc × rust_age_h at the params rust_pressure, for
#   several combined effective_susceptibility levels. A grey reference
#   line shows the background hazard (weibull + senescence) so the
#   rust component can be read against the non-rust mortality context.
#   Calls weibull_hazard() and senescence_hazard() from R/mortality.R.
#
# PANEL C (one per fire_age) -- Fire × rust: the resprout bonus in time
#   For a plant that burned at fire_age, shows two lines through the
#   resprout window and beyond at a reference effective_susceptibility:
#     "No fire (counterfactual)" -- p_rust = eff_susc × rust_age_h only,
#       following the plant's age as it would without fire. Shown as a
#       dashed declining curve continuing through all years.
#     "Resprouting / recovered" -- elevated by eff_susc × rust_resprout_h
#       during the window (years 0 to resprout_yrs_base - 1), then drops
#       back to the counterfactual line exactly at recovery.
#   The vertical gap between the two lines during the window is constant
#   (= eff_susc × rust_resprout_h) and is the quantity directly
#   calibrated by rust_dose_response$resprout$max_effect. The step-down
#   at recovery is its visual signature.
#   Calls dose_response() from R/mortality.R.
#
# NOTE: rust_start_year does not gate these plots -- all curves show
# "what the hazard would be if rust were active." If you want to check
# the actual runtime gating, look at nominate_deaths() directly.
#
# USAGE (from project root):
#   source("diagnostics/diag_rust_age.R")
#   diag_rust_age()
#   diag_rust_age(save_path = "outputs/rust_age.png")
#
#   # Explore effect of tuning the resprout bonus
#   p <- get_default_params()
#   p$rust_dose_response$resprout$max_effect <- 0.8
#   diag_rust_age(params = p)

if (file.exists("params.R")) {
  .root <- "."
} else if (file.exists("../params.R")) {
  .root <- ".."
} else {
  stop("Source from project root or diagnostics/ directory.")
}
source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))   # hill_weight(), dose_response(),
                                              # weibull_hazard(), senescence_hazard()

# --- Internal helpers: call model functions directly -------------------------

# Age-decaying rust hazard (before disease triangle multiplier)
# Calls hill_weight() from R/mortality.R
.rust_age_h <- function(age, rust_pressure, params) {
  peak    <- rust_pressure * params$rust_dose_response$age_peak
  floor_h <- params$rust_dose_response$age_floor
  floor_h + (peak - floor_h) * hill_weight(
    age,
    params$rust_dose_response$age_half_sat,
    params$rust_dose_response$age_hill
  )
}

# Resprout bonus hazard (flat across ages; calls dose_response from R/mortality.R)
.rust_resprout_h <- function(rust_pressure, params) {
  dose_response(rust_pressure, params$rust_dose_response$resprout)
}

# Background (non-rust) hazard contribution at a given age
# Calls weibull_hazard() and senescence_hazard() from R/mortality.R
.background_h <- function(age, params) {
  weibull_hazard(age, params) + senescence_hazard(age, params)
}

# =============================================================================
# diag_rust_age()
# =============================================================================
diag_rust_age <- function(
    params         = get_default_params(),
    ages           = 0:50,
    rust_pressures = c(0.5, 1.0, 2.0),       # Panel A: curves by pressure
    eff_susc_vals  = c(0.25, 0.50, 0.75, 1.00), # Panel B: curves by eff_susc
    fire_ages      = c(5, 25),                # Panel C: one sub-panel each
    eff_susc_ref   = 0.7,                     # Panel C: reference eff_susc
    window_extra   = 3L,                      # Panel C: years to show post-recovery
    save_path      = NULL,
    width_in       = 14,
    height_in      = 10,
    res            = 150
) {
  wlen           <- params$resprout_yrs_base
  rp             <- params$rust_pressure
  resprout_bonus <- .rust_resprout_h(rp, params)
  age_hs         <- params$rust_dose_response$age_half_sat
  age_fl         <- params$rust_dose_response$age_floor

  # --- Panel A: rust_age_h vs age ------------------------------------------
  A_df <- do.call(rbind, lapply(rust_pressures, function(p_val) {
    data.frame(
      age      = ages,
      hazard   = sapply(ages, .rust_age_h, rust_pressure = p_val, params = params),
      label    = paste0("rust_pressure = ", p_val),
      rp       = p_val,
      stringsAsFactors = FALSE
    )
  }))

  # --- Panel B: p_rust vs age, multiple eff_susc ---------------------------
  B_rust <- do.call(rbind, lapply(eff_susc_vals, function(es) {
    data.frame(
      age    = ages,
      hazard = es * sapply(ages, .rust_age_h, rust_pressure = rp, params = params),
      label  = paste0("eff_susc = ", es),
      es     = es,
      stringsAsFactors = FALSE
    )
  }))
  B_bg <- data.frame(
    age    = ages,
    hazard = sapply(ages, .background_h, params = params)
  )

  # --- Panel C: counterfactual vs resprouting/recovered --------------------
  C_df <- do.call(rbind, lapply(fire_ages, function(fa) {
    years <- seq(-1L, wlen + window_extra)
    do.call(rbind, lapply(years, function(yr) {
      age         <- max(fa + 1L + yr, 0L)
      rust_age    <- .rust_age_h(age, rp, params)
      is_resprout <- (yr >= 0L && yr < wlen)
      bonus       <- if (is_resprout) resprout_bonus else 0
      rbind(
        data.frame(fire_age = fa, year = yr,
                   p_rust   = eff_susc_ref * rust_age,
                   line     = "No fire (counterfactual)",
                   stringsAsFactors = FALSE),
        data.frame(fire_age = fa, year = yr,
                   p_rust   = eff_susc_ref * (rust_age + bonus),
                   line     = "Resprouting / recovered",
                   stringsAsFactors = FALSE)
      )
    }))
  }))

  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_rust_age(A_df, B_rust, B_bg, C_df, params, fire_ages, wlen,
                  eff_susc_ref, resprout_bonus, age_hs, age_fl, rp,
                  rust_pressures)
  } else {
    .base_rust_age(A_df, B_rust, B_bg, C_df, params, fire_ages, wlen,
                   eff_susc_ref, resprout_bonus, age_hs, age_fl, rp,
                   rust_pressures)
  }

  invisible(list(panel_A = A_df, panel_B = B_rust, bg = B_bg, panel_C = C_df))
}

# =============================================================================
# ggplot2 rendering
# =============================================================================
.gg_rust_age <- function(A_df, B_rust, B_bg, C_df, params, fire_ages, wlen,
                          eff_susc_ref, resprout_bonus, age_hs, age_fl, rp,
                          rust_pressures) {
  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes

  # Palettes
  n_rp <- length(rust_pressures)
  n_es <- length(unique(B_rust$es))
  rp_pal <- grDevices::colorRampPalette(c("#b5d4f4","#1a478f"))(n_rp)
  es_pal <- grDevices::colorRampPalette(c("#fdd0a2","#d94801"))(n_es)

  peak_ref <- params$rust_pressure * params$rust_dose_response$age_peak

  # -- Panel A: base rust hazard --------------------------------------------
  A_df$label <- factor(A_df$label, levels = unique(A_df$label))
  pA <- gg(A_df, aes(x = age, y = hazard, colour = label)) +
    ggplot2::geom_hline(yintercept = age_fl, linetype = "dashed",
                         colour = "#888780", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = age_hs, linetype = "dotted",
                         colour = "#888780", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = rp_pal, name = NULL) +
    ggplot2::annotate("text", x = age_hs, y = peak_ref * 1.05,
                       label = paste0("age_half_sat=", age_hs),
                       hjust = -0.05, size = 2.8, colour = "#888780") +
    ggplot2::annotate("text", x = max(A_df$age) * 0.98, y = age_fl * 1.15,
                       label = paste0("age_floor=", age_fl),
                       hjust = 1, size = 2.8, colour = "#888780") +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .07))) +
    ggplot2::labs(x = "Age (years)", y = "rust_age_h (hazard)",
                  title = "Rust age-hazard (before disease triangle)",
                  subtitle = paste0("age_peak=", params$rust_dose_response$age_peak,
                                     "  age_floor=", age_fl,
                                     "  age_half_sat=", age_hs,
                                     "  age_hill=", params$rust_dose_response$age_hill)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # -- Panel B: p_rust by eff_susc, with background -------------------------
  B_rust$label <- factor(B_rust$label, levels = unique(B_rust$label))
  pB <- gg(B_rust, aes(x = age, y = hazard, colour = label)) +
    ggplot2::geom_line(data = B_bg, aes(x = age, y = hazard),
                        colour = "#cccccc", linewidth = 1.2,
                        linetype = "solid", inherit.aes = FALSE) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = es_pal, name = NULL) +
    ggplot2::annotate("text", x = 1, y = B_bg$hazard[2] * 1.05,
                       label = "background\n(weibull+senescence)",
                       hjust = 0, size = 2.5, colour = "#aaaaaa") +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .07))) +
    ggplot2::labs(x = "Age (years)", y = "p_rust = eff_susc × rust_age_h",
                  title = paste0("Rust probability by age  (rust_pressure=", rp, ")"),
                  subtitle = "Grey = background hazard reference") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # -- Panel C: one sub-panel per fire_age ----------------------------------
  C_panels <- lapply(fire_ages, function(fa) {
    sub  <- C_df[C_df$fire_age == fa, ]
    sub$line <- factor(sub$line, levels = c("No fire (counterfactual)",
                                             "Resprouting / recovered"))
    ymax <- max(sub$p_rust, na.rm = TRUE) * 1.08

    strip_df <- data.frame(xmin = 0, xmax = wlen, ymin = 0, ymax = ymax)

    gg(sub, aes(x = year, y = p_rust, colour = line, linetype = line)) +
      ggplot2::geom_rect(data = strip_df,
                          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                          fill = "#d0e4f4", alpha = 0.4, colour = NA,
                          inherit.aes = FALSE) +
      ggplot2::geom_vline(xintercept = 0,    linetype = "dotted",
                           colour = "#888", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = wlen, linetype = "dashed",
                           colour = "#333", linewidth = 0.6) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_colour_manual(
        values = c("No fire (counterfactual)" = "#888780",
                   "Resprouting / recovered"  = "#e34948"),
        name = NULL) +
      ggplot2::scale_linetype_manual(
        values = c("No fire (counterfactual)" = "dashed",
                   "Resprouting / recovered"  = "solid"),
        name = NULL) +
      ggplot2::annotate("text", x = wlen + 0.1, y = ymax * 0.5,
                         label = paste0("recovery\n(yr ", wlen, ")"),
                         hjust = 0, size = 2.6, colour = "#333") +
      ggplot2::annotate("text", x = 0.2, y = ymax * 0.97,
                         label = paste0("gap = ", round(eff_susc_ref * resprout_bonus, 3)),
                         hjust = 0, size = 2.6, colour = "#e34948") +
      ggplot2::scale_y_continuous(limits = c(0, ymax),
                                    expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(x = "Year since fire  (shaded = resprout window)",
                    y = "p_rust",
                    title = paste0("Fire × rust: burned at age ", fa),
                    subtitle = paste0("eff_susc=", eff_susc_ref,
                                       "  resprout bonus=",
                                       round(resprout_bonus, 3))) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
  })

  # -- Assemble ------------------------------------------------------------
  all_panels <- c(list(pA, pB), C_panels)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(all_panels, ncol = 2))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    do.call(gridExtra::grid.arrange, c(all_panels, ncol = 2))
  } else {
    for (p in all_panels) { print(p); readline("Enter for next ...") }
  }
}

# =============================================================================
# base R fallback
# =============================================================================
.base_rust_age <- function(A_df, B_rust, B_bg, C_df, params, fire_ages, wlen,
                             eff_susc_ref, resprout_bonus, age_hs, age_fl, rp,
                             rust_pressures) {

  n_panels <- 2L + length(fire_ages)
  op <- par(mfrow = c(2L, ceiling(n_panels / 2L)),
             mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  n_rp <- length(rust_pressures)
  n_es <- length(unique(B_rust$es))
  rp_pal <- grDevices::colorRampPalette(c("#b5d4f4","#1a478f"))(n_rp)
  es_pal <- grDevices::colorRampPalette(c("#fdd0a2","#d94801"))(n_es)

  # Panel A
  rp_labels <- unique(A_df$label)
  ymax_A <- max(A_df$hazard, na.rm = TRUE) * 1.08
  plot(NA, xlim = range(A_df$age), ylim = c(0, ymax_A),
       xlab = "Age (years)", ylab = "rust_age_h",
       main = "Rust age-hazard (before disease triangle)", las = 1)
  abline(h = age_fl, lty = 2, col = "#888780")
  abline(v = age_hs, lty = 3, col = "#888780")
  for (j in seq_along(rp_labels)) {
    sub <- A_df[A_df$label == rp_labels[j], ]
    lines(sub$age, sub$hazard, col = rp_pal[j], lwd = 2)
  }
  legend("topright", bty = "n", cex = 0.75, legend = rp_labels,
         col = rp_pal, lwd = 2)

  # Panel B
  es_labels <- unique(B_rust$label)
  ymax_B <- max(c(B_rust$hazard, B_bg$hazard), na.rm = TRUE) * 1.08
  plot(NA, xlim = range(B_rust$age), ylim = c(0, ymax_B),
       xlab = "Age (years)", ylab = "p_rust",
       main = paste0("Rust probability by age (rust_pressure=", rp, ")"), las = 1)
  lines(B_bg$age, B_bg$hazard, col = "#cccccc", lwd = 2)
  for (j in seq_along(es_labels)) {
    sub <- B_rust[B_rust$label == es_labels[j], ]
    lines(sub$age, sub$hazard, col = es_pal[j], lwd = 2)
  }
  legend("topright", bty = "n", cex = 0.75,
         legend = c("background", es_labels),
         col = c("#cccccc", es_pal), lwd = 2)

  # Panel C (one per fire_age)
  for (fa in fire_ages) {
    sub  <- C_df[C_df$fire_age == fa, ]
    xr   <- range(sub$year)
    ymax <- max(sub$p_rust, na.rm = TRUE) * 1.08

    plot(NA, xlim = xr, ylim = c(0, ymax),
         xlab = "Year since fire", ylab = "p_rust",
         main = paste0("Fire × rust: burned at age ", fa), las = 1)
    rect(0, 0, wlen, ymax,
         col = grDevices::adjustcolor("#2a78d6", 0.08), border = NA)
    abline(v = 0,    lty = 3, col = "#888")
    abline(v = wlen, lty = 2, col = "#333")

    cf  <- sub[sub$line == "No fire (counterfactual)", ]
    rsp <- sub[sub$line == "Resprouting / recovered", ]
    cf  <- cf[order(cf$year), ]
    rsp <- rsp[order(rsp$year), ]
    lines(cf$year,  cf$p_rust,  col = "#888780", lty = 2, lwd = 2)
    lines(rsp$year, rsp$p_rust, col = "#e34948", lty = 1, lwd = 2)

    legend("topright", bty = "n", cex = 0.75,
           legend = c("No fire (counterfactual)", "Resprouting / recovered"),
           col = c("#888780","#e34948"), lty = c(2,1), lwd = 2)
  }
}
