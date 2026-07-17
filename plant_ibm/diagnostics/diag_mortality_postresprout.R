# diagnostics/diag_mortality_postresprout.R
#
# Diagnostic: mortality of resprouting plants after fire, decomposed into
# fire-independent (background), rust age-based, and rust resprout-bonus
# components -- the last being what new literature helps calibrate.
#
# PLOT A (left panel)
#   Mortality in the first post-fire year, as a function of age at fire.
#   Lines: background (no fire), resprouting at each rust intensity.
#   The gap between env_t=0 and env_t>0 resprouting lines is driven by
#   rust_dose_response$resprout (the bonus) + residual age-based rust.
#
# PLOT B (one panel per entry in fire_ages)
#   Annual mortality year-by-year across the resprout window, for reference
#   ages at time of fire. Shows the mortality step-down when recovery occurs.
#   Observed data points can be overlaid for direct calibration.
#
# Calls model functions from R/mortality.R directly; the union formula and
# rust hazard assembly are reproduced inline (they live in nominate_deaths())
# with explicit comments marking where each mirrors the model source.
#
# NOTE on window duration: the model draws resprout_yrs_remain from
# rpois(1, resprout_yrs_base) -- this diagnostic uses resprout_yrs_base
# as the expected (mean) duration for a single representative trajectory.
#
# USAGE (from project root):
#   source("diagnostics/diag_mortality_postresprout.R")
#   diag_mortality_postresprout()
#   diag_mortality_postresprout(save_path = "outputs/resprout_mort.png")
#
#   # Overlay observed mortality data on Plot B panels
#   # obs_data: named or positional list, one data.frame per fire_age entry.
#   # Columns required: year_since_fire (integer), mortality_rate (numeric).
#   # Optional: mort_lo, mort_hi for error bars (e.g. 95% CI).
#   diag_mortality_postresprout(
#     fire_ages = c(5, 25),
#     obs_data  = list(
#       data.frame(year_since_fire = 0:3, mortality_rate = c(0.30,0.18,0.10,0.08)),
#       data.frame(year_since_fire = 0:3, mortality_rate = c(0.12,0.08,0.06,0.05))
#     )
#   )

# --- locate project root -----------------------------------------------------
if (file.exists("params.R")) {
  .root <- "."
} else if (file.exists("../params.R")) {
  .root <- ".."
} else {
  stop("Source from project root or diagnostics/ directory.")
}

source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))   # all model hazard functions live here

# --- internal: mortality components for one (age, resprout-state) point ------
#
# Reproduces the per-individual block inside nominate_deaths() using the
# same model functions, for scalar inputs. Comments mark which model
# function each line calls. The only piece not a single callable is the
# rust hazard assembly and the union formula -- both are inline in
# nominate_deaths() -- so they're reproduced here with matching logic.
#
# annual_env_t is used directly as the pathogen vertex, regardless of
# rust_start_year -- this always computes "what would happen at this
# rust level" for the purpose of the diagnostic.
.resprout_components <- function(age, N_total, annual_env_t, params,
                                  is_resprouting,
                                  resist_score = 0, env_susc = 1) {
  age <- max(age, 0L)

  # Background hazards -- direct model function calls
  p_age <- weibull_hazard(age, params)                # R/mortality.R
  p_sen <- senescence_hazard(age, params)             # R/mortality.R
  p_dd  <- dd_hazard(N_total, params) *               # R/mortality.R
           dd_age_weight(age, params)                  # R/mortality.R

  # juv_decline: gated by !resprout in nominate_deaths; resprouting
  # individuals do NOT receive this term regardless of age. N_canopy
  # approximated as 40% of N_total (typical adult fraction at equilibrium).
  p_jd <- if (!is_resprouting && age < params$age_first_flower_mean)
    juv_decline_hazard(age, N_total * 0.40, params)   # R/mortality.R
  else 0

  # Rust -- disease triangle: host x environment x pathogen
  eff_susc   <- (1 - resist_score) * env_susc * annual_env_t

  # Age-decaying rust component (applies to all individuals when rust active)
  rust_decay <- hill_weight(age,                      # R/mortality.R
                             params$rust_dose_response$age_half_sat,
                             params$rust_dose_response$age_hill)
  peak       <- params$rust_pressure * params$rust_dose_response$age_peak
  rust_age_h <- params$rust_dose_response$age_floor +
                (peak - params$rust_dose_response$age_floor) * rust_decay
  p_rust_age <- eff_susc * rust_age_h

  # Resprout-state bonus: extra vulnerability of actively resprouting tissue,
  # on top of whatever age-based susceptibility the plant already carries.
  # This is the term the new article data will most directly calibrate.
  # Calls dose_response() from R/mortality.R.
  resprout_h      <- if (is_resprouting)
    dose_response(params$rust_pressure, params$rust_dose_response$resprout)
  else 0
  p_rust_resprout <- eff_susc * resprout_h

  # Total rust hazard: components summed before union, mirroring
  # nominate_deaths()'s: p_rust <- eff_susc * (rust_age_h + resprout_h)
  p_rust <- p_rust_age + p_rust_resprout

  # Probability union -- mirrors nominate_deaths() exactly
  p_total      <- 1 - (1-p_age)*(1-p_sen)*(1-p_dd)*(1-p_jd)*(1-p_rust)
  p_background <- 1 - (1-p_age)*(1-p_sen)*(1-p_dd)*(1-p_jd)  # no rust

  list(p_age = p_age, p_sen = p_sen, p_dd = p_dd, p_jd = p_jd,
       p_rust_age = p_rust_age, p_rust_resprout = p_rust_resprout,
       p_background = p_background, p_total = p_total)
}

# =============================================================================
# diag_mortality_postresprout()
# =============================================================================
diag_mortality_postresprout <- function(
    params       = get_default_params(),
    N_total      = 2000,
    rust_env_t   = c(0, 0.5, 1.0),
    fire_ages    = c(4, 12),
    age_range_A  = 3:60,
    window_extra = 3,
    obs_data     = NULL,   # list of data.frames (one per fire_ages entry)
    save_path    = NULL,
    width_in     = 14,
    height_in    = 5,
    res          = 150
) {

  wlen <- params$resprout_yrs_base   # expected (Poisson mean) window length

  # --- Plot A: first post-fire year, by age at fire -------------------------
  A_df <- do.call(rbind, lapply(age_range_A, function(fa) {
    age <- fa + 1L  # age at mortality check: aging precedes mortality in loop

    bg <- .resprout_components(age, N_total, 0, params, is_resprouting = FALSE)
    rows <- list(data.frame(age_at_fire = fa, p = bg$p_background,
                             label = "Background (no fire/rust)",
                             env_t = NA_real_, stringsAsFactors = FALSE))

    for (et in rust_env_t) {
      cr  <- .resprout_components(age, N_total, et, params, is_resprouting = TRUE)
      lbl <- if (et == 0) "Resprout, no rust" else paste0("Resprout env_t=", et)
      rows[[length(rows) + 1]] <- data.frame(
        age_at_fire = fa, p = cr$p_total,
        label = lbl, env_t = et, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }))

  # --- Plot B: trajectory across resprout window, one df per fire_age ------
  B_list <- lapply(fire_ages, function(fa) {
    years <- seq(-1L, wlen + window_extra)
    do.call(rbind, lapply(rust_env_t, function(et) {
      do.call(rbind, lapply(years, function(yr) {
        age     <- fa + 1L + yr   # age at mortality check each year
        is_resp <- yr >= 0L && yr < wlen
        cr      <- .resprout_components(age, N_total, et, params, is_resp)
        data.frame(
          fire_age        = fa,
          year_since_fire = yr,
          env_t           = et,
          is_resprouting  = is_resp,
          p_total         = cr$p_total,
          p_background    = cr$p_background,
          p_rust_age      = cr$p_rust_age,
          p_rust_resprout = cr$p_rust_resprout,
          label           = if (et == 0) "No rust" else paste0("Rust env_t=", et),
          stringsAsFactors = FALSE
        )
      }))
    }))
  })
  names(B_list) <- paste0("age_", fire_ages)

  # --- flags for annotation -------------------------------------------------
  bonus_on <- params$rust_dose_response$resprout$max_effect > 0

  # --- render ---------------------------------------------------------------
  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_resprout(A_df, B_list, wlen, fire_ages, rust_env_t, obs_data, bonus_on)
  } else {
    .base_resprout(A_df, B_list, wlen, fire_ages, rust_env_t, obs_data, bonus_on)
  }

  invisible(list(plot_A = A_df, plot_B = B_list))
}

# =============================================================================
# ggplot2 rendering
# =============================================================================
.gg_resprout <- function(A_df, B_list, wlen, fire_ages, rust_env_t,
                          obs_data, bonus_on) {
  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes

  # colour palette: one colour per env_t value
  n_et   <- length(rust_env_t)
  et_pal <- grDevices::colorRampPalette(c("#2a78d6", "#e34948"))(n_et)
  names(et_pal) <- as.character(rust_env_t)

  bonus_note <- if (!bonus_on) "  [rust resprout bonus off in params]" else ""

  # --- Plot A ---------------------------------------------------------------
  A_bg   <- A_df[is.na(A_df$env_t), ]
  A_resp <- A_df[!is.na(A_df$env_t), ]
  A_resp$env_t_ch <- as.character(A_resp$env_t)

  pA <- gg(A_resp, aes(x = age_at_fire, y = p, colour = env_t_ch)) +
    ggplot2::geom_line(data = A_bg, aes(x = age_at_fire, y = p),
                        colour = "#888780", linetype = "dashed",
                        linewidth = 0.9, inherit.aes = FALSE) +
    ggplot2::annotate("text",
                       x = max(A_df$age_at_fire) * 0.97, y = tail(A_bg$p, 1),
                       label = "background\n(no fire)",
                       hjust = 1, vjust = -0.3, colour = "#888780", size = 2.8) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::scale_colour_manual(
      values = et_pal,
      labels = ifelse(rust_env_t == 0, "Resprout, no rust",
                       paste0("Resprout, env_t=", rust_env_t)),
      name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .07))) +
    ggplot2::labs(x = "Age at time of fire (years)",
                  y = "Annual mortality probability",
                  title = paste0("First post-fire year", bonus_note)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # --- Plot B ---------------------------------------------------------------
  B_all <- do.call(rbind, B_list)
  B_all$panel  <- paste0("Burned at age ", B_all$fire_age)
  B_all$env_ch <- as.character(B_all$env_t)

  # grey background strip for the resprout window
  strip_df <- data.frame(
    panel = paste0("Burned at age ", fire_ages),
    xmin  = 0, xmax = wlen
  )

  pB <- gg(B_all, aes(x = year_since_fire, y = p_total, colour = env_ch)) +
    ggplot2::geom_rect(data = strip_df,
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                        fill = "#d0e4f4", alpha = 0.4, colour = NA,
                        inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 0,    linetype = "dotted", colour = "#888") +
    ggplot2::geom_vline(xintercept = wlen, linetype = "dashed", colour = "#333",
                         linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::facet_wrap(~ panel, ncol = length(fire_ages), scales = "fixed") +
    ggplot2::scale_colour_manual(
      values = et_pal,
      labels = ifelse(rust_env_t == 0, "No rust",
                       paste0("Rust env_t=", rust_env_t)),
      name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .07))) +
    ggplot2::annotate("text", x = wlen + 0.15, y = Inf,
                       label = paste0("recovery\n(yr ", wlen, ")"),
                       hjust = 0, vjust = 1.3, size = 2.6, colour = "#333") +
    ggplot2::labs(x = "Year since fire  (0 = fire year,  shaded = resprout window)",
                  y = "Annual mortality probability",
                  title = paste0("Resprout window trajectory", bonus_note)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"))

  # optional observed data overlay
  if (!is.null(obs_data)) {
    for (i in seq_along(fire_ages)) {
      if (i > length(obs_data) || is.null(obs_data[[i]])) next
      od       <- obs_data[[i]]
      od$panel <- paste0("Burned at age ", fire_ages[i])
      pB <- pB +
        ggplot2::geom_point(data = od, inherit.aes = FALSE,
                             aes(x = year_since_fire, y = mortality_rate),
                             colour = "black", size = 2.5, shape = 16)
      if (all(c("mort_lo","mort_hi") %in% names(od))) {
        pB <- pB +
          ggplot2::geom_errorbar(data = od, inherit.aes = FALSE,
                                  aes(x = year_since_fire,
                                      ymin = mort_lo, ymax = mort_hi),
                                  width = 0.25, colour = "black")
      }
    }
  }

  # combine
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, widths = c(1.1, 1.9)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, ncol = 2, widths = c(1.1, 1.9))
  } else {
    print(pA); readline("Enter for Plot B ... "); print(pB)
  }
}

# =============================================================================
# base R fallback
# =============================================================================
.base_resprout <- function(A_df, B_list, wlen, fire_ages, rust_env_t,
                            obs_data, bonus_on) {
  n_panels <- 1L + length(fire_ages)
  op <- par(mfrow = c(1L, n_panels), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  n_et   <- length(rust_env_t)
  et_col <- grDevices::colorRampPalette(c("#2a78d6", "#e34948"))(n_et)
  bonus_note <- if (!bonus_on) "\n[resprout bonus off]" else ""

  # --- Plot A ---------------------------------------------------------------
  yA <- max(A_df$p, na.rm = TRUE) * 1.08
  plot(NA, xlim = range(A_df$age_at_fire), ylim = c(0, yA),
       xlab = "Age at fire (years)", ylab = "Annual mortality probability",
       main = paste0("Year 1 post-fire mortality", bonus_note), las = 1)

  bg <- A_df[is.na(A_df$env_t), ]
  lines(bg$age_at_fire, bg$p, col = "#888780", lty = 2, lwd = 1.5)

  for (j in seq_len(n_et)) {
    sub <- A_df[!is.na(A_df$env_t) & A_df$env_t == rust_env_t[j], ]
    lines(sub$age_at_fire, sub$p, col = et_col[j], lwd = 2)
  }

  legend("topleft", bty = "n", cex = 0.75,
         legend = c("Background (no fire)",
                    ifelse(rust_env_t == 0, "Resprout, no rust",
                            paste0("Resprout, env_t=", rust_env_t))),
         col = c("#888780", et_col),
         lty = c(2L, rep(1L, n_et)), lwd = c(1.5, rep(2, n_et)))

  # --- Plot B (one per fire_age) -------------------------------------------
  B_all <- do.call(rbind, B_list)
  yB    <- max(B_all$p_total, na.rm = TRUE) * 1.08

  for (k in seq_along(fire_ages)) {
    fa  <- fire_ages[k]
    sub <- B_all[B_all$fire_age == fa, ]
    xr  <- range(sub$year_since_fire)

    plot(NA, xlim = xr, ylim = c(0, yB),
         xlab = "Year since fire", ylab = "Annual mortality probability",
         main = paste0("Resprout: burned age ", fa, bonus_note), las = 1)

    rect(0, 0, wlen, yB,
         col = grDevices::adjustcolor("#2a78d6", 0.08), border = NA)
    abline(v = 0,    lty = 3, col = "#888")
    abline(v = wlen, lty = 2, col = "#333")
    mtext(paste0("recovery yr ", wlen), at = wlen, side = 3, cex = 0.65)

    for (j in seq_len(n_et)) {
      et_sub <- sub[sub$env_t == rust_env_t[j], ]
      et_sub <- et_sub[order(et_sub$year_since_fire), ]
      lines(et_sub$year_since_fire, et_sub$p_total, col = et_col[j], lwd = 2)
    }

    if (!is.null(obs_data) && k <= length(obs_data) && !is.null(obs_data[[k]])) {
      od <- obs_data[[k]]
      points(od$year_since_fire, od$mortality_rate, pch = 16, cex = 1.3)
      if (all(c("mort_lo","mort_hi") %in% names(od)))
        arrows(od$year_since_fire, od$mort_lo, od$year_since_fire, od$mort_hi,
               angle = 90, code = 3, length = 0.06)
    }

    legend("topright", bty = "n", cex = 0.75,
           legend = ifelse(rust_env_t == 0, "No rust", paste0("env_t=", rust_env_t)),
           col = et_col, lwd = 2)
  }
}
