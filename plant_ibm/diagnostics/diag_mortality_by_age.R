# diagnostics/diag_mortality_by_age.R
#
# Diagnostic: annual mortality probability by age (0:max_age), for a
# range of population sizes and rust intensities.
#
# CALLS THE ACTUAL MODEL FUNCTIONS from R/mortality.R -- this is
# simultaneously a calibration aid and a code check. Any change to
# weibull_hazard(), senescence_hazard(), juv_decline_hazard(), etc.
# is immediately reflected here without touching this script.
#
# USAGE (from project root):
#   source("diagnostics/diag_mortality_by_age.R")
#   diag_mortality_by_age()
#
# Or override params for a specific scenario:
#   p <- get_default_params()
#   p$juv_decline_dose_response$max_effect <- 0.5
#   diag_mortality_by_age(params = p)
#
# Save to file:
#   diag_mortality_by_age(save_path = "outputs/mortality_check.png")
#
# ACHIEVING A BATHTUB SHAPE (high at 0-1, plateau, rising again):
#   p$weibull_k                          <- 1      # flat background, not rising
#   p$juv_decline_dose_response          <- list(  # high juvenile mortality under canopy
#       max_effect = 0.5, form = "saturating",
#       half_sat = 800, hill = 1)
#   p$juv_decline_age_half_sat           <- 2      # halved by age 2
#   p$juv_decline_age_hill               <- 2      # fairly sharp
#   p$senescence_dose_response           <- list(  # rising late mortality
#       max_effect = 0.25, form = "sigmoid",
#       half_sat = 35, hill = 3)
#   diag_mortality_by_age(params = p)

# --- locate project root (works whether sourced from root or diagnostics/) ---
.root <- if (file.exists("params.R")) {
  "."
} else if (file.exists("../params.R")) {
  ".."
} else {
  stop("Source this script from the project root or diagnostics/ directory.")
}
source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))  # weibull_hazard, senescence_hazard,
                                            # juv_decline_hazard, hill_weight,
                                            # dose_response (all model functions)

# --- internal: compute all mortality components for one (age, N_canopy) pair -
# Mirrors nominate_deaths() logic, calling the same model functions with
# scalar inputs. N_canopy is the canopy density driving juv_decline_hazard;
# it is the only way population size feeds into mortality now that the
# separate density-dependent mortality term has been removed.
#
# resist_score = 0, env_susc = 1: fully susceptible plant, worst case.
.compute_components <- function(age, N_canopy, annual_env_t, params,
                                 resist_score = 0, env_susc = 1) {

  # Background Weibull hazard
  p_age <- weibull_hazard(age, params)

  # Senescence hazard: rising Hill function (off by default, max_effect=0)
  p_sen <- senescence_hazard(age, params)

  # Juvenile-decline hazard: only for age < AFR (non-resprouting).
  # N_canopy is the canopy density from canopy_density(pop); this is
  # how population density enters mortality -- via shade on juveniles,
  # not a separate density-dependent hazard term.
  p_jd <- if (age < params$age_first_flower_mean) {
    juv_decline_hazard(age, N_canopy, params)
  } else {
    0
  }

  # Rust hazard: disease triangle for a non-resprouting individual.
  p_rust <- 0
  if (!is.infinite(params$rust_start_year) && annual_env_t > 0) {
    eff_susc   <- (1 - resist_score) * env_susc * annual_env_t
    rust_decay <- hill_weight(age,
                              params$rust_dose_response$age_half_sat,
                              params$rust_dose_response$age_hill)
    peak       <- params$rust_pressure * params$rust_dose_response$age_peak
    rust_h     <- params$rust_dose_response$age_floor +
                  (peak - params$rust_dose_response$age_floor) * rust_decay
    p_rust     <- eff_susc * rust_h
  }

  # Probability union -- mirrors nominate_deaths() exactly
  p_total <- 1 - (1 - p_age) * (1 - p_sen) * (1 - p_jd) * (1 - p_rust)

  c(p_age = p_age, p_sen = p_sen, p_jd = p_jd, p_rust = p_rust, p_total = p_total)
}

# ============================================================================
# diag_mortality_by_age()
#
# Produces a 2-panel figure:
#   Left  -- total mortality by age for N_canopy_values, with and without
#             rust. Lines differ only in the juv_decline component (canopy-
#             dependent height), making the canopy-shading effect on juvenile
#             mortality directly visible. After the juvenile period, all lines
#             converge (senescence and rust are canopy-independent).
#   Right -- component breakdown at a single reference canopy density.
# ============================================================================
diag_mortality_by_age <- function(
    params           = get_default_params(),
    N_canopy_values  = c(200, 800, 2000),  # canopy densities to compare
    rust_env_t       = 0.7,
    N_canopy_ref     = NULL,  # canopy for component panel (NULL: middle value)
    ages             = 0:60,
    save_path        = NULL,
    width_in         = 12,
    height_in        = 5,
    res              = 150
) {

  if (is.null(N_canopy_ref)) N_canopy_ref <- N_canopy_values[ceiling(length(N_canopy_values) / 2)]

  # ---- build data frame: total mortality for all N_canopy values -----------
  total_df <- do.call(rbind, lapply(seq_along(N_canopy_values), function(i) {
    nc <- N_canopy_values[i]
    no_rust   <- sapply(ages, function(a) .compute_components(a, nc, 0,          params)["p_total"])
    with_rust <- sapply(ages, function(a) .compute_components(a, nc, rust_env_t, params)["p_total"])
    rbind(
      data.frame(age = ages, p = no_rust,   rust = "no rust",
                 N_canopy = paste0("N_canopy = ", format(nc, big.mark=",")),
                 nc_val = nc, stringsAsFactors = FALSE),
      data.frame(age = ages, p = with_rust, rust = paste0("rust (env_t=", rust_env_t, ")"),
                 N_canopy = paste0("N_canopy = ", format(nc, big.mark=",")),
                 nc_val = nc, stringsAsFactors = FALSE)
    )
  }))

  # ---- build data frame: component breakdown at N_canopy_ref ---------------
  comp_raw <- lapply(ages, function(a) .compute_components(a, N_canopy_ref, rust_env_t, params))
  comp_df <- data.frame(
    age        = ages,
    Weibull    = sapply(comp_raw, `[[`, "p_age"),
    Senescence = sapply(comp_raw, `[[`, "p_sen"),
    Juv_decline = sapply(comp_raw, `[[`, "p_jd"),
    Rust       = sapply(comp_raw, `[[`, "p_rust"),
    Total      = sapply(comp_raw, `[[`, "p_total")
  )

  sen_on  <- params$senescence_dose_response$max_effect > 0
  jd_on   <- params$juv_decline_dose_response$max_effect > 0
  rust_on <- !is.infinite(params$rust_start_year) && rust_env_t > 0

  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .plot_gg(total_df, comp_df, N_canopy_values, N_canopy_ref,
             rust_env_t, sen_on, jd_on, rust_on, params)
  } else {
    .plot_base(total_df, comp_df, N_canopy_values, N_canopy_ref,
               rust_env_t, sen_on, jd_on, rust_on, params)
  }

  invisible(list(total = total_df, components = comp_df))
}

# ============================================================================
# ggplot2 rendering
# ============================================================================
.plot_gg <- function(total_df, comp_df, N_canopy_values, N_canopy_ref, rust_env_t,
                      sen_on, jd_on, rust_on, params) {

  gg <- ggplot2::ggplot

  # colour ramp: one colour per N_canopy value
  pal <- c("#2a78d6", "#008300", "#e34948", "#eda100", "#4a3aa7")[seq_along(N_canopy_values)]
  names(pal) <- paste0("N_canopy = ", format(N_canopy_values, big.mark=","))

  # -- panel 1: total mortality by canopy density, with/without rust ---------
  p1 <- gg(total_df, ggplot2::aes(x = age, y = p, colour = N_canopy,
                                    linetype = rust)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = pal, name = "Canopy density") +
    ggplot2::scale_linetype_manual(
      values = c("no rust" = "solid",
                 setNames("dashed", paste0("rust (env_t=", rust_env_t, ")"))),
      name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult=c(0,.05))) +
    ggplot2::labs(x = "Age (years)", y = "Annual mortality probability",
                  title = "Total mortality by age and canopy density") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.spacing.y = ggplot2::unit(1, "pt"))

  # -- panel 2: component breakdown at N_canopy_ref -------------------------
  comp_long <- utils::stack(comp_df[, c("Weibull","Senescence","Juv_decline","Rust","Total")])
  comp_long$age <- rep(comp_df$age, 5)
  names(comp_long) <- c("hazard","component","age")
  comp_long$component <- as.character(comp_long$component)

  # label inactive terms
  comp_long$label <- comp_long$component
  comp_long$label[comp_long$component == "Senescence"  & !sen_on]  <- "Senescence (off)"
  comp_long$label[comp_long$component == "Juv_decline" & !jd_on]   <- "Juv decline (off)"
  comp_long$label[comp_long$component == "Rust"        & !rust_on] <- "Rust (off)"

  comp_pal <- c(
    "Weibull"           = "#888780",
    "Senescence"        = "#e87ba4",
    "Senescence (off)"  = "#e1e0d9",
    "Juv_decline"       = "#2a78d6",
    "Juv decline (off)" = "#b5d4f4",
    "Rust"              = "#e34948",
    "Rust (off)"        = "#f7c1c1",
    "Total"             = "#0b0b0b"
  )

  p2 <- gg(comp_long, ggplot2::aes(x = age, y = hazard, colour = label,
                                     linewidth = component == "Total")) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = comp_pal, name = "Component") +
    ggplot2::scale_linewidth_manual(values = c("FALSE" = 0.75, "TRUE" = 1.4),
                                     guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult=c(0,.05))) +
    ggplot2::labs(x = "Age (years)", y = "Hazard / probability",
                  title = paste0("Components at N_canopy = ",
                                  format(N_canopy_ref, big.mark=","),
                                 "  (rust env_t = ", rust_env_t, ")")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(18, "pt"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(p1, p2, ncol = 2))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  } else {
    print(p1)
    readline("press Enter for component panel ...")
    print(p2)
  }
}

# ============================================================================
# base R fallback
# ============================================================================
.plot_base <- function(total_df, comp_df, N_canopy_values, N_canopy_ref, rust_env_t,
                        sen_on, jd_on, rust_on, params) {

  pal <- c("#2a78d6","#008300","#e34948","#eda100","#4a3aa7")[seq_along(N_canopy_values)]
  ages <- comp_df$age

  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(old_par), add = TRUE)

  # -- panel 1 ---------------------------------------------------------------
  ymax <- max(total_df$p, na.rm = TRUE) * 1.05
  plot(NA, xlim = range(ages), ylim = c(0, ymax),
       xlab = "Age (years)", ylab = "Annual mortality probability",
       main = "Total mortality by age", las = 1)

  for (i in seq_along(N_canopy_values)) {
    nc_lab <- paste0("N_canopy = ", format(N_canopy_values[i], big.mark=","))
    sub_nr <- total_df[total_df$N_canopy == nc_lab & total_df$rust == "no rust", ]
    sub_r  <- total_df[total_df$N_canopy == nc_lab & grepl("^rust", total_df$rust), ]
    lines(sub_nr$age, sub_nr$p, col = pal[i], lwd = 2, lty = 1)
    lines(sub_r$age,  sub_r$p,  col = pal[i], lwd = 1.5, lty = 2)
  }

  legend("topleft", bty = "n", cex = 0.8,
         legend = c(paste0("N_canopy=", format(N_canopy_values, big.mark=","), " no rust"),
                    paste0("N_canopy=", format(N_canopy_values, big.mark=","), " rust")),
         col    = c(pal, pal),
         lty    = c(rep(1, length(N_canopy_values)), rep(2, length(N_canopy_values))),
         lwd    = c(rep(2, length(N_canopy_values)), rep(1.5, length(N_canopy_values))))

  # -- panel 2 ---------------------------------------------------------------
  comp_cols <- c(Weibull     = "#888780", Senescence  = "#e87ba4",
                 Juv_decline = "#2a78d6", Rust        = "#e34948",
                 Total       = "#0b0b0b")

  components <- c("Weibull","Senescence","Juv_decline","Rust","Total")
  ymax2 <- max(unlist(comp_df[, components]), na.rm = TRUE) * 1.05

  plot(NA, xlim = range(ages), ylim = c(0, ymax2),
       xlab = "Age (years)", ylab = "Hazard / probability",
       main = paste0("Components: N_canopy=", format(N_canopy_ref, big.mark=","),
                     "  rust env_t=", rust_env_t), las = 1)

  labels <- c("Weibull",
              if (!sen_on) "Senescence (off)" else "Senescence",
              if (!jd_on)  "Juv decline (off)" else "Juv decline",
              if (!rust_on)"Rust (off)" else paste0("Rust (env_t=", rust_env_t, ")"),
              "Total")
  lwds <- c(1.2, 1.2, 1.2, 1.2, 2)

  for (j in seq_along(components)) {
    lines(ages, comp_df[[components[j]]], col = comp_cols[components[j]], lwd = lwds[j])
  }
  legend("topleft", bty = "n", cex = 0.8, legend = labels,
         col = comp_cols, lwd = lwds)
}
