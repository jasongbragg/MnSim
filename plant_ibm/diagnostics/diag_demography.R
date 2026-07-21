# diagnostics/diag_demography.R
#
# Demographic diagnostic: population trajectory and age structure
# contrasting fire alone vs fire + myrtle rust, with rust status
# consistent throughout (spin-up AND post-fire). The two trajectories
# diverge from the very start of the pre-fire window, so you can see:
#
#   pre-fire gap:   effect of rust on a standing, unburned population
#   post-fire gap:  combined rust + fire vs fire alone
#   within-scenario drop at year 0: fire effect independent of rust
#
# FOUR PHASES:
#   1a. SPIN-UP, rust OFF: 500 yr, fire off, rust off → rust-free equilibrium
#   1b. SPIN-UP, rust ON:  500 yr, fire off, rust on (max constant pressure,
#       annual_env_t ≈ 1.0) → rust-impacted equilibrium
#   2a. POST-FIRE, rust OFF: single fire at year 0, rust stays off, 100 yr
#   2b. POST-FIRE, rust ON:  single fire at year 0, rust stays on, 100 yr
#
# PANEL A -- Population trajectory, year -prefire_years to +postfire_years
#   Two lines through the whole window (both diverge pre-fire):
#     Blue:  rust OFF -- fire effect only
#     Red:   rust ON  -- fire + max rust pressure
#   Bold: N_IUCN (= reproductive adults, n_flowering) -- the IUCN metric.
#   Thin: N_alive -- total population.
#   Vertical dashed line at year 0 (fire). Horizontal dashed line at
#   rust-OFF pre-fire N_IUCN (undisturbed, rust-free reference).
#
# PANEL B -- Age structure: four snapshots via plot_age_structure_compare()
#   Pre-fire equil rust-OFF | Pre-fire equil rust-ON |
#   Year 100 rust-OFF       | Year 100 rust-ON
#   Shared y-axis. Reads top-row: what does rust do to an undisturbed stand?
#   Reads bottom-row: what does fire do to each? Compare columns for rust.
#
# USAGE (from project root):
#   source("diagnostics/diag_demography.R")
#   diag_demography()
#   diag_demography(save_path = "outputs/demography.png")
#
#   # Faster test run
#   diag_demography(spinup_years = 150L, postfire_years = 60L)
#
#   # Add resistance buffering
#   p <- get_default_params()
#   p$resist_freq0 <- c(0.05)
#   diag_demography(params = p)

if (file.exists("params.R")) {
  .root <- "."
} else if (file.exists("../params.R")) {
  .root <- ".."
} else {
  stop("Source from project root or diagnostics/ directory.")
}
source(file.path(.root, "params.R"))
source(file.path(.root, "R/mortality.R"))
source(file.path(.root, "R/genetics.R"))
source(file.path(.root, "R/individuals.R"))
source(file.path(.root, "R/recruitment.R"))
source(file.path(.root, "R/fire.R"))
source(file.path(.root, "R/census.R"))
source(file.path(.root, "R/simulate.R"))
source(file.path(.root, "diagnostics/diag_age_structure.R"))

# =============================================================================
# diag_demography()
# =============================================================================
diag_demography <- function(
    params          = get_default_params(),
    spinup_years    = 500L,
    prefire_years   = 20L,
    postfire_years  = 100L,
    rust_env_alpha  = 1e4,   # Beta(1e4, 1): annual_env_t ≈ 1.0 every year
    rust_env_beta   = 1,     # maximum constant pathogen pressure
    save_path       = NULL,
    width_in        = 14,
    height_in       = 10,
    res             = 150
) {

  # shared helper: apply rust-ON settings to a params list
  .rust_on <- function(p) {
    p$rust_start_year   <- 1L
    p$annual_env_alpha  <- rust_env_alpha
    p$annual_env_beta   <- rust_env_beta
    p
  }

  # --- Phase 1a: spin-up, rust OFF -----------------------------------------
  message(sprintf("Phase 1a: spin-up rust-OFF (%d yr)...", spinup_years))
  p_spin_cfact                  <- params
  p_spin_cfact$fire_years       <- integer(0)
  p_spin_cfact$fire_prob_annual <- 0
  p_spin_cfact$rust_start_year  <- Inf
  p_spin_cfact$n_years          <- spinup_years

  res_spin_cfact <- run_simulation(p_spin_cfact, year0 = 1L, verbose = FALSE)
  init_cfact <- list(individuals = res_spin_cfact$individuals,
                     resist_gt   = res_spin_cfact$resist_gt)

  # --- Phase 1b: spin-up, rust ON ------------------------------------------
  message(sprintf("Phase 1b: spin-up rust-ON  (%d yr)...", spinup_years))
  p_spin_rust                  <- .rust_on(params)
  p_spin_rust$fire_years       <- integer(0)
  p_spin_rust$fire_prob_annual <- 0
  p_spin_rust$n_years          <- spinup_years

  res_spin_rust <- run_simulation(p_spin_rust, year0 = 1L, verbose = FALSE)
  init_rust <- list(individuals = res_spin_rust$individuals,
                    resist_gt   = res_spin_rust$resist_gt)

  # --- Phase 2a: post-fire, rust OFF (from cfact init) ---------------------
  message("Phase 2a: post-fire rust-OFF...")
  p_cfact                  <- params
  p_cfact$fire_years       <- c(1L)
  p_cfact$fire_prob_annual <- 0
  p_cfact$rust_start_year  <- Inf
  p_cfact$n_years          <- postfire_years + 1L  # +1 so rel_year reaches postfire_years

  res_cfact <- run_simulation(p_cfact, init_state = init_cfact,
                               year0 = 1L, verbose = FALSE)

  # --- Phase 2b: post-fire, rust ON (from rust init) -----------------------
  message("Phase 2b: post-fire rust-ON...")
  p_rust         <- .rust_on(params)
  p_rust$fire_years       <- c(1L)
  p_rust$fire_prob_annual <- 0
  p_rust$n_years          <- postfire_years + 1L

  res_rust <- run_simulation(p_rust, init_state = init_rust,
                              year0 = 1L, verbose = FALSE)

  # --- Build trajectory data frame -----------------------------------------
  keep <- c("rel_year", "N_alive", "N_IUCN", "scenario")

  # Pre-fire: last prefire_years of each spin-up census, labelled -N..-1
  .pre <- function(census, sc) {
    df           <- tail(census, prefire_years)
    df$rel_year  <- seq_len(prefire_years) - prefire_years - 1L  # -20..-1
    df$scenario  <- sc
    df[, keep]
  }

  # Post-fire: year 1 (fire year) → rel_year 0, year 2 → rel_year 1, ...
  .post <- function(census, sc) {
    df           <- census
    df$rel_year  <- df$year - 1L
    df$scenario  <- sc
    df[df$rel_year <= postfire_years, keep]
  }

  traj <- rbind(
    .pre(res_spin_cfact$census, "Rust OFF"),
    .pre(res_spin_rust$census,  "Rust ON"),
    .post(res_cfact$census, "Rust OFF"),
    .post(res_rust$census,  "Rust ON")
  )

  traj$scenario <- factor(traj$scenario, levels = c("Rust OFF", "Rust ON"))

  # Reference: rust-free pre-fire N_IUCN
  ref_NIUCN <- tail(res_spin_cfact$census$N_IUCN, 1L)

  # --- Render --------------------------------------------------------------
  if (!is.null(save_path)) {
    grDevices::png(save_path, width = width_in, height = height_in,
                   units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  yr <- paste0("Yr ", postfire_years)
  age_pops <- setNames(
    list(res_spin_cfact$individuals,
         res_spin_rust$individuals,
         res_cfact$individuals,
         res_rust$individuals),
    c("Pre-fire, Rust OFF",
      "Pre-fire, Rust ON",
      paste0(yr, ", Rust OFF"),
      paste0(yr, ", Rust ON"))
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .gg_demography(traj, age_pops, ref_NIUCN, prefire_years, postfire_years)
  } else {
    .base_demography(traj, age_pops, ref_NIUCN, prefire_years, postfire_years)
  }

  invisible(list(
    trajectory    = traj,
    age_snapshots = age_pops,
    ref_NIUCN     = ref_NIUCN
  ))
}

# =============================================================================
# ggplot2 rendering
# =============================================================================
.gg_demography <- function(traj, age_pops, ref_NIUCN,
                            prefire_years, postfire_years) {
  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes

  pal <- c("Rust OFF" = "#2a78d6", "Rust ON" = "#e34948")

  # -- Panel A: trajectory --------------------------------------------------
  pA <- gg(traj, aes(x = rel_year, colour = scenario)) +
    ggplot2::geom_line(aes(y = N_alive), linewidth = 0.5, alpha = 0.35) +
    ggplot2::geom_line(aes(y = N_IUCN),  linewidth = 1.1) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::geom_hline(yintercept = ref_NIUCN, linetype = "dashed",
                         colour = "#888780", linewidth = 0.6) +
    ggplot2::annotate("text",
                       x = -prefire_years + 0.5, y = ref_NIUCN,
                       label = paste0("Pre-fire N_IUCN (rust-free) = ", ref_NIUCN),
                       hjust = 0, vjust = -0.4, size = 2.8, colour = "#888780") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         colour = "#333", linewidth = 0.7) +
    ggplot2::annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.3,
                       label = "Fire", size = 3, colour = "#333") +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                  expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::scale_x_continuous(
      breaks = seq(-prefire_years, postfire_years, by = 10)) +
    ggplot2::labs(
      x = "Year relative to fire",
      y = "Count",
      title = "Demographic trajectory: fire with and without myrtle rust",
      subtitle = paste0(
        "Bold: N_IUCN (reproductive adults)  ",
        "Thin: N_alive  |  ",
        "Rust ON: annual_env_t \u2248 1 (max constant pressure), ",
        "active throughout spin-up and post-fire"
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())

  # -- Panel B: 4-snapshot age structure ------------------------------------
  pB <- plot_age_structure_compare(age_pops, shared_y = TRUE,
    title = paste0(
      "Age structure: pre-fire equilibrium (left pair) vs year ",
      postfire_years, " post-fire (right pair)"
    )
  )

  # -- assemble -------------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, ncol = 1, heights = c(1, 1))
  } else {
    print(pA)
    readline("Enter for age structure panel ... ")
    print(pB)
  }
}

# =============================================================================
# base R fallback
# =============================================================================
.base_demography <- function(traj, age_pops, ref_NIUCN,
                              prefire_years, postfire_years) {

  op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.7, 0))
  on.exit(par(op), add = TRUE)

  pal <- c("Rust OFF" = "#2a78d6", "Rust ON" = "#e34948")
  ymax <- max(traj$N_alive, na.rm = TRUE) * 1.05
  xr   <- range(traj$rel_year)

  plot(NA, xlim = xr, ylim = c(0, ymax), las = 1,
       xlab = "Year relative to fire", ylab = "Count",
       main = "Demographic trajectory: fire with and without myrtle rust")
  abline(h = ref_NIUCN, lty = 2, col = "#888780")
  abline(v = 0, lty = 2, col = "#333")
  mtext("Fire", at = 0, side = 3, cex = 0.8)

  for (sc in levels(traj$scenario)) {
    sub <- traj[traj$scenario == sc, ]
    sub <- sub[order(sub$rel_year), ]
    lines(sub$rel_year, sub$N_alive, col = pal[sc], lwd = 1, lty = 1)
    lines(sub$rel_year, sub$N_IUCN,  col = pal[sc], lwd = 2)
  }
  legend("topright", bty = "n", cex = 0.8,
         legend = names(pal), col = pal, lwd = 2)

  plot_age_structure_compare(age_pops, shared_y = TRUE,
    title = paste0("Age structure: pre-fire (left pair) vs yr ",
                   postfire_years, " (right pair)"))
}
