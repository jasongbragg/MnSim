# talk/make_talk_figures.R
#
# Figures for a preliminary talk on the M. nodosa IUCN A4 assessment.
# Re-run after any model change to refresh every figure from the current
# model code and params.R.
#
#   cd plant_ibm
#   /opt/R/4.3.1/bin/Rscript talk/make_talk_figures.R            # everything
#   /opt/R/4.3.1/bin/Rscript talk/make_talk_figures.R pegg       # one stage
#   stages: pegg | range | figures | all (default)
#
# OUTPUTS (talk/figures/)
#   fig1_model_diagram.png      static schematic (from talk/model_diagram.svg)
#   fig2_pegg_fit.png           model vs Pegg et al. (2025) dead resprouters over time
#   fig3_trajectories.png       range-wide ensemble: N relative to pre-rust reference
#   fig4_outcomes.png           distribution of A4 declines and IUCN category shares
#   fig5_map.png                median decline by grid cell
#   talk_numbers.txt            the numbers quoted on the slides
#
# CACHING
#   Every simulation is cached as one .rds per run under RESULTS_DIR, keyed by
#   a hash of its full parameter set and scenario. Re-running reuses results
#   whose parameters are unchanged, so after a model change point RESULTS_DIR
#   at a new folder (or delete the old one) to force fresh runs.
#
# ENSEMBLE DESIGN (range-wide)
#   run i: own RNG stream (MASTER_SEED + i), so adding runs never changes
#   earlier ones. Each run draws
#     - a grid cell, with probability proportional to abundance (weight_n_raw);
#       cell-specific fire frequency/extent and rust pressure
#     - a simulation seed (fire timing, annual rust years, demography)
#     - values for the parameters in PARAM_DRAWS (model uncertainty)
#   and gives one A4 decline. The ensemble of declines is the result.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ============================================================================
# Configuration
# ============================================================================
RESULTS_DIR <- "results/talk_v1/"
FIG_DIR     <- "talk/figures/"
N_CORES     <- max(1L, parallel::detectCores() - 1L)
N_PEGG      <- 100L     # Pegg-scenario runs
N_RANGE     <- 300L     # range-wide A4 runs
MASTER_SEED <- 20261022L

# Model uncertainty: parameters drawn independently per run, uniform on
# [min, max]. 'path' is the location in the params list. Set to list() to
# switch parameter uncertainty off (site + seed variation only).
PARAM_DRAWS <- list(
  resprout_rust_mortality = list(
    path  = c("rust_dose_response", "resprout", "max_effect"),
    min   = 0.35, max = 0.60,
    label = "Extra rust mortality of resprouting plants (max_effect)",
    note  = "range consistent with Pegg et al. 2025 (24/50 dead, 95% CI 0.34-0.63)")
)

PEGG_CSV     <- "talk/data/pegg2025_yarringully_fig3a.csv"
DIAGRAM_SVG  <- "talk/model_diagram.svg"
QLD_LAT      <- -28.16   # approx. NSW/Qld border on the coast

# ============================================================================
# Model code
# ============================================================================
source("params.R")
for (f in c("individuals", "genetics", "mortality", "fire", "recruitment",
            "census", "simulate"))
  source(file.path("R", paste0(f, ".R")))
source("runner/params_utils.R")
source("runner/iucn_a4.R")
source("spatial/03_range_wide_iucn.R")   # build_site_params(), ASSESS, YEAR_*

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,     recursive = TRUE, showWarnings = FALSE)

TALK_VERSION <- "talk-v1"

# Apply this run's parameter draws (uses the current RNG stream).
draw_params <- function(p) {
  vals <- list()
  for (nm in names(PARAM_DRAWS)) {
    d <- PARAM_DRAWS[[nm]]
    v <- stats::runif(1, d$min, d$max)
    p[[d$path]] <- v
    vals[[nm]]  <- v
  }
  list(params = p, values = vals)
}

base_value <- function(d) get_default_params()[[d$path]]

# ============================================================================
# Stage 1: Pegg scenario -- cumulative death of resprouters through time
# ============================================================================
# Same scenario as runner/pegg_validation.R (spin-up 1910, rust 2010, one
# fire 2020 burning every tree), continued to the end of PEGG_LAST_YEAR.
# Cohort: plants alive at the end of 2019 that were not killed by the fire,
# i.e. the resprouters. cum_dead[y] = share of them dead (any non-fire cause)
# by the end of model year y, plotted at 1 Jan of y + 1.
PEGG_SPINUP <- 1910L; PEGG_RUST <- 2010L; PEGG_FIRE <- 2020L
PEGG_LAST_YEAR <- 2025L

pegg_trajectory <- function(i) {
  set.seed(MASTER_SEED + 100000L + i)
  dr <- draw_params(get_default_params())
  p  <- dr$params
  p$seed             <- 5000L + i
  p$rust_start_year  <- PEGG_RUST
  p$fire_years       <- PEGG_FIRE
  p$fire_prob_annual <- 0
  p$fire_p_fimp      <- 1
  p$n_years          <- PEGG_LAST_YEAR - PEGG_SPINUP + 1L

  id  <- make_run_id(list(params = p, version = TALK_VERSION, kind = "pegg"))
  out <- file.path(RESULTS_DIR, "pegg", paste0(id, ".rds"))
  if (file.exists(out)) return(readRDS(out))

  pop <- run_simulation(p, year0 = PEGG_SPINUP)$individuals
  cal_birth <- ifelse(pop$birth_year < PEGG_SPINUP,
                      PEGG_SPINUP + pop$birth_year, pop$birth_year)
  prefire <- cal_birth <= PEGG_FIRE - 1L &
             (pop$alive | (!is.na(pop$death_year) & pop$death_year >= PEGG_FIRE))
  resprouters <- prefire & !(pop$death_cause %in% "fire")
  n_resp <- sum(resprouters)
  yrs <- PEGG_FIRE:PEGG_LAST_YEAR
  cum <- vapply(yrs, function(y)
    sum(resprouters & pop$death_cause %in% "mortality" & pop$death_year <= y) / n_resp,
    numeric(1))
  res <- data.frame(run = i, model_year = yrs,
                    date = as.Date(sprintf("%d-01-01", yrs + 1L)),
                    cum_dead = cum, n_resprouters = n_resp,
                    frac_fire_killed = sum(prefire & pop$death_cause %in% "fire") / sum(prefire))
  for (nm in names(dr$values)) res[[nm]] <- dr$values[[nm]]
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, out)
  res
}

run_pegg_stage <- function() {
  message(sprintf("Pegg scenario: %d runs on %d cores", N_PEGG, N_CORES))
  res <- parallel::mclapply(seq_len(N_PEGG), function(i)
    tryCatch(pegg_trajectory(i), error = function(e) {
      message("  pegg run ", i, ": ", conditionMessage(e)); NULL }),
    mc.cores = N_CORES)
  df <- do.call(rbind, res)
  saveRDS(df, file.path(RESULTS_DIR, "pegg_trajectories.rds"))
  df
}

# ============================================================================
# Stage 2: range-wide A4 ensemble
# ============================================================================
load_sites <- function() {
  s <- read.csv(SITE_PARAMS_PATH)
  ok <- !is.na(s$annual_env_alpha) & !is.na(s$annual_env_beta) &
        !is.na(s$fire_prob_annual) & !is.na(s$fire_p_fimp)
  attr(s, "n_dropped")      <- sum(!ok)
  attr(s, "weight_dropped") <- sum(s$weight_n_raw[!ok]) / sum(s$weight_n_raw)
  s <- s[ok, ]
  s$region <- ifelse(s$lat_centre > QLD_LAT, "Qld", "NSW")
  s
}

plan_range_run <- function(i, sites, base) {
  set.seed(MASTER_SEED + i)
  k  <- sample.int(nrow(sites), 1L, prob = sites$weight_n_raw)
  sd <- sample.int(99999L, 1L)
  dr <- draw_params(base)
  p  <- build_site_params(dr$params, sites[k, ])
  p$seed <- sd
  row <- data.frame(run = i, cell_row = k, cell_id = sites$cell_id[k],
                    lon = sites$lon_centre[k], lat = sites$lat_centre[k],
                    region = sites$region[k], sim_seed = sd,
                    run_id = iucn_run_id(p, ASSESS))
  for (nm in names(dr$values)) row[[nm]] <- dr$values[[nm]]
  list(params = p, row = row)
}

run_range_stage <- function() {
  sites <- load_sites()
  base  <- get_default_params()
  base$fire_years <- integer(0)
  message(sprintf("Range-wide: %d runs on %d cores (%d cells; %d cells without rust/fire values dropped, %.1f%% of abundance weight)",
                  N_RANGE, N_CORES, nrow(sites), attr(sites, "n_dropped"),
                  100 * attr(sites, "weight_dropped")))
  plans <- lapply(seq_len(N_RANGE), plan_range_run, sites = sites, base = base)
  plan  <- do.call(rbind, lapply(plans, `[[`, "row"))
  write.csv(plan, file.path(RESULTS_DIR, "range_plan.csv"), row.names = FALSE)

  rdir <- file.path(RESULTS_DIR, "range")
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  parallel::mclapply(seq_len(N_RANGE), function(i) {
    out <- file.path(rdir, paste0(plans[[i]]$row$run_id, ".rds"))
    if (file.exists(out)) return(invisible(NULL))
    tryCatch(
      suppressWarnings(run_iucn_a4(
        plans[[i]]$params,
        year_spinup = YEAR_SPINUP, year_gentime_start = YEAR_GENTIME_START,
        year_gentime_end = YEAR_GENTIME_END, year_rust = YEAR_RUST,
        year_assess = YEAR_ASSESS, year_end = YEAR_END, avg_years = AVG_YEARS,
        save_rds = out, plot = FALSE, verbose = FALSE)),
      error = function(e) message("  range run ", i, ": ", conditionMessage(e)))
    invisible(NULL)
  }, mc.cores = N_CORES)
  invisible(plan)
}

# Gather summaries and N_IUCN / N_ref trajectories for the planned runs.
TRAJ_YEARS <- 1950:2130
collect_range <- function() {
  plan <- read.csv(file.path(RESULTS_DIR, "range_plan.csv"))
  plan <- plan[plan$run <= N_RANGE, ]
  rdir <- file.path(RESULTS_DIR, "range")
  rows <- list(); traj <- list()
  for (j in seq_len(nrow(plan))) {
    f <- file.path(rdir, paste0(plan$run_id[j], ".rds"))
    if (!file.exists(f)) next
    r <- readRDS(f)
    s <- summarise_iucn_a4(r)
    rows[[length(rows) + 1L]] <- cbind(plan[j, setdiff(names(plan), "run_id")], s)
    cen <- r$census
    N <- cen$N_IUCN[match(TRAJ_YEARS, cen$year)]
    if (isTRUE(cen$extinct[nrow(cen)])) N[TRAJ_YEARS > max(cen$year)] <- 0
    traj[[length(traj) + 1L]] <- data.frame(run = plan$run[j], year = TRAJ_YEARS,
                                            rel = N / s$N_ref)
  }
  if (length(rows) == 0L) stop("No range-wide results found in ", rdir)
  df <- do.call(rbind, rows)
  df$category <- factor(vapply(df$pct_decline, iucn_category, ""),
                        levels = c("LC/NT", "VU", "EN", "CR"))
  list(runs = df, traj = do.call(rbind, traj))
}

# ============================================================================
# Stage 3: figures
# ============================================================================
CAT_COL <- c("LC/NT" = "#8fb3cf", VU = "#eda100", EN = "#d6452f", CR = "#7a0f12")
INK     <- "#2b2b2b"
theme_talk <- function(base = 15) {
  theme_minimal(base_size = base) +
    theme(plot.title    = element_text(face = "bold", colour = INK),
          plot.subtitle = element_text(colour = "#555555"),
          plot.caption  = element_text(colour = "#777777", size = base * 0.62, hjust = 0),
          panel.grid.minor = element_blank(),
          axis.title    = element_text(colour = INK))
}
wrap <- function(x, width = 165) paste(strwrap(x, width), collapse = "\n")
save_fig <- function(p, name, w = 12, h = 6.75) {
  ggsave(file.path(FIG_DIR, name), p, width = w, height = h, dpi = 200, bg = "white")
  message("  wrote ", file.path(FIG_DIR, name))
}
draws_caption <- function() {
  if (length(PARAM_DRAWS) == 0L) return("No parameter draws (site and seed variation only).")
  paste(vapply(PARAM_DRAWS, function(d)
    sprintf("%s ~ U(%.2f, %.2f) [%s]", d$label, d$min, d$max, d$note), ""),
    collapse = "; ")
}

fig_diagram <- function() {
  if (!file.exists(DIAGRAM_SVG)) { message("  no ", DIAGRAM_SVG); return(invisible()) }
  rsvg::rsvg_png(DIAGRAM_SVG, file.path(FIG_DIR, "fig1_model_diagram.png"),
                 width = 2400)
  message("  wrote ", file.path(FIG_DIR, "fig1_model_diagram.png"))
}

fig_pegg <- function(traj) {
  obs <- read.csv(PEGG_CSV, comment.char = "#")
  obs$date <- as.Date(obs$date)
  obs$k  <- round(obs$pct_dead / 100 * obs$n_trees)
  ci <- t(mapply(function(k, n) stats::binom.test(k, n)$conf.int, obs$k, obs$n_trees))
  obs$lo <- ci[, 1]; obs$hi <- ci[, 2]; obs$p <- obs$k / obs$n_trees
  obs$kind <- ifelse(obs$exact, "Pegg et al. (text: 24 of 50)",
                     "Pegg et al. (read from Fig. 3a)")

  start <- data.frame(run = unique(traj$run), date = as.Date("2020-01-01"), cum_dead = 0)
  tr <- rbind(start, traj[, c("run", "date", "cum_dead")])
  q <- do.call(rbind, lapply(split(tr, tr$date), function(d) data.frame(
    date = d$date[1], med = median(d$cum_dead),
    q05 = quantile(d$cum_dead, .05), q25 = quantile(d$cum_dead, .25),
    q75 = quantile(d$cum_dead, .75), q95 = quantile(d$cum_dead, .95))))
  at22 <- traj$cum_dead[traj$model_year == 2021]

  ggplot() +
    annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-01-31"),
             ymin = -Inf, ymax = Inf, fill = "#f2efe8") +
    annotate("text", x = as.Date("2020-03-20"), y = 0.97, hjust = 0, vjust = 1,
             label = "field monitoring", colour = "#8a7f6a", size = 4.2) +
    geom_ribbon(data = q, aes(date, ymin = q05, ymax = q95), fill = "#2a78d6", alpha = 0.15) +
    geom_ribbon(data = q, aes(date, ymin = q25, ymax = q75), fill = "#2a78d6", alpha = 0.28) +
    geom_line(data = q, aes(date, med), colour = "#1d5fae", linewidth = 1.1) +
    geom_point(data = q, aes(date, med), colour = "#1d5fae", size = 2.4) +
    geom_errorbar(data = obs, aes(date, ymin = lo, ymax = hi), width = 18,
                  colour = "#b0b0b0", linewidth = 0.5) +
    geom_point(data = obs, aes(date, p, shape = kind, fill = kind), size = 3.2,
               colour = INK, stroke = 0.6) +
    scale_shape_manual(values = c(21, 23), name = NULL) +
    scale_fill_manual(values = c("white", "#d6452f"), name = NULL) +
    geom_vline(xintercept = as.Date("2020-01-01"), colour = "#e08000",
               linewidth = 0.9, linetype = "dashed") +
    annotate("text", x = as.Date("2020-01-20"), y = 0.86, hjust = 0,
             label = "fire", colour = "#e08000", size = 4.5) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                 limits = as.Date(c("2019-10-01", "2026-03-01"))) +
    labs(x = NULL, y = "Resprouting trees dead",
         title = "Post-fire death of resprouting M. nodosa under myrtle rust",
         subtitle = sprintf(
           "Model: median, 50%% and 90%% of %d runs. Jan 2022: model %.0f%%, observed 48%% (95%% CI 34-63%%)",
           length(unique(traj$run)), 100 * median(at22)),
         caption = wrap(paste0(
           "Field data: Pegg et al. 2025, Austral Ecology, Yarringully NR, 50 trees resprouting when monitoring began; error bars 95% binomial CI. ",
           "Monthly values read off the published figure (approximate).\n",
           "Model: annual time step (deaths shown at the end of each model year); every tree burned (2020); rust present since 2010. ",
           draws_caption()))) +
    theme_talk() +
    theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = NA))
}

fig_trajectories <- function(R) {
  tr <- R$traj; runs <- R$runs
  q <- do.call(rbind, lapply(split(tr, tr$year), function(d) data.frame(
    year = d$year[1], med = median(d$rel, na.rm = TRUE),
    q05 = quantile(d$rel, .05, na.rm = TRUE), q25 = quantile(d$rel, .25, na.rm = TRUE),
    q75 = quantile(d$rel, .75, na.rm = TRUE), q95 = quantile(d$rel, .95, na.rm = TRUE))))
  set.seed(1L)
  show <- sample(unique(tr$run), min(25L, length(unique(tr$run))))
  ymax <- 1.8
  ggplot() +
    annotate("rect", xmin = 1999, xmax = 2009, ymin = 0, ymax = ymax,
             fill = "#a1d99b", alpha = 0.45) +
    annotate("rect", xmin = min(runs$cmp_start), xmax = max(runs$cmp_end),
             ymin = 0, ymax = ymax, fill = "#fdae6b", alpha = 0.30) +
    annotate("text", x = 2004, y = ymax * 0.97, label = "reference\n(pre-rust)",
             vjust = 1, size = 4, colour = "#2f6b2a") +
    annotate("text", x = (min(runs$cmp_start) + max(runs$cmp_end)) / 2, y = ymax * 0.97,
             label = "comparison\n(3 generations on;\nvaries by run)", vjust = 1, size = 4,
             colour = "#9a5400") +
    geom_line(data = tr[tr$run %in% show, ], aes(year, pmin(rel, ymax), group = run),
              colour = "#7f7f7f", linewidth = 0.25, alpha = 0.6) +
    geom_ribbon(data = q, aes(year, ymin = q05, ymax = pmin(q95, ymax)), fill = "#2a78d6", alpha = 0.13) +
    geom_ribbon(data = q, aes(year, ymin = q25, ymax = q75), fill = "#2a78d6", alpha = 0.25) +
    geom_line(data = q, aes(year, med), colour = "#1d5fae", linewidth = 1.2) +
    geom_hline(yintercept = 1, colour = "#2f6b2a", linetype = "dotted") +
    geom_vline(xintercept = YEAR_RUST, colour = "#d6452f", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = YEAR_ASSESS, colour = "#4a3aa7", linewidth = 0.8) +
    annotate("text", x = YEAR_RUST - 1, y = 0.08, label = "rust\n2010", hjust = 1,
             colour = "#d6452f", size = 4) +
    annotate("text", x = YEAR_ASSESS + 1, y = 0.08, label = "now", hjust = 0,
             colour = "#4a3aa7", size = 4) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0),
                       labels = function(x) paste0(round(100 * x), "%")) +
    scale_x_continuous(breaks = seq(1960, 2120, 20), expand = c(0, 0)) +
    labs(x = NULL, y = "Mature individuals\n(% of pre-rust reference)",
         title = "Each run: one population, one future",
         subtitle = sprintf(
           "%d runs across the range. Grey: %d individual runs (sharp drops are fires). Blue: median, 50%% and 90%% bands",
           length(unique(tr$run)), length(show)),
         caption = wrap(paste0("Mature individuals = alive plants of reproductive age, including those not flowering after fire or rust. ",
                          "Each run's comparison window ends 3 x its own generation time after 2009.\n",
                          draws_caption()))) +
    theme_talk()
}

fig_outcomes <- function(R) {
  df <- R$runs
  med <- median(df$pct_decline, na.rm = TRUE)
  shares <- function(x) as.numeric(prop.table(table(x)))
  sh <- shares(df$category)
  lo <- min(0, 10 * floor(min(df$pct_decline, na.rm = TRUE) / 10))
  ytop <- max(graphics::hist(df$pct_decline, breaks = seq(lo - 2.5, 102.5, 2.5),
                             plot = FALSE)$counts)

  pA <- ggplot(df, aes(pct_decline, fill = category)) +
    geom_histogram(binwidth = 2.5, boundary = 0, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = c(30, 50, 80), colour = "#9a9a9a", linetype = "dashed") +
    annotate("text", x = c(30, 50, 80), y = Inf, label = c("VU", "EN", "CR"),
             hjust = -0.2, vjust = 1.4, colour = CAT_COL[c("VU", "EN", "CR")],
             fontface = "bold", size = 5) +
    geom_vline(xintercept = med, colour = INK, linewidth = 1) +
    annotate("text", x = med + 1.2, y = ytop * 1.12, hjust = 0,
             label = sprintf("median %.0f%%", med), size = 4.3, colour = INK) +
    scale_fill_manual(values = CAT_COL, drop = FALSE, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_x_continuous(limits = c(lo, 100), labels = function(x) paste0(x, "%")) +
    labs(x = "Projected decline in mature individuals over three generations",
         y = "Runs") +
    theme_talk()

  mk <- function(g, x) if (length(x)) data.frame(group = g, category = x) else NULL
  grp <- rbind(
    mk("All runs",   df$category),
    mk("NSW cells",  df$category[df$region == "NSW"]),
    mk("Qld cells*", df$category[df$region == "Qld"]))
  grp$group <- factor(grp$group, levels = rev(intersect(c("All runs", "NSW cells", "Qld cells*"), grp$group)))
  cnt <- as.data.frame(table(grp$group, grp$category))
  names(cnt) <- c("group", "category", "n")
  cnt$share <- cnt$n / ave(cnt$n, cnt$group, FUN = sum)
  ntot <- tapply(cnt$n, cnt$group, sum)
  levels(cnt$group) <- sprintf("%s\n(n = %d)", levels(cnt$group), ntot[levels(cnt$group)])

  pB <- ggplot(cnt, aes(share, group, fill = category)) +
    geom_col(width = 0.65, colour = "white", position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = ifelse(share >= 0.06, sprintf("%.0f%%", 100 * share), "")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              colour = "white", fontface = "bold", size = 4.3) +
    scale_fill_manual(values = CAT_COL, drop = FALSE, name = NULL) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = "Share of runs", y = NULL) +
    theme_talk() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  (pA | pB) + plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
      title = sprintf("Range-wide outcome: %.0f%% of runs threatened  (VU %.0f%%, EN %.0f%%, CR %.0f%%)",
                      100 * mean(df$pct_decline >= 30, na.rm = TRUE),
                      100 * sh[2], 100 * sh[3], 100 * sh[4]),
      subtitle = sprintf("IUCN criterion A4, %d runs: grid cells drawn in proportion to abundance x simulation seed x parameter draws",
                         nrow(df)),
      caption = wrap(paste0("* Qld cells: no Qld fire history yet; fire frequency and extent there come from the prior. ",
                       "Rust pressure from the pathogen SDM everywhere. No resistance yet (fully susceptible).\n",
                       draws_caption()), 190),
      theme = theme_talk())
}

fig_map <- function(R) {
  df <- R$runs
  sites <- load_sites()
  cell <- do.call(rbind, lapply(split(df, df$cell_id), function(d) data.frame(
    cell_id = d$cell_id[1], lon = d$lon[1], lat = d$lat[1], n = nrow(d),
    med = median(d$pct_decline, na.rm = TRUE))))
  aus <- ggplot2::map_data("world", region = "Australia")
  ggplot() +
    geom_polygon(data = aus, aes(long, lat, group = group), fill = "#f1efe9", colour = "#bdb6a6",
                 linewidth = 0.3) +
    geom_hline(yintercept = QLD_LAT, colour = "#9a9a9a", linetype = "dashed", linewidth = 0.4) +
    annotate("text", x = 148.3, y = QLD_LAT + 0.25, label = "Qld (fire: prior only)",
             hjust = 0, size = 3.8, colour = "#777777") +
    annotate("text", x = 148.3, y = QLD_LAT - 0.25, label = "NSW", hjust = 0, size = 3.8,
             colour = "#777777") +
    geom_point(data = sites, aes(lon_centre, lat_centre), colour = "#c9c9c9", size = 0.8) +
    geom_point(data = cell, aes(lon, lat, colour = med, size = n), alpha = 0.9) +
    scale_colour_gradientn(colours = c("#8fb3cf", "#eda100", "#d6452f", "#7a0f12"),
                           values = scales::rescale(c(0, 30, 50, 80, 100)),
                           limits = c(min(0, min(cell$med)), 100),
                           name = "Median\ndecline (%)") +
    scale_size_area(max_size = 6, name = "Runs in cell") +
    coord_quickmap(xlim = c(148, 154.2), ylim = c(-37, -24)) +
    labs(x = NULL, y = NULL, title = "Where the runs come from",
         subtitle = "Cells sampled in proportion to abundance;\ncolour = median decline of the runs in that cell",
         caption = "Grey dots: all occupied cells.") +
    theme_talk() + theme(panel.grid = element_blank(), axis.text = element_blank())
}

write_numbers <- function(R, traj) {
  df <- R$runs
  sh <- prop.table(table(df$category))
  at22 <- traj$cum_dead[traj$model_year == 2021]
  lines <- c(
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    sprintf("Model commit: %s", tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                                         error = function(e) "unknown")),
    sprintf("Parameter draws: %s", draws_caption()),
    "",
    sprintf("PEGG: %d runs; resprouters dead by Jan 2022: median %.1f%% (50%%: %.1f-%.1f; 90%%: %.1f-%.1f); observed 48%% (24/50, 95%% CI 33.7-62.6%%)",
            length(unique(traj$run)), 100 * median(at22),
            100 * quantile(at22, .25), 100 * quantile(at22, .75),
            100 * quantile(at22, .05), 100 * quantile(at22, .95)),
    sprintf("PEGG: trees killed outright by the fire (model): median %.1f%%",
            100 * median(tapply(traj$frac_fire_killed, traj$run, `[`, 1))),
    "",
    sprintf("RANGE: %d runs over %d distinct cells", nrow(df), length(unique(df$cell_id))),
    sprintf("RANGE: decline median %.1f%%, IQR %.1f-%.1f%%, 90%% %.1f-%.1f%%",
            median(df$pct_decline, na.rm = TRUE),
            quantile(df$pct_decline, .25, na.rm = TRUE), quantile(df$pct_decline, .75, na.rm = TRUE),
            quantile(df$pct_decline, .05, na.rm = TRUE), quantile(df$pct_decline, .95, na.rm = TRUE)),
    sprintf("RANGE: category shares  LC/NT %.0f%%  VU %.0f%%  EN %.0f%%  CR %.0f%%",
            100 * sh[["LC/NT"]], 100 * sh[["VU"]], 100 * sh[["EN"]], 100 * sh[["CR"]]),
    sprintf("RANGE: threatened (>= 30%%): %.0f%% of runs", 100 * mean(df$pct_decline >= 30, na.rm = TRUE)),
    sprintf("RANGE: generation time median %.1f yr; window median %.0f yr",
            median(df$T_g, na.rm = TRUE), median(df$window_len, na.rm = TRUE)),
    sprintf("RANGE: flowering-only decline median %.1f%%", median(df$pct_decline_flowering, na.rm = TRUE)),
    sprintf("RANGE: Qld-cell runs: %d (fire from prior)", sum(df$region == "Qld"))
  )
  for (r in intersect(c("NSW", "Qld"), df$region)) {
    s <- prop.table(table(df$category[df$region == r]))
    lines <- c(lines, sprintf("RANGE %s: LC/NT %.0f%%  VU %.0f%%  EN %.0f%%  CR %.0f%%  (median %.1f%%)",
                              r, 100 * s[["LC/NT"]], 100 * s[["VU"]], 100 * s[["EN"]], 100 * s[["CR"]],
                              median(df$pct_decline[df$region == r], na.rm = TRUE)))
  }
  writeLines(lines, file.path(FIG_DIR, "talk_numbers.txt"))
  cat(lines, sep = "\n")
}

make_figures <- function() {
  traj <- readRDS(file.path(RESULTS_DIR, "pegg_trajectories.rds"))
  R    <- collect_range()
  saveRDS(R, file.path(RESULTS_DIR, "range_collected.rds"))
  fig_diagram()
  save_fig(fig_pegg(traj),       "fig2_pegg_fit.png")
  save_fig(fig_trajectories(R),  "fig3_trajectories.png")
  save_fig(fig_outcomes(R),      "fig4_outcomes.png", w = 13, h = 6.75)
  save_fig(fig_map(R),           "fig5_map.png", w = 7.5, h = 8.5)
  write_numbers(R, traj)
}

# ============================================================================
# Main
# ============================================================================
if (sys.nframe() == 0L) {
  stage <- commandArgs(trailingOnly = TRUE)
  stage <- if (length(stage)) stage[1] else "all"
  if (stage %in% c("pegg", "all"))    run_pegg_stage()
  if (stage %in% c("range", "all"))   run_range_stage()
  if (stage %in% c("figures", "all")) make_figures()
}
