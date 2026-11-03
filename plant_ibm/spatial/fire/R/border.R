# spatial/fire/R/border.R -- do the NSW and Qld records agree where they meet?
#
# Cells within BORDER_BAND degrees of latitude of the border, each side,
# share climate and vegetation, so their fire regimes should be similar if
# the two records measure the same thing. Compares:
#   - pooled fire-year rate, sum(k) / sum(T), with a Poisson rate-ratio test
#   - share of cells with at least one fire year
#   - burnt fraction in fire years
#   - year by year: share of band cells with a fire year (shows changes in a
#     record over time, e.g. Qld's automated-mapping years before 2013)

suppressPackageStartupMessages(library(ggplot2))

border_compare <- function(reg, cells, cy, years_by_state) {
  band <- border_band(cells)
  reg <- reg[, setdiff(names(reg), c("state", "lat_centre", "lon_centre"))]
  d <- merge(reg, sf::st_drop_geometry(cells[, c("cell_id", "state", "lat_centre", "lon_centre")]),
             by = "cell_id")
  d$band <- d$cell_id %in% cells$cell_id[band]
  b <- d[d$band, ]
  k <- tapply(b$n_wildfire_years_sig, b$state, sum)
  T <- tapply(b$T_obs, b$state, sum)
  pt <- if (all(c("NSW", "Qld") %in% names(k)) && all(T > 0))
    stats::poisson.test(c(k[["Qld"]], k[["NSW"]]), c(T[["Qld"]], T[["NSW"]])) else NULL
  summ <- b |> dplyr::group_by(state) |> dplyr::summarise(
    cells = dplyr::n(), cell_years = sum(T_obs), fire_years = sum(n_wildfire_years_sig),
    rate_per_cell_year = fire_years / cell_years,
    share_cells_with_fire = mean(n_wildfire_years_sig > 0),
    median_frac_burned = stats::median(mean_frac_burned, na.rm = TRUE),
    median_FRI = stats::median(FRI_years), .groups = "drop")

  # year by year
  yy <- do.call(rbind, lapply(c("NSW", "Qld"), function(s) {
    ids <- b$cell_id[b$state == s]; yrs <- years_by_state[[s]]
    if (!length(ids)) return(NULL)
    sig <- cy[cy$cell_id %in% ids & cy$frac_burned >= AREA_THRESH, ]
    data.frame(state = s, year = yrs,
               share = vapply(yrs, function(y) length(unique(sig$cell_id[sig$year == y])) / length(ids), 0))
  }))
  list(summary = summ, rate_ratio = pt, by_year = yy, cells = d)
}

plot_border <- function(bc, scenario, path) {
  d <- bc$cells
  rr <- bc$rate_ratio
  sub <- if (is.null(rr)) "" else sprintf("Qld/NSW fire-year rate in the band: %.2f (95%% CI %.2f-%.2f)",
                                          rr$estimate, rr$conf.int[1], rr$conf.int[2])
  pA <- ggplot(d, aes(lat_centre, n_wildfire_years_sig / T_obs, colour = state)) +
    annotate("rect", xmin = QLD_BORDER_LAT - BORDER_BAND, xmax = QLD_BORDER_LAT + BORDER_BAND,
             ymin = -Inf, ymax = Inf, fill = "#eeeeee") +
    geom_vline(xintercept = QLD_BORDER_LAT, linetype = "dashed", colour = "#888888") +
    geom_jitter(width = 0.02, height = 0, alpha = 0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.6) +
    scale_colour_manual(values = c(NSW = "#1d5fae", Qld = "#c0392b")) +
    labs(x = "Latitude (cell centre)", y = "Fire years per year observed (k / T)",
         title = "Fire-year rate by latitude", subtitle = sub) + theme_minimal(base_size = 12)
  pB <- ggplot(bc$by_year, aes(year, share, colour = state)) +
    geom_line() + geom_point(size = 1) +
    scale_colour_manual(values = c(NSW = "#1d5fae", Qld = "#c0392b")) +
    labs(x = NULL, y = "Share of band cells with a fire year",
         title = sprintf("Border band (+/- %.1f deg), year by year", BORDER_BAND)) +
    theme_minimal(base_size = 12)
  p <- patchwork::wrap_plots(pA, pB, ncol = 1) +
    patchwork::plot_annotation(title = paste("Border check:", scenario))
  ggsave(path, p, width = 10, height = 9, dpi = 160, bg = "white")
}
