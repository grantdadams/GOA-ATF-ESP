# Phase-plane diagram: spawning biomass against fishing mortality, both scaled
# to their Tier 3 reference points (Guidelines 4.9.2, figure 5).
#
# The track is the assessment's own history, with the two projected years from
# SPM on the end, as the Guidelines ask for. Reference points come from SPM so
# the whole figure is on one basis.
#
# Both axes are ratios, which makes them scale-free: Rceattle's F is a
# multiplier on a selectivity that peaks above 1, but F and F35% share that
# scale, so it cancels. No need to normalise anything.

library(ggplot2)

# The Amendment 56 Tier 3 control rule, drawn on these axes.
#
# FOFL  = F35% * ramp   and   max FABC = F40% * ramp, where
#   ramp = min(1, (SSB/B40% - alpha) / (1 - alpha))
#
# The x axis is SSB/B35%, so SSB/B40% = x * (B35%/B40%) = x * 0.875. Both lines
# therefore share a rising limb and differ only in where they flatten: FOFL at
# 1, max FABC at F40%/F35%, both at SSB = B40%, which is x = 1.143.
control_rule <- function(x, f_ratio, b35_over_b40, alpha = 0.05) {
  ramp <- pmin(1, pmax(0, (x * b35_over_b40 - alpha) / (1 - alpha)))
  f_ratio * ramp
}


# fit          a fitted Rceattle model (supplies the historical track)
# dir          the SPM run folder (supplies reference points and the projection)
# alt          SPM harvest scenario for the projected years; 2 is the author's F
# label_years  which years to label. Default is decadal plus the ends and the
#              projection. "all" labels every year, which is unreadable for a
#              track as tightly clustered as arrowtooth's.
phase_plane <- function(fit, dir, alt = 2, label_years = NULL,
                        stock = "GOA arrowtooth flounder") {

  q <- fit$quantities
  d <- fit$data_list
  hind <- d$styr:d$endyr

  # Reference points, from SPM's summary.
  ref <- read.csv(file.path(dir, "spm_summary.csv"))
  rv  <- function(v) as.numeric(ref$value[match(v, ref$variable)])
  b35 <- rv("SSB_ofl"); b40 <- rv("SSB_40")
  f35 <- rv("F_ofl");   f40 <- rv("F_abc")

  # Historical track from the assessment.
  track <- data.frame(
    year   = hind,
    b_b35  = q$ssb[1, as.character(hind)] / b35,
    f_f35  = q$F_spp[1, as.character(hind)] / f35,
    period = "Assessment"
  )

  # The two projected years from SPM, averaged over its simulations.
  detail <- read.csv(file.path(dir, "spm_detail.csv"))
  proj_years <- c(d$endyr + 1, d$endyr + 2)
  proj <- do.call(rbind, lapply(proj_years, function(y) {
    rows <- detail$Alt == alt & detail$Year == y
    data.frame(year = y,
               b_b35 = mean(detail$SSB[rows]) / b35,
               f_f35 = mean(detail$F[rows]) / f35,
               period = "Projected")
  }))

  dat <- rbind(track, proj)
  if (is.null(label_years)) {
    label_years <- c(min(hind), seq(1980, d$endyr, by = 10), d$endyr, proj_years)
  } else if (identical(label_years, "all")) {
    label_years <- dat$year
  }

  # The control rule, across the whole x range.
  x <- seq(0, max(dat$b_b35) * 1.05, length.out = 400)
  rule <- rbind(
    data.frame(x = x, y = control_rule(x, 1,         b35 / b40), rule = "FOFL = F35%"),
    data.frame(x = x, y = control_rule(x, f40 / f35, b35 / b40), rule = "max FABC = F40%")
  )

  ggplot(dat, aes(b_b35, f_f35)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_line(data = rule, aes(x, y, colour = rule), linewidth = 1.1) +
    geom_path(colour = "grey35", linewidth = 0.3) +
    geom_point(aes(shape = period, fill = period), size = 1.9, colour = "grey20") +
    ggrepel::geom_text_repel(
      data = subset(dat, year %in% label_years), aes(label = year),
      size = 2.7, colour = "grey20", min.segment.length = 0,
      segment.size = 0.2, segment.colour = "grey60", max.overlaps = Inf) +
    scale_colour_manual(values = c("FOFL = F35%"     = "#d1615d",
                                   "max FABC = F40%" = "#5aa9b5"), name = NULL) +
    scale_shape_manual(values = c(Assessment = 21, Projected = 24), name = NULL) +
    scale_fill_manual(values = c(Assessment = "grey25", Projected = "white"), name = NULL) +
    expand_limits(x = 0, y = 0) +
    labs(x = expression(SSB / B[35 * "%"]), y = expression(F / F[35 * "%"]),
         title = "Phase-plane diagram of management trajectory", subtitle = stock) +
    theme_light(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}
