# =============================================================================
# Internal (Rceattle / TMB::oneStepPredict) vs external (afscOSA / compResidual)
# OSA residuals -- GOA arrowtooth, model 25.0
#
# The two are different estimators of the same quantity, so this script separates
# differences that are BOOKKEEPING (must be zero, or the comparison is
# meaningless) from ones that are STRUCTURAL (expected, and are the point).
#
# ---- BOOKKEEPING: get these right or nothing else means anything -------------
#  1. Bin order. Rceattle lays a joint-sex (Sex == 3) composition out as ONE
#     multinomial of 2*nbin bins: females in bins 1..nbin, males in
#     nbin+1..2*nbin. Both ATF comp fleets are Sex == 3, so the survey age comp
#     is a single 42-bin multinomial and the fishery length comp a single 58-bin
#     one -- NOT two per-sex multinomials.
#  2. Which bin is dropped. compResidual::resMulti() drops the LAST row of the
#     matrix it is handed; Rceattle drops the bin flagged `is_last_bin`, also the
#     last. For joint sex that is the MALE plus group only. Splitting by sex
#     instead gives 2*(nbin-1) residuals from a different decomposition
#     (40 vs 41 here) -- see the ordering checks at the bottom.
#  3. Residual labelling. resMulti() residual k is the conditional for bin k
#     given bins 1..k-1, so labels are index[1:(n-1)], the FIRST n-1 bins.
#     afscOSA::run_osa() does this correctly. NOTE:
#     R/Functions/Plot osa function.R uses `index[-1]` (the LAST n-1), which
#     shifts every bubble one bin -- fix that before trusting those figures.
#  4. Sample size. Rceattle's OSA counts use RAW Sample_size (comp_n[, 2]), NOT
#     Sample_size * Comp_weights: the OSA rebuild uses unweighted densities. So
#     pass raw Sample_size as afscOSA's `N` (0.106-0.119 weights here would
#     otherwise shrink N ~9x).
#  5. Row alignment. comp_obs / comp_hat / comp_ctl are all in comp_data row
#     order, so they join 1:1. Read them from mod$obj$env$data (the rearranged
#     TMB data) -- mod$data_list is the PRE-rearrange input and has no comp_obs.
#  6. Likelihood. Comp_loglike defaults to MultinomialAFSC here, and the OSA
#     rebuild residualizes the AFSC pseudo-likelihood under the full multinomial
#     -- so resMulti() (not resDirM()/`theta`) is the right external analogue.
#
# ---- STRUCTURAL: expected to differ -----------------------------------------
#  7. Randomization. resMulti() residuals are randomized-quantile (the counts are
#     integers), so they are STOCHASTIC. With N=200 spread over 42 bins many
#     cells hold 0-2 fish, the binomial CDF jumps are large, and the
#     randomization dominates. This is the single biggest source of
#     internal-vs-external disagreement, and it is noise in the EXTERNAL
#     estimator, not error in either. Quantified below.
#  8. Random effects. With random_rec = TRUE observations are not independent
#     given the fixed effects; internal OSA re-integrates the random effects as
#     each observation is added. External OSA conditions on the fitted comp_hat
#     and treats every composition as independent.
#  9. Cross-source conditioning. Internal residualizes index (then catch) first,
#     then compositions ordered by YEAR, then fleet, then bin -- so the two comp
#     fleets are interleaved chronologically and each is conditional on the
#     survey index. External runs each fleet in isolation.
# 10. Counts. Internal: (proportion + comp_offset) * N, continuous
#     (discrete = FALSE). External: round(N * proportion), integer. Negligible
#     in practice (checked below: cor > 0.99).
# =============================================================================

library(Rceattle)
library(afscOSA)       # remotes::install_github('noaa-afsc/afscOSA')
library(compResidual)  # remotes::install_github('fishfollower/compResidual/compResidual')
library(dplyr)
library(tidyr)
library(ggplot2)

out_dir <- "Results/Figures/Diagnostics/OSA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito, assigned in fixed order and never cycled
col_int <- "#0072B2"; col_ext <- "#D55E00"; col_bar <- "#009E73"


# ---- Fit --------------------------------------------------------------------
data_2025 <- Rceattle::read_data(
  file = "Data/GOA_25.1.1_arrowtooth_single_species_1977-2024_new_data_clean_LA_WA_AE.xlsx")
data_2025$estDynamics <- 0
data_2025$index_data$Log_sd <- data_2025$index_data$Log_sd / data_2025$index_data$Observation

mod_25_0 <- Rceattle::fit_mod(data_list = data_2025, estimateMode = 0,
                              random_rec = TRUE, msmMode = 0, initMode = 1,
                              fit_control = fit_control(verbose = 1, phase = TRUE))

# Fixed-effect twin. With no random effects the observations ARE independent
# given the parameters, so (8) and (9) drop out and internal vs external should
# agree up to (7). This is the CONTROL that tells you whether disagreement in
# mod_25_0 is real signal or a bookkeeping mistake.
mod_25_0_fe <- Rceattle::fit_mod(data_list = data_2025, estimateMode = 0,
                                 random_rec = FALSE, msmMode = 0, initMode = 1,
                                 fit_control = fit_control(verbose = 1, phase = TRUE))


# ---- Pull each composition block in Rceattle's exact bin order ---------------
osa_block <- function(mod, fleet_code, sex, comp_type) {
  tmb  <- mod$obj$env$data                # rearranged TMB data
  ctl  <- tmb$comp_ctl                    # Fleet_code, Species, Sex, Age0_Length1, Year
  rows <- which(ctl[, 1] == fleet_code & ctl[, 3] == sex & ctl[, 4] == comp_type &
                  tmb$comp_n[, 2] > 0 & ctl[, 5] > 0 & ctl[, 5] <= tmb$endyr)
  sp     <- ctl[rows[1], 2]
  nbin   <- if (comp_type == 0) tmb$nages[sp] else tmb$nlengths[sp]
  n_comp <- nbin * if (sex == 3) 2L else 1L      # joint sex = ONE 2*nbin multinomial
  fc <- mod$data_list$fleet_control
  list(fleet_code = fleet_code, sex = sex, comp_type = comp_type,
       fleet_name = fc$Fleet_name[match(fleet_code, fc$Fleet_code)],
       nbin = nbin, n_comp = n_comp,
       index_label = if (comp_type == 0) "age" else "length",
       years  = ctl[rows, 5],
       N      = tmb$comp_n[rows, 2],                                  # RAW Sample_size
       obs    = tmb$comp_obs[rows, seq_len(n_comp), drop = FALSE],    # row-normalised
       exp    = mod$quantities$comp_hat[rows, seq_len(n_comp), drop = FALSE],
       offset = tmb$comp_offset)
}

comp_blocks <- function(mod) {
  tmb <- mod$obj$env$data
  u <- unique(tmb$comp_ctl[tmb$comp_n[, 2] > 0, c(1, 3, 4), drop = FALSE])
  lapply(seq_len(nrow(u)), function(i) osa_block(mod, u[i, 1], u[i, 2], u[i, 3]))
}

# (year, bin) key for a block's residual matrix, in resMulti()'s output order
block_key <- function(b) data.frame(
  fleet = b$fleet_code, sex = b$sex,
  year = rep(b$years, times = b$n_comp - 1L),
  age_length_bin = rep(seq_len(b$n_comp - 1L), each = length(b$years)))


# ---- External OSA, two flavours ---------------------------------------------
# (a) stock afscOSA::run_osa() with correctly ordered input
ext_afscosa <- function(b, seed = 99801) {
  o <- afscOSA::run_osa(obs = b$obs, exp = b$exp, N = b$N, fleet = b$fleet_name,
                        index = seq_len(b$n_comp), years = b$years,
                        index_label = b$index_label, seed = seed)
  if (is.null(o)) return(NULL)
  o$res %>% transmute(fleet = b$fleet_code, sex = b$sex, year,
                      age_length_bin = index, ext_afscosa = resid)
}

# (b) resMulti() fed Rceattle's own obsvec counts, (proportion + offset) * N, to
#     isolate afscOSA's round()/no-offset input construction (10) from the
#     conditioning differences (8, 9).
ext_matched <- function(b, seed = 99801) {
  o <- round((b$obs + b$offset) * b$N, 0)
  p <- b$exp + b$offset; p <- p / rowSums(p)
  set.seed(seed); res <- compResidual::resMulti(t(o), t(p))
  if (!all(is.finite(res))) { warning("non-finite matched OSA: ", b$fleet_name); return(NULL) }
  cbind(block_key(b),
        ext_matched = as.vector(t(matrix(res, nrow = nrow(res), ncol = ncol(res)))))
}

# (c) THE NOISE FLOOR. resMulti() residuals are randomized-quantile, so a single
#     call is one draw from a distribution. Average over `S` seeds to get the
#     expected external residual (`ext_bar`), and keep the draws so the
#     seed-to-seed correlation can be used as the yardstick for "how well could
#     ANY two OSA calculations agree on these data?".
ext_draws <- function(b, S = 30) {
  o <- round(b$N * b$obs / rowSums(b$obs), 0)
  p <- b$exp / rowSums(b$exp)
  sapply(seq_len(S), function(s) {
    set.seed(s); r <- compResidual::resMulti(t(o), t(p))
    as.vector(t(matrix(r, nrow = nrow(r), ncol = ncol(r))))
  })
}


# ---- Assemble ---------------------------------------------------------------
compare_osa <- function(mod, label, S = 30) {
  blocks <- comp_blocks(mod)
  int <- Rceattle::osa_residuals(mod, parallel = FALSE) %>%
    filter(source == "comp") %>%
    transmute(fleet, sex, year, age_length_bin, internal = residual)

  bar <- bind_rows(lapply(blocks, function(b) {
    E <- ext_draws(b, S)
    cbind(block_key(b), ext_bar = rowMeans(E), ext_sd_seed = apply(E, 1, sd))
  }))
  meta <- bind_rows(lapply(blocks, function(b) data.frame(
    fleet = b$fleet_code, sex = b$sex, fleet_name = b$fleet_name,
    nbin = b$nbin, index_label = b$index_label)))

  key <- c("fleet", "sex", "year", "age_length_bin")
  int %>%
    full_join(bind_rows(lapply(blocks, ext_afscosa)), by = key) %>%
    full_join(bind_rows(lapply(blocks, ext_matched)), by = key) %>%
    full_join(bar, by = key) %>%
    left_join(meta, by = c("fleet", "sex")) %>%
    mutate(model   = label,
           # split joint-sex bins back onto one age/length axis for plotting
           sex_lab = ifelse(sex == 3 & age_length_bin > nbin, "Male", "Female"),
           bin     = ifelse(sex == 3 & age_length_bin > nbin,
                            age_length_bin - nbin, age_length_bin)) %>%
    arrange(fleet, year, age_length_bin)
}

cmp <- bind_rows(
  compare_osa(mod_25_0_fe, "FE (random_rec = FALSE)"),
  compare_osa(mod_25_0,    "RE (random_rec = TRUE)"))


# ---- Agreement, against the right yardstick ---------------------------------
# Read `cor_int_bar` against `cor_draw_bar`: if the internal residual tracks the
# seed-averaged external residual AT LEAST as well as a single external draw
# does, the two estimators agree as well as the external one agrees with itself.
noise_floor <- function(mod, label, S = 30) {
  int_all <- Rceattle::osa_residuals(mod, source = "comp", parallel = FALSE)
  bind_rows(lapply(comp_blocks(mod), function(b) {
    E <- ext_draws(b, S); k <- block_key(b)
    d <- int_all %>% filter(fleet == b$fleet_code) %>%
      transmute(year, age_length_bin, int = residual) %>%
      merge(cbind(k, ebar = rowMeans(E), e1 = E[, 1]), by = c("year", "age_length_bin"))
    Em <- E[match(paste(d$year, d$age_length_bin), paste(k$year, k$age_length_bin)), ]
    sc <- cor(Em)[upper.tri(cor(Em))]
    data.frame(model = label, fleet = b$fleet_name, n = nrow(d),
               cor_int_draw = cor(d$int, d$e1),        # internal vs ONE external draw
               cor_draw_draw = mean(sc),               # two external draws (noise floor)
               floor_lo = quantile(sc, 0.025), floor_hi = quantile(sc, 0.975),
               cor_int_bar  = cor(d$int, d$ebar),      # internal vs expected external
               cor_draw_bar = mean(apply(Em, 2, cor, y = d$ebar)),
               sdnr_int = sd(d$int), sdnr_draw = sd(d$e1), sdnr_bar = sd(d$ebar))
  }))
}
osa_noise <- bind_rows(noise_floor(mod_25_0_fe, "FE (random_rec = FALSE)"),
                       noise_floor(mod_25_0,    "RE (random_rec = TRUE)"))
rownames(osa_noise) <- NULL
print(as.data.frame(osa_noise), digits = 3)

# Simple pairwise agreement, incl. the afscOSA-vs-matched check that (10) is
# negligible.
agreement <- function(d, a, b) {
  ok <- is.finite(d[[a]]) & is.finite(d[[b]]); x <- d[[a]][ok]; y <- d[[b]][ok]
  data.frame(n = length(x), pearson = cor(x, y), spearman = cor(x, y, method = "spearman"),
             sign_agree = mean(sign(x) == sign(y)),
             mean_abs_diff = mean(abs(x - y)), max_abs_diff = max(abs(x - y)))
}
osa_agreement <- cmp %>% group_by(model, fleet_name) %>%
  group_modify(~ bind_rows(
    cbind(pair = "internal vs afscOSA",    agreement(.x, "internal", "ext_afscosa")),
    cbind(pair = "internal vs ext_bar",    agreement(.x, "internal", "ext_bar")),
    cbind(pair = "afscOSA vs matched(10)", agreement(.x, "ext_afscosa", "ext_matched")))) %>%
  ungroup()
print(as.data.frame(osa_agreement), digits = 3)

# SDNR: the internal comp residuals are deterministic and carry NO randomization
# variance, so their null SDNR is below 1. Do NOT read internal comp SDNR
# against the N(0,1) chi-square interval osa_diagnostics() prints -- it will
# look under-dispersed for that reason alone.
sdnr_tbl <- cmp %>%
  pivot_longer(c(internal, ext_afscosa, ext_bar), names_to = "method", values_to = "resid") %>%
  filter(is.finite(resid)) %>%
  group_by(model, fleet_name, method) %>%
  summarise(n = n(), sdnr = sd(resid), mean = mean(resid),
            lower = quantile(resid, .025), upper = quantile(resid, .975), .groups = "drop") %>%
  mutate(sdnr_lo = sqrt(qchisq(.025, n - 1) / (n - 1)),
         sdnr_hi = sqrt(qchisq(.975, n - 1) / (n - 1)))
print(as.data.frame(sdnr_tbl), digits = 3)


# ---- Plots ------------------------------------------------------------------
p_scatter <- cmp %>%
  pivot_longer(c(ext_afscosa, ext_bar), names_to = "external", values_to = "ext") %>%
  filter(is.finite(internal), is.finite(ext)) %>%
  mutate(external = recode(external, ext_afscosa = "afscOSA (single seed)",
                           ext_bar = "external, averaged over 30 seeds")) %>%
  ggplot(aes(ext, internal)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(alpha = 0.3, size = 1, colour = col_int) +
  facet_grid(model ~ fleet_name + external) +
  labs(x = "External OSA residual", y = "Internal OSA residual",
       title = "Internal (Rceattle) vs external (afscOSA) OSA residuals",
       subtitle = "Grey line is 1:1. Scatter against a single seed is dominated by resMulti() randomization.") +
  coord_equal() + theme_classic(base_size = 9)
ggsave(file.path(out_dir, "osa_internal_vs_external_scatter.png"), p_scatter,
       width = 10, height = 6, dpi = 300)

p_qq <- cmp %>%
  pivot_longer(c(internal, ext_afscosa, ext_bar), names_to = "method", values_to = "resid") %>%
  filter(is.finite(resid)) %>%
  ggplot(aes(sample = resid, colour = method)) +
  stat_qq(size = .7, alpha = .45) + stat_qq_line(colour = "grey40", linewidth = .4) +
  facet_grid(model ~ fleet_name) +
  scale_colour_manual(values = c(internal = col_int, ext_afscosa = col_ext, ext_bar = col_bar),
                      name = NULL) +
  labs(x = "Theoretical N(0,1) quantile", y = "Sample quantile",
       title = "OSA residual normal QQ plots by estimator") +
  theme_classic(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(out_dir, "osa_internal_vs_external_qq.png"), p_qq,
       width = 9, height = 6, dpi = 300)

p_bub <- cmp %>%
  filter(model == "RE (random_rec = TRUE)") %>%
  pivot_longer(c(internal, ext_afscosa), names_to = "method", values_to = "resid") %>%
  filter(is.finite(resid)) %>%
  ggplot(aes(year, bin, size = abs(resid), colour = resid > 0)) +
  geom_point(alpha = .7) +
  facet_grid(fleet_name + sex_lab ~ method, scales = "free_y") +
  scale_colour_manual(values = c(`FALSE` = col_int, `TRUE` = col_ext),
                      labels = c("negative", "positive"), name = NULL) +
  scale_size_area(max_size = 4, name = "|residual|") +
  labs(x = "Year", y = "Age / length bin",
       title = "OSA residual bubbles: internal vs afscOSA (model 25.0)") +
  theme_classic(base_size = 9)
ggsave(file.path(out_dir, "osa_internal_vs_external_bubbles.png"), p_bub,
       width = 9, height = 9, dpi = 300)


# =============================================================================
# Ordering checks -- what actually moves the numbers
# =============================================================================
b1 <- comp_blocks(mod_25_0)[[1]]
o1 <- round(b1$N * b1$obs / rowSums(b1$obs), 0); p1 <- b1$exp / rowSums(b1$exp)
# Compare on SEED-AVERAGED residuals. A single resMulti() call draws from one
# RNG stream consumed in column order, so permuting years re-deals the
# randomization and a raw comparison measures the RNG, not the conditional
# structure. Averaging removes it. `S = 200` here; the MC reference below is the
# yardstick for "indistinguishable".
bar_ <- function(o, p, S = 200) Reduce(`+`, lapply(seq_len(S), function(s) {
  set.seed(s); as.matrix(compResidual::resMulti(t(o), t(p))) })) / S

# (i) YEAR order is irrelevant externally: resMulti() treats each column (year)
#     as an independent multinomial, so residuals just permute with the years.
perm  <- sample(seq_along(b1$years))
b_ord <- bar_(o1, p1); b_prm <- bar_(o1[perm, ], p1[perm, ])
cat("year-order invariance (external): cor =", cor(as.vector(b_ord[, perm]), as.vector(b_prm)),
    " mean|diff| =", mean(abs(b_ord[, perm] - b_prm)), "\n")
# MC reference: two independent seed-averages of the SAME ordering. If the
# year-order mean|diff| above is <= this, the orderings are indistinguishable.
mc <- function(seeds) Reduce(`+`, lapply(seeds, function(s) {
  set.seed(s); as.matrix(compResidual::resMulti(t(o1), t(p1))) })) / length(seeds)
cat("  MC reference (same order, independent seed sets): mean|diff| =",
    mean(abs(mc(1:100) - mc(201:300))), "\n")

# (ii) BIN order DOES matter -- it defines the conditional sequence and which bin
#      is dropped.
rv <- rev(seq_len(b1$n_comp))
cat("bin-order sensitivity (external), cor(fwd, rev) =",
    cor(as.vector(b_ord), rev(as.vector(bar_(o1[, rv], p1[, rv])))), "\n")

# (iii) The joint-sex trap: splitting a Sex == 3 composition into two
#       multinomials is a DIFFERENT decomposition and yields 2*(nbin-1)
#       residuals, not 2*nbin-1.
fi <- seq_len(b1$nbin); mi <- b1$nbin + fi
n_split <- sum(sapply(list(fi, mi), function(ii) {
  set.seed(1)
  length(compResidual::resMulti(
    t(round(b1$N * b1$obs[, ii] / rowSums(b1$obs[, ii]), 0)),
    t(b1$exp[, ii] / rowSums(b1$exp[, ii]))))
}))
cat("joint-sex: joint n =", length(b_ord), " vs per-sex split n =", n_split, "\n")

# (iv) Internal conditioning order. osa_residuals() hardcodes
#      source -> year -> fleet -> bin (WHAM-style). Re-run oneStepPredict() with
#      a chronological order to see how much the RE conditioning sequence moves
#      individual residuals. Slow -- run when you care.
osa_reorder <- function(mod, ord = c("wham", "chronological")) {
  ord <- match.arg(ord)
  osa_dat <- Rceattle:::build_osa_data(mod$obj$env$data, build_osa = TRUE)
  sel <- osa_dat$obs_ctl
  sel <- sel[sel$source %in% c("index", "catch", "comp") & !sel$is_last_bin, ]
  src <- c("index", "catch", "comp")
  sel <- if (ord == "wham") sel[order(match(sel$source, src), sel$year, sel$fleet_code, sel$bin_index), ]
         else               sel[order(sel$year, match(sel$source, src), sel$fleet_code, sel$bin_index), ]
  r <- TMB::oneStepPredict(Rceattle:::.osa_build_obj(mod, osa_dat),
                           observation.name = "obsvec", data.term.indicator = "keep",
                           method = "oneStepGaussianOffMode", subset = sel$obs_pos + 1L,
                           discrete = FALSE, parallel = FALSE, seed = 123)
  data.frame(obs_pos = sel$obs_pos, source = sel$source, fleet = sel$fleet_code,
             year = sel$year, age_length_bin = sel$age_length_bin, residual = r$residual)
}
ord_cmp <- merge(osa_reorder(mod_25_0, "wham"), osa_reorder(mod_25_0, "chronological"),
                 by = "obs_pos", suffixes = c("_wham", "_chron"))
cat("internal conditioning order (RE model): cor =",
    with(ord_cmp, cor(residual_wham, residual_chron)), " max|diff| =",
    with(ord_cmp, max(abs(residual_wham - residual_chron))), "\n")
