
# https://grantdadams.github.io/Rceattle/articles/model-diagnostics.html
library(Rceattle)
library(ggplot2)

# Data ----
data_2025_c_la_wa_ae <- Rceattle::read_data( file = "Data/GOA_25.1.1_arrowtooth_single_species_1977-2024_new_data_clean_LA_WA_AE.xlsx")
data_2025_c_la_wa_ae$estDynamics = 0
data_2025_c_la_wa_ae$index_data$Log_sd <- data_2025_c_la_wa_ae$index_data$Log_sd/data_2025_c_la_wa_ae$index_data$Observation

# Fit model -----
# ceattle_ss_2025_RE_c_la_wa_ae
mod_25_0 <- Rceattle::fit_mod(data_list = data_2025_c_la_wa_ae,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1)

# Diagnostics ----
# * Summaries -----
summary(mod_25_0)
convergence_diagnostics(mod_25_0)

# * Plots ----
plot_index(mod_25_0)
plot_logindex(mod_25_0)
plot_indexresidual(mod_25_0)
plot_comp(mod_25_0)
plot_catch(mod_25_0)

# * OSAs ----
osa <- osa_residuals(mod_25_0)
head(osa)

# Statistical diagnostics (Stewart & Monnahan 2025): SDNR and lower/upper tail
osa_diagnostics(osa)

# Q-Q plot (with SDNR / tail annotation) + residual-by-year
plot(osa)


# * Retrospectives ----
mod_25_0_retro <- retrospective(Rceattle = mod_25_0, peels = 5)
mod_25_0_retro$mohns # Mohn's rho for each quantity

# Plot historical trajectories across peels
plot_biomass(mod_25_0_retro$Rceattle_list)


# * Jitters ----
mod_25_0_jitters <- jitter(Rceattle = mod_25_0, njitter = 100, phase = TRUE)

# Histogram of NLL differences relative to the best run
hist(log(mod_25_0_jitters$nll - min(mod_25_0_jitters$nll)),
     main = "Jitter NLL spread (log scale)",
     xlab = "log(NLL - min NLL)")

# Overlay biomass trajectories — tight overlap indicates a stable optimum
plot_biomass(mod_25_0_jitters$Rceattle_list) + theme(legend.position="none")


# * Self-test ----
mod_25_0_sims <- self_test(mod_25_0, nsim = 100, start = "estimated")

# Number of simulations that converged (non-converged runs are dropped)
length(mod_25_0_sims)

# Overlay biomass / SSB trajectories across simulations — the original fit's
# trajectory should sit inside the spread of the refits.
# One grey per simulation KEPT, not per simulation requested: self_test() drops
# the ones that did not converge, so a hardcoded 100 leaves the black `1` past
# the end of the list and draws the fit itself grey, among the greys.
plot_biomass(c(mod_25_0_sims, list(mod_25_0)),
             line_col = c(rep("grey", length(mod_25_0_sims)), 1)) +
  theme(legend.position = "none")

# * Likelihood profiles ----
# sigmaR
prof_sigmaR <- profile(
  fitted = mod_25_0,
  param    = "sigmaR",
  slots    = list(1),
  values   = list(seq(0.1, 1.5, by = 0.05))
)

plot(prof_sigmaR$grid$slot_1,
     prof_sigmaR$nll - min(prof_sigmaR$nll, na.rm = TRUE),
     type = "l", xlab = "sigmaR", ylab = "dNLL")


# 2-D cross-profile: M1 across sex
prof_M_sex <- profile(
  fitted = mod_25_0,
  param    = "M1",
  slots    = list(c(1, 1, 1), c(1, 2, 1)),  # males and females
  values   = list(seq(0.10, 0.40, by = 0.05),
                  seq(0.10, 0.40, by = 0.05))
)
