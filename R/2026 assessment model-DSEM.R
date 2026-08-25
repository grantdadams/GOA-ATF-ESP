# Libraries ----
remotes::install_version("dsem", version = "3.0.0") # Need this version
remotes::install_github("grantdadams/Rceattle@dsem-v5-integration")
library(Rceattle) # dsem-v5-integration
library(dsem)
library(dplyr)
library(ggplot2)

# Data ----
data_2026 <- Rceattle::read_data( file = "2026/data/GOA_arrowtooth_2026.xlsx")
data_2026$estDynamics = 0
data_2026$index_data$Log_sd <- data_2026$index_data$Log_sd/data_2026$index_data$Observation
# Return to previous comp weights
#data_2026$fleet_control$Comp_weights <- c(1,1,0.25)

plot_data(data_2026)

# SAFE (Hindcast) Model ----
# Base 25.0 = single-species TMB based Rceattle model with data improvements that fixes sex-specific M (females = 0.2 and males = 0.35) and treats annual recruitment as random effects.
# Alternative 25.1 = Model 25.0 but estimates sex-specific M instead of fixing M.

# * Fit Base 25.0 -----
# Fixed M
mod_25_0 <- Rceattle::fit_mod(data_list = data_2026,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1)


# * Fit Alternative 25.1 -----
# Estimate sex-specific M
mod_25_1 <- Rceattle::fit_mod(data_list = data_2026,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1,
                              M1Fun = build_M1(M1_model = 2))


# Diagnostics ----
# * M25.0 Diags ----
# * Summaries -----
summary(mod_25_0)
convergence_diagnostics(mod_25_0)

# * Plots ----
plot_index(mod_25_0)
plot_index(mod_25_0, log = TRUE)
plot_indexresidual(mod_25_0)
plot_comp(mod_25_0)
plot_catch(mod_25_0)

# * OSAs ----
osa_25_0 <- osa_residuals(mod_25_0)
head(osa_25_0)

# Statistical diagnostics (Stewart & Monnahan 2025): SDNR and lower/upper tail
osa_diagnostics(osa_25_0)

# Q-Q plot (with SDNR / tail annotation) + residual-by-year
plot(osa_25_0)

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
plot_biomass(c(mod_25_0_sims, list(mod_25_0)), line_col = c(rep("grey", 100), 1)) + theme(legend.position="none")

# * Likelihood profiles ----
# sigmaR
prof_sigmaR_25_0 <- profile(
  fitted = mod_25_0,
  param    = "sigmaR",
  slots    = list(1),
  values   = list(seq(0.1, 1.5, by = 0.05))
)

plot(prof_sigmaR_25_0$grid$slot_1,
     prof_sigmaR_25_0$nll - min(prof_sigmaR_25_0$nll, na.rm = TRUE),
     type = "l", xlab = "sigmaR", ylab = "dNLL")


# 2-D cross-profile: M1 across sex
prof_M_sex_25_0 <- profile(
  fitted = mod_25_0,
  param    = "M1",
  slots    = list(c(1, 1, 1), c(1, 2, 1)),  # males and females
  values   = list(seq(0.10, 0.40, by = 0.05),
                  seq(0.10, 0.40, by = 0.05))
)

#not sure about this one
plot(prof_M_sex_25_0$grid$slot_1,
     prof_M_sex_25_0$nll - min(prof_M_sex_25_0$nll, na.rm = TRUE),
     type = "l", xlab = "M", ylab = "dNLL")

# * M25.1 Diagnostics ----
# * Summaries -----
summary(mod_25_1)
convergence_diagnostics(mod_25_1)

# * Plots ----
plot_index(mod_25_1)
plot_index(mod_25_1, log = TRUE)
plot_indexresidual(mod_25_1)
plot_comp(mod_25_1)
plot_catch(mod_25_1)

# * OSAs ----
osa_25_1 <- osa_residuals(mod_25_1)
head(osa_25_1)

# Statistical diagnostics (Stewart & Monnahan 2025): SDNR and lower/upper tail
osa_diagnostics(osa_25_1)

# Q-Q plot (with SDNR / tail annotation) + residual-by-year
plot(osa_25_1)

# * Retrospectives ----
mod_25_1_retro <- retrospective(Rceattle = mod_25_1, peels = 5)
mod_25_1_retro$mohns # Mohn's rho for each quantity

# Plot historical trajectories across peels
plot_biomass(mod_25_1_retro$Rceattle_list)
# plot_biomass(mod_25_1_retro$Rceattle_list, model_names = paste("Arrowtooth; mohns =", round(mod_25_1_retro$mohns[1,2], 3)))


# * Jitters ----
mod_25_1_jitters <- jitter(Rceattle = mod_25_1, njitter = 100, phase = TRUE)

# Histogram of NLL differences relative to the best run
hist(log(mod_25_1_jitters$nll - min(mod_25_1_jitters$nll)),
     main = "Jitter NLL spread (log scale)",
     xlab = "log(NLL - min NLL)")

# Overlay biomass trajectories — tight overlap indicates a stable optimum
plot_biomass(mod_25_1_jitters$Rceattle_list) + theme(legend.position="none")


# * Self-test ----
mod_25_1_sims <- self_test(mod_25_1, nsim = 100, start = "estimated")

# Number of simulations that converged (non-converged runs are dropped)
length(mod_25_1_sims)

# Overlay biomass / SSB trajectories across simulations — the original fit's
# trajectory should sit inside the spread of the refits.
plot_biomass(c(mod_25_1_sims, list(mod_25_1)), line_col = c(rep("grey", 100), 1)) + theme(legend.position="none")

# * Likelihood profiles ----
# sigmaR
prof_sigmaR_25_1 <- profile(
  fitted = mod_25_1,
  param    = "sigmaR",
  slots    = list(1),
  values   = list(seq(0.1, 1.5, by = 0.05))
)

plot(prof_sigmaR_25_1$grid$slot_1,
     prof_sigmaR_25_1$nll - min(prof_sigmaR_25_1$nll, na.rm = TRUE),
     type = "l", xlab = "sigmaR", ylab = "dNLL")


# 2-D cross-profile: M1 across sex
prof_M_sex_25_1 <- profile(
  fitted = mod_25_1,
  param    = "M1",
  slots    = list(c(1, 1, 1), c(1, 2, 1)),  # males and females
  values   = list(seq(0.10, 0.40, by = 0.05),
                  seq(0.10, 0.40, by = 0.05))
)

# Model Comparison ----

plot_selectivity(mod_25_0)
plot_selectivity(mod_25_1)
plot_f(list(mod_25_0,mod_25_1), model_names = c("Model 25.0", "Model 25.1"), file = "Results/Model comparison", legend.pos = "bottomleft")
plot_biomass(list(mod_25_0,mod_25_1), model_names = c("Model 25.0", "Model 25.1"), incl_proj = TRUE, add_ci = TRUE, file = "Results/Model comparison", legend.pos = "bottomleft")
plot_ssb(list(mod_25_0,mod_25_1), model_names = c("Model 25.0", "Model 25.1"), incl_proj = TRUE, add_ci = TRUE, file = "Results/Model comparison", legend.pos = "bottomleft")
plot_recruitment(list(mod_25_0,mod_25_1), model_names = c("Model 25.0", "Model 25.1"), incl_proj = TRUE, add_ci = TRUE, file = "Results/Model comparison")
plot_index(list(mod_25_0,mod_25_1), model_names = c("Model 25.0", "Model 25.1"), file = "Results/Model comparison")

mod_25_0$sdrep
mod_25_0$opt
mod_25_0$quantities

# Projection Models ----

# * Fit Base 25.0 -----
# Fixed M
mod_25_0 <- Rceattle::fit_mod(data_list = data_2026,
                              estimateMode = 2, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1)


# * Fit Alternative 25.1 -----
# Estimate sex-specific M
mod_25_1 <- Rceattle::fit_mod(data_list = data_2026,
                              estimateMode = 2, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1,
                              M1Fun = build_M1(M1_model = 2))



# Research Models ----
# Cannibalism Multispecies Model



# ESP-DSEM Integrated Model ----

# SEM ----
#FIXME @kalei Need to update

# * IID sem ----
# - For model comparison
atf_iid_sem = "
  # source                  link  target,                    lag param_name        start
  # ------------------------------------------------------------------------------------
  # --- AR1 processes (one per environmental variable) ---
  EnvIndex1           ->  EnvIndex1,            1,  EnvIndex1_AR1,            0
  EnvIndex2           ->  EnvIndex2,            1,  EnvIndex2_AR2,            0
  EnvIndex3           ->  EnvIndex3,            1,  EnvIndex3_AR3,            0

  # --- Recruitment variance ---
  recdevs1 <-> recdevs1,                        0,  sigmaR1,                  1
"

# * Full SEM ----
atf_sem = "
  # source                  link  target,                    lag param_name        start
  # ------------------------------------------------------------------------------------
  # --- AR1 processes (one per environmental variable) ---
  EnvIndex1           ->  EnvIndex1,            1,  EnvIndex1_AR1,            0
  EnvIndex2           ->  EnvIndex2,            1,  EnvIndex2_AR2,            0
  EnvIndex3           ->  EnvIndex3,            1,  EnvIndex3_AR3,            0

  # --- Recruitment ---
  EnvIndex1           ->  recdevs1,             1,  EnvIndex1_to_R,           0
  EnvIndex2           ->  recdevs1,             1,  EnvIndex2_to_R,           0
  EnvIndex3           ->  recdevs1,             1,  EnvIndex3_to_R,           0

  # --- Recruitment variance ---
  recdevs1 <-> recdevs1,                        0,  sigmaR1,                  1
"

# Run DSEM model ----
# * IID SEM ----
atf_iid_mod <- Rceattle::fit_mod(data_list = data_2026,
                                  estimateMode = 0, # Estimate
                                  random_rec = TRUE,
                                  dsem = build_DSEM(
                                    sem = atf_iid_sem,
                                    family = "fixed",
                                    sigmaR_prior_sd = 0.5 # Need sd prior to get it to converge if missing a lot of years of env data (assumes mean as initial value)
                                  ),
                                  msmMode = 0, # Single species mode
                                  initMode = 1,
                                  fit_control = fit_control(
                                    verbose = 1,
                                    phase = TRUE))

summary(atf_iid_mod)
AIC(atf_iid_mod)

# * Full SEM ----
atf_dsem_mod <- Rceattle::fit_mod(data_list = data_2026,
                                  estimateMode = 0, # Estimate
                                  random_rec = TRUE,
                                  dsem = build_DSEM(
                                    sem = atf_sem,
                                    family = "fixed",
                                    sigmaR_prior_sd = 0.5 # Need sd prior to get it to converge if missing a lot of years of env data (assumes mean as initial value)
                                  ),
                                  msmMode = 0, # Single species mode
                                  initMode = 1,
                                  fit_control = fit_control(
                                    verbose = 1,
                                    phase = TRUE))

summary(atf_dsem_mod)
AIC(atf_dsem_mod)
