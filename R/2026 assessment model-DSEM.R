# Libraries ----
remotes::install_version("dsem", version = "3.0.0") # Need this version
remotes::install_github("grantdadams/Rceattle@dsem-v5-integration")
library(Rceattle) # dsem-v5-integration
library(dsem)
library(dplyr)
library(ggplot2)
library(here)

# Globals ----
yr<-2026 #year for output, usually current year
dir.create(here(yr, "results"))

# Data ----
# 2025 Model
data_2025 <- Rceattle::read_data( file = "2026/data/GOA_arrowtooth_2025.xlsx")
data_2025$estDynamics = 0
data_2025$index_data$Log_sd <- data_2025$index_data$Log_sd/data_2025$index_data$Observation

# Excel is storing the environmental columns as text, so they arrive as
# character and TMB rejects them. Turn them back into numbers.
env_columns <- setdiff(names(data_2025$env_data), "Year")
data_2025$env_data[env_columns] <- lapply(data_2025$env_data[env_columns], as.numeric)

# 2026 Model
data_2026 <- Rceattle::read_data( file = "2026/data/GOA_arrowtooth_2026.xlsx")
data_2026$estDynamics = 0
data_2026$index_data$Log_sd <- data_2026$index_data$Log_sd/data_2026$index_data$Observation
# Return to previous comp weights
#data_2026$fleet_control$Comp_weights <- c(1,1,0.25)

# Excel is storing the environmental columns as text, so they arrive as
# character and TMB rejects them. Turn them back into numbers.
env_columns <- setdiff(names(data_2026$env_data), "Year")
data_2026$env_data[env_columns] <- lapply(data_2026$env_data[env_columns], as.numeric)

plot_data(data_2025)
plot_data(data_2026)

# SAFE (Hindcast) Model ----
# Base 25.0 = single-species TMB based Rceattle model with data improvements that fixes sex-specific M (females = 0.2 and males = 0.35) and treats annual recruitment as random effects.
# Alternative 25.1 = Model 25.0 but estimates sex-specific M instead of fixing M.

# * Fit Base 25.0 -----
# Fixed M Old
mod_25_old <- Rceattle::fit_mod(data_list = data_2025,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1)
summary(mod_25_old)
convergence_diagnostics(mod_25_old)

# Fixed M
mod_25_0 <- Rceattle::fit_mod(data_list = data_2026,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = TRUE),
                              initMode = 1)

summary(mod_25_0)
convergence_diagnostics(mod_25_0)

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

summary(mod_25_1)
convergence_diagnostics(mod_25_1)

# * Fit Alternative 25.2 -----
# - Cannibalism
mod_25_2 <- Rceattle::fit_mod(data_list = data_2026,
                              inits = mod_25_1$estimated_params,
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # Random recruitment
                              msmMode = 1, # Multi species mode
                              fit_control = fit_control(
                                verbose = 1,
                                phase = FALSE), # No phasing
                              initMode = 1,
                              M1Fun = build_M1(M1_model = 2))

summary(mod_25_2)
convergence_diagnostics(mod_25_2)

# Make a list of models to run for rest of script
mod_list<-c("mod_25_old","mod_25_0","mod_25_1","mod_25_2")

# Set i for your loop
i <- 2
print(mod_list[i])
dir.create(here(yr, "results", mod_list[i]))

# Fetch the actual model object using get()
mod <- get(mod_list[i])

# Diagmod# Diagnostics ----
# * Summaries -----
summary(mod)
convergence_diagnostics(mod)

# * Plots ----
plot_index(mod)
plot_index(mod, log = TRUE)
plot_indexresidual(mod)
plot_biomass(mod, incl_proj = TRUE, add_ci = TRUE)
plot_ssb(mod, incl_proj = TRUE, add_ci = TRUE)
plot_recruitment(mod, incl_proj = TRUE, add_ci = TRUE)
plot_comp(mod)
plot_catch(mod)
plot_selectivity(mod)

# * OSAs ----
osa <- osa_residuals(mod)
head(osa)

# Statistical diagnostics (Stewart & Monnahan 2025): SDNR and lower/upper tail
osa_diags<-osa_diagnostics(osa)
readr::write_csv(osa_diags, here::here(yr, "results", mod_list[i],"osa_diagnostics.csv"))

# Q-Q plot (with SDNR / tail annotation) + residual-by-year
osa_plots<-plot(osa)

ggsave(filename = here::here(yr, "results", mod_list[i],"osa_aggregate.png"),
       plot = osa_plots$aggregate, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

ggsave(here::here(yr, "results", mod_list[i],"osa_composition.png"),
       plot = osa_plots$composition, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

# * Retrospectives ----
retro <- retrospective(Rceattle = mod, peels = 10)
mohns <- retro$mohns # Mohn's rho for each quantity
readr::write_csv(mohns, here::here(yr, "results", mod_list[i],"mohns_diagnostics.csv"))

# Plot historical trajectories across peels
retro_plot_biomass<-plot_biomass(retro$Rceattle_list, add_ci = TRUE)
retro_plot_ssb<-plot_ssb(retro$Rceattle_list, add_ci = TRUE)
retro_plot_recruitment<-plot_recruitment(retro$Rceattle_list, add_ci = TRUE)

ggsave(filename = here::here(yr, "results", mod_list[i],"biomass_retro.png"),
       plot = retro_plot_biomass, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)
ggsave(filename = here::here(yr, "results", mod_list[i],"ssb_retro.png"),
       plot = retro_plot_ssb, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)
ggsave(filename = here::here(yr, "results", mod_list[i],"recruitment_retro.png"),
       plot = retro_plot_recruitment, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)


# * Jitters ----
mod_jitters <- jitter(Rceattle = mod, njitter = 100, phase = TRUE)

# Histogram of NLL differences relative to the best run
jitter_data <- data.frame(
  nll_diff = log(mod_jitters$nll - min(mod_jitters$nll))
)

jitter_hist <- ggplot(data = jitter_data, aes(x = nll_diff)) +
  geom_histogram(fill = "grey70", color = "black", binwidth = 1, boundary = 0) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  labs(title = "Jitter NLL spread (log scale)",
       x = "log(NLL - min NLL)",
       y = "Frequency") +
  theme_bw() # Or whatever theme you prefer

ggsave(filename = here::here(yr, "results", mod_list[i],"jitter_hist.png"),
       plot = jitter_hist, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

# Overlay biomass trajectories — tight overlap indicates a stable optimum
biomass_jitters <- plot_biomass(mod_jitters$Rceattle_list) + theme(legend.position="none")

ggsave(filename = here::here(yr, "results", mod_list[i],"biomass_jitters.png"),
       plot = biomass_jitters, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

# * Self-test ----
mod_sims <- self_test(mod, nsim = 100, start = "estimated")

# Number of simulations that converged (non-converged runs are dropped)
num_sims<- length(mod_sims)
num_sims_df <- data.frame(converged_simulations = num_sims)
readr::write_csv(num_sims_df, here::here(yr, "results", mod_list[i],"num_sims_converged.csv"))

# Overlay biomass / SSB trajectories across simulations — the original fit's
# trajectory should sit inside the spread of the refits.
biomass_sims <- plot_biomass(c(mod_sims, list(mod)), line_col = c(rep("grey", 100), 1)) + theme(legend.position="none")

ggsave(filename = here::here(yr, "results", mod_list[i],"biomass_sims.png"),
       plot = biomass_sims, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

# * Likelihood profiles ----
# sigmaR
prof_sigmaR <- profile(
  fitted = mod,
  param    = "sigmaR",
  slots    = list(1),
  values   = list(seq(0.1, 1.5, by = 0.05))
)
prof_sigmaR
minimum?

plot(prof_sigmaR$grid$slot_1,
     prof_sigmaR$nll - min(prof_sigmaR$nll, na.rm = TRUE),
     type = "l", xlab = "sigmaR", ylab = "dNLL")

plot_profile(prof_sigmaR, xlab = "sigmaR")


# 2-D cross-profile: M1 across sex
prof_M_sex <- profile(
  fitted = mod,
  param    = "M1",
  slots    = list(c(1, 1, 1), c(1, 2, 1)),  # males and females
  values   = list(seq(0.10, 0.40, by = 0.05),
                  seq(0.10, 0.40, by = 0.05))
)

#not sure about this one
plot(prof_M_sex$grid$slot_1,
     prof_M_sex$nll - min(prof_M_sex$nll, na.rm = TRUE),
     type = "l", xlab = "M", ylab = "dNLL")

# Model Comparison ----

f_compare <- plot_f(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"))

ggsave(filename = here::here(yr, "results","f_comparison.png"),
       plot = f_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

biom_compare <- plot_biomass(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), incl_proj = TRUE, add_ci = TRUE)

ggsave(filename = here::here(yr, "results", "biomass_comparison.png"),
       plot = biom_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

ssb_compare <- plot_ssb(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), incl_proj = TRUE, add_ci = TRUE, file = "Results/Model comparison", legend.pos = "bottomleft")

ggsave(filename = here::here(yr, "results", "ssb_comparison.png"),
       plot = ssb_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

recruit_compare <- plot_recruitment(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), incl_proj = TRUE, add_ci = TRUE, file = "Results/Model comparison")

ggsave(filename = here::here(yr, "results", "recruit_comparison.png"),
       plot = recruit_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

index_compare <- plot_index(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), file = "Results/Model comparison")

ggsave(filename = here::here(yr, "results", "index_comparison.png"),
       plot = index_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

sel_compare <- plot_selectivity(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), file = "Results/Model comparison")

ggsave(filename = here::here(yr, "results", "selectivity_comparison.png"),
       plot = sel_compare, units = 'in', bg = 'white', height = 8, width = 11, dpi = 300)

#plot_mortality(list(mod_25_old,mod_25_0,mod_25_1,mod_25_2), model_names = c("Model 25.0 Fix M","Model 25.0 Fix M New Data", "Model 25.1 Estimate M", "Model 25.1 Multispecies"), file = "Results/Model comparison")

# Projection Models ----


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
