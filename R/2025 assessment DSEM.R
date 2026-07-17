# Libraries ----
# remotes::install_version("dsem", version = "3.0.0") # Need this version
# remotes::install_github("grantdadams/Rceattle@dev-DSEM")
library(Rceattle) # dev-DSEM
library(dplyr)

# Data ----
data_2025_c_la_wa_ae <- Rceattle::read_data( file = "Data/GOA_25.1.1_arrowtooth_single_species_1977-2024_new_data_clean_LA_WA_AE.xlsx")
data_2025_c_la_wa_ae$estDynamics = 0
data_2025_c_la_wa_ae$index_data$Log_sd <- data_2025_c_la_wa_ae$index_data$Log_sd/data_2025_c_la_wa_ae$index_data$Observation

#FIXME update to 2026 data


# Add DSEM environmental data to Rceattle data object ----
# @kalei - TODO
env_data <- data.frame(Year = 1977:2025, EnvIndex1 = rnorm(49), EnvIndex2 = rnorm(49), EnvIndex3 = rnorm(49)) # Random (will be commented out)

data_2025_c_la_wa_ae$env_data <- data_2025_c_la_wa_ae$env_data %>%
  dplyr::full_join(env_data)

plot_data(data_2025_c_la_wa_ae)


# SEM ----
#FIXME Need to update

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
atf_iid_mod <- Rceattle::fit_mod(data_list = data_2025_c_la_wa_ae,
                                  estimateMode = 0, # Estimate
                                  random_rec = TRUE,
                                  dsem = build_DSEM(
                                    sem = atf_iid_sem,
                                    family = "fixed",
                                    sigmaR_prior_sd = 0.5 # Need sd prior to get it to converge
                                  ),
                                  msmMode = 0, # Single species mode
                                  initMode = 1,
                                  fit_control = fit_control(
                                    verbose = 1,
                                    phase = TRUE))

summary(atf_iid_mod)
AIC(atf_iid_mod)

# * Full SEM ----
atf_dsem_mod <- Rceattle::fit_mod(data_list = data_2025_c_la_wa_ae,
                                  estimateMode = 0, # Estimate
                                  random_rec = TRUE,
                                  dsem = build_DSEM(
                                    sem = atf_sem,
                                    family = "fixed",
                                    sigmaR_prior_sd = 0.5 # Need sd prior to get it to converge
                                  ),
                                  msmMode = 0, # Single species mode
                                  initMode = 1,
                                  fit_control = fit_control(
                                    verbose = 1,
                                    phase = TRUE))

summary(atf_dsem_mod)
AIC(atf_dsem_mod)
