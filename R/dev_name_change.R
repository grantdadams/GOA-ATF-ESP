# Code used to fit the 2023 arrowtooth flounder assessment in CEATTLE
# uses the "dev_srr" version of Rceattle

library(Rceattle)
library(readxl)
library(dplyr)
library(TMB)

# Load data ----
mydata_atf <- Rceattle::read_data( file = "Data/GOA_23.1.1_arrowtooth_single_species_1977-2023.xlsx")
mydata_atf$estDynamics = 0
mydata_atf$index_data$Log_sd <- mydata_atf$index_data$Log_sd/mydata_atf$index_data$Observation
mydata_atf$fleet_control$proj_F_prop <- c(1,1,1)
mydata_atf$fleet_control$Index_sd_prior <- mydata_atf$fleet_control$Survey_sd_prior
mydata_atf$fleet_control$Estimate_index_sd <- mydata_atf$fleet_control$Estimate_survey_sd
mydata_atf$fleet_control$Comp_loglike <- -1
mydata_atf$fleet_control$Age_max_selected <- NA

inits <- build_params(mydata_atf)
for(i in 1:length(mod_objects$estimated_params)){
  parname = names(mod_objects$estimated_params)[i]
  if(parname %in% names(inits)){
    inits[parname] = mod_objects$estimated_params[i]
  }
}

inits$ln_F = mod_objects$estimated_params$ln_mean_F + mod_objects$estimated_params$F_dev

# Single-species models ----
# - Fit single-species models and no fishing
# * Fix M ----
ceattle_ss <- Rceattle::fit_mod(data_list = mydata_atf,
                                inits = inits, # Initial parameters = 0
                                file = NULL , #"Models/ss", # Don't save
                                estimateMode = 4, # Estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single species mode
                                verbose = 1,
                                phase = FALSE,
                                initMode = 2)
plot_biomass(list(mod_objects, ceattle_ss))
ceattle_ss$quantities$jnll_comp
mod_objects$quantities$jnll_comp
