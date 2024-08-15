# Code used to fit the 2023 arrowtooth flounder assessment in CEATTLE
# uses the "dev_srr" version of Rceattle

library(Rceattle)
library(readxl)
library(dplyr)
library(TMB)

# Load data ----
mydata_atf <- Rceattle::read_data( file = "Data/GOA_23.1.1_arrowtooth_single_species_1977-2023.xlsx")
mydata_atf$estDynamics = 0
mydata_atf$srv_biom$Log_sd <- mydata_atf$srv_biom$Log_sd/mydata_atf$srv_biom$Observation
mydata_atf$fleet_control$proj_F_prop <- c(1,1,1)

# * Adjust ageing error ----
agerr_mat <- as.matrix(mydata_atf$age_error[1:21,3:23])
mydata_atf$age_error[1:21,3:23] <- diag(1, 21) # Remove ageing error from predicted

mydata_atf$comp_data[1:16,9:29] <- as.matrix(mydata_atf$comp_data[1:16,9:29]) %*% agerr_mat # Females
mydata_atf$comp_data[1:16,30:50] <- as.matrix(mydata_atf$comp_data[1:16,30:50]) %*% agerr_mat # Males

for(i in 1:nrow(mydata_atf$comp_data)){
  if(mydata_atf$comp_data$Age0_Length1[i] == 0 & mydata_atf$comp_data$Fleet_code[i] == 1){
    mydata_atf$comp_data[i,9:50] <- mydata_atf$comp_data[i,9:50]/sum(mydata_atf$comp_data[i,9:50])
  }
}


# Single-species models ----
# - Fit single-species models with tier-3 HCR using CEATTLE likelihoods
bridging_model_3 <- Rceattle::fit_mod(data_list = mydata_atf,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      msmMode = 0, # Single species mode
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                      DynamicHCR = FALSE, # equilibrium reference points
                                                      FsprTarget = 0.4, # F40%
                                                      FsprLimit = 0.35, # F35%
                                                      Plimit = 0,
                                                      Alpha = 0.05))

# - No fishing
bridging_model_3_nof <- Rceattle::fit_mod(data_list = mydata_atf,
                                          inits = NULL, # Initial parameters = 0
                                          file = NULL, # Don't save
                                          estimateMode = 0, # Estimate
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 0, # Single species mode
                                          verbose = 1,
                                          phase = "default",
                                          initMode = 1)


# * Model w/ M ----
# - Fit single-species models with tier-3 HCR using CEATTLE likelihoods, but estimate Ms
bridging_model_4 <- Rceattle::fit_mod(data_list = mydata_atf,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      msmMode = 0, # Single species mode
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                      DynamicHCR = FALSE, # equilibrium reference points
                                                      FsprTarget = 0.4, # F40%
                                                      FsprLimit = 0.35, # F35%
                                                      Plimit = 0,
                                                      Alpha = 0.05),
                                      M1Fun = build_M1(M1_model = 2) # Estimate M
)



# Multi-species model ----
# (cannibalism)
ms_model <- Rceattle::fit_mod(data_list = mydata_atf,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              verbose = 1,
                              phase = "default",
                              suit_meanyr = 2018,
                              initMode = 1,
                              msmMode = 1, # Multi-species model
                              M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)




# Compare models ----
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)

SAFE2023_comp <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 4)

SAFE2023_mod <- bridging_model_3
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000

# -- Comp data
comp_temp <- SAFE2023_mod$data_list$comp_data %>%
  dplyr::select(-paste0("Comp_", 1:117))

comp_temp <- comp_temp %>%
  left_join(SAFE2023_comp)

SAFE2023_mod$quantities$comp_hat <- comp_temp %>%
  dplyr::select(paste0("Comp_", 1:117)) %>%
  as.matrix()



# Plot final ----
MPcols <- rev(oce::oce.colorsViridis(3))
line_col <- c("grey60", MPcols[2:3])
model_list <- list(SAFE2023_mod, bridging_model_3, ms_model)
model_names = c("ADMB", "TMB single-spp", "TMB multi-spp")

plot_biomass(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_ssb(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_recruitment(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

line_col <- MPcols[2:3]
model_list <- list(bridging_model_3, ms_model)
model_names = c("TMB single-spp", "TMB multi-spp")
plot_catch(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_b_eaten(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

# * Plot diagnostics ----
plot_pearson(ms_model, file = "Results/Figures/Diagnostics/Final_ms_")
plot_pearson(bridging_model_3, file = "Results/Figures/Diagnostics/Final_ss_")
plot_pearson(SAFE2023_mod, file = "Results/Figures/Diagnostics/Final_ADMB_")

plot_pearson(ms_model, file = "Results/Figures/Diagnostics/Final_ms_6_", max_z = 6)
plot_pearson(bridging_model_3, file = "Results/Figures/Diagnostics/Final_ss_6_", max_z = 6)

plot_pearson(ms_model, file = "Results/Figures/Diagnostics/Final_ms_51_", max_z = 51)
plot_pearson(bridging_model_3, file = "Results/Figures/Diagnostics/Final_ss_51_", max_z = 51)




