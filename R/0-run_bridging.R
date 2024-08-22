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
mydata_atf_ageerror <- mydata_atf # Save data with "correct" ageing error set-up

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
# * Model 2 ----
# - Fix n-at-age and parameters to 2023 SAFE-multinomial selectivity values
mydata_atf_fixed <- mydata_atf
mydata_atf_fixed$wt[,6:ncol(mydata_atf_fixed$wt)] <- mydata_atf_fixed$wt[,6:ncol(mydata_atf_fixed$wt)] * 0.001 # Convert wt to tonnes
mydata_atf_fixed$estDynamics = 1 # Do not estimate
mydata_atf_fixed$NByageFixed[,5:25] <- mydata_atf_fixed$NByageFixed[,5:25] / 1000 # Numbers at age in CEATTLE is in 1,000's

# Make a parameters object to fill in with SAFE model estimates
mydata_atf_fixed$srr_prior_mean <- 9
mydata_atf_fixed$initMode <- 1
inits <- build_params(mydata_atf_fixed)

# Fishery selectivity - Non-parametric
inits$sel_coff[3,1,1:19] <- c(-4.51694290833, -3.70259980774, -2.90228183839, -2.14347044956, -1.45768792710, -0.868760941137, -0.390282722764, -0.0312224931896, 0.201213866945, 0.317107821839, 0.359031061254, 0.408570118597, 0.442205388565, 0.424347346174, 0.398254677665, 0.374330633827, 0.354849669770, 0.340369152572, 0.330663722726
)
inits$sel_coff[3,2,1:19] <- c(-4.56241988161, -3.65989030018, -2.77330837325, -1.93535643447, -1.18687503496, -0.566563860235, -0.100731040158, 0.205993405883, 0.370842111120, 0.437969489649, 0.441362712492, 0.421323478258, 0.387411477262, 0.349174984741, 0.310334510039, 0.273480355699, 0.240778867356, 0.214146429961, 0.196024113277
)

# Fishing mortality
inits$ln_mean_F[3] <- -4.10928686426
inits$F_dev[3,1:47] <- c(-0.177978417724, -0.285519674842, -0.379484797513, -0.333640700870, -0.377434679601, -0.842175623138, -0.528428968796, -1.13519569584, -1.95144949500, -2.19209537185, -0.797343656032, -0.771884491581, -1.48027433576, -0.422772922152, 0.349095106282, 0.535988212551, 0.351552383253, 0.494724755078, 0.236107198591, 0.419042326678, 0.0859306576161, -0.156970897161, 0.0606876147991, 0.474036993652, 0.289689815949, 0.354974406565, 0.704866515809, 0.0312545284044, 0.232903215652, 0.527265591631, 0.418580218424, 0.543538896733, 0.382034931302, 0.373121242209, 0.624386187662, 0.241829259196, 0.316947376678, 0.887487145387, 0.297833007717, 0.383252388915, 0.737419566196, 0.423042648172, 0.724895475427, 0.605539138364, -0.119803917246, 0.0458048197550, -0.201377979530
)

# Recruitment (recruitment in the SAFE model is mean_rec * exp(rec_devs) for each sex, so needs to be multiplied by 2 for CEATTLE)
rec_devs <- c(-0.193891196100, -0.102380215922, -0.119560997407, -0.0882990076250, -0.129117933501, -0.146055003125, -0.141335283066, -0.127634637257, -0.289192661784, -0.286125646152, -0.267628324752, -0.284329728819, -0.183537052124, -0.0160034490241, -0.147406775790, -0.194670626260, -0.207234256988, -0.271954367590, -0.288244675028, -0.291678767211, -0.171145345964, -0.105222115041, -0.00918346547033, -0.130785774246, -0.267456683472, -0.309763287996, -0.189087371857, 0.134817799530, 0.515084336470, 0.459390984524, 0.420757842134, 0.461291534029, 0.380781306864, 0.245113183221, 0.285292200803, 0.281316964420, 0.0434715972116, 0.0297848013545, 0.102022090218, 0.0750247855189, 0.265584065442, 0.518193299520, 0.595705388860, 0.841138163072, 0.410989072363, 0.306116057813, 0.216589245502, 0.251695966675, 0.287792199671, 0.210945726508, -0.0932783021390, -0.132174199071, -0.420630472653, -0.443081298889, -0.198021413041, -0.0341137832990, -0.0169380428935, -0.220630297025, -0.343067418290, -0.209164725275, -0.118916372151, 0.183713058619, -0.249015985084, -0.0752681408134, -0.00786989171491, -0.00333465900830, 0.00181798057919)
mean_rec <- 19.6917382653

inits$rec_pars[1,1] =  mean(log(exp(mean_rec+rec_devs)*2))
inits$rec_dev[1,1:47] <- rec_devs[21:67]
inits$init_dev[1,1:20] <- rev(rec_devs[1:20])

# Survey 1 (bottom trawl) - Logistic, q = 1
inits$ln_sel_slp[1,1,1] <- log(1.71587657721) # Fem
inits$ln_sel_slp[1,1,2] <- log(1.70293207654) # Males
inits$sel_inf[1,1,1] <- 2.78883897273 # Fem
inits$sel_inf[1,1,2] <- 2.97249336208 # Males

# Survey 2 (survey 1 length comp) - Logistic, q = 1
inits$ln_sel_slp[1,2,1] <- log(1.71587657721) # Fem
inits$ln_sel_slp[1,2,2] <- log(1.70293207654) # Males
inits$sel_inf[1,2,1] <- 2.78883897273 # Fem
inits$sel_inf[1,2,2] <- 2.97249336208 # Males

mydata_atf_fixed$estDynamics <- 0
bridging_model_2 <- Rceattle::fit_mod(data_list = mydata_atf_fixed,
                                       inits = inits, # Initial parameters = 0
                                       file = NULL, # Don't save
                                       estimateMode = 4, # Estimate
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 0, # Single species mode
                                       verbose = 1,
                                       phase = NULL,
                                       initMode = 1)


# * Model 3 ----
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
                                      # HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      #                 DynamicHCR = FALSE, # equilibrium reference points
                                      #                 FsprTarget = 0.4, # F40%
                                      #                 FsprLimit = 0.35, # F35%
                                      #                 Plimit = 0,
                                      #                 Alpha = 0.05)
                                      )

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


# * Model 4 ----
# - Fit single-species models with tier-3 HCR using CEATTLE likelihoods
# - AND correct ageing error
bridging_model_4 <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      msmMode = 0, # Single species mode
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      # HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      #                 DynamicHCR = FALSE, # equilibrium reference points
                                      #                 FsprTarget = 0.4, # F40%
                                      #                 FsprLimit = 0.35, # F35%
                                      #                 Plimit = 0,
                                      #                 Alpha = 0.05)
                                      )

# -- No fishing
bridging_model_4_nof <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                                          inits = NULL, # Initial parameters = 0
                                          file = NULL, # Don't save
                                          estimateMode = 0, # Estimate
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 0, # Single species mode
                                          verbose = 1,
                                          phase = "default",
                                          initMode = 1)

# * Model 5 ----
# - Model 4 with rec as random effects
bridging_model_5 <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                                          inits = bridging_model_4_nof$estimated_params, # Initial parameters = 0
                                          file = NULL, # Don't save
                                          estimateMode = 0, # Estimate
                                          random_rec = TRUE, # No random recruitment
                                          msmMode = 0, # Single species mode
                                          verbose = 1,
                                          phase = NULL,
                                          initMode = 1)


# * Model 6 ----
# - Fit single-species models with tier-3 HCR using CEATTLE likelihoods, but estimate Ms
bridging_model_6 <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      msmMode = 0, # Single species mode
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      # HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      #                 DynamicHCR = FALSE, # equilibrium reference points
                                      #                 FsprTarget = 0.4, # F40%
                                      #                 FsprLimit = 0.35, # F35%
                                      #                 Plimit = 0,
                                      #                 Alpha = 0.05),
                                      M1Fun = build_M1(M1_model = 2) # Estimate M
)



# Multi-species model ----
# (cannibalism)
ms_model <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                              inits = bridging_model_6$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = FALSE, # No random recruitment
                              verbose = 1,
                              phase = NULL,
                              suit_meanyr = 2015,
                              initMode = 1,
                              msmMode = 1, # Multi-species model
                              M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)

ms_model_RE <- Rceattle::fit_mod(data_list = mydata_atf_ageerror,
                              inits = ms_model$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 0, # Estimate
                              random_rec = TRUE, # No random recruitment
                              verbose = 1,
                              phase = NULL,
                              suit_meanyr = 2015,
                              initMode = 1,
                              msmMode = 1, # Multi-species model
                              M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)


# # * Input catch from single-species Tier-3 and project ----
# # ** Model specifications ----
# hind_yrs <- mydata_atf$styr: mydata_atf$endyr
# hind_nyrs <- length(hind_yrs)
# proj_yrs <- (mydata_atf$endyr + 1) : mydata_atf$projyr
# proj_nyrs <- length(proj_yrs)
# nflts = nrow(ms_model$data_list$fleet_control)
# new_years <- proj_yrs
#
#
# # ** Update rec devs (mean recruitment, not R0)
# ms_model_proj <- ms_model
# ms_model_proj$data_list$endyr <- 2050
# rec_dev <- log(mean(ms_model_proj$quantities$R[1,1:hind_nyrs]))  - log(ms_model_proj$quantities$R0[1]) #
# ms_model_proj$estimated_params$rec_dev[1,proj_yrs - ms_model_proj$data_list$styr + 1] <- replace(
#   ms_model_proj$estimated_params$rec_dev[1,proj_yrs - ms_model_proj$data_list$styr + 1],
#   values =  rec_dev)
#
#
# # -- wt
# proj_wt <- ms_model_proj$data_list$wt %>%
#   group_by(Wt_index , Sex) %>%
#   slice(rep(n(),  proj_nyrs)) %>%
#   mutate(Year = proj_yrs)
# ms_model_proj$data_list$wt  <- rbind(ms_model_proj$data_list$wt, proj_wt)
# ms_model_proj$data_list$wt <- dplyr::arrange(ms_model_proj$data_list$wt, Wt_index, Year)
#
# # -- Pyrs
# proj_Pyrs <- ms_model_proj$data_list$Pyrs %>%
#   group_by(Species, Sex) %>%
#   slice(rep(n(),  proj_nyrs)) %>%
#   mutate(Year = proj_yrs)
# ms_model_proj$data_list$Pyrs  <- rbind(ms_model_proj$data_list$Pyrs, proj_Pyrs)
# ms_model_proj$data_list$Pyrs <- dplyr::arrange(ms_model_proj$data_list$Pyrs, Species, Year)
#
#
# # ** GET RECOMMENDED TAC FROM Single-species model ----
# # - Get projected catch data from single-species model
# new_catch_data <- bridging_model_4$data_list$fsh_biom #FIXME: final single-species model
# dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
# new_catch_data$Catch[dat_fill_ind] <- bridging_model_4$quantities$fsh_bio_hat[dat_fill_ind]
# ms_model_proj$data_list$fsh_biom <- new_catch_data
#
#
# # ** Update parameters ----
# # -- F_dev
# ms_model_proj$estimated_params$F_dev <- cbind(ms_model_proj$estimated_params$F_dev, matrix(0, nrow= nrow(ms_model_proj$estimated_params$F_dev), ncol = length(new_years)))
#
# # -- Time-varing survey catachbilitiy - Assume last year - filled by columns
# ms_model_proj$estimated_params$ln_srv_q_dev <- cbind(ms_model_proj$estimated_params$ln_srv_q_dev, matrix(ms_model_proj$estimated_params$ln_srv_q_dev[,ncol(ms_model_proj$estimated_params$ln_srv_q_dev)], nrow= nrow(ms_model_proj$estimated_params$ln_srv_q_dev), ncol = length(new_years)))
#
# # -- Time-varing selectivity - Assume last year - filled by columns
# ln_sel_slp_dev = array(0, dim = c(2, nflts, 2, hind_nyrs + length(new_years)))  # selectivity deviations paramaters for logistic
# sel_inf_dev = array(0, dim = c(2, nflts, 2, hind_nyrs + length(new_years)))  # selectivity deviations paramaters for logistic
# # sel_coff_dev = array(0, dim = c(nflts, 2, nselages_om, hind_nyrs + length(new_years)))  # selectivity deviations paramaters for non-parameteric
#
# ln_sel_slp_dev[,,,1:hind_nyrs] <- ms_model_proj$estimated_params$ln_sel_slp_dev
# sel_inf_dev[,,,1:hind_nyrs] <- ms_model_proj$estimated_params$sel_inf_dev
# # sel_coff_dev[,,,1:hind_nyrs] <- ms_model_proj$estimated_params$# sel_coff_dev
#
# ln_sel_slp_dev[,,,(hind_nyrs + 1):(hind_nyrs + length(new_years))] <- ln_sel_slp_dev[,,,hind_nyrs]
# sel_inf_dev[,,,(hind_nyrs + 1):(hind_nyrs + length(new_years))] <- sel_inf_dev[,,,hind_nyrs]
# # sel_coff_dev[,,,(hind_nyrs + 1):(hind_nyrs + length(new_years))] <- # sel_coff_dev[,,,hind_nyrs]
#
# ms_model_proj$estimated_params$ln_sel_slp_dev <- ln_sel_slp_dev
# ms_model_proj$estimated_params$sel_inf_dev <- sel_inf_dev
# # ms_model_proj$estimated_params$# sel_coff_dev <- # sel_coff_dev
#
#
# # ** Update map ----
# # - (Only new parameter we are estimating is the F_dev of the new years)
# ms_model_proj$map <- build_map(
#   data_list = ms_model_proj$data_list,
#   params = ms_model_proj$estimated_params,
#   debug = TRUE,
#   random_rec = ms_model_proj$data_list$random_rec)
# ms_model_proj$map$mapFactor$dummy <- as.factor(NA); ms_model_proj$map$mapList$dummy <- NA
#
#
# # -- Estimate terminal F for catch
# new_f_yrs <- (ncol(ms_model_proj$map$mapList$F_dev) - length(new_years) + 1) : ncol(ms_model_proj$map$mapList$F_dev) # - Years of new F
# f_fleets <- ms_model_proj$data_list$fleet_control$Fleet_code[which(ms_model_proj$data_list$fleet_control$Fleet_type == 1)] # Fleet rows for F
# ms_model_proj$map$mapList$F_dev[f_fleets,new_f_yrs] <- replace(ms_model_proj$map$mapList$F_dev[f_fleets,new_f_yrs], values = 1:length(ms_model_proj$map$mapList$F_dev[f_fleets,new_f_yrs]))
# ms_model_proj$map$mapFactor$F_dev <- factor(ms_model_proj$map$mapList$F_dev)
#
# # ** Estimate F ----
# ms_model_proj <- Rceattle::fit_mod(data_list = ms_model_proj$data_list,
#                                    inits = ms_model_proj$estimated_params,
#                                    map = ms_model_proj$map,
#                                    bounds = NULL,
#                                    file = NULL, # Don't save
#                                    estimateMode = 1, # Estimate
#                                    random_rec = FALSE, # No random recruitment
#                                    verbose = 1,
#                                    phase = NULL,
#                                    suit_meanyr = 2015,
#                                    initMode = 1,
#                                    msmMode = 1, # Multi-species model
#                                    getsd = FALSE,
#                                    M1Fun = build_M1(M1_model = 2) # Estimate residual M
# )
#
# plot_ssb(ms_model_proj, incl_proj = T)



# Compare models ----
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)

SAFE2023_mod <- bridging_model_3
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000


# - SAFE model with fixed multinomial and selectivity penalties set to the same for sexes (lines 470-472 of dat file)
SAFE2023multisel <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 3)
SAFE2023multisel_mod <- bridging_model_3
SAFE2023multisel_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023multisel$Biomass
SAFE2023multisel_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023multisel$SSB
SAFE2023multisel_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023multisel$Recruitment/1000
SAFE2023multisel_mod$quantities$fsh_bio_hat <- SAFE2023multisel$Catch

# CEATTLE with fixed parameters
# bridging_model_2$quantities$biomass <- bridging_model_2$quantities$biomass/1000
# bridging_model_2$quantities$biomassSSB <- bridging_model_2$quantities$biomassSSB/1000
bridging_model_2$quantities$R <- bridging_model_2$quantities$R/1000
# bridging_model_2$quantities$fsh_bio_hat <- bridging_model_2$quantities$fsh_bio_hat/1000
# bridging_model_2$quantities$srv_bio_hat <- bridging_model_2$quantities$srv_bio_hat/1000


# Plot bridging ----
model_list <- list(SAFE2023_mod, SAFE2023multisel_mod, bridging_model_2, bridging_model_3, bridging_model_4, bridging_model_5, ms_model, ms_model_RE)
model_names = c("Base - ADMB", "Model 1 - ADMB", "Model 2 - SS", "Model 3 - SS", "Model 4 - SS", "Model 5 - SS", "Model 6 - MS", "Model 7 - MS")

plot_biomass(model_list, model_names = NULL, file = "Results/Figures/Bridging_", width = 6, height = 3)
plot_ssb(model_list, model_names = NULL, file = "Results/Figures/Bridging_", width = 6, height = 3)
plot_recruitment(model_list, model_names = model_names, file = "Results/Figures/Bridging_", width = 6, height = 3)


# JNLL ----
source("R/Functions/likelihood comparisons.R", echo=TRUE)

model_list <- list(bridging_model_2, bridging_model_3, bridging_model_4, bridging_model_5, ms_model, ms_model_RE)
model_names = c("Model 2", "Model 3", "Model 4", "Model 5", "MS", "MS RE")

ss_ll <- get_atf_ll(bridging_model_2) %>%
  rename("Model 2" = Value)

for(i in 2:length(model_list)){
  ss_ll <- cbind(ss_ll,
                 get_atf_ll(model_list[[i]]) %>%
                   dplyr::select(Value)
                 )
}
colnames(ss_ll)[2:ncol(ss_ll)] <- model_names

write.csv(ss_ll, file = "Results/Bridging_jnll.csv")


