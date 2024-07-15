# Code used to fit the 2023 arrowtooth flounder assessment in CEATTLE
# uses the "dev" version of Rceattle

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

# From ADMB:
# Obs_srv_agecomp_fem
# 7.90666e-05 0.00313498 0.0379455 0.0652694 0.108601 0.130278 0.0905713 0.0523955 0.0289754 0.0175692 0.0132002 0.0113169 0.0109475 0.0115989 0.0108947 0.00733653 0.00320759 0.000875116 0.000146982 1.45321e-05 0
# 0.000760599 0.0298807 0.125971 0.10811 0.0813278 0.0790064 0.0683911 0.0589958 0.0451463 0.0323879 0.0223119 0.0157865 0.0121851 0.00876744 0.00583284 0.00395838 0.00247833 0.00119619 0.000408521 9.80151e-05 1.83925e-05
# 0.00107533 0.0307406 0.0552633 0.0727462 0.106164 0.100667 0.0758981 0.0623449 0.045969 0.0333777 0.0284382 0.0223733 0.0150648 0.0103278 0.00826743 0.00735492 0.00638912 0.00493648 0.00326254 0.00188223 0.00189843
# 0.00296019 0.046004 0.0638507 0.0617656 0.0668547 0.0760286 0.0806364 0.0816902 0.0745241 0.0523471 0.0331168 0.0234664 0.0166097 0.0106142 0.0063218 0.00376536 0.0023633 0.00148951 0.00085809 0.00044536 0.000401546


# Single-species models ----
# * Model 2a ----
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
inits$sel_coff[3,1,1:19] <- c(-5.13814131542, -4.06310922873, -3.02882477819, -2.09692893332, -1.33902047192, -0.762445774915, -0.301102813587, 0.0780622242856, 0.316610831580, 0.416969496168, 0.417508794464, 0.474223235902, 0.474223235928, 0.465280753866, 0.459188238711, 0.434205258057, 0.411996409676, 0.374425006662, 0.341596785184
)
inits$sel_coff[3,2,1:19] <- c(-5.48089549544, -4.22269598437, -2.99993774945, -1.88559271870, -0.962287534289, -0.289490042924, 0.118187123366, 0.316005878529, 0.437744989536, 0.448964602862, 0.454400002470, 0.450705441666, 0.420118603503, 0.385269934062, 0.360668710119, 0.331192360965, 0.304399771159, 0.282100759059, 0.265129831272
)

# Fishing mortality
inits$ln_mean_F[3] <- -4.13075343786
inits$F_dev[3,1:47] <- c(-0.0808287792366, -0.204706554083, -0.309917330618, -0.275603570708, -0.330809565602, -0.809494282234, -0.511742035688, -1.13268323793, -1.96104418150, -2.21162816731, -0.826052881746, -0.806797459240, -1.52140553557, -0.467207213704, 0.304221183497, 0.494104768845, 0.315559775223, 0.466578678181, 0.215177708587, 0.402822471394, 0.0745989206053, -0.164523277968, 0.0558179032618, 0.472208865479, 0.292168734909, 0.360970086620, 0.713959022873, 0.0436683210806, 0.248932951976, 0.547528589388, 0.442140860959, 0.568035920374, 0.405968176094, 0.395137256747, 0.642066417697, 0.255902527385, 0.326243717213, 0.893422845978, 0.301651309143, 0.383553329146, 0.732483837447, 0.411130450083, 0.702908815236, 0.572908383992, -0.162884487843, -0.00586266412835, -0.258680604317
)

# Recruitment (recruitment in the SAFE model is mean_rec * exp(rec_devs) for each sex, so needs to be multiplied by 2 for CEATTLE)
rec_devs <- c(-0.673372487520, -0.116252298134, -0.133223988152, -0.129888274265, -0.249476455968, -0.289065973436, -0.290813438052, -0.471304426971, -0.312284029440, -0.479091188985, -0.0744103146349, -0.522980830938, -0.00528333424906, -0.0808681374866, -0.0443303648467, -0.205927499849, -0.316502674685, -0.330242269043, -0.167119850975, -0.138387489507, -0.154680498622, -0.00939497093961, 0.0709556908117, -0.0272725461674, -0.210078063538, -0.241330199244, -0.101859593400, 0.213211748235, 0.553890618015, 0.418194967182, 0.362402871467, 0.438995361391, 0.383109733098, 0.239611704633, 0.318413474217, 0.293654145026, 0.0256320677247, 0.0202395208230, 0.0935169864899, 0.0632086077868, 0.227092193199, 0.457561246651, 0.538672444691, 0.796239505883, 0.386488212512, 0.311326388851, 0.229338281171, 0.291395300562, 0.340151008589, 0.263389698693, -0.0579317859696, -0.0854817459168, -0.413476514312, -0.440241525736, -0.155148618732, 0.0588964949996, 0.107627986227, -0.0927434556432, -0.242949701617, -0.102087606264, -0.0243900515327, 0.328816035412, -0.344789475147, -0.0867691297293, -0.0240277026218, 0.0163069396890, -0.00286072176250)
mean_rec <- 19.6756606490

inits$rec_pars[1,1] =  mean(log(exp(mean_rec+rec_devs)*2))
inits$rec_dev[1,1:47] <- rec_devs[21:67]
inits$init_dev[1,1:20] <- rev(rec_devs[1:20])

# Survey 1 - Logistic, q = 1
inits$ln_sel_slp[1,1,1] <- log(1.64676409761) # Fem
inits$ln_sel_slp[1,1,2] <- log(1.64315096373) # Males
inits$sel_inf[1,1,1] <- 2.83302919170 # Fem
inits$sel_inf[1,1,2] <- 3.01664650825 # Males

# Survey 2 - Logistic, q = 1
inits$ln_sel_slp[1,2,1] <- log(1.64676409761) # Fem
inits$ln_sel_slp[1,2,2] <- log(1.64315096373) # Males
inits$sel_inf[1,2,1] <- 2.83302919170 # Fem
inits$sel_inf[1,2,2] <- 3.01664650825 # Males

atf_fixed_numbers <- Rceattle::fit_mod(data_list = mydata_atf_fixed,
                                       inits = inits, # Initial parameters = 0
                                       file = NULL, # Don't save
                                       estimateMode = 4, # Estimate
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 0, # Single species mode
                                       verbose = 1,
                                       phase = "default",
                                       initMode = 1)

mydata_atf_fixed$estDynamics <- 0
bridging_model_2a <- Rceattle::fit_mod(data_list = mydata_atf_fixed,
                             inits = inits, # Initial parameters = 0
                             file = NULL, # Don't save
                             estimateMode = 4, # Estimate
                             random_rec = FALSE, # No random recruitment
                             msmMode = 0, # Single species mode
                             verbose = 1,
                             phase = "default",
                             initMode = 1,
                             TMBfilename = "src/ceattle_v01_10_atf_ll")

# Check likelihoods
bridging_model_2a$quantities$jnll_comp

check <- bridging_model_2a$data_list$comp_data
check[,9:ncol(check)] <- bridging_model_2a$quantities$true_age_comp_hat
sum(check$Sample_size[1:16] * (bridging_model_2a$data_list$comp_data[1:16,9:50] + 0.000001) * log((check[1:16,9:50] + 0.000001)))
bridging_model_2a$quantities$jnll_comp[3,]/bridging_model_2a$data_list$fleet_control$Comp_weights

# From ADMB:
# Pred_srv_agecomp_fem
# 0.0130388 0.0338264 0.067657 0.0877003 0.0960432 0.0885791 0.0670319 0.0472871 0.0390952 0.0308334 0.0212311 0.0174427 0.0157855 0.0150555 0.0117996 0.0103732 0.00505266 0.00647735 0.00353869 0.00342746 0.0148114
# 0.0123368 0.0463766 0.122149 0.10918 0.072712 0.0529302 0.044865 0.0440604 0.0396867 0.0298525 0.0210064 0.0173447 0.0136712 0.00941238 0.00773074 0.00699647 0.00667363 0.00523103 0.00459975 0.00224102 0.0125365
# 0.00913608 0.0374923 0.0911718 0.106113 0.102383 0.0981026 0.0572967 0.034182 0.0242661 0.0204108 0.0199658 0.0179464 0.0134877 0.00948712 0.00783056 0.00617262 0.00425055 0.00349216 0.00316187 0.00301748 0.0111335
# 0.00726879 0.0337822 0.0796151 0.0924068 0.097198 0.0857543 0.0649369 0.0557884 0.0517306 0.0297498 0.0175838 0.0124335 0.010433 0.0101932 0.00915924 0.0068861 0.00484779 0.0040059 0.00316285 0.00218177 0.0106922
#
#
# bridging_model_2a$quantities$sel[,1,,1]
# bridging_model_2a$quantities$sel[,2,,1]
#
# head(bridging_model_2a$data_list$NByageFixed[,5:25])
# head(t(bridging_model_2a$quantities$NByage[1,1,1:21,1:47]))
#
#
# t(bridging_model_2a$quantities$NByage[1,1,1:21,1:47])/1000 - bridging_model_2a$data_list$NByageFixed[1:47,5:25]


# * Model 2b ----
# - Fit single-species models with tier-3 HCR using ADMB likelihoods
bridging_model_2b <- Rceattle::fit_mod(data_list = mydata_atf,
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
                                       TMBfilename = "src/ceattle_v01_10_atf_ll")


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
                                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                      DynamicHCR = FALSE, # equilibrium reference points
                                                      FsprTarget = 0.4, # F40%
                                                      FsprLimit = 0.35, # F35%
                                                      Plimit = 0,
                                                      Alpha = 0.05))



# Multi-species model ----
# (cannibalism)
ms_model <- Rceattle::fit_mod(data_list = mydata_atf,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      msmMode = 1, # Multi-species model
                                      M1Fun = build_M1(M1_model = 1) # Estimate residual M
                              )


# * Input catch from single-species Tier-3
mydata_atf_ms <- mydata_atf
mydata_atf_ms$endyr <- 2050

proj_catch <- data.frame(Fleet_name = "ATF_total",
                         Fleet_code = 3, Species = 1,
                         Year = 2024:2050,
                         Month = 0, Selectivity_block = 1,
                         Catch = bridging_model_3$quantities$fsh_bio_hat[48:74],
                         Log_sd = 0.05)

mydata_atf_ms$fsh_biom <- rbind(
  mydata_atf_ms$fsh_biom,
  proj_catch)


bridging_model_3 <- Rceattle::fit_mod(data_list = mydata_atf_ms,
                                      inits = NULL, # Initial parameters = 0
                                      file = NULL, # Don't save
                                      estimateMode = 0, # Estimate
                                      random_rec = FALSE, # No random recruitment
                                      verbose = 1,
                                      phase = "default",
                                      initMode = 1,
                                      msmMode = 1, # Multi-species model
                                      M1Fun = build_M1(M1_model = 1) # Estimate residual M
                                      )


################################################
# Compare models
################################################
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)

SAFE2023_mod <- atf_base
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000


# - SAFE model with fixed multinomial and selectivity penalties set to the same for sexes (lines 470-472 of dat file)
SAFE2023multisel <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 3)

SAFE2023multisel_mod <- atf_base
SAFE2023multisel_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023multisel$Biomass
SAFE2023multisel_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023multisel$SSB
SAFE2023multisel_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023multisel$Recruitment/1000
SAFE2023multisel_mod$quantities$fsh_bio_hat <- SAFE2023multisel$Catch

# CEATTLE with fixed parameters
# bridging_model_2a$quantities$biomass <- bridging_model_2a$quantities$biomass/1000
# bridging_model_2a$quantities$biomassSSB <- bridging_model_2a$quantities$biomassSSB/1000
bridging_model_2a$quantities$R <- bridging_model_2a$quantities$R/1000
# bridging_model_2a$quantities$fsh_bio_hat <- bridging_model_2a$quantities$fsh_bio_hat/1000
# bridging_model_2a$quantities$srv_bio_hat <- bridging_model_2a$quantities$srv_bio_hat/1000

#
model_list <- list(atf_base, SAFE2023_mod, SAFE2023multi_mod, SAFE2023multisel_mod, bridging_model_2a)
model_names = c("CEATTLE", "SAFE", "SAFE-multi", "SAFE-multi-sel", "Fixed")


plot_biomass(list(atf_fixed_numbers, bridging_model_2a, SAFE2023multisel_mod))

plot_biomass(model_list, model_names = model_names)
plot_ssb(model_list, model_names = model_names)
plot_recruitment(model_list, model_names = model_names)
plot_catch(list(SAFE2023multisel_mod, bridging_model_2a, atf_base), model_names = c("SAFE", "CEATTLE", "FIXED"))
plot_index(list(bridging_model_2a, atf_base), model_names = c("CEATTLE", "FIXED"))

