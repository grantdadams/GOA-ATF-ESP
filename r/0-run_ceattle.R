library(Rceattle)
library(readxl)
library(dplyr)


################################################
# Single-species mode
################################################
mydata_atf <- Rceattle::read_data( file = "Data/GOA_23.1.1_arrowtooth_single_species_1977-2023.xlsx")
mydata_atf$estDynamics = 0
mydata_atf$srv_biom$Log_sd <- mydata_atf$srv_biom$Log_sd/mydata_atf$srv_biom$Observation
mydata_atf$fleet_control$proj_F_prop <- c(1,1,1)

# - Fit single-species models with tier-3 HCR
atf_base <- Rceattle::fit_mod(data_list = mydata_atf,
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



#########################################################
# Model 2 - Fix n-at-age and parameters to 2023 SAFE-multinomial selectivity values
#########################################################
mydata_atf_fixed <- mydata_atf
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
inits$init_dev[1,1:20] <- rec_devs[1:20]

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
atf_fixed_params <- Rceattle::fit_mod(data_list = mydata_atf_fixed,
                                       inits = inits, # Initial parameters = 0
                                       file = NULL, # Don't save
                                       estimateMode = 4, # Estimate
                                       random_rec = FALSE, # No random recruitment
                                       msmMode = 0, # Single species mode
                                       verbose = 1,
                                       phase = "default",
                                       initMode = 1)

t(atf_fixed_params$quantities$NByage[1,1,,1:10])/1000
head(mydata_atf_fixed$NByageFixed)


################################################
# Multi-species mode (cannibalism)
################################################


################################################
# Compare models
################################################
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)

SAFE2023_mod <- atf_base
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000

# - SAFE model with fixed multinomial
SAFE2023multi <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 2)

SAFE2023multi_mod <- atf_base
SAFE2023multi_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023multi$Biomass
SAFE2023multi_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023multi$SSB
SAFE2023multi_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023multi$Recruitment/1000

# - SAFE model with fixed multinomial and selectivity penalties set to the same for sexes (lines 470-472 of dat file)
SAFE2023multisel <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 3)

SAFE2023multisel_mod <- atf_base
SAFE2023multisel_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023multisel$Biomass
SAFE2023multisel_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023multisel$SSB
SAFE2023multisel_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023multisel$Recruitment/1000

# CEATTLE with fixed parameters
atf_fixed_params$quantities$biomass <- atf_fixed_params$quantities$biomass/1000
atf_fixed_params$quantities$biomassSSB <- atf_fixed_params$quantities$biomassSSB/1000
atf_fixed_params$quantities$R <- atf_fixed_params$quantities$R/1000
atf_fixed_params$quantities$fsh_bio_hat <- atf_fixed_params$quantities$fsh_bio_hat/1000
atf_fixed_params$quantities$srv_bio_hat <- atf_fixed_params$quantities$srv_bio_hat/1000

#
model_list <- list(atf_base, SAFE2023_mod, SAFE2023multi_mod, SAFE2023multisel_mod, atf_fixed_params)
model_names = c("CEATTLE", "SAFE", "SAFE-multi", "SAFE-multi-sel", "Fixed")

plot_biomass(model_list, model_names = model_names)
plot_ssb(model_list, model_names = model_names)
plot_recruitment(model_list, model_names = model_names)
plot_catch(model_list, model_names = model_names)
plot_index(model_list, model_names = model_names)

