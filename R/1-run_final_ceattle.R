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


# Single-species models ----
# - Fit single-species models and no fishing
# * Fix M ----
ceattle_ss <- Rceattle::fit_mod(data_list = mydata_atf,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single species mode
                                verbose = 1,
                                phase = "default",
                                initMode = 1)

# * Rec as random effects
ceattle_ss_RE <- Rceattle::fit_mod(data_list = mydata_atf,
                                   inits = ceattle_ss$estimated_params, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate
                                   random_rec = TRUE, # Random recruitment
                                   msmMode = 0, # Single species mode
                                   verbose = 1,
                                   phase = NULL,
                                   initMode = 1)



# Multi-species model ----
# (cannibalism)
ceattle_ms <- Rceattle::fit_mod(data_list = mydata_atf,
                                inits = ceattle_ss$estimated_params, # Initial parameters = 0
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

# * Rec as random effects
ceattle_ms_RE <- Rceattle::fit_mod(data_list = mydata_atf,
                                   inits = ceattle_ss_RE$estimated_params, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate
                                   random_rec = TRUE, # Random recruitment
                                   verbose = 1,
                                   phase = NULL,
                                   suit_meanyr = 2015,
                                   initMode = 1,
                                   msmMode = 1, # Multi-species model
                                   M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)


# Compare likelihoods ----
source("R/Functions/likelihood comparisons.R", echo=TRUE)

ss_ll <- get_atf_ll(ceattle_ss_RE) %>%
  rename(SS = Value)
ms_ll <- get_atf_ll(ceattle_ms_RE) %>%
  rename(MS = Value)

write.csv(cbind(ss_ll, ms_ll[,2]), file = "Results/Final_jnll.csv")



## ABC ----
abc_list <- lapply(quantities_list, function(x) abc_calc_tmb(x, datlist = dat))


# Compare models ----
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)
SAFE2023_comp <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 4)
SAFE2023_index <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 5)

SAFE2023_mod <- ceattle_ss
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000

# -- Comp data
comp_temp <- SAFE2023_mod$data_list$comp_data %>%
  dplyr::select(-paste0("Comp_", 1:117))

comp_temp <- comp_temp %>%
  left_join(SAFE2023_comp) # Exclude years not in obs

SAFE2023_mod$quantities$comp_hat <- comp_temp %>%
  dplyr::select(paste0("Comp_", 1:117)) %>%
  as.matrix()

# -- Index
SAFE2023_mod$quantities$srv_bio_hat <- SAFE2023_index %>%
  dplyr::filter(Year %in% SAFE2023_mod$data_list$srv_biom$Year) %>%
  dplyr::pull(Pred)


# Plot final ----
MPcols <- rev(oce::oce.colorsViridis(3))
line_col <- c("grey60", MPcols[2:3])
model_list <- list(SAFE2023_mod, ceattle_ss_RE, ceattle_ms_RE)
model_names = c("ADMB", "TMB single-spp", "TMB multi-spp")

plot_biomass(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_ssb(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_recruitment(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

line_col <- MPcols[2:3]
model_list <- list(ceattle_ss_RE, ceattle_ms_RE)
model_names = c("TMB single-spp", "TMB multi-spp")
plot_catch(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

source("R/Functions/Plot_b_eatent_1spp function.R", echo=TRUE)
plot_b_eaten(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)


# Plot diagnostics ----
# * Comp data ----
dev.off()
source("R/Functions/Pearson plot function.R", echo=TRUE)
plot_pearson(SAFE2023_mod, file = "Results/Figures/Diagnostics/Final_ADMB_")
plot_pearson(ceattle_ss_RE, file = "Results/Figures/Diagnostics/Final_ss_")
plot_pearson(ceattle_ms_RE, file = "Results/Figures/Diagnostics/Final_ms_")

# * Index fits ----
line_col <- c("grey60", MPcols[2:3])
model_list <- list(SAFE2023_mod, ceattle_ss_RE, ceattle_ms_RE)
model_names = c("ADMB", "TMB single-spp", "TMB multi-spp")
plot_index(model_list, model_names = model_names, file = "Results/Figures/Diagnostics/Final_", width = 6, height = 3, line_col = line_col)
plot_logindex(model_list, model_names = model_names, file = "Results/Figures/Diagnostics/Final_", width = 6, height = 3, line_col = line_col)


# * Run retrospectives ----
ss_retro <- retrospective(ceattle_ss, peels = 10)
ss_RE_retro <- retrospective(ceattle_ss_RE, peels = 10)

ms_retro <- retrospective(ceattle_ms, peels = 10)
ms_RE_retro <- retrospective(ceattle_ms_RE, peels = 10)

plot_biomass(ss_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_", width = 6, height = 3)
plot_biomass(ss_RE_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_RE_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_RE_", width = 6, height = 3)

plot_biomass(ms_retro$Rceattle_list, model_names = paste("Mohns =", round(ms_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/MS_", width = 6, height = 3)
plot_biomass(ms_RE_retro$Rceattle_list, model_names = paste("Mohns =", round(ms_RE_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/MS_RE_", width = 6, height = 3)


# * Profile sigmaR ----
rsigma_vec <- seq(from = 0.05, to = 2, by = 0.05)

profile_rsigma <- function(model = NULL, rsigma_vec = NULL, species = NULL){
  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = detectCores() - 6
  registerDoParallel(cores)

  # Loop through Rsigma
  profile_list <- foreach(i = 1:length(rsigma_vec)) %dopar% {
    library(Rceattle)
    library(dplyr)

    # Update sigmaR
    inits <- model$estimated_params
    inits$ln_rec_sigma[species] <- log(rsigma_vec[i])

    # Build map
    data_list <- model$data_list
    # data_list$estDynamics <- rep(1, data_list$nspp)
    # data_list$estDynamics[species] <- 0
    map <- Rceattle::build_map(data_list, params = inits, debug = FALSE, random_rec = FALSE)

    # Estimate
    mod_prof <- fit_mod(
      data_list = data_list,
      inits = inits,
      map =  map,
      bounds = NULL,
      file = NULL,
      estimateMode = 1,
      HCR = build_hcr(HCR = model$data_list$HCR, # Tier3 HCR
                      DynamicHCR = model$data_list$DynamicHCR,
                      FsprTarget = model$data_list$FsprTarget,
                      FsprLimit = model$data_list$FsprLimit,
                      Ptarget = model$data_list$Ptarget,
                      Plimit = model$data_list$Plimit,
                      Alpha = model$data_list$Alpha,
                      Pstar = model$data_list$Pstar,
                      Sigma = model$data_list$Sigma,
                      Fmult = model$data_list$Fmult,
                      HCRorder = model$data_list$HCRorder
      ),
      recFun = build_srr(srr_fun = model$data_list$srr_fun,
                         srr_pred_fun  = model$data_list$srr_pred_fun ,
                         proj_mean_rec  = model$data_list$proj_mean_rec ,
                         srr_meanyr = model$data_list$srr_meanyr,
                         R_hat_yr = model$data_list$R_hat_yr,
                         srr_est_mode  = model$data_list$srr_est_mode ,
                         srr_prior_mean  = model$data_list$srr_prior_mean,
                         srr_prior_sd   = model$data_list$srr_prior_sd,
                         Bmsy_lim = model$data_list$Bmsy_lim,
                         srr_env_indices = model$data_list$srr_env_indices),
      M1Fun = build_M1(M1_model= model$data_list$M1_model,
                       updateM1 = FALSE,
                       M1_use_prior = model$data_list$M1_use_prior,
                       M2_use_prior = model$data_list$M2_use_prior,
                       M1_prior_mean = model$data_list$M1_prior_mean,
                       M1_prior_sd = model$data_list$M1_prior_sd),
      random_rec = model$data_list$random_rec,
      niter = model$data_list$niter,
      msmMode = model$data_list$msmMode,
      avgnMode = model$data_list$avgnMode,
      suitMode = model$data_list$suitMode,
      suit_meanyr = model$data_list$suit_meanyr,
      initMode = model$data_list$initMode,
      phase = NULL,
      loopnum = 1,
      getsd = FALSE,
      verbose = 0)

    mod_prof
  }

  closeAllConnections()
  gc()

  return(profile_list)
}


# - Run profile
profile1 <- profile_rsigma(model = ceattle_ss, rsigma_vec, species = 1)
profile2 <- profile_rsigma(model = ceattle_ss_RE, rsigma_vec, species = 1)

profile3 <- profile_rsigma(model = ceattle_ms, rsigma_vec, species = 1)
profile4 <- profile_rsigma(model = ceattle_ms_RE, rsigma_vec, species = 1)

# - Combine
pml_models <- list(profile1, profile3)
mml_models <- list(profile2, profile4)


# - Plot
par(mfrow = c(1,2))

for(i in 1:2){

  # -- Plot penalized likelihood profile
  y = sapply(pml_models[[i]], function(x) x$opt$objective)
  y = y-min(y)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = c("Single-spp", "Multi-spp")[i], col = "red", ylim = c(0,10))

  # -- Plot marginalized maximum likelihood profile
  y = sapply(mml_models[[i]], function(x) x$opt$objective)
  y = y-min(y)
  lines(y = y, x = rsigma_vec, col = 1)

  # -- Plot MLE
  abline(v = exp(list(ceattle_ss_RE, ceattle_ms_RE)[[i]]$estimated_params$ln_rec_sigma[1]), lty = 2)
}

legend("topright", c("Penalized likelihood", "Random effects", "Minima"), col = c(2,1,1), lty = c(1,1,2), bty = "n")


