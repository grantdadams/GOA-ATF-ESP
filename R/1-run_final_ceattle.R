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

# * Estimate M ----
ceattle_ss_M <- Rceattle::fit_mod(data_list = mydata_atf,
                                  inits = NULL, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  verbose = 1,
                                  phase = "default",
                                  initMode = 1,
                                  M1Fun = build_M1(M1_model = 2) # Estimate M (sex-specific)
)

# * Rec as random effects
ceattle_ss_M_RE <- Rceattle::fit_mod(data_list = mydata_atf,
                                     inits = ceattle_ss_M$estimated_params, # Initial parameters = 0
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = TRUE, # Random recruitment
                                     msmMode = 0, # Single species mode
                                     verbose = 1,
                                     phase = NULL,
                                     initMode = 1,
                                     M1Fun = build_M1(M1_model = 2) # Estimate M (sex-specific)
)



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

ceattle_ms <- Rceattle::fit_mod(data_list = mydata_atf,
                                inits = ceattle_ms$estimated_params, # Initial parameters = 0
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
                                   loopnum = 3,
                                   msmMode = 1, # Multi-species model
                                   M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)

ceattle_ms_RE <- Rceattle::fit_mod(data_list = mydata_atf,
                                   inits = ceattle_ms_RE$estimated_params, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate
                                   random_rec = TRUE, # Random recruitment
                                   verbose = 1,
                                   phase = NULL,
                                   suit_meanyr = 2015,
                                   initMode = 1,
                                   loopnum = 3,
                                   msmMode = 1, # Multi-species model
                                   M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)

# Save outputs for ESP ----
# Biomass eaten due to cannibalism, annual ration across all ages, and M at age 1-5
years <- ceattle_ms_RE$data_list$styr:ceattle_ms_RE$data_list$endyr
nyrs <- length(years)


weighted_ration <- function(Rceattle, spp = 1, minage = 1, maxage = max(Rceattle$data_list$nages)){
  yrs <- Rceattle$data_list$styr:Rceattle$data_list$endyr
  spp_wt_index <- Rceattle$data_list$pop_wt_index

  wt_f <- Rceattle$data_list$wt %>%
    dplyr::filter(Wt_index == spp_wt_index &
                    Year %in% yrs &
                    Sex == 1) %>%
    dplyr::select(paste0("Age", minage:maxage))

  wt_m <- Rceattle$data_list$wt %>%
    dplyr::filter(Wt_index == spp_wt_index &
                    Year %in% yrs &
                    Sex == 2) %>%
    dplyr::select(paste0("Age", minage:maxage))

  females <- apply(Rceattle$quantities$ration[spp,1,minage:maxage,1:length(yrs)] * Rceattle$quantities$NByage[spp,1,minage:maxage,1:length(yrs)] * wt_f, 1, sum)

  males <- apply(Rceattle$quantities$ration[spp,2,minage:maxage,1:length(yrs)] * Rceattle$quantities$NByage[spp,2,minage:maxage,1:length(yrs)] * wt_m, 1, sum)

  return(
    males + females
    )
}

ESP_data <- data.frame(
  Year = years,
  Beaten_tonnes = apply(ceattle_ms_RE$quantities$B_eaten_as_prey[,,,1:nyrs], 3, sum),
  Ration_tonnes = weighted_ration(ceattle_ms_RE),

  # - M-at-age/sex
  M_f_age1 = ceattle_ms_RE$quantities$M[1,1,1,1:nyrs],
  M_f_age2 = ceattle_ms_RE$quantities$M[1,1,2,1:nyrs],
  M_f_age3 = ceattle_ms_RE$quantities$M[1,1,3,1:nyrs],
  M_f_age4 = ceattle_ms_RE$quantities$M[1,1,4,1:nyrs],
  M_f_age5 = ceattle_ms_RE$quantities$M[1,1,5,1:nyrs],

  M_m_age1 = ceattle_ms_RE$quantities$M[1,2,1,1:nyrs],
  M_m_age2 = ceattle_ms_RE$quantities$M[1,2,2,1:nyrs],
  M_m_age3 = ceattle_ms_RE$quantities$M[1,2,3,1:nyrs],
  M_m_age4 = ceattle_ms_RE$quantities$M[1,2,4,1:nyrs],
  M_m_age5 = ceattle_ms_RE$quantities$M[1,2,5,1:nyrs]
)
write.csv(ESP_data, file = "Results/MS CEATTLE ESP indicators.csv")


# Compare likelihoods ----
source("R/Functions/likelihood comparisons.R", echo=TRUE)

ss_ll <- get_atf_ll(ceattle_ss_RE) %>%
  rename(SS = Value)
ss_M_ll <- get_atf_ll(ceattle_ss_M_RE) %>%
  rename(MS = Value)
ms_ll <- get_atf_ll(ceattle_ms_RE) %>%
  rename(MS = Value)

write.csv(cbind(ss_ll, ss_M_ll[,2], ms_ll[,2]), file = "Results/Final_jnll.csv")



## ABC ----
# abc_list <- lapply(quantities_list, function(x) abc_calc_tmb(x, datlist = dat))


# Compare models ----
# - SAFE model
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)
SAFE2023_comp <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 4)
SAFE2023_index <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 5)
SAFE2023_sel <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 6)

SAFE2023_mod <- ceattle_ss
SAFE2023_mod$quantities$biomass[1,1:length(1977:2023)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$biomassSSB[1,1:length(1977:2023)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:length(1977:2023)] <- SAFE2023$Recruitment/1000
SAFE2023_mod$quantities$sel[1,1,1:21,] <- SAFE2023_sel$Survsel_fem
SAFE2023_mod$quantities$sel[1,2,1:21,] <- SAFE2023_sel$Survsel_mal
SAFE2023_mod$quantities$sel[2,1,1:21,] <- SAFE2023_sel$Survsel_fem
SAFE2023_mod$quantities$sel[2,2,1:21,] <- SAFE2023_sel$Survsel_mal
SAFE2023_mod$quantities$sel[3,1,1:21,] <- SAFE2023_sel$Fishsel_fem
SAFE2023_mod$quantities$sel[3,2,1:21,] <- SAFE2023_sel$Fishsel_mal


SAFE2023_mod$quantities$fsh_bio_hat[1:length(1977:2023)] <- SAFE2023_index$`Pred catch`

# -- Comp data
comp_temp <- SAFE2023_mod$data_list$comp_data %>%
  dplyr::select(-paste0("Comp_", 1:117))

comp_temp <- comp_temp %>%
  dplyr::left_join(SAFE2023_comp %>%
                     dplyr::select(-Sample_size),
                   by = c("Fleet_name", "Fleet_code", "Species", "Sex", "Age0_Length1", "Year", "Month")) # Exclude years not in obs

SAFE2023_mod$quantities$comp_hat <- comp_temp %>%
  dplyr::select(paste0("Comp_", 1:117)) %>%
  as.matrix()

# -- Index
SAFE2023_mod$quantities$srv_bio_hat <- SAFE2023_index %>%
  dplyr::filter(Year %in% SAFE2023_mod$data_list$srv_biom$Year) %>%
  dplyr::pull(Pred)


# - Index model as biomass
index_mod <- ceattle_ss
yrs_srv <- ceattle_ss$data_list$srv_biom$Year - 1977 + 1
index_mod$quantities$biomass[1,] <- 1
index_mod$quantities$biomass[1,yrs_srv] <- ceattle_ss$data_list$srv_biom$Observation


# Plot final ----
MPcols <- rev(oce::oce.colorsViridis(4))
line_col <- c("grey60", MPcols[2:4])
model_list <- list(SAFE2023_mod, ceattle_ss_RE, ceattle_ss_M_RE, ceattle_ms_RE)
model_names = c("ADMB", "TMB single-spp (fix M)", "TMB single-spp (est M)", "TMB multi-spp")

# * Time-series ----
plot_biomass(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_ssb(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_recruitment(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

plot_catch(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

source("R/Functions/Plot_b_eatent_1spp function.R", echo=TRUE)
plot_b_eaten(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

# * Mortality ----
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 2, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 3, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 4, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 5, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 6, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 7, model_names = model_names[2:4], line_col = line_col[2:4])
plot_m_at_age(model_list[2:4], file = "Results/Figures/Final_", width = 5, height = 5, age = 8, model_names = model_names[2:4], line_col = line_col[2:4])

# - M1
round(ceattle_ms_RE$quantities$M1[1,1:2,1], 3)

# - Females average M
round(mean(ceattle_ms_RE$quantities$M[1,1,1,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,1,2,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,1,3,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,1,4,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,1,5,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,1,6,1:length(1977:2023)]), 3)

# - Males average M
round(mean(ceattle_ms_RE$quantities$M[1,2,1,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,2,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,3,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,4,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,5,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,6,1:length(1977:2023)]), 3)
round(mean(ceattle_ms_RE$quantities$M[1,2,7,1:length(1977:2023)]), 3)


# Plot diagnostics ----
# * Comp data ----
dev.off()
plot_comp(SAFE2023_mod, file = "Results/Figures/Diagnostics/Comp/Final_ADMB_"); legend("topright", legend = substitute(paste(bold('ADMB'))), bty = "n")
plot_comp(ceattle_ss_RE, file = "Results/Figures/Diagnostics/Comp/Final_ss_"); legend("topright", legend = substitute(paste(bold('CEATTLE single-spp (fix M)'))), bty = "n")
plot_comp(ceattle_ss_M_RE, file = "Results/Figures/Diagnostics/Comp/Final_ss_M_"); legend("topright", legend = substitute(paste(bold('CEATTLE single-spp (est M)'))), bty = "n")
plot_comp(ceattle_ms_RE, file = "Results/Figures/Diagnostics/Comp/Final_ms_"); legend("topright", legend = substitute(paste(bold('CEATTLE multi-spp'))), bty = "n")

source("R/Functions/Pearson plot function.R", echo=TRUE)
plot_pearson(SAFE2023_mod, file = "Results/Figures/Diagnostics/Comp/Final_ADMB_")
plot_pearson(ceattle_ss_RE, file = "Results/Figures/Diagnostics/Comp/Final_ss_")
plot_pearson(ceattle_ss_M_RE, file = "Results/Figures/Diagnostics/Comp/Final_ss_M_")
plot_pearson(ceattle_ms_RE, file = "Results/Figures/Diagnostics/Comp/Final_ms_")




# * OSA plots ---
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
## devtools::install_github("fishfollower/compResidual/compResidual")
library(ggplot2)
source("R/Functions/Plot osa function.R", echo=TRUE)
plot_osa_comps(SAFE2023_mod, file = "Results/Figures/Diagnostics/OSA/Final_ADMB_", model_name = "ADMB")
plot_osa_comps(ceattle_ss_RE, file = "Results/Figures/Diagnostics/OSA/Final_ss_", model_name = "Single-spp fix M")
plot_osa_comps(ceattle_ss_M_RE, file = "Results/Figures/Diagnostics/OSA/Final_ss_M_", model_name = "Single-spp est M")
plot_osa_comps(ceattle_ms_RE, file = "Results/Figures/Diagnostics/OSA/Final_ms_", model_name = "Multi-spp")

# o <- round(Neff*obs/rowSums(obs),0); p=exp/rowSums(exp)
# ## default output
# res <-  compResidual::resMulti(t(o), t(p))

# * Index fits ----
plot_index(model_list, model_names = model_names, file = "Results/Figures/Diagnostics/Final_", width = 6, height = 3, line_col = line_col)
plot_logindex(model_list, model_names = model_names, file = "Results/Figures/Diagnostics/Final_", width = 6, height = 3, line_col = line_col)


# * Selectivity ----
plot_selectivity(SAFE2023_mod); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('ADMB'))), bty = "n")
plot_selectivity(ceattle_ss_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE single-spp (fix M)'))), bty = "n")
plot_selectivity(ceattle_ss_M_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE single-spp (est M)'))), bty = "n")
plot_selectivity(ceattle_ms_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE multi-spp'))), bty = "n")


# * Retrospectives ----

# ** Single-spp ----
ss_retro <- retrospective(ceattle_ss, peels = 10)
ss_RE_retro <- retrospective(ceattle_ss_RE, peels = 10)

plot_biomass(ss_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_", width = 6, height = 3)
plot_biomass(ss_RE_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_RE_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_RE_", width = 6, height = 3)

# ** Single-spp est M ----
ss_M_retro <- retrospective(ceattle_ss_M, peels = 10)
ss_M_RE_retro <- retrospective(ceattle_ss_M_RE, peels = 10)

plot_biomass(ss_M_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_M_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_M_", width = 6, height = 3)
plot_biomass(ss_M_RE_retro$Rceattle_list, model_names = paste("Mohns =", round(ss_M_RE_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/SS_M_RE_", width = 6, height = 3)

# - M retro plot
endyr <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$data_list$endyr)
M_fem <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$quantities$M[1,1,1,1])
M_males <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$quantities$M[1,2,1,1])

plot(x = endyr, y = M_males, ylab = "M", xlab = "Terminal year", ylim = range(c(M_fem, M_males, 0, 0.4)), type = "l", lwd = 2)
lines(x = endyr, y = M_fem, lwd = 2, col = "blue")
abline(h = 0.2, lty = 2, col = "blue", lwd = 2)
abline(h = 0.35, lty = 2, col = 1, lwd = 2)
legend("bottomleft", legend = c("Females", "Males", "Est", "Fix"), col = c("blue", 1,1,1), lty = c(1,1,1,2), lwd = 2, bty = "n")

# ** Multi-spp ----
ms_retro <- retrospective(ceattle_ms, peels = 10)
ms_RE_retro <- retrospective(ceattle_ms_RE, peels = 10)

plot_biomass(ms_retro$Rceattle_list, model_names = paste("Mohns =", round(ms_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/MS_", width = 6, height = 3)
plot_biomass(ms_RE_retro$Rceattle_list, model_names = paste("Mohns =", round(ms_RE_retro$mohns[1,2], 2)), file = "Results/Figures/Diagnostics/MS_RE_", width = 6, height = 3)

# - M retro plot
endyr <- sapply(ms_RE_retro$Rceattle_list, function(x) x$data_list$endyr)
M_fem <- sapply(ms_RE_retro$Rceattle_list, function(x) x$quantities$M1[1,1,1])
M_males <- sapply(ms_RE_retro$Rceattle_list, function(x) x$quantities$M1[1,2,1])

plot(x = endyr, y = M_males, ylab = "M1", xlab = "Terminal year", ylim = range(c(M_fem, M_males, 0, 0.4)), type = "l", lwd = 2)
lines(x = endyr, y = M_fem, lwd = 2, col = "blue")
abline(h = 0.2, lty = 2, col = "blue", lwd = 2)
abline(h = 0.35, lty = 2, col = 1, lwd = 2)
legend("bottomleft", legend = c("Females", "Males", "Est", "Fix"), col = c("blue", 1,1,1), lty = c(1,1,1,2), lwd = 2, bty = "n")


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
    mod_prof <- tryCatch(
      fit_mod(
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
        verbose = 0),
      error = function(ex) {
        return(NULL) }
    )

    mod_prof
  }

  closeAllConnections()
  gc()

  return(profile_list)
}


# -- Run profile
profile1 <- profile_rsigma(model = ceattle_ss, rsigma_vec, species = 1)
profile2 <- profile_rsigma(model = ceattle_ss_RE, rsigma_vec, species = 1)

profile3 <- profile_rsigma(model = ceattle_ss_M, rsigma_vec, species = 1)
profile4 <- profile_rsigma(model = ceattle_ss_M_RE, rsigma_vec, species = 1)

profile5 <- profile_rsigma(model = ceattle_ms, rsigma_vec, species = 1)
profile6 <- profile_rsigma(model = ceattle_ms_RE, rsigma_vec, species = 1)

# -- Combine
pml_models <- list(profile1, profile3, profile5)
mml_models <- list(profile2, profile4, profile6)


# -- Plot
par(mfrow = c(1,3))

for(i in 1:3){

  # -- Plot penalized likelihood profile
  y = sapply(pml_models[[i]], function(x) x$opt$objective)
  y = y-min(y)

  plot(y = y, x = rsigma_vec, ylab = "dNLL", xlab = "sigmaR", type = "l", main = c("Single-spp (fix M)", "Single-spp (est M)", "Multi-spp")[i], col = "red", ylim = c(0,10))

  # -- Plot marginalized maximum likelihood profile
  null_opt <- sapply(mml_models[[i]], function(x) !is.null(x))
  y = sapply(mml_models[[i]][null_opt], function(x) x$opt$objective)
  y = y-min(y)
  lines(y = y, x = rsigma_vec[null_opt], col = 1)

  # -- Plot MLE
  abline(v = exp(list(ceattle_ss_RE, ceattle_ss_M_RE, ceattle_ms_RE)[[i]]$estimated_params$ln_rec_sigma[1]), lty = 2)
}

legend("bottomright", c("Penalized likelihood", "Random effects", "Minima"), col = c(2,1,1), lty = c(1,1,2), bty = "n")


