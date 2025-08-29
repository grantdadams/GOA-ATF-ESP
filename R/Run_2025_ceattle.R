# Code used to fit the 2025 arrowtooth flounder assessment in CEATTLE
# uses the "dev-name-change" version of Rceattle

library(Rceattle)
library(readxl)
library(dplyr)
library(TMB)

# Load data ----
# * 2025 data ----
data_2025 <- Rceattle::read_data( file = "Data/GOA_25.1.1_arrowtooth_single_species_1977-2024.xlsx")
data_2025$estDynamics = 0
data_2025$index_data$Log_sd <- data_2025$index_data$Log_sd/data_2025$index_data$Observation
data_2025$fleet_control$proj_F_prop <- c(1,1,1)

# * 2023 data ----
data_2023 <- Rceattle::read_data( file = "Data/GOA_23.1.1_arrowtooth_single_species_1977-2023.xlsx")
data_2023$estDynamics = 0
data_2023$index_data$Log_sd <- data_2023$index_data$Log_sd/data_2023$index_data$Observation
data_2023$fleet_control$proj_F_prop <- c(1,1,1)


# 1) 2023 single-species models ----
# - Fit single-species models and no fishing and fix M
ceattle_ss_2023 <- Rceattle::fit_mod(data_list = data_2023,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = TRUE, # No random recruitment
                                msmMode = 0, # Single species mode
                                verbose = 1,
                                phase = FALSE,
                                initMode = 1)
plot_biomass(list(ceattle_ss_2023),
             model_names = c("2023 Rceattle"),
             file = "Results/Model comparison")

# 2) 2025 new data/old weighting ----
# * Penalized likelihood ----
ceattle_ss_2025 <- Rceattle::fit_mod(data_list = data_2025,
                                     inits = NULL, # Initial parameters = 0
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = FALSE, # No random recruitment
                                     msmMode = 0, # Single species mode
                                     verbose = 1,
                                     phase = FALSE,
                                     initMode = 1)

# * Random recruitment ----
ceattle_ss_2025_RE <- Rceattle::fit_mod(data_list = data_2025,
                                   inits = ceattle_ss_2025$estimated_params, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   estimateMode = 0, # Estimate
                                   random_rec = TRUE, # Random recruitment
                                   msmMode = 0, # Single species mode
                                   verbose = 1,
                                   phase = FALSE,
                                   initMode = 1)

plot_biomass(list(ceattle_ss_2023, ceattle_ss_2025_RE),
             model_names = c("2023 Rceattle","2025 Rceattle"),
             file = "Results/Model comparison")


# 3) 2025 Dirichlet multinomial ----
data_2025$fleet_control$Comp_weights <- c(1,1,1)
data_2025$fleet_control$Comp_loglike <- c(1,1,1) # DM

# * Penalized likelihood ----
ceattle_ss_2025_DM <- Rceattle::fit_mod(data_list = data_2025,
                                     inits = NULL, # Initial parameters = 0
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = FALSE, # No random recruitment
                                     msmMode = 0, # Single species mode
                                     verbose = 1,
                                     phase = FALSE,
                                     initMode = 1)

# * Random recruitment ----
ceattle_ss_2025_DM_RE <- Rceattle::fit_mod(data_list = data_2025,
                                        inits = ceattle_ss_2025_DM$estimated_params, # Initial parameters = 0
                                        file = NULL, # Don't save
                                        estimateMode = 0, # Estimate
                                        random_rec = TRUE, # Random recruitment
                                        msmMode = 0, # Single species mode
                                        verbose = 1,
                                        phase = FALSE,
                                        initMode = 1)

plot_biomass(list(ceattle_ss_2023, ceattle_ss_2025_RE, ceattle_ss_2025_DM_RE),
             model_names = c("2023 Rceattle","2025 Rceattle", "2025 Rceattle DM"),
             file = "Results/Model comparison")


# 4) 2025 Estimate M ----
# * Penalized likelihood ----
ceattle_ss_2025_M <- Rceattle::fit_mod(data_list = data_2025,
                                  inits = ceattle_ss_2025_DM_RE$estimated_params, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  estimateMode = 0, # Estimate
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 0, # Single species mode
                                  verbose = 1,
                                  phase = FALSE,
                                  initMode = 1,
                                  M1Fun = build_M1(M1_model = 2) # Estimate M (sex-specific)
)

# * Random recruitment ----
ceattle_ss_2025_M_RE <- Rceattle::fit_mod(data_list = data_2025,
                                     inits = ceattle_ss_2025_DM_RE$estimated_params, # Initial parameters = 0
                                     file = NULL, # Don't save
                                     estimateMode = 0, # Estimate
                                     random_rec = TRUE, # Random recruitment
                                     msmMode = 0, # Single species mode
                                     verbose = 1,
                                     phase = FALSE,
                                     initMode = 1,
                                     M1Fun = build_M1(M1_model = 2) # Estimate M (sex-specific)
)

plot_biomass(list(ceattle_ss_2023, ceattle_ss_2025_RE, ceattle_ss_2025_DM_RE, ceattle_ss_2025_M_RE),
             model_names = c("2023 Rceattle","2025 Rceattle", "2025 Rceattle DM", "2025 Rceattle M & DM"),
             file = "Results/Model comparison")



# 5) 2025 multi-species model ----
# (cannibalism)
# * Penalized likelihood ----
ceattle_ms_2025 <- Rceattle::fit_mod(data_list = data_2025,
                                inits = ceattle_ss_2025$estimated_params, # Initial parameters = 0
                                file = NULL, # Don't save
                                estimateMode = 0, # Estimate
                                random_rec = FALSE, # No random recruitment
                                verbose = 1,
                                phase = FALSE,
                                suit_styr = 1990,
                                suit_endyr = 2015,
                                initMode = 1,
                                msmMode = 1, # Multi-species model
                                M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)

# * Random recruitment ----
ceattle_ms_2025_RE <- Rceattle::fit_mod(data_list = data_2025,
                                   inits = ceattle_ms_2025$estimated_params, # Initial parameters = 0
                                   file = NULL,       # Don't save
                                   estimateMode = 0,  # Estimate
                                   random_rec = TRUE, # Random recruitment
                                   verbose = 1,
                                   phase = FALSE,
                                   suit_styr = 1990,
                                   suit_endyr = 2015,
                                   initMode = 1,
                                   loopnum = 3,
                                   msmMode = 1, # Multi-species model
                                   M1Fun = build_M1(M1_model = 2) # Estimate residual M (sex-specific)
)


# Compare model updates ----
years <- ceattle_ss$data_list$styr:ceattle_ss$data_list$endyr
nyrs <- length(years)


# * Get old ADMB output ----
SAFE2023 <- read_excel("Data/2023_SAFE_biomass_estimate.xlsx", sheet = 1)
SAFE2023_mod <- ceattle_ss_2023
SAFE2023_mod$quantities$biomass[1,1:(nyrs-2)] <- SAFE2023$Biomass
SAFE2023_mod$quantities$ssb[1,1:(nyrs-2)] <- SAFE2023$SSB
SAFE2023_mod$quantities$R[1,1:(nyrs-2)] <- SAFE2023$Recruitment/1000


plot_biomass(list(ceattle_ss, ceattle_ss_RE, SAFE2023_mod, ceattle_ss_2023, ceattle_ss_2023_dw, ceattle_ss_2023_dw_data, ceattle_ss_2023_endyr),
             model_names = c("2025 Rceattle", "2025 Rceattle RE", "2023 ADMB","2023 Rceattle", "2023 DW Rceattle", "2023 DW no old data Rceattle", "2023 w/ 2021 age Rceattle"),
             file = "Results/Model comparison")


# Save outputs for ESP ----
# Biomass eaten due to cannibalism, annual ration across all ages, and M at age 1-5


weighted_ration <- function(Rceattle, spp = 1, minage = 1, maxage = max(Rceattle$data_list$nages)){
  yrs <- Rceattle$data_list$styr:Rceattle$data_list$endyr
  spp_wt_index <- Rceattle$data_list$pop_wt_index

  wt_f <- Rceattle$data_list$weight %>%
    dplyr::filter(Wt_index == spp_wt_index &
                    Sex == 1) %>%
    dplyr::select(-Year) %>%
    dplyr::cross_join(data.frame(Year = yrs)) %>%
    dplyr::select(paste0("Age", minage:maxage))

  wt_m <- Rceattle$data_list$weight %>%
    dplyr::filter(Wt_index == spp_wt_index &
                    Sex == 2)  %>%
    dplyr::select(-Year) %>%
    dplyr::cross_join(data.frame(Year = yrs)) %>%
    dplyr::select(paste0("Age", minage:maxage))

  females <- apply(Rceattle$quantities$ration[spp,1,minage:maxage,1:length(yrs)] * Rceattle$quantities$N_at_age[spp,1,minage:maxage,1:length(yrs)] * wt_f, 1, sum)

  males <- apply(Rceattle$quantities$ration[spp,2,minage:maxage,1:length(yrs)] * Rceattle$quantities$N_at_age[spp,2,minage:maxage,1:length(yrs)] * wt_m, 1, sum)

  return(
    males + females
    )
}

ESP_data <- data.frame(
  Year = years,
  Beaten_tonnes = apply(ceattle_ms_RE$quantities$B_eaten_as_prey[,,,1:nyrs], 3, sum),
  Ration_tonnes = weighted_ration(ceattle_ms_RE),

  # - M-at-age/sex
  M_f_age1 = ceattle_ms_RE$quantities$M_at_age[1,1,1,1:nyrs],
  M_f_age2 = ceattle_ms_RE$quantities$M_at_age[1,1,2,1:nyrs],
  M_f_age3 = ceattle_ms_RE$quantities$M_at_age[1,1,3,1:nyrs],
  M_f_age4 = ceattle_ms_RE$quantities$M_at_age[1,1,4,1:nyrs],
  M_f_age5 = ceattle_ms_RE$quantities$M_at_age[1,1,5,1:nyrs],

  M_m_age1 = ceattle_ms_RE$quantities$M_at_age[1,2,1,1:nyrs],
  M_m_age2 = ceattle_ms_RE$quantities$M_at_age[1,2,2,1:nyrs],
  M_m_age3 = ceattle_ms_RE$quantities$M_at_age[1,2,3,1:nyrs],
  M_m_age4 = ceattle_ms_RE$quantities$M_at_age[1,2,4,1:nyrs],
  M_m_age5 = ceattle_ms_RE$quantities$M_at_age[1,2,5,1:nyrs]
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


# Plot final ----
MPcols <- rev(oce::oce.colorsViridis(4))
line_col <- MPcols[2:4]
model_list <- list(ceattle_ss_RE, ceattle_ss_M_RE, ceattle_ms_RE)
model_names = c("TMB single-spp (fix M)", "TMB single-spp (est M)", "TMB multi-spp")

# * Time-series ----
plot_biomass(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_ssb(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)
plot_recruitment(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

plot_catch(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

plot_b_eaten(model_list, model_names = model_names, file = "Results/Figures/Final_", width = 6, height = 3, line_col = line_col)

# * Mortality ----
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 2, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 3, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 4, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 5, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 6, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 7, model_names = model_names, line_col = line_col[2:4])
plot_m_at_age(model_list, file = "Results/Figures/Final_", width = 5, height = 5, age = 8, model_names = model_names, line_col = line_col[2:4])

# - M1
round(ceattle_ms_RE$quantities$M1_at_age[1,1:2,1,1], 3)

# - Females average M
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,1,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,2,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,3,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,4,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,5,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,1,6,1:nyrs]), 3)

# - Males average M
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,1,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,2,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,3,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,4,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,5,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,6,1:nyrs]), 3)
round(mean(ceattle_ms_RE$quantities$M_at_age[1,2,7,1:nyrs]), 3)


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
plot_selectivity(ceattle_ss_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE single-spp (fix M)'))), bty = "n")
plot_selectivity(ceattle_ss_M_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE single-spp (est M)'))), bty = "n")
plot_selectivity(ceattle_ms_RE); legend(y = 0.15, x = 12.5, legend = substitute(paste(bold('CEATTLE multi-spp'))), bty = "n")


# Diagnostics ----
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


# * M retro plot ----
endyr <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$data_list$endyr)
M_fem <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$quantities$M_at_age[1,1,1,1])
M_males <- sapply(ss_M_RE_retro$Rceattle_list, function(x) x$quantities$M_at_age[1,2,1,1])

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
M_fem <- sapply(ms_RE_retro$Rceattle_list, function(x) x$quantities$M1_at_age[1,1,1])
M_males <- sapply(ms_RE_retro$Rceattle_list, function(x) x$quantities$M1_at_age[1,2,1])

plot(x = endyr, y = M_males, ylab = "M1", xlab = "Terminal year", ylim = range(c(M_fem, M_males, 0, 0.4)), type = "l", lwd = 2)
lines(x = endyr, y = M_fem, lwd = 2, col = "blue")
abline(h = 0.2, lty = 2, col = "blue", lwd = 2)
abline(h = 0.35, lty = 2, col = 1, lwd = 2)
legend("bottomleft", legend = c("Females", "Males", "Est", "Fix"), col = c("blue", 1,1,1), lty = c(1,1,1,2), lwd = 2, bty = "n")


# * Profile sigmaR ----
source("R/Functions/Profile functions.R", echo=TRUE)
rsigma_vec <- seq(from = 0.05, to = 2, by = 0.05)

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


# * Profile M ----
library(ggplot2)
m_vec <- seq(from = 0.1, to = 0.8, by = 0.02)
m_profile <- profile_m(model = ceattle_ss_M_RE, m_vec = m_vec, species = 1)

null_vec <- sapply(m_profile, is.null)
ll_vec <- sapply(m_profile, function(x) ifelse(is.null(x), NA, x$opt$objective))
ll_mat <- matrix(ll_vec, length(m_vec), length(m_vec))
contour(x = m_vec, y = m_vec, z = ll_mat, xlab = "Male M" , ylab = "Female M")


ll_df <- data.frame(F_Mort = sapply(m_profile, function(x) ifelse(is.null(x), NA, x$quantities$M_at_age[1,1,1,1])),
                    M_mort = sapply(m_profile, function(x) ifelse(is.null(x), NA, x$quantities$M_at_age[1,2,1,1])),
                    LL = ll_vec)

p1 <- ggplot(ll_df, aes(y = F_Mort, x = M_mort, z = LL)) +
  geom_contour( bins = 15) +
  geom_contour_filled( bins = 15) +
  theme_classic() +
  # xlim(0.1, 0.4) +
  theme(legend.position="none")

p1

ggsave(filename = "Results/Figures/M profile.png", width = 6, height = 6, units = "in")


