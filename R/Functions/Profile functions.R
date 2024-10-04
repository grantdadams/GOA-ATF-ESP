# Functions to profile Rceattle parameters for a single-species


profile_rsigma <- function(model = NULL, rsigma_vec = NULL, species = NULL, filename = NULL){
  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = detectCores() - 6
  registerDoParallel(cores)

  # Loop through Rsigma
  profile_list <- foreach(i = 1:length(rsigma_vec),
                          .packages = c("Rceattle", "dplyr")) %dopar% {

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
                                file = ifelse(is.null(filename), NULL, paste0(filename, rsigma_vec[i])),
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



profile_m <- function(model = NULL, m_vec = NULL, species = NULL, filename = NULL){
  ### Set up parallel processing
  library(foreach)
  library(doParallel)

  cores = detectCores() - 6
  registerDoParallel(cores)
  nsex <- model$data_list$nsex[species]

  # Pointers for M1 in two-sex models used below
  m1_pointer <- 1:(length(m_vec))
  if(nsex == 2){
    m1_pointer <- matrix(1:(length(m_vec) ^ nsex), length(m_vec), length(m_vec))
  }

  # Loop through Rsigma
  profile_list <- foreach(i = 1:(length(m_vec) ^ nsex),
                          .packages = c("Rceattle", "dplyr"),
                          .export = c("m1_pointer", "nsex")) %dopar% {

                            # Update M1
                            inits <- model$estimated_params
                            inits$ln_M1[species, 1,] <- log(m_vec[i])

                            if(nsex == 2){
                              inits$ln_M1[species, 1,] <- log(m_vec[which(m1_pointer == i, arr.ind = TRUE)[1]])
                              inits$ln_M1[species, 2,] <- log(m_vec[which(m1_pointer == i, arr.ind = TRUE)[2]])
                            }

                            # Map out M!
                            model$map$mapList$ln_M1[] <- NA
                            model$map$mapFactor$ln_M1 <- factor(model$map$mapList$ln_M1)

                            # Estimate
                            mod_prof <- tryCatch(
                              fit_mod(
                                data_list = model$data_list,
                                inits = inits,
                                map =  model$map,
                                bounds = NULL,
                                # file = ifelse(is.null(filename), NULL, paste0(filename, rsigma_vec[i])),
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
                                loopnum = 3,
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
