#' Plot biomass eaten
#'
#' @description Function the plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include time series mean as horizontal line
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_b_eaten <-  function(Rceattle,
                          file = NULL,
                          model_names = NULL,
                          line_col = NULL,
                          species = NULL,
                          spnames = NULL,
                          add_ci = FALSE,
                          lwd = 3,
                          save = FALSE,
                          right_adj = 0,
                          width = 7,
                          height = 6.5,
                          minyr = NULL,
                          incl_proj = FALSE,
                          mod_cex = 1,
                          alpha = 0.4,
                          mod_avg = rep(FALSE, length(Rceattle)),
                          mse = FALSE,
                          OM = TRUE) {

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(years, max)))
  if(is.null(minyr)){minyr <- min((sapply(years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  maxage <- max(Rceattle[[1]]$data_list$nages)
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(0, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))

  for (i in 1:length(Rceattle)) {

    # - Get quantities
    quantity[,1:nyrs_vec[i] , i] <- apply(Rceattle[[i]]$quantities$B_eaten_as_prey[,,,1:nyrs_vec[i]], 3, sum)

    # # Get SD of quantity
    # # NOTE: No uncertainty estimates currently
    # if(add_ci & !mse){
    #   sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "B_eaten_as_prey")
    #   sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
    #   quantity_sd[,  1:nyrs_vec[i], i] <-
    #     replace(quantity_sd[,,,, i], values = sd_temp)
    # }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,1:nyrs_vec[i], i] <- apply(
        apply(Rceattle[[i]]$asymptotic_samples$B_eaten_as_prey[,,,1:nyrs_vec[i],], c(3), function(x) sum), # Sum across age-sex
        c(1,2), sd(as.vector(log(x)))) # SD across samples
      log_quantity_mu[,1:nyrs_vec[i], i] <- apply(
        apply(Rceattle[[i]]$asymptotic_samples$B_eaten_as_prey[,,,1:nyrs_vec[i],], c(3), function(x) sum), # Sum across age-sex
        c(1,2), mean(as.vector(log(x)))) # Mean across samples
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.25) )
    quantity <- apply( quantity, c(1,2), mean ) # Get mean quantity

    # -- Put back in array for indexing below
    quantity <- array(quantity, dim = c(nspp, nyrs,  1))
    quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
    quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
    quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
    quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
  }

  # - Model Average
  for (i in 1:length(Rceattle)) {
    if(mod_avg[i]){
      quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
    }
  }

  # - Rescale
  quantity <- quantity
  quantity_upper95 <- quantity_upper95
  quantity_lower95 <- quantity_lower95
  quantity_upper50 <- quantity_upper50
  quantity_lower50 <- quantity_lower50


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_b_eaten_as_prey_trajectory", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_b_eaten_as_prey_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = "Biomass consumed (mt)",
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$endyr, lwd  = lwd, col = "grey", lty = 2)
      }

      # Legends
      legend("topleft",
             legend = spnames[spp[j]],
             bty = "n",
             cex = 1)

      if (spp[j] == 1) {
        if(!is.null(model_names)){
          legend(
            "topright",
            legend = model_names,
            lty = rep(1, length(line_col)),
            lwd = lwd,
            col = line_col,
            bty = "n",
            cex = mod_cex
          )
        }
      }


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(years[[k]], rev(years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(years[[k]], rev(years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = years[[k]],
          y = quantity[spp[j], 1:length(years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}
