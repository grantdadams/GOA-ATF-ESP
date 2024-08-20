#' Plot time series of survey comp data
#'
#' @description Function the plots the survey comp data as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param species Species names for legend
#' @param cex Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#'
#' @return Returns and saves a figure with the catch trajectory.
#' @export
plot_pearson <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           species = NULL,
           cex = 3,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           max_z = NULL,
           incl_proj = FALSE) {

    # Make sure we are using only one model
    if(class(Rceattle) != "Rceattle"){
      stop("Please only use one Rceattle model")
    }

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }


    # Get colors
    colvec=c("red", "blue", "black")

    # Extract data objects
    # Get observed
    comp_data <- Rceattle$data_list$comp_data
    # Normalize
    comp_data[,grep("Comp_", colnames(comp_data))] = comp_data[,grep("Comp_", colnames(comp_data))]/rowSums(comp_data[,grep("Comp_", colnames(comp_data))], na.rm = TRUE)

    fleet_control <- Rceattle$data_list$fleet_control


    # Get estimated
    comp_hat <- Rceattle$data_list$comp_data
    comp_hat[,grep("Comp_", colnames(comp_data))] = Rceattle$quantities$comp_hat

    #############################################
    # Plot comps Type 2 - Pearson residual
    #############################################
    for(comp_type in c(0,1)){ # Age0, Length 1

      srv <- unique(comp_data$Fleet_code[which(comp_data$Age0_Length1 == comp_type)])
      nsrv <- length(srv)


      if(nsrv > 0){
        loops <- ifelse(is.null(file), 1, 2)
        for (i in 1:loops) {

          # Plot configuration
          par(
            mar = c(3, 3 , 1 , 1) ,
            oma = c(0 , 0 , 0 , 0),
            tcl = -0.35,
            mgp = c(1.75, 0.5, 0)
          )

          # Loop around fleets with comp data
          for (j in 1:nsrv) {

            #  Plot name
            if (i == 2) {
              filename <- paste0(file, paste0(c("_comps_pearson_residual_", "_comps_pearson_residual_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
              png(
                file = filename ,
                width = 7,# 169 / 25.4,
                height = 6.5,# 150 / 25.4,

                units = "in",
                res = 300
              )
            }


            # Extract comps
            srv_ind <- which(comp_data$Fleet_code == srv[j] & comp_data$Age0_Length1 == comp_type)
            comp_tmp <- comp_data[srv_ind,]
            comp_hat_tmp <- comp_hat[srv_ind, ]

            # Reorganize and clean
            comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
            comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))

            comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
            comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))

            # Calculate pearson residual
            comp_tmp$comp_hat <- comp_hat_tmp$comp
            comp_tmp <- comp_tmp[which(comp_tmp$comp > 0),]
            comp_tmp <- comp_tmp[which(comp_tmp$comp_hat > 0),]

            # Calculate pearson
            comp_tmp$pearson <- (comp_tmp$comp - comp_tmp$comp_hat) / sqrt( ( comp_tmp$comp_hat * (1 - comp_tmp$comp_hat)) / comp_tmp$Sample_size)
            if(is.null(max_z)){
              max_pearson <- max(abs(comp_tmp$pearson), na.rm = TRUE)
            } else{
              max_pearson = max_z
            }


            sp <- fleet_control$Species[which(fleet_control$Fleet_code == srv[j])]
            nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

            plot(
              y = NA,
              x = NA,
              ylim = c(0, nages[sp] * 1.25),
              xlim = c(min(abs(comp_tmp$Year), na.rm = TRUE), max(abs(comp_tmp$Year), na.rm = TRUE) + right_adj),
              xlab = "Year",
              ylab = c("Pearson residual: Age comp", "Pearson residual: Length comp")[comp_type + 1],
              xaxt = "s"
            )


            # Legends
            legend("topleft", as.character(fleet_control$Fleet_name[srv[j]]), bty = "n", cex = 1)



            ############################
            # Legend
            #############################
            round = 0
            if(sum(seq(from = 1, to = max_pearson, length.out = 3) < 3) == 3){
              round = 1
            }

            # Positive
            x_loc <- c(max(abs(comp_tmp$Year), na.rm = T) - 6, max(abs(comp_tmp$Year), na.rm = T) - 3.5, max(abs(comp_tmp$Year), na.rm = T) - 1)
            symbols( x = x_loc , y = rep(nages[sp] * 1.1, 3) , circle = round(seq(from = 1, to = max_pearson, length.out = 3) , round), inches=0.20,add=T, bg = colvec[3])
            text(x = x_loc, y = rep(nages[sp] * 1.23, 3), labels = round(seq(from = 1, to = max_pearson, length.out = 3) , round))

            # Negative
            x_loc <- c(max(abs(comp_tmp$Year), na.rm = T) - 8.5, max(abs(comp_tmp$Year), na.rm = T) - 11, max(abs(comp_tmp$Year), na.rm = T) - 13.5)
            symbols( x = x_loc , y = rep(nages[sp] * 1.1, 3) , circle = -round(seq(from = -1, to = -max_pearson, length.out = 3) , round), inches=0.20,add=T, bg = NA)
            text(x = x_loc, y = rep(nages[sp] * 1.23, 3), labels = round(seq(from = -1, to = -max_pearson, length.out = 3) , round) )


            #############################
            # Plot comp pearson residuals
            #############################
            if(nrow(comp_tmp) > 0){

              # Get colors
              comp_tmp$colors <- NA
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 0, colvec[3], comp_tmp$colors) # Combined sex / 1 sex model
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 1, colvec[1], comp_tmp$colors) # Females
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 2, colvec[2], comp_tmp$colors) # Males

              # Joint sex
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age <= nages[sp], colvec[1], comp_tmp$colors) # Females
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age > nages[sp], colvec[2], comp_tmp$colors) # Males

              # Adjust years for joint
              comp_tmp$Year <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age > nages[sp], comp_tmp$Year + 0.2, comp_tmp$Year) # Males
              comp_tmp$Year <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age <= nages[sp], comp_tmp$Year - 0.2, comp_tmp$Year) # Females

              # Adjust age for joint data
              comp_tmp$age <- ifelse(comp_tmp$age > nages[sp], comp_tmp$age - nages[sp], comp_tmp$age)

              # Background colors
              comp_tmp$bg_colors <- ifelse(comp_tmp$pearson > 0, comp_tmp$colors, NA)

              # Plot
              symbols( x = c(0, comp_tmp$Year) , y = c(0, comp_tmp$age), circle = c(max_pearson, abs(comp_tmp$pearson)), inches=0.2,add=T, bg = comp_tmp$bg_colors, fg = comp_tmp$colors)
            }



            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }
  }
