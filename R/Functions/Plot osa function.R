#' Explore OSA residuals for multinomial composition data and
#' compare to Pearson
#' @param obs,exp,pearson the observed, expected and Pearson
#'   residual matrices with rows as years and columns as ages (or
#'   lengths)
#' @param index,years vectors giving the index of ages (or length) and years
#' @param index_label character value indicating 'age' or 'length bin' depending on comp type
#' @param stock,survey characters given the stock and survey,
#' used to create filename stock_survey.pdf
#' @param outpath folder name for output, e.g. 'figs'
#' @return returns nothing but creates a PDF file in the working
#' directory
#'
plot_osa_comps <- function(model, file = '', model_name = NULL){

  Fleets <- model$data_list$comp_data %>%
    dplyr::filter(Sample_size > 0) %>%
    dplyr::distinct(Fleet_name, Fleet_code, Age0_Length1, Sex)

  comp_obs <- model$data_list$comp_data
  comp_hat <- comp_obs %>%
    dplyr::select(-grep("Comp_", colnames(comp_data)))
  comp_hat <- cbind(comp_hat, model$quantities$comp_hat)


  for(flt in 1:nrow(Fleets)){

    # - Get comps
    comp_obs_temp <- comp_obs %>%
      dplyr::filter(Fleet_code == Fleets$Fleet_code[flt] &
                      Age0_Length1 == Fleets$Age0_Length1[flt] &
                      Sex == Fleets$Sex[flt])
    comp_hat_temp <- comp_hat %>%
      dplyr::filter(Fleet_code == Fleets$Fleet_code[flt] &
                      Age0_Length1 == Fleets$Age0_Length1[flt] &
                      Sex == Fleets$Sex[flt])
    n_eff = comp_obs_temp %>%
      dplyr::pull(Sample_size)

    cols_keep <- grep("Comp_", colnames(comp_obs_temp))
    comp_obs_info = comp_obs_temp %>%
      dplyr::select(-cols_keep)
    comp_obs_temp = comp_obs_temp %>%
      dplyr::select(cols_keep)
    comp_hat_temp = comp_hat_temp %>%
      dplyr::select(cols_keep)

    o <- round(n_eff*comp_obs_temp/rowSums(comp_obs_temp, na.rm = TRUE),0)
    p <- comp_hat_temp/rowSums(comp_hat_temp, na.rm = TRUE)


    # - Get info
    year <- comp_obs_info$Year
    if(Fleets$Age0_Length1[flt] == 0){
      index <- model$data_list$minage:(model$data_list$minage+model$data_list$nages-1)
      sex <- rep(Fleets$Sex, length(index))

      if(Fleets$Sex[flt] == 3){
        sex <- rep(c("Females", "Males"), each = length(index))
        index <- c(index, index)
      }
    }

    if(Fleets$Age0_Length1[flt] == 1){
      index <- 1:(model$data_list$nlengths)
      sex <- rep(Fleets$Sex, length(index))

      if(Fleets$Sex[flt] == 3){
        sex <- rep(c("Females", "Males"), each = length(index))
        index <- c(index, index)
      }
    }
    index <- paste0(sex, "_", index)

    # - Remove columns with no obs
    o <- o[, 1:length(index)]
    p <- p[, 1:length(index)]

    res <-  compResidual::resMulti(t(o), t(p))

    if(!all(is.finite(res))){
      warning("failed to calculate OSA residuals for ", Fleets$Fleet_name[flt])
    }else{

      # * Plot full OSA ---
      filename <- paste0(file, paste0("_osa_full_plot_fleet_code_", Fleets$Fleet_code[flt], ".png"))
      png(
        file = filename ,
        width = 7,# 169 / 25.4,
        height = 6.5,# 150 / 25.4,

        units = "in",
        res = 300
      )
      plot(res)
      dev.off()


      # * Plot OSA bubbles ----
      mat <- t(matrix(res, nrow=nrow(res), ncol=ncol(res)))
      dimnames(mat) <- list(year=year, index=index[-1])
      reslong <- reshape2::melt(mat, value.name='resid')
      reslong$sex <- as.character(sapply(reslong$index, function(x) strsplit(as.character(x), "_")[[1]])[1,])
      reslong$index <- as.numeric(sapply(reslong$index, function(x) strsplit(as.character(x), "_")[[1]])[2,])

      g1 <- ggplot(reslong, aes(year, index, size=abs(resid),
                                color=resid>0)) + geom_point() +
        ggtitle(paste0(model_name, ': OSA w/o index 1')) +
        facet_wrap(~sex) +  ylab(Fleets$Fleet_name[flt]) +
        theme_classic()
      print(g1)

      filename <- paste0(file, paste0("_osa_bubbles_fleet_code_", Fleets$Fleet_code[flt], ".png"))
      ggsave(filename = filename ,
             g1,
             width = 7,# 169 / 25.4,
             height = 6.5,# 150 / 25.4,
             units = "in",
             dpi = 300)

      # ## ind is age/len bin to drop
      # for(ind in 1:length(index)){
      #   ## assumes **last column** dropped so put it there
      #   index2 <- index[-ind]
      #   o2 <- cbind(o[,-ind], o[,ind])
      #   p2 <- cbind(p[,-ind], p[,ind])
      #   res <- resMulti(t(o2), t(p2))
      #   ## not sure why these fail sometimes?
      #   if(!all(is.finite(res))) {warning('failed when ind=',ind); break}
      #   mat <- t(matrix(res, nrow=nrow(res), ncol=ncol(res)))
      #   dimnames(mat) <- list(year=years, index=index2)
      #   reslong <- reshape2::melt(mat, value.name='resid')
      #   g <- ggplot(reslong, aes(year, index, size=abs(resid),
      #                            color=resid>0)) + geom_point()+
      #     ggtitle(paste0('OSA w/o ', index_label, ' ', index[ind])) + ylim(range(index)) +
      #     ylab(index_label)
      #   print(g)
      # }
    }
  }
}
