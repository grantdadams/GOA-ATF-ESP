## Script to run the projections for various models
library(tidyverse)
source("~/GitHub/GOApollock/R/proj_fns.R", echo=TRUE)

## Final 2022 model
abc_calc_tmb <- function(replist, datlist){
  ## New automated way using the new 'spm' package. Confirmed the
  ## same as manually doing it with 'proj' below.
  write_spm_inputs_tmb(replist, datlist, path = 'R/Functions/proj')
  prev_dir <- getwd()
  setwd('R/Functions/proj')
  system("spm")
  setwd(prev_dir)
  bf <- read.csv('R/Functions/proj/spm_detail.csv')
  exec_table <- get_exec_table_tmb(replist, datlist, bf)
  exec_tableF <- format_exec_table(exec_table)

  proj_scens <- bf %>% select(-Spp) %>%
    group_by(Alternative,Yr)%>%
    summarize_all('mean') %>% arrange(Alternative,Yr) %>% ungroup %>%
    as.data.frame()

  return(list(exec_table = exec_table, exec_tableF = exec_tableF, proj_scens = proj_scens))
}


write_spm_inputs_tmb <- function(replist, datlist, path=getwd()){
  ## delete old files in case it fails it'll stop there
  ff <- file.path(path, c("goa_atf.txt", "spm.dat", "spp_catch.dat"))
  trash <- lapply(ff, function(x) if(file.exists(x)) file.remove(x))
  ayr <- datlist$endyr
  write_spm_setup(ayr, path=path, catch=tail(datlist$cattot,1))
  write_spm_dynamics_tmb(replist, datlist, path)
  ## old unused file, just need dummy values so it runs
  ## tacpar <- readLines(tacpar.dat'); dput(tacpar)
  write.table(file=file.path(path,'tacpar.dat'),
              x=c("7", "6", " 2.60197 0.417 0.372 0.361 0.296 0.2733 0.125",
                  " -0.300856 -0.741664 -1.26797 -1.78614 -2.15198 -2.39595 -2.59939",
                  " -0.00703968 -0.956592 -1.48909 -2.09122 -2.26005 -2.31812 -2.45321",
                  " 0.372276 -1.35066 -1.76076 -2.22163 -2.31726 -2.28107 -2.34758",
                  " 0.505535 -1.19926 -1.68812 -2.39543 -2.40752 -2.5904 -2.58726",
                  " 0.368722 -1.59025 -1.79548 -2.67422 -2.61358 -2.41092 -2.70204",
                  " 0.308248 -1.78052 -2.23264 -2.85605 -3.39375 -3.05209 -2.94219",
                  " 0.180676 -1.7586 -1.80222 -3.09402 -2.2618 -3.43215 -3.06417"))
}



write_spm_dynamics_tmb <- function(replist, datlist, path){
  hindyrs <- dat_list$styr:dat_list$endyr
  ayr = dat_list$endyr

  x <- c()
  x <-
    c(x,paste("# proj input file written by write_spm_inputs on", Sys.time()))
  x <- c(x, "GOA_ATF")
  x <- c(x, "0 # Flag to tell if this is a SSL forage species")
  x <- c(x, "0 # Flag to Dorn's version of a constant buffer")
  x <- c(x, "1 # number of fisheries")
  x <- c(x, "2 # number of sexes")
  x <- c(x, paste(mean(tail(replist$F_spp[1,1:length(hindyrs)],5)), " # average F over last 5 years"))
  x <- c(x, "1 # Author's F multiplier (to MaxPermissible)")
  x <- c(x, "0.4 # ABC SPR")
  x <- c(x, "0.35 # OFL SPR")
  x <- c(x, "2 # Month of spawning")
  x <- c(x, "21 # Number of ages")
  x <- c(x, "1 # Fratio")

  x <- c(x, "# Natural mortality females")
  x <- c(x, paste(rowMeans(replist$M[1,1,,(length(hindyrs)-5):length(hindyrs)]), collapse=' '))
  x <- c(x, "# Natural mortality males")
  x <- c(x, paste(rowMeans(replist$M[1,2,,(length(hindyrs)-5):length(hindyrs)]), collapse=' '))
  x <- c(x, "# Maturity")
  x <- c(x, paste(dat_list$pmature[1,2:22], collapse=' '))
  x <- c(x, "# Spawning WAA")
  x <- c(x, paste(dat_list$wt %>%
                    dplyr::filter(Year <= dat_list$endyr & Year > dat_list$endyr-5 & Sex == 1) %>%
                    dplyr::select(paste0("Age", 1:21)) %>%
                    colMeans(), collapse=' '))
  x <- c(x, "# Fishery WAA")
  x <- c(x, paste(dat_list$wt %>%
                    dplyr::filter(Year <= dat_list$endyr & Year > dat_list$endyr-5 & Sex == 1) %>%
                    dplyr::select(paste0("Age", 1:21)) %>%
                    colMeans(), collapse=' '))

  x <- c(x, "# Fishery selex females, averaged over last 5 years")
  x <- c(x, paste(replist$sel[3,1,,1], collapse=' ')) # 2023 (nyrs+1) not time-varying
  x <- c(x, "# Fishery selex males, averaged over last 5 years")
  x <- c(x, paste(replist$sel[3,2,,1], collapse=' ')) # 2023 (nyrs+1) not time-varying

  x <- c(x, "# current year starting NAA females")
  x <- c(x, paste(replist$NByage[1,1,,length(hindyrs)], collapse=' '))
  x <- c(x, "# current year starting NAA males")
  x <- c(x, paste(replist$NByage[1,2,,length(hindyrs)], collapse=' '))

  ind <- which(hindyrs %in% 1978:(ayr-1))
  recs <- replist$R[,ind]
  ssb <- replist$biomassSSB[,ind-1]
  x <- c(x, paste(length(recs), " # num of recruits"))
  x <- c(x, "# Recruits from 1978 to this year -1 due to uncertain last year")
  x <- c(x, paste(recs, collapse=' '))
  x <- c(x, '# SSB from 1977 to this year-2 to match recruits')
  x <- c(x, paste(ssb, collapse=' '))
  writeLines(x, con=file.path(path, 'goa_atf.txt'))
}


#' Get table from projection output
#' @param replist A list as read in by \link{read_pk_rep}
#' @param bigfile A data frame of full projection scenario
#'   outputs
#' @export
get_exec_table_tmb <- function(replist, datlist, bigfile, maxABCratio=1){
  ayr <- tail(datlist$yrs,1)
  M <- c(.3,.3)
  means <- bigfile %>% filter(Yr > ayr & Yr <= ayr+2) %>%
    group_by(Alternative, Yr, Spp, SpNo) %>%
    summarize_all(mean) %>% ungroup
  sumbio <- replist[[113]][1:2]*1e6
  ssb <- filter(means,  Alternative==1) %>% pull(SSB) %>%
    round(0)
  if('B0' %in% names(means))
    bfrac <- filter(means, Alternative==1) %>% select(B0,B40,B35) %>%
    round(-3) %>% t
  else {warning("no B0 col"); bfrac=matrix(NA,3,2)}
  fofl <- filter(means,  Alternative==2) %>% pull(FOFL) %>% round(3)
  fabc <-filter(means, Alternative==1) %>% pull(F) %>% round(3)
  maxfabc <- fabc
  ofl <- filter(means,  Alternative==2) %>% pull(OFL) %>% round(0)
  maxabc <- filter(means,  Alternative==1) %>% pull(Catch) %>% round(0)
  abc <- maxabc*maxABCratio
  tab <- rbind(M, sumbio, ssb, bfrac, fofl, maxfabc, fabc,ofl,maxabc,abc)
  tab <- as.data.frame(tab) %>% setNames(c(ayr+1,ayr+2)) %>%
    cbind(name=row.names(tab),.)
  row.names(tab) <- NULL
  return(tab)
}

write_spm_setup <- function(ayr, path, catch){
  x <- c()
  x <- c(x,paste("# proj input file written by write_proj_inputs on",
                 Sys.time()))
  x <- c(x, " std # run name")
  x <- c(x, " 3 # tier")
  x <- c(x, "7 # number of alt scenarios")
  x <- c(x, paste(1:7, collapse=' '))
  x <- c(x, "1 # flag to set TAC=ABC (1 is true)")
  x <- c(x, "1 # SR type (1 ricker, 2 bholt")
  x <- c(x, "1 # proj rec form (default: 1= use obs mean/SD")
  x <- c(x, "0 # SR conditioning (0 means no)")
  x <- c(x, "0 # turn off prior")
  x <- c(x, "1 # flag to write bigfile")
  x <- c(x, "14 # number of proj years")
  x <- c(x, "1000 # number of simulations")
  x <- c(x, paste(ayr, " # begin year"))
  x <- c(x, "5 # number of years with specified catch")
  x <- c(x, "1 # numberof species")
  x <- c(x, "0 # OY min")
  x <- c(x, "2e6 # OY max")
  x <- c(x, "goa_atf.txt # input file name for biology")
  x <- c(x, "1 # ABC multipliers")
  x <- c(x, "1000 # pop scalars")
  x <- c(x, "0.75 # new alt 4 Fabc SPRs (unused?)")
  x <- c(x, "1 # num of TAC model categories")
  x <- c(x, "1 # TAC model indices")
  x <- c(x, paste(ayr,catch,"# Catch in each year starting w/ begining NAA")) #FIXME: add projected catch * yieldratio of 0.147
  writeLines(x, con=file.path(path, 'spm.dat'))
}
