
# Install dependencies ----
install.packages("pacman")
install.packages("TMB", type = "source")
install.packages("Matrix", type = "source")
pacman::p_load(dplyr,
               ggplot2,
               MASS,
               oce,
               readxl,
               TMB,
               devtools,
               writexl,
               reshape2,
               gplots,
               tidyr,
               remotes,
               testthat,
               foreach,
               R.utils,
               rmarkdown,
               knitr,
               doParallel)
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# NOTE: both DSEM and dev-name-change will work, but use the dev-name-change for now
remotes::install_github("grantdadams/Rceattle", ref = "dev-DSEM") # Master branch was the branch previously used for the ATF assessment
remotes::install_github("grantdadams/Rceattle") # Master branch was the branch previously used for the ATF assessment
