
# Install dependencies ----
install.packages("remotes")
remotes::install_github("grantdadams/Rceattle")

# NOTE: both DSEM and main will work
remotes::install_github("grantdadams/Rceattle", ref = "dev-DSEM")
remotes::install_github("grantdadams/Rceattle")
