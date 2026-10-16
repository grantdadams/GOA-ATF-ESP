
# Install dependencies ----
install.packages("remotes")
remotes::install_github("afsc-assessments/Rceattle@5.64.0")

# NOTE: both DSEM and main will work
remotes::install_github("grantdadams/Rceattle", ref = "dev-DSEM")
remotes::install_github("afsc-assessments/Rceattle@5.64.0")
