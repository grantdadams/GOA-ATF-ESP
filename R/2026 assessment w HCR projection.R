# GOA arrowtooth flounder -- fit the two 2026 SAFE models for projection.
#
#   Model 26.0  Base. Sex-specific M fixed (females 0.2, males 0.35 per year).
#   Model 26.1  Alternative. Sex-specific M estimated.
#
# Both are fitted under the Tier 3 harvest control rule and saved to
# 2026/models/, where "2026 assessment SPM projections.R" picks them up. Takes a
# few minutes.

library(Rceattle)

dir.create("2026/models", showWarnings = FALSE, recursive = TRUE)

# Data ----
data_2026 <- Rceattle::read_data(file = "2026/data/GOA_arrowtooth_2026.xlsx")
data_2026$estDynamics <- 0

# The likelihood reads Log_sd as a CV, but the workbook stores it on the same
# scale as the observation, so divide through.
data_2026$index_data$Log_sd <- data_2026$index_data$Log_sd / data_2026$index_data$Observation

# Excel is storing the environmental columns as text, so they arrive as
# character and TMB rejects them. Turn them back into numbers.
env_columns <- setdiff(names(data_2026$env_data), "Year")
data_2026$env_data[env_columns] <- lapply(data_2026$env_data[env_columns], as.numeric)

# The workbook has Proj_F_proportion at 0 for every fleet, which leaves no
# fishing mortality to allocate in the projection. Arrowtooth has one fishery,
# so it takes all of it.
is_fishery <- data_2026$fleet_control$Fleet_type %in% c(1, "Fishery")
data_2026$fleet_control$Proj_F_proportion <- ifelse(is_fishery, 1, 0)

# Harvest control rule ----
# Amendment-56 Tier 3: FOFL = F35%, max FABC = F40%, both ramped down below
# B40% with alpha = 0.05, and no directed fishing below B20%.
tier3 <- Rceattle::build_hcr(
  HCR        = "NPFMC",
  DynamicHCR = FALSE,
  Ftarget    = 0.40,   # F40% -- max FABC
  Flimit     = 0.35,   # F35% -- FOFL
  Plimit     = 0,      # directed fishing prohibited below BX%
  Alpha      = 0.05
)

# * Model 26.0 -- fixed sex-specific M ----
mod_26_0 <- Rceattle::fit_mod(
  data_list    = data_2026,
  estimateMode = "Estimate",     # hindcast + HCR projection
  random_rec   = TRUE,           # recruitment deviations as random effects
  msmMode      = 0,              # single species
  initMode     = 1,
  HCR          = tier3,
  fit_control  = Rceattle::fit_control(verbose = 1, phase = TRUE)
)
saveRDS(mod_26_0, "2026/models/mod_26_0.RDS")

# * Model 26.1 -- estimated sex-specific M ----
mod_26_1 <- Rceattle::fit_mod(
  data_list    = data_2026,
  estimateMode = "Estimate",
  random_rec   = TRUE,
  msmMode      = 0,
  initMode     = 1,
  M1Fun        = Rceattle::build_M1(M1_model = 2),
  HCR          = tier3,
  fit_control  = Rceattle::fit_control(verbose = 1, phase = TRUE)
)
saveRDS(mod_26_1, "2026/models/mod_26_1.RDS")
