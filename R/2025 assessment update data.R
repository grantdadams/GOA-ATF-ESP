# Update the 2025 dataset with new data
# See vignette fro creating a data-set from scratch
# https://grantdadams.github.io/Rceattle/articles/data-without-excel.html
# https://github.com/grantdadams/Rceattle/blob/main/examples/Build_data_without_excel.R

library(Rceattle)
library(dplyr)

# Data ----
data_2025_c_la_wa_ae <- Rceattle::read_data( file = "Data/GOA_25.1.1_arrowtooth_single_species_1977-2024_new_data_clean_LA_WA_AE.xlsx")
data_2025_c_la_wa_ae$estDynamics = 0
data_2025_c_la_wa_ae$index_data$Log_sd <- data_2025_c_la_wa_ae$index_data$Log_sd/data_2025_c_la_wa_ae$index_data$Observation


# Update data ----
# * Controls ----
data_2025_c_la_wa_ae$endyr <- 2026 # New end year


# * Catch ----
# - Add in 2026
catch_data <- data_2025_c_la_wa_ae$catch_data[1,]
catch_data$Year <- data_2025_c_la_wa_ae$endyr
catch_data$Catch <- 99999 # @kalei add in 2026 catch
data_2025_c_la_wa_ae$catch_data <- rbind(data_2025_c_la_wa_ae$catch_data , catch_data)


# * Comp ----
# - Age
comp_data <- data_2025_c_la_wa_ae$comp_data %>%
  dplyr::group_by(Fleet_code, Age0_Length1) %>%
  slice(n()) %>%
  dplyr::mutate(Year = data_2025_c_la_wa_ae$endyr) %>%
  as.data.frame()

# - For length set Age0_Length1 = 1

# - Combine
data_2025_c_la_wa_ae$comp_data <- rbind(data_2025_c_la_wa_ae$comp_data, comp_data) %>%
  dplyr::arrange(Fleet_code, Age0_Length1, Year)


# * Index data ----
index_data <- data.frame(
  Fleet_name = "ATF_bottom_trawl", Fleet_code = 1,
  Species = 1,
  Year = 2025, Month = 0,
  Selectivity_block = 1, Q_block = 1,
  Observation = 99999, Log_sd = 99999)

data_2025_c_la_wa_ae$index_data <- rbind(data_2025_c_la_wa_ae$index_data, index_data)


# Save data ----
write_data(data_2025_c_la_wa_ae, "Data/GOA_26_arrowtooth_new_data_clean_LA_WA_AE.xlsx")

