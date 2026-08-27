# GOA arrowtooth flounder -- 2026 harvest projections.
#
# Catch advice for 2027 and 2028, worked out two ways for each of the two SAFE
# models:
#
#   1. Rceattle's own Tier 3 projection, run when the model was fitted.
#   2. SPM, the standard AFSC projection model, which also gives the seven
#      harvest scenarios.
#
# Run "2026 fit projection models.R" first; this script reads the saved fits.
#
# The two should agree closely in 2027. They drift apart later because SPM
# redraws recruitment 1000 times and reports the mean, while Rceattle carries a
# single trajectory forward.

# Install spmR from your local cloned folder
devtools::install("C:/Users/kalei.shotwell/Work/GitHub/spmR")

library(Rceattle)
library(spmR)     # remotes::install_github("afsc-assessments/spmR")

source("R/Functions/spm_bridge.R")

out_dir <- "2026/projections"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fits <- list("26.0" = readRDS("2026/models/mod_26_0.RDS"),
             "26.1" = readRDS("2026/models/mod_26_1.RDS"))
endyr <- fits[["26.0"]]$data_list$endyr        # 2026; advice is for 2027 and 2028

# Compile SPM once. Needs ADMB installed.
#if (!file.exists(file.path(out_dir, "spm"))) build_spm(out_dir)

# Project both models ----------------------------------------------------

for (model in names(fits)) {
  fit <- fits[[model]]

  # Each model gets its own folder of SPM inputs and outputs.
  run_dir <- file.path(out_dir, paste0("mod_", sub(".", "_", model, fixed = TRUE)))
  # Trick the spmR file.exists() check by creating an extensionless copy
  file.copy(
    from = file.path(run_dir, "spm.exe"),
    to = file.path(run_dir, "spm"),
    overwrite = TRUE
  )
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(out_dir, "spm"), file.path(run_dir, "spm"), overwrite = TRUE)
  Sys.chmod(file.path(run_dir, "spm"), "0755")

  write_spm_inputs(fit, run_dir)
  detail <- run_spm(run_dir)

  from_rceattle <- rceattle_exec_table(fit)
  from_spm      <- spm_exec_table(run_dir, detail, endyr)
  scenarios     <- tier3_scenario_table(detail, years = endyr + 1:2)

  # Put the two routes side by side, year by year.
  comparison <- data.frame(
    Quantity = from_rceattle$Quantity,
    from_rceattle[[2]], from_spm[[2]],
    from_rceattle[[3]], from_spm[[3]],
    check.names = FALSE
  )
  names(comparison) <- c("Quantity",
                         paste("Rceattle", endyr + 1), paste("SPM", endyr + 1),
                         paste("Rceattle", endyr + 2), paste("SPM", endyr + 2))

  write.csv(comparison, file.path(run_dir, "comparison.csv"), row.names = FALSE)
  write.csv(scenarios,  file.path(run_dir, "scenarios.csv"),  row.names = FALSE)

  cat("\n==== Model", model, "====\n")
  print(comparison, row.names = FALSE)
  cat("\nSeven harvest scenarios (SPM):\n")
  print(as.data.frame(scenarios), row.names = FALSE)
}
