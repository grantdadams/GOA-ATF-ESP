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

library(Rceattle)
library(spmR)     # remotes::install_github("afsc-assessments/spmR")
                  # or, from a local clone: devtools::install("path/to/spmR")

source("R/Functions/spm_bridge.R")

out_dir <- "2026/projections"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fits <- list("26.0" = readRDS("2026/models/mod_26_0.RDS"),
             "26.1" = readRDS("2026/models/mod_26_1.RDS"))
endyr <- fits[["26.0"]]$data_list$endyr        # 2026; advice is for 2027 and 2028

# Catch specified for the projection years, in tons. The fishery takes only a
# small share of the ABC, so the projection is run on the catch actually
# expected rather than the full ABC: the ABC in force (119,985 t) times the
# five-year average yield ratio of 0.1034184856. Both years take the same value
# because the second is a rollover.
#
# Without this the projection removes the whole ABC in 2027, which understates
# the 2028 population and so the 2028 ABC and OFL. The 2027 column is not
# affected either way.
#
# The yield ratio is the assessment author's, kept as given. For the record: it
# matches the mean catch for 2021-2025 (12,475 t) to within 0.5%, and sits 10%
# below the mean for 2022-2026 (13,844 t), so it appears to predate the 2026
# catch entering the workbook. Recomputing it on the newer window would raise
# it, and would need the ABC series, since the ratio is the mean of catch over
# ABC year by year rather than mean catch over one ABC.
specified_catch <- c("2027" = 12408.667, "2028" = 12408.667)

# The OFL specified for 2025 by the last assessment (2026 was 143,347 t). Used
# for Flimit, which the SARA file and SIS need: the F that would have taken that
# OFL, given where this model now puts the 2025 stock.
ofl_last_year <- 142832

# The SPM program. On a Mac or Linux this compiles it, which needs ADMB
# installed; on Windows put a prebuilt spm.exe in 2026/projections/ instead.
# run_spm() picks up whichever of the two is there, so nothing below changes.
#if (!file.exists(file.path(out_dir, spm_program()))) build_spm(out_dir)


# Project both models ----------------------------------------------------

for (model in names(fits)) {
  fit <- fits[[model]]

  # Each model gets its own folder of SPM inputs and outputs.
  run_dir <- file.path(out_dir, paste0("mod_", sub(".", "_", model, fixed = TRUE)))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  write_spm_inputs(fit, run_dir, fixed_catch = specified_catch)
  detail <- run_spm(run_dir)

  from_rceattle <- rceattle_exec_table(fit, fixed_catch = specified_catch)
  from_spm      <- spm_exec_table(run_dir, detail, endyr)
  scenarios     <- scenario_years(detail, fit)   # assessment year + all projection years

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
  cat("scenarios.csv covers", min(scenarios$Year), "to", max(scenarios$Year), "\n")
  print(comparison, row.names = FALSE)
  cat("\nSeven harvest scenarios (SPM):\n")
  print(head(scenarios, 14), row.names = FALSE)

  # Flimit for the SARA file and SIS.
  if (is.na(ofl_last_year)) {
    cat("\nFlimit: not calculated -- set ofl_last_year to the OFL specified for",
        endyr - 1, "\n")
  } else {
    cat("\nFlimit for", endyr - 1, ":", round(flimit(fit, ofl_last_year), 4),
        " (F35% is", round(fit$quantities$Flimit[1], 4), ")\n")
  }
}
