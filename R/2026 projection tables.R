# GOA arrowtooth flounder -- 2026 harvest projections.
#
# Two routes to the same Amendment-56 Tier 3 advice, for models 26.0 and 26.1:
#
#   1. Rceattle's internal NPFMC Tier 3 harvest control rule. One deterministic
#      trajectory forward from the fitted model, using the model's own SPR
#      reference points, selectivity and projected recruitment.
#   2. The standard projection model (SPM, afsc-assessments/spmR). The
#      operational AFSC projection: 1000 simulated recruitment trajectories and
#      the seven standard harvest scenarios, following the approach Cole
#      Monnahan uses for GOA pollock.
#
# Run "2026 fit projection models.R" first; this script reads the saved fits.
#
# The two routes should agree closely on the reference points and on the first
# projection year. They diverge later because SPM redraws recruitment from the
# historical mean and variance while Rceattle carries its own projection
# forward, and because SPM reports means over simulations.

library(Rceattle)
library(dplyr)
library(spmR)     # remotes::install_github("afsc-assessments/spmR")

source("R/Functions/spm_bridge.R")

out_dir <- "2026/projections"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mods <- list("26.0" = readRDS("2026/models/mod_26_0.RDS"),
             "26.1" = readRDS("2026/models/mod_26_1.RDS"))
ayr <- mods[["26.0"]]$data_list$endyr        # 2026; advice is for 2027 and 2028

# SPM executable ---------------------------------------------------------
# spmR ships the ADMB template but no macOS binary, so build it once. Needs
# ADMB on the PATH.
if (!file.exists(file.path(out_dir, "spm"))) build_spm(out_dir)


# Projections ------------------------------------------------------------

internal <- list()
spm_tab  <- list()
spm_alts <- list()

for (nm in names(mods)) {
  m <- mods[[nm]]

  # * Route 1: Rceattle's internal Tier 3 HCR ----
  internal[[nm]] <- rceattle_exec_table(m)

  # * Route 2: SPM ----
  run_dir <- file.path(out_dir, paste0("mod_", sub("\\.", "_", nm)))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(out_dir, "spm"), file.path(run_dir, "spm"), overwrite = TRUE)
  Sys.chmod(file.path(run_dir, "spm"), "0755")

  write_spm_inputs(m, path = run_dir, spp_file = "goa_atf.dat",
                   run_name = "GOA_ATF", nsims = 1000, nproj = 14)
  detail <- run_spm(run_dir)

  spm_tab[[nm]]  <- spm_exec_table(run_dir, detail, endyr = ayr)
  spm_alts[[nm]] <- tier3_scenario_table(detail, years = ayr + 1:2)

  write.csv(spm_tab[[nm]],  file.path(run_dir, "exec_table_spm.csv"), row.names = FALSE)
  write.csv(internal[[nm]], file.path(run_dir, "exec_table_rceattle.csv"), row.names = FALSE)
  write.csv(spm_alts[[nm]], file.path(run_dir, "scenario_table.csv"), row.names = FALSE)
}


# Comparison -------------------------------------------------------------
# The SPR reference points are the check that the bridge is faithful: SPM
# recomputes F40% and F35% from the M, maturity, weight-at-age, selectivity and
# spawning timing written into its input file, so agreement with Rceattle's own
# Ftarget and Flimit means all of those crossed over correctly.

compare <- function(nm) {
  a <- internal[[nm]]; b <- spm_tab[[nm]]
  data.frame(Quantity = a$Quantity,
             Rceattle = a[[2]], SPM = b[[2]],
             Rceattle_yr2 = a[[3]], SPM_yr2 = b[[3]],
             check.names = FALSE) |>
    setNames(c("Quantity",
               paste0(c("Rceattle ", "SPM "), ayr + 1),
               paste0(c("Rceattle ", "SPM "), ayr + 2)))
}

for (nm in names(mods)) {
  cat("\n==== Model", nm, "====\n")
  print(compare(nm), row.names = FALSE)
  cat("\nSeven harvest scenarios (SPM):\n")
  print(as.data.frame(spm_alts[[nm]]), row.names = FALSE)
}

write.csv(compare("26.0"), file.path(out_dir, "comparison_26_0.csv"), row.names = FALSE)
write.csv(compare("26.1"), file.path(out_dir, "comparison_26_1.csv"), row.names = FALSE)
