# Run the standard projection model (SPM) from a fitted Rceattle model.
#
# SPM is the ADMB projection model behind the Tier 3 executive summary table and
# the seven harvest scenarios (github.com/afsc-assessments/spmR). Three
# functions here are all you need:
#
#   write_spm_inputs(fit, dir)   write SPM's three input files from a fit
#   run_spm(dir)                 run SPM and return its results
#   spm_exec_table(...)          the Tier 3 table from those results
#
# plus rceattle_exec_table(fit) for the same table out of Rceattle's own
# projection, and build_spm(dir) to compile SPM once.
#
# Four things about SPM's file format are easy to get wrong and produce no error
# if you do. All four were checked against SPM's source and its bundled GOA
# arrowtooth example:
#
#   1. Recruitment is FEMALE recruits, not total. Passing total recruitment
#      doubles B100%, B40% and B35%.
#   2. Spawning month is a month number, and SPM uses (month - 1) / 12 as the
#      fraction of the year before spawning. Rceattle's spawn_month is already
#      that fraction times 12, so SPM's month = Rceattle's + 1.
#   3. Spawning weight-at-age is one row (females). Mortality, maturity, fishery
#      weight, selectivity and numbers-at-age are two rows each, females then
#      males. Drop a row and every field after it shifts.
#   4. Numbers, weights, catch, recruitment and SSB are written in Rceattle's
#      own units (thousands of fish, kg, tons) with SPM's scalar set to 1, so
#      SPM's output comes back in tons.


# Write SPM's input files ------------------------------------------------
#
# fit         a fitted Rceattle model
# dir         directory to write into
# stock_name  name SPM prints on its output
# rec_years   years of recruitment SPM draws future recruitment from.
#             Default 1978 to endyr - 1; the terminal year is not yet informed
#             by data. Also sets B100%, which scales with mean recruitment.
# nproj       projection years
# nsims       simulated recruitment trajectories
# alt4_spr    SPR rate for harvest scenario 4. Check this against how your SAFE
#             chapter words scenario 4.
# author_f    author's F as a fraction of max FABC; 1 is max permissible
# avg_years   terminal years to average M, weight-at-age and selectivity over

write_spm_inputs <- function(fit,
                             dir,
                             stock_name = "GOA_ATF",
                             rec_years  = NULL,
                             nproj      = 14,
                             nsims      = 1000,
                             alt4_spr   = 0.75,
                             author_f   = 1,
                             avg_years  = 5) {

  q <- fit$quantities
  d <- fit$data_list
  endyr <- d$endyr
  nages <- d$nages[1]
  ages  <- paste0("Age", 1:nages)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # This code assumes a two-sex stock, ages starting at 1, and one fishery.
  fishery <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type %in% c(1, "Fishery")]
  if (d$nsex[1] != 2)   stop("Written for a two-sex stock; nsex is ", d$nsex[1])
  if (d$minage[1] != 1) stop("Written for ages starting at 1; minage is ", d$minage[1])
  if (length(fishery) != 1) stop("Expected one fishery, found ", length(fishery))

  # Average the last few years of anything that could vary through time.
  recent <- as.character((endyr - avg_years + 1):endyr)
  avg <- function(x) rowMeans(matrix(x, nrow = nages))

  m_female   <- avg(q$M_at_age[1, 1, , recent])
  m_male     <- avg(q$M_at_age[1, 2, , recent])
  sel_female <- avg(q$sel_at_age[fishery, 1, , recent])
  sel_male   <- avg(q$sel_at_age[fishery, 2, , recent])
  avg_f      <- mean(q$F_spp[1, recent])   # drives harvest scenario 3

  # Selectivity goes across on Rceattle's scale rather than rescaled to peak at
  # 1. SPM's F is then the same number as Rceattle's Ftarget and Flimit, and
  # every catch and biomass is unaffected either way.

  # Weight-at-age, in kg. Rceattle can hold several series; take the one the
  # spawning biomass and the fishery each use.
  weight <- function(index, sex) {
    w <- d$weight[d$weight$Wt_index == index & d$weight$Sex == sex, ]
    if (length(unique(w$Year)) > 1) w <- w[w$Year %in% (endyr - avg_years + 1):endyr, ]
    colMeans(w[, ages, drop = FALSE])
  }
  wt_index  <- d$ssb_wt_index[1]
  fsh_index <- d$fleet_control$Weight_index[d$fleet_control$Fleet_code == fishery]

  # Maturity. SPM reads a male ogive as well but never uses it -- spawning
  # biomass is females only -- so the female ogive is written twice.
  maturity <- as.numeric(d$maturity[1, ages])

  # Recruitment, and the spawning biomass that produced it. Age-1 fish in year
  # y came from spawning in year y - 1.
  if (is.null(rec_years)) rec_years <- 1978:(endyr - 1)
  if (min(rec_years) <= d$styr) stop("rec_years must start after ", d$styr)
  if (max(rec_years) > endyr)   stop("rec_years must end by ", endyr)
  recruits <- q$N_at_age[1, 1, 1, as.character(rec_years)]      # females only
  ssb      <- q$ssb[1, as.character(rec_years - 1)]

  # Numbers at age at the start of the first projection year. The age-1 slot is
  # set to mean recruitment so that year is drawn the same way as later years.
  n_female <- q$N_at_age[1, 1, , as.character(endyr + 1)]
  n_male   <- q$N_at_age[1, 2, , as.character(endyr + 1)]
  n_female[1] <- n_male[1] <- mean(recruits)

  row <- function(x) paste(format(as.numeric(x), scientific = FALSE, trim = TRUE),
                           collapse = " ")

  # * Biology file ----
  writeLines(c(
    paste("# written from an Rceattle fit on", Sys.Date()),
    "#runname",        stock_name,
    "#ssl_spp",        "0",
    "#Dorn_buffer",    "0",
    "#nfsh",           "1",
    "#nsex",           "2",
    "#avgF5yr",        row(avg_f),
    "#F40_mult",       row(author_f),
    "#spr_abc",        "0.4",              # Tier 3: max FABC is F40%
    "#spr_msy",        "0.35",             # Tier 3: FOFL is F35%
    "#sp_mo",          row(d$spawn_month[1] + 1),
    "#nages",          row(nages),
    "#Frat",           "1",
    "#M females",      row(m_female),
    "#M males",        row(m_male),
    "#pmat females",   row(maturity),
    "#pmat males",     row(maturity),
    "#wtage_sp",       row(weight(wt_index, 1)),
    "#wtage_fsh females", row(weight(fsh_index, 1)),
    "#wtage_fsh males",   row(weight(fsh_index, 2)),
    "#sel females",    row(sel_female),
    "#sel males",      row(sel_male),
    "#N females",      row(n_female),
    "#N males",        row(n_male),
    "#nrec",           row(length(recruits)),
    paste0("#R female recruits ", min(rec_years), "-", max(rec_years)), row(recruits),
    paste0("#SSB ", min(rec_years) - 1, "-", max(rec_years) - 1),       row(ssb)
  ), file.path(dir, "goa_atf.dat"))

  # * Setup file ----
  writeLines(c(
    paste("# written from an Rceattle fit on", Sys.Date()),
    "#rn",               "std",
    "#Tier",             "3",
    "#nalts",            "7",
    "#alts",             paste(1:7, collapse = " "),
    "#tac_flag",         "1",
    "#srr_type",         "1",
    # Draw future recruitment from the observed mean and variance, so the
    # stock-recruit curve above does not enter the projection.
    "#srr_form",         "1",
    "#srr_conditioning", "0",
    "#srr_reserved",     "0",
    "#spm_detail_flag",  "1",
    "#nprj_yrs",         as.character(nproj),
    "#nsims",            as.character(nsims),
    # Start the year after the assessment: the terminal year's catch is fitted.
    "#beg_yr",           as.character(endyr + 1),
    "#nyrs_fixed_catch", "0",
    "#nspp",             "1",
    "#OY_min",           "0",
    # Optimum yield cap, set high so it cannot bind on a single stock.
    "#OY_max",           "2000000",
    "#data_files",       "goa_atf.dat",
    "#ABC_mults",        "1",
    "#scalars",          "1",
    "#alt4_spr",         as.character(alt4_spr),
    "#nTAC_cat",         "1",
    "#TACind",           "1",
    "#fixed_catch",      paste(endyr, 0, "# placeholder, not read")
  ), file.path(dir, "spm.dat"))

  # * ABC-to-TAC file ----
  # Only used when several stocks share a TAC. One stock takes TAC = ABC, so
  # these are placeholders that let SPM read its input.
  writeLines(c("7", "6", "1 1 1 1 1 1 1", rep("0 0 0 0 0 0 0", 7)),
             file.path(dir, "tacpar.dat"))

  invisible(dir)
}


# Compile and run SPM ----------------------------------------------------

# Compile SPM from the template in the spmR package. Needs ADMB installed.
# Only has to be done once.
build_spm <- function(dir) {
  template <- system.file("admb", "spm.tpl", package = "spmR")
  if (template == "") {
    stop("Install spmR first: remotes::install_github('afsc-assessments/spmR')")
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir <- normalizePath(dir)

  build <- file.path(tempdir(), "spmbuild")
  dir.create(build, showWarnings = FALSE)
  file.copy(template, file.path(build, "spm.tpl"), overwrite = TRUE)

  here <- setwd(build)
  on.exit(setwd(here))
  message <- system2("admb", "spm", stdout = TRUE, stderr = TRUE)
  if (!file.exists("spm")) stop("ADMB build failed:\n", paste(message, collapse = "\n"))

  file.copy("spm", file.path(dir, "spm"), overwrite = TRUE)
  Sys.chmod(file.path(dir, "spm"), "0755")
  invisible(file.path(dir, "spm"))
}


# Run SPM on the input files in `dir` and return one row per simulation, year
# and harvest scenario.
run_spm <- function(dir) {
  if (!file.exists(file.path(dir, "spm"))) {
    stop("No spm program in ", dir, ". Run build_spm('", dir, "') first.")
  }
  here <- setwd(dir)
  on.exit(setwd(here))

  # Delete the last run's output so a failure this time cannot be read as a
  # success. SPM exits with an error code even when it works, so the output
  # file appearing is what tells us the run finished.
  unlink(c("spm_detail.csv", "spm_summary.csv"))
  system2("./spm", stdout = "spm_run.log", stderr = "spm_run.log")
  if (!file.exists("spm_detail.csv")) {
    stop("SPM failed; see ", file.path(dir, "spm_run.log"))
  }
  read.csv("spm_detail.csv")
}


# Executive summary tables -----------------------------------------------

exec_rows <- c("Total biomass (t)", "Female spawning biomass (t)",
               "B100% (t)", "B40% (t)", "B35% (t)",
               "FOFL", "max FABC", "FABC",
               "OFL (t)", "max ABC (t)", "ABC (t)")

# Stack the eleven rows into a table, rounding F rates to 4 decimals and
# biomasses and catches to whole tons.
exec_table <- function(rows, years) {
  x <- do.call(rbind, rows)
  is_rate <- exec_rows %in% c("FOFL", "max FABC", "FABC")
  x[is_rate, ]  <- round(x[is_rate, ], 4)
  x[!is_rate, ] <- round(x[!is_rate, ])
  out <- data.frame(Quantity = exec_rows, x, check.names = FALSE)
  setNames(out, c("Quantity", years))
}


# Tier 3 table from SPM. `alt` is the harvest scenario: 2 is the author's F,
# which is the same as scenario 1 when author_f = 1.
spm_exec_table <- function(dir, detail, endyr, alt = 2) {
  years <- c(endyr + 1, endyr + 2)

  # Reference points, computed once for the whole projection.
  ref <- read.csv(file.path(dir, "spm_summary.csv"))
  refpt <- function(name) {
    value <- as.numeric(ref$value[match(name, ref$variable)])
    c(value, value)                     # same in both year columns
  }

  # Everything else is a mean over the simulations, one value per year.
  yearly <- function(column) {
    sapply(years, function(y) mean(detail[[column]][detail$Alt == alt & detail$Year == y]))
  }

  exec_table(list(
    yearly("Tot_biom"),
    yearly("SSB"),
    refpt("SSB_100"),
    refpt("SSB_40"),
    refpt("SSB_ofl"),
    refpt("F_ofl"),
    refpt("F_abc"),
    refpt("F_abc"),
    yearly("OFL"),
    yearly("ABC"),
    yearly("ABC")
  ), years)
}


# Tier 3 table from Rceattle's own projection under build_hcr(HCR = "NPFMC").
# One trajectory rather than a mean over simulations, so it answers a slightly
# different question than the SPM table, but the first year should agree.
rceattle_exec_table <- function(fit) {
  q <- fit$quantities
  d <- fit$data_list
  years <- c(d$endyr + 1, d$endyr + 2)

  fishery <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type %in% c(1, "Fishery")]
  wt_index <- d$fleet_control$Weight_index[d$fleet_control$Fleet_code == fishery]

  sb0  <- q$SB0[1, as.character(years[1])]   # unfished female spawning biomass
  fabc <- q$Ftarget[1]                       # F40%
  fofl <- q$Flimit[1]                        # F35%

  # Catch the projected population takes at a given F, summed over sex and age.
  # At Ftarget this reproduces Rceattle's own projected catch exactly, which is
  # what makes the OFL row comparable to the ABC row.
  catch_at <- function(year, f) {
    year <- as.character(year)
    f_at_age <- f * q$sel_at_age[fishery, , , year]
    z <- f_at_age + q$M_at_age[1, , , year]
    sum(q$weight_hat[wt_index, , , year] * f_at_age / z * (1 - exp(-z)) *
          q$N_at_age[1, , , year])
  }

  by_year <- function(f) sapply(years, f)

  exec_table(list(
    by_year(function(y) q$biomass[1, as.character(y)]),
    by_year(function(y) q$ssb[1, as.character(y)]),
    c(sb0, sb0),
    c(0.40 * sb0, 0.40 * sb0),
    c(0.35 * sb0, 0.35 * sb0),
    c(fofl, fofl),
    c(fabc, fabc),
    c(fabc, fabc),
    by_year(function(y) catch_at(y, fofl)),
    by_year(function(y) catch_at(y, fabc)),
    by_year(function(y) catch_at(y, fabc))
  ), years)
}
