# Run the standard projection model (SPM) from a fitted Rceattle model.
#
# SPM is the ADMB projection model behind the Tier 3 executive summary table and
# the seven harvest scenarios (github.com/afsc-assessments/spmR).
#
# In the order you use them:
#
#   build_spm(dir)                    compile SPM. Once, and only on Mac or
#                                     Linux -- on Windows drop in an spm.exe
#   write_spm_inputs(fit, dir, ...)   write SPM's three input files from a fit
#   run_spm(dir)                      run SPM and return its results
#   spm_exec_table(dir, ...)          the Tier 3 executive summary table
#   scenario_years(detail, fit)       every scenario and year, for the
#                                     projection tables
#   rceattle_exec_table(fit, ...)     the same Tier 3 table, but from
#                                     Rceattle's own projection instead
#
# The last one is the cross-check: the two models share no code, so when they
# agree the translation below is right.
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


# A fixed_catch must be named with the years it applies to, running
# consecutively from the first projection year. Both routes below check it:
# without the names the SPM file loses its year column, and the Rceattle table
# would quietly go back to removing the whole ABC. Returns the years.
check_fixed_catch <- function(fixed_catch, endyr) {
  if (is.null(names(fixed_catch))) {
    stop("fixed_catch needs to be named with its years, as in ",
         "c(\"", endyr + 1, "\" = 12408.667, \"", endyr + 2, "\" = 12408.667)")
  }
  years <- suppressWarnings(as.numeric(names(fixed_catch)))
  expected <- endyr + seq_along(fixed_catch)
  if (anyNA(years) || !all(years == expected)) {
    stop("fixed_catch years must run consecutively from ", endyr + 1,
         "; got ", paste(names(fixed_catch), collapse = ", "))
  }
  years
}


# Write SPM's input files ------------------------------------------------
#
# fit           a fitted Rceattle model
# dir           directory to write into
# stock_name    name SPM prints on its output
# fixed_catch   catch in tons to force in the first projection years, named by
#               year, e.g. c("2027" = 12408.667, "2028" = 12408.667). Use this
#               when the fishery is not expected to take the full ABC. SPM
#               applies it to every scenario except 6 and 7, which have to use
#               FOFL by definition. NULL lets every scenario take its own ABC.
# rec_years     years of recruitment SPM draws future recruitment from.
#               Default 1978 to endyr - 1; the terminal year is not yet informed
#               by data. Also sets B100%, which scales with mean recruitment.
# nproj         projection years
# nsims         simulated recruitment trajectories
# alt4_spr      SPR rate for harvest scenario 4. Check this against how your
#               SAFE chapter words scenario 4.
# author_f      author's F as a fraction of max FABC; 1 is max permissible
# avg_years     terminal years to average M, weight-at-age and selectivity over

write_spm_inputs <- function(fit,
                             dir,
                             stock_name  = "GOA_ATF",
                             fixed_catch = NULL,
                             rec_years   = NULL,
                             nproj       = 14,
                             nsims       = 1000,
                             alt4_spr    = 0.75,
                             author_f    = 1,
                             avg_years   = 5) {

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

  # Average the last few years of any biology that could vary through time.
  # Time-invariant for arrowtooth, so the window makes no difference here.
  recent <- as.character((endyr - avg_years + 1):endyr)
  avg <- function(x) rowMeans(matrix(x, nrow = nages))

  m_female   <- avg(q$M_at_age[1, 1, , recent])
  m_male     <- avg(q$M_at_age[1, 2, , recent])
  sel_female <- avg(q$sel_at_age[fishery, 1, , recent])
  sel_male   <- avg(q$sel_at_age[fishery, 2, , recent])

  # Average F for harvest scenario 3, on the same scale as the selectivity
  # below. Guidelines 4.11.3 define this window as assessment year - 5 through
  # assessment year - 1, so it stops the year BEFORE endyr: the assessment
  # year's own catch is an extrapolation to 31 December, not a landed total.
  f_years <- as.character((endyr - 5):(endyr - 1))
  avg_f   <- mean(q$F_spp[1, f_years])

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
  # Specified catch: SPM reads one "year catch" line per year, starting at the
  # first projection year and running consecutively.
  if (is.null(fixed_catch)) {
    n_fixed <- 0
    catch_lines <- paste(endyr, 0, "# placeholder, not read")
  } else {
    n_fixed <- length(fixed_catch)
    catch_lines <- paste(check_fixed_catch(fixed_catch, endyr), fixed_catch)
  }

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
    "#nyrs_fixed_catch", as.character(n_fixed),
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
    "#fixed_catch",      catch_lines
  ), file.path(dir, "spm.dat"))

  # * ABC-to-TAC file ----
  # Only used when several stocks share a TAC. One stock takes TAC = ABC, so
  # these are placeholders that let SPM read its input.
  writeLines(c("7", "6", "1 1 1 1 1 1 1", rep("0 0 0 0 0 0 0", 7)),
             file.path(dir, "tacpar.dat"))

  invisible(dir)
}


# Compile and run SPM ----------------------------------------------------

# What the SPM program is called here: spm.exe on Windows, spm elsewhere.
spm_program <- function() {
  if (.Platform$OS.type == "windows") "spm.exe" else "spm"
}

# Find the SPM program: in the run folder if it is there, otherwise in the
# projections folder above it.
find_spm <- function(dir) {
  here  <- file.path(dir, spm_program())
  above <- file.path(dirname(dir), spm_program())
  if (file.exists(here)) return(here)
  if (file.exists(above)) return(above)
  stop("No ", spm_program(), " in ", dir, " or ", dirname(dir),
       ". On a Mac or Linux run build_spm() to compile it; on Windows copy in ",
       "a prebuilt spm.exe.")
}

# Compile SPM from the template in the spmR package. Needs ADMB installed.
# Only has to be done once. On Windows, use a prebuilt spm.exe instead.
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
  program <- spm_program()
  if (!file.exists(program)) {
    stop("ADMB build failed:\n", paste(message, collapse = "\n"))
  }

  file.copy(program, file.path(dir, program), overwrite = TRUE)
  Sys.chmod(file.path(dir, program), "0755")
  invisible(file.path(dir, program))
}


# Run SPM on the input files in `dir` and return one row per simulation, year
# and harvest scenario.
run_spm <- function(dir) {
  program <- normalizePath(find_spm(dir))
  here <- setwd(dir)
  on.exit(setwd(here))

  # Delete the last run's output so a failure this time cannot be read as a
  # success. SPM exits with an error code even when it works, so the output
  # file appearing is what tells us the run finished.
  unlink(c("spm_detail.csv", "spm_summary.csv"))
  system2(program, stdout = "spm_run.log", stderr = "spm_run.log")
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

# Build the table. `rows` is a list of eleven pairs, one per row of exec_rows
# above and in the same order, each holding the first year's value then the
# second. F rates are rounded to 4 decimals, biomasses and catches to tons.
exec_table <- function(rows, years) {
  x <- do.call(rbind, rows)     # eleven pairs -> an 11 x 2 table
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

  # Everything else is a mean over the 1000 simulations, one value per year.
  yearly <- function(column) {
    sapply(years, function(y) {
      rows <- detail$Alt == alt & detail$Year == y
      mean(detail[[column]][rows])
    })
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



scenario_names <- c("Maximum permissible ABC", "Author-specified ABC",
                    "Average recent F", "Alternative SPR rate", "No fishing",
                    "OFL threshold determination", "Status-determination ramp")

scenario_columns <- c("Alt", "Scenario", "Year", "Catch", "ABC", "OFL",
                      "SSB", "B_B35", "Tot_biom", "Rec", "F")

# Every projection year of every harvest scenario, averaged over the
# simulations, plus the assessment year itself from the fit. One row per
# scenario and year, so the projection tables in the SAFE chapter can be
# pivoted out of it however they need to be laid out.
scenario_years <- function(detail, fit) {
  columns <- c("Catch", "ABC", "OFL", "SSB", "Tot_biom", "Rec", "F", "B35")

  out <- aggregate(detail[columns], by = detail[c("Alt", "Year")], FUN = mean)
  out$Scenario <- scenario_names[out$Alt]
  out$B_B35 <- out$SSB / out$B35            # B35 is the same in every simulation
  out <- out[, scenario_columns]

  # The assessment year uses SPM's B35 too, so the whole B/B35% column is on one
  # basis. Rceattle's own B35 differs by about 0.2%, because it averages
  # recruitment over a slightly different set of years.
  out <- rbind(assessment_year(fit, b35 = unique(detail$B35)[1]), out)
  out <- out[order(out$Alt, out$Year), ]
  rownames(out) <- NULL
  out
}


# The assessment year, taken from the fit. SPM starts the year after, so this
# row is history: the same under every scenario, with the catch already landed
# and no ABC or OFL of its own.
assessment_year <- function(fit, b35) {
  q <- fit$quantities
  d <- fit$data_list
  year <- d$endyr
  chr <- as.character(year)

  # The catch actually landed, not the model's fit to it. The two agree to
  # about 0.1%, but a projection table should report the observed number.
  fleet <- fishery_fleet(fit)
  catch <- d$catch_data$Catch[d$catch_data$Fleet_code == fleet &
                                d$catch_data$Year == year]
  if (length(catch) != 1) {
    stop("Expected one ", year, " catch row for the fishery, found ", length(catch))
  }

  data.frame(
    Alt      = seq_along(scenario_names),
    Scenario = scenario_names,
    Year     = year,
    Catch    = catch,
    ABC      = NA_real_,
    OFL      = NA_real_,
    SSB      = q$ssb[1, chr],
    B_B35    = q$ssb[1, chr] / b35,
    Tot_biom = q$biomass[1, chr],
    Rec      = q$R[1, chr],
    F        = q$F_spp[1, chr]
  )
}


# Tier 3 table from Rceattle's own projection under build_hcr(HCR = "NPFMC").
# One trajectory rather than a mean over simulations, so it answers a slightly
# different question than the SPM table, but the first year should agree.
#
# fixed_catch works as it does for SPM: catch in tons named by year. Rceattle's
# projection always takes the full ABC, so when a catch is specified the second
# year's population has to be worked out again from the first year's. The
# reference points and the F rates do not change -- the Tier 3 ramp keys on the
# previous year's spawning biomass, and arrowtooth spawn before the fishery
# opens, so 2027 spawning biomass is already set whatever the 2027 catch is.
rceattle_exec_table <- function(fit, fixed_catch = NULL) {
  q <- fit$quantities
  d <- fit$data_list
  years <- c(d$endyr + 1, d$endyr + 2)

  sb0  <- q$SB0[1, as.character(years[1])]   # unfished female spawning biomass
  fabc <- q$Ftarget[1]                       # F40%
  fofl <- q$Flimit[1]                        # F35%

  # Numbers at age, sex by age, in each of the two years.
  numbers <- list(q$N_at_age[1, , , as.character(years[1])],
                  q$N_at_age[1, , , as.character(years[2])])

  # If a catch is specified for the first year, redo the second year's numbers.
  # Recruitment carries over untouched: it is set by the first year's spawning
  # biomass, which the first year's catch cannot change.
  if (!is.null(fixed_catch)) {
    check_fixed_catch(fixed_catch, d$endyr)
    numbers[[2]] <- age_one_year(fit, years[1], catch = fixed_catch[[1]],
                                 recruits = numbers[[2]][, 1])
  }

  # Each quantity, for the first year then the second. Both years are worked
  # out the same way, from that year's population.
  biomass <- c(total_biomass(fit, years[1], numbers[[1]]),
               total_biomass(fit, years[2], numbers[[2]]))
  ssb     <- c(spawning_biomass(fit, years[1], numbers[[1]]),
               spawning_biomass(fit, years[2], numbers[[2]]))
  ofl     <- c(catch_at_f(fit, years[1], fofl, numbers[[1]]),
               catch_at_f(fit, years[2], fofl, numbers[[2]]))
  abc     <- c(catch_at_f(fit, years[1], fabc, numbers[[1]]),
               catch_at_f(fit, years[2], fabc, numbers[[2]]))

  both <- function(x) c(x, x)     # a reference point: the same in both years

  exec_table(list(
    biomass,
    ssb,
    both(sb0),
    both(0.40 * sb0),
    both(0.35 * sb0),
    both(fofl),
    both(fabc),
    both(fabc),
    ofl,
    abc,
    abc
  ), years)
}


# The pieces the table above is built from. Each takes numbers at age as a
# 2 x nages matrix, females then males, so a population other than the one
# Rceattle projected can be passed in.

# Which fleet is the fishery, and which weight-at-age series it uses.
fishery_fleet <- function(fit) {
  fc <- fit$data_list$fleet_control
  fc$Fleet_code[fc$Fleet_type %in% c(1, "Fishery")]
}

# Catch in tons taken at fishing mortality `f`, by the Baranov equation. At
# Rceattle's own projected F this reproduces its projected catch exactly.
catch_at_f <- function(fit, year, f, numbers) {
  q <- fit$quantities
  year <- as.character(year)
  fleet <- fishery_fleet(fit)
  wt_index <- fit$data_list$fleet_control$Weight_index[
    fit$data_list$fleet_control$Fleet_code == fleet]

  f_at_age <- f * q$sel_at_age[fleet, , , year]
  z <- f_at_age + q$M_at_age[1, , , year]
  sum(q$weight_hat[wt_index, , , year] * f_at_age / z * (1 - exp(-z)) * numbers)
}

# The fishing mortality that takes `catch` tons.
f_for_catch <- function(fit, year, catch, numbers) {
  most <- catch_at_f(fit, year, 5, numbers)
  if (catch > most) {
    stop("A catch of ", round(catch), " t is more than the ", year,
         " population could yield at any plausible F (at most ", round(most), " t).")
  }
  uniroot(function(f) catch_at_f(fit, year, f, numbers) - catch,
          interval = c(0, 5), tol = 1e-10)$root
}

total_biomass <- function(fit, year, numbers) {
  wt <- fit$quantities$weight_hat[fit$data_list$pop_wt_index[1], , , as.character(year)]
  sum(wt * numbers)
}

# Females only. Arrowtooth spawn at the start of the year (spawn_month 0), so
# no mortality is applied before spawning.
spawning_biomass <- function(fit, year, numbers) {
  d <- fit$data_list
  if (d$spawn_month[1] != 0) {
    stop("Written for a stock that spawns at the start of the year.")
  }
  wt <- fit$quantities$weight_hat[d$ssb_wt_index[1], , , as.character(year)]
  maturity <- as.numeric(d$maturity[1, paste0("Age", 1:d$nages[1])])
  sum(numbers[1, ] * maturity * wt[1, ])
}

# Flimit, for the SARA file and SIS (Guidelines section 4.15): the F that would
# have produced a catch equal to last year's OFL, given where this model now
# thinks the stock was in that year.
#
# `ofl` is the OFL that was specified for that year by the previous assessment,
# so it has to be supplied -- it is not something this model produces. Give it
# in tons, like everything else here; the spreadsheet the Guidelines link to
# wants kt instead. The year defaults to endyr - 1, the most recent complete
# year, which is the year the overfishing determination is made for.
#
# Reading it: Flimit below F35% means this model now sees a bigger stock in that
# year than the assessment that set the OFL did, so less F was needed to take
# it. Above F35% means the opposite.
flimit <- function(fit, ofl, year = fit$data_list$endyr - 1) {
  numbers <- fit$quantities$N_at_age[1, , , as.character(year)]
  f_for_catch(fit, year, ofl, numbers)
}


# Numbers at age one year on, if `catch` tons are taken in `year`. Survivors
# age by one, the oldest two ages pool into the plus group, and `recruits`
# fills age 1.
age_one_year <- function(fit, year, catch, recruits) {
  q <- fit$quantities
  nages <- fit$data_list$nages[1]
  year <- as.character(year)

  f <- f_for_catch(fit, year, catch, q$N_at_age[1, , , year])
  z <- f * q$sel_at_age[fishery_fleet(fit), , , year] + q$M_at_age[1, , , year]
  survivors <- q$N_at_age[1, , , year] * exp(-z)

  next_year <- matrix(0, nrow = 2, ncol = nages)
  next_year[, 2:nages] <- survivors[, 1:(nages - 1)]
  next_year[, nages]   <- next_year[, nages] + survivors[, nages]
  next_year[, 1]       <- recruits
  next_year
}
