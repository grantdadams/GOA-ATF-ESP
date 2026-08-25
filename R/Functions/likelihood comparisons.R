

# Likelihoods ----
#
# UNVERIFIED ROW NUMBERS -- check before this table goes in a document.
# The Value column below addresses `jnll_comp` by bare row number. Rceattle now
# fixes those rows with the JnllRow enum (ceattle.cpp), and the order has moved
# since these numbers were written: rows 1-6 are now Index data, Catch data,
# Composition data, CAAL data, Non-parametric selectivity, Selectivity deviates.
# On that order "Catch" is reading Composition data and "Survey" is reading
# Catch data. Index by `rownames(model$quantities$jnll_comp)` instead, once it
# is settled which component each label is meant to name -- "Slectivity" could
# be either of the two selectivity rows.
#
# `index_ln_q` / `R_ln_sd` were renamed to `index_log_q` / `R_log_sd`; until
# that was fixed this function errored on exp(NULL).
get_atf_ll <- function(model){

  years <- model$data_list$styr:model$data_list$endyr
  nyrs <- length(years)

  projyears <- model$data_list$styr:model$data_list$projyr
  nprojyrs <- length(projyears)

  atf_ll <- data.frame( Component = c(
    "Catch",
    "Fishery length",
    "Survey",
    "Survey age",
    "PENALTIES",
    "Recruitment deviations",
    "Slectivity",
    "Joint nll",
    "Marginal nll",
    "N parm",
    "q",
    "Mean R",
    "sigmaR",
    "F40%",
    "2023 Biomass",
    "2023 SSB",
    "B0",
    "B40",
    "SB0",
    "SB40",
    "ABC"
  ),
  Value = c(
    model$quantities$jnll_comp[3,3],
    model$quantities$jnll_comp[4,3],
    model$quantities$jnll_comp[2,1],
    model$quantities$jnll_comp[4,1],
    NA, # Penalties
    model$quantities$jnll_comp[11,1] + model$quantities$jnll_comp[12,1],
    model$quantities$jnll_comp[5,3],
    sum(model$quantities$jnll_comp),
    ifelse(is.null(model$opt$objective),  model$quantities$jnll, model$opt$objective),
    length(model$obj$par),
    exp(model$estimated_params$index_log_q)[1],
    mean(model$quantities$R[,1:nyrs]),
    exp(model$estimated_params$R_log_sd)[1],
    NA,
    model$quantities$biomass[1,nyrs],
    model$quantities$ssb[1,nyrs],
    model$quantities$biomass[1,nprojyrs],
    model$quantities$biomass[1,nprojyrs] * 0.4,
    model$quantities$ssb[1,nprojyrs],
    model$quantities$ssb[1,nprojyrs] * 0.4,
    NA
  )
  )

}
