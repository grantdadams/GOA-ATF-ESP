

# Likelihoods ----
get_atf_ll <- function(model){

  atf_ll <- data.frame( Component = c(
    "Catch",
    "Fishery length",
    "Survey",
    "Survey age",
    "PENALTIES",
    "Recruitment deviations",
    "Slectivity",
    "F dev",
    "Joint nll",
    "Marginal nll",
    "N parm",
    "q",
    "Mean R",
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
    model$quantities$jnll_comp[2,3],
    model$quantities$jnll_comp[3,3],
    model$quantities$jnll_comp[1,1],
    model$quantities$jnll_comp[3,1],
    NA, # Penalties
    model$quantities$jnll_comp[11,1] + model$quantities$jnll_comp[12,1],
    model$quantities$jnll_comp[5,3],
    model$quantities$jnll_comp[13,3],
    model$quantities$jnll,
    ifelse(is.null(model$opt$objective),  model$quantities$jnll, model$opt$objective),
    length(model$obj$par),
    exp(model$estimated_params$ln_srv_q)[1],
    mean(model$quantities$R[,length(1977:2023)]),
    NA,
    model$quantities$biomass[1,length(1977:2023)],
    model$quantities$biomassSSB[1,length(1977:2023)],
    model$quantities$biomass[1,length(1977:2050)],
    model$quantities$biomass[1,length(1977:2050)] * 0.4,
    model$quantities$biomassSSB[1,length(1977:2050)],
    model$quantities$biomassSSB[1,length(1977:2050)] * 0.4,
    NA
  )
  )

}
