recdevs <- c(rev(bridging_model_3$estimated_params$init_dev[1,1:20]), bridging_model_3$estimated_params$rec_dev[1,1:47])
R0 = bridging_model_3$estimated_params$rec_pars[1,1]
m=mean(recdevs)
R0=R0+m
recdevs=recdevs-m

R0
recdevs

bridging_model_3$quantities$R[1,1:47]

fdev <- bridging_model_3$estimated_params$F_dev[3,]
m=mean(fdev)
Fhat=bridging_model_3$estimated_params$ln_mean_F[3]+m
fdev=fdev-m

Fhat
fdev
