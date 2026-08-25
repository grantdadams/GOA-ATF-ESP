#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <goa_atf.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  pad_evalout = new ofstream("atf_gen2.mcmc.out");       //wasgen;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_fut.allocate("styr_fut");
  endyr_fut.allocate("endyr_fut");
  assess.allocate("assess");
  nsurv.allocate("nsurv");
cout<<"nsurv"<<nsurv<<std::endl;
  median_rec.allocate("median_rec");
  first_age.allocate("first_age");
  last_age.allocate("last_age");
  first_length.allocate("first_length");
  last_length.allocate("last_length");
  nages_read.allocate("nages_read");
cout<<"nages_read"<<nages_read<<std::endl;
 nages=last_age-first_age+1;          // # of ages in the model  
cout<<"nages"<<nages<<std::endl;
  nsurv_aged.allocate("nsurv_aged");
cout<<"nsurv_aged"<<nsurv_aged<<std::endl;    
  aa.allocate(1,nages);
  yy.allocate(styr,endyr);
 for (i=styr;  i<=endyr; i++) yy(i)=i;
 for (i=1; i<=nages;i++) aa(i)=i;
  phase_F40.allocate("phase_F40");
  phase_logistic_sel.allocate(1,2*nsurv,"phase_logistic_sel");
  phase_fishery_sel.allocate(1,2,"phase_fishery_sel");
  phase_alphabeta.allocate("phase_alphabeta");
  q_Phase.allocate(1,nsurv,"q_Phase");
  phase_selcoffs.allocate("phase_selcoffs");
cout<<"phase_selcoffs"<<phase_selcoffs<<std::endl;         
  nselages.allocate("nselages");
  nselages_srv.allocate(1,nsurv,"nselages_srv");
  fishsel_LB_f.allocate(1,2,"fishsel_LB_f");
  fishsel_LB_m.allocate(1,2,"fishsel_LB_m");
  fishsel_UB_f.allocate(1,2,"fishsel_UB_f");
  fishsel_UB_m.allocate(1,2,"fishsel_UB_m");
  fishsel_prior_f.allocate(1,2,"fishsel_prior_f");
  fishsel_prior_m.allocate(1,2,"fishsel_prior_m");
  nsel_params.allocate(1,nsurv,"nsel_params");
cout<<"nsel_params"<<nsel_params<<std::endl;
  sel_prior_f.allocate(1,2*nsurv,"sel_prior_f");
  sel_prior_m.allocate(1,2*nsurv,"sel_prior_m");
cout<<"sel_prior_f"<<sel_prior_f<<std::endl;  
  sel_LB_f.allocate(1,2*nsurv,"sel_LB_f");
  sel_LB_m.allocate(1,2*nsurv,"sel_LB_m");
  sel_UB_f.allocate(1,2*nsurv,"sel_UB_f");
  sel_UB_m.allocate(1,2*nsurv,"sel_UB_m");
  sel1_desc_prior_f.allocate(1,2,"sel1_desc_prior_f");
  sel1_desc_prior_m.allocate(1,2,"sel1_desc_prior_m");
  sel1_desc_LB_f.allocate(1,2,"sel1_desc_LB_f");
  sel1_desc_LB_m.allocate(1,2,"sel1_desc_LB_m");
  sel1_desc_UB_f.allocate(1,2,"sel1_desc_UB_f");
  sel1_desc_UB_m.allocate(1,2,"sel1_desc_UB_m");
  nlen.allocate("nlen");
  nobs_fish.allocate("nobs_fish");
  yrs_fish.allocate(1,nobs_fish,"yrs_fish");
cout<<"yrs_fish"<<yrs_fish<<std::endl; 
  nsamples_fish.allocate(1,2,1,nobs_fish,"nsamples_fish");
  nobs_srv.allocate(1,nsurv,"nobs_srv");
  yrs_srv.allocate(1,nsurv,1,nobs_srv,"yrs_srv");
  nobs_srv_length.allocate(1,nsurv,"nobs_srv_length");
cout<<"nobs_srv_length"<<nobs_srv_length<<std::endl; 
  yrs_srv_length.allocate(1,nsurv,1,nobs_srv_length,"yrs_srv_length");
cout<<"yrs_srv_length"<<yrs_srv_length<<std::endl;
  nsamples_srv_length_fem.allocate(1,nsurv,1,nobs_srv_length,"nsamples_srv_length_fem");
  nsamples_srv_length_mal.allocate(1,nsurv,1,nobs_srv_length,"nsamples_srv_length_mal");
  obs_p_fish.allocate(1,2,1,nobs_fish,1,nlen,"obs_p_fish");
cout<<"obs_p_fish"<<obs_p_fish(1)<<std::endl;  
  obs_p_srv_length_fem.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"obs_p_srv_length_fem");
  obs_p_srv_length_mal.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"obs_p_srv_length_mal");
  catch_bio.allocate(styr,endyr,"catch_bio");
cout<<"catch_bio"<<catch_bio<<std::endl;
  obs_srv.allocate(1,nsurv,1,nobs_srv,"obs_srv");
cout<<"obs_srv"<<obs_srv<<std::endl;
  obs_srv_sd.allocate(1,nsurv,1,nobs_srv,"obs_srv_sd");
  wt.allocate(1,2,1,nages,"wt");
cout<<"wt"<<wt<<std::endl;  
  maturity.allocate(1,nages,"maturity");
cout<<"maturity"<<maturity<<std::endl;
  lenage.allocate(1,2,1,nages,1,nlen,"lenage");
cout<<"lenage"<<lenage<<std::endl; 
if(assess==1){nyrs_temps = nobs_srv(1);}  
if(assess==2) {nyrs_temps = 33;}//nobs_srv(1); change back for BSAI assessment
if(assess==3){nyrs_temps = nobs_srv(1);}
cout<<"nyrs_temps"<<nyrs_temps<<std::endl; 
  bottom_temps.allocate(1,nyrs_temps,"bottom_temps");
cout<<"nyrs_temps"<<nyrs_temps<<std::endl;
  monot_sel.allocate("monot_sel");
cout<<"monot_sel"<<monot_sel<<std::endl;
  wt_like.allocate(1,8,"wt_like");
cout<<"wt_like"<<wt_like<<std::endl;               
  nobs_srv_age.allocate(1,nsurv_aged,"nobs_srv_age");
  yrs_srv_age.allocate(1,nsurv_aged,1,nobs_srv_age,"yrs_srv_age");
  nsamples_srv_age.allocate(1,nsurv_aged,1,2,1,nobs_srv_age,"nsamples_srv_age");
  obs_p_srv_age_fem.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"obs_p_srv_age_fem");
  obs_p_srv_age_mal.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"obs_p_srv_age_mal");
  M.allocate(1,2,"M");
cout<<"M"<<M<<std::endl;            
  M_f.allocate(1,nages);
  M_m.allocate(1,nages);
 for (i=1; i<= nages; i++) { M_f(i)=M(1); M_m(i)=M(2); }
  offset_const.allocate("offset_const");
  q_Lower_bound.allocate(1,nsurv,"q_Lower_bound");
  q_Upper_bound.allocate(1,nsurv,"q_Upper_bound");
  q_surv_prior_mean.allocate(1,nsurv,"q_surv_prior_mean");
  nparams_srv.allocate(1,nsurv,"nparams_srv");
cout<<"q_surv_prior_mean"<<q_surv_prior_mean<<std::endl;  
  mean_log_rec_prior.allocate("mean_log_rec_prior");
  log_avg_fmort_prior.allocate("log_avg_fmort_prior");
 cout<<"the fmort prior "<<log_avg_fmort_prior<<endl;
  like_wght.allocate(1,7,"like_wght");
  fpen_mult.allocate(1,2,"fpen_mult");
  catch_err.allocate("catch_err");
  fmort_boundL.allocate("fmort_boundL");
  fmort_boundH.allocate("fmort_boundH");
 cout<<"fmort bound H is  "<<fmort_boundH<<endl;
  agerr_matrix.allocate(1,nages_read,1,nages_read,"agerr_matrix");
  obs_p_fish_norm.allocate(1,2,1,nobs_fish,1,nlen);
  eof_marker.allocate("eof_marker");
 cout<<"the universal answer is "<<eof_marker<<endl;
  cv_srv.allocate(1,nsurv,1,nobs_srv);
   styr_rec=styr-nages+1;
   if(nselages>nages) 
   {nselages=nages;  
   cout<<"Warning selectivity: is set to be estimated on more ages than are in the model."<<std::endl;  }
   for (i=1; i<= nsurv; i++){
   if(nselages_srv(i)>nages) nselages_srv(i)=nages;
   }
   //calculate cv for surveys
   for (int j=1;j<=nsurv;j++){
   for (i=1;i<=nobs_srv(j);i++){ 
   cv_srv(j,i)=obs_srv_sd(j,i)/(double)obs_srv(j,i); }} 
   //change weights to tons
   wt=wt*.001;
  catch_err_like=.5/(catch_err*catch_err);
  obs_sexr.allocate(1,nobs_fish);
  obs_sexr_srv_2.allocate(1,nsurv,1,nobs_srv_length);
  pred_sexr.allocate(styr,endyr);
}

void model_parameters::initializationfunction(void)
{
  F40.set_initial_value(.20);
  F35.set_initial_value(.21);
  F30.set_initial_value(.23);
  mean_log_rec.set_initial_value(mean_log_rec_prior);
  log_avg_fmort.set_initial_value(log_avg_fmort_prior);
  q_surv.set_initial_value(q_surv_prior_mean);
  fmort_dev.set_initial_value(0.00001);
  fishsel_params_f.set_initial_value(fishsel_prior_f);
  fishsel_params_m.set_initial_value(fishsel_prior_m);
  srv_params_f.set_initial_value(sel_prior_f);
  srv_params_m.set_initial_value(sel_prior_m);
  srv1desc_params_f.set_initial_value(sel1_desc_prior_f);
  srv1desc_params_m.set_initial_value(sel1_desc_prior_m);
  alpha.set_initial_value(1.);
  beta.set_initial_value(0.);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  q_surv.allocate(1,nsurv,q_Lower_bound,q_Upper_bound,q_Phase,"q_surv");
  fishsel_params_f.allocate(1,2,fishsel_LB_f,fishsel_UB_f,phase_fishery_sel,"fishsel_params_f");
  fishsel_params_m.allocate(1,2,fishsel_LB_m,fishsel_UB_m,phase_fishery_sel,"fishsel_params_m");
  srv_params_f.allocate(1,2*nsurv,sel_LB_f,sel_UB_f,phase_logistic_sel,"srv_params_f");
  srv_params_m.allocate(1,2*nsurv,sel_LB_m,sel_UB_m,phase_logistic_sel,"srv_params_m");
  srv1desc_params_f.allocate(1,2,sel1_desc_LB_f,sel1_desc_UB_f,phase_fishery_sel,"srv1desc_params_f");
  srv1desc_params_m.allocate(1,2,sel1_desc_LB_m,sel1_desc_UB_m,phase_fishery_sel,"srv1desc_params_m");
  alpha.allocate(phase_alphabeta,"alpha");
  beta.allocate(phase_alphabeta,"beta");
  mean_log_rec.allocate(1,"mean_log_rec");
  rec_dev.allocate(styr_rec,endyr,-15,15,2,"rec_dev");
  log_avg_fmort.allocate(1,"log_avg_fmort");
  fmort_dev.allocate(styr,endyr,fmort_boundL,fmort_boundH,1,"fmort_dev");
  log_selcoffs_fish.allocate(1,2,1,nselages,phase_selcoffs,"log_selcoffs_fish");
  sexr_param_fish.allocate(1.0,1.0,-5,"sexr_param_fish");
  F40.allocate(0.01,1.,phase_F40,"F40");
  F35.allocate(0.01,1.,phase_F40,"F35");
  F30.allocate(0.01,1.,phase_F40,"F30");
  log_sel_fish.allocate(1,2,1,nages,"log_sel_fish");
  #ifndef NO_AD_INITIALIZE
    log_sel_fish.initialize();
  #endif
  sel.allocate(1,2,1,nages,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  sel_srv.allocate(1,2,1,nsurv,1,nages,"sel_srv");
  #ifndef NO_AD_INITIALIZE
    sel_srv.initialize();
  #endif
  avgsel_fish.allocate(1,2,"avgsel_fish");
  #ifndef NO_AD_INITIALIZE
    avgsel_fish.initialize();
  #endif
  popn.allocate(1,2,styr,endyr,"popn");
  #ifndef NO_AD_INITIALIZE
    popn.initialize();
  #endif
  totn_srv.allocate(1,nsurv,1,2,styr,endyr,"totn_srv");
  #ifndef NO_AD_INITIALIZE
    totn_srv.initialize();
  #endif
  temp1.allocate(1,nages,"temp1");
  #ifndef NO_AD_INITIALIZE
    temp1.initialize();
  #endif
  temp2.allocate(1,nages,"temp2");
  #ifndef NO_AD_INITIALIZE
    temp2.initialize();
  #endif
  explbiom.allocate(styr,endyr,"explbiom");
  #ifndef NO_AD_INITIALIZE
    explbiom.initialize();
  #endif
  pred_bio.allocate(styr,endyr,"pred_bio");
  #ifndef NO_AD_INITIALIZE
    pred_bio.initialize();
  #endif
  fspbio.allocate(styr,endyr,"fspbio");
  pred_srv.allocate(1,nsurv,styr,endyr,"pred_srv");
  #ifndef NO_AD_INITIALIZE
    pred_srv.initialize();
  #endif
  pred_p_fish.allocate(1,2,styr,endyr,1,nlen,"pred_p_fish");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish.initialize();
  #endif
  pred_p_srv_age_fem.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"pred_p_srv_age_fem");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_age_fem.initialize();
  #endif
  pred_p_srv_age_mal.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"pred_p_srv_age_mal");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_age_mal.initialize();
  #endif
  pred_p_srv_len_fem.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"pred_p_srv_len_fem");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_len_fem.initialize();
  #endif
  pred_p_srv_len_mal.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"pred_p_srv_len_mal");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_len_mal.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  natage.allocate(1,2,styr,endyr,1,nages,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  totalbiomass.allocate(styr,endyr,"totalbiomass");
  catage.allocate(1,2,styr,endyr,1,nages,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  Z.allocate(1,2,styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,2,styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,2,styr,endyr,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  fmort.allocate(styr,endyr,"fmort");
  #ifndef NO_AD_INITIALIZE
    fmort.initialize();
  #endif
  rbar.allocate("rbar");
  #ifndef NO_AD_INITIALIZE
  rbar.initialize();
  #endif
  surv.allocate(1,2,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  offset.allocate(1,7,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  catch_like.allocate("catch_like");
  #ifndef NO_AD_INITIALIZE
  catch_like.initialize();
  #endif
  sexr_like.allocate("sexr_like");
  #ifndef NO_AD_INITIALIZE
  sexr_like.initialize();
  #endif
  age_like.allocate(1,nsurv_aged,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  length_like.allocate(1,nsurv+1,"length_like");
  #ifndef NO_AD_INITIALIZE
    length_like.initialize();
  #endif
  sel_like.allocate(1,4,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  surv_like.allocate(1,nsurv,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  endbiom.allocate("endbiom");
  depletion.allocate("depletion");
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  tmp.allocate("tmp");
  #ifndef NO_AD_INITIALIZE
  tmp.initialize();
  #endif
  pred_sexr.allocate(styr,endyr,"pred_sexr");
  #ifndef NO_AD_INITIALIZE
    pred_sexr.initialize();
  #endif
  sigmar.allocate("sigmar");
  #ifndef NO_AD_INITIALIZE
  sigmar.initialize();
  #endif
  ftmp.allocate("ftmp");
  #ifndef NO_AD_INITIALIZE
  ftmp.initialize();
  #endif
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  SBF30.allocate("SBF30");
  #ifndef NO_AD_INITIALIZE
  SBF30.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,4,1,nages,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  nage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"nage_future");
  #ifndef NO_AD_INITIALIZE
    nage_future.initialize();
  #endif
  fspbiom_fut.allocate(1,4,styr_fut,endyr_fut,"fspbiom_fut");
  F_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"F_future");
  #ifndef NO_AD_INITIALIZE
    F_future.initialize();
  #endif
  Z_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"Z_future");
  #ifndef NO_AD_INITIALIZE
    Z_future.initialize();
  #endif
  S_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"S_future");
  #ifndef NO_AD_INITIALIZE
    S_future.initialize();
  #endif
  catage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"catage_future");
  #ifndef NO_AD_INITIALIZE
    catage_future.initialize();
  #endif
  avg_rec_dev_future.allocate("avg_rec_dev_future");
  #ifndef NO_AD_INITIALIZE
  avg_rec_dev_future.initialize();
  #endif
  avg_F_future.allocate(1,4,"avg_F_future");
  #ifndef NO_AD_INITIALIZE
    avg_F_future.initialize();
  #endif
  catch_future.allocate(1,3,styr_fut,endyr_fut,"catch_future");
  future_biomass.allocate(1,4,styr_fut,endyr_fut,"future_biomass");
  explbiom_fut.allocate(styr_fut,endyr_fut,"explbiom_fut");
  #ifndef NO_AD_INITIALIZE
    explbiom_fut.initialize();
  #endif
  maxsel_fish.allocate("maxsel_fish");
  #ifndef NO_AD_INITIALIZE
  maxsel_fish.initialize();
  #endif
  maxsel_srv.allocate(1,nsurv,"maxsel_srv");
  #ifndef NO_AD_INITIALIZE
    maxsel_srv.initialize();
  #endif
  mlike.allocate("mlike");
  #ifndef NO_AD_INITIALIZE
  mlike.initialize();
  #endif
  qlike.allocate("qlike");
  #ifndef NO_AD_INITIALIZE
  qlike.initialize();
  #endif
  flike.allocate("flike");
  #ifndef NO_AD_INITIALIZE
  flike.initialize();
  #endif
  qtime.allocate(styr,endyr,"qtime");
  #ifndef NO_AD_INITIALIZE
    qtime.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  obs_mean_sexr=0.34;  //initial value for avg proportion of male population estimated from shelf surveys; calculated below
  obs_SD_sexr=0.0485;  //initial value for standard deviation of mean male population proportion: calculated below
  
  for(i=1; i<=nobs_fish;i++)
  {
    obs_sexr(i) = sum(obs_p_fish(1,i))/sum(obs_p_fish(1,i) + obs_p_fish(2,i)); 
  }  
  for(i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{    
    	obs_sexr_srv_2(i,j)=sum(obs_p_srv_length_mal(i,j)/
    	      (sum(obs_p_srv_length_mal(i,j))+sum(obs_p_srv_length_fem(i,j))));
	}
  }
  obs_mean_sexr=mean(obs_sexr_srv_2); //previously was just estimated from shelf survey data so kept that here.
  obs_SD_sexr=std_dev(obs_sexr_srv_2(1)); 
 //Compute offset for multinomial and length bin proportions
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_
  offset.initialize();
  if(assess<3){
  for (i=1; i <= nobs_fish; i++)
  {
  double sumtot ;  
  sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i)); 
  obs_p_fish_norm(1,i) = obs_p_fish(1,i) / sumtot; 
  obs_p_fish_norm(2,i) = obs_p_fish(2,i) / sumtot; 
  for(k=1; k<=2;k++) 
  {
    offset(1) -= nsamples_fish(k,i)*obs_p_fish_norm(k,i) * log(obs_p_fish_norm(k,i)+offset_const); //this multiplies elements together then sums them up.
  } 
  } 
  }
  if(assess==3)
  {
	for (i=1; i <= nobs_fish; i++)
  {
  double sumtot ;
  sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i));
  obs_p_fish_norm(1,i) /= sum(obs_p_fish(1,i));
  obs_p_fish_norm(2,i) /= sum(obs_p_fish(2,i));
  for(k=1; k<=2;k++)
    offset(1) -= nsamples_fish(k,i)*obs_p_fish_norm(k,i) * log(obs_p_fish_norm(k,i)+.0001);
  }
  }   
 
  //this loops over all surveys and makes sure all proportions sum to 1.
  for(i=1;i<=nsurv;i++){
	for(j=1;j<=nobs_srv_length(i);j++){    
		double sumtot;
		sumtot=sum(obs_p_srv_length_fem(i,j)+obs_p_srv_length_mal(i,j));
        obs_p_srv_length_mal(i,j)=obs_p_srv_length_mal(i,j)/sumtot;  //changing these to proportions rather than numbers
        obs_p_srv_length_fem(i,j)=obs_p_srv_length_fem(i,j)/sumtot;
        offset(i+1)-= nsamples_srv_length_fem(i,j)*obs_p_srv_length_fem(i,j)*log(obs_p_srv_length_fem(i,j)+offset_const)
                   +nsamples_srv_length_mal(i,j)*obs_p_srv_length_mal(i,j)*log(obs_p_srv_length_mal(i,j)+offset_const); 
	}
  } 
  //survey age offsets
  for (i=1;i<=nsurv_aged;i++)
  {
    for (j=1;j<=nobs_srv_age(i);j++)
    {
	double sumtot;
	sumtot=sum(obs_p_srv_age_fem(i,j)+obs_p_srv_age_mal(i,j));
	obs_p_srv_age_fem(i,j)=(obs_p_srv_age_fem(i,j)/sumtot)*agerr_matrix;
	obs_p_srv_age_mal(i,j)=(obs_p_srv_age_mal(i,j)/sumtot)*agerr_matrix;
	offset(i+nsurv+1)-=nsamples_srv_age(i,1,j)*obs_p_srv_age_fem(i,j)*log(obs_p_srv_age_fem(i,j)+offset_const)+
	             nsamples_srv_age(i,2,j)*obs_p_srv_age_mal(i,j)*log(obs_p_srv_age_mal(i,j)+offset_const);
    }  
  } 
                                                 
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& evalout= *pad_evalout;
   get_selectivity();     
   get_mortality();  
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));   
   get_numbers_at_age();
   get_catch_at_age();
  if (active(F40))
    compute_spr_rates();
  if (last_phase())
  {
    Future_projections();      
  }
  if (sd_phase() || mceval_phase()) 
   { 
	Do_depend();
    if (mceval_phase())
    {
      evalout << obj_fun << " " ;
      for (i=styr;i<=endyr;i++)
      evalout<<  fspbio(i) << " ";
      for (i=styr;i<=endyr;i++)    
      evalout<< natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" ";
      for (i=styr;i<=endyr;i++)
      evalout << 2*natage(1,i,1) <<" ";
      evalout <<  endl;
    }
  }
    evaluate_the_objective_function();  
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::get_selectivity(void)
{
  ofstream& evalout= *pad_evalout;
  if(active(log_selcoffs_fish))// 
 {           
    for(k=1;k<=2;k++)
    {
     for (j=1;j<=nselages;j++)
      {
        log_sel_fish(k,j)=log_selcoffs_fish(k,j);
      }
    }
   //sets selectivity of ages older than nselages to selectivity at nselages  
   for(k=1;k<=2;k++)
   {
     for (j=nselages+1;j<=nages;j++)
     {
       log_sel_fish(k,j)=log_sel_fish(k,j-1);
     }
   }
 for(k=1;k<=2;k++)
   {
     avgsel_fish(k)=log(mean(mfexp(log_selcoffs_fish(k))));
   }
 //vector=vector-scalar same as  vector-=scalar  
 //scaling selectivities by subracting the mean so exp(mean(s))=1.   
 //selectivities can be greater than 1 but mean is 1. 
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
    }     
 }//  end if(active(log_selcoffs_fish))
  else
    { 
	   sel(1)=get_sel(fishsel_params_f(1),fishsel_params_f(2));  
       sel(2)=get_sel(fishsel_params_m(1),fishsel_params_m(2));  
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
            if(j<=nselages)
             {
               sel(1,j)=1./(1.+mfexp(-1.*fishsel_params_f(1)*(double(j)-fishsel_params_f(2))));
               sel(2,j)=1./(1.+mfexp(-1.*fishsel_params_m(1)*(double(j)-fishsel_params_m(2))));
             }
            else
            {
             sel(1,j)=sel(1,j-1);
             sel(2,j)=sel(2,j-1);
            }    				
          } 
     }
  for(i=1;i<=nsurv;i++)   
  {
     if(nsel_params(i)==4)
     {              
     sel_srv(1,i) = get_sel(srv_params_f(i),srv_params_f(i+1),srv1desc_params_f(1),srv1desc_params_f(2));
     sel_srv(2,i) = get_sel(srv_params_m(i),srv_params_m(i+1),srv1desc_params_m(1),srv1desc_params_m(2));         
     }
     if(nsel_params(i)==2)
       {
       sel_srv(1,i) = get_sel(srv_params_f((2*i)-1),srv_params_f(2*i));
       sel_srv(2,i) = get_sel(srv_params_m((2*i)-1),srv_params_m(2*i)); 
       }        
  }  
  for (j=1;j<=nages;j++)   
  {           
	for (i=1;i<=nsurv;i++) 
	{
    if (j>nselages_srv(i))
    {
	sel_srv(1,i,j)=sel_srv(1,i,j-1); 
	sel_srv(2,i,j)=sel_srv(2,i,j-1); 
    } 
    } 
  }
}

dvar_vector model_parameters::get_sel(const dvariable& slp, const dvariable& a50)
{
  ofstream& evalout= *pad_evalout;
   {
	dvar_vector sel_tmp(1,nages);
    for (j=1;j<=nages;j++)  //this is selectivity for the surveys
    {  
    sel_tmp(j)=1./(1.+mfexp(-1.*slp*(double(j)-a50))); 
    }          
    return(sel_tmp);
   }          
}

dvar_vector model_parameters::get_sel(const dvariable& slp, const dvariable& a50, const dvariable& dslp, const dvariable& d50)
{
  ofstream& evalout= *pad_evalout;
   {
	dvar_vector sel_tmp(1,nages);
   for (j=1;j<=nages;j++)  //this is selectivity for the surveys         
   {
	  sel_tmp(j) = 1./(1.+mfexp(-1.*slp*(double(j)-a50)));           
      sel_tmp(j) *= 1./(1.+mfexp(dslp*(double(j)-d50)));
   }
 	return(sel_tmp);
  }          
}

void model_parameters::get_mortality(void)
{
  ofstream& evalout= *pad_evalout;
  maxsel_fish=max(sel(1));     //1 is females
  if(maxsel_fish<max(sel(2)))  //if highest female selectivity is > male selectivity, make maxsel_fish=male high selectivity
      maxsel_fish=max(sel(2));
  	  if (assess==3)   //added to match Kamchatka assessment - consider revising this
	  {
	  maxsel_fish=1.0;
      }
  fmort = mfexp(log_avg_fmort+fmort_dev); 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i<=endyr;i++)
    {
      F(k,i)=(sel(k)/maxsel_fish)*fmort(i);
      Z(k,i)=F(k,i) + M(k); 
    }
  } 
  S = mfexp(-1.0*Z); 
}

void model_parameters::get_numbers_at_age(void)
{
  ofstream& evalout= *pad_evalout;
  for(i=1;i<=nsurv;i++)
  {
   maxsel_srv(i)=max(sel_srv(1,i));
   if(maxsel_srv(i)<max(sel_srv(2,i)))
   maxsel_srv(i)=max(sel_srv(2,i));
  }     
    if (assess==3)   //added to match Kamchatka assessment - consider revising this
  {    
  for(i=2;i<=nsurv;i++)  
  maxsel_srv(i)=1.0;  
  }   
  if (assess==3)   //added to match Kamchatka assessment - consider removing this
  {
  maxsel_srv(1)=max(sel_srv(1,1));
  sel_srv(1,1) /= maxsel_srv(1); // shelf survey normalized by 4th age class 
  maxsel_srv(1)=max(sel_srv(2,1));
  sel_srv(2,1) /= maxsel_srv(1); // shelf survey normalized by 4th age class 
  }
  int itmp;
 //calc initial population  
  for (j=1;j<nages;j++)
    {
      itmp=styr+1-j;
      natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
      natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));
    }
    itmp=styr+1-nages;
  //last age    
    natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
    natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));
 // Now do for next several years----------------------------------
  if(assess<3)
  {
  for (i=styr+1;i<=endyr;i++)
  {
    //for age 1 recruits in the last year use value read in from data file
    if(i<=(endyr-1))
    {
      natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
      natage(2,i,1)=natage(1,i,1);
    }
    else
    {
      natage(1,i,1)=median_rec;
      natage(2,i,1)=natage(1,i,1);
    }
  }
  }
  if(assess==3)
  {
   for (i=styr+1;i<=endyr;i++)
   {
  //for age 1 recruits in the last year use value read in from data file
   natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
   natage(2,i,1)=natage(1,i,1);
   }  
  }
 //numbers at age 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i< endyr;i++)
    {
      for(j=1;j<nages;j++)
      {
        natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
      }
      natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
      popn(k,i)= natage(k,i)*sel(k);
    }
    popn(k,endyr)=natage(k,endyr)*sel(k);
  }
  for (i=styr;i<=endyr;i++)
  {
      pred_sexr(i)=sum(natage(2,i))/(sum((natage(1,i)+natage(2,i))));  //calculation of prop. of males in pred. population 
  }
  //predicted survey values
  fspbio.initialize(); 
  qtime=q_surv(1); 
  for (j=1;j<=nsurv;j++)
 {
    for (i=styr;i<=endyr;i++)
   {
  fspbio(i) = natage(1,i)*elem_prod(wt(1),maturity);
  explbiom(i)=0.;
  pred_bio(i)=0.; 
  pred_srv(j,i)=0.;
  //catchability calculation for survey years
  if (assess==1 && (i>=1982) && (i-1981 <= nobs_srv(1)))      //JNI catchability calculation for survey years    
   {   
   qtime(i)=q_surv(1)*mfexp(-alpha+beta*bottom_temps(i-1981));
   }
  for(k=1;k<=2;k++)
    {
    if (j==1 && assess==1)
      {             
    pred_srv(j,i) += qtime(i)*(natage(k,i)*elem_prod(sel_srv(k,j)/maxsel_srv(j),wt(k)));maxsel_srv(j);   //shelf survey, dividing by the maxsel constrains female selectivity to be 1.0
      } 
    else 
      {
    pred_srv(j,i) += q_surv(j)*(natage(k,i)*elem_prod(sel_srv(k,j),wt(k)));///maxsel_srv(j);         //slope survey JNI  do not need to divide by maxsel_srv if it is logistic but does not hurt
      } 
       //Aleutian Islands survey JNI
    //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
    // are set different in the tpl file the program will take to value from the bin file and use that 
    explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
    pred_bio(i)+=natage(k,i)*wt(k);
      }
   }  
 }     
  //Fit survey length compositions
  for (i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{
			ii=yrs_srv_length(i,j);
			pred_p_srv_len_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii))*lenage(1);
			pred_p_srv_len_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii))*lenage(2);
			dvariable sum_tot=sum(pred_p_srv_len_fem(i,j)+pred_p_srv_len_mal(i,j));
			pred_p_srv_len_fem(i,j)/=sum_tot;
			pred_p_srv_len_mal(i,j)/=sum_tot;
	 }
   }
  //Fit survey age composition
   for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
			ii=yrs_srv_age(i,j);
			pred_p_srv_age_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii));
			pred_p_srv_age_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii));
			dvariable sum_tot=sum(pred_p_srv_age_fem(i,j)+pred_p_srv_age_mal(i,j));
			pred_p_srv_age_fem(i,j)/=sum_tot;
			pred_p_srv_age_mal(i,j)/=sum_tot;
	 }
   }
  depletion=pred_bio(endyr)/pred_bio(styr);
  endbiom=pred_bio(endyr);
}

void model_parameters::get_catch_at_age(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr; i<=endyr; i++)
  {
    pred_catch(i)=0.;
    for(k=1;k<=2;k++)
    {      
      //--Baranov's equation here-----------------------------------
      for (j = 1 ; j<= nages; j++)
      {
        catage(k,i,j) = natage(k,i,j)*F(k,i,j)*(1.-S(k,i,j))/Z(k,i,j);
        pred_catch(i) += catage(k,i,j)*wt(k,j);
      }
      pred_p_fish(k,i)=elem_prod(sel(k),natage(k,i))*lenage(k)/(popn(1,i)+popn(2,i));
    }
  }   
}

void model_parameters::Future_projections(void)
{
  ofstream& evalout= *pad_evalout;
  for(k=1;k<=2;k++)
  {
    nage_future(k,styr_fut)(2,nages)=++elem_prod(natage(k,endyr)(1,nages-1),S(k,endyr)(1,nages-1));
    nage_future(k,styr_fut,nages)+=natage(k,endyr,nages)*S(k,endyr,nages);
   }
    future_biomass.initialize();
    catch_future.initialize();
    for (int l=1;l<=4;l++)
    {
      switch (l)
      {
        case 1:
          ftmp=F40;
          break;
        case 2:
          ftmp=F35;
          break;
        case 3:
          ftmp=F30;
          break;
        case 4:
          ftmp.initialize();
          break;
      }
     for(k=1;k<=2;k++)
     {
      for (i=endyr+1;i<=endyr_fut;i++)
      {
        for (j=1;j<=nages;j++)
        {
          F_future(k,i,j) = (sel(k,j)/maxsel_fish)*ftmp;
          Z_future(k,i,j) = F_future(k,i,j)+M(k);
          S_future(k,i,j) = exp(-1.*Z_future(k,i,j));
        }
      }
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(k,i,1)  = median_rec;
       // Now graduate for the next year....
        nage_future(k,i+1)(2,nages) = ++elem_prod(nage_future(k,i)(1,nages-1),S_future(k,i)(1,nages-1));
        nage_future(k,i+1,nages)   += nage_future(k,i,nages)*S_future(k,i,nages);
      }
      nage_future(k,endyr_fut,1)  = median_rec;
      // Now get catch at future ages
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        for (j = 1 ; j<= nages; j++)
        {
          catage_future(k,i,j) = nage_future(k,i,j) * F_future(k,i,j) * ( 1.- S_future(k,i,j) ) / Z_future(k,i,j);
         if(k==1)
          {
          fspbiom_fut(l,i) += nage_future(1,i,j)*wt(1,j)*maturity(j);
          }
        }
        if (l!=4) catch_future(l,i)   += catage_future(k,i)*wt(k);
        future_biomass(l,i) += nage_future(k,i)*wt(k);
      }   //end loop over future years
     }   //end loop over sex
     fspbiom_fut(l)=0.;
     for(i=styr_fut;i<=endyr_fut;i++)
       fspbiom_fut(l,i) = elem_prod(nage_future(1,i),wt(1)) * maturity;
    }   //End of loop over F's
}

void model_parameters::compute_spr_rates(void)
{
  ofstream& evalout= *pad_evalout;
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0.initialize();
  SBF40.initialize();
  SBF35.initialize();
  SBF30.initialize();
  // Initialize the recruit (1) for each F  (F40 etc)
  for (i=1;i<=3;i++)
    Nspr(i,1)=1.;
  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*exp(-1.*M(1));
    Nspr(2,j)=Nspr(2,j-1)*exp(-1.*(M(1)+F40*sel(1,j-1)/maxsel_fish));
    Nspr(3,j)=Nspr(3,j-1)*exp(-1.*(M(1)+F35*sel(1,j-1)/maxsel_fish));
    Nspr(4,j)=Nspr(4,j-1)*exp(-1.*(M(1)+F30*sel(1,j-1)/maxsel_fish));
  }
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*exp(-1.*M(1))/(1.-exp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*exp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*exp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  Nspr(4,nages)=Nspr(4,nages-1)*exp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));
  for (j=1;j<=nages;j++)
  {
   // Kill them off till april (0.25) atf spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
  	//##################################Below the spawning fraction is made so that there is no spawn_fract 
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*exp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
    SBF30  += Nspr(4,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F30*sel(1,j)/maxsel_fish));
  }
  sprpen    = 200.*square((SBF40/SB0)-0.4);
  sprpen   += 200.*square((SBF35/SB0)-0.35);
  sprpen   += 200.*square((SBF30/SB0)-0.30);
}

void model_parameters::Do_depend(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr;  i<=endyr;  i++) 
  totalbiomass(i)=natage(1,i)*wt(1) + natage(2,i)*wt(2);
  obj_fun += 1.*sexr_like;             // male proportion prior, emphasis factor = 1
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& evalout= *pad_evalout;
  length_like.initialize();
  age_like.initialize();   
  sel_like.initialize(); 
  fpen.initialize();
  rec_like.initialize();
  surv_like.initialize();
  catch_like.initialize();
  sexr_like.initialize();
  obj_fun.initialize(); 
  if (active(rec_dev))
  {
  length_like.initialize();
  int ii;
    //recruitment likelihood - norm2 is sum of square values   
    rec_like = norm2(rec_dev);
    //length likelihood
    for(k=1;k<=2;k++)
    {
      for (i=1; i <= nobs_fish; i++)
      {
        ii=yrs_fish(i);
        //fishery length likelihood fitting
          length_like(1) -= nsamples_fish(k,i)*(offset_const+obs_p_fish_norm(k,i))*log(pred_p_fish(k,ii)+offset_const);
      }
    }
    //add the offset to the likelihood   
    length_like(1)-=offset(1);
  //survey length composition fitting
   for (i=1;i<=nsurv;i++)
   {
     for (j=1;j<=nobs_srv_length(i);j++) 
     {    
	   length_like(i+1)-=((nsamples_srv_length_fem(i,j)*(offset_const+obs_p_srv_length_fem(i,j))*log(pred_p_srv_len_fem(i,j)+offset_const))
	                      +(nsamples_srv_length_mal(i,j)*(offset_const+obs_p_srv_length_mal(i,j))*log(pred_p_srv_len_mal(i,j)+offset_const)));
	 } 
	  length_like(i+1)-=offset(i+1); 
    } 
  for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
   age_like(i)-=nsamples_srv_age(i,1,j)*(offset_const+obs_p_srv_age_fem(i,j))*log(pred_p_srv_age_fem(i,j)+offset_const)+
                  nsamples_srv_age(i,2,j)*(offset_const+obs_p_srv_age_mal(i,j))*log(pred_p_srv_age_mal(i,j)+offset_const); 	
	}	
	age_like(i)-=offset(i+nsurv+1);
  }
  //end of if(active (rec_dev))
  }
  // Fit to indices (lognormal)  
  //weight each years estimate by 1/(2*variance) - use cv as an approx to s.d. of log(biomass) 
  for (i=1;i<=nsurv;i++)
  {   
  surv_like(i) = norm2(elem_div(log(obs_srv(i))-log(pred_srv(i)(yrs_srv(i))),sqrt(2)*cv_srv(i)));   //yrs_srv(i) is an index here.
  }  
   catch_like=norm2(log(catch_bio+offset_const)-log(pred_catch+offset_const));
   // sex ratio likelihood
   sexr_like=0.5*norm2((obs_mean_sexr-pred_sexr)/obs_SD_sexr); 
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences  
  if(active(log_selcoffs_fish))
  {  
    sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1)))); //fishery females
    sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2)))); //fishery males 
   for (j=1;j<nages;j++)
   {
    if(monot_sel==1)
    { 
        if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
        sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));   //monotonicity constraint fishery females
        if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
        sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));  //monotonicity constraing fishery males
    }
   }  
    obj_fun+=1.*sum(sel_like);    
    obj_fun+=1.*square(avgsel_fish(1));
    obj_fun+=1.*square(avgsel_fish(2));    
  } //end if active(log_selcoffs_fish)
  // Phases less than 3, penalize High F's
    if (current_phase()<2)
    {
       //F's are low for arrowtooth changed the value to compare from .2 to .001
       //don't know if makes any difference since the penalty is reduced at the end
       fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-offset_const);
    }
    else
    {
      fpen=fpen_mult(1)*norm2(mfexp(fmort_dev+log_avg_fmort)-offset_const);      //0.001 for BSAI and 0.01 for GOA see fpen_mult values  
    }
    if (active(fmort_dev))
    {
      fpen+=fpen_mult(2)*norm2(fmort_dev);    //0.01 for BSAI and 1 for GOA
    }
  obj_fun += rec_like;
  obj_fun += like_wght(1)*length_like(1);    //this is fishery length likelihood
  for(i=2;i<=nsurv+1;i++)
   {
   obj_fun+=like_wght(6)*length_like(i);     //survey likelihood
   } 
  obj_fun+= like_wght(7)*sum(age_like); 
  for(i=1;i<=nsurv;i++)
  {
  obj_fun += like_wght(i+1)*surv_like(i);  
  }     
  obj_fun += catch_err_like*like_wght(5)*catch_like;        // large emphasis to fit observed catch 
  //obj_fun += fpen;   
  obj_fun += sprpen;     
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Styr" << endl<<styr<<endl;  
  report << "Endyr" << endl<<endyr<<endl; 
  report << "nages_read" <<endl<<nages_read<<endl;  
  report << "nobs_srv" <<endl<<nobs_srv<<endl; 
  report << "nobs_fish" <<endl<<nobs_fish<<endl; 
  report << "yrs_fish" <<endl<<yrs_fish<<endl; 
  report << "yrs_srv" <<endl<<yrs_srv<<endl; 
  report << "nobs_srv_length" <<endl<<nobs_srv_length<<endl; 
  report << "nobs_srv_age"<<endl<<nobs_srv_age<<endl;
  report << "assessment 1 is BSAI ATF, 2 is GOA ATF, 3 is Kamchatka"<<endl<<assess<<endl; 
  report << "Numbers_fem" << endl<<natage(1)<<endl;
  report << "Numbers_mal" << endl<<natage(2)<<endl;
  report << "Catch_est_fem"<<endl<<catage(1)<<endl;   
  report << "Catch_est_mal"<<endl<<catage(2)<<endl; 
  report << "Mort_est_fem"<<endl<<F(1)<<endl; 
  report << "Mort_est_mal"<<endl<<F(2)<<endl; 
  report << "Fishsel_fem"<<endl<<sel(1)/maxsel_fish<<endl; 
  report << "Fishsel_mal"<<endl<<sel(2)/maxsel_fish<<endl; 
  report << "Survsel_fem"<<endl<<sel_srv(1)/maxsel_srv(1)<<endl; 
  report << "Survsel_mal"<<endl<<sel_srv(2)/maxsel_srv(1)<<endl; 
  report << "Obs_srv_biomass"<<endl<<obs_srv<<endl; 
  report << "Pred_srv_biomass"<<endl<<pred_srv<<endl;
  report << "obs_srv_sd" <<endl <<obs_srv_sd<<endl;  
  report << "Obs_srv_lengthcomp_fem"<<endl<<obs_p_srv_length_fem<<endl; 
  report << "Pred_srv_lengthcomp_fem"<<endl<<pred_p_srv_len_fem<<endl; 
  report << "Obs_srv_lengthcomp_mal"<<endl<<obs_p_srv_length_mal<<endl; 
  report << "Pred_srv_lengthcomp_mal"<<endl<<pred_p_srv_len_mal<<endl;
  report << "Yrs_srv_lengthcomp"<<endl<<yrs_srv_length<<endl; 
  report << "Yrs_srv_age"<<endl<<yrs_srv_age<<endl;
  report << "Obs_srv_agecomp_fem"<<endl<<obs_p_srv_age_fem<<endl; 
  report << "Pred_srv_agecomp_fem"<<endl<<pred_p_srv_age_fem<<endl; 
  report << "Obs_srv_agecomp_mal"<<endl<<obs_p_srv_age_mal<<endl;
  report << "Pred_srv_agecomp_mal"<<endl<<pred_p_srv_age_mal<<endl;
  report << "obs_p_fish" << endl<<obs_p_fish << endl;
  report << "obs_p_fish_norm" << endl<<obs_p_fish_norm << endl;
  report << "pred_p_fish" << endl<<pred_p_fish << endl;
  report << "Obs_catch"<<endl<<catch_bio<<endl; 
  report << "Pred_catch"<<endl<<pred_catch<<endl; 
  report << "FSB"<<endl<<fspbio<<endl; 
  report << "fmort"<<endl<<fmort<<endl;
  report << "F40"<<endl<<F40<<endl;
  report << "F35"<<endl<<F35<<endl;
  report << "F30"<<endl<<F30<<endl;
  report << "SBF40"<<endl<<SBF40<<endl;
  report << "SBF35"<<endl<<SBF35<<endl;
  report << "SBF30"<<endl<<SBF30<<endl;
  report << "SB0"<<endl<<SB0<<endl; 
  report << "surv_like"<<endl<<surv_like<<endl;
  report << "length_like_fishery"<<endl<<length_like(1)<<endl; 
  report << "length_like_survey"<<endl<<length_like(2)<<endl; 
  report << "age_like_survey1"<<endl<<age_like<<endl;  
  report << "fpen" <<endl<<fpen<<endl;   
  report << "catch_like"<<endl<<catch_like<<endl;   
  report << "rec_like"<<endl<<rec_like<<endl; 
  report << "total_likelihood"<<endl<<obj_fun<<endl; 
  report << "sexr_like"<<endl<<sexr_like<<endl;  
  report << "sel_like" <<endl<<sel_like<<endl;
  report <<"offset"<<endl<<offset<<endl; 
  report <<"explbiom"<<endl<<explbiom<<endl; 
  report <<"fspbio"<<endl<<fspbio<<endl; 
  report <<"pred_bio"<<endl<<pred_bio<<endl; 
  report <<"lenage"<<endl<<lenage<<endl;  
  report <<"fishmort"<<endl<<mfexp(log_avg_fmort+fmort_dev)<<endl;
  report <<"wt"<<endl<<wt<<endl; 
  report << "sprpen"<<endl<<sprpen<<endl;
  report << "nsamples_srv_length_fem"<<endl<<nsamples_srv_length_fem<<endl; 
  report << "nsamples_srv_length_mal"<<endl<<nsamples_srv_length_mal<<endl; 
  report << "fspbiom_fut"<<endl<<fspbiom_fut<<endl;
  report << "q_surv"<<endl<<q_surv<<endl;
  report << "M"<<endl<<M<<endl;
  report << "q_time_catchability"<<endl<<qtime<<endl;
  report << "pred_sexr"<<endl<<pred_sexr<<endl; 
  report << "obs_mean_sexr"<<endl<<obs_mean_sexr<<endl;
  report << "obs_SD_sexr"<<endl<<obs_SD_sexr<<endl; 
  report << "alpha"<<endl<<alpha<<endl; 
  report << "beta"<<endl<<beta<<endl; 
  report << "obs_srv_sd"<<endl<<obs_srv_sd<<endl; 
  report << "mean_log_rec"<<endl<<mean_log_rec<<endl; 
  report << "rec_dev"<<endl<<rec_dev<<endl; 
  report<<"Num_parameters_Estimated "<<endl<<initial_params::nvarcalc()<<endl;
    report << "obj_fun"<<endl<<obj_fun<<endl;
  report << "like_wght"<<endl<<like_wght<<endl;
  report << "nsamples_fish"<<endl<<nsamples_fish<<endl;
  report<<"nsamples_srv_age"<<endl<<nsamples_srv_age<<endl; 
  report << "SARA form for Angie Grieg" << endl;
  report << "ATF        # stock  " << endl;
  report << "BSAI       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
  report << "2013       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
  report << "3a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
  report << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
  report << "partial    # UPDATE (new benchmark full partial)" << endl;
  report << "2          # LIFE_HIST - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "2          # ASSES_FREQ - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "5          # ASSES_LEV - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "5          # CATCH_DAT - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "3          # ABUND_DAT - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "567000     # Minimum B  Lower 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "665000     # Maximum B  Upper 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "202138     # BMSY  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
  report << "ADMB       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
  report << "           # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
  report << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
  report << "1          # number of fisheries" << endl;
  report << "1          # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
  report << "1          # recruitment age used by model or size" << endl;
  report << "1          # age+ or mmCW+ used for biomass estimate" << endl;
  report << "Single age        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
  report << "Age model         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
  report << "Age of maximum F  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
  report << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
  report << "ALL" << endl; 
  report <<"#FISHERYYEAR - list years used in the model " << endl;
    for (i=styr;  i<=endyr; i++)
       report << i << "	";
       report<<endl;  
  report<<"#AGE - list of ages used in the model"<<endl;
    for (i=1; i<=21;i++)
       report << i << "	";
       report<<endl;    
  report <<"#RECRUITMENT - Number of recruits by year " << endl;
    for (i=styr;  i<=endyr;  i++)
	   report  << 2*natage(1,i,1) << "	";
	   report<<endl;     
  report <<"#SPAWNBIOMASS - Spawning biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << fspbio(i) << "	";
       report<<endl;  
  report <<"#TOTALBIOMASS - Total biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << natage(1,i)*wt(1) + natage(2,i)*wt(2) << "	";
       report<<endl;
  report <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
	for (i=styr;  i<=endyr;  i++)
	   report  << (F(1,i,13)+ F(2,i,13))/2<< "	";
	   report<<endl;
  report <<"#TOTALCATCH - Total catch by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << catch_bio(i) << "	";
       report<<endl;
 report <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
       for (i=1;  i<=13;  i++) 
       report  << maturity(i) <<"	";
       report<< endl; 
 report <<"#SPAWNWT - Average spawning weight (in kg) by age"<< endl; 
       report <<"0.019	0.041	0.113	0.224	0.376	0.566	0.784	1.028	1.292	1.569	1.855	2.142	2.417	2.667	2.881	3.057	3.198	3.308	3.393"<<endl;                              
       report<<endl;
 report <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
 for (i=1;  i<=13;  i++) 
 report  << 0.2 <<"	";
 report<< endl;   
 for (i=1;  i<=13;  i++) 
 report  << 0.35 <<"	";
 report<< endl;
 report << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
   for (i=styr; i<=endyr;i++)
     report <<natage(1,i)<< "	";
     report<<endl;
   for (i=styr; i<=endyr;i++)
     report <<natage(2,i)<< "	";
     report<<endl;
 report <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
    report <<wt(1)*1000  << "	";
    report<<endl; //1 is females        
    report <<wt(2)*1000  << "	";
    report<<endl; //2 is males
  report << "#SELECTIVITY - Estimated fishery selectivity for females (first) males (second) at age " << endl;
    for (j=1; j<=nages;j++)
      report <<" " <<sel(1,j)<< "	";
      report<<endl;
    for (j=1; j<=nages;j++)
      report <<" "  <<sel(2,j)<< "	";
      report<<endl;
 report << "#SURVEYDESC"<<endl;
 report<<"EBS_trawl_survey BS_slope_trawl_survey AI_trawl_survey"<<endl;
 report<<"SURVEYMULT"<<endl;
 report<<"1 1 1"<<endl;
 report << "#GOA_trawl_survey - Gulf of Alaska survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nsurv;i++)
      report << yrs_srv(i) << "	";
      report<<endl;
    for (i=1; i<=nsurv;i++) 
      report<< obs_srv(i)<< "	";
      report<< endl;
    for (i=1; i<=nsurv;i++){     
    for (j=1; j<=nobs_srv(i);j++){
      report<< pred_srv(i,yrs_srv(i,j))<< "	"; }}
      report<< endl;
 report<<"#STOCKNOTES"<<endl;
 report<<"SAFE report indicates that this stock was not subjected to overfishing in 2012 and is neither overfished nor approaching a condition of being overfished in 2013."<<endl;
 write_newproj();
}

void model_parameters::write_newproj(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream newproj("proj.dat");
 newproj <<"#Species name here:"<<endl;
 newproj <<"GOA_ATF"<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"2"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 newproj << mean(fmort(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << "2"<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<M_f<<endl;
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<M_m<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_sdivide_by_2"<<aa<<endl<<maturity<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt(1)*1000<< endl;
 newproj <<"#_Wt_at_age_fishery1 female" <<aa<<endl<<wt(1)*1000<< endl;
 newproj <<"#_Wt_at_age_spawners1 male"<<aa<<endl<<wt(2)*1000<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<sel(1)/max(sel(1))<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<sel(2)/max(sel(2))<< endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage(1,endyr)/1000 <<endl<<natage(2,endyr)/1000 << endl;
 newproj <<"#_N_recruitment_years"<<endl<<"38"<<endl; //endyr-1978-1<< endl;
 // newproj <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<endl<<mfexp(mean_log_rec+rec_dev(1977,2015))/1000*2<< endl;
 newproj <<"#_Spawners per recruitment (starting at 1977)"<<endl<<fspbio(1977,endyr-3)<< endl;
 newproj.close();
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{4000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3 1e-4 1e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(300);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
