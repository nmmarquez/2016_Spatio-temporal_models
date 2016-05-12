
# working stuff
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework solution" )
library(TMB)

# Load simulator
source( "C:/Users/James.Thorson/Desktop/Project_git/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework/Homework_Week_5_simulator.R" )

# Compile
Version = "HW5"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Results to save
Results = array(NA, dim=c(100,2,4), dimnames=list(NULL,c("trend","no_trend"),c("RMSE","sigma_spatial","sigma_obs","beta_hat")) )

# loop through replicates
for( repI in 10:dim(Results)[1]){
  # Simulate data
  DF = Sim_Fn()
  # Sort from south to north
  DF = DF[ order(DF[,'y_i']), ]

  # Make inputs
  Data = list("a_i"=DF[,'a_i'], "l_i"=DF[,'l_i'], "loc_i"=DF[,'y_i']-mean(DF[,'y_i']))  # mean-center the predictor variable
  Params = list( "log_Linf_mean"=log(10), "log_kappa"=log(0.2), "log_Lzero"=log(1), "log_rho"=log(0.2), "beta"=0, "ln_sigma_spatial"=1, "ln_sigma_measurement"=1, "epsilon_i"=rep(0,nrow(DF)) )
  Random = c("epsilon_i") #, "beta", "log_Linf_mean", "log_Lzero")  # Using REML, so only treating variance parameters as fixed effects

  # Loop through estimation models
  for( estI in 1:2){
    # Settings
    Map = list()
    if( estI == 2 ){
      Map[["beta"]] = factor(NA)
    }

    # Build object
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map )
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
    Report = Obj$report()
    # SD = sdreport( Obj )

    # Record estimates
    Results[repI,estI,"RMSE"] = sqrt( mean((DF$Linf_i - Report$Linf_i)^2) )
    Results[repI,estI,c("sigma_spatial","sigma_obs","beta_hat")] = unlist(Report[c("sigma_spatial","sigma_measurement","beta")])
  }
}

##############
# Check RMSE
##############

# Should show that model with trend is has much lower SD for spatial process, but essentially no difference in the RMSE when predicting Linf_i
# The lesson is that data are often sufficient to estimate spatial variation in parameters even when the model is missing covariates
# However, including covariates will generally decrease the SD of unexplained spatial variation, thus you can partition variance into a component explained by covariates
apply( Results, MARGIN=2:3, FUN=mean, na.rm=TRUE )



