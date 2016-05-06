rm(list=ls())
set.seed(123)
library(TMB)
# library(knitr)
library(RandomFields)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework/")

model_name <- "hw5"
if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

M <- 100

Sim_Fn <- function( n_i=1000, Scale=2, Sigma2=1, logSD_spatial=0.1, L0=10, 
                   Linf_0=100, beta_y=0.02, growth_rate=0.1, print=NULL, plot=F,
                   mortality_rate=growth_rate*1.6, logSD_resid=0.05){
    if(!is.null(print)){
        print(print)  
    }
    y_i = runif( n=n_i, min=32, max=49 )
    
    # Simulate spatial process
    RMmodel = RMgauss(var=logSD_spatial^2, scale=Scale)
    Linf_i = Linf_0 * exp( RFsimulate(model=RMmodel, x=rep(0,n_i), y=y_i)@data[,1] - 
                               Sigma2/2 ) * exp( beta_y*(y_i-40.5) )
    if(plot){
        plot( y=Linf_i, x=y_i )
    }
    
    # Simulate ages of samples
    Survival_a = exp( -mortality_rate * 1:100 )
    a_i = sample( size=n_i, x=1:100, prob=Survival_a, replace=TRUE)
    l_i = Linf_i - (Linf_i-L0) * exp( -growth_rate * a_i ) * 
        exp( rnorm(n_i, mean=-logSD_resid^2/2, sd=logSD_resid) )
    
    # Bundle and return stuff
    DF = data.frame( "y_i"=y_i, "a_i"=a_i, "l_i"=l_i, "Linf_i"=Linf_i)
    DF <- DF[order(DF$y_i),]
    DF
}

if (!file.exists("./data.rda")){
    df_list <- lapply(1:M, function(x) Sim_Fn(print=x, plot = T, print=x))
    save(df_list, file="./data.rda")
}
load("./data.rda")

run_model <- function(df, use_site=F, print=NULL, return_reports=F, re_option=0){
    if(!is.null(print)){
        print(print)  
    }
    N <- nrow(df)
    covs <- use_site + 1
    X <- df[, c("Linf_i", "y_i")]; X[,1] <- 1; X[,2] <- X[,2] - 40.5
    random <- c("pi")
    Map <- list()
    if (re_option > 2){
        random <- c()
        Map[["pi"]] <- factor(rep(NA, N))
        Map[["log_sigma_p"]] <- factor(NA)
    }
    dyn.load(dynlib(model_name))
    Params <- list(pi=rep(0, N), beta=rep(0, covs), log_l0=0, log_k=0,
                   logit_rho=0,  log_sigma_p=0, log_sigma_l=0)
    Data <- list(Y=df$l_i, A=df$a_i, P=df$y_i, X=as.matrix(X[,1:covs]),
                 re_option=re_option)
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, 
                     random=random, silent=T, map=Map)
    Obj$env$tracemgc <- FALSE
    Obj$env$inner.control$trace <- FALSE
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                  control=list(iter.max=1000, eval.max=500))
    Report <- Obj$report()
    dyn.unload(dynlib(model_name))
    if(Opt$convergence != 0){
        print("Model had partial or non_convergence!")
    }
    if(return_reports){
        return(Report)
    }
    else{
        return(Report$l_inf)
    }
}



results <- lapply(c(FALSE, TRUE), function(bool) sapply(1:M, function(i)
    run_model(df_list[[i]], use_site=bool, print=i)))

rmse <- function(true_vec, est_vec){
    (sum((true_vec - est_vec)**2) / length(true_vec))**.5
}

b0_rmse <- sapply(1:M, function(i) rmse(df_list[[i]]$Linf_i, results[[1]][,i]))
b1_rmse <- sapply(1:M, function(i) rmse(df_list[[i]]$Linf_i, results[[2]][,i]))
sum(sign(b0_rmse - b1_rmse))

Report <- run_model(df_list[[100]], use_site=T, return_reports = T)
plot(df_list[[100]]$y_i, Report$l_inf)
plot(Report$pi)
