rm(list=ls())
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 2 -- mixed-effects/Homework/")
set.seed(123)
library(TMB)
library(ggplot2)

true_mu <- 2

sim_clams <- function(mu=true_mu, sigma_s=sqrt(1), sigma_y=sqrt(.5), site_num=10, site_obs=10){
    lambda_site <- rlnorm(site_num, mu, sigma_s)
    # this doesnt seem right...
    y_mean <- sapply(lambda_site, function(x) rlnorm(site_obs, log(x), sigma_y))
    y_obs <- apply(y_mean, 2, function(x) rpois(site_obs, x))
    y_obs
}

M <- 1000

simulated_clam_data <- lapply(1:M, function(x) sim_clams())

model_name <- "hw2"
if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

run_model <- function(data, site_var=FALSE, ind_var=FALSE, use_REML=TRUE){
    random <- c()
    if(use_REML & (site_var | ind_var)) random <- union(random, "x0")
    if(site_var) random <- union(random, "z")
    if(ind_var) random <- union(random, "e")
    dyn.load(dynlib("hw2"))
    Params <- list(x0=0, z=rep(0, ncol(data) * site_var), 
                   log_sdz=rep(0, site_var), log_sde=rep(0, ind_var),
                   e=matrix(0, nrow=nrow(data)*ind_var, ncol=ncol(data)*ind_var))
    Data <- list(y=data, site_var=as.numeric(site_var), ind_var=as.numeric(ind_var))
    Obj <- MakeADFun(data=Data, parameters=Params, DLL="hw2", random=random)
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
    Report <- Obj$report()
    SD <- sdreport(Obj)
    beta_int <- summary(SD)["x0",]
    beta_int["lower_bound"] <- beta_int["Estimate"] - 1.96 * beta_int["Std. Error"]
    beta_int["upper_bound"] <- beta_int["Estimate"] + 1.96 * beta_int["Std. Error"]
    dyn.unload(dynlib("hw2"))
    beta_int
}

model_spec <- list(glm=list(site_var=FALSE, ind_var=FALSE),
                   glmm_site=list(site_var=TRUE, ind_var=FALSE),
                   glmm_ind=list(site_var=FALSE, ind_var=TRUE),
                   glmm_both=list(site_var=TRUE, ind_var=TRUE))

modeled_sims <- lapply(model_spec, function(x) 
    as.data.frame(t(sapply(simulated_clam_data, function(y)
        run_model(data=y, site_var=x$site_var, ind_var=x$ind_var)))))

for(models in names(modeled_sims)){
    modeled_sims[[models]]$sim_num <- 1:M
    modeled_sims[[models]]$std_err <- modeled_sims[[models]][,"Std. Error"]
}

plot_results <- function(results_df){
    p <- ggplot(results_df, aes(Estimate, sim_num))
    p + geom_point() + xlim(-1, 4) + 
        geom_vline(xintercept=true_mu, colour="red") +
        geom_errorbarh(aes(xmax = upper_bound, xmin = lower_bound))
}

plot_results(modeled_sims$glm)
plot_results(modeled_sims$glmm_site)
plot_results(modeled_sims$glmm_ind)
plot_results(modeled_sims$glmm_both)

sum(modeled_sims$glmm_both$upper_bound >= true_mu & 
    modeled_sims$glmm_both$lower_bound <= true_mu) / M
