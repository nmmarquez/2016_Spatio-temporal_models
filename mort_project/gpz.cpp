#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_VECTOR(log_rate_mort);
    DATA_VECTOR(age);
    DATA_IVECTOR(age_group);
    DATA_IVECTOR(time_group);
    DATA_INTEGER(option);
    printf("%s\n", "Data loaded");
    
    // Parameters
    PARAMETER_VECTOR(log_gpz);
    PARAMETER(log_sigma_obs);
    PARAMETER_VECTOR(epsilon_age);
    PARAMETER_VECTOR(epsilon_time);
    PARAMETER(logit_rho_age);
    PARAMETER(logit_rho_time);
    printf("%s\n", "Parameters loaded");
    
    // Transforms
    vector <Type> gpz = exp(log_gpz);
    Type sigma_obs = exp(log_sigma_obs);
    Type rho_age = 1 / (1 + exp(-logit_rho_age));
    Type rho_time = 1 / (1 + exp(-logit_rho_time));
    int N = log_rate_mort.size();
    printf("%s\n", "Data Transformed");
    
    // Objective funcction
    Type nll = 0.;
    
    // Probability of random effects
    vector<Type> age_effect(N);
    vector<Type> time_effect(N);

    using namespace density;
    if (option >= 2){
        nll += AR1(rho_age)(epsilon_age);
    }
    
    if (option >= 3){
        nll += AR1(rho_time)(epsilon_time);
    }

    printf("%s\n", "Setting age effects per observation");
    for(int n=0; n<N; n++){
        age_effect[n] = epsilon_age[age_group[n]];
        time_effect[n] = epsilon_time[time_group[n]];
    }
    
    printf("%s\n", "make_predictions");
    Type alpha = gpz[0];
    Type beta = gpz[1];
    Type gamma = gpz[2];
    Type delta = gpz[3];
    Type zeta = gpz[4];
    log_rate_mort_hat = log(alpha * exp(-1. * beta * age) + gamma +
        delta * exp(zeta * age)) + age_effect + time_effect;

    printf("%s\n", "evaluate data likelihood");
    for(int n=0; n<N; n++){
        nll -= dnorm(log_rate_mort[n], log_rate_mort_hat[n], sigma_obs, true);
    }

    printf("%s\n", "Report values");
    // Reporting
    REPORT(alpha);
    REPORT(beta);
    REPORT(gamma);
    REPORT(delta);
    REPORT(zeta);
    REPORT(sigma_obs);
    REPORT(nll);
    REPORT(epsilon_age);
    REPORT(epsilon_time);
    REPORT(age_effect);
    REPORT(age_group);
    REPORT(log_rate_mort_hat);
    return nll;
    }
