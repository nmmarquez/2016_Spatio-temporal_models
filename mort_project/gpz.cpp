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
    PARAMETER_VECTOR(gpz);
    PARAMETER(log_sigma_obs);
    PARAMETER_VECTOR(epsilon_age);
    PARAMETER_VECTOR(epsilon_time);
    PARAMETER_ARRAY(epsilon_age_time);
    PARAMETER(logit_rho_age);
    PARAMETER(logit_rho_time);
    PARAMETER(logit_rho_age2);
    PARAMETER(logit_rho_time2);
    printf("%s\n", "Parameters loaded");
    
    // Transforms
    // vector <Type> gpz = exp(log_gpz);
    Type sigma_obs = exp(log_sigma_obs);
    Type rho_age = 1 / (1 + exp(-logit_rho_age));
    Type rho_time = 1 / (1 + exp(-logit_rho_time));
    Type rho_age2 = 1 / (1 + exp(-logit_rho_age2));
    Type rho_time2 = 1 / (1 + exp(-logit_rho_time2));
    int N = log_rate_mort.size();
    printf("%s\n", "Data Transformed");
    
    // Objective funcction
    Type nll = 0.;
    
    // Probability of random effects
    vector<Type> age_effect(N);
    vector<Type> time_effect(N);
    vector<Type> age_time_effect(N);

    using namespace density;
    if (option >= 2){
        nll += AR1(rho_age)(epsilon_age);
    }
    
    if (option >= 3){
        nll += AR1(rho_time)(epsilon_time);
    }
    if (option >= 4){
        nll += SEPARABLE(AR1(rho_age2), AR1(rho_time2))(epsilon_age_time);
    }

    printf("%s\n", "Setting age effects per observation");
    for(int n=0; n<N; n++){
        age_effect[n] = epsilon_age[age_group[n]];
        time_effect[n] = epsilon_time[time_group[n]];
        age_time_effect[n] = epsilon_age_time(age_group[n], time_group[n]);
    }
    
    printf("%s\n", "make_predictions");
    Type alpha = exp(gpz[0]);
    Type beta = gpz[1];
    Type gamma = exp(gpz[2]);
    Type delta = exp(gpz[3]);
    Type zeta = gpz[4];
    vector<Type> log_rate_mort_hat = log(alpha * exp(-1. * beta * age) + gamma +
        delta * exp(zeta * age)) + age_effect + time_effect + age_time_effect;

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
