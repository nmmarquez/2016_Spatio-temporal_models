#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_VECTOR(log_rate_mort);
    DATA_VECTOR(age);
    DATA_IVECTOR(age_group);
    DATA_INTEGER(option);
    
    // Parameters
    PARAMETER_VECTOR(log_gpz);
    PARAMETER(log_sigma_obs);
    PARAMETER_VECTOR(age_fixed);
    
    // Transforms
    vector <Type> gpz = exp(log_gpz);
    Type sigma_obs = exp(log_sigma_obs);
    int N = log_rate_mort.size();
    
    // Objective funcction
    Type nll = 0.;
    
    // Probability of random effects
    // Comprabale to an AR1
    vector<Type> log_rate_mort_hat(N);
    Type alpha;
    Type beta;
    Type gamma;
    Type delta;
    Type zeta;
    
    vector<Type> age_effect(N);
    for(int n=0; n<N; n++){
        age_effect[n] = age_fixed[age_group[n]];
    }
    
    
    using namespace density;
    if(option == 0){
        alpha = gpz[0];
        beta = gpz[1];
        gamma = gpz[2];
        log_rate_mort_hat = log(alpha * exp(beta * age) + gamma);
    }
    
    if(option == 1){
        alpha = gpz[0];
        beta = gpz[1];
        gamma = gpz[2];
        delta = gpz[3];
        zeta = gpz[4];
        log_rate_mort_hat = log(alpha * exp(-1. * beta * age) + gamma + 
            delta * exp(zeta * age)) + age_effect;
    }
    
    for(int n=0; n<N; n++){
        nll -= dnorm(log_rate_mort[n], log_rate_mort_hat[n], sigma_obs, true);
    }
        
    // Reporting
    REPORT(alpha);
    REPORT(beta);
    REPORT(gamma);
    REPORT(delta);
    REPORT(zeta);
    REPORT(sigma_obs);
    REPORT(nll);
    REPORT(age_fixed);
    
    return nll;
    }
