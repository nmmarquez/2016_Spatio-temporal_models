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
    DATA_VECTOR(unique_age);
    printf("%s\n", "Data loaded");
    
    // Parameters
    PARAMETER_VECTOR(log_gpz);
    PARAMETER(log_sigma_obs);
    PARAMETER_VECTOR(epsilon_age);
    PARAMETER(log_sigma_age);
    PARAMETER(logit_rho);
    printf("%s\n", "Parameters loaded");
    
    // Transforms
    vector <Type> gpz = exp(log_gpz);
    Type sigma_obs = exp(log_sigma_obs);
    Type sigma_age = exp(log_sigma_age);
    Type sigma_age2 = pow(log_sigma_age, 2.);
    Type rho = 1 / (1 + exp(-logit_rho));
    int N = log_rate_mort.size();
    int A = epsilon_age.size();
    printf("%s\n", "Data Transformed");
    
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
    /*for(int a=0; a<A; a++){
        Type dist;
        if((a == 0) & (option == 1)){
            nll -= dnorm(epsilon_age[a], Type(0.), sigma_age);
        }
        else if (option == 1){
            dist = unique_age[a]-unique_age[a-1];
            nll -= dnorm(epsilon_age[a], pow(rho,dist)*epsilon_age[a-1], 
                         pow(sigma_age2*(1.-pow(rho,2.*dist)), 0.5), true);
        }
    }*/

    using namespace density;
    if (option == 1){
        nll += SCALE( AR1(rho), pow(sigma_age2 / (1-pow(rho,2)),0.5))( epsilon_age );
    }

    printf("%s\n", "Setting age effects per observation");
    for(int n=0; n<N; n++){
        age_effect[n] = epsilon_age[age_group[n]];
    }
    
    printf("%s\n", "make_predictions");
    alpha = gpz[0];
    beta = gpz[1];
    gamma = gpz[2];
    delta = gpz[3];
    zeta = gpz[4];
    log_rate_mort_hat = log(alpha * exp(-1. * beta * age) + gamma +
        delta * exp(zeta * age)) + age_effect;

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
    REPORT(age_effect);
    REPORT(age_group);
    return nll;
    }
