#include <TMB.hpp>

template<class Type>
Type logit_scaled(Type x, Type switch_point, Type scale){
    return 1. / (1. + exp(scale * (switch_point - x)));
}

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
    //Type N0 = exp(gpz[0]);
    //Type lambda = -1. * exp(gpz[1]);
    Type c = gpz[0];
    Type eta = exp(gpz[1]);
    //Type eta = 16.;
    Type ceiling = exp(gpz[2]);
    Type scale = exp(gpz[3]);
    Type m = gpz[4];
    Type b = gpz[5];
    Type rho = eta + exp(gpz[6]);
    // Type rho = 26.;
    Type inf_term;
    Type ya_term;
    Type sns_prob;
    vector<Type> log_rate_mort_hat(N);
    for(int n=0; n<N; n++){
        // N0 * exp(lambda * age[n])
        log_rate_mort_hat[n] = c +  // infant term
            ceiling * logit_scaled(age[n], eta, scale) + //+ // young adult term
            logit_scaled(age[n], rho, Type(3.)) * (m * age[n] + b); // sns term
    }

    printf("%s\n", "evaluate data likelihood");
    for(int n=0; n<N; n++){
        nll -= dnorm(log_rate_mort[n], log_rate_mort_hat[n], sigma_obs, true);
    }

    printf("%s\n", "Report values");
    // Reporting
//    REPORT(N0);
//    REPORT(lambda);
    REPORT(c);
    REPORT(eta);
    REPORT(ceiling);
    REPORT(m);
    REPORT(b);
    REPORT(scale);
    REPORT(rho);
    REPORT(sigma_obs);
    REPORT(nll);
    REPORT(epsilon_age);
    REPORT(epsilon_time);
    REPORT(age_effect);
    REPORT(age_group);
    REPORT(log_rate_mort_hat);
    return nll;
}
