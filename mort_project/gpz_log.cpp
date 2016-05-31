#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;

template<class Type>
Type logit_scaled(Type x, Type switch_point, Type scale){
    return 1. / (1. + exp(scale * (switch_point - x)));
}

template<class Type>
SparseMatrix<Type> car_Q(SparseMatrix<Type> graph, Type rho, Type sigma){
    SparseMatrix<Type> Q = rho * graph * Type(-1.);
    for (int i = 0; i < Q.rows(); i++){
       Q.coeffRef(i,i) = graph.col(i).sum() / sigma;
    }
}

template<class Type>
SparseMatrix<Type> ar_Q(int N, Type rho, Type sigma) {
    SparseMatrix<Type> Q(N,N);
    Q.insert(0,0) = (1.) / pow(sigma, 2.);
    for (size_t n = 1; n < N; n++) {
        Q.insert(n,n) = (1. + pow(rho, 2.)) / pow(sigma, 2.);
        Q.insert(n-1,n) = (-1. * rho) / pow(sigma, 2.);
        Q.insert(n,n-1) = (-1. * rho) / pow(sigma, 2.);
    }
    Q.coeffRef(N-1,N-1) = (1.) / pow(sigma, 2.);
    return Q; 
}

// template<class Type>
// SparseMatrix<Type> i_Q(int N, Type anything) { // that is a phat bug
//     SparseMatrix<Type> Q(N,N);
//     for (size_t n = 0; n < N; n++) {
//         Q.insert(n,n) = 1.;
//     }
//     return Q;
// }

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_VECTOR(log_rate_mort);
    DATA_VECTOR(age);
    DATA_IVECTOR(age_group);
    DATA_IVECTOR(time_group);
    DATA_IVECTOR(loc_group);
    DATA_INTEGER(option);
    DATA_SPARSE_MATRIX(graph); // queens matrix of location graph
    printf("%s\n", "Data loaded");

    // Parameters
    PARAMETER_VECTOR(gpz);
    PARAMETER(log_sigma_obs);
    PARAMETER_ARRAY(phi);
    PARAMETER_ARRAY(epsilon);
    PARAMETER(logit_rho_age);
    PARAMETER(logit_rho_time);
    PARAMETER(logit_rho_loc);
    PARAMETER(log_sigma_age);
    PARAMETER(log_sigma_time);
    PARAMETER(log_sigma_loc);
    PARAMETER(secular);
    printf("%s\n", "Parameters loaded");

    // Transforms
    Type sigma_obs = exp(log_sigma_obs);
    Type sigma_age = exp(log_sigma_age);
    Type sigma_time = exp(log_sigma_time);
    Type sigma_loc = exp(log_sigma_loc);
    Type rho_age = 1 / (1 + exp(-logit_rho_age));
    Type rho_time = 1 / (1 + exp(-logit_rho_time));
    Type rho_loc = 1 / (1 + exp(-logit_rho_loc));
    int N = log_rate_mort.size();
    int L = loc_group.maxCoeff() + 1;
    int A = age_group.maxCoeff() + 1;
    int T = time_group.maxCoeff() + 1;
    printf("%s\n", "Data Transformed");

    // Objective funcction
    Type nll = 0.;

    // Identity Matrix
    SparseMatrix<Type> Q_loc2(L,L);
    for (size_t n = 0; n < L; n++) {
        Q_loc2.insert(n,n) = 1.;
    }
    
    // Probability of random effects
    //SparseMatrix<Type> Q_loc = car_Q(graph, rho_loc, sigma_loc);
    SparseMatrix<Type> Q_age = ar_Q(A, rho_age, sigma_age);
    SparseMatrix<Type> Q_time = ar_Q(T, rho_time, sigma_time);
    
    if(option >= 1){
        nll += SEPARABLE(GMRF(Q_time), SEPARABLE(GMRF(Q_age), GMRF(Q_loc2)))(phi);
    }
    //if(option >= 2){
    //    nll += GMRF(Q_loc)(epsilon);
    //}

    printf("%s\n", "make_predictions");
    Type N0 = exp(gpz[0]);
    Type lambda = -1. * exp(gpz[1]);
    Type c = gpz[2];
    Type m = gpz[3];
    Type b = gpz[4];
    Type rho = exp(gpz[5]);
    printf("%s\n", "parameters gpz loaded");
    Type sw_term;
    Type inf_term;
    Type sns_term;
    Type re_term;
    vector<Type> log_rate_mort_hat(N);
    for(int n=0; n<N; n++){
        sw_term = logit_scaled(age[n], rho, Type(3.));
        inf_term = (N0 * exp(lambda * age[n]) + c);
        sns_term = (m * age[n] + b);
        re_term = phi(loc_group[n], age_group[n], time_group[n]); //+ epsilon(loc_group[n]);
        log_rate_mort_hat[n] = inf_term * (1 - sw_term) + sw_term * sns_term + 
            secular * time_group[n] + re_term; 
    }

    printf("%s\n", "evaluate data likelihood");
    for(int n=0; n<N; n++){
        nll -= dnorm(log_rate_mort[n], log_rate_mort_hat[n], sigma_obs, true);
    }

    printf("%s\n", "Report values");
    // Reporting
    REPORT(N0);
    REPORT(lambda);
    REPORT(c);
    REPORT(secular);
    REPORT(m);
    REPORT(b);
    REPORT(sigma_obs);
    REPORT(sigma_age);
    REPORT(nll);
    REPORT(age_group);
    REPORT(log_rate_mort_hat);
    //REPORT(Q_loc);
    REPORT(Q_age);
    REPORT(Q_time);
    REPORT(phi);
    REPORT(rho_age);
    return nll;
}
