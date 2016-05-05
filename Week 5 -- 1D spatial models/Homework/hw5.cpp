// State-space dynamic linear model
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
    Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
    if(give_log) return logres; else return exp(logres);
}


template<class Type>
Type objective_function<Type>::operator() () {
    
    DATA_VECTOR(Y);
    DATA_VECTOR(A);
    DATA_VECTOR(P);
    DATA_MATRIX(X);

    PARAMETER(logit_rho);
    PARAMETER_VECTOR(beta);
    PARAMETER_VECTOR(pi);
    PARAMETER_VECTOR(epsilon);
    PARAMETER(log_l0);
    PARAMETER(log_k);
    PARAMETER(log_sigma_e);
    PARAMETER(log_sigma_p2);
    PARAMETER(log_sigma_l);

    // Bookkeeping variables
    int N = Y.size();
    int K = X.cols();

    // exponentiate variables
    Type l0 = exp(log_l0);
    Type k = exp(log_k);
    Type sigma_p2 = exp(log_sigma_p2);
    Type sigma_e = exp(log_sigma_e);
    Type sigma_l = exp(log_sigma_l);
    Type rho = 1 / (1 + exp(-logit_rho));

    vector<Type> log_l_inf(N);
    vector<Type> l_inf(N);
    vector<Type> l_inf_obs(N);
    vector<Type> y_pred(N);

    // linear prediction
    for(int n = 0; n < N; n++){
        log_l_inf[n] = pi[n] + epsilon[n];
        for(int k = 0; k < K; k++){
            log_l_inf[n] += beta[k] * X(n,k);
        }
    }

    l_inf = exp(log_l_inf);

    for(int n = 0; n < N; n++){
        y_pred[n] = l_inf[n] - ((l_inf[n] - l0) * exp(Type(-1.) * k * A[n]));
    }

    // record nll
    Type nll = 0.;

    // random effects structure
    Type dist;
    for(int n = 0; n < N; n++){
        // nll -= dnorm(epsilon[n], Type(0.), sigma_e);
        if(n == 0){
            nll -= dnorm(pi[n], Type(0.), pow(sigma_p2, .5));
        }
        else{
            dist = P[n]-P[n-1];
            nll -= dnorm(pi[n], pow(rho,dist)*pi[n-1], pow(sigma_p2*(1-pow(rho,2*dist)), 0.5), true);
        }
    }

    // data likelihood
    for(int n = 0; n < N; n++){
        nll -= dnorm(log(y_pred[n]), log(Y[n]), sigma_l, true);
    }

    // report the parameters I need
    REPORT(y_pred);
    REPORT(l_inf);
    REPORT(epsilon);
    REPORT(sigma_l);
    REPORT(sigma_p2);
    REPORT(rho);

    return nll;
}
