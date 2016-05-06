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
    
    DATA_VECTOR(Y); // response
    DATA_VECTOR(A); // age
    DATA_VECTOR(P); // pointwise location
    DATA_MATRIX(X); // covariates
    DATA_INTEGER(re_option); // how to account for re structure

    PARAMETER(logit_rho); // 
    PARAMETER_VECTOR(beta); // weights
    PARAMETER_VECTOR(pi); // spatial errors
    PARAMETER(log_l0); // 
    PARAMETER(log_k); // 

    PARAMETER(log_sigma_p);
    PARAMETER(log_sigma_l);

    // Bookkeeping variables
    int N = Y.size();
    int K = X.cols();

    // exponentiate variables
    Type l0 = exp(log_l0);
    Type k = exp(log_k);
    Type sigma_p = exp(log_sigma_p);
    Type sigma_p2 = pow(sigma_p, 2.);
    Type sigma_l = exp(log_sigma_l);
    Type rho = 1 / (1 + exp(-logit_rho));

    vector<Type> log_l_inf(N);
    vector<Type> l_inf(N);
    vector<Type> l_inf_obs(N);
    vector<Type> y_pred(N);

    // record nll
    Type nll = 0.;

    // random effects structure
    Type dist;

    if(re_option == 0){
        for(int n = 0; n < N; n++){
            if(n == 0){
                nll -= dnorm(pi[n], Type(0.), sigma_p);
            }
            else{
                dist = P[n]-P[n-1];
                nll -= dnorm(pi[n], pow(rho,dist)*pi[n - 1], pow(sigma_p2*(1.-pow(rho,2.*dist)), 0.5), true);
            }
        }
    }

    if(re_option == 1){
        using namespace density;
        matrix<Type> cov_mat(N,N);
        for(int n1 = 0; n1 < N; n1++){
            for(int n2 = n1; n2 < N; n2++){
                dist = abs(P[n2]-P[n1]);
                cov_mat(n1,n2) = sigma_p2 * pow(rho, dist);
                if(n1!=n2) cov_mat(n2,n1) = cov_mat(n1,n2);
            }
        }
    nll += MVNORM(cov_mat)(pi);
    }
    
    if(re_option == 2){
        for(int n = 0; n < N; n++){
            if(n == 0){
                nll -= dnorm(pi[n], Type(0.), sigma_p);
            }
            else{
                dist = P[n]-P[n-1];
                nll -= dnorm(pi[n], rho*pi[n - 1], sigma_p, true);
            }
        }
    }
    
    // linear prediction
    for(int n = 0; n < N; n++){
        log_l_inf[n] = Type(0.); // initiate to zero
        log_l_inf[n] += pi[n]; // spatial random effect
        for(int k = 0; k < K; k++){
            log_l_inf[n] += beta[k] * X(n,k); // add in covariates separately
        }
    }
    
    l_inf = exp(log_l_inf);
    
    for(int n = 0; n < N; n++){
        y_pred[n] = l_inf[n] - ((l_inf[n] - l0) * exp(Type(-1.) * k * A[n]));
    }

    // data likelihood
    for(int n = 0; n < N; n++){
        nll -= dnorm(log(y_pred[n]), log(Y[n]), sigma_l, true);
    }

    // report the parameters I need
    REPORT(y_pred);
    REPORT(l_inf);
    REPORT(pi);
    REPORT(l0);
    REPORT(sigma_l);
    REPORT(sigma_p);
    REPORT(rho);
    REPORT(k);
    REPORT(beta);

    return nll;
}
