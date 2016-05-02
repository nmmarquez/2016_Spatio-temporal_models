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
    
    PARAMETER_VECTOR(l_inf);
    PARAMETER(l0);
    PARAMETER(k);
    PARAMETER(log_sigma);
    
    int N y.size();
    
    yhat vector<Type> (N);
    for(int i = 0; i < n; i++){
        if (l_inf.size() == 1){
            yhat[i] = exp(l_inf[0] - (l_inf[0] - l0) * exp(-k * A[i]));
        }
        else{
            yhat[i] = exp(l_inf[i] - (l_inf[i] - l0) * exp(-k * A[i]));
        }
    }
    
}
