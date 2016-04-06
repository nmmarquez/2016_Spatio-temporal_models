#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR(pred_bool);
  DATA_INTEGER(use_gamma);

  // Parameters
  PARAMETER_VECTOR(beta);
  PARAMETER(theta);
  PARAMETER(log_sigma);

  // Objective funcction
  Type zero_prob = 1 / (1 + exp(-theta));
  Type sigma = exp(log_sigma);
  int n_data = y.size();
  vector<Type> jnll_vec(n_data);
  Type jnll = 0;
  Type pred_jnll = 0;

  // Linear predictor
  vector<Type> exp_val = exp(X * beta);

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    if(y(i)==0) jnll_vec(i) -= log(zero_prob);
    if(y(i)!=0){
        if(use_gamma){
            Type shape = exp_val(i);
            jnll_vec(i) -= log(1-zero_prob) + dgamma(y(i), shape, sigma, true);
        }
        else{
            Type meanlog = log(exp_val(i));
            jnll_vec(i) -= log(1-zero_prob) + dlnorm(y(i), meanlog, sigma, true);
        }
    } 
    // Running counter
    if(pred_bool(i)==0) jnll += jnll_vec(i);
    if(pred_bool(i)==1) pred_jnll += jnll_vec(i);
  }

  // Reporting
  REPORT(zero_prob);
  REPORT(sigma);
  REPORT(exp_val);
  REPORT(pred_jnll);
  REPORT(jnll_vec);
  REPORT(theta);
  REPORT(beta);
  return jnll;
}
