
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( l_i );
  DATA_VECTOR( a_i );
  DATA_VECTOR( loc_i );

  // Parameters
  PARAMETER( log_Linf_mean );
  PARAMETER( log_kappa );
  PARAMETER( log_Lzero );
  PARAMETER( log_rho );
  PARAMETER( beta ); // Spatial trend
  PARAMETER( ln_sigma_spatial );
  PARAMETER( ln_sigma_measurement );

  // Random effects
  PARAMETER_VECTOR( epsilon_i );

  // Objective funcction
  int n_i = l_i.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  Type sigma_measurement = exp(ln_sigma_measurement);
  Type sigma_spatial = exp( ln_sigma_spatial );
  Type rho = exp( log_rho );
  vector<Type> dist_i(n_i);
  Type kappa = exp(log_kappa);
  Type Linf_mean = exp(log_Linf_mean );
  Type Lzero = exp(log_Lzero );

  // Probability of random effects
  dist_i(0) = NA_REAL; // reminder that I don't need this slot
  jnll_comp(1) -= dnorm( epsilon_i(0), Type(0.0), sigma_spatial, true );
  for(int i=1; i<n_i; i++){
    dist_i(i) = loc_i(i) - loc_i(i-1);
    jnll_comp(1) -= dnorm( epsilon_i(i), pow(rho, dist_i(i))*epsilon_i(i-1), sigma_spatial*pow(1-pow(rho,2*dist_i(i)),0.5), true );
  }

  // Probability of data conditional on random effects
  vector<Type> lpred_i(n_i);
  vector<Type> Linf_i(n_i);
  for( int i=0; i<n_i; i++){
    Linf_i(i) = Linf_mean * exp(epsilon_i(i)) * exp(beta*loc_i(i));
    lpred_i(i) = Linf_i(i) - (Linf_i(i)-Lzero) * exp(-kappa*a_i(i));
    jnll_comp(0) -= dlognorm( l_i(i), log(lpred_i(i)), sigma_measurement, true );
  }

  // Objective function
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( beta );
  REPORT( epsilon_i );
  REPORT( dist_i );
  REPORT( lpred_i );
  REPORT( Linf_i );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( rho );
  REPORT( Linf_mean );
  REPORT( kappa );
  REPORT( Lzero );
  REPORT( sigma_spatial );
  REPORT( sigma_measurement );

  return jnll;
}
