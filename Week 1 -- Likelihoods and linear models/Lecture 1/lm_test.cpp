
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);

  // Parameters
  PARAMETER(mean);
  PARAMETER(log_sd);

  // Objective funcction
  Type sd = exp(log_sd);

  // Probability of data conditional on fixed effect values
  Type jnll = -1 * sum(dnorm(y, mean, sd, true));
  
  // Reporting
  return jnll;
}

