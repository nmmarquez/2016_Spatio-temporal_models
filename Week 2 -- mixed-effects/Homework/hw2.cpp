#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
    printf ("%s \n", "Testing the print function inside of TMB C++");
    // Data
    DATA_ARRAY(y);
    DATA_INTEGER(site_var);
    DATA_INTEGER(ind_var);
    
    // Parameters
    PARAMETER(x0);
    PARAMETER_VECTOR(log_sdz);
    PARAMETER_VECTOR(log_sde);
    PARAMETER_VECTOR(z);
    PARAMETER_ARRAY(e);
    
    printf ("%s \n", "Parameters loaded");
    
    // Objective funcction
    Type jnll = 0;
    
    // dims 
    int site_obs = y.dim(0);        // number of observations per site
    int num_sites = y.dim(1);        // number of sites
    
    printf ("%s \n", "dimensions found");
    
    // Probability of data conditional on fixed and random effect values
    array<Type> ypred(site_obs, num_sites);
    for(int i=0; i<site_obs; i++){
        for(int j=0; j<num_sites; j++){
            if(site_var and ind_var){
                ypred(i,j) = exp(x0 + z(j) + e(i,j));
            }
            if(site_var and not ind_var){
                ypred(i,j) = exp(x0 + z(j));
            }
            if(not site_var and ind_var){
                ypred(i,j) = exp(x0 + e(i,j));
            }
            if(not site_var and not ind_var){
                ypred(i,j) = exp(x0);
            }
        jnll -= dpois(y(i,j), ypred(i,j), true);
        }
    }
    
    // Probability of random coefficients
    Type sdz;
    if(site_var){
        for( int j=0; j<num_sites; j++){
            jnll -= dnorm(z(j), Type(0.0), exp(log_sdz[0]), true);
        }
        sdz = exp(log_sdz[0]);
    }
    Type sde;
    if(ind_var){
        for(int i=0; i<site_obs; i++){
            for( int j=0; j<num_sites; j++){
                jnll -= dnorm(e(i,j), Type(0.0), exp(log_sde[0]), true);
            }
        }
        sde = exp(log_sde[0]);
    }
    
    // Reporting
    REPORT(sdz);
    REPORT(z);
    REPORT(e);
    ADREPORT(e);
    ADREPORT(sdz);
    REPORT(sde);
    ADREPORT(sde);
    ADREPORT(z);
    REPORT(x0);
    ADREPORT(x0);
    
    return jnll;
}
