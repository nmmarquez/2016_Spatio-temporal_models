# HW\#2 Generalized Linear Mixed Models with TMB  
#### Neal Marquez

In order to test the effects of mis-specifying a model on a parameters 
standard error estimate we will simulate data that resembles collecting data 
on ten counts of a species at 10 different sites where sites have some 
biological difference in habitat suitability. The data follows the following 
distribution.

$$
log(\gamma_{s}) \sim \mathcal{N}(\mu, 1)
$$
$$
log(\lambda_{s,c}) \sim \mathcal{N}(log(\gamma_{s}), .5)
$$
$$
y_{s,c} \sim Poisson(\lambda_{s,c})
$$

$\mu$ is the log mean of the expected counts and has a set value of $2$.  
$log(\gamma_{s})$ is the log mean for site $s$.  
$\lambda_{s,c}$ is the expected count observation for site $s$ count $c$.  

We will attempt to estimate $\mu$ using four different models. 

1. a generalized linear model, with only an intercept term  __glm__
2. a GLMM with only among-site variability __glmm_site__
3. a GLMM with only overdispersion __glmm_ind__
4. a GLMM with both among-site variability and overdispersion __glmm_both__  

In order to get an idea of how often our confidence intervals for $\mu$ cover 
the true value we ill simulate 1000 data sets using the specifications above 
and observe (1) the range of the confidence intervals of our parameter of 
concern and (2) how often the confidence interval covers the true value.