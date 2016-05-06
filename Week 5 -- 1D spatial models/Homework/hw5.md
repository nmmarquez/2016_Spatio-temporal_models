# HW\#2 1D spatial models using TMB 
#### Neal Marquez

## Simulating Fish growth rates

We believe that fish growth patterns follow a similar patterns to that 
described by the von Bertalanffy growth curve

$$
\hat{l}(a) = l_{\infty} - (l_{\infty} - l_{0})exp(-ka) 
$$
$$
l(a) \sim \mathcal{N}(\hat{l}(a), \sigma_{l})
$$

Where $a$ is the age of the fish, $l_{\infty}$ is the maximum length of a fish, 
l_{0} is the age of a fish at birth, $k$ is the growth rate of the fish along 
the curve, and $l(a)$ is the observed length of a fish given its age. In 
addition $l_{\infty}$ varies along the coastline with a general trend and via 
some unobserved spatial variation due to unmeasured environmental factors.

## Modeling unobserved $l_{\infty}$ 

I order to predict $l_{\infty}$ which is unobserved we will use observed fish 
lengths in order to estimate the value. We assume that $l_{\infty}$ follows a
1-D spatial process described be either 

(1) $l_{\infty} \sim \mathcal{N}(\beta_{0}, \Sigma)$

or 

(2) $l_{\infty} \sim \mathcal{N}(\beta_{0} + \beta_{1}s, \Sigma)$

Where $\beta_{0}$ is an intercept term and $\beta_{1}$ is the spatial trend 
across the coastline. $\Sigma$ is characterized by two parameters $\rho$ and
$\sigma_{p}$ which capture both the pointwise variation and the geostatistical 
range.

## Simulation process 

1000 observations of fish along the coast were simulated using a Random Gaussian 
field in order to simulate geospatial random error to which weights were then 
applied via an intercept and location ($\beta_0$ and $\beta_1$). From this 
$\hat{l}(a)$ was calculated using a distribution of ges that resemble a fish
mortality and population patterns. $l(a)$ was drawn from a normal distribution 
with mean $\hat{l}(a)$ and standard deviation $\sigma_{l}$. In the estimation 
process we will use simulated values of $l(a)$ in order to try and recover the 
true values of $l_{\infty}$ that were used in order to generate those values.
This process was repeated 100 times for validation.

![Simulated values of $l_{\infty}$](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework/Rplot.png  ""){#id .class width=495 height=295px} 

## Evaluation
RMSE was used to evaluate the two models. Estimated values with model (2).

![Estimated values of $l_{\infty}$](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework/Rplot01.png  ""){#id .class width=495 height=295px}\

When using covariates to inform the model for estimates of $l_{\infty}$ the model was able to better recover 
$l_{\infty}$ when the model was able to successfully converge. Both models had difficulty converging in some positions
which may have to do with the simulation process not directly resembling the simulation process. This was true even 
when applying an equidistant AR model. This seems like it my be because the model simulation produces greater 
uncertainty farther away from the origin. In future simulations it may be advantages to explicitly write out the 
simulation.