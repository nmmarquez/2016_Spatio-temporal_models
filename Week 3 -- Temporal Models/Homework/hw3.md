# HW\#3 Temporal Models
### Neal Marquez

In order to understand the underlying process of modeling a temporally 
auto-regressive phenomena we will simulate data using the following 
structure.

$$
S_{t} \sim \mathcal{N}(\beta S_{t-1}, \sigma^{2}_{s})
$$
$$
Y_{t} \sim \mathcal{N}(S_{t}, \sigma^{2}_{y})
$$

Three parameterizations of $\beta$, $\beta = \{-.5, 0, .5\}$, and 
$\sigma_{y}$, $\sigma_{y} = \{.2, .4, .8\}$, will be used resulting in nine 
unique simulation sets where each set will consist of 100 simulation runs.
From these 100 simulation runs a mean, sd, 2.5, and 97.5 quantiles can be 
calculated and histograms of the parameter estimates can be observed.
Note that the temporal component of this model comes from the normal distribution 
around the prior state $S_{t-1}$. The observation that we see is a function of this 
past state, a growth or decay factor $\beta$ and some measurement error captured 
by $\sigma_{y}$.


## Parameter estimates from 9 unique simulations  

|    b| sigma_proc| sigma.obs|  b_hat| sigma_proc_hat| sigma_obs_hat|
|----:|----------:|---------:|------:|--------------:|-------------:|
| -0.5|        0.4|       0.2| -0.501|          0.388|         0.170|
|  0.0|        0.4|       0.2| -0.006|          0.270|         0.200|
|  0.5|        0.4|       0.2|  0.503|          0.379|         0.198|
| -0.5|        0.4|       0.4| -0.502|          0.377|         0.391|
|  0.0|        0.4|       0.4|  0.002|          0.310|         0.285|
|  0.5|        0.4|       0.4|  0.500|          0.386|         0.382|
| -0.5|        0.4|       0.8| -0.493|          0.386|         0.774|
|  0.0|        0.4|       0.8|  0.006|          0.511|         0.428|
|  0.5|        0.4|       0.8|  0.500|          0.301|         0.808|

|    b| sigma_proc| sigma.obs|  b_sd| sigma_proc_sd| sigma_obs_sd|
|----:|----------:|---------:|-----:|-------------:|------------:|
| -0.5|        0.4|       0.2| 0.036|         0.062|        0.111|
|  0.0|        0.4|       0.2| 0.051|         0.204|        0.209|
|  0.5|        0.4|       0.2| 0.037|         0.069|        0.100|
| -0.5|        0.4|       0.4| 0.044|         0.113|        0.119|
|  0.0|        0.4|       0.4| 0.055|         0.273|        0.268|
|  0.5|        0.4|       0.4| 0.043|         0.100|        0.111|
| -0.5|        0.4|       0.8| 0.071|         0.195|        0.123|
|  0.0|        0.4|       0.8| 0.085|         0.410|        0.424|
|  0.5|        0.4|       0.8| 0.060|         0.201|        0.104|

|    b| sigma_proc| sigma.obs|  b_q25| sigma_proc_q25| sigma_obs_q25|
|----:|----------:|---------:|------:|--------------:|-------------:|
| -0.5|        0.4|       0.2| -0.555|          0.272|         0.000|
|  0.0|        0.4|       0.2| -0.098|          0.000|         0.000|
|  0.5|        0.4|       0.2|  0.442|          0.226|         0.000|
| -0.5|        0.4|       0.4| -0.586|          0.156|         0.071|
|  0.0|        0.4|       0.4| -0.098|          0.000|         0.000|
|  0.5|        0.4|       0.4|  0.409|          0.174|         0.141|
| -0.5|        0.4|       0.8| -0.621|          0.000|         0.479|
|  0.0|        0.4|       0.8| -0.140|          0.000|         0.000|
|  0.5|        0.4|       0.8|  0.397|          0.000|         0.617|

|    b| sigma_proc| sigma.obs| b_q975| sigma_proc_q975| sigma_obs_q975|
|----:|----------:|---------:|------:|---------------:|--------------:|
| -0.5|        0.4|       0.2| -0.422|           0.479|          0.338|
|  0.0|        0.4|       0.2|  0.096|           0.485|          0.482|
|  0.5|        0.4|       0.2|  0.572|           0.483|          0.365|
| -0.5|        0.4|       0.4| -0.420|           0.582|          0.560|
|  0.0|        0.4|       0.4|  0.091|           0.632|          0.620|
|  0.5|        0.4|       0.4|  0.570|           0.556|          0.543|
| -0.5|        0.4|       0.8| -0.341|           0.699|          0.973|
|  0.0|        0.4|       0.8|  0.147|           0.966|          0.985|
|  0.5|        0.4|       0.8|  0.613|           0.664|          0.966|


The values of the 95% confidence intervals cover the true parameter in all 
simulation parameterizations. When $\beta$ is 0 the estimates of either $\sigma$ 
has a bimodal distribution. The histograms below show this effect

### Estimates of $\sigma_{s}$
![Estimates of $\sigma_{s}$](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week\ 3\ --\ Temporal\ Models/Homework/plot9.png  ""){#id .class width=495 height=295px}\  

### Estimates of $\sigma_{y}$
![Estimates of $\sigma_{y}$](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week\ 3\ --\ Temporal\ Models/Homework/plot9_sigma_proc.png  ""){#id .class width=495 height=295px}\  

This is likely do to the inability of the model to distinguish where the 
variance is coming from when $\beta$ is zero since there is no 
auto-regression. 

## Gompertz model simulation  

A Gompertz model can be defined as the following  
$$
S_{t} \sim \mathcal{N}(\alpha + \beta S_{t-1}, \sigma^{2}_{s})
$$
$$
Y_{t} \sim \mathcal{N}(S_{t}, \sigma^{2}_{y})
$$

The primary difference is the inclusion of the intercept term $\alpha$ which affects the trend of the population. 
In this model a population will converge on a stable state. We will simulate 
data following this structure with $\alpha = 1.2$, $\beta = .7$, 
$\sigma_{s} = .2$, $\sigma_{y} = .5$. We will start with an initial value of 
4 for the first observation $y$ which leads to a difficult estimation of the 
variance terms $\sigma_{s}$ & $\sigma_{y}$ and even the means of $\alpha$ and 
$\beta$. This is shown in the very non-normal distributions of the parameter 
estimates and lack of convergence to the mean.
 
### Estimates of Gompertz parameters when start values are at stable point 
![Estimates of Gompertz parameters](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week\ 3\ --\ Temporal\ Models/Homework/plot4.png  ""){#id .class width=495 height=295px}\  

This effect is attenuated when we change the start value. So that the model may more readily identify 
the parameters within 100 time points.

### Estimates of Gompertz parameters when start values are not at stable point  
![Estimates of Gompertz parameters](/home/neal/Documents/Classes/2016_Spatio-temporal_models/Week\ 3\ --\ Temporal\ Models/Homework/plot4b.png  ""){#id .class width=495 height=295px}\  