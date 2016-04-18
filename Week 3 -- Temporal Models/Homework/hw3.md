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