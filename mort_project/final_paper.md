# Forecasting US Mortality 
### Neal Marquez

## Modeling Mortality

Modeling human mortality has had a long history of attempts and revisions 
since Gompertz first made his claim of human mortality patterns in the 1880s.
The ability to accurately model and describe rates of mortality can have 
numerous effects on how social institutions and policies can be structured 
in order to meet the needs of its populations. For example knowing the age
pattern of mortality allows a society to build the medical infrastructure
to cope with the needs of different aged individuals. Knowing the root causes
of mortality allow for the change of policy in order to reduce specific 
drivers of mortality.

Another aspect of mortality that is often of concern is the future rates
and how they differ from the past. Over the past century there has been 
dramatic declines in mortality across all locations and ages. While that
progress has stagnated to some degree we still see dramatic decreases in 
mortality among developing countries where medical infrastructure is 
expanding and improving. As this happens we would expect there to be a
direct effects on the population, production and needs of a society.

## History of Descriptive Models

Descriptive models of mortality have an inherent age structure to them 
that allow for a simple descriptive explanation of how death effects 
individuals across their life-span. Gompertz first claim about human
mortality was that there is a linear increase in log rate mortality
as one aged. Below is an example of this phenomenon where we graph the 
log rate mortality of Women in Alabama between the ages of 22 and 77.

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/al_female_sns_mort.png ""){#id .class width=400 height=220px}\


While this model holds well for older ages, those who are part of the
senescence group, younger individuals experience quiet different 
mortality patterns as they age. From birth the rate of mortality
decreases until the senescence period is started at which time the 
mortality rates will then again increase. Below is the log rate 
mortality for males in Texas which shows this pattern that is nearly a human
universal. 

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/tx_male_all.png ""){#id .class width=400 height=220px}\

In order to capture this in 1983 Siler proposed a model to capture this switch
in mortality patterns by decompsing rate of mortality into three terms, an 
infant mortality term, a constant risk term and a senescence term.

$$
m_x = \alpha exp(\beta x) + c + \delta exp(\eta x)
$$

This model of mortality captures the change across ages well but still only 
applies to a single location of interest for a single snapshot in time. Others
have proposed more complex models, see Heligman-Pollard 1980, however as more
terms have been added to try and capture subtle nuances of mortality the 
generalizable usefulness of the model tends to suffer.

## Forecasting and Lee Carter

These models all offer a descriptive frame work for mortality over age however the do 
not offer a good solution for mortality over time or across regions. In 2000 
Lee & Carter developed a model that is till this day widely used for forecasting 
mortality at the all_cause and cause specific level. Abandoning the traditional
framework of looking at age patterns the model argues that better forecasts can be
made by assuming that ages are largely independent in their level and only similar in 
their rate of change. In doing so the correlation between ages is lost and as the 
forecasted model continues outward it will propagate differences between age groups.
The model is as follows

$$
m_{at} \sim \mathcal{N}(\mu_{at}, \sigma^{2})
$$
$$
\mu_{at} = \beta{a}\gamma{t}
$$
$$
\gamma_{t} = \gamma_{t-1} + \theta + \epsilon_{t}
$$
$$
\epsilon_{t} \sim \mathcal{N}(0, \sigma^{2}_{rw})
$$

While this model has preformed well when tested on US mortality data from 1970 
to 2005 within the US it has performed lackluster in other environments. 
Because of the lack of age structure th model produces nonsensical results 
where adjacent age groups have differing and sometimes opposite rates of 
change log rate mortality. in 2006 Girosi and King wrote a response to the model
showing where the model works well and the many times that it doesnt and 
criticizing the approach for not pooling information across age and 
geography.

## Modified GeoTemporal Siler

In order to fit all the dimensions of concern while still maintaining a coherent 
age structure this project attempts to use a modification of the Siler model 
which accounts for deviations away from the expected value due to relatedness 
across multiple dimensions while also including temporal change. The model can 
be broken down into three familiar components. The first component is a 
modification of the Siler model specified above and follows the form 

$$
S_x = (N0exp(\lambda x) + c) \times (1 - pr_x) + pr_x \times (mx + b)
$$
where
$$
pr_x = 1 / (1 + exp(3 * (\kappa-x)))
$$

The first term in the first equation corresponds to infant mortality term in
the Siler model and the second term is the senescence component. The model also
has a temporal trend added to it which accounts for some technological 
innovation over time.

$$
temporal = \beta \times time
$$

This component has the ability incorporate within it the effects of covariates
the mark medical and technological innovations that effect mortality however
we will simply use time in this exercise as a proxy for their effect on mortality.
That is to say that we expect that as time passes mortality will decrease as 
innovations happen. 

The last component is the structured random error component which captures
relatedness across three dimensions of error, age, space and time. The 
structure is as follows.

$$
\phi \sim \mathcal{N}(0, Q^{-1})
$$

$Q$ acts as the precision matrix for this model which is composed of three
independent portions.

$$
Q = Q_{loc} \otimes Q_{age} \otimes Q_{time}
$$

The age and the time precision components follow an AR1 process where 
the precision matrix is as follows

$$
\begin{aligned}
Q^{AR}_{i,j} = 
\begin{cases}
    \frac{1}{\sigma^2} ,& \text{if  } i = j = 0 | i = j = max(i) \\
    \frac{1 + \rho^2}{\sigma^2} ,& \text{else if  } i = j \\
    \frac{-\rho}{\sigma^2},  & \text{else if  } i \sim j \\
    0, & \text{otherwise}  
\end{cases}
 \end{aligned}
$$

For modeling purposes we will be using discrete units for bot age and time
where the ages reflect Global Burden of Disease age Groups and years are 
single year groups. This allows for a trivial application for the AR1 precision
matrix to be applied to each dimension.

The precision matrix for a the geospatial portion follows a conditional 
auto-regressive form where elements are autoregressive if they are 
considered neighbors. Because we are using an ariel approach we determine
neighbors as those geographical units that share a border with one another.
The matrix is defined as follows

$$
\begin{aligned}
Q^{CAR}_{i,j} = 
\begin{cases}
    \frac{1}{\sigma^2} ,& \text{if  } i = j \\
    \frac{-\rho}{\sigma},  & \text{else if  } i \sim j \\
    0, & \text{otherwise}  
\end{cases}
 \end{aligned}
$$

Restrictions to this precision matrix require that $\sigma$
be greater than zero and the matrix be diagonal dominant. That 
is that any diagonal element for any row must be greater than 
the sum of the absolute values of the off diagonal for that row.
This may be achieved, by making the maximum value of $pho$ be
slightly less than 1 over the maximum numbers of neighbors any one
location has. This may be achieved by a logit transform with a 
scale.

The full model is then seen as 

$$
log(\mu_{lat}) = S_x + temporal_t + \phi_{lat}
$$

with a probability distribution

$$
log(m_{lat}) \sim \mathcal{N}(\mu_{lat}, \sigma_{obs})
$$

## Fitting the model 

In order to test the model we will evaluate the model as fitted on 
US state mortality data between the years of 1990 to 2015. US data 
has the luxury of being from a near complete and comprehensive 
vital registration system that tracks most deaths at the state level 
for any given year well while documenting age. This data also has the 
benefit of being tested on many times such that comprehensive 
benchmarks exist for forecasting. 

To evaluate the model performance the model will be fit on the years
spanning between 1990 and 2005 and then error metrics will be 
calculated using holdout data from 2012 to 2015. After this full 
forecasts will be made out to 2030 using the entire data set. 
Separate models will be run for males and females.

## Results

Below is the fit of data between the years 1990 and 2005 
predicting from then on forward for the state of California. 
The second plot switches the age and time dimensions to show 
that as the model forecasts it maintains the gompertz age 
structure. The second plot is fitted using the full data set and
forecasted out to 2030.

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/ca_forecast_age.png ""){#id .class width=400 height=220px}\

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/ca_forecast_time.png ""){#id .class width=400 height=220px}\

Below is the error map shown with model for out of sample results
between 2013 and 2015.

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/log_rate_error_map.png ""){#id .class width=400 height=220px}\

And below is the difference in out of sample error within age groups.

![0](/home/neal/Documents/Classes/2016_Spatio-temporal_models/mort_project/forecasted_error_by_age.png ""){#id .class width=400 height=220px}\ 


## Discussion

While the model fits the data relatively well future changes to the exact parameterizations 
of the model will need to be done to account for the current biases in the model. 
Currently the model systematically is biased upwards for young age groups relative
to old age groups. This can be seen in the error by age group graph above for out 
of sample data. This could mean that changes to log rate mortality are changing more 
rapidly for young ages compared to older ages. 

In addition the map above shows that there may be some regional effects that are not 
captured well in this model. In order to alleviate this a random effect that just 
captures geography effect could be added. 

While this model was unable to come to the level of predictive validity of either the 
Girosi-King model or the Lee-Carter model, it retention of age pattern offers hope that 
it can be used in the future. In these runs of US dat the temporal trend for the model
was shown to be relatively weak however the US is a relatively stagnant country 
when it comes to decreases in log rate mortality over time and testing on other countries
is a must to get the full scope of the generalizability of this model. 

## References

Girosi, F., King, G.  (2006).  Demographic  Forecasting.  Cambridge  University  Press,  Cambridge.

Siler, W. (1983). Parameters of mortality in human populations with widely varying life spans. Statistics in Medicine, 2, 373-380. 

Lee,R.D., Carter, L.R.  (1992).  Modelling  and  forecasting  U.S.  mortality.  Journal  of  the American Statistical Association, 87(419), 659-671.

Heligman, L.,  Pollard, J.H.  (1980).  The  age  pattern  of  mortality.  Journal  of  the  Institute of Actuaries, 107(1, No 434), 49-80