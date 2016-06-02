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
across multiple dimensions while also including temporal change.

