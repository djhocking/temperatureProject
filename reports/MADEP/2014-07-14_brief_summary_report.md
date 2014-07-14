Massachusetts Daily Stream Temperature Modeling
===================

-> **Daniel Hocking, Ben Letcher, Kyle O'Neil** <-

-> *Silvio O. Conte Anadromous Fish Research Center, U.S. Geological Survey* <-

### Brief Summary Report: 14 July 2014

Our goal is to model daily stream temperature over large geographic areas as a function of air temperature, precipitation, and landscape characteristics. We used an index of air-water temperature synchrony to determine the part of the year where air temperature was directly affecting water temperature (i.e. minimal effects of ice-cover, phase change, and snow melt). We analyzed the synchronized period of each year at each site with a linear mixed effects model implemented in a Bayesian framework to be flexible and scaleable to large data sets. 

The model included air temperature, 1-day lagged air temperature, 2-day lagged air temperature, amount of precipitation that day and in each of the pervious two days (2 lags), drainage area, percent forest cover in the catchment, elevation of the stream reach, surficial coarseness of the catchment (how much sand, gravel, and rocks), percentage of the catchment that is comprised of wetland, area of stream impoundment, snow-water equivalent, latitude, longitude, and a cubic function of day of the year. We used year and each measurement site as random effects to account for correlation not explained by the other predictor variables and to adjust for unequal length time series at different sites and years. We modeled Massachusetts stream temperature data acquired from the MA Department of Environmental Protection, MA Division of Fisheries and Wildlife, and the U.S. Geological Survey.

Our model performed well with a root mean squared error of 0.89 <sup>o</sup>C, suggesting a typical accuracy within 1 <sup>o</sup>C. We found air temperature, forest cover, elevation, impoundments, and day of the year to be the most important predictors of stream temperature. There was also more correlation within sites than within years.

From this model, we will be able to predict daily stream temperature across time and space. This will allow us to calculate additional derived metrics of interest such as the average maximum summer temperature, number of days over a threshold temperature, frequency of thermally impaired days, resiliency to climate warming, etc. We can also compare the predicted stream temperatures for potential management actions (e.g. increased or decreased forest cover). The hierarchical mixed effects approach allows the incorporation of data from short monitoring periods and should provide accurate estimates even for sites that were monitored primarily during relatively extreme events (e.g. a once a decade heat wave). The model accuracy and generality will continue to improve as more data is incorporated. 

We are also investigating further enhancements to the model and incorporation of additional predictor variables. These include riparian forest cover, riparian impervious surfaces, distance to the nearest upstream dam, and various interactions among predictor variables.

### Figures

Examples of observed (blue) and predicted (red) stream temperature compared with air temperature (black).


![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 




