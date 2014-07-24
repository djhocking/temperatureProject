rm(list=ls())

#library(devtools)
# devtools::install_github("hadley/dplyr") # not working

library(ggplot2)
library(dplyr)
library(DataCombine) # for the slide function

#setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
#baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
baseDir <- 'C:/Users/dhocking/Documents/temperatureProject/'
setwd(baseDir)

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
dataLocalDir <- paste0(baseDir, 'localData/')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))
source(paste0(baseDir, 'code/functions/dataIndexingFunctions.R'))

load((paste0(dataOutDir, 'modSummary.RData')))
load(paste0(dataLocalDir, 'daymetFullRecordObservedMASites.RData'))
load(paste0(dataOutDir, 'tempDataSync.RData'))

load(paste0(dataOutDir, 'springFallBreakpoints.RData'))

#Northeast
CTDEP  <- T
MAFW   <- T
MAUSGS <- T
MADEP  <- T 
NHFG   <- F
NHDES  <- F
USFS   <- F
VTFWS  <- F
MEDMR  <- F

#Montana
MTUSGSYellowstone <- F
MTUSGSGlacier <- F

sourceChoice <- list( CTDEP,   MAFW,   MAUSGS, MADEP,   NHFG,   NHDES,   MEDMR,   USFS,   VTFWS,    MTUSGSYellowstone,   MTUSGSGlacier)
sourceNames  <- c   ('CTDEP', 'MAFW', 'MAUSGS', 'MADEP', 'NHFG', 'NHDES', 'MEDMR', 'USFS', 'VTFWS',  'MTUSGSYellowstone', 'MTUSGSGlacier')

dataSource <- sourceNames[sourceChoice == T]

fields <- c("agency", "date", "AgencyID", "year", "site", "date", "dOY", "temp", "airTemp", "prcp", "srad", "dayl", "swe")

covariateData <- readStreamTempData(timeSeries=FALSE, covariates=TRUE, dataSourceList=dataSource, fieldListTS=fields, fieldListCD='ALL', directory=dataInDir)
springFallBPs$site <- as.character(springFallBPs$site)


########## How to add BP for years without data and clip data to the sync period ??? #######
# Join with break points
covariateDataBP <- left_join(covariateData, springFallBPs, by=c('site', 'year'))
# rm(covariateData)

# temp hack
climateData$site <- as.character(climateData$site)
tempFullSync <- left_join(climateData, covariateData, by=c('site'))

# Clip to syncronized season
# tempFullSync <- filter(tempDataBP, dOY >= finalSpringBP & dOY <= finalFallBP)

# temp hack
tempFullSync <- filter(tempFullSync, dOY >= 50 & dOY <= 350)
tempFullSync$Latitude <- tempFullSync$Latitude.x
tempFullSync$Longitude <- tempFullSync$Longitude.x
##################

# Order by group and date
tempFullSync <- tempFullSync[order(tempFullSync$site,tempFullSync$year,tempFullSync$dOY),]

# For checking the order of tempFullSync
tempFullSync$count <- 1:length(tempFullSync$year)

tempFullSync <- tempFullSync[order(tempFullSync$count),] # just to make sure tempFullSync is ordered for the slide function

# airTemp
tempFullSync <- slide(tempFullSync, Var = "airTemp", GroupVar = "site", slideBy = -1, NewVar='airTempLagged1')
tempFullSync <- slide(tempFullSync, Var = "airTemp", GroupVar = "site", slideBy = -2, NewVar='airTempLagged2')

# prcp
tempFullSync <- slide(tempFullSync, Var = "prcp", GroupVar = "site", slideBy = -1, NewVar='prcpLagged1')
tempFullSync <- slide(tempFullSync, Var = "prcp", GroupVar = "site", slideBy = -2, NewVar='prcpLagged2')
tempFullSync <- slide(tempFullSync, Var = "prcp", GroupVar = "site", slideBy = -3, NewVar='prcpLagged3')

tempFullSync <- left_join(tempFullSync, tempDataSync[ , c("site", "date", "temp")], by = c("site", "date"))

# Make dataframe with just variables for modeling and order before standardizing
tempFullSync <- tempFullSync[ , c("year", "site", "date",  "FEATUREID", "HUC4", "HUC8", "HUC12", "temp", "Latitude", "Longitude", "airTemp", "airTempLagged1", "airTempLagged2", "prcp", "prcpLagged1", "prcpLagged2", "prcpLagged3", "dOY", "Forest", "Herbacious", "Agriculture", "Developed", "TotDASqKM", "ReachElevationM", "ImpoundmentsAllSqKM", "HydrologicGroupAB", "SurficialCoarseC", "CONUSWetland", "ReachSlopePCNT", "srad", "dayl", "swe")] #  "finalSpringBP", "finalFallBP", "agency", ""date","fsite", "fyear", "AgencyID" 

summary(tempFullSync)
dim(tempFullSync)
#tempFullSync <- na.omit(tempFullSync) ####### Change this so don't take out NA in stream temperature
dim(tempFullSync)

# Standardize for Analysis

stdCovs <- function(x, y, varNames){
  xStd <- as.data.frame(matrix(NA, dim(x)[1], length(varNames)))
  names(xStd) <- varNames
  for(i in 1:length(varNames)){
    xStd[ , varNames[i]] <- (x[ , varNames[i]] - mean(y[ , varNames[i]], na.rm=T)) / sd(y[ , varNames[i]], na.rm=T)
  }
  return(xStd)
}

varNames1 <- names(tempFullSync[ ,9:dim(tempFullSync)[2]])

tempFullStd <- stdCovs(tempFullSync, tempDataSync, varNames1)
tempFullSyncS <- cbind(tempFullSync[ ,c(1:8)], tempFullStd)

summary(tempFullSyncS)
tempFullSyncS[is.na(tempFullSyncS)] <- 0

fixEf <- modSummary@fixEf[ ,"Mean"]
names(fixEf) <- row.names(modSummary@fixEf)

#tempFullSync <- tempFullSync[which(tempFullSync$site %in% unique(tempDataSync$site)), ]
#tempFullSyncS <- tempFullSync[which(tempFullSyncS$site %in% unique(tempDataSync$site)), ]
sites <- unique(tempFullSync$site)
BSite <- modSummary@BSite
BYear <- modSummary@BYear

tempFullSyncS$cYear <- as.character(tempFullSyncS$year)


# Split data by site-year then do predictions for those with observed stream temperature data and those without, then recombine. The problem is that sites outside of the years observed won't get the site-specific values and years with data but at different sites won't get the site-specific data.
tempFullSyncS$siteYear <- paste0(tempFullSyncS$site, tempFullSyncS$year)
tempDataSyncS$siteYear <- paste0(tempDataSyncS$site, tempDataSyncS$year)

tempFullSiteYearS <- tempFullSyncS[which(tempFullSyncS$siteYear %in% unique(tempDataSyncS$siteYear)), ]
tempFullMeanS <- subset(tempFullSyncS, !(tempFullSyncS$siteYear %in% unique(tempDataSyncS$siteYear)))

# this will work fo MA because not predicting to any completely new sites
tempFullSiteYearS <- filter(tempFullSyncS, filter = year %in% unique(tempDataSyncS$year))
tempFullSiteS <- filter(tempFullSyncS, filter = !(year %in% unique(tempDataSyncS$year)))

tempFullSiteYearS$tempPredicted <- modSummary@fixEf["intercept", "Mean"] +
  BSite[tempFullSiteYearS$site, "intercept.site"] + 
  BYear[tempFullSiteYearS$cYear, "intercept.year"] + 
  modSummary@fixEf["lat", "Mean"]*tempFullSiteYearS$Latitude + 
  modSummary@fixEf["lon", "Mean"]*tempFullSiteYearS$Longitude + 
  BSite[tempFullSiteYearS$site, "airTemp"]*tempFullSiteYearS$airTemp + 
  BSite[tempFullSiteYearS$site, "airTempLag1"]*tempFullSiteYearS$airTempLagged1 + 
  BSite[tempFullSiteYearS$site, "airTempLag2"]*tempFullSiteYearS$airTempLagged2 + 
  BSite[tempFullSiteYearS$site, "precip"]*tempFullSiteYearS$prcp + 
  BSite[tempFullSiteYearS$site, "precipLag1"]*tempFullSiteYearS$prcpLagged1 + 
  BSite[tempFullSiteYearS$site, "precipLag2"]*tempFullSiteYearS$prcpLagged2 + 
  BSite[tempFullSiteYearS$site, "drainage"]*tempFullSiteYearS$TotDASqKM + 
  BSite[tempFullSiteYearS$site, "forest"]*tempFullSiteYearS$Forest + 
  BSite[tempFullSiteYearS$site, "elevation"]*tempFullSiteYearS$ReachElevationM + 
  BSite[tempFullSiteYearS$site, "coarseness"]*tempFullSiteYearS$SurficialCoarseC + 
  BSite[tempFullSiteYearS$site, "wetland"]*tempFullSiteYearS$CONUSWetland + 
  BSite[tempFullSiteYearS$site, "impoundments"]*tempFullSiteYearS$ImpoundmentsAllSqKM + 
  BSite[tempFullSiteYearS$site, "swe"]*tempFullSiteYearS$swe + 
  BYear[tempFullSiteYearS$cYear, "dOY"]*tempFullSiteYearS$dOY + 
  BYear[tempFullSiteYearS$cYear, "dOY2"]*((tempFullSiteYearS$dOY)^2) + 
  BYear[tempFullSiteYearS$cYear, "dOY3"]*((tempFullSiteYearS$dOY)^3)

tempFullSiteS$tempPredicted <- modSummary@fixEf["intercept", "Mean"] +
  BSite[tempFullSiteS$site, "intercept.site"] + 
  modSummary@fixEf["lat", "Mean"]*tempFullSiteS$Latitude + 
  modSummary@fixEf["lon", "Mean"]*tempFullSiteS$Longitude + 
  BSite[tempFullSiteS$site, "airTemp"]*tempFullSiteS$airTemp + 
  BSite[tempFullSiteS$site, "airTempLag1"]*tempFullSiteS$airTempLagged1 + 
  BSite[tempFullSiteS$site, "airTempLag2"]*tempFullSiteS$airTempLagged2 + 
  BSite[tempFullSiteS$site, "precip"]*tempFullSiteS$prcp + 
  BSite[tempFullSiteS$site, "precipLag1"]*tempFullSiteS$prcpLagged1 + 
  BSite[tempFullSiteS$site, "precipLag2"]*tempFullSiteS$prcpLagged2 + 
  BSite[tempFullSiteS$site, "drainage"]*tempFullSiteS$TotDASqKM + 
  BSite[tempFullSiteS$site, "forest"]*tempFullSiteS$Forest + 
  BSite[tempFullSiteS$site, "elevation"]*tempFullSiteS$ReachElevationM + 
  BSite[tempFullSiteS$site, "coarseness"]*tempFullSiteS$SurficialCoarseC + 
  BSite[tempFullSiteS$site, "wetland"]*tempFullSiteS$CONUSWetland + 
  BSite[tempFullSiteS$site, "impoundments"]*tempFullSiteS$ImpoundmentsAllSqKM + 
  BSite[tempFullSiteS$site, "swe"]*tempFullSiteS$swe + 
  modSummary@fixEf["dOY", "Mean"]*tempFullSiteS$dOY + 
  modSummary@fixEf["dOY2", "Mean"]*((tempFullSiteS$dOY)^2) + 
  modSummary@fixEf["dOY3", "Mean"]*((tempFullSiteS$dOY)^3)


tempFullS <- rbind(tempFullSiteYearS, tempFullSiteS)

tempFull <- left_join(tempFullSync, tempFullS[ , c("year", "site", "date", "tempPredicted")], by = c("year", "site", "date")) 


# plot observed and predicte vs day of the year for all sites in all years
sites <- unique(as.character(tempFull$site))

for(i in 1:length(unique(tempFull$site))){
  dataSite <- filter(tempFull, filter = site == sites[i])
  dataSiteObs <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, tempPredicted)) + 
    coord_cartesian(xlim = c(100, 300), ylim = c(0, 35)) + 
    geom_point(data=dataSiteObs, aes(dOY, temp), colour='blue') +
    geom_point(colour = 'red', size=1) + 
    geom_line(colour = 'red', size=0.1) + 
    geom_point(aes(dOY, airTemp), size=1) + 
    ggtitle(dataSite$site[i]) + 
    facet_wrap(~year) + 
    xlab(label = 'Day of the year') + ylab('Temperature (C)') + 
    theme(axis.text.x = element_text(angle = 45))
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/fullRecord/', dataSite$site[i], '.png'), plot=foo, dpi=300 , width=12,height=8, units='in' )
} # surprisingly fast but wouldn't do for all catchments

yearPredict <- filter(tempFull, site == "MADEP_W0989_T1", year == "2005")
dataSiteObs <- filter(tempDataSync, filter = site == "MADEP_W0989_T1")
foo <- ggplot(yearPredict, aes(dOY, tempPredicted)) + 
  coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + 
  geom_point(data=dataSiteObs, aes(dOY, temp), colour='blue') +
  geom_point(colour = 'red', size=1) + 
  geom_line(colour = 'red', size=0.1) + 
  geom_point(aes(dOY, airTemp), size=1) + 
  #ggtitle(dataSite$site[i]) + 
  facet_wrap(~year) + 
  xlab(label = 'Day of the year') + ylab('Temperature (C)') + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(filename="C:/Users/dhocking/Documents/temperatureProject/presentations/yearPredict.png", plot=foo, dpi=300, width=12, height=8, units="in")


yearPredict <- filter(tempFull, site == "MADEP_W0989_T1", year > 2001 & year <= 2013)
dataSiteObs <- filter(tempDataSync, filter = site == "MADEP_W0989_T1")
foo <- ggplot(yearPredict, aes(dOY, tempPredicted)) + 
  coord_cartesian(xlim = c(100, 300), ylim = c(0, 30)) + 
  geom_point(data=dataSiteObs, aes(dOY, temp), colour='blue') +
  geom_point(colour = 'red', size=1) + 
  geom_line(colour = 'red', size=0.1) + 
  geom_point(aes(dOY, airTemp), size=1) + 
  #ggtitle(dataSite$site[i]) + 
  facet_wrap(~year) + 
  xlab(label = 'Day of the year') + ylab('Temperature (C)') + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(filename="C:/Users/dhocking/Documents/temperatureProject/presentations/multiYearPredict.png", plot=foo, dpi=300, width=12, height=8, units="in")




# plot observed and predicte vs day of the year for all sites
sites <- unique(tempDataSync$site)

for(i in 1:length(unique(tempDataSync$site))){
  dataSiteObs <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSiteObs, aes(dOY, temp)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + geom_point(aes(dOY, streamTempPred), colour = 'red', size=1) + geom_line(aes(dOY, streamTempPred), colour = 'red', size=0.1) + geom_point(aes(dOY, airTemp), colour='black', size=1) + ggtitle(unique(tempDataSync$fsite)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast

rmse(tempFull[which(!is.na(tempFull$temp)), "temp"] - tempFull[!is.na(tempFull$temp), "tempPredicted"])


############## Derived metrics ##########

# Mean maximum daily mean temperature by site (over years)
bySite <- group_by(tempFull, site)
bySiteYear <- group_by(bySite, year, add = TRUE)
maxTemp <- filter(bySite, tempPredicted == max(tempPredicted))
maxTempSite <- summarise(maxTemp, mean(tempPredicted)) # not needed - already max.t
#summarise(by.site.year, sd(mean(tempPredicted))) # not working based on filter or grouping

(maxTempSiteYear <- summarise(bySiteYear, max(tempPredicted)))
names(maxTempSiteYear) <- c("site", "year", "maxTempPredicted")
derivedSiteMetrics <- summarise(maxTempSiteYear, meanMaxTemp = mean(maxTempPredicted))
# maxTempSiteYear1 <- left_join(as.data.frame(maxTempSiteYear), tempFull, by=c("site", "tempPredicted"))

# Maximum max daily mean temperature
maxMaxTemp <- bySiteYear %>%
  summarise(maxTemp = max(tempPredicted)) %>%
  summarise(maxMaxTemp = max(maxTemp))

derivedSiteMetrics <- left_join(derivedSiteMetrics, maxMaxTemp, by = "site")

# ggplot(tempFull, aes(dOY, temp)) + geom_point(size=1, colour='black') + geom_point(aes(dOY, tempPredicted), colour = 'red', size=0.75) + ylab(label="Stream temperature (C)") + xlab("Day of the year") + geom_point(data=maxTempSiteYear1, aes(dOY, tempPredicted), colour = "green") + facet_grid(site ~ year) # max temp points all replicated on every panel

# Number of days with stream temp > 18C
meanDays18 <- bySiteYear %>%
  filter(tempPredicted > 18) %>%
  summarise(days18 = n()) %>%
  summarise(meanDays18 = mean(days18))

derivedSiteMetrics <- left_join(derivedSiteMetrics, meanDays18, by = "site")

# Number of years with mean maximum over 18 C
yearsMaxTemp18 <- summarise(
  filter(summarise(bySiteYear, maxTemp = max(tempPredicted)), maxTemp > 18),
  yearsMaxTemp18 = n()
)
derivedSiteMetrics <- left_join(derivedSiteMetrics, yearsMaxTemp18, by = "site")

# frequency of years with a mean max over 18 C
derivedSiteMetrics <- mutate(derivedSiteMetrics, freqMax18 = yearsMaxTemp18/length(unique(bySiteYear$year)))

# Resistance to peak air temperature
## Need to think of a way to make more general rather than by specific dOY (60 day max moving window air temp?)
meanResist <- bySiteYear %>%
  filter(dOY >= 145 & dOY <= 275) %>%
  mutate(absResid = abs(airTemp - tempPredicted)) %>%
  summarise(resistance = sum(absResid)) %>%
  summarise(meanResist = mean(resistance))

derivedSiteMetrics <- left_join(derivedSiteMetrics, meanResist, by = "site")


WB.2011.summer <- tempFull[which(tempFull$site == "MAUSGS_WEST_BROOK" & tempFull$year == 2011 & tempFull$dOY >=145 & tempFull$dOY <= 275), ]
sum(WB.2011.summer$airTemp - WB.2011.summer$tempPredicted)

ggplot(tempFull[which(tempFull$site == "MAUSGS_WEST_BROOK" & tempFull$year == 2011), ], aes(dOY, tempPredicted)) + 
  geom_point(size=2, colour = "red") + geom_line(colour = 'red') +
  geom_point(data=tempFull[which(tempFull$site == "MAUSGS_WEST_BROOK" & tempFull$year == 2011), ], aes(dOY, airTemp), colour = "black", size=2) + 
  geom_line(data=tempFull[which(tempFull$site == "MAUSGS_WEST_BROOK" & tempFull$year == 2011), ], aes(dOY, airTemp), colour = "black") + 
  geom_ribbon(data = tempFull[which(tempFull$site == "MAUSGS_WEST_BROOK" & tempFull$year == 2011 & tempFull$dOY >=145 & tempFull$dOY <= 275), ], aes(x=dOY, ymin=tempPredicted, ymax=airTemp), fill="dark grey", alpha=.5) +
  xlab("Day of the year") +
  ylab("Temperature (C)") #+ theme_classic()

ggplot(tempFull[which(tempFull$site == "WB OBEAR" & tempFull$year == 2010), ], aes(dOY.real, tempPredicted)) + 
  geom_point(size=2, colour = "black") + geom_line(colour = 'black') +
  geom_abline(intercept = 18, slope=0, colour='red') +
  geom_point(data = tempFull[which(tempFull$site == "WB OBEAR" & tempFull$year == 2010 & tempFull$tempPredicted >= 18), ], aes(dOY.real, tempPredicted), colour='red') +
  xlab("Day of the year") +
  ylab("Stream temperature (C)") #+ theme_classic()

# Reset ggplot2 theme default to gray
theme_set(theme_gray())



# Air-Water Resiliency

# RMSE for each site (flag highest)
meanRMSE <- bySiteYear %>%
  filter(!(is.na(temp))) %>%
  mutate(error2 = (temp - tempPredicted)^2) %>%
  summarise(RMSE = sqrt(mean(error2))) %>%
  summarise(meanRMSE = mean(RMSE))

derivedSiteMetrics <- left_join(derivedSiteMetrics, meanRMSE, by = "site")

# total observations (days with data) per site
totObs <- bySiteYear %>%
  filter(!is.na(temp)) %>%
  summarise(Obs = n()) %>%
  summarise(totObs = sum(Obs))

derivedSiteMetrics <- left_join(derivedSiteMetrics, totObs, by = "site")

# Flag based on RMSE > 90%
derivedSiteMetrics <- mutate(derivedSiteMetrics, flag = ifelse(meanRMSE > quantile(derivedSiteMetrics$meanRMSE, probs = c(0.9), na.rm=TRUE), "Flag", ""))

summary(derivedSiteMetrics)

derivedSiteMetricsClean <- na.omit(derivedSiteMetrics)
write.table(derivedSiteMetricsClean, file = 'reports/MADEP/derivedSiteMetrics.csv', sep=',', row.names = F)



