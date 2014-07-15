rm(list=ls())

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

load(paste0(dataOutDir, 'springFallBreakpoints.RData'))

#Northeast
CTDEP  <- F
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


# Make dataframe with just variables for modeling and order before standardizing
tempFullSync <- tempFullSync[ , c("year", "site", "date",  "FEATUREID", "HUC4", "HUC8", "HUC12", "Latitude", "Longitude", "airTemp", "airTempLagged1", "airTempLagged2", "prcp", "prcpLagged1", "prcpLagged2", "prcpLagged3", "dOY", "Forest", "Herbacious", "Agriculture", "Developed", "TotDASqKM", "ReachElevationM", "ImpoundmentsAllSqKM", "HydrologicGroupAB", "SurficialCoarseC", "CONUSWetland", "ReachSlopePCNT", "srad", "dayl", "swe")] #  "finalSpringBP", "finalFallBP", "agency", ""date","fsite", "fyear", "AgencyID","temp", 

summary(tempFullSync)
dim(tempFullSync)
#tempFullSync <- na.omit(tempFullSync) ####### Change this so don't take out NA in stream temperature
dim(tempFullSync)


# Standardize for Analysis

tempFullSyncS <- cbind(tempFullSync[ ,c(1:9)],
                       apply(X = tempFullSync[ ,10:dim(tempFullSync)[2]], MARGIN=2,
                             FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))


tempFullSyncS[is.na(tempFullSyncS)] <- 0

fixEf <- modSummary@fixEf[ ,"Mean"]
names(fixEf) <- row.names(modSummary@fixEf)

tempFullSync$tempPredicted <- modSummary@fixEf["intercept", "Mean"] + 
  modSummary@fixEf["lat", "Mean"]*tempFullSyncS$Latitude + 
  modSummary@fixEf["lon", "Mean"]*tempFullSyncS$Longitude + 
  modSummary@fixEf["airTemp", "Mean"]*tempFullSyncS$airTemp + 
  modSummary@fixEf["airTempLag1", "Mean"]*tempFullSyncS$airTempLagged1 + 
  modSummary@fixEf["airTempLag2", "Mean"]*tempFullSyncS$airTempLagged2 + 
  modSummary@fixEf["precip", "Mean"]*tempFullSyncS$prcp + 
  modSummary@fixEf["precipLag1", "Mean"]*tempFullSyncS$prcpLagged1 + 
  modSummary@fixEf["precipLag2", "Mean"]*tempFullSyncS$prcpLagged2 + 
  modSummary@fixEf["drainage", "Mean"]*tempFullSyncS$TotDASqKM + 
  modSummary@fixEf["forest", "Mean"]*tempFullSyncS$Forest + 
  modSummary@fixEf["elevation", "Mean"]*tempFullSyncS$ReachElevationM + 
  modSummary@fixEf["coarseness", "Mean"]*tempFullSyncS$SurficialCoarseC + 
  modSummary@fixEf["wetland", "Mean"]*tempFullSyncS$CONUSWetland + 
  modSummary@fixEf["impoundments", "Mean"]*tempFullSyncS$ImpoundmentsAllSqKM + 
  modSummary@fixEf["swe", "Mean"]*tempFullSyncS$swe + 
  modSummary@fixEf["dOY", "Mean"]*tempFullSyncS$dOY + 
  modSummary@fixEf["dOY2", "Mean"]*tempFullSyncS$dOY^2 + 
  modSummary@fixEf["dOY3", "Mean"]*tempFullSyncS$dOY^3




# plot observed and predicte vs day of the year for all sites
sites <- unique(tempDataSync$site)

for(i in 1:length(unique(tempDataSync$site))){
  dataSite <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, temp)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + geom_point(aes(dOY, streamTempPred), colour = 'red', size=1) + geom_line(aes(dOY, streamTempPred), colour = 'red', size=0.1) + geom_point(aes(dOY, airTemp), colour='black', size=1) + ggtitle(unique(tempDataSync$fsite)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast


# plot observed and predicte vs day of the year for all sites
sites <- unique(tempFullSync$site)

for(i in 1:length(unique(tempFullSync$site))){
  dataSite <- filter(tempFullSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, tempPredicted)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + ggtitle(unique(tempFullSync$site)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast





