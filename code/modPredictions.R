rm(list=ls())

library(ggplot2)
library(dplyr)
library(DataCombine) # for the slide function

#setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
#baseDir <- 'C:/Users/dhocking/Documents/temperatureProject/'
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

summary(tempFullSyncS)
tempFullSyncS[is.na(tempFullSyncS)] <- 0

fixEf <- modSummary@fixEf[ ,"Mean"]
names(fixEf) <- row.names(modSummary@fixEf)

load(paste0(dataOutDir, 'tempDataSync.RData'))
#tempFullSync <- tempFullSync[which(tempFullSync$site %in% unique(tempDataSync$site)), ]
#tempFullSyncS <- tempFullSync[which(tempFullSyncS$site %in% unique(tempDataSync$site)), ]
sites <- unique(tempFullSync$site)
BSite <- modSummary@BSite
BYear <- modSummary@BYear


# Split data by site-year then do predictions for those with observed stream temperature data and those without, then recombine. The problem is that sites outside of the years observed won't get the site-specific values and years with data but at different sites won't get the site-specific data.
tempFullSyncS$siteYear <- paste0(tempFullSyncS$site, tempFullSyncS$year)
tempDataSyncS$siteYear <- paste0(tempDataSyncS$site, tempDataSyncS$year)

tempFullSiteYearS <- tempFullSyncS[which(tempFullSyncS$siteYear %in% unique(tempDataSyncS$siteYear)), ]
tempFullMeanS <- subset(tempFullSyncS, !(tempFullSyncS$siteYear %in% unique(tempDataSyncS$siteYear)))

# this will work fo MA because not predicting to any completely new sites
tempFullSiteYearS <- filter(tempFullSyncS, filter = year %in% unique(tempDataSyncS$year))
tempFullSiteS <- filter(tempFullSyncS, filter = !(year %in% unique(tempDataSyncS$year)))


system.time(tempFullSiteYearS$tempPredicted <- modSummary@fixEf["intercept", "Mean"] +
  BSite[tempFullSiteYearS$site, "intercept.site"] + # check
  BYear[tempFullSiteYearS$year, "intercept.year"] + # check
  modSummary@fixEf["lat", "Mean"]*tempFullSiteYearS$Latitude + 
  modSummary@fixEf["lon", "Mean"]*tempFullSiteYearS$Longitude + 
  (BSite[tempFullSiteYearS$site, "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSiteYearS$airTemp + 
  (BSite[tempFullSiteYearS$site, "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSiteYearS$airTempLagged1 + 
  (BSite[tempFullSiteYearS$site, "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSiteYearS$airTempLagged2 + 
  (BSite[tempFullSiteYearS$site, "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSiteYearS$prcp + 
  (BSite[tempFullSiteYearS$site, "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSiteYearS$prcpLagged1 + 
  (BSite[tempFullSiteYearS$site, "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSiteYearS$prcpLagged2 + 
  (BSite[tempFullSiteYearS$site, "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSiteYearS$TotDASqKM + 
  (BSite[tempFullSiteYearS$site, "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSiteYearS$Forest + 
  (BSite[tempFullSiteYearS$site, "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSiteYearS$ReachElevationM + 
  (BSite[tempFullSiteYearS$site, "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSiteYearS$SurficialCoarseC + 
  (BSite[tempFullSiteYearS$site, "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSiteYearS$CONUSWetland + 
  (BSite[tempFullSiteYearS$site, "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSiteYearS$ImpoundmentsAllSqKM + 
  (BSite[tempFullSiteYearS$site, "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSiteYearS$swe + 
  (BYear[tempFullSiteYearS$year, "dOY"] + modSummary@fixEf["dOY", "Mean"])*tempFullSiteYearS$dOY + 
  (BYear[tempFullSiteYearS$year, "dOY2"] + modSummary@fixEf["dOY2", "Mean"])*tempFullSiteYearS$dOY^2 + 
  (BYear[tempFullSiteYearS$year, "dOY3"] + modSummary@fixEf["dOY3", "Mean"])*tempFullSiteYearS$dOY^3)






tempFullSiteYearS$tempPredicted <- NA
pb <- txtProgressBar(min = 0, max = dim(tempFullSyncS)[1], style = 3)
for(i in 1:dim(tempFullSyncS)[1]){
  tempFullSiteYearS$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] +
      BSite[tempFullSiteYearS$site[i], "intercept.site"] + # check
      BYear[tempFullSiteYearS$year[i], "intercept.year"] + # check
      modSummary@fixEf["lat", "Mean"]*tempFullSiteYearS$Latitude[i] + 
      modSummary@fixEf["lon", "Mean"]*tempFullSiteYearS$Longitude[i] + 
      (BSite[tempFullSiteYearS$site[i], "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSiteYearS$airTemp[i] + 
      (BSite[tempFullSiteYearS$site[i], "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSiteYearS$airTempLagged1[i] + 
      (BSite[tempFullSiteYearS$site[i], "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSiteYearS$airTempLagged2[i] + 
      (BSite[tempFullSiteYearS$site[i], "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSiteYearS$prcp[i] + 
      (BSite[tempFullSiteYearS$site[i], "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSiteYearS$prcpLagged1[i] + 
      (BSite[tempFullSiteYearS$site[i], "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSiteYearS$prcpLagged2[i] + 
      (BSite[tempFullSiteYearS$site[i], "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSiteYearS$TotDASqKM[i] + 
      (BSite[tempFullSiteYearS$site[i], "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSiteYearS$Forest[i] + 
      (BSite[tempFullSiteYearS$site[i], "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSiteYearS$ReachElevationM[i] + 
      (BSite[tempFullSiteYearS$site[i], "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSiteYearS$SurficialCoarseC[i] + 
      (BSite[tempFullSiteYearS$site[i], "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSiteYearS$CONUSWetland[i] + 
      (BSite[tempFullSiteYearS$site[i], "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSiteYearS$ImpoundmentsAllSqKM[i] + 
      (BSite[tempFullSiteYearS$site[i], "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSiteYearS$swe[i] + 
      (BYear[tempFullSiteYearS$year[i], "dOY"] + modSummary@fixEf["dOY", "Mean"])*tempFullSiteYearS$dOY[i] + 
      (BYear[tempFullSiteYearS$year[i], "dOY2"] + modSummary@fixEf["dOY2", "Mean"])*tempFullSiteYearS$dOY[i]^2 + 
      (BYear[tempFullSiteYearS$year[i], "dOY3"] + modSummary@fixEf["dOY3", "Mean"])*tempFullSiteYearS$dOY[i]^3
  
  setTxtProgressBar(pb, i)
}
close(pb)


tempFullSiteS$tempPredicted <- NA
pb <- txtProgressBar(min = 0, max = dim(tempFullSyncS)[1], style = 3)
for(i in 1:dim(tempFullSyncS)[1]){
  tempFullSiteS$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] +
    BSite[tempFullSiteS$site[i], "intercept.site"] + # check
    BYear[tempFullSiteS$year[i], "intercept.year"] + # check
    modSummary@fixEf["lat", "Mean"]*tempFullSiteS$Latitude[i] + 
    modSummary@fixEf["lon", "Mean"]*tempFullSiteS$Longitude[i] + 
    (BSite[tempFullSiteS$site[i], "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSiteS$airTemp[i] + 
    (BSite[tempFullSiteS$site[i], "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSiteS$airTempLagged1[i] + 
    (BSite[tempFullSiteS$site[i], "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSiteS$airTempLagged2[i] + 
    (BSite[tempFullSiteS$site[i], "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSiteS$prcp[i] + 
    (BSite[tempFullSiteS$site[i], "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSiteS$prcpLagged1[i] + 
    (BSite[tempFullSiteS$site[i], "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSiteS$prcpLagged2[i] + 
    (BSite[tempFullSiteS$site[i], "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSiteS$TotDASqKM[i] + 
    (BSite[tempFullSiteS$site[i], "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSiteS$Forest[i] + 
    (BSite[tempFullSiteS$site[i], "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSiteS$ReachElevationM[i] + 
    (BSite[tempFullSiteS$site[i], "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSiteS$SurficialCoarseC[i] + 
    (BSite[tempFullSiteS$site[i], "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSiteS$CONUSWetland[i] + 
    (BSite[tempFullSiteS$site[i], "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSiteS$ImpoundmentsAllSqKM[i] + 
    (BSite[tempFullSiteS$site[i], "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSiteS$swe[i] + 
    modSummary@fixEf["dOY", "Mean"]*tempFullSiteS$dOY[i] + 
    modSummary@fixEf["dOY2", "Mean"]*tempFullSiteS$dOY[i]^2 + 
    modSummary@fixEf["dOY3", "Mean"]*tempFullSiteS$dOY[i]^3
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Make predictions
# this will be incredibly slow but I'm not sure of a better way to do this and incorporate the random effects unless it's done in JAGS but then it will have to be done for every iteration which will be time and memory intensive
tempFullSync$tempPredicted <- NA

pb <- txtProgressBar(min = 0, max = dim(tempFullSyncS)[1], style = 3)
for(i in 1:dim(tempFullSyncS)[1]){
  if(tempFullSyncS$site[i] %in% unique(tempDataSyncS$site) & tempFullSyncS$year[i] %in% unique(tempDataSyncS$year)){ # if from a site and year with data: use random site and year effects
      tempFullSync$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] +
        BSite[tempFullSyncS$site[i], "intercept.site"] + # check
        BYear[tempFullSyncS$year[i], "intercept.year"] + # check
        modSummary@fixEf["lat", "Mean"]*tempFullSyncS$Latitude[i] + 
        modSummary@fixEf["lon", "Mean"]*tempFullSyncS$Longitude[i] + 
        (BSite[tempFullSyncS$site[i], "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSyncS$airTemp[i] + 
        (BSite[tempFullSyncS$site[i], "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSyncS$airTempLagged1[i] + 
        (BSite[tempFullSyncS$site[i], "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSyncS$airTempLagged2[i] + 
        (BSite[tempFullSyncS$site[i], "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSyncS$prcp[i] + 
        (BSite[tempFullSyncS$site[i], "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSyncS$prcpLagged1[i] + 
        (BSite[tempFullSyncS$site[i], "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSyncS$prcpLagged2[i] + 
        (BSite[tempFullSyncS$site[i], "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSyncS$TotDASqKM[i] + 
        (BSite[tempFullSyncS$site[i], "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSyncS$Forest[i] + 
        (BSite[tempFullSyncS$site[i], "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSyncS$ReachElevationM[i] + 
        (BSite[tempFullSyncS$site[i], "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSyncS$SurficialCoarseC[i] + 
        (BSite[tempFullSyncS$site[i], "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSyncS$CONUSWetland[i] + 
        (BSite[tempFullSyncS$site[i], "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSyncS$ImpoundmentsAllSqKM[i] + 
        (BSite[tempFullSyncS$site[i], "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSyncS$swe[i] + 
        (BYear[tempFullSyncS$year[i], "dOY"] + modSummary@fixEf["dOY", "Mean"])*tempFullSyncS$dOY[i] + 
        (BYear[tempFullSyncS$year[i], "dOY2"] + modSummary@fixEf["dOY2", "Mean"])*tempFullSyncS$dOY[i]^2 + 
        (BYear[tempFullSyncS$year[i], "dOY3"] + modSummary@fixEf["dOY3", "Mean"])*tempFullSyncS$dOY[i]^3
    } else {
        if(tempFullSyncS$site[i] %in% unique(tempDataSyncS$site)){# if from a site with data but not from a year with data just use the random site intercept and slopes and mean year effects
      tempFullSync$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] +
        BSite[tempFullSyncS$site[i], "intercept.site"] + # check
        modSummary@fixEf["lat", "Mean"]*tempFullSyncS$Latitude[i] + 
        modSummary@fixEf["lon", "Mean"]*tempFullSyncS$Longitude[i] + 
        (BSite[tempFullSyncS$site[i], "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSyncS$airTemp[i] + 
        (BSite[tempFullSyncS$site[i], "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSyncS$airTempLagged1[i] + 
        (BSite[tempFullSyncS$site[i], "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSyncS$airTempLagged2[i] + 
        (BSite[tempFullSyncS$site[i], "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSyncS$prcp[i] + 
        (BSite[tempFullSyncS$site[i], "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSyncS$prcpLagged1[i] + 
        (BSite[tempFullSyncS$site[i], "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSyncS$prcpLagged2[i] + 
        (BSite[tempFullSyncS$site[i], "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSyncS$TotDASqKM[i] + 
        (BSite[tempFullSyncS$site[i], "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSyncS$Forest[i] + 
        (BSite[tempFullSyncS$site[i], "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSyncS$ReachElevationM[i] + 
        (BSite[tempFullSyncS$site[i], "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSyncS$SurficialCoarseC[i] + 
        (BSite[tempFullSyncS$site[i], "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSyncS$CONUSWetland[i] + 
        (BSite[tempFullSyncS$site[i], "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSyncS$ImpoundmentsAllSqKM[i] + 
        (BSite[tempFullSyncS$site[i], "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSyncS$swe[i] + 
        modSummary@fixEf["dOY", "Mean"]*tempFullSyncS$dOY[i] + 
        modSummary@fixEf["dOY2", "Mean"]*tempFullSyncS$dOY[i]^2 + 
        modSummary@fixEf["dOY3", "Mean"]*tempFullSyncS$dOY[i]^3
    } else {
      if(tempFullSyncS$year[i] %in% unique(tempDataSyncS$year)){ # if not from a site with data but from a year with data use the random year effects only
        tempFullSync$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] + 
          BSite[tempFullSyncS$site[i], "intercept.site"] + # check
          modSummary@fixEf["lat", "Mean"]*tempFullSyncS$Latitude[i] + 
          modSummary@fixEf["lon", "Mean"]*tempFullSyncS$Longitude[i] + 
          modSummary@fixEf["airTemp", "Mean"]*tempFullSyncS$airTemp[i] + 
          modSummary@fixEf["airTempLag1", "Mean"]*tempFullSyncS$airTempLagged1[i] + 
          modSummary@fixEf["airTempLag2", "Mean"]*tempFullSyncS$airTempLagged2[i] + 
          modSummary@fixEf["precip", "Mean"]*tempFullSyncS$prcp[i] + 
          modSummary@fixEf["precipLag1", "Mean"]*tempFullSyncS$prcpLagged1[i] + 
          modSummary@fixEf["precipLag2", "Mean"]*tempFullSyncS$prcpLagged2[i] + 
          modSummary@fixEf["drainage", "Mean"]*tempFullSyncS$TotDASqKM[i] + 
          modSummary@fixEf["forest", "Mean"]*tempFullSyncS$Forest[i] + 
          modSummary@fixEf["elevation", "Mean"]*tempFullSyncS$ReachElevationM[i] + 
          modSummary@fixEf["coarseness", "Mean"]*tempFullSyncS$SurficialCoarseC[i] + 
          modSummary@fixEf["wetland", "Mean"]*tempFullSyncS$CONUSWetland[i] + 
          modSummary@fixEf["impoundments", "Mean"]*tempFullSyncS$ImpoundmentsAllSqKM[i] + 
          modSummary@fixEf["swe", "Mean"]*tempFullSyncS$swe[i] + 
          (BYear[tempFullSyncS$year[i], "dOY"] + modSummary@fixEf["dOY", "Mean"])*tempFullSyncS$dOY[i] + 
          (BYear[tempFullSyncS$year[i], "dOY2"] + modSummary@fixEf["dOY2", "Mean"])*tempFullSyncS$dOY[i]^2 + 
          (BYear[tempFullSyncS$year[i], "dOY3"] + modSummary@fixEf["dOY3", "Mean"])*tempFullSyncS$dOY[i]^3
    } else { # if not from a site or year with observed data just predict using the mean fixed effects
      tempFullSync$tempPredicted[i] <- modSummary@fixEf["intercept", "Mean"] + 
        modSummary@fixEf["lat", "Mean"]*tempFullSyncS$Latitude[i] + 
        modSummary@fixEf["lon", "Mean"]*tempFullSyncS$Longitude[i] + 
        modSummary@fixEf["airTemp", "Mean"]*tempFullSyncS$airTemp[i] + 
        modSummary@fixEf["airTempLag1", "Mean"]*tempFullSyncS$airTempLagged1[i] + 
        modSummary@fixEf["airTempLag2", "Mean"]*tempFullSyncS$airTempLagged2[i] + 
        modSummary@fixEf["precip", "Mean"]*tempFullSyncS$prcp[i] + 
        modSummary@fixEf["precipLag1", "Mean"]*tempFullSyncS$prcpLagged1[i] + 
        modSummary@fixEf["precipLag2", "Mean"]*tempFullSyncS$prcpLagged2[i] + 
        modSummary@fixEf["drainage", "Mean"]*tempFullSyncS$TotDASqKM[i] + 
        modSummary@fixEf["forest", "Mean"]*tempFullSyncS$Forest[i] + 
        modSummary@fixEf["elevation", "Mean"]*tempFullSyncS$ReachElevationM[i] + 
        modSummary@fixEf["coarseness", "Mean"]*tempFullSyncS$SurficialCoarseC[i] + 
        modSummary@fixEf["wetland", "Mean"]*tempFullSyncS$CONUSWetland[i] + 
        modSummary@fixEf["impoundments", "Mean"]*tempFullSyncS$ImpoundmentsAllSqKM[i] + 
        modSummary@fixEf["swe", "Mean"]*tempFullSyncS$swe[i] + 
        modSummary@fixEf["dOY", "Mean"]*tempFullSyncS$dOY[i] + 
        modSummary@fixEf["dOY2", "Mean"]*tempFullSyncS$dOY[i]^2 + 
        modSummary@fixEf["dOY3", "Mean"]*tempFullSyncS$dOY[i]^3
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)
}

# plot observed and predicte vs day of the year for all sites
sites <- unique(as.character(tempFullSync$site))

for(i in 1:length(unique(tempFullSync$site))){
  dataSite <- filter(tempFullSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, tempPredicted)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + ggtitle(unique(tempFullSync$site)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/fullRecord/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast


# plot observed and predicte vs day of the year for all sites
sites <- unique(tempDataSync$site)

for(i in 1:length(unique(tempDataSync$site))){
  dataSite <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, temp)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + geom_point(aes(dOY, streamTempPred), colour = 'red', size=1) + geom_line(aes(dOY, streamTempPred), colour = 'red', size=0.1) + geom_point(aes(dOY, airTemp), colour='black', size=1) + ggtitle(unique(tempDataSync$fsite)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast



