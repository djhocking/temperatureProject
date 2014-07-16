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

stdCovs <- function(x, y, varNames){
  xStd <- as.data.frame(matrix(NA, dim(x)[1], length(varNames)))
  names(xStd) <- varNames
  for(i in 1:length(varNames)){
    xStd[ , varNames[i]] <- (x[ , varNames[i]] - mean(y[ , varNames[i]], na.rm=T)) / sd(y[ , varNames[i]], na.rm=T)
  }
  return(xStd)
}

load(paste0(dataOutDir, 'tempDataSync.RData'))
varNames1 <- names(tempFullSync[ ,8:dim(tempFullSync)[2]])

tempFullStd <- stdCovs(tempFullSync, tempDataSync, varNames1)
tempFullSyncS <- cbind(tempFullSync[ ,c(1:7)], tempFullStd)

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


system.time(tempFullSiteYearS$tempPredicted <- modSummary@fixEf["intercept", "Mean"] +
  BSite[tempFullSiteYearS$site, "intercept.site"] + 
  BYear[tempFullSiteYearS$cYear, "intercept.year"] + 
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
  (BYear[tempFullSiteYearS$cYear, "dOY"] + modSummary@fixEf["dOY", "Mean"])*tempFullSiteYearS$dOY + 
  (BYear[tempFullSiteYearS$cYear, "dOY2"] + modSummary@fixEf["dOY2", "Mean"])*((tempFullSiteYearS$dOY)^2) + 
  (BYear[tempFullSiteYearS$cYear, "dOY3"] + modSummary@fixEf["dOY3", "Mean"])*((tempFullSiteYearS$dOY)^3))


tempFullSiteS$tempPredicted <- modSummary@fixEf["intercept", "Mean"] +
  BSite[tempFullSiteS$site, "intercept.site"] + 
  modSummary@fixEf["lat", "Mean"]*tempFullSiteS$Latitude + 
  modSummary@fixEf["lon", "Mean"]*tempFullSiteS$Longitude + 
  (BSite[tempFullSiteS$site, "airTemp"] + modSummary@fixEf["airTemp", "Mean"])*tempFullSiteS$airTemp + 
  (BSite[tempFullSiteS$site, "airTempLag1"] + modSummary@fixEf["airTempLag1", "Mean"])*tempFullSiteS$airTempLagged1 + 
  (BSite[tempFullSiteS$site, "airTempLag2"] + modSummary@fixEf["airTempLag2", "Mean"])*tempFullSiteS$airTempLagged2 + 
  (BSite[tempFullSiteS$site, "precip"] + modSummary@fixEf["precip", "Mean"])*tempFullSiteS$prcp + 
  (BSite[tempFullSiteS$site, "precipLag1"] + modSummary@fixEf["precipLag1", "Mean"])*tempFullSiteS$prcpLagged1 + 
  (BSite[tempFullSiteS$site, "precipLag2"] + modSummary@fixEf["precipLag2", "Mean"])*tempFullSiteS$prcpLagged2 + 
  (BSite[tempFullSiteS$site, "drainage"] + modSummary@fixEf["drainage", "Mean"])*tempFullSiteS$TotDASqKM + 
  (BSite[tempFullSiteS$site, "forest"] + modSummary@fixEf["forest", "Mean"])*tempFullSiteS$Forest + 
  (BSite[tempFullSiteS$site, "elevation"] + modSummary@fixEf["elevation", "Mean"])*tempFullSiteS$ReachElevationM + 
  (BSite[tempFullSiteS$site, "coarseness"] + modSummary@fixEf["coarseness", "Mean"])*tempFullSiteS$SurficialCoarseC + 
  (BSite[tempFullSiteS$site, "wetland"] + modSummary@fixEf["wetland", "Mean"])*tempFullSiteS$CONUSWetland + 
  (BSite[tempFullSiteS$site, "impoundments"] + modSummary@fixEf["impoundments", "Mean"])*tempFullSiteS$ImpoundmentsAllSqKM + 
  (BSite[tempFullSiteS$site, "swe"] + modSummary@fixEf["swe", "Mean"])*tempFullSiteS$swe + 
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
    ggtitle(unique(tempFullSync$site)[i]) + 
    facet_wrap(~year) + 
    xlab(label = 'Day of the year') + ylab('Temperature (C)') + 
    theme(axis.text.x = element_text(angle = 45))
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/fullRecord/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=12,height=8, units='in' )
} # surprisingly fast but wouldn't do for all catchments


# plot observed and predicte vs day of the year for all sites
sites <- unique(tempDataSync$site)

for(i in 1:length(unique(tempDataSync$site))){
  dataSiteObs <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSiteObs, aes(dOY, temp)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + geom_point(aes(dOY, streamTempPred), colour = 'red', size=1) + geom_line(aes(dOY, streamTempPred), colour = 'red', size=0.1) + geom_point(aes(dOY, airTemp), colour='black', size=1) + ggtitle(unique(tempDataSync$fsite)[i]) + facet_wrap(~year) + xlab(label = 'Day of the year') + ylab('Temperature (C)')
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=6,height=4, units='in' )
} # surprisingly fast

rmse(tempFullSync[!is.na(tempDataSync$temp), "temp"] - tempFullSync[!is.na(tempDataSync$temp), "tempPredicted"])
############## Derived metrics ##########

# Mean maximum daily mean temperature by site (over years)
library(dplyr)
by.site <- group_by(tempFullSync, site)
by.site.year <- group_by(by.site, year, add = TRUE)
max.t <- filter(by.site, streamTemp == max(streamTemp))
summarise(max.t, mean(streamTemp)) # not needed - already max.t
summarise(by.site.year, sd(mean(streamTemp))) # not working based on filter or grouping

(max.t.site.year <- summarise(by.site.year, max(streamTemp)))
names(max.t.site.year) <- c("site", "year", "streamTemp")
max.t.site.year1 <- merge(as.data.frame(max.t.site.year), pred.t, by=c("site", "streamTemp"), all.x=T, all.y=F)

ggplot(pred.t, aes(dOY.real, temp)) + geom_point(size=1, colour='black') + geom_point(aes(dOY.real, streamTemp), colour = 'red', size=0.75) + ylab(label="Stream temperature (C)") + xlab("Day of the year") + geom_point(data=max.t.site.year1, aes(dOY.real, streamTemp), colour = "green") + facet_grid(site ~ year) # max temp points all replicated on every panel

# Number of days with stream temp > 20C
days.20 <- summarise(by.site.year, days.20 = length(streamTemp >= 20))
summarise(days.20, mean(days.20))

ggplot(pred.t[which(pred.t$site == "WB OBEAR" & pred.t$year == 2010), ], aes(dOY.real, streamTemp)) + 
  geom_point(size=2, colour = "black") + geom_line(colour = 'black') +
  geom_abline(intercept = 18, slope=0, colour='red') +
  geom_point(data = pred.t[which(pred.t$site == "WB OBEAR" & pred.t$year == 2010 & pred.t$streamTemp >= 18), ], aes(dOY.real, streamTemp), colour='red') +
  xlab("Day of the year") +
  ylab("Stream temperature (C)") #+ theme_classic()

# Resistance to peak air temperature
WB.2011.summer <- pred.t[which(pred.t$site == "WEST BROOK" & pred.t$year == 2011 & pred.t$dOY.real >=145 & pred.t$dOY.real <= 275), ]
sum(WB.2011.summer$airTemp - WB.2011.summer$streamTemp)

ggplot(pred.t[which(pred.t$site == "WEST BROOK" & pred.t$year == 2011), ], aes(dOY.real, streamTemp)) + 
  geom_point(size=2, colour = "black") + geom_line(colour = 'black') +
  geom_point(data=et2[which(et2$site == "WEST BROOK" & et2$year == 2011), ], aes(dOY, airTemp), colour = "red", size=2) + 
  geom_line(data=et2[which(et2$site == "WEST BROOK" & et2$year == 2011), ], aes(dOY, airTemp), colour = "red") + 
  geom_ribbon(data = pred.t[which(pred.t$site == "WEST BROOK" & pred.t$year == 2011 & pred.t$dOY.real >=145 & pred.t$dOY.real <= 275), ], aes(x=dOY.real, ymin=streamTemp, ymax=airTemp), fill="dark grey", alpha=.5) +
  xlab("Day of the year") +
  ylab("Temperature (C)") #+ theme_classic()

# Reset ggplot2 theme default to gray
theme_set(theme_gray())


