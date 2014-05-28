rm(list=ls())

library(DataCombine) # for the slide function
library(ggplot2)
library(relaimpo)
#library(lme4)
library(plyr)
library(reshape)
#library(ggmap)
library(foreign)
library(maptools)
#library(gridExtra)
library(nlme)

setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')

baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

WB <- T

if(WB) {
  load("dataOut/etWB.RData")
} else {
  load(paste0(dataOutDir, 'et.RData'))
  et$date <- et$date.x
}

et <- et[order(et$count),] # just to make sure et is ordered for the slide function

et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -1, NewVar='airTempLagged1')
et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -2, NewVar='airTempLagged2')

et <- slide(et, Var = "prcp", GroupVar = "site", slideBy = -1, NewVar='prcpLagged1')
et <- slide(et, Var = "prcp", GroupVar = "site", slideBy = -2, NewVar='prcpLagged2')

# Log scale swe before splitting into segments:
et$sweL <- log(et$swe + 0.001)

# Scale daymet vaiables in segment 2:
et$daylS <- (et$dayl - mean(et$dayl,na.rm=T)) / sd(et$dayl, na.rm=T)
et$sradS <- (et$srad - mean(et$srad,na.rm=T)) / sd(et$srad, na.rm=T)
et$sweLS <- (et$sweL - mean(et$sweL,na.rm=T)) / sd(et$sweL, na.rm=T)


# Function that returns Root Mean Squared Error
rmse <- function(error) {
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error) {
  mean(abs(error))
}


et$day <- as.numeric(et$date)

if(WB) {
  et2 <- et[ , c("siteYear", "year", "site", "date", "season", "temp", "airTemp", "airTempLagged1", "airTempLagged2", "dOY", "day", "prcp", "prcpLagged1", "prcpLagged2", "srad", "dayl", "swe")]
  
  et2$temp[which(et2$temp < 0 & et2$temp > -1)] <- 0
  et2$temp[which(et2$temp <= -1)] <- NA
  et2 <- na.omit(et2)
  
  # Scale if necessary or desired
  # Zuur p. 485
  # log.dams isn't standardized because not continuous but could be
  etS <- cbind(et2[ ,c(1:6)],
               apply(X = et2[ ,7:dim(et2)[2]], MARGIN=2,
                     FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))
  
} else {
  et2 <- et[ , c("siteYear", "year", "site", "HUC_4", "HUC_8", "HUC_12", "TNC_DamCount", "date", "dOY", "day", "temp", "airTemp", "prcp", "prcpLagged1", "prcpLagged2", "srad", "dayl", "swe", "airTempLagged1", "airTempLagged2", "Latitude", "Longitude", "Forest", "Agriculture", "Herbacious", "Developed", "Wetland", "Water", "Impervious", "BasinSlopeDEG", "HydrologicGroupAB", "SurficialCoarseC", "TotDASqKM", "BasinElevationM",  "CONUSWetland", "ImpoundmentsOpenSqKM", "ReachSlopePCNT")] # flow, "season",
  
  et2$temp[which(et2$temp < 0 & et2$temp > -1)] <- 0
  et2$temp[which(et2$temp <= -1)] <- NA
  et2 <- et2[which(et2$TotDASqKM <= 2000), ]
  et2 <- et2[which(et2$ImpoundmentsOpenSqKM <= 40), ]
  et2 <- na.omit(et2)
  
  # Scale if necessary or desired
  # Zuur p. 485
  # log.dams isn't standardized because so far from normal
  etS <- cbind(et2[ ,c(1:11)],
               apply(X = et2[ ,12:dim(et2)[2]], MARGIN=2,
                     FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))
}


et2$jDay <- as.numeric(as.Date(et2$date, format = "%m%d%Y")) + 2440588
etS$jDayS <- (et2$jDay - mean(et2$jDay, na.rm = T)) / (sd(et2$jDay, na.rm = T))
etS$jDayC <- (et2$jDay - mean(et2$jDay, na.rm = T))

etS$sinDOY <- sin(2*pi/360*etS$dOY)
etS$cosDOY <- cos(2*pi/360*etS$dOY)
etS$scDOY <- sin(2*pi/360*etS$dOY) + cos(2*pi/360*etS$dOY)

etS$month <- format(etS$date, '%m')

etS$season <- NA

## Need better way to create season
#for(i in 1:length(etS$month)){
#if(etS$month[i] == '12') { etS$season[i]  <- 1 }
#if(etS$month[i]  == '01') { etS$season[i] <- 1 }
#if(etS$month[i]  == '02') { etS$season[i] <- 1 }
#if(etS$month[i] == '03') { etS$season[i] <- 2 }
#if(etS$month[i] == '04') { etS$season[i] <- 2 }
#if(etS$month[i] == '05') { etS$season[i] <- 2 }
#if(etS$month[i] == '06') { etS$season[i] <- 3 }
#if(etS$month[i] == '07') { etS$season[i] <- 3 }
#if(etS$month[i] == '08') { etS$season[i] <- 3 }
#if(etS$month[i] == '09') { etS$season[i] <- 4 }
#if(etS$month[i] == '10') { etS$season[i] <- 4 }
#if(etS$month[i] == '11') { etS$season[i] <- 4 }
#}
etS$fSeason <- as.factor(etS$season)

etS$fyear <- as.factor(etS$year)
etS$fsite <- as.factor(etS$site)

#ggplot(data = et2, aes(date, temp)) + geom_point(size = 1) + facet_wrap(~site)

#-----------Check Data Relationships---------------

## put histograms on the diagonal
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  }

panel.cor <- function(x, y, method = "pearson", use = "pairwise.complete.obs", digits = 2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  x = x
  y = y
  Cor <- switch(method,
                spearman = cor(x, y, method = "spearman", use = use),
                pearson = cor(x, y, method = "pearson", use = use),
                kendell = cor(x, y, method = "kendall", use = use),
                stop("\nThe type of correlation must be spearman, pearson, or kendall\n"))
  txt <- format(c(Cor, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep="")
  #if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)
  #text(0.9, 0.9, Cor.Type, cex = 1)
}

panel.smooth.big <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) {
  x <- sample(x = x, size=1000, replace = F)
  y <- sample(x = y, size = 1000, replace = F)
    #x <- rnorm(500, mean(x, na.rm = T), sd(x, na.rm = T)) # alpha bleeding
    #y <- rnorm(500, mean(y, na.rm = T), sd(y, na.rm = T)) # alpha bleeding
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}

Pairs1 <- etS[c(runif(1000, 1, length(etS$temp))) , c("temp", "date.x", "dOY", "airTemp", "airTempLagged1", "airTempLagged2", "sinDOY", "cosDOY", "prcp", "srad", "dayl", "swe")]

pairs(Pairs1, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)

Pairs2 <- etS[c(runif(1000, 1, length(etS$temp))) , c("temp", "Forest", "Agriculture", "HydrologicGroupAB", "SurficialCoarseC", "TotDASqKM", "BasinElevationM",  "CONUSWetland", "ImpoundmentsOpenSqKM", "ReachSlopePCNT")]

pairs(Pairs2, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)

par(mfrow = c(1,1))

#-------------Analysis---------------------
# STEP 1: Try full, beyond optimal model to get a working autocorrelation

#system.time(lmeFull <- lme(temp ~ airTemp+airTempLagged1+airTempLagged2+
#                 Latitude+Longitude+
#                 Forest+ Agriculture+
#                 BasinElevationM+ ReachSlopePCNT+ 
#                 CONUSWetland+ SurficialCoarseC+
#                 dayl + srad + swe +
#                 sin(2*pi/360*dOY) + cos(2*pi/360*dOY), random = list(site = ~ airTemp + airTempLagged1 + airTempLagged2 + Forest, year = ~sin(2*pi/360*dOY) + cos(2*pi/360*dOY)), data = etS, na.action = "na.omit")) # |year/site ??? TotDASqKM+ ImpoundmentsOpenSqKM+ 



################## Cubic Day of the Year ##################
library(nlme)

# consider replacing lags with 5 and/or 10 day average airT and total precip. Also need drainage area
lme1 <- lme(temp ~ airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 +  airTemp*prcp + dOY + I(dOY^2) + I(dOY^3), random = list(year = ~1, site = ~airTemp ), data = etS1)
summary(lme1)

plot(lme1)
plot(etS1$date, etS1$temp, cex = 0.1, pch = 16, ylim = c(-5, 25)) 
I <- order(etS1$date)
points(etS1$date[I], fitted(lme1)[I], col = 'red', cex = 0.1)
abline(h = 0)

plot(fitted(lme1), resid(lme1))
abline(h = 0, col = 'red')

acf(resid(lme1), lag = 1000) 

plot(etS1$dOY, resid(lme1))

rmse(resid(lme1)) # 0.886


# without dOY
lme2 <- lme(temp ~ airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 +  airTemp*prcp, random = list(year = ~1, site = ~airTemp), data = etS1)
summary(lme2)

plot(lme2)
plot(etS1$date, etS1$temp, cex = 0.1, pch = 16, ylim = c(-5, 25)) 
I <- order(etS1$date)
points(etS1$date[I], fitted(lme2)[I], col = 'red', cex = 0.1)
abline(h = 0)

plot(fitted(lme2), resid(lme2))
abline(h = 0, col = 'red')

acf(resid(lme2), lag = 1000) # autocorrelation worse without dOY

plot(etS1$dOY, resid(lme2))

rmse(resid(lme2)) # 1.23


#library(gamm4) - better if don't account for correlation structure

# bam - for large data


library(mgcv)

?gamm
?mgcv
?smooth.terms
?adaptive.smooth
vignette('mgcv')
vignette('gamm')

# gamm4 is faster and has a better, more accuracte (no PQL) estimator but can't handle autocorrelation. Try first with gamm4. If full model has no autocorrelation of the residuals, continue model selection with gamm4. Otherwise switch to mgcv for the autocorrelation functions

# use bam instead because of size of data

# use bam within mgcv
library(mgcv)
system.time(bamWinter <- bam(temp ~ airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 +  airTemp*prcp + s(dOY, by = fsite) + s(fsite, bs = 're') + s(fyear, bs = 're'), data = etS))
summary(bamWinter) # r2 = 0.974
plot(bamWinter)
plot(etS$date, etS$temp, cex = 0.1, pch = 16, ylim = c(-5, 25)) 
I <- order(etS$date)
points(etS$date[I], fitted(bamWinter)[I], col = 'red', cex = 0.1)
abline(h = 0) # does seem to improve the winter estimates

rmse(resid(bamWinter)) # 1.005

# Analysis of Non-Winter Data (Alt: could use winter break points)
etS1 <- etS[which(etS$dOY > 83 & etS$dOY < 358), ]

# bam

system.time(bam1 <- bam(temp ~ airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 +  airTemp*prcp + s(dOY, by = fsite) + s(fsite, bs = 're') + s(fyear, bs = 're'), data = etS1))
summary(bam1) # r2 = 0.971
plot(bam1)
plot(etS1$date, etS1$temp, cex = 0.1, pch = 16, ylim = c(-5, 25)) 
I <- order(etS1$date)
points(etS1$date[I], fitted(bam1)[I], col = 'red', cex = 0.1)
abline(h = 0) # does seem to improve the winter estimates

plot(fitted(bam1), resid(bam1))
abline(h = 0, col = 'red')

plot(etS1$temp, resid(bam1))
abline(h = 0, col = 'red')

par(mfrow = c(1, 2))
hist(resid(bam1)); boxplot(resid(bam1))
par(mfrow = c(1,1))

acf(resid(bam1), lag = 1000) # moderate correlation problem

rmse(resid(bam1)) # 0.886
mae(resid(bam1)) # 0675

vis.gam(bam1,theta=30,ticktype="detailed", view=c("airTemp", "dOY"))
vis.gam(bam1,theta=-45,ticktype="detailed",se=2,  view=c("airTemp", "dOY"))
vis.gam(bam1,plot.type="contour",  view=c("airTemp", "dOY"))

vis.gam(bam1,theta=30,ticktype="detailed", view=c("prcp", "dOY"))
vis.gam(bam1,theta=-45,ticktype="detailed",se=2,  view=c("prcp", "dOY"))
vis.gam(bam1,plot.type="contour",  view=c("prcp", "dOY"))

vis.gam(bam1,theta=30,ticktype="detailed", view=c("airTemp", "site"))
vis.gam(bam1,theta=-45,ticktype="detailed",se=2,  view=c("airTemp", "site"))
vis.gam(bam1,plot.type="contour",  view=c("airTemp", "site"))



# gamm4
library(gamm4)
?gamm4

system.time(gamm4Full <- gamm4(temp ~ airTemp + airTempLagged1+airTempLagged2+ prcp + prcpLagged1 + prcpLagged2 +Latitude+Longitude+Forest+ Agriculture+BasinElevationM+ ReachSlopePCNT+CONUSWetland+ SurficialCoarseC + s(dOY) + prcp*airTemp, random = ~(1| site) + (1 | year), data = etS)) 

# Error in crossprod(root.phi %*% Zt) : 
# Cholmod error 'problem too large' at file ../Core/cholmod_aat.c, line 173
# In addition: Warning message:
#   In optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp),  :
#               convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded


system.time(gamm4Full <- gamm4(temp ~ airTemp + airTempLagged1+airTempLagged2+ prcp + prcpLagged1 + prcpLagged2 + s(dOY) + prcp*airTemp, random = ~(1| site), data = etS))


summary(gamm4Full)
plot(fitted(gamm4Full), resid(gamm4Full))

vis.gam(gamm4Full,theta=30,ticktype="detailed")
vis.gam(gamm4Full,theta=-45,ticktype="detailed",se=2)
vis.gam(gamm4Full,plot.type="contour")

#### mgcv

system.time(gammFull <- gamm(temp ~ airTemp + airTempLagged1+airTempLagged2+ prcp + prcpLagged1 + prcpLagged2 +Latitude+Longitude+Forest+ Agriculture+BasinElevationM+ ReachSlopePCNT+CONUSWetland+ SurficialCoarseC + s(dOY) + prcp*airTemp, random = list(site =~airTemp + airTempLagged1+airTempLagged2+ prcp + prcpLagged1 + prcpLagged2, year = ~ 1), data = etS)) # consider sum of 3 day precip * airTemp for the dampening effect of flow

vis.gam(ct1,theta=30,ticktype="detailed")
vis.gam(ct1,theta=-45,ticktype="detailed",se=2)
vis.gam(ct1,plot.type="contour")

# Test if can do random smoothers for each site
gammTest <- gamm(temp ~ airTemp + s(dOY, by = site), random = list(site = ~airTemp, year = ~1), data = etS)
summary(gammTest)

vis.gam(ct1,theta=30,ticktype="detailed")
vis.gam(ct1,theta=-45,ticktype="detailed",se=2)
vis.gam(ct1,plot.type="contour")

system.time(m1Gamm <- gamm(temp ~ s(jDayS), data = etS))
summary(m1Gamm$gam, cor = FALSE)
plot(m1Gamm$gam, seWithMean = TRUE)
summary(m1Gamm$lme)
acf(resid(m1Gamm$lme, type = "normalised"))

system.time(m2Gamm <- gamm(temp ~ s(jDayC), data = etS))
summary(m2Gamm$gam, cor = FALSE)
plot(m2Gamm$gam)
summary(m2Gamm$lme)

m2GammPred <- predict(m2Gamm$gam, se = FALSE, type = "response")
plot(etS$date.x, etS$temp, cex = 0.1, pch = 16) 
I <- order(etS$date.x)
lines(etS$date.x[I], m2GammPred[I], col = 'red')
points(etS$date.x[I], fitted(m2Gamm$gam)[I], col = 'red')

plot(fitted(m2Gamm$lme), residuals(m2Gamm$lme))
# plot(fitted(m2Gamm$gam), residuals(m2Gamm$gam)) # same as above thankfully


system.time(m2Gammb <- gamm(temp ~ airTemp + s(jDayC, k = 32), data = etS))
summary(m2Gammb$gam, cor = FALSE)
plot(m2Gammb$gam)
summary(m2Gammb$lme)

m2GammbPred <- predict(m2Gammb$gam, se = FALSE, type = "response")
plot(etS$date.x, etS$temp, cex = 0.1, pch = 16, ylim = c(-5, 30)) # adding many knots allows jDay to fit well which could alleviate need for year effect
I <- order(etS$date.x)
abline(h = 0)
points(etS$date.x[I], fitted(m2Gammb$gam)[I], col = 'red', cex=0.1) # problem is still predicts water temps well below 0
# lines(etS$date.x[I], m2GammbPred[I], col = 'blue', lwd = 0.2) 

# Add site effects and more knots
system.time(m2Gammc <- gamm(temp ~ airTemp + s(jDayC, k = 100), random = list(site = ~1), data = etS))
summary(m2Gammc$gam, cor = FALSE)
plot(m2Gammc$gam)
summary(m2Gammc$lme)


############
# GAM smoothers by site ???
############

system.time(m3Gamm <- gamm(temp ~  s(year) + s(dOY), random = list(site =~1), data = etS)) # 13 min
summary(m3Gamm$gam)
summary(m3Gamm$lme)
plot(m3Gamm$gam)

m3GammPred <- predict(m3Gamm$gam, se = FALSE, type = "response")
plot(etS$date.x, etS$temp, cex = 0.1, pch = 16) # some water temps well below 0
I <- order(etS$date.x)
#lines(etS$date.x[I], m3GammPred[I], col = 'red')
points(etS$date.x[I], fitted(m3Gamm$gam)[I], col = 'red')

plot(fitted(m3Gamm$lme), residuals(m3Gamm$lme))

rmse(resid(m3Gamm$lme))

system.time(m4Gamm <- gamm(temp ~ airTemp+airTempLagged1+airTempLagged2+ ForestS + s(year) + s(dOY), random = list(site =~1), data = etS))
summary(m4Gamm$gam)
summary(m4Gamm$lme)
plot(m4Gamm$gam)

system.time(m5Gamm <- gamm(temp ~ airTemp+airTempLagged1+airTempLagged2+ LatitudeS+LongitudeS+ ForestS+ BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ WetSlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ daylS + sradS + sweLS + s(year) + s(dOY), random = list(site =~1 + airTemp), data = etS)) # 26 min
summary(m5Gamm$gam)
summary(m5Gamm$lme)
plot(m5Gamm$gam)

system.time(m5GammIdent <- gamm(temp ~ airTemp+airTempLagged1+airTempLagged2+ LatitudeS+LongitudeS+ ForestS+ BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ WetSlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ daylS + sradS + sweLS + s(year) + s(dOY), random = list(site =~1 + airTemp), data = etS, weights = varIdent(form =~1 | site))) # 10 min
summary(m5GammIdent$gam)
summary(m5GammIdent$lme)
plot(m5GammIdent$gam)

system.time(m5GammPow <- gamm(temp ~ airTemp+airTempLagged1+airTempLagged2+ LatitudeS+LongitudeS+ ForestS+ BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ WetSlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ daylS + sradS + sweLS + s(year) + s(dOY), random = list(site =~1 + airTemp), data = etS, weights = varPower(form =~dOY))) # 37 min
summary(m5GammPow$gam)
summary(m5GammPow$lme)
plot(m5GammPow$gam)

system.time(m6Gamm <- gamm(temp ~ airTemp+airTempLagged1+airTempLagged2+ LatitudeS+LongitudeS+ ForestS+ BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ WetSlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ daylS + sradS + sweLS + s(year) + s(dOY), random = list(site =~1 + airTemp + s(dOY)), data = etS)) #
summary(m6Gamm$gam)
summary(m6Gamm$lme)
plot(m6Gamm$gam)


# try separate smoothers for HUC 4 or HUC 8 or based on latitudinal groupings

