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
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
#baseDir <- 'D:/projects/temperatureProject/'

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

load(paste0(dataOutDir, 'etWB.RData'))

et$day <- as.numeric(et$date)
et <- et[order(et$count),] # just to make sure et is ordered for the slide function
# airTemp
et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -1, NewVar='airTempLagged1')
et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -2, NewVar='airTempLagged2')

# replace missing spring and fall BP with means for clipping the data to the synchronized season
et1 <- et
et1[is.na(et$springBP), "springBP"] <- mean(et$springBP, na.rm=T)
et1[is.na(et$fallBP), "fallBP"] <- mean(et$fallBP, na.rm=T)
et1 <- et[which(et1$dOY >= et1$springBP & et1$dOY <= et1$fallBP), ]

# Make dataframe with just variables for modeling and order before standardizing
et2 <- et1[ , c("year", "site", "date", "season", "temp", "airTemp", "airTempLagged1", "airTempLagged2", "dOY", "day", "flow", "srad", "dayl", "swe")]

# Scale if necessary or desired
# Zuur p. 485
# log.dams isn't standardized because not continuous but could be
etS <- cbind(et2[ ,c(1:5)],
             apply(X = et2[ ,6:dim(et2)[2]], MARGIN=2,
                   FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))

etS$fyear <- as.factor(etS$year)
etS$fsite <- as.factor(etS$site)

etS <- na.omit(etS)
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

Pairs1 <- etS[c(runif(1000, 1, length(etS$temp))) , c("temp", "date", "dOY", "airTemp", "airTempLagged1", "airTempLagged2", "flow", "srad", "dayl", "swe")]

pairs(Pairs1, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)
par(mfrow = c(1,1))

#-------------Analysis---------------------
# STEP 1: Try full, beyond optimal model to get a working autocorrelatio

################## Cubic Day of the Year ##################

###### Independent and Crossed Random Effects of year and site########
library(lme4)

# Independent Random Intercepts
lmer1 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (1|year) + (1|site), data = etS)
summary(lmer1)
ranef(lmer1)

# Crossed Random Intercepts??????????
lmer2 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (1|year/site), data = etS)
summary(lmer2)
ranef(lmer2)

# Add random slopes
lmer3 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2|year/site), data = etS)
summary(lmer3)
ranef(lmer3)

# Add random slopes for dOY
lmer4 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 + dOY + I(dOY^2) + I(dOY^3)|year/site), data = etS) # convergence warnings

# Slopes for nested ...
lmer5 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 |year:site) + (dOY + I(dOY^2) + I(dOY^3)|year), data = etS)
summary(lmer5) # works and might make the most sense

# Random slopes for site within year only
lmer6 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 + dOY + I(dOY^2) + I(dOY^3) |year:site) + (1|year), data = etS)
summary(lmer6)
plot(lmer6)
plot(acf(resid(lmer6))) # not sure if this works with sites and years

# Random slopes for crossed
lmer7 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 |site) + (dOY + I(dOY^2) + I(dOY^3)|year), data = etS)
summary(lmer7) # works and potentially easier to define resilient sites


############# Nested Random Effects #################
library(nlme)

# consider replacing lags with 5 and/or 10 day average airT and total precip. Also need drainage area
lme2 <- lme(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3), random = list(year = ~1, site = ~1 ), data = etS)
summary(lme2)

# compare to lme4
cbind(fixef(lmer2), fixef(lme1)) # virtually identical
cbind(ranef(lmer2)$year, ranef(lme1)[1]) # essentially identical
data.frame(siteYearLMER=row.names(ranef(lmer2)$site), siteYearLME=row.names(ranef(lme1)$site), ranef(lmer2)$site, ranef(lme1)[2]) # different order but 

plot(lme1)
plot(etS1$date, etS1$temp, cex = 0.1, pch = 16, ylim = c(-5, 25)) 
I <- order(etS1$date)
points(etS1$date[I], fitted(lme1)[I], col = 'red', cex = 0.1)
abline(h = 0)

plot(fitted(lme1), resid(lme1))
abline(h = 0, col = 'red')

acf(resid(lme1), lag = 1000) 

plot(etS1$dOY, resid(lme1))

rmse(resid(lme1)) # 0.943

ranef(lme1)

plot(etS1$dOY, fitted(lme1))

etS1$lme1 <- fitted(lme1)
ggplot(data=etS1[etS1$year >= 2007, ], aes(dOY, lme1)) + geom_point(aes(dOY, temp)) + geom_line(colour='red') + facet_grid(year ~ site) 

# Add random slopes
system.time(lme3 <- lme(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3), random = list(year = ~1, site = ~airTemp + airTempLagged1 + airTempLagged2 + dOY + I(dOY^2) + I(dOY^3)), data = etS, control=c(maxIter=1000))) # no convergence

# Random slopes by year
system.time(lme4 <- lme(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3), random = list(year = ~ dOY + I(dOY^2) + I(dOY^3), site = ~airTemp + airTempLagged1 + airTempLagged2), data = etS, control=c(maxIter=1000))) # no convergence

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

etS1$lme1 <- fitted(lme1)
ggplot(data=etS1[etS1$year >= 2007, ], aes(dOY, lme1)) + geom_point(aes(dOY, temp)) + geom_line(colour='red') + facet_grid(year ~ site) 

plot(bam1)


system.time(bam1 <- bam(temp ~ airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 +  airTemp*prcp + s(dOY, by = fsite:fyear) + s(fsite, bs = 're') + s(fyear, bs = 're'), data = etS1))
summary(bam1) # r2 = 0.971


# gamm4
library(gamm4)
?gamm4

system.time(gamm4Full <- gamm4(temp ~ airTemp + airTempLagged1+airTempLagged2+ prcp + prcpLagged1 + prcpLagged2 +Latitude+Longitude+Forest+ Agriculture+BasinElevationM+ ReachSlopePCNT+CONUSWetland+ SurficialCoarseC + s(dOY) + prcp*airTemp, random = ~(1| site) + (1 | year), data = etS1)) 

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



################# JAGS #####################

plot(etS$airTemp, etS$temp)
points(etS$airTemp[10111], etS$temp[10111], col='red')

sink("code/cubicDay.txt")
cat("
    model{
    # Priors
    alpha ~ dnorm(0, 0.01)
    b.air ~ dnorm(0, 0.01)
    b.air1 ~ dnorm(0, 0.01)
    b.air2 ~ dnorm(0, 0.01)
    b.prcp ~ dnorm(0, 0.01)
    b.prcp1 ~ dnorm(0, 0.01)
    b.prcp2 ~ dnorm(0, 0.01)
    b.lat ~ dnorm(0, 0.01)
    b.lon ~ dnorm(0, 0.01)
    b.forest ~ dnorm(0, 0.01)
    b.elev ~ dnorm(0, 0.01)
    b.slope ~ dnorm(0, 0.01)
    b.wetland ~ dnorm(0, 0.01)
    b.coarse ~ dnorm(0, 0.01)
    b.drain ~ dnorm(0, 0.01)
    b.air.drain ~ dnorm(0, 0.01)
    sigma ~ dunif(0.00001, 0.99999)
    tau <- 1/(sigma*sigma)
    
    # Hyperpriors for Random Effects
    sigma.site ~ dunif(0, 5)
    tau.site <- 1 / sigma.site * sigma.site
    #tau.site ~ dgamma(0.01, 0.01)    
    #sigma.site <- pow(1/tau.site, 0.5) # sqrt of 1/tau
    
    sigma.year ~ dunif(0, 5)
    tau.year <- 1 / sigma.year * sigma.year
    
    for(i in 1:n.sites) {
    b.site[i] ~ dnorm(0, tau.site)
    }
    
    for(i in 1:n.years) {
    b.year[i] ~ dnorm(0, tau.year)
    }
    
    # Autocorrelation structure
    
    # Linear Model
    for(i in 1:n.obs) {
    stream.mu[i] <- alpha + b.site[site[i]] + b.year[year[i]] + b.air*airTemp[i] + b.air1*airTempLagged1[i] + b.air2*airTempLagged2[i] + b.prcp*prcp[i] + b.prcp1*prcpLagged1[i] + b.prcp2*prcpLagged2[i] + b.lat*Latitude[i] + b.lon*Longitude[i] + b.forest*forest[i] + b.elev*BasinElevationM[i] + b.slope*ReachSlopePCNT[i] + b.wetland*CONUSWetland[i] + b.coarse*SurficialCoarseC[i] + b.drain*Drainage[i] + b.air.drain*airTemp[i]*Drainage[i]
    temp[i] ~ dnorm(stream.mu[i], tau)T(0, 80) # prevent stream temperatures below 0
    }
    }
    ", fill = TRUE)
sink()

n.obs <- dim(etS)[1]
n.sites <- length(unique(etS$site))

inits1 <- function() {
  list(b.air = rnorm(1, 3.5, 0.5))
}

params1 <- c(    "alpha",
                 "b.air",
                 "b.air1",
                 "b.air2",
                 "b.prcp",
                 "b.prcp1",
                 "b.prcp2",
                 "b.lat",
                 "b.lon",
                 "b.forest",
                 "b.elev",
                 "b.slope",
                 "b.wetland",
                 "b.coarse",
                 "b.drain",
                 "b.air.drain",
                 "sigma.site",
                 "sigma.year",
                 "sigma")

data1 <- list(temp = etS$temp,
              airTemp = etS$airTemp,
              airTempLagged1 = etS$airTempLagged1,
              airTempLagged2 = etS$airTempLagged2,
              prcp = etS$prcp,
              prcpLagged1 = etS$prcpLagged1,
              prcpLagged2 = etS$prcpLagged2,
              Latitude = etS$Latitude,
              Longitude = etS$Longitude,
              site = as.factor(etS$site),
              year = as.factor(etS$year),
              forest = etS$Forest,
              BasinElevationM = etS$BasinElevationM,
              ReachSlopePCNT = etS$ReachSlopePCNT,
              CONUSWetland = etS$CONUSWetland,
              SurficialCoarseC = etS$SurficialCoarseC,
              Drainage = etS$TotDASqKM,
              n.obs = n.obs,
              n.years = length(unique(etS$year)),
              n.sites = n.sites)



n.burn = 100
n.it = 3000
n.thin = 6

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data1", "inits1", "params1", "n.obs", "n.sites", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/cubicDay.txt", data1, inits1, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params1, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
})

out32 <- mcmc.list(out)
stopCluster(CL)

rm(out)

#pdf(file="C:/Users/dhocking/Dropbox/out32.pdf", width=10, height=10)
plot(out32[ , c("alpha",
                "b.air",
                "b.air1",
                "b.air2",
                "b.prcp",
                "b.prcp1",
                "b.prcp2",
                "b.lat",
                "b.lon",
                "b.forest",
                "b.elev",
                "b.slope",
                "b.wetland",
                "b.coarse",
                "b.drain",
                "b.air.drain",
                "sigma.site",
                "sigma.year",
                "sigma")])
dev.off()
summary(out3[ , c("alpha",
                  "b.air",
                  "b.air1",
                  "b.air2",
                  "b.prcp",
                  "b.prcp1",
                  "b.prcp2",
                  "b.lat",
                  "b.lon",
                  "b.forest",
                  "b.elev",
                  "b.slope",
                  "b.wetland",
                  "b.coarse",
                  "b.drain",
                  "b.air.drain",
                  "sigma.site",
                  "sigma.year",
                  "sigma")])

library(ggmcmc)
#ggmcmc(ggs(out3[ , c("B.air", "B.air1", "B.air2", "gam", "alpha", "beta", "sigma")]))

ss <- summary(out3[ , -c(1:5)])
streamTmeans <- ss$statistics[ , "Mean"]

Coefs <- as.list(summary(out3[ , c(1:6)])$statistics[ , "Mean"])
names(Coefs)

si <- Coefs$mu + (Coefs$alpha - Coefs$mu) / (1 + exp(Coefs$b0.g*(Coefs$beta - etS$temp)))
#si2 <- sapply(si, FUN = rnorm, sd=Coefs$sigma, simplify = TRUE)
lsi.mu <- log(si)
err <- NA
for(i in 1:length(si)) { err[i] <- rnorm(1, lsi.mu[i], Coefs$sigma)}
si2 <- exp(lsi)

plot(etS$temp, si2)
abline(0, 1, col = 'red')
#m1 <- apply(out3[ , "stream.mu"])
Resids <- etS$temp - si
rmse(Resids)

plot(etS$temp, streamTmeans)
abline(0, 1, col = 'red')
#m1 <- apply(out3[ , "stream.mu"])
Resids <- etS$temp - streamTmeans
rmse(Resids)




Tstream <- matrix(NA, length(etS$temp), 1)
for(i in 1:length(etS$temp)){
  Tstream[i] <- apply(as.matrix(out3[ , c(paste("stream.mu[", i, "]", sep =""))]), 2, FUN = quantile, probs = c(0.5))
}

plot(etS$temp, Tstream)
abline(0, 1, col = 'red')

rmse((etS$temp - Tstream))

txt1 <- "done"
write.table(txt1, file = "C:/Users/dhocking/Dropbox/done.txt")

