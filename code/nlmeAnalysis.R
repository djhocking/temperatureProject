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

#load("dataOut/etWB.RData")
load(paste0(dataInDir, 'nlme.RData'))

et <- et[order(et$count),] # just to make sure et is ordered for the slide function

et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -1, NewVar='airTempLagged1')
et <- slide(et, Var = "airTemp", GroupVar = "site", slideBy = -2, NewVar='airTempLagged2')

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

et2 <- et[ , c("siteYear", "year", "site", "HUC_4", "HUC_8", "HUC_12", "TNC_DamCount", "date.x", "dOY", "temp", "airTemp", "prcp", "srad", "dayl", "swe", "airTempLagged1", "airTempLagged2", "Latitude", "Longitude", "Forest", "Agriculture", "Herbacious", "Developed", "Wetland", "Water", "Impervious", "BasinSlopeDEG", "HydrologicGroupAB", "SurficialCoarseC", "TotDASqKM", "BasinElevationM",  "CONUSWetland", "ImpoundmentsOpenSqKM", "ReachSlopePCNT")] # flow, "season",
et2 <- na.omit(et2)

et2$day <- as.numeric(et2$date)

# Scale if necessary or desired
# Zuur p. 485
# log.dams isn't standardized because so far from normal
etS <- cbind(et2[ ,c(1:9)],
                      apply(X = et2[ ,10:dim(et2)[2]], MARGIN=2,
                            FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))

#ggplot(data = et2, aes(date, temp)) + geom_point(size = 1) + facet_wrap(~site)

#-----------Check Correlations---------------

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

etS$sinDOY <- sin(2*pi/360*etS$dOY)
etS$cosDOY <- cos(2*pi/360*etS$dOY)
etS$scDOY <- sin(2*pi/360*etS$dOY) + cos(2*pi/360*etS$dOY)

Pairs1 <- etS[c(runif(1000, 1, length(etS$temp))) , c("temp", "date.x", "dOY", "airTemp", "prcp", "srad", "dayl", "swe", "airTempLagged1", "airTempLagged2", "scDOY")]

pairs(Pairs1, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)

Pairs2 <- etS[c(runif(1000, 1, length(etS$temp))) , c("temp", "Forest", "Agriculture", "HydrologicGroupAB", "SurficialCoarseC", "TotDASqKM", "BasinElevationM",  "CONUSWetland", "ImpoundmentsOpenSqKM", "ReachSlopePCNT")]

pairs(Pairs2, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)


#-------------Analysis---------------------
# STEP 1: Try full, beyond optimal model to get a working autocorrelation

# run fixed effect only to get starting values
lmFull <- lm(temp ~ airTemp+airTempLagged1+airTempLagged2+
                 Latitude+Longitude+
                 Forest+ Agriculture+
                 BasinElevationM+ ReachSlopePCNT+  
                 CONUSWetland+ SurficialCoarseC+
                 dayl + srad + swe +
                 sin(2*pi/360*dOY) + cos(2*pi/360*dOY), data = etS) # TotDASqKM+ImpoundmentsOpenSqKM+ 
                 

system.time(lmeFull <- lme(temp ~ airTemp+airTempLagged1+airTempLagged2+
                 Latitude+Longitude+
                 Forest+ Agriculture+
                 BasinElevationM+ ReachSlopePCNT+ 
                 CONUSWetland+ SurficialCoarseC+
                 dayl + srad + swe +
                 sin(2*pi/360*dOY) + cos(2*pi/360*dOY), random = list(site = ~ airTemp + airTempLagged1 + airTempLagged2 + Forest, year = ~sin(2*pi/360*dOY) + cos(2*pi/360*dOY)), data = etS, na.action = "na.omit")) # |year/site ??? TotDASqKM+ ImpoundmentsOpenSqKM+ 
summary(lmeFull) # didn't converge


system.time(lmeFull <- lme(temp ~ airTemp+airTempLagged1+airTempLagged2+
                 Latitude+Longitude+
                 Forest+ Agriculture+
                 BasinElevationM+ ReachSlopePCNT+ TotDASqKM+ 
                 CONUSWetland+ SurficialCoarseC+ImpoundmentsOpenSqKM+ 
                 dayl + srad + swe, random = list(site = ~ airTemp, year = ~1), data = etS, na.action = "na.omit", control = lmeControl(msMaxIter =90, maxIter=200)
)) # |year/site ???
summary(lmeFull) 

system.time(lmeFullAR1 <- update(lmeFull, correlation = corAR1(form = ~day | site) ))



system.time(lmeFull <- lme(temp ~ airTemp+airTempLagged1+airTempLagged2+
                 Latitude+Longitude+ prcp+
                 Forest+ Agriculture+
                 BasinElevationM+ ReachSlopePCNT+ TotDASqKM+ 
                 CONUSWetland+ SurficialCoarseC+ImpoundmentsOpenSqKM+ 
                 dayl + srad + swe, random = list(site = ~ airTemp + airTempLagged1 + airTempLagged2, year = ~1), data = etS, na.action = "na.omit", control = lmeControl(msMaxIter =90, maxIter=200)
)) # |year/site ???
summary(lmeFull) 

system.time(lmeFullAR1 <- update(lmeFull, correlation = corAR1() )) # 14.8 hours
summary(lmeFullAR1)

system.time(lmeFullAR2 <- update(lmeFull, correlation = corARMA(p = 2, q = 0) )) # 17.5 hours
summary(lmeFullAR2)

cbind(fixef(lmeFull), fixef(lmeFullAR1), fixef(lmeFullAR2))

plot(ACF(lmeFull, maxLag = 400))
plot(ACF(lmeFullAR1, maxLag = 400))
plot(ACF(lmeFullAR2, maxLag = 400))


# consider adding season and/or 10 day average (or 10-day average min) temperature and interaction to account for winter.

# consider adding 5 day precip

# how to account for snowmelt?








pdDiag(~sin(2*pi*dayS)))

lmeTest6b <- lme(temp ~ sin(2*pi/360*dayS) + cos(2*pi/360*dayS), data = etS, random = ~1| site)





system.time(lmeTest1 <- lme(temp ~ airTemp, random = ~ 1|site, data = et2, method = "ML", na.action = "na.omit")) # |year/site ??? 186s (3 min)
summary(lmeTest1)
rmse(resid(lmeTest1))
qplot(fitted(lmeTest1), resid(lmeTest1))
df1 <- data.frame(fit = fitted(lmeTest1), res = resid(lmeTest1), et2)
ggplot(data = df1, aes(date, res)) + geom_point(size = 1) + facet_wrap(~site) #+ stat_smooth(method = gam, formula = y ~ s(x))

system.time(lmeTest1Cor <- lme(temp ~ airTemp, random = ~ 1|site, correlation = corARMA(p = 1, q = 0, form = ~date|site), data = et2, method = "ML", na.action = na.exclude))
#### Took 2 hours but used all physical memory (8 GB) and 65 GB of virtual memory with JUST West Brook data. Causing RStudio to crash because R would temporarily freeze. Was able to run in R. This seems wrong for a 15065 x 14 data frame ###
summary(lmeTest1Cor)
rmse(resid(lmeTest1Cor))
qplot(fitted(lmeTest1Cor), resid(lmeTest1Cor))
df1 <- data.frame(fit = fitted(lmeTest1Cor), res = resid(lmeTest1Cor), et2)
ggplot(data = df1, aes(date, res)) + geom_point(size = 1) + facet_wrap(~site)
anova(lmeTest1, lmeTest1Cor) # doesn't improve model fit

system.time(lmeTest1Cor2 <- lme(temp ~ airTemp + jDay, random = ~ 1|site, correlation = corARMA(p = 1, q = 0, form = ~jDay|site), data = et2, method = "REML", na.action = na.exclude))

system.time(lmeTest2 <- lme(temp ~ airTemp, random = list(site = ~ 1, year = ~1), data = et, method = "ML", na.action = "na.omit")) # |year/site ??? 186s (3 min)
summary(lmeTest2)
rmse(resid(lmeTest2))

system.time(lmeTest5 <- lme(temp~airTemp+airTempLagged1+airTempLagged2+ daylS + sradS + sweLS, random = ~ 1 + airTemp+airTempLagged1+airTempLagged2|site, data = et, method = "ML", na.action = "na.omit", correlation = corCAR1(form = ~date|site)))
summary(lmeTest5)

rmse(resid(lmeTest5))
mae(resid(lmeTest5))

system.time(lmeTest5 <- lme(temp~airTemp+airTempLagged1+airTempLagged2+ daylS + sradS + sweLS, random = ~ 1 + airTemp+airTempLagged1+airTempLagged2|site, data = et, method = "ML", na.action = "na.omit", correlation = corCAR1(form = ~date|site)))
summary(lmeTest5)

rmse(resid(lmeTest5))
mae(resid(lmeTest5))


# Pinheiro and Bates p240
ggplot(data = Ovary, aes(Time, follicles)) + geom_point(size = 1) + facet_wrap(~Mare)
fm10var.lme <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary, random = pdDiag(~sin(2*pi*Time)))
summary(fm10var.lme)

# ** Need to scale date or dOY or convert 2*pi to radians (divide by 360)** 
et2$day <- as.numeric(et2$date)
et2$dayS <- (et2$day - mean(et2$day)) / sd(et2$day)
etGdf <- groupedData(temp ~ sin(2*pi*dayS)|site, data = et2)

lmeTest6a <- lme(temp ~ sin(2*pi*dayS) + cos(2*pi*dayS), data = etGdf, random = pdDiag(~sin(2*pi*dayS)))

lmeTest6b <- lme(temp ~ sin(2*pi*dayS) + cos(2*pi*dayS), data = et2, random = ~sin(2*pi*dayS)| site)

cbind(coef(lmeTest6a), coef(lmeTest6b))
cbind(fixef(lmeTest6a), fixef(lmeTest6b)) # not sure why these are different



summary(lmeTest6a)
summary(lmeTest6b)

lmeTest6AR1 <- update(lmeTest6a, correlation = corAR1() )
summary(lmeTest6AR1)

anova(lmeTest6a, lmeTest6AR1)

qplot(fitted(lmeTest6AR1), resid(lmeTest6AR1))


# Example of invocation of functions
rmse(resid(lmeTest3))
mae(resid(lmeTest3))



