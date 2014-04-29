rm(list=ls())

library(ggplot2)
#library(relaimpo)
#library(lme4)
#library(plyr)
#library(reshape)
#library(ggmap)
#library(foreign)
#library(maptools)
#library(gridExtra)
library(DataCombine) # for the slide function
library(nlme)

setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')

baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

load("dataOut/etWB.RData")

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

et2 <- et[ , c("year", "site", "date", "dOY", "temp", "flow", "airTemp", "prcp", "sradS", "daylS", "sweLS", "season", "airTempLagged1", "airTempLagged2")]
et2 <- na.omit(et2)

system.time(lmeTest1 <- lme(temp ~ airTemp, random = ~ 1|site, data = et2, method = "ML", na.action = "na.omit")) # |year/site ??? 186s (3 min)
summary(lmeTest1)
rmse(resid(lmeTest1))
qplot(fitted(lmeTest1), resid(lmeTest1))
df1 <- data.frame(fit = fitted(lmeTest1), res = resid(lmeTest1), et2)
ggplot(data = df1, aes(date, res)) + geom_point(size = 1) + facet_wrap(~site) #+ stat_smooth(method = gam, formula = y ~ s(x))

system.time(lmeTest1Cor <- lme(temp ~ airTemp + date, random = ~ 1|site, correlation = corARMA(p = 1, q = 0, form = ~date|site), data = et2, method = "REML", na.action = na.exclude))
summary(lmeTest1Cor)
rmse(resid(lmeTest1Cor))
qplot(fitted(lmeTest1Cor), resid(lmeTest1Cor))
df1 <- data.frame(fit = fitted(lmeTest1Cor), res = resid(lmeTest1Cor), et2)
ggplot(data = df1, aes(date, res)) + geom_point(size = 1) + facet_wrap(~site)

system.time(lmeTest1Cor <- lme(temp ~ airTemp + jDay, random = ~ 1|site, correlation = corARMA(p = 1, q = 0, form = ~jDay|site), data = et2, method = "REML", na.action = na.exclude))

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



# Example of invocation of functions
rmse(resid(lmeTest3))
mae(resid(lmeTest3))



