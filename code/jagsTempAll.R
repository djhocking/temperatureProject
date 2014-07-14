rm(list=ls())

library(ggplot2)
library(dplyr)
library(nlme)

#setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
#baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
baseDir <- 'C:/Users/dhocking/Documents/temperatureProject/'
setwd(baseDir)

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
dataLocalDir <- paste0(baseDir, 'localData')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

loadLocalData <- F

if(loadLocalData) {
  load(paste0(dataLocalDir, 'etSCT'))
} else {
  load(paste0(dataOutDir, 'tempDataSync.RData'))
}


########## JAGS Model ##############
sink("code/correlatedSlopes.txt")
cat("
    model{
    # Likelihood
    for(i in 1:n){ # n observations
    temp[i] ~ dnorm(stream.mu[i], tau)
    stream.mu[i] <- inprod(B.0[], X.0[i, ]) + inprod(B.site[site[i], ], X.site[i, ]) + inprod(B.year[year[i], ], X.year[i, ]) #  
    }
    
    # prior for model variance
    sigma ~ dunif(0, 100)
    tau <- pow(sigma, -2)
    
    for(k in 1:K.0){
    B.0[k] ~ dnorm(0, 0.001) # priors coefs for fixed effect predictors
    }
    
    # Priors for random effects of site
    for(j in 1:J){ # J sites
    B.site[j, 1:K] ~ dmnorm(mu.site[ ], tau.B.site[ , ])
    }
    mu.site[1] <- 0
    for(k in 2:K){
    mu.site[k] ~ dnorm(0, 0.0001)
    }
    
    # Prior on multivariate normal std deviation
    tau.B.site[1:K, 1:K] ~ dwish(W.site[ , ], df.site)
    df.site <- K + 1
    sigma.B.site[1:K, 1:K] <- inverse(tau.B.site[ , ])
    for(k in 1:K){
    for(k.prime in 1:K){
    rho.B.site[k, k.prime] <- sigma.B.site[k, k.prime]/sqrt(sigma.B.site[k, k]*sigma.B.site[k.prime, k.prime])
    }
    sigma.b.site[k] <- sqrt(sigma.B.site[k, k])
    }
    
    # YEAR EFFECTS
    # Priors for random effects of year
    for(t in 1:Ti){ # Ti years
    B.year[t, 1:L] ~ dmnorm(mu.year[ ], tau.B.year[ , ])
    }
    mu.year[1] <- 0
    for(l in 2:L){
    mu.year[l] ~ dnorm(0, 0.0001)
    }
    
    # Prior on multivariate normal std deviation
    tau.B.year[1:L, 1:L] ~ dwish(W.year[ , ], df.year)
    df.year <- L + 1
    sigma.B.year[1:L, 1:L] <- inverse(tau.B.year[ , ])
    for(l in 1:L){
    for(l.prime in 1:L){
    rho.B.year[l, l.prime] <- sigma.B.year[l, l.prime]/sqrt(sigma.B.year[l, l]*sigma.B.year[l.prime, l.prime])
    }
    sigma.b.year[l] <- sqrt(sigma.B.year[l, l])
    }

    }
    ", fill = TRUE)
sink()

# Fixed effects
#variables.fixed <- c("intercept",  
                     #"drainage", 
                     #"forest",
                     #"elevation")
#K.0 <- length(variables.fixed)
X.0 <- data.frame(intercept = 1,
                  lat = tempDataSyncS$Latitude,
                  lon = tempDataSyncS$Longitude)
variables.fixed <- names(X.0)
K.0 <- length(variables.fixed)


# Random site effects
#variables.site <- c("Intercept-site",
 #                   "Air Temperature",
  #                  "Air Temp Lag1",
   #                 "Air Temp Lag2",
    #                "Precip",
     #               "Precip Lag1",
      #              "Precip Lag2")

# Slope, Aspect, Dams/Impoundments, Agriculture, Wetland, Coarseness, dayl, srad, swe

X.site <- data.frame(intercept.site = 1, 
                     airTemp = tempDataSyncS$airTemp, 
                     airTempLag1 = tempDataSyncS$airTempLagged1,
                     airTempLag2 = tempDataSyncS$airTempLagged2,
                     precip = tempDataSyncS$prcp,
                     precipLag1 = tempDataSyncS$prcpLagged1,
                     precipLag2 = tempDataSyncS$prcpLagged3,
                     drainage = tempDataSyncS$TotDASqKM,
                     forest = tempDataSyncS$Forest,
                     elevation = tempDataSyncS$ReachElevationM,
                     coarseness = tempDataSyncS$SurficialCoarseC,
                     wetland = tempDataSyncS$CONUSWetland,
                     impoundments = tempDataSyncS$ImpoundmentsAllSqKM,
                     swe = tempDataSyncS$swe)
variables.site <- names(X.site)
J <- length(unique(tempDataSyncS$site))
K <- length(variables.site)
n <- dim(tempDataSyncS)[1]
W.site <- diag(K)

# Random Year effects
#variables.year <- c("Intercept-year",
  #                  "dOY",
   #                 "dOY2",
    #                "dOY3")

X.year <- data.frame(intercept.year = 1, 
                     dOY = tempDataSyncS$dOY, 
                     dOY2 = tempDataSyncS$dOY^2,
                     dOY3 = tempDataSyncS$dOY^3)
variables.year <- names(X.year)
Ti <- length(unique(tempDataSyncS$year))
L <- length(variables.year)
W.year <- diag(L)

data <- list(n = n, 
             J = J, 
             K = K, 
             Ti = Ti,
             L = L,
             K.0 = K.0,
             X.0 = X.0,
             W.site = W.site,
             W.year = W.year,
             temp = tempDataSyncS$temp,
             X.site = X.site, #as.matrix(X.site),
             X.year = as.matrix(X.year),
             site = as.factor(tempDataSyncS$site),
             year = as.factor(tempDataSyncS$year))

inits <- function(){
  list(#B.raw = array(rnorm(J*K), c(J,K)), 
    #mu.site.raw = rnorm(K),
    sigma = runif(1),
    #tau.B.site.raw = rwish(K + 1, diag(K)),
    xi = runif(K))
}

params <- c("sigma",
            "B.0",
            "B.site",
            "rho.B.site",
            "mu.site",
            "sigma.b.site",
            "B.year",
            "rho.B.year",
            "mu.year",
            "sigma.b.year",
            "stream.mu")

#M1 <- bugs(tempDataSyncS, )

n.burn = 5000
n.it = 3000
n.thin = 3

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data", "inits", "params", "K", "J", "Ti", "L", "n", "W.site", "W.year", "X.site", "X.year", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

system.time(out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/correlatedSlopes.txt", data, inits, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
}))

M3 <- mcmc.list(out)
stopCluster(CL)

#pdf("/Users/Dan/Dropbox/correlatedSlopes.pdf")
pdf("C:/Users/dhocking/Dropbox/correlatedSlopes.pdf")
plot(M3[ , 1:50])
dev.off()

rm(out)

pairs(as.matrix(M3[ , c(1:8, 17:20)]))

memory.limit(size = 1e6)

summary.stats <- summary(M3)
summary.stats[1:1000, 1:2]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, K.0+K+L, 4))
names(fix.ef) <- c("Mean", "Std. Error", "LCI", "UCI")
row.names(fix.ef) <- c(variables.fixed, variables.site, variables.year)
for(k in 1:K.0){
  fix.ef[k, 1:2] <- summary.stats$statistics[paste0('B.0[',k,']') , c("Mean", "SD")]
  fix.ef[k, 3:4] <- summary.stats$quantiles[paste0('B.0[',k,']') , c("2.5%", "97.5%")]
}
for(k in 1:K){
  fix.ef[k+K.0, 1:2] <- summary.stats$statistics[paste0('mu.site[',k,']') , c("Mean", "SD")]
  fix.ef[k+K.0, 3:4] <- summary.stats$quantiles[paste0('mu.site[',k,']') , c("2.5%", "97.5%")]
}
for(l in 1:L){
  fix.ef[l+K.0+K, 1:2] <- summary.stats$statistics[paste0('mu.year[',l,']') , c("Mean", "SD")]
  fix.ef[l+K.0+K, 3:4] <- summary.stats$quantiles[paste0('mu.year[',l,']') , c("2.5%", "97.5%")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef.site <- as.data.frame(matrix(NA, K, 2))
names(ran.ef.site) <- c("Variance", "SD")
row.names(ran.ef.site) <- variables.site
for(k in 1:K){
  ran.ef.site[k, 2] <- summary.stats$statistics[paste0('sigma.b.site[',k,']') , c("Mean")]
  ran.ef.site[k, 1] <- ran.ef.site[k, 2] ^ 2
}
ran.ef.site

# Make Random Effects Output like summary(lmer)
ran.ef.year <- as.data.frame(matrix(NA, L, 2))
names(ran.ef.year) <- c("Variance", "SD")
row.names(ran.ef.year) <- variables.year
for(l in 1:L){
  ran.ef.year[l, 2] <- summary.stats$statistics[paste0('sigma.b.year[',l,']') , c("Mean")]
  ran.ef.year[l, 1] <- ran.ef.year[l, 2] ^ 2
}
ran.ef.year

# Make correlation matrix of random site effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- summary.stats$statistics[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.site <- round(cor.site, digits=3)
cor.site[upper.tri(cor.site, diag=TRUE)] <- ''
cor.site

# Make correlation matrix of random year effects
cor.year <- as.data.frame(matrix(NA, L, L))
names(cor.year) <- variables.year
row.names(cor.year) <- variables.year
for(l in 1:L){
  for(l.prime in 1:L){
    cor.year[l, l.prime] <- summary.stats$statistics[paste('rho.B.year[',l,',',l.prime,']', sep=""), "Mean"]
  }
}
cor.year <- round(cor.year, digits=3)
cor.year[upper.tri(cor.year, diag=TRUE)] <- ''
cor.year

# combine model summary results into an S4 Object
setClass("jagsSummary",
         representation(fixEf="data.frame",
                        ranEf="list",
                        ranCor="list"))

modSummary <- new("jagsSummary")
modSummary@fixEf <- fix.ef
modSummary@ranEf <- list(ranSite=ran.ef.site, ranYear=ran.ef.year)
modSummary@ranCor <- list(corSite=cor.site, corYear=cor.year)

modSummary
str(modSummary)

# Add predicted stream temperatures to dataframe
tempDataSync$streamTempPred <- NA
tempDataSync$streamTempPredLCI <- NA
tempDataSync$streamTempPredUCI <- NA
for(i in 1:n){
  tempDataSync$streamTempPred[i] <- summary.stats$statistics[paste0('stream.mu[',i,']') , c("Mean")]
  tempDataSync$streamTempPredLCI[i] <- summary.stats$quantiles[paste0('stream.mu[',i,']') , c("2.5%")]
  tempDataSync$streamTempPredUCI[i] <- summary.stats$quantiles[paste0('stream.mu[',i,']') , c("97.5%")]
}

########## Check model fit #############
tempDataSync$err <- tempDataSync$streamTempPred - tempDataSync$temp
rmse(tempDataSync$err)

# Add fit for validation data

########################################

############ Plots ###############
ggplot(tempDataSync, aes(temp, streamTempPred)) + geom_point() +geom_abline(slope=1, intercept=0, colour='red')

ggplot(tempDataSync, aes(airTemp, streamTempPred)) + geom_point(aes(colour = dOY, size=0.5, alpha=0.5)) + facet_grid(. ~ year)

ggplot(tempDataSync, aes(dOY, streamTempPred)) + geom_point(size=0.75) + geom_point(data=tempDataSync, aes(dOY, temp), colour = "red", size=0.75) + facet_grid(. ~ year)

ggplot(tempDataSync, aes(dOY, streamTempPred, colour = year)) + geom_line(size=0.5)


ggplot(tempDataSync, aes(dOY, streamTempPred, colour = year)) + geom_point(aes(group = year))

ggplot(tempDataSync, aes(airTemp, streamTempPred)) + geom_line(size = 0.5, alpha = 0.3) + theme_bw() + theme(legend.position="none")
ggplot(tempDataSync, aes(airTemp, streamTempPred, colour = site)) + geom_point(size = 0.5, alpha = 0.8) + theme_bw() + theme(legend.position="none")

ggplot(tempDataSync, aes(airTemp, streamTempPred, group = site)) + 
  geom_line(stat='smooth', method='lm', se=F, alpha=0.2, size=0.5, colour='black') + 
  theme_bw() +
  theme(legend.position="none") + 
  xlab('Air Temperature (C)') + 
  ylab('Predicted Stream Temperature (C)')

# plot observed and predicte vs day of the year for all sites
sites <- unique(tempDataSync$site)

for(i in 1:length(unique(tempDataSync$site))){
  dataSite <- filter(tempDataSync, filter = site == sites[i])
  foo <- ggplot(dataSite, aes(dOY, temp)) + coord_cartesian(xlim = c(50, 350), ylim = c(0, 30)) + geom_point(colour = 'blue') + geom_line(colour = 'blue') + geom_point(aes(dOY, streamTempPred), colour = 'red', size=1) + geom_line(aes(dOY, streamTempPred), colour = 'red', size=0.1) + geom_point(aes(dOY, airTemp), colour='black', size=0.1) + ggtitle(unique(tempDataSync$fsite)[i]) + facet_wrap(~year)
  ggsave(filename=paste0(dataLocalDir,'/', 'plots/', unique(tempDataSync$fsite)[i], '.png'), plot=foo, dpi=300 , width=8,height=5, units='in' )
} # surprisingly fast



###### Consider: Slope, Aspect, Dams/Impoundments, Agriculture, Wetland, Coarseness, dayl, srad, swe

####### Add correlation structure for Lat-Lon and/or HUC ########

##### Add autoregressive component to the model? #########

##### Add interactions such as airTemp*drainage





