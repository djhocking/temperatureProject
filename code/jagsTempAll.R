rm(list=ls())

setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
#baseDir <- 'D:/projects/temperatureProject/'

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
dataLocalDir <- paste0(baseDir, 'localData')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

loadLocalData <- F
loadetS <- F

if(loadLocalData) {
  load(paste0(dataLocalDir, 'etSCT'))
} else {
  if(loadetS) {
    load(paste0(dataOutDir, 'etS.RData'))
  } else {
    load(paste0(dataOutDir, 'et.RData'))

str(et)
summary(et)

# replace missing spring and fall BP with means for clipping the data to the synchronized season
et1 <- et
et1[which(is.na(et$springBP)), "springBP"] <- mean(et$springBP, na.rm=T)
et1[is.na(et$fallBP), "fallBP"] <- mean(et$fallBP, na.rm=T)
et1 <- et1[which(et1$dOY >= et1$springBP & et1$dOY <= et1$fallBP), ]

# Make dataframe with just variables for modeling and order before standardizing
et2 <- et1[ , c("agency", "date", "AgencyID", "siteYear", "year", "site", "date", "springBP", "fallBP", "FEATUREID", "HUC_4", "HUC_8", "HUC_12", "temp", "Latitude", "Longitude", "airTemp", "airTempLagged1", "airTempLagged2", "prcp", "prcpLagged1", "prcpLagged2", "prcpLagged3", "dOY", "Forest", "Herbacious", "Agriculture", "Developed", "TotDASqKM", "ReachElevationM", "ImpoundmentsAllSqKM", "HydrologicGroupAB", "SurficialCoarseC", "CONUSWetland", "ReachSlopePCNT", "srad", "dayl", "swe")] #  

dim(et2)
summary(et2) # strange there are some temp <0 and lots of stream temp NA
et2 <- na.omit(et2)
dim(et2)

#### temp just use CT - later add ability to pick state ######
# et2 <- et2[et2$agency == "CTDEP", ]
# dim(et2)
#############

etS <- cbind(et2[ ,c(1:14)],
             apply(X = et2[ ,15:dim(et2)[2]], MARGIN=2,
                   FUN = function(x){(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}))

etS$fyear <- as.factor(etS$year)
etS$fsite <- as.factor(etS$site)

save(etS, et2, file = paste0(dataOutDir, 'etS.RData'))
}
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
    
    # YEAR EFFECTiS
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
                  lat = etS$Latitude,
                  lon = etS$Longitude)
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
                     airTemp = etS$airTemp, 
                     airTempLag1 = etS$airTempLagged1,
                     airTempLag2 = etS$airTempLagged2,
                     precip = etS$prcp,
                     precipLag1 = etS$prcpLagged1,
                     precipLag2 = etS$prcpLagged3,
                     drainage = etS$TotDASqKM,
                     forest = etS$Forest,
                     elevation = etS$ReachElevationM,
                     coarseness = etS$SurficialCoarseC,
                     wetland = etS$CONUSWetland,
                     impoundments = etS$ImpoundmentsAllSqKM,
                     swe = etS$swe)
variables.site <- names(X.site)
J <- length(unique(etS$site))
K <- length(variables.site)
n <- dim(etS)[1]
W.site <- diag(K)

# Random Year effects
#variables.year <- c("Intercept-year",
  #                  "dOY",
   #                 "dOY2",
    #                "dOY3")

X.year <- data.frame(intercept.year = 1, 
                     dOY = etS$dOY, 
                     dOY2 = etS$dOY^2,
                     dOY3 = etS$dOY^3)
variables.year <- names(X.year)
Ti <- length(unique(etS$year))
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
             temp = etS$temp,
             X.site = X.site, #as.matrix(X.site),
             X.year = as.matrix(X.year),
             site = as.factor(etS$site),
             year = as.factor(etS$year))

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

#M1 <- bugs(etS, )

n.burn = 3000
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
pdf("/Users/Dan/Dropbox/correlatedSlopes.pdf")
plot(M3[ , 1:50])
dev.off()

summary(M3[ , 1:50])

summary(M3)$statistics[ , "Mean"]

rm(out)

pairs(as.matrix(M3[ , c(1:8, 17:20)]))

#summary(M3)$statistics[ , "Mean"]

summary.stats <- summary(M3)$statistics
summary.stats[1:1000, 1:2]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, K.0+K+L, 2))
names(fix.ef) <- c("Mean", "Std. Error")
row.names(fix.ef) <- c(variables.fixed, variables.site, variables.year)
for(k in 1:K.0){
  fix.ef[k, ] <- summary.stats[paste0('B.0[',k,']') , c("Mean", "SD")]
}
for(k in 1:K){
  fix.ef[k+K.0, ] <- summary.stats[paste0('mu.site[',k,']') , c("Mean", "SD")]
}
for(l in 1:L){
  fix.ef[l+K.0+K, ] <- summary.stats[paste0('mu.year[',l,']') , c("Mean", "SD")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef.site <- as.data.frame(matrix(NA, K, 2))
names(ran.ef.site) <- c("Variance", "Std. Dev.")
row.names(ran.ef.site) <- variables.site
for(k in 1:K){
  ran.ef.site[k, 2] <- summary.stats[paste0('sigma.b.site[',k,']') , c("Mean")]
  ran.ef.site[k, 1] <- ran.ef.site[k, 2] ^ 2
}
ran.ef.site

# Make Random Effects Output like summary(lmer)
ran.ef.year <- as.data.frame(matrix(NA, L, 2))
names(ran.ef.year) <- c("Variance", "Std. Dev.")
row.names(ran.ef.year) <- variables.year
for(k in 1:L){
  ran.ef.year[k, 2] <- summary.stats[paste0('sigma.b.year[',k,']') , c("Mean")]
  ran.ef.year[k, 1] <- ran.ef.year[k, 2] ^ 2
}
ran.ef.year

# Make correlation matrix of random site effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- summary.stats[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.site <- round(cor.site, digits=3)
cor.site[upper.tri(cor.site, diag=TRUE)] <- ''
cor.site

# Make correlation matrix of random year effects
cor.year <- as.data.frame(matrix(NA, L, L))
names(cor.year) <- variables.year
row.names(cor.year) <- variables.year
for(k in 1:L){
  for(k.prime in 1:L){
    cor.year[k, k.prime] <- summary.stats[paste('rho.B.year[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.year <- round(cor.year, digits=3)
cor.year[upper.tri(cor.year, diag=TRUE)] <- ''
cor.year

pred.t <- as.data.frame(matrix(NA, n, 2))
for(i in 1:n){
  pred.t[i, 1] <- summary.stats[paste0('stream.mu[',i,']') , c("Mean")]
  print(i)
}
pred.t[ , 2] <- etS$site
pred.t$airTemp <- (etS$airTemp * sd(et2$airTemp)) + mean(et2$airTemp)
pred.t$year <- etS$year
pred.t$dOY <- etS$dOY
names(pred.t) <- c("streamTemp", "site", "airTemp", "year", "dOY")

all.data <- data.frame(et2, pred.t)

ggplot(all.data, aes(airTemp, streamTemp)) + geom_point(aes(colour = dOY, size=0.5, alpha=0.5)) + facet_grid(site ~ year)

ggplot(all.data, aes(dOY, streamTemp)) + geom_point(size=0.75) + geom_point(data=etS, aes(dOY, temp), colour = "red", size=0.75) + facet_grid(site ~ year)

ggplot(all.data, aes(dOY, streamTemp, colour = year)) + geom_line(size=0.5)


ggplot(all.data, aes(dOY, streamTemp, colour = year)) + geom_point(aes(group = year))

ggplot(all.data, aes(airTemp, streamTemp)) + geom_line(size = 0.5, alpha = 0.5) + theme_bw() + theme(legend.position="none")
ggplot(all.data, aes(airTemp, streamTemp, colour = site)) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() + theme(legend.position="none")

ggplot(all.data, aes(airTemp, streamTemp, group = site)) + 
  geom_line(stat='smooth', method='lm', se=F, alpha=0.1, size=0.5, colour='black') + 
  theme_bw() +
  theme(legend.position="none") + 
  xlab('Air Temperature (C)') + 
  ylab('Predicted Stream Temperature (C)')

# dOY.eff.site <- 

err <- pred.t$streamTemp - etS$temp
rmse(err)

# compare with lme4
# Random slopes for crossed
system.time(lmer7 <- lmer(temp ~ TotDASqKM + Forest + ReachElevationM + airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 + prcp + prcpLagged1 + prcpLagged2 | site) + (dOY + I(dOY^2) + I(dOY^3)|year), data = etS))
summary(lmer7)
ranef(lmer7)

for(j in 1:J){
  air.eff.site[j] <- sum(ranef(lmer7)$site[j, 2:4])
}
air.eff.site
names(air.eff.site) <- row.names(rand.site)

rand.site <- ranef(lmer7)$site
rand.site$site <- row.names(rand.site)

library(dplyr)
by.site <- group_by(rand.site, site)
summarise(by.site, sum())



###### Consider: Slope, Aspect, Dams/Impoundments, Agriculture, Wetland, Coarseness, dayl, srad, swe

####### Add correlation structure for Lat-Lon and/or HUC ########





