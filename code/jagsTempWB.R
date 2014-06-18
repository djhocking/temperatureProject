rm(list=ls())

setwd('/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/')
#setwd('C:/Users/dhocking/Documents/temperatureProject/')

#baseDir <- 'C:/KPONEIL/GitHub/projects/temperatureProject/'
baseDir <- '/Users/Dan/Documents/Research/Stream_Climate_Change/temperatureProject/'
#baseDir <- 'D:/projects/temperatureProject/'

dataInDir <- paste0(baseDir, 'dataIn/')
dataOutDir <- paste0(baseDir, 'dataOut/')
graphsDir <- paste0(baseDir, 'graphs/')

source(paste0(baseDir, 'code/functions/temperatureModelingFunctions.R'))

load(paste0(dataOutDir, 'etSWB.RData'))

######### Non-nested model with random intercepts #########
sink("code/cubicDay.txt")
cat("
    model{
    # Priors
    alpha ~ dnorm(0, 0.01)
    b.air ~ dnorm(0, 0.01)
    b.air1 ~ dnorm(0, 0.01)
    b.air2 ~ dnorm(0, 0.01)
    b.flow ~ dnorm(0, 0.01)
    b.dOY ~ dnorm(0, 0.01)
    b.dOY2 ~ dnorm(0, 0.01)
    b.dOY3 ~ dnorm(0, 0.01)

    sigma ~ dunif(0, 100)
    tau <- 1/(sigma*sigma)
    
    # Hyperpriors for Random Effects
    sigma.site ~ dunif(0, 100)
    tau.site <- 1 / (sigma.site * sigma.site)
    #tau.site ~ dgamma(0.01, 0.01)    
    #sigma.site <- pow(1/tau.site, 0.5) # sqrt of 1/tau
    
    sigma.year ~ dunif(0, 100)
    tau.year <- 1 / (sigma.year * sigma.year)
    
    for(i in 1:n.sites) {
      b.site[i] ~ dnorm(0, tau.site)
    }
    
    for(i in 1:n.years) {
      b.year[i] ~ dnorm(0, tau.year)
    }
    
    # Autocorrelation structure
    
    # Linear Model
    for(i in 1:n.obs) {
      stream.mu[i] <- alpha + b.site[site[i]] + 
                              b.year[year[i]] + 
                              b.air*airTemp[i] + b.air1*airTempLagged1[i] + b.air2*airTempLagged2[i] + 
                              b.flow*flow[i] +
                              b.dOY*dOY[i] + b.dOY2*dOY2[i] + b.dOY3*dOY3[i] 
                         
      temp[i] ~ dnorm(stream.mu[i], tau)
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
                 "b.year",
                 "b.site",
                 "b.air",
                 "b.air1",
                 "b.air2",
                 "b.flow",
                 "b.dOY",
                 "b.dOY2",
                 "b.dOY3",
                 "sigma.site",
                 "sigma.year",
                 "sigma")

data1 <- list(temp = etS$temp,
              airTemp = etS$airTemp,
              airTempLagged1 = etS$airTempLagged1,
              airTempLagged2 = etS$airTempLagged2,
              flow = etS$flow,
              dOY = etS$dOY,
              dOY2 = etS$dOY^2,
              dOY3 = etS$dOY^3,
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



n.burn = 500
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

o <- mcmc.list(out)
stopCluster(CL)
plot(o)

rm(out)


######## Non-nested model with random intercepts and slopes #########
sink("code/cubicDaySlopes.txt")
cat("
    model{
    # Priors
    alpha ~ dnorm(0, 0.01)
    b.flow ~ dnorm(0, 0.01)
    b.dOY ~ dnorm(0, 0.01)
    b.dOY2 ~ dnorm(0, 0.01)
    b.dOY3 ~ dnorm(0, 0.01)
    
    sigma ~ dunif(0, 100)
    tau <- 1/(sigma*sigma)
    
    # Hyperpriors for Random Effects
    sigma.site ~ dunif(0, 100)
    tau.site <- 1 / (sigma.site * sigma.site)

    sigma.air ~ dunif(0, 100)
    tau.air <- 1 / (sigma.air * sigma.air)

    sigma.air1 ~ dunif(0, 100)
    tau.air1 <- 1 / (sigma.air1 * sigma.air1)

    sigma.air2 ~ dunif(0, 100)
    tau.air2 <- 1 / (sigma.air2 * sigma.air2)
    
    sigma.year ~ dunif(0, 100)
    tau.year <- 1 / (sigma.year * sigma.year)
    
    for(j in 1:n.sites) {
    b.site[j] ~ dnorm(0, tau.site)
    b.air[j] ~ dnorm(0, tau.air)
    b.air1[j] ~ dnorm(0, tau.air1)
    b.air2[j] ~ dnorm(0, tau.air2)
    }
    
    for(t in 1:n.years) {
    b.year[t] ~ dnorm(0, tau.year)
    }
    
    # Autocorrelation structure
    
    # Linear Model
    for(i in 1:n.obs) {
    stream.mu[i] <- alpha + b.site[site[i]] + b.year[year[i]] + b.air[site[i]]*airTemp[i] + b.air1[site[i]]*airTempLagged1[i] + b.air2[site[i]]*airTempLagged2[i] + b.flow*flow[i] + b.dOY*dOY[i] + b.dOY2*dOY2[i] + b.dOY3*dOY3[i] 
    
    temp[i] ~ dnorm(stream.mu[i], tau)
    }
  }
    ", fill = TRUE)
sink()

n.obs <- dim(etS)[1]
n.sites <- length(unique(etS$site))

inits1 <- function() {
  list(b.air = rnorm(n.sites, 3.5, 0.5))
}

params1 <- c(    "alpha",
                 "b.year",
                 "b.site",
                 "b.air",
                 "b.air1",
                 "b.air2",
                 "b.flow",
                 "b.dOY",
                 "b.dOY2",
                 "b.dOY3",
                 "sigma.site",
                 "sigma.year",
                 "sigma.air",
                 "sigma.air1",
                 "sigma.air2",
                 "sigma")

data1 <- list(temp = etS$temp,
              airTemp = etS$airTemp,
              airTempLagged1 = etS$airTempLagged1,
              airTempLagged2 = etS$airTempLagged2,
              flow = etS$flow,
              dOY = etS$dOY,
              dOY2 = etS$dOY^2,
              dOY3 = etS$dOY^3,
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



n.burn = 1000
n.it = 3000
n.thin = 3

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data1", "inits1", "params1", "n.obs", "n.sites", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/cubicDaySlopes.txt", data1, inits1, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params1, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
})

o <- mcmc.list(out)
stopCluster(CL)
pdf("/Users/Dan/Dropbox/cubicSlopes.pdf")
plot(o)
dev.off()

rm(out)


######### Non-nested model with random intercepts and slopes ########
sink("code/cubicDaySlopes.txt")
cat("
    model{
    # Priors
    alpha ~ dnorm(0, 0.01)
    b.flow ~ dnorm(0, 0.01)
    
    sigma ~ dunif(0, 100)
    tau <- 1/(sigma*sigma)
    
    # Hyperpriors for Random Effects
    sigma.site ~ dunif(0, 100)
    tau.site <- 1 / (sigma.site * sigma.site)
    
    mu.air ~ dnorm(0, 0.001)
    sigma.air ~ dunif(0, 100)
    tau.air <- 1 / (sigma.air * sigma.air)
    
    mu.air1 ~ dnorm(0, 0.001)
    sigma.air1 ~ dunif(0, 100)
    tau.air1 <- 1 / (sigma.air1 * sigma.air1)
    
    mu.air2 ~ dnorm(0, 0.001)
    sigma.air2 ~ dunif(0, 100)
    tau.air2 <- 1 / (sigma.air2 * sigma.air2)
    
    sigma.year ~ dunif(0, 100)
    tau.year <- 1 / (sigma.year * sigma.year)
    
    mu.dOY ~ dnorm(0, 0.001)
    sigma.dOY ~ dunif(0, 100)
    tau.dOY <- 1 / (sigma.dOY * sigma.dOY)
    
    mu.dOY2 ~ dnorm(0, 0.001)
    sigma.dOY2 ~ dunif(0, 100)
    tau.dOY2 <- 1 / (sigma.dOY2 * sigma.dOY2)
    
    mu.dOY3 ~ dnorm(0, 0.001)
    sigma.dOY3 ~ dunif(0, 100)
    tau.dOY3 <- 1 / (sigma.dOY3 * sigma.dOY3)
    
    for(j in 1:n.sites) {
    b.site[j] ~ dnorm(0, tau.site)
    b.air[j] ~ dnorm(mu.air, tau.air)
    b.air1[j] ~ dnorm(mu.air1, tau.air1)
    b.air2[j] ~ dnorm(mu.air2, tau.air2)
    }
    
    for(t in 1:n.years) {
    b.year[t] ~ dnorm(0, tau.year)
    b.dOY[t] ~ dnorm(mu.dOY, tau.dOY)
    b.dOY2[t] ~ dnorm(mu.dOY2, tau.dOY2)
    b.dOY3[t] ~ dnorm(mu.dOY3, tau.dOY3)
    }
    
    # Autocorrelation structure
    
    # Linear Model
    for(i in 1:n.obs) {
    stream.mu[i] <- alpha + b.site[site[i]] + b.year[year[i]] + b.air[site[i]]*airTemp[i] + b.air1[site[i]]*airTempLagged1[i] + b.air2[site[i]]*airTempLagged2[i] + b.flow*flow[i] + b.dOY[year[i]]*dOY[i] + b.dOY2[year[i]]*dOY2[i] + b.dOY3[year[i]]*dOY3[i] 
    
    temp[i] ~ dnorm(stream.mu[i], tau)
    }
    }
    ", fill = TRUE)
sink()

n.obs <- dim(etS)[1]
n.sites <- length(unique(etS$site))

inits1 <- function() {
  list(b.air = rnorm(n.sites, 2.5, 0.5))
}

params1 <- c(    "alpha",
                 "b.year",
                 "b.site",
                 "b.air",
                 "b.air1",
                 "b.air2",
                 "b.flow",
                 "b.dOY",
                 "b.dOY2",
                 "b.dOY3",
                 "sigma.site",
                 "sigma.year",
                 "mu.dOY",
                 "mu.dOY2",
                 "mu.dOY3",
                 "sigma.dOY",
                 "sigma.dOY2",
                 "sigma.dOY3",
                 "mu.air",
                 "mu.air1",
                 "mu.air2",
                 "sigma.air",
                 "sigma.air1",
                 "sigma.air2",
                 "sigma")

data1 <- list(temp = etS$temp,
              airTemp = etS$airTemp,
              airTempLagged1 = etS$airTempLagged1,
              airTempLagged2 = etS$airTempLagged2,
              flow = etS$flow,
              dOY = etS$dOY,
              dOY2 = etS$dOY^2,
              dOY3 = etS$dOY^3,
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
n.it = 1000
n.thin = 1

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data1", "inits1", "params1", "n.obs", "n.sites", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/cubicDaySlopes.txt", data1, inits1, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params1, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
})

o <- mcmc.list(out)
stopCluster(CL)
pdf("/Users/Dan/Dropbox/cubicSlopes.pdf")
plot(o)
dev.off()

rm(out)


############# Correlated Random Intercepts and Slopes #############
# Inverse Wishart and Multivariate normal (Gelman and Hill p378)

# Random Site Effects
sink("code/correlatedSlopesSites.txt")
cat("
    model{
    # Likelihood
    for(i in 1:n){ # n observations
    temp[i] ~ dnorm(stream.mu[i], tau)
    stream.mu[i] <- inprod(B.site[site[i], ], X.site[i, ]) #+ inprod(B.year[year[i], ], X.year[i, ]) # inprod(b.0[], X.0[i, ]) + 
    }
    
    # prior for model variance
    sigma ~ dunif(0, 100)
    tau <- pow(sigma, -2)
    
    # priors coefs for fixed effect predictors
    #for(k in 1:K.0){
    # b.0[k] ~ dnorm(0, 0.001) 
    #}
    
    # Priors for random effects of site
    for(j in 1:J){ # J sites
    for(k in 1:K){ # K random site effects
    B.site[j, k] <- xi[k]*B.site.raw[j, k]
    }
    B.site.raw[j, 1:K] ~ dmnorm(mu.site.raw[ ], tau.B.site.raw[ , ])
    }
    for(k in 1:K){
    mu.site[k] <- xi[k]*mu.site.raw[k]
    mu.site.raw[k] ~ dnorm(0, 0.001)
    xi[k] ~ dunif(0, 100)
    }
    
    # Prior on multivariate normal std deviation
    tau.B.site.raw[1:K, 1:K] ~ dwish(W.site[ , ], df.site)
    df.site <- K + 1
    sigma.B.site.raw[1:K, 1:K] <- inverse(tau.B.site.raw[ , ])
    for(k in 1:K){
    for(k.prime in 1:K){
    rho.B.site[k, k.prime] <- sigma.B.site.raw[k, k.prime]/sqrt(sigma.B.site.raw[k, k]*sigma.B.site.raw[k.prime, k.prime])
    }
    sigma.B.site[k] <- abs(xi[k])*sqrt(sigma.B.site.raw[k, k])
    }
    }
    ", fill = TRUE)
sink()

variables.site <- c("Intercept-site",
                    "Air Temperature",
                    "Air Temp Lag1",
                    "Air Temp Lag2")
J <- length(unique(etS$site))
K <- length(variables.site)
n <- dim(etS)[1]
W.site <- diag(K)
X.site <- data.frame(int = 1, 
                     airT = etS$airTemp, 
                     airT1 = etS$airTempLagged1,
                     airT2 = etS$airTempLagged2)

data <- list(n = n, 
             J = J, 
             K = K, 
             Ti = Ti,
             L = L,
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
            "B.site",
            "rho.B.site",
            "mu.site",
            "sigma.B.site")

#M1 <- bugs(etS, )

n.burn = 2000
n.it = 3000
n.thin = 3

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data", "inits", "params", "K", "J", "Ti", "L", "n", "W.site", "X.site", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

system.time(out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/correlatedSlopesSites.txt", data, inits, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
}))

M1 <- mcmc.list(out)
stopCluster(CL)
pdf("/Users/Dan/Dropbox/correlatedSlopesSites.pdf")
plot(M1)
dev.off()

summary(M1)

summary(M1)$statistics[ , "Mean"]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, k, 2))
names(fix.ef) <- c("Mean", "Std. Error")
row.names(fix.ef) <- variables.site
for(k in 1:K){
  fix.ef[k, ] <- summary(M1)$statistics[paste0('mu.site[',k,']') , c("Mean", "SD")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef <- as.data.frame(matrix(NA, k, 2))
names(ran.ef) <- c("Variance", "Std. Dev.")
row.names(ran.ef) <- variables.site
for(k in 1:K){
  ran.ef[k, 2] <- summary(M1)$statistics[paste0('sigma.B.site[',k,']') , c("Mean")]
  ran.ef[k, 1] <- ran.ef[k, 2] ^ 2
}
ran.ef

# Make correlation matrix of random effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- summary(M1)$statistics[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.site <- round(cor.site, digits=3)
cor.site[upper.tri(cor.site, diag=TRUE)] <- ''
cor.site

rm(out)

# compare with lme4
# Random slopes for crossed
lmer6 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + (airTemp + airTempLagged1 + airTempLagged2 |site), data = etS)
summary(lmer6)
ranef(lmer6)



########## Add fixed flow to Random Site Effects ########
sink("code/correlatedSlopesSites.txt")
cat("
    model{
    # Likelihood
    for(i in 1:n){ # n observations
    temp[i] ~ dnorm(stream.mu[i], tau)
    stream.mu[i] <- inprod(b.0[], X.0[i, ]) + inprod(B.site[site[i], ], X.site[i, ]) #+ inprod(B.year[year[i], ], X.year[i, ]) #  
    }
    
    # prior for model variance
    sigma ~ dunif(0, 100)
    tau <- pow(sigma, -2)
    
    # priors coefs for fixed effect predictors
    for(k in 1:K.0){
     b.0[k] ~ dnorm(0, 0.001) 
    }
    
    # Priors for random effects of site
    for(j in 1:J){ # J sites
    for(k in 1:K){ # K random site effects
    B.site[j, k] <- xi[k]*B.site.raw[j, k]
    }
    B.site.raw[j, 1:K] ~ dmnorm(mu.site.raw[ ], tau.B.site.raw[ , ])
    }

    for(k in 1:K){
    mu.site[k] <- xi[k]*mu.site.raw[k]
    mu.site.raw[k] ~ dnorm(0, 0.001)
    xi[k] ~ dunif(0, 100)
    }
    
    # Prior on multivariate normal std deviation
    tau.B.site.raw[1:K, 1:K] ~ dwish(W.site[ , ], df.site)
    df.site <- K + 1
    sigma.B.site.raw[1:K, 1:K] <- inverse(tau.B.site.raw[ , ])
    for(k in 1:K){
    for(k.prime in 1:K){
    rho.B.site[k, k.prime] <- sigma.B.site.raw[k, k.prime]/sqrt(sigma.B.site.raw[k, k]*sigma.B.site.raw[k.prime, k.prime])
    }
    sigma.B.site[k] <- abs(xi[k])*sqrt(sigma.B.site.raw[k, k])
    }
    }
    ", fill = TRUE)
sink()

variables.site <- c("Intercept-site",
                    "Air Temperature",
                    "Air Temp Lag1",
                    "Air Temp Lag2")
J <- length(unique(etS$site))
K <- length(variables.site)
n <- dim(etS)[1]
W.site <- diag(K)
X.site <- data.frame(int.site = 1, 
                     airT = etS$airTemp, 
                     airT1 = etS$airTempLagged1,
                     airT2 = etS$airTempLagged2)

variables.fixed <- c("flow")
K.0 <- length(variables.fixed)
X.0 <- data.frame(flow = etS$flow)

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
            "B.site",
            "b.0",
            "rho.B.site",
            "mu.site",
            "sigma.B.site")

#M1 <- bugs(etS, )

n.burn = 100
n.it = 1000
n.thin = 1

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data", "inits", "params", "K", "J", "Ti", "L", "n", "W.site", "X.site", "n.burn", "n.it", "n.thin", "K.0", "X.0"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

system.time(out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/correlatedSlopesSites.txt", data, inits, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
}))

M2 <- mcmc.list(out)
stopCluster(CL)
pdf("/Users/Dan/Dropbox/correlatedSlopesSites.pdf")
plot(M2)
dev.off()

summary(M2)

summary(21)$statistics[ , "Mean"]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, k, 2))
names(fix.ef) <- c("Mean", "Std. Error")
row.names(fix.ef) <- variables.site
for(k in 1:K){
  fix.ef[k, ] <- summary(M1)$statistics[paste0('mu.site[',k,']') , c("Mean", "SD")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef <- as.data.frame(matrix(NA, k, 2))
names(ran.ef) <- c("Variance", "Std. Dev.")
row.names(ran.ef) <- variables.site
for(k in 1:K){
  ran.ef[k, 2] <- summary(M1)$statistics[paste0('sigma.B.site[',k,']') , c("Mean")]
  ran.ef[k, 1] <- ran.ef[k, 2] ^ 2
}
ran.ef

# Make correlation matrix of random effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- summary(M1)$statistics[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.site <- round(cor.site, digits=3)
cor.site[upper.tri(cor.site, diag=TRUE)] <- ''
cor.site

rm(out)

# compare with lme4
# Random slopes for crossed
lmer2 <- lmer(temp ~ flow + airTemp + airTempLagged1 + airTempLagged2 + (airTemp + airTempLagged1 + airTempLagged2 |site), data = etS)
summary(lmer2)
ranef(lmer2)


############# Site and Year Effects ###############
sink("code/correlatedSlopes.txt")
cat("
    model{
      # Likelihood
      for(i in 1:n){ # n observations
        temp[i] ~ dnorm(stream.mu[i], tau)
        stream.mu[i] <- inprod(B.site[site[i], ], X.site[i, ]) + inprod(B.year[year[i], ], X.year[i, ]) # inprod(b.0[], X.0[i, ]) + 
      }
      
      # Priors for fixed effects


      # prior for model variance
      sigma ~ dunif(0, 100)
      tau <- pow(sigma, -2)
      
      #for(k in 1:K.0){
       # b.0[k] ~ dnorm(0, 0.001) # priors coefs for fixed effect predictors
      #}
      
      # Priors for random effects of site
      for(j in 1:J){ # J sites
        for(k in 1:K){ # K random site effects
          B.site[j, k] <- xi[k]*B.site.raw[j, k]
        }
        B.site.raw[j, 1:K] ~ dmnorm(mu.site.raw[ ], tau.B.site.raw[ , ])
      }
      for(k in 2:K){
        mu.site[k] <- xi[k]*mu.site.raw[k]
        mu.site.raw[k] ~ dnorm(0, 0.0001)
        xi[k] ~ dunif(0, 100)
      }
        mu.site[1] <- 0 #xi[1]*mu.site.raw[1]
      #  mu.site.raw[1] <- 0
        #xi[1] <- 1

      
      # Prior on multivariate normal std deviation
      tau.B.site.raw[1:K, 1:K] ~ dwish(W.site[ , ], df.site)
      df.site <- K + 1
      sigma.B.site.raw[1:K, 1:K] <- inverse(tau.B.site.raw[ , ])
      for(k in 1:K){
        for(k.prime in 1:K){
          rho.B.site[k, k.prime] <- sigma.B.site.raw[k, k.prime]/sqrt(sigma.B.site.raw[k, k]*sigma.B.site.raw[k.prime, k.prime])
        }
        sigma.B.site[k] <- abs(xi[k])*sqrt(sigma.B.site.raw[k, k])
      }
      
      # YEAR EFFECTiS
      # Priors for random effects of year
      for(t in 1:Ti){ # Ti years
        for(l in 1:L){ # L random year effects
          B.year[t, l] <- xt[l]*B.year.raw[t, l]
        }
        B.year.raw[t, 1:L] ~ dmnorm(mu.year.raw[ ], tau.B.year.raw[ , ])
      }
      for(l in 1:L){
        mu.year[l] <- xt[l]*mu.year.raw[l]
        mu.year.raw[l] ~ dnorm(0, 0.0001)
        xt[l] ~ dunif(0, 100)
      }
      
      # Prior on multivariate normal std deviation
      tau.B.year.raw[1:L, 1:L] ~ dwish(W.year[ , ], df.year)
      df.year <- L + 1
      sigma.B.year.raw[1:L, 1:L] <- inverse(tau.B.year.raw[ , ])
      for(l in 1:L){
        for(l.prime in 1:L){
          rho.B.year[l, l.prime] <- sigma.B.year.raw[l, l.prime]/sqrt(sigma.B.year.raw[l, l]*sigma.B.year.raw[l.prime, l.prime])
        }
        sigma.B.year[l] <- abs(xt[l])*sqrt(sigma.B.year.raw[l, l])
      }
    }
    ", fill = TRUE)
sink()

variables.site <- c("Intercept-site",
                    "Air Temperature",
                    "Air Temp Lag1",
                    "Air Temp Lag2")
J <- length(unique(etS$site))
K <- length(variables.site)
n <- dim(etS)[1]
W.site <- diag(K)
X.site <- data.frame(int = 1, 
                airT = etS$airTemp, 
                airT1 = etS$airTempLagged1,
                airT2 = etS$airTempLagged2)

variables.year <- c("Intercept-year",
                    "dOY",
                    "dOY2",
                    "dOY3")
Ti <- length(unique(etS$year))
L <- length(variables.year)
W.year <- diag(L)
X.year <- data.frame(int=1, 
                dOY = etS$dOY, 
                dOY2 = etS$dOY^2,
                dOY3 = etS$dOY^3)


data <- list(n = n, 
             J = J, 
             K = K, 
             Ti = Ti,
             L = L,
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
            "B.site",
            "rho.B.site",
            "mu.site",
            "sigma.B.site",
            "B.year",
            "rho.B.year",
            "mu.year",
            "sigma.B.year")

#M1 <- bugs(etS, )

n.burn = 100
n.it = 500
n.thin = 1

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
plot(M3)
dev.off()

summary(M3)

summary(M3)$statistics[ , "Mean"]

rm(out)

pairs(as.matrix(M3[ , c(1:4, 17:20)]))

# compare with lme4
# Random slopes for crossed
lmer7 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 |site) + (dOY + I(dOY^2) + I(dOY^3)|year), data = etS)
summary(lmer7)
ranef(lmer7)


#################################
# remove xi then make mu.site and mu.year = 0
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

variables.site <- c("Intercept-site",
                    "Air Temperature",
                    "Air Temp Lag1",
                    "Air Temp Lag2")
J <- length(unique(etS$site))
K <- length(variables.site)
n <- dim(etS)[1]
W.site <- diag(K)
X.site <- data.frame(int = 1, 
                     airT = etS$airTemp, 
                     airT1 = etS$airTempLagged1,
                     airT2 = etS$airTempLagged2)

variables.fixed <- c("intercept", "flow")
K.0 <- length(variables.fixed)
X.0 <- data.frame(int = 1,
                  flow = etS$flow)

variables.year <- c("Intercept-year",
                    "dOY",
                    "dOY2",
                    "dOY3")
Ti <- length(unique(etS$year))
L <- length(variables.year)
W.year <- diag(L)
X.year <- data.frame(int=1, 
                     dOY = etS$dOY, 
                     dOY2 = etS$dOY^2,
                     dOY3 = etS$dOY^3)


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

n.burn = 1000
n.it = 1000
n.thin = 1

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
plot(M3)
dev.off()

summary(M3)

summary(M3)$statistics[ , "Mean"]

rm(out)

pairs(as.matrix(M3[ , c(1:8, 17:20)]))

#summary(M3)$statistics[ , "Mean"]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, K.0+k, 2))
names(fix.ef) <- c("Mean", "Std. Error")
row.names(fix.ef) <- c(variables.fixed, variables.site)
for(k in 1:K.0){
  fix.ef[k, ] <- summary(M3)$statistics[paste0('B.0[',k,']') , c("Mean", "SD")]
}
for(k in 1:K){
  fix.ef[k+K.0, ] <- summary(M3)$statistics[paste0('mu.site[',k,']') , c("Mean", "SD")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef <- as.data.frame(matrix(NA, k, 2))
names(ran.ef) <- c("Variance", "Std. Dev.")
row.names(ran.ef) <- variables.site
for(k in 1:K){
  ran.ef[k, 2] <- summary(M3)$statistics[paste0('sigma.b.site[',k,']') , c("Mean")]
  ran.ef[k, 1] <- ran.ef[k, 2] ^ 2
}
ran.ef

# Make Random Effects Output like summary(lmer)
ran.ef2 <- as.data.frame(matrix(NA, k, 2))
names(ran.ef2) <- c("Variance", "Std. Dev.")
row.names(ran.ef2) <- variables.year
for(k in 1:K){
  ran.ef2[k, 2] <- summary(M3)$statistics[paste0('sigma.b.year[',k,']') , c("Mean")]
  ran.ef2[k, 1] <- ran.ef2[k, 2] ^ 2
}
ran.ef2

# Make correlation matrix of random effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- summary(M3)$statistics[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.site <- round(cor.site, digits=3)
cor.site[upper.tri(cor.site, diag=TRUE)] <- ''
cor.site

# Make correlation matrix of random effects
cor.year <- as.data.frame(matrix(NA, K, K))
names(cor.year) <- variables.year
row.names(cor.year) <- variables.year
for(k in 1:K){
  for(k.prime in 1:K){
    cor.year[k, k.prime] <- summary(M3)$statistics[paste('rho.B.year[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.year <- round(cor.year, digits=3)
cor.year[upper.tri(cor.year, diag=TRUE)] <- ''
cor.year

sum.pred <- summary(M3)$statistics
pred.t <- as.data.frame(matrix(NA, n, 2))
for(i in 1:n){
  pred.t[i, 1] <- sum.pred[paste0('stream.mu[',i,']') , c("Mean")]
  print(i)
}
pred.t[ , 2] <- etS$site
pred.t$airTemp <- (etS$airTemp * sd(et2$airTemp)) + mean(et2$airTemp)
pred.t$year <- etS$year
pred.t$dOY <- etS$dOY
names(pred.t) <- c("streamTemp", "site", "airTemp", "year", "dOY")

ggplot(pred.t, aes(airTemp, streamTemp)) + geom_point(aes(colour = dOY)) + facet_grid(site ~ year)

ggplot(pred.t, aes(dOY, streamTemp)) + geom_point(size=0.75) + geom_point(data=etS, aes(dOY, temp), colour = "red", size=0.75) + facet_grid(site ~ year)

err <- pred.t$streamTemp - etS$temp
rmse(err)

# compare with lme4
# Random slopes for crossed
lmer7 <- lmer(temp ~ airTemp + airTempLagged1 + airTempLagged2 + flow + dOY + I(dOY^2) + I(dOY^3) + (airTemp + airTempLagged1 + airTempLagged2 | site) + (dOY + I(dOY^2) + I(dOY^3)|year), data = etS)
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

###############








###### Simulate Data #########
N <- 10000
k <- 4
x <- 1:N
f <- rep(rnorm(k, 0, 4), each = N/k)
e <- rnorm(N)
y <- x + f + e

fac <- gl(k, N/k)
library(lme4)
fm1 <- lmer(y ~ x + (1|fac))
summary(fm1)
