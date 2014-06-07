################# JAGS #####################

#plot(etS$airTemp, etS$temp)
#points(etS$airTemp[10111], etS$temp[10111], col='red')

# To have more than 1 randomly varying slope it is best to use the scaled inverse-Wishart model (Gelman and Hill p377) - only for nested effects or correlated intercepts and slopes?


etS <- etS[!is.na(etS$flow),]

# Non-nested model with random intercepts
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


# Non-nested model with random intercepts and slopes
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


# Non-nested model with random intercepts and slopes
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
    
    sigma.air ~ dunif(0, 100)
    tau.air <- 1 / (sigma.air * sigma.air)
    
    sigma.air1 ~ dunif(0, 100)
    tau.air1 <- 1 / (sigma.air1 * sigma.air1)
    
    sigma.air2 ~ dunif(0, 100)
    tau.air2 <- 1 / (sigma.air2 * sigma.air2)
    
    sigma.year ~ dunif(0, 100)
    tau.year <- 1 / (sigma.year * sigma.year)
    
    sigma.dOY ~ dunif(0, 100)
    tau.dOY <- 1 / (sigma.dOY * sigma.dOY)
    
    sigma.dOY2 ~ dunif(0, 100)
    tau.dOY2 <- 1 / (sigma.dOY2 * sigma.dOY2)
    
    sigma.dOY3 ~ dunif(0, 100)
    tau.dOY3 <- 1 / (sigma.dOY3 * sigma.dOY3)
    
    for(j in 1:n.sites) {
    b.site[j] ~ dnorm(0, tau.site)
    b.air[j] ~ dnorm(0, tau.air)
    b.air1[j] ~ dnorm(0, tau.air1)
    b.air2[j] ~ dnorm(0, tau.air2)
    }
    
    for(t in 1:n.years) {
    b.year[t] ~ dnorm(0, tau.year)
    b.dOY[t] ~ dnorm(0, tau.dOY)
    b.dOY2[t] ~ dnorm(0, tau.dOY2)
    b.dOY3[t] ~ dnorm(0, tau.dOY3)
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
                 "sigma.dOY",
                 "sigma.dOY3",
                 "sigma.dOY3",
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


