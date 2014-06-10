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

sink("code/correlatedSlopes.txt")
cat("
    model{
      # Likelihood
      for(i in 1:n){ # n observations
        temp[i] ~ dnorm(stream.mu[i], tau)
        stream.mu[i] <- inprod(B.site[site[i], ], X.site[i, ]) #+ inprod(B.year[year[i], X.year[i, ]]) # inprod(b.0[], X.0[i, ]) + 
      }
      
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
      for(k in 1:K){
        mu.site[k] <- xi[k]*mu.site.raw[k]
        mu.site.raw[k] ~ dnorm(0, 0.0001)
        xi[k] ~ dunif(0, 100)
      }
      
      # Prior on multivariate normal std deviation
      tau.B.site.raw[1:K, 1:K] ~ dwish(W[ , ], df.site)
      df.site <- K + 1
      sigma.B.site.raw[1:K, 1:K] <- inverse(tau.B.site.raw[ , ])
      for(k in 1:K){
        for(k.prime in 1:K){
          rho.B.site[k, k.prime] <- sigma.B.site.raw[k, k.prime]/sqrt(sigma.B.site.raw[k, k]*sigma.B.site.raw[k.prime, k.prime])
        }
        sigma.B.site[k] <- abs(xi[k])*sqrt(sigma.B.site.raw[k, k])
      }

    # 
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
W <- diag(K)
X.site <- data.frame(int=1, 
                airT=etS$airTemp, 
                airT1=etS$airTempLagged1,
                airT2=etS$airTempLagged2)
data <- list(n = n, 
             J = J, 
             K = K, 
             W = W,
             temp = etS$temp,
             X.site = as.matrix(X.site),
             site = etS$site)

inits <- function(){
  list(#B.raw = array(rnorm(J*K), c(J,K)), 
       #mu.site.raw = rnorm(K),
       sigma = runif(1),
       #tau.B.site.raw = rwish(K + 1, diag(K)),
       xi = runif(K))
}

params <- c("B.site",
            "rho.B.site",
            "mu.site",
            #"stream.mu",
            "sigma",
            "sigma.B.site")

#M1 <- bugs(etS, )

n.burn = 1000
n.it = 1000
n.thin = 1

library(parallel)
library(rjags)

CL <- makeCluster(3)
clusterExport(cl=CL, list("data", "inits", "params", "K", "J", "n", "W", "X.site", "n.burn", "n.it", "n.thin"), envir = environment())
clusterSetRNGStream(cl=CL, iseed = 2345642)

system.time(out <- clusterEvalQ(CL, {
  library(rjags)
  load.module('glm')
  jm <- jags.model("code/correlatedSlopes.txt", data, inits, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
}))

M1 <- mcmc.list(out)
stopCluster(CL)
pdf("/Users/Dan/Dropbox/correlatedSlopes.pdf")
plot(M1)
dev.off()

summary(M1)

rm(out)


