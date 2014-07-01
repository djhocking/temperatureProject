rm(list=ls())

library(ggplot2)
library(rjags)

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

n.burn = 5000
n.it = 5000
n.thin = 5

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
#plot(M3)
#dev.off()

Summary.M3 <- summary(M3)

Summary.M3$statistics[ , "Mean"]

rm(out)

pairs(as.matrix(M3[ , c(1:8, 17:20)]))

#summary(M3)$statistics[ , "Mean"]

# Make "Fixed Effects" Output like summary(lmer)
fix.ef <- as.data.frame(matrix(NA, K.0+k, 2))
names(fix.ef) <- c("Mean", "Std. Error")
row.names(fix.ef) <- c(variables.fixed, variables.site)
for(k in 1:K.0){
  fix.ef[k, ] <- Summary.M3$statistics[paste0('B.0[',k,']') , c("Mean", "SD")]
}
for(k in 1:K){
  fix.ef[k+K.0, ] <- Summary.M3$statistics[paste0('mu.site[',k,']') , c("Mean", "SD")]
}
fix.ef

# Make Random Effects Output like summary(lmer)
ran.ef <- as.data.frame(matrix(NA, k, 2))
names(ran.ef) <- c("Variance", "Std. Dev.")
row.names(ran.ef) <- variables.site
for(k in 1:K){
  ran.ef[k, 2] <- Summary.M3$statistics[paste0('sigma.b.site[',k,']') , c("Mean")]
  ran.ef[k, 1] <- ran.ef[k, 2] ^ 2
}
ran.ef

# Make Random Effects Output like summary(lmer)
ran.ef2 <- as.data.frame(matrix(NA, k, 2))
names(ran.ef2) <- c("Variance", "Std. Dev.")
row.names(ran.ef2) <- variables.year
for(k in 1:K){
  ran.ef2[k, 2] <- Summary.M3$statistics[paste0('sigma.b.year[',k,']') , c("Mean")]
  ran.ef2[k, 1] <- ran.ef2[k, 2] ^ 2
}
ran.ef2

# Make correlation matrix of random effects
cor.site <- as.data.frame(matrix(NA, K, K))
names(cor.site) <- variables.site
row.names(cor.site) <- variables.site
for(k in 1:K){
  for(k.prime in 1:K){
    cor.site[k, k.prime] <- Summary.M3$statistics[paste('rho.B.site[',k,',',k.prime,']', sep=""), "Mean"]
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
    cor.year[k, k.prime] <- Summary.M3$statistics[paste('rho.B.year[',k,',',k.prime,']', sep=""), "Mean"]
  }
}
cor.year <- round(cor.year, digits=3)
cor.year[upper.tri(cor.year, diag=TRUE)] <- ''
cor.year

sum.pred <- Summary.M3$statistics
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

pred.t$dOY.real <- etS$dOY*(sd(et2$dOY)) + mean(et2$dOY)
pred.t$temp <- etS$temp

# Are figures for presentations?
presentation <- T
if(presentation){
  theme_bw_present <- theme_set(theme_bw())
  theme_bw_present <- theme_update(axis.text=element_text(size=12), 
                                   axis.title=element_text(size=14,face="bold"))
}


ggplot(pred.t, aes(airTemp, streamTemp)) + geom_point(aes(colour = dOY.real)) + facet_grid(site ~ year)

# Observed vs. Predicted
ggplot(pred.t, aes(dOY.real, temp)) + geom_point(size=1, colour='black') + geom_point(aes(dOY.real, streamTemp), colour = 'red', size=0.75) + facet_grid(site ~ year) + ylab(label="Stream temperature (C)") + xlab("Day of the year")

############## Derived metrics ##########

# Mean maximum daily mean temperature by site (over years)
library(dplyr)
by.site <- group_by(pred.t, site)
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

######### Model fit and predictive ability ###########
err <- pred.t$streamTemp - etS$temp
rmse(err)

########### compare with lme4 ############
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
