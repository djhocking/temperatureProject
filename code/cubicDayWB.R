sink("code/cubicDayWB.txt")
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
    stream.mu[i] <- alpha + b.site[site[i]] + b.year[year[i]] + b.air*airTemp[i] + b.air1*airTempLagged1[i] + b.air2*airTempLagged2[i] + b.prcp*prcp[i] + b.prcp1*prcpLagged1[i] + b.prcp2*prcpLagged2[i]
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
  jm <- jags.model("code/cubicDayWB.txt", data1, inits1, n.adapt=n.burn, n.chains=1)
  fm <- coda.samples(jm, params1, n.iter = n.it, thin = n.thin)
  return(as.mcmc(fm))
})

out32 <- mcmc.list(out)
stopCluster(CL)

rm(out)

pdf(file="C:/Users/dhocking/Dropbox/cubicDayWB.pdf", width=10, height=10)
plot(out32[ , c("alpha",
                "b.air",
                "b.air1",
                "b.air2",
                "b.prcp",
                "b.prcp1",
                "b.prcp2",
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

