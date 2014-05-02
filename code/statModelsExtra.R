

#Model 0: Air temp & lags only
#---------------------------------------------------------------------
m0 <- temp~airTemp+airTempLagged1+airTempLagged2

#Model 1: All main effects
#---------------------------------------------------------------------
m1 <- temp~airTemp+airTempLagged1+airTempLagged2+
  LatitudeS+LongitudeS+
  ForestS+ AgricultureLS+ 
  BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ 
  WetlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ 
  daylS + sradS + sweLS


m0S2 <- lm(m0, data=et2)
#-----------------------
if(validate){
  vm0S2 <- validateModel( m0S2, et2, valData2 )
  if(createBiasMaps) {makeBiasMap(vm0S2)}
}

m1S2 <- lm(m1, data=et2)
#-----------------------
if(validate){
  vm1S2 <- validateModel( m1S2, et2, valData2 )
  if(createBiasMaps) { makeBiasMap(vm1S2) }
}


AIC(m0S2)
AIC(m0S2,m1S2)$AIC




m0a <- paste0("temp~airTemp+airTempLagged1+airTempLagged2")
m1a <- paste0("temp~airTemp+airTempLagged1+airTempLagged2+
  LatitudeS+LongitudeS+
  ForestS+ AgricultureLS+ 
  BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ 
  WetlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ 
  daylS + sradS + sweLS")

segModels <- list(m0a, m1a)

mS2  <- list()
vmS2 <- list()

#Make this a funciton
for (i in 1:length(segModels)){
  
  mS2[[i]] <- lm(as.formula(segModels[[i]]), data=et2)  #tried to make this into a function, but the "validateModel" function can't handle it.
  
  if (validate){
    vmS2[[i]] <- validateModel( mS2[[i]], et2, valData2 )
    if(createBiasMaps) {makeBiasMap(vmS2[[i]])}
  }
}



AIC(mS2[[i]])#, m1S2)$df
AIC(m0S2,m1S2)$AIC




































m1S2 <- lm(m1, data=et2, na.action="na.exclude")

rez <- resid(m1S2)


plot(et2$prcp, rez)


hist(et2$prcp[which(et2$prcp == 0)], rez[which(et2$prcp == 0)])



# Predict stream temperatures:
et2[,c('pred', 'lowr', 'upr')] 


temp <- predict (m1S2, newdata=et2, interval = "confidence") #default confidence interval (CI) = 95%