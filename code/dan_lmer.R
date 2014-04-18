

ggplot(data = et[which(et$site == "CTDEP_1"), ], aes(date.y, airTemp)) +
  geom_line()


m1 <- temp~airTemp+airTempLagged1+airTempLagged2+
  LatitudeS+LongitudeS+
  ForestS+ AgricultureLS+ 
  BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ 
  WetlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ 
  daylS + sradS + sweLS

m1S2 <- lm(m1, data=et2)
summary(m1S2)

lmm1 <- temp~airTemp+airTempLagged1+airTempLagged2+
  LatitudeS+LongitudeS+
  ForestS+ AgricultureLS+ 
  BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ 
  WetlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS+ 
  daylS + sradS + sweLS + (1|site)

et2$fsite <- as.factor(et2$site)
lmm1S2 <- lmer(lmm1, data = et2)
summary(lmm1S2)

cbind(coef(m1S2), fixef(lmm1S2)) # shows that with lm sites with lots of data are exaggerating the effect of lat and lon and other variables that are constant over the time of the data (through pseudorelication)
