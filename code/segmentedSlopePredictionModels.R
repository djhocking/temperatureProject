
#---------------------------------------------------------------------

#                            OTHER MODELS

#---------------------------------------------------------------------




#Model 3: Selected interactions
#---------------------------------------------------------------------
m3 <- 'temp~airTemp + airTempLagged1 + airTempLagged2 + LatitudeS + LongitudeS + BasinElevationMS + ReachSlopePCNTS + 
          TotDASqKMLS + SurficialCoarseCLS + ForestS + AgricultureLS + WetlandOrWaterS + ImpoundmentsOpenSqKMLS +
          daylS + sradS + sweLS + 
          airTemp*ReachSlopePCNTS + airTemp*TotDASqKMLS + airTemp*WetlandOrWaterS + 
          airTemp*ImpoundmentsOpenSqKMLS + airTemp*sweLS + 
          airTempLagged1 + airTempLagged1*ReachSlopePCNTS + airTempLagged1*TotDASqKMLS + airTempLagged1*WetlandOrWaterS + 
          airTempLagged1*ImpoundmentsOpenSqKMLS + airTempLagged1*sweLS + 
          airTempLagged2 + airTempLagged2*ReachSlopePCNTS + airTempLagged2*TotDASqKMLS + airTempLagged2*WetlandOrWaterS + 
          airTempLagged2*ImpoundmentsOpenSqKMLS + airTempLagged2*sweLS + 
          ReachSlopePCNTS*ImpoundmentsOpenSqKMLS + ReachSlopePCNTS*daylS + ReachSlopePCNTS*sradS +
          TotDASqKMLS*ImpoundmentsOpenSqKMLS + TotDASqKMLS*daylS + TotDASqKMLS*sradS +
          SurficialCoarseCLS*ForestS +   
          ForestS*daylS + ForestS*sradS + ForestS*sweLS +
          WetlandOrWaterS*daylS + WetlandOrWaterS*sradS + 
          ImpoundmentsOpenSqKMLS*daylS + ImpoundmentsOpenSqKMLS*sradS + 
          daylS*sweLS + 
          sradS*sweLS'

#Model 4: Two-way interactions  without Daymet data.
#---------------------------------------------------------------------
m4 <- 'temp~(airTemp+airTempLagged1+airTempLagged2+
          #segment+
          LatitudeS+LongitudeS+
          ForestS+ AgricultureLS+ 
          BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMLS+ 
          WetlandOrWaterS+ SurficialCoarseCLS+ImpoundmentsOpenSqKMLS 
          )^2'



#Model 7: Selected interactions with prcp
#---------------------------------------------------------------------
m7 <- 'temp~airTemp + airTempLagged1 + airTempLagged2 + LatitudeS + LongitudeS + BasinElevationMS + ReachSlopePCNTS + 
          TotDASqKMLS + SurficialCoarseCLS + ForestS + AgricultureLS + WetlandOrWaterS + ImpoundmentsOpenSqKMLS +
          daylS + sradS + sweLS + prcpLS +
          airTemp*ReachSlopePCNTS + airTemp*TotDASqKMLS + airTemp*WetlandOrWaterS + 
          airTemp*ImpoundmentsOpenSqKMLS + airTemp*sweLS + 
          airTempLagged1 + airTempLagged1*ReachSlopePCNTS + airTempLagged1*TotDASqKMLS + airTempLagged1*WetlandOrWaterS + 
          airTempLagged1*ImpoundmentsOpenSqKMLS + airTempLagged1*sweLS + 
          airTempLagged2 + airTempLagged2*ReachSlopePCNTS + airTempLagged2*TotDASqKMLS + airTempLagged2*WetlandOrWaterS + 
          airTempLagged2*ImpoundmentsOpenSqKMLS + airTempLagged2*sweLS + 
          ReachSlopePCNTS*ImpoundmentsOpenSqKMLS + ReachSlopePCNTS*daylS + ReachSlopePCNTS*sradS +
          TotDASqKMLS*ImpoundmentsOpenSqKMLS + TotDASqKMLS*daylS + TotDASqKMLS*sradS +
          SurficialCoarseCLS*ForestS +   
          ForestS*daylS + ForestS*sradS + ForestS*sweLS +
          WetlandOrWaterS*daylS + WetlandOrWaterS*sradS + 
          ImpoundmentsOpenSqKMLS*daylS + ImpoundmentsOpenSqKMLS*sradS + 
          daylS*sweLS + 
          sradS*sweLS +
          prcpLS*BasinElevationMS'
