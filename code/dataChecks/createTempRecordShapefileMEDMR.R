rm(list = ls())

library(foreign)

load('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/MEDMR/streamTempSitesObservedClimateData_MEDMR.RData')

e <- masterData

e1 <- e[!is.na(e$temp),c('site', 'AgencyID', 'date', 'temp', 'airTemp', 'dayl', 'srad', 'swe', 'vp', 'prcp', 'Latitude', 'Longitude')]


load('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/MEDMR/covariateData_MEDMR.RData')

f <- covariateData


f1 <- f[,c('site', 'Forest', 'Agriculture', 'BasinElevationM', 'ReachSlopePCNT', 'TotDASqKM', 'WetlandOrWater', 'SurficialCoarseC', 'ImpoundmentsOpenSqKM')]

et <- merge(e1,f1, by = 'site', all.x = T, sort = F)

et <- et[,-1]

et <- replace(et, is.na(et), -9999)


names(et) <- c('site', 'date', 'streamTemp', 'airTemp', 'dayl', 'srad', 'swe', 'vp', 'prcp', 'Latitude', 'Longitude', 'Forest', 'Agro', 'BasnElevM', 'RchSlpPCT', 'TotDASqKM', 'WetOrWater', 'SurfCrseC', 'ImpOpnSqKM')

head(et)


write.dbf(et, file = 'C:/KPONEIL/temporary/tempRecordMEDMR.dbf')


