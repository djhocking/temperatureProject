library(plyr)

rm(list=ls())


agency <- c('CTDEP', 'MAFW', 'MAUSGS', 'MEDMR', 'NHDES', 'NHFG', 'USFS', 'VTFWS')

for (i in 1:length(agency)){
  
  load(paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', agency[i], '/streamTempSitesObservedClimateData_', agency[i],'.RData'))
  
  e <- masterData
  e$siteYear <- paste(e$site,e$year,sep='_')
  
  d <- ddply( e, .(siteYear), summarise, All = all(is.na(temp)))
  d1 <- d$siteYear[d$All == TRUE]
  
  e1 <- e[!e$siteYear %in% d1,]
  
  masterData <- e1[,-which(colnames(e1) == 'siteYear')]
  
  save(masterData, file = paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', agency[i], '/streamTempSitesObservedClimateData_', agency[i],'.RData'))
  
  print(paste0( ' Agency # ', i, '.   Removed ', length(which(is.na(e$temp))) - length(which(is.na(e1$temp))), ' records of NA.'))

}
