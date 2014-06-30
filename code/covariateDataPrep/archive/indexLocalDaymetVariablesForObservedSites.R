#Author:
#  Kyle O'Neil
#Created:
#  12/15/13
#Last Updated:
#  01/10/14
#Language:
#	 R
#Description:
#  This code pulls daily Daymet climate data from a group of NetCDF files for a given list of sites.
#  This version works off of an existing dataframe of stream temperature data to index climate data, 
#  but the code can easily be modified to index any location/time period in the current range of
#  NetCDF files which is VT, NH, MA, CT, and RI over 1980-2012.
#========================================================================================================

####Goes through a list of basins (GAGES II) and aggregate met data
rm(list=ls())
library(sp)
library(rgdal)
library(rgeos)
library(maptools)     # loads sp library too
library(chron)
library(ncdf)

#Load the function that indexes daymet tiles based on a lat/lon point:
source("C:/KPONEIL/GitHub/projects/temperatureProject/code/functions/indexDaymetTileByLatLon.R")

#==================================================================================================================
#                             Read in Site data
#==================================================================================================================

Source <- "MADEP"

setwd(paste0("C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/", Source))

load(paste0("dailyStreamTemp", Source, ".RData"))

masterData$agency <- Source
masterData$AgencyID <- masterData$site
masterData$site <- paste0(Source, '_', masterData$AgencyID)


#Write out a .csv file for the creation of a shapefile:
#------------------------------------------------------
load(paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', Source, '/streamTempData_', Source, '.RData'))
shapefile <- masterData[, c('site', 'Latitude', 'Longitude')]
write.csv(shapefile, file = paste0('C:/KPONEIL/USGS/Stream Temperature/Shapefiles/siteLocations/', Source, '.csv'))


#load("NewCovariateData_865sites.RData")
#load("NewCovariateData_MEDMRsites.RData")

#Daymet variables you want:
Variables <- c("dayl", "srad", "swe", "tmax", "tmin", "vp") #"prcp", 

stats <- data.frame(site = NA, SiteLat = NA, SiteLon = NA, VarLat = NA, VarLon = NA, MinDist = NA)

#==================================================================================================================
#                             Loop through NetCDF daily data files
#==================================================================================================================

Sites <- unique(masterData$site)

start.time <- proc.time()[3]
for ( i in 1:length(Sites)){

  print(i)
 
  CurRecord <- masterData[which(masterData$site %in% Sites[i]),]
  CurRecord <- CurRecord[which(CurRecord$year < 2014),]
  
  SiteLon <- unique(CurRecord$Longitude)#covariate.data$Longitude[which(covariate.data$site == Sites[i])]
  SiteLat <- unique(CurRecord$Latitude)#covariate.data$Latitude [which(covariate.data$site == Sites[i])]
  
  #Index the tile by site location:
  Tile <- indexDaymetTileByLatLon(SiteLat,SiteLon)
  
  BegYr <- min(CurRecord$year)
  EndYr <- max(CurRecord$year)
    
  for (j in 1:length(Variables)){
    
    for ( year in BegYr:EndYr ){
            
      #Now open the NetCDF with the known location:
      #--------------------------------------------
      NCDF <- open.ncdf(paste0("F:/KPONEIL/SourceData/climate/DAYMET/unzipped/Daily/", Tile, "_", year,"/", Variables[j], ".nc"))    #netcdf
      #Dimension limits of each of the variables we'll use:
      #----------------------------------------------------
      start1 = c(1,1)
      latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
      loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
      YDcount  <- NCDF$var$yearday$varsize      
      
      start2 = c(1, 1, 1)
      varcount = c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2], NCDF$var$yearday$varsize)

      #Read in variables:
      #------------------
      lat = get.var.ncdf ( nc=NCDF, varid="lat",                 start = start1, count = latcount )
      lon = get.var.ncdf ( nc=NCDF, varid="lon",                 start = start1, count = loncount )
      dOY = get.var.ncdf ( nc=NCDF, varid="yearday",             start = 1,      count = YDcount  )
      var = get.var.ncdf ( nc=NCDF, varid= paste0(Variables[j]), start = start2, count = varcount )
      
      close.ncdf(NCDF)

      dOY <- dOY + 1  #Daymet doy starts at 0.
      
      if (year == BegYr){
        coords <- cbind( as.vector(lon), as.vector(lat))
      
        distances <- spDistsN1(coords, c(SiteLon, SiteLat), longlat = TRUE)
  
        MinDist <- min(distances)
  
        distpos <- which(distances == MinDist)[1]
      
        coords <- as.data.frame(coords) # for indexing
        VarLon  <- coords[distpos, 1]
        VarLat  <- coords[distpos, 2]

        position <- which(lat == VarLat & lon == VarLon, arr.in = TRUE)  
      }
      
      temp.var <- data.frame(year, dOY, var[position[1], position[2], 1:365])
      names(temp.var) <- c("year", "dOY", paste0(Variables[j]))
      
      if ( year == BegYr ) ( main.var <- temp.var)
      if ( year >  BegYr ) ( main.var <- rbind(main.var, temp.var))

      rm(temp.var)
      
      newstats <- data.frame(Sites[i], SiteLat, SiteLon, VarLat, VarLon, MinDist)
      names(newstats)[1] <- "site"
      stats <- rbind(stats, newstats)
    }

    if (j == 1) {all.vars <- main.var} else {all.vars <- merge(all.vars, main.var, by = c('year','dOY'), all.x = T)}  
    
  }
  all.vars$site <- Sites[i]
  TempRecord <- merge(CurRecord, all.vars, by = c("site", "year", "dOY"), all.x = T)

  if (i == 1) {FullRecord <- TempRecord} else {FullRecord <- rbind(FullRecord, TempRecord)}
}

end.time   <- proc.time()[3]
print(paste0((end.time-start.time)/3600, " hours"))

FullRecord$airTemp <- (FullRecord$tmin + FullRecord$tmax)/2


#Save paired stream and air temp:
#--------------------------------
masterData <- FullRecord[,c('site', 'year', 'dOY', 'date', 'AgencyID', 'agency', 'temp', 'airTemp', 'Latitude', 'Longitude')]
save(masterData, file = paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', Source, '/streamTempData_', Source, '.RData'))


#Save all climate data:
#----------------------
masterData <- FullRecord[order(FullRecord$site,FullRecord$year,FullRecord$dOY),]
save(masterData, file = paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', Source, '/streamTempSitesObservedClimateData_', Source, '_NeedPrcp.RData'))






#
#
# Stop here.
#
#















#Look at a couple of metrics:

length(which(is.na(FullRecord)))/length(which(!is.na(FullRecord)))*100
#length(which(is.na(FullRecord$prcp)))/length(FullRecord$prcp)*100





#Check to make sure things ran properly before running next piece

#master.data <- FullRecord
#save(FullRecord, stats, file = "C:/KPONEIL/USGS/Stream Temperature/StreamTempsWithDayMet_865.RData")
#save(master.data, file = "C:/KPONEIL/USGS/Stream Temperature/TempsWithDayMet_865FINAL.RData")



#write.csv(FullRecord, file = "C:/KPONEIL/USGS/Stream Temperature/test.csv")


