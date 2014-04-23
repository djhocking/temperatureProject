#Author:
#  Kyle O'Neil
#Created:
#  12/15/13
#Last Updated:
#  01/15/14
#Language:
#	 R
#Description:
#  This code is a modification of "DayMetNetCDF.R" that calculates the spatial average of upstream climate
#   climate data for the sites of interest.
#========================================================================================================

rm(list=setdiff(ls(), c("Catchments")))

library(sp)
library(rgdal)
library(rgeos)
library(maptools)     # loads sp library too
library(chron)
library(ncdf)

baseDir <- 'C:/KPONEIL/GitHub/projects/'


#Load the function that indexes daymet tiles based on a lat/lon point:
source(paste0(baseDir, 'temperatureProject/code/functions/indexDaymetTileByLatLon.R'))
#==================================================================================================================
#                             Read in spatial data
#==================================================================================================================

proj4.NHD  <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#proj4.Daymet <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=lonlat +no_defs "
#proj4.NHDLambert <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83  +units=m +no_defs "

#Catchments <- readShapePoly ( "C:/KPONEIL/USGS/NHDPlusV2/Modified Versions/NENY_NHDCatchment_Lambert.shp", proj4string=CRS(proj4.NHDLambert))
Catchments <- readShapePoly ( "C:/KPONEIL/USGS/NHDPlusV2/Modified Versions/NENY_NHDCatchment.shp", proj4string=CRS(proj4.NHD))

dataSource <- "VTFWS"

load(paste0(baseDir, 'temperatureProject/dataIn/delineatedCatchments/DelineatedCatchments_NHDPlus_NENY.RData'))
DelineatedCatchmentsMaster <- NENYDelineatedCatchments
MasterLength <- length(DelineatedCatchmentsMaster)

#==================================================================================================================
#                             Read in Site data
#==================================================================================================================

setwd(paste0(baseDir, 'temperatureProject/dataIn/', dataSource))

load(paste0("streamTempData_",   dataSource, ".RData"))
load(paste0("covariateData_", dataSource, ".RData"))

Sites <- unique(masterData$site)

#Daymet variables you want:
Variables <- c("prcp")  #, "dayl", "srad", "swe", "tmax", "tmin", "vp")

stats <- data.frame(site = NA, SiteLat = NA, SiteLon = NA, VarLat = NA, VarLon = NA, MinDist = NA)

#==================================================================================================================
#                               Get a master list of Daymet coordinates
#==================================================================================================================

Tiles <- c(11754, 11755, 11756, 11934, 11935, 11936, 12114, 12115, 12116, 12117, 12295, 12296, 12297)
numTiles <- length(Tiles)

#Create a master list of the daymet coords.
#------------------------------------------
for (i in 1:numTiles){
  print(i)
  
  #Now open the NetCDF with the known location:
  #--------------------------------------------
  NCDF <- open.ncdf(paste0("F:/KPONEIL/SourceData/Projected/DAYMET/Daily/", Tiles[i], "_2010/prcp.nc"))    #netcdf   
         
  #Dimension limits of each of the variables we'll use:
  #----------------------------------------------------
  start1 = c(1,1)
  latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
  loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
 
  #Read in variables:
  #------------------
  lat = get.var.ncdf ( nc=NCDF, varid="lat",                 start = start1, count = latcount )
  lon = get.var.ncdf ( nc=NCDF, varid="lon",                 start = start1, count = loncount )
  
  close.ncdf(NCDF)
  TempCoords <- cbind( as.vector(lon), as.vector(lat))
  colnames(TempCoords) <- c("Longitude", "Latitude")
  
  if (i ==1) {MasterCoords <- TempCoords}
  if (i > 1) {MasterCoords <- rbind(MasterCoords, TempCoords)}
}

MasterCoordsMatrix <- MasterCoords
MasterCoords <- as.data.frame(MasterCoords)

#==================================================================================================================
#                               Loop through Sites and NetCDFs, getting data.
#==================================================================================================================

start.time <- proc.time()[3]
for ( i in 1:length(Sites)){
    
  print(i)
  
  #----------------------------------------------------------------------------------
  #Get the catchment polygon:
  #----------------------------------------------------------------------------------
  featureID <- covariateData$FEATUREID[covariateData$site %in% Sites[i]]
  features <- DelineatedCatchmentsMaster[[which(sapply(c(1:MasterLength),FUN = function(x){DelineatedCatchmentsMaster[[x]][1]==featureID})==TRUE)]]
  
  CatchmentShape <- Catchments[Catchments$FEATUREID %in% features,]
  BasinShape     <- gUnaryUnion(CatchmentShape) #dissolve individual catchments
  a <- SpatialPoints(MasterCoords, proj4string=CRS(proj4.NHD))
  
  inside <- as.data.frame(a[!is.na(over(a, BasinShape)),])
  #----------------------------------------------------------------------------------
  
  
  #----------------------------------------------------------------------------------
  #If no point falls within the catchment, find the nearest one.
  #----------------------------------------------------------------------------------
  if(nrow(inside) == 0 ){
    
    TempLat <- covariateData$Latitude [covariateData$site %in% Sites[i]]
    TempLon <- covariateData$Longitude[covariateData$site %in% Sites[i]]
    
    distances <- spDistsN1(MasterCoordsMatrix, c(TempLon, TempLat), longlat = TRUE)
    MinDist <- min(distances)
    distpos <- which(distances == MinDist)[1]
    
    NearLon  <- MasterCoords[distpos, 1]
    NearLat  <- MasterCoords[distpos, 2]
    
    inside[1,1] <- NearLon
    inside[1,2] <- NearLat
  }
  #----------------------------------------------------------------------------------
  
  
  #----------------------------------------------------------------------------------
  #Pull the Tiles for the points within the catchment:
  #----------------------------------------------------------------------------------
  for(k in 1:length(inside[,1])){
    
    SiteLon <- inside[k,1]
    SiteLat <- inside[k,2]
    
    #Index the tile by site location:
    Tile <- indexDaymetTileByLatLon(SiteLat,SiteLon)
    
    temp <- data.frame('Longitude' = SiteLon, 'Latitude' = SiteLat, 'Tile' = Tile)
    
    if (k ==1) SpatialLocs <- temp
    if (k > 1) SpatialLocs <- rbind(SpatialLocs, temp)
  }
  rm(SiteLat, SiteLon)
  
  SpatialLocs <- SpatialLocs[ order(SpatialLocs$Tile), ]
  SubTiles <- unique(SpatialLocs$Tile)
  #----------------------------------------------------------------------------------
  
  
  #----------------------------------------------------------------------------------
  #Link record with catchment area
  #----------------------------------------------------------------------------------
  CurRecord <- masterData[which(masterData$site == Sites[i]),]
  CurRecord <- CurRecord[which(CurRecord$year < 2013),]
  
  BegYr <- min(CurRecord$year)
  EndYr <- max(CurRecord$year)
    
  for (j in 1:length(Variables)){
    
    for ( year in BegYr:EndYr ){
      
      for ( t in 1:length(SubTiles)){
  
        NCDF <- open.ncdf(paste0("F:/KPONEIL/SourceData/Projected/DAYMET/Daily/", SubTiles[t], "_", year,"/", Variables[j], ".nc"))    #netcdf
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
      
        TileCoords <- as.data.frame(cbind( as.vector(lon), as.vector(lat)))
        names(TileCoords) <- c('Lon', 'Lat')
        
        xx <- SpatialLocs[which(SpatialLocs$Tile == SubTiles[t]),]
        
        for (m in 1:length(xx[,1])){
          position <- which(lon == xx$Longitude[m] & lat == xx$Latitude[m], arr.in = TRUE)
          
          if ( t == 1 & m == 1) {temp.var <- data.frame(year, dOY, var[position[1], position[2], 1:365])} else(temp.var <- cbind(temp.var, var[position[1], position[2], 1:365]))
          
        }
      } 
      
      ifelse( ncol(temp.var) > 3, R <- rowMeans(temp.var[,-c(1,2)], na.rm = FALSE, dims = 1),  R <- temp.var[,-c(1,2)] )
      
      temp.var <- data.frame(temp.var[,c(1,2)], R)
      names(temp.var) <- c("year", "dOY", paste0(Variables[j]))
      
      if ( year == BegYr ) ( main.var <- temp.var)
      if ( year >  BegYr ) ( main.var <- rbind(main.var, temp.var))

      rm(temp.var, R)
    }

    if (j == 1) {all.vars <- main.var} else {all.vars <- merge(all.vars, main.var, by = c('year','dOY'), all.x = T)}  
  }

  #Add data into main dataframe
  all.vars$site <- Sites[i]
  TempRecord <- merge(CurRecord, all.vars, by = c("site", "year", "dOY"), all.x = T)
  if (i == 1) {FullRecord <- TempRecord} else {FullRecord <- rbind(FullRecord, TempRecord)}
}
end.time   <- proc.time()[3]
print(paste0((end.time-start.time)/3600, " hours"))



#Merge data in with other daymet and stream temp dataframes:
#===========================================================

MergeCols <- c('site', 'year', 'dOY', 'prcp')

UpstreamVariables <- FullRecord[,names(FullRecord) %in% MergeCols]

load(paste0(baseDir, 'temperatureProject/dataIn/', dataSource, '/streamTempSitesObservedClimateData_',dataSource, '_NeedPrcp.RData'))

masterData <- merge(masterData, UpstreamVariables, by = c('site', 'year', 'dOY'), all.x = T)

#Check for reasonable NA count. (Don't include stream temperature column.)
length(which(is.na(masterData[,-7])))

save(masterData, file = paste0(baseDir, 'temperatureProject/dataIn/', dataSource, '/streamTempSitesObservedClimateData_',dataSource, '.RData'))




