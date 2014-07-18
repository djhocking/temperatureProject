Retired: 7/17/14

#=========================================================================================================
# This function pulls the Daymet variables for a time series at one location nearest to the site location.
# It takes the stream temperature record and the string of variables to pull from Daymet.
#     At minimum, the record needs a "site", "Latitude", "Longitude", "year", and "dOY" columns.
# It returns the original dataframe with new columns for the Daymet variables.
#=========================================================================================================

indexLocalDaymetVariablesForObservedSites <- function(record, variables, daymetDirectory){
  
  library(ncdf)
  library(sp)
  
  start.time <- proc.time()[3]
  
  sites <- unique(record$site)
  
  for ( i in 1:length(sites)){
    
    print(paste0(round(i/length(sites), digits = 3)*100, '% done.'))
    
    # Select site
    curRecord <- record[which(record$site %in% sites[i]),]
    curRecord <- curRecord[which(curRecord$year < 2014),]
    
    # Site coordinates
    siteLon <- unique(curRecord$Longitude)
    siteLat <- unique(curRecord$Latitude)
    
    #Index the tile by site location:
    tile <- indexDaymetTileByLatLon(siteLat,siteLon)
    
    # Site time range
    begYr <- min(curRecord$year)
    
    # Loop through the variables and years in NetCDF files
    for (j in 1:length(variables)){
      
      for ( year in unique(curRecord$year) ){
        
        #Open the NetCDF with the known location:
        #--------------------------------------------
        NCDF <- open.ncdf(paste0(daymetDirectory, tile, '_', year,'/', variables[j], '.nc'))    #netcdf
        
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
        var = get.var.ncdf ( nc=NCDF, varid= paste0(variables[j]), start = start2, count = varcount )
        
        close.ncdf(NCDF)
        
        dOY <- dOY + 1  #Daymet doy starts at 0.
        
        if (year == begYr){
          coords <- cbind( as.vector(lon), as.vector(lat))
          
          distances <- spDistsN1(coords, c(siteLon, siteLat), longlat = TRUE)
          
          minDist <- min(distances)
          
          distpos <- which(distances == minDist)[1]
          
          coords <- as.data.frame(coords) # for indexing
          varLon  <- coords[distpos, 1]
          varLat  <- coords[distpos, 2]
          
          position <- which(lat == varLat & lon == varLon, arr.in = TRUE)  
        }
        
        tempVar <- data.frame(year, dOY, var[position[1], position[2], 1:365])
        names(tempVar) <- c("year", "dOY", paste0(variables[j]))
        
        if ( year == begYr ) ( mainVar <- tempVar)
        if ( year >  begYr ) ( mainVar <- rbind(mainVar, tempVar))
        
        rm(tempVar)
        
      }# End year loop
      
      # Join the variables into one dataframe
      if (j == 1) {allVars <- mainVar} else {allVars <- merge(allVars, mainVar, by = c('year','dOY'), all.x = T)}  
      
    }# End variable loop
    allVars$site <- sites[i]
    tempRecord <- merge(curRecord, allVars, by = c("site", "year", "dOY"), all.x = T, all.y = F, sort = F)
    
    if (i == 1) {fullRecord <- tempRecord} else {fullRecord <- rbind(fullRecord, tempRecord)}
  }
  
  # How long it takes to run
  end.time   <- proc.time()[3]
  print(paste0((end.time-start.time)/3600, " hours"))
  
  fullRecord$airTemp <- (fullRecord$tmin + fullRecord$tmax)/2
  
  return(fullRecord)
}
#=========================================================================================================
