library(foreign)


DaymetTiles <- c(11754, 11755, 11756, 11934, 11935, 11936, 12114, 12115, 12116, 12117, 12295, 12296, 12297)
Year <- 2010

for ( i in 1:length(DaymetTiles)){
  
  print(i)
  setwd("C:/KPONEIL/USGS/Stream Temperature/data/temperature/fromKyle/BP_Analysis/BP_Analysis/DaymetClimateData")
  
  load(paste0('NHD_DaymetTile_' , DaymetTiles[i], '_', Year, '.RData'))
  
  print(paste0('NHD_DaymetTile_' , DaymetTiles[i], '_', Year, '.RData'))
  
  test <- data.frame(unique(FullRecord$FEATUREID), 1)
  
  if ( i == 1 ) {TestPreds <- test}  else (TestPreds <- rbind(TestPreds, test))
  
}

names(TestPreds) <- c('FEATUREID', 'test')
write.dbf(TestPreds, file = ('C:/KPONEIL/USGS/Stream Temperature/data/temperature/fromKyle/BP_Predictions/TESTALL.dbf'))


length(unique(TestPreds$FEATUREID))

length(TestPreds$FEATUREID)






#===================================
#Creating a shapefile of the Points:
#===================================
Tiles <- c(11754, 11755, 11756, 11934, 11935, 11936, 12114, 12115, 12116, 12117, 12295, 12296, 12297)


for (t in 1:length(Tiles)){
  
  NCDF <- open.ncdf(paste0("C:/KPONEIL/SourceData/Projected/DAYMET/Daily/", Tiles[t], "_2000/prcp.nc"))    #netcdf
  
  #Dimension limits of each of the variables we'll use:
  start1 = c(1,1)
  latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
  loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
  
  
  start2 = c(1, 1, 1)
  varcount = c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2], NCDF$var$yearday$varsize)
  
  #Read in variables:
  #------------------
  lat = get.var.ncdf ( nc=NCDF, varid="lat",                 start = start1, count = latcount )
  lon = get.var.ncdf ( nc=NCDF, varid="lon",                 start = start1, count = loncount )
  dOY = get.var.ncdf ( nc=NCDF, varid="yearday",             start = 1,      count = YDcount  )
  var = get.var.ncdf ( nc=NCDF, varid= paste0(Variables[j]), start = start2, count = varcount )
  
  
  close.ncdf(NCDF)
  
  DataPres<- length(which(!is.na(var[,,1:365])))/365
  
  templist <- cbind( as.vector(lon), as.vector(lat), Tiles[t], as.vector(var[,,DataPres]))
  
  ifelse (t == 1, coords <- templist, coords <- rbind(coords, templist))
  
  rm(lat, lon, templist)
}

colnames(coords) <- c("Lon", "Lat", "Tile")

distances <- spDistsN1(coords[,1:2], c(SiteLon, SiteLat), longlat = TRUE)

coords <- as.data.frame(coords) # for indexing

