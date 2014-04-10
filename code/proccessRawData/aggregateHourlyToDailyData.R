#Steps for prepping Access database for R code that aggregates hourly to daily data:
#  1) Copy/Paste record from Access into Excel
#  2) Remove spaces within cells
#  3) Date:
#     Format cells -> custom: Yyyy-mm-dd
#  4) Time:
#     Format cells -> time: 24-hr time
#  5) Save as .CSV file

rm(list=ls())

setwd('C:/KPONEIL/SourceData/Raw/streamTemperature/VT')

TEMP1      <- read.csv("VTFWSdata1.csv", header=T)
TEMP1$temp <- as.numeric(as.character(TEMP1$temp)) #Because TEMP gets read in as factor

TEMP2      <- read.csv("VTFWSdata2.csv", header=T)
TEMP2$temp <- as.numeric(as.character(TEMP2$temp))

TEMP       <- rbind(TEMP1,TEMP2)
  #colnames(TEMP) <- c('site', 'date', 'time', 'temp')

sites      <- read.csv('VTFWSsites.csv', header = T)
  #colnames(sites) <- c('site', 'stream', 'Latitude', 'Longitude')

TEMP <- TEMP[TEMP$site %in% sites$site,]

temp.date <- as.Date(TEMP$date)
#temp.temp <- TEMP$temp
temp.temp <- (TEMP$temp - 32)*(5/9)  # <---- Convert from Fahrenheit to Celsius

temp.id <- TEMP$site

#Number of sites
uni <- unique(TEMP$site)

#Mean and max:
for (i in 1:length(uni)) {
  print(i/length(uni))
  
  a <- which(temp.id == uni[i])
  uni.date <- unique(temp.date[a])
  temp.1 <- cbind(temp.date[a],temp.temp[a])
  temp.2 <- tapply(temp.1[,2],temp.1[,1],mean)
  
  out <- data.frame(uni[i], uni.date, temp.2[])
  colnames(out) <- c('site', 'date', 'temp')
  
  if( i == 1) { finalData <- out} else ( finalData <- rbind(finalData, out) )
}

finalData$year <- as.numeric(strftime(finalData$date, '%Y'))
finalData$dOY <- as.numeric(strftime(finalData$date, '%j'))




#Loop through all sites and fill gaps with NAs
Sites <- unique(finalData$site)
tempFrame <- finalData

for ( i in 1:length(Sites) ){
  
  expData <- tempFrame[tempFrame$site == Sites[i],]
  
  start.date <- min(expData$date)
  end.date   <- max(expData$date)  
  
  Record <- data.frame(seq(from=as.Date(start.date),to=as.Date(end.date),by="day"))
  names(Record) <- "date"
  Record$year <- as.numeric(strftime(Record$date, '%Y'))
  Record$dOY <- as.numeric(strftime(Record$date, '%j'))
  
  newRecord <- merge(Record, expData, by = c('date', 'dOY', 'year'), all.x = T, all.y = F, sort = F)
  
  #Fill in blanks:
  newRecord$site       <- expData$site[1]
  newRecord$AgencyID   <- expData$AgencyID[1]
  newRecord$agency     <- expData$agency[1]
  #newRecord$Latitude   <- expData$Latitude[1]
  #newRecord$Longitude  <- expData$Longitude[1]
  
  if(i == 1) { finalData <- newRecord} else ( finalData <- rbind(finalData, newRecord))
  
}


masterData <- merge(finalData, sites, by = 'site', all.x = T, sort = F)

save(masterData,file="dailyStreamTempVTFWS.Rdata")
