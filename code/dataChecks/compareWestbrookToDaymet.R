rm(list = ls())

library(ggplot2)
library(gridExtra)
library(reshape2)

# Prep the Smith data
# ===================
#Met Station Catchment: 6780769

load('F:/KPONEIL/SourceData/climate/UMASS (Westbrook)/WestbrookMetData.RData')

smith <- wbMetData

smith$year <- as.numeric(strftime(smith$date, '%Y'))
smith$dOY <- as.numeric(strftime(smith$date, '%j'))

smith <- smith[smith$year %in% c(2011,2012),c('year', 'dOY', 'date', 'meanAirTemp', 'srad', 'prcp')]
names(smith) <- c('year', 'dOY', 'date', 'airTempSmith', 'sradSmith', 'prcpSmith')

smith <- replace(smith, smith == 99999, NA)


# Prep the DayMet data:
# =====================

load('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/westbrookAirTemp/WestbrookDaymet2011.RData')
FR11 <- FullRecord

load('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/westbrookAirTemp/WestbrookDaymet2012.RData')
FR12 <- FullRecord


FR11 <- FR11[,c('FEATUREID', 'year', 'dOY', 'srad', 'airTemp', 'prcp')]
FR12 <- FR12[,c('FEATUREID', 'year', 'dOY', 'srad', 'airTemp', 'prcp')]

names(FR11) <- c('FEATUREID', 'year', 'dOY', 'sradDay', 'airTempDay', 'prcpDay')
names(FR12) <- c('FEATUREID', 'year', 'dOY', 'sradDay', 'airTempDay', 'prcpDay')



# Make plots:
# ===========
features <- unique(FR12$FEATUREID)

for( i in 1:length(features)){
  
  # Index individual catchments
  Day <- rbind ( FR11[FR11$FEATUREID == features[i],],FR12[FR12$FEATUREID == features[i],] )

  compare <- merge(smith, Day, by = c('year', 'dOY'), all.x = T, sort = F)

  compare$season <- ifelse( compare$dOY<80 ,'winter',
                    ifelse( compare$dOY<172,'spring',
                    ifelse( compare$dOY<264,'summer', 
                    ifelse( compare$dOY<355,'fall'  ,'winter'  ))))
  compare$season <- ifelse( compare$dOY>355,'winter',compare$season )
  
  # Metrics
  rmseAir  <- round(sqrt(sum((compare$airTempDay-compare$airTempSmith)^2, na.rm=T)/nrow(compare)), digits = 3)
  rmsePrcp <- round(sqrt(sum((compare$prcpDay-compare$prcpSmith)^2, na.rm=T)/nrow(compare)), digits = 3)
  rmseSrad <- round(sqrt(sum((compare$sradDay-compare$sradSmith)^2, na.rm=T)/nrow(compare)), digits = 3)
  
  
  # GGploting (need to melt first, then plot):
  # ==========================================
  airMelt <- melt(compare[,c('date','airTempSmith','airTempDay')],id.vars='date')
  air <- ggplot( airMelt, aes(x = date, y = value, colour = variable)) + 
          scale_color_manual(values=c('red', 'blue')) +
          geom_line() +
          ylab('deg C') +
          ggtitle(paste0('Air Temperature Comparison.   RMSE: ', rmseAir))
  
  prcpMelt <- melt(compare[,c('date','prcpSmith','prcpDay')],id.vars='date')
  prcp <- ggplot( prcpMelt, aes(x = date, y = value, colour = variable)) + 
            scale_color_manual(values=c('red', 'blue')) +
            geom_line() +
            ylab('mm') +
            ggtitle(paste0('Precipitation Comparison.   RMSE: ', rmsePrcp))
  
  
  sradMelt <- melt(compare[,c("date","sradSmith","sradDay")],id.vars="date")
  srad <- ggplot( sradMelt, aes(x = date, y = value, colour = variable)) + 
            scale_color_manual(values=c("red", "blue")) +
            geom_line() +
            ylab('W/m^2') +
            ggtitle(paste0('Solar Radiation Comparison.   RMSE: ', rmseSrad))
  
  
  gOut <- arrangeGrob( air, prcp, srad, ncol=1 )

  ggsave(plot=gOut,file=paste('F:/KPONEIL/SourceData/streamTemperature/MA/westbrookAirTemp/sourceCompare/',features[i],".png",sep=''),dpi=300,width=6,height=8, units='in', scale=2)

  # GGploting Daymet v Smith:
  # ==========================================
  airRegSeason <- ggplot( compare, aes(x = airTempDay, y = airTempSmith, colour = season, size = 1.5)) + 
    geom_point() +
    geom_abline() +
    ggtitle(paste0('Daymet v Smith Comparison (Air Temperature by Season)'))

  airRegDOY <- ggplot( compare, aes(x = airTempDay, y = airTempSmith, colour = dOY, size = 1.5)) + 
    scale_color_gradient(low="red", high="blue") +
    geom_point() +
    geom_abline() +
    ggtitle(paste0('Daymet v Smith Comparison (Air Temperature by dOY)'))
  
  gAirOut <- arrangeGrob( airRegSeason, airRegDOY, ncol=1 )
  
  ggsave(plot=gAirOut,file=paste0('F:/KPONEIL/SourceData/streamTemperature/MA/westbrookAirTemp/sourceCompare/AirRegression',features[i],".png"),dpi=300,width=6,height=8, units='in', scale=2)
  
  prcpRegSeason <- ggplot( compare, aes(x = prcpDay, y = prcpSmith, colour = season)) + 
    geom_point() +
    geom_abline() +
    ggtitle(paste0('Daymet v Smith Comparison (Precipitation by Season)'))
  
  prcpRegDOY <- ggplot( compare, aes(x = prcpDay, y = prcpSmith, colour = dOY)) + 
    scale_color_gradient(low="red", high="blue") +
    geom_point() +
    geom_abline() +
    ggtitle(paste0('Daymet v Smith Comparison (Precipitation by dOY)'))
  
  gPrcpOut <- arrangeGrob( prcpRegSeason, prcpRegDOY, ncol=1 )
  
  ggsave(plot=gPrcpOut,file=paste0('F:/KPONEIL/SourceData/streamTemperature/MA/westbrookAirTemp/sourceCompare/PrcpRegression',features[i],".png"),dpi=300,width=6,height=8, units='in', scale=2)
  
}

