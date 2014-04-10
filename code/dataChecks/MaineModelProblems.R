

#=====================================================================
# NE Model validated over ME data with sites removed
meanRMSE meanBiasMean meanBiasSD AIC.df      AIC
vm0S2 3.156222    -2.574663   1.684789      5 475451.8
vm1S2 1.822638  -0.09863293   1.553108     18 429583.9
vm2S2 1.410757   0.08682173   1.115775    138 394389.1
vm3S2 1.358385   0.03003508   1.182799     49 412413.3
vm4S2 1.514816   0.09325699     1.2402     93 403976.1

#Removed Sites:
MEDMR  MEDMR_3WESTBR7.70	2001	RMSE
MEDMR	MEDMR_3WESTBR7.70	2002	RMSE
MEDMR	MEDMR_4BEAVER1.46	2002	RMSE
MEDMR	MEDMR_4WESTER0.54	2005	RMSE
MEDMR	MEDMR_5MAINST64.25	2003	RMSE
MEDMR	MEDMR_5MAINST64.25	2004	RMSE
MEDMR	MEDMR_5MAINST64.25	2005	RMSE
MEDMR	MEDMR_5MOPANG1.78	2001	RMSE
MEDMR	MEDMR_5MOPANG34.81	2006	RMSE
MEDMR	MEDMR_5MOPANG34.81	2010	RMSE
MEDMR	MEDMR_6BEAVER2.13	2002	RMSE
MEDMR	MEDMR_6BEAVER2.13	2004	RMSE
MEDMR	MEDMR_6BEAVER2.13	2005	RMSE
MEDMR	MEDMR_6CHASEM2.17	2002	RMSE
MEDMR	MEDMR_6CHASEM2.17	2004	RMSE
MEDMR	MEDMR_6CHASEM2.17	2005	RMSE

#=====================================================================


#=====================================================================
# NE Model validated over ME data. (Bad RMSE sites included.)

# Segment 2
meanRMSE meanBiasMean meanBiasSD AIC.df      AIC
vm0S2 3.110084    -2.511081   1.671053      5 491790.2
vm1S2 1.833838   -0.1123063   1.537582     18 444929.0
vm2S2 1.457744   0.07091599   1.121146    138 408307.8
vm3S2 1.385221   0.01583692   1.188425     49 427448.0
vm4S2 1.556295   0.07324081   1.244442     93 418194.7
  
# Segment 3
meanRMSE meanBiasMean meanBiasSD AIC.df      AIC
vm0S3  2.61666    -1.971768     1.4931      5 456480.8
vm1S3 1.519299  -0.04706604   1.218885     18 396294.8
vm2S3 1.249932  -0.03156526  0.8711015    132 361776.6
vm3S3 1.220526   0.01925172  0.9267518     49 377510.6
#=====================================================================





















et2[et2$site == 'MEDMR_4BEAVER1.46',]

s ='MEDMR_4BEAVER1.46'


vm2S2$v






#Log of:
swe





et2rem <- et2[et2$site != 'MEDMR_4BEAVER1.46',]


m2S2rem <- lm(temp~(airTemp+airTempLagged1+airTempLagged2+
                   #segment+
                   LatitudeS+LongitudeS+
                   ForestS+ AgricultureS+ 
                   BasinElevationMS+ ReachSlopePCNTS+ TotDASqKMS+ 
                   WetlandOrWaterS+ SurficialCoarseCS+log(ImpoundmentsOpenSqKM+0.001)+ 
                   daylS + sradS + sweS)^2,
           data=et2rem)

vm1S2rem <- validateModel( m2S2rem,et2rem )



#------------------------------------------------------
#Plot some stuff...
#------------------------------------------------------
valRes <- vm2S2$v

valRes$site <- valRes$s
valRes$year <- valRes$y


compare <- merge(valRes, siteData, by = c('site', 'year'), all.x = T, sort = F)



#Need to repeat the scaling in the object that gets plotted...
compare$LatitudeS <- (compare$Latitude-mean(compare$Latitude,na.rm=T))/sd(compare$Latitude,na.rm=T)
compare$LongitudeS <- (compare$Longitude-mean(compare$Longitude,na.rm=T))/sd(compare$Longitude,na.rm=T)
compare$ForestS <-   (compare$Forest-mean(compare$Forest,na.rm=T))/sd(compare$Forest,na.rm=T)

compare$AgricultureLS <-   (log(compare$Agriculture + 0.001) -mean(log(compare$Agriculture + 0.001),na.rm=T))/sd(log(compare$Agriculture + 0.001),na.rm=T)

compare$BasinElevationMS <-   (compare$BasinElevationM-mean(compare$BasinElevationM,na.rm=T))/sd(compare$BasinElevationM,na.rm=T)
compare$ReachSlopePCNTS <-   (compare$ReachSlopePCNT-mean(compare$ReachSlopePCNT,na.rm=T))/sd(compare$ReachSlopePCNT,na.rm=T)

compare$TotDASqKMLS <-   (log(compare$TotDASqKM+0.001) -mean(log(compare$TotDASqKM+0.001),na.rm=T))/sd(log(compare$TotDASqKM+0.001),na.rm=T)

compare$WetlandOrWaterS <-   (compare$WetlandOrWater-mean(compare$WetlandOrWater,na.rm=T))/sd(compare$WetlandOrWater,na.rm=T)
compare$SurficialCoarseCLS <-   (log(compare$SurficialCoarseC + 1) -mean(log(compare$SurficialCoarseC + 1),na.rm=T))/sd(log(compare$SurficialCoarseC + 1),na.rm=T)

compare$ImpoundmentsOpenSqKMLS <-   (log(compare$ImpoundmentsOpenSqKM+1) -mean(log(compare$ImpoundmentsOpenSqKM+1),na.rm=T))/sd(log(compare$ImpoundmentsOpenSqKM+1),na.rm=T)





d <- compare

unique(d[,c('site','year', 'rmse')])

g <- list()
vList <- c(#'airTemp','airTempLagged1','airTempLagged2',
  'LatitudeS','LongitudeS',
  'ForestS', 'AgricultureLS', 
  'BasinElevationMS', 'ReachSlopePCNTS', 'TotDASqKMLS', 
  'WetlandOrWaterS', 'SurficialCoarseCLS'
  #  'daylS' , 'sradS' , 'sweS' 
)

i=0
for(v in vList){
  i=i+1
  print(c(i,v))
  g[[i]] <- ggplot(d, aes_string(x=v,y='rmse')) +
    geom_point() +
    geom_text( aes(label = site))+
    scale_y_log10()
}



























save(compare, file = 'C:/KPONEIL/m2S2rmse.RData')


plot(compare$Forest, compare$rmse)
plot(compare$Agriculture, compare$rmse)
plot(compare$BasinElevationM, compare$rmse)
plot(compare$Agriculture, compare$rmse)
plot(compare$Agriculture, compare$rmse)





model = m2S2

d = valData
  
s = 'MEDMR_3MAINST1.84'

y = 2004

#s = 'MEDMR_4BEAVER1.46'

# leave out one site validation
validateModel <- function ( model, d, vd ) {
  f=model$call
  siteList <- unique(vd$site)
  v <- data.frame(s=NA,y=NA,n=NA,rmse=NA,biasMean=NA,biasSD=NA,
                  sigma=NA,lat=NA,lon=NA)
  g <- list()
  i=0
  for(s in siteList){
    
    # run model without site = s
    withOut <- d[d$site != s,]    
    m <- lm( f, data=withOut )
    
    yearList <- unique( vd$year[vd$site == s ] )
    
    
    summary(m)
    
    for(y in yearList){
      
      # data for the site/year combo
      with <- vd[vd$site == s & vd$year == y,]
      # predict temp for missing site/year based on all the other sites
      p <- cbind(with,pred=predict(m,with))
      
      i=i+1
      g[[i]] <- ggplot(p[,c('temp','pred')],aes(temp,pred))+
        geom_point()+
        geom_smooth(method='lm') +
        geom_abline(intercept=0,slope=1) +
        scale_x_continuous('Observed temperature') +
        scale_y_continuous('Predicted temperature') +
        ggtitle(paste(as.character(s), as.character(y),sep='_'))
      
      p$bias <- p$pred-p$temp
      # #      v <- rbind(v,c(as.character(s), as.character(y), nrow(p), 
      #                      sqrt(sum(p$bias^2, na.rm=T)/nrow(p)),
      #                      mean(p$bias,na.rm=T),
      #                      sd(p$bias,na.rm=T),
      #                      summary(m)$sigma,
      #                      unique(p$Latitude),unique(p$Longitude))) 
      v <- rbind(v,data.frame(s=as.character(s), 
                              y=y, n=nrow(p), 
                              rmse=sqrt(sum(p$bias^2, na.rm=T)/nrow(p)),
                              biasMean=mean(p$bias,na.rm=T),
                              biasSD=sd(p$bias,na.rm=T),
                              sigma=summary(m)$sigma,
                              lat=unique(p$Latitude),lon=unique(p$Longitude)))
      
    }
  }
  v <- v[-1,] #get rid of first dummy row
  v$biasMeanAbs <- abs(v$biasMean)
  v$biasMeanDir <- ifelse( v$biasMean>0,1,-1 )
  
  meanRMSE <- mean((v$rmse),na.rm=T)
  meanBiasMean <- mean((v$biasMean),na.rm=T)
  meanBiasSD <- mean((v$biasSD),na.rm=T)
  
  return( list(v=v,g=g,means=list(meanRMSE=meanRMSE,meanBiasMean=meanBiasMean,meanBiasSD=meanBiasSD) ) )
}
















#===================================================
#Look at errors:
#===================================================

et2[,c('ImpoundmentsOpenSqKM', 'WetlandsOpenSqKM')]


et2Pairs <- pairs(~SurficialCoarseC+log(ImpoundmentsOpenSqKM+0.001)+ log(WetlandsOpenSqKM +0.001),data=et2)



Latitude+Longitude+Forest+ Impervious+ Agriculture+ BasinElevationM+ ReachSlopePCNT+ TotDASqKM+ WetlandOrWater+ 
  
  
  
  
  
  png(filename=paste0(subGraphsDir, '/summerBP/summerBP_',e1$site[1],'_',year,'.png'),width=1000, height=600, bg="white") 

pairs(et2[airTemp+airTempLagged1+airTempLagged2+
            #segment+
            Latitude+Longitude+
            Forest+ Agriculture+ 
            BasinElevationM+ ReachSlopePCNT+ TotDASqKM+ 
            WetlandOrWater+ SurficialCoarseC+log(ImpoundmentsOpenSqKM+0.001)+ log(WetlandsOpenSqKM +0.001)+
            daylS + sradS + sweS)^2,






