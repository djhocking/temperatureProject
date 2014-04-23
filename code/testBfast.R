
library(bfast)

###################################################
#  test effect of h on bps for one site/year combo
###################################################

et1 <- et[et$site %in% 'WB MITCHELL'  & et$year == 2011 & !is.na(et$temp),]#c('date','tempIndex')]
#et1 <- et[et$site %in% 'CTDEP_241'  & !is.na(et$temp),]#c('date','tempIndex')]

#et1 <- data.frame(date=1:length(et1),temp=et1)
plot(et1$date,et1$tempIndex)
ets <- ts(et1$tempIndex,start=c(2010,1),frequency=7)

ets <- ifelse(abs(ets)>2,2,ets)
#ets <- bfastts(et1$temp,et1$date, type='irregular')

s <- seq(0.00825,0.05,0.005)



b <- list(); firstBP <- list(); secondBP <- list()
i <- 0
for (pts in 2:6){ #},c(0.01,0.025,0.05,0.1)){
  i=i+1
  print(h)
#eb <- bfast(ets,h=.0085,max.iter=1,season="none")
eb <- getBfastBP1Or3(ets,2011,1,pts)
b[[i]] <- eb$output[[1]]$Vt.bp

# pull out breakpoints
#biggestDiffIndex <- order(diff(b[[i]]))[length(b[[i]]) - 1] 
#firstBP[[i]] <-  b[[i]][biggestDiffIndex]
#secondBP[[i]] <- b[[i]][biggestDiffIndex + 1]

#b[[i]]$h <- h
plot(eb,main=h)
}

###################################################
#  test summer BP with diff h for one site
###################################################


#et2 <- et1[et1$dOY >= firstBP[[1]] & et1$dOY <= secondBP[[1]],]
#et2 <- et1[et1$dOY >= et1$springBP & et1$dOY <= et1$fallBP,]

et2 <- et[et$site %in% 'WEST BROOK'  
          & !is.na(et$temp),]#c('date','tempIndex')]
#et1 <- data.frame(date=1:length(et1),temp=et1)
table(et2$year,et2$summerBP)

maxNumBPs <- 4
h=1/(maxNumBPs+1)
years=2007:2009 #unique(et2$year)

for( y in years ){
  print(y)
  e <- et2[et2$year == y ,]
    #       et2$dOY >= et2$springBP & et2$dOY <= et2$fallBP,]
  #plot(e$date,e$tempIndex)
  es <- ts(e$temp,start=c(2010,1),frequency=7)

  eb <- bfast(es,h=h,max.iter=1,season="none")
  plot(eb,main=paste(y,maxNumBPs))
}

et2 <- et[et$site %in% 'WEST BROOK' 
          & et$year %in% 2002:2004 
          & !is.na(et$temp),]#c('date','tempIndex')]

maxNumBPs <- 4*3
h=1/(maxNumBPs+1)


es <- ts(et2$temp,start=c(2010,1),frequency=365)

eb <- bfast(es,h=h,max.iter=1000,season="dummy",hpc='foreach')
plot(eb,main=paste(y,maxNumBPs))
