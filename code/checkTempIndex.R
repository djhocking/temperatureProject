ggplot( et[et$year %in% 2003 & et$site %in% 'WEST BROOK',], aes(dOY,temp))+
  geom_point(size=2) +
  geom_point( aes(dOY,airTemp), colour='red', size=2) +
 # geom_point( aes(dOY,(waterDelta-airDelta)/waterDelta), colour='green', size=2) +
  geom_line() +
  geom_line( aes(dOY,airTemp), colour='red') +
#  geom_line( aes(dOY,(waterDelta-airDelta)/waterDelta), colour='green') +
  geom_point( aes(dOY,tempIndex), colour='blue', size=2) +
  geom_point( aes(dOY,movingMean), colour='orange', size=2) +
  geom_hline( aes(yintercept=quantileLo), colour='black') +
  geom_hline( aes(yintercept=quantileHi), colour='black') +
  geom_vline( aes(xintercept=as.numeric(springBP)),size=1.25) +
  geom_vline( aes(xintercept=as.numeric(fallBP)),size=1.25) +
  geom_vline( aes(xintercept=as.numeric(summerBP)),size=1.25) +
  theme_bw() +
  #  theme(panel.grid.major = element_blank(),
  #        panel.grid.minor = element_blank())+
  scale_y_continuous(expression(paste("Water temperature (",degree, "C)", sep = "")),lim=c(-20,30))+ 
  scale_x_continuous('Day of year')+#, limit=c(100,150)) +
  facet_wrap(~year)



ggplot( et[et$year %in% 2003 & et$site %in% 'WEST BROOK',], aes(dOY,temp))+
  geom_point(size=2) +
  geom_point( aes(dOY,airTemp), colour='red', size=2) +
  geom_point( aes(dOY,(waterDelta-airDelta)), colour='green', size=2) +
  geom_line() +
  geom_line( aes(dOY,airTemp), colour='red') +
  geom_line( aes(dOY,rollapply((waterDelta-airDelta), width=window, fill=NA, mean)), colour='green',size=2) +
  geom_point( aes(dOY,airDelta), colour='red', size=2) +
  geom_point( aes(dOY,waterDelta), colour='black', size=2) +
  geom_hline( aes(yintercept=quantileLo), colour='black') +
  geom_hline( aes(yintercept=quantileHi), colour='black') +
  geom_vline( aes(xintercept=as.numeric(springBP)),size=1.25) +
  geom_vline( aes(xintercept=as.numeric(fallBP)),size=1.25) +
  geom_vline( aes(xintercept=as.numeric(summerBP)),size=1.25) +
  theme_bw(base_size=gBaseSize) +
  #  theme(panel.grid.major = element_blank(),
  #        panel.grid.minor = element_blank())+
  scale_y_continuous(expression(paste("Water temperature (",degree, "C)", sep = "")),lim=c(-20,30))+ 
  scale_x_continuous('Day of year')+#, limit=c(100,150)) +
  facet_wrap(~year)


et <- slide(et, Var = "waterDelta", GroupVar = "siteYear", slideBy = 0)


ggplot(et[et$year %in% 2002:2005 & et$site %in% c('WEST BROOK'),], aes(airDelta,waterDelta))+
  geom_point(aes(color=factor(segment)),size=2) +
  geom_smooth(aes(color=factor(segment)),method='lm')+
  geom_abline(intercept=0,slope=0.1)+
  facet_grid(segment~year)


ggplot(et[et$year %in% 2003 & et$site %in% 'WEST BROOK',], aes(dOY,waterDelta))+
  geom_point(size=2) +  
  geom_line(size=1) + 
  geom_point( aes(dOY,airDelta), colour='green', size=2)  +
  geom_line( aes(dOY,airDelta), colour='green', size=1)



tI <- expand.grid(a=10:20,w=seq(5,20,0.1))
tI$tI <- (tI$w-tI$a)/(tI$w + .00000001)
tI$tI2 <- (tI$w/tI$a)
tI$tI3 <- (tI$w-tI$a)

ggplot(tI,aes(a,(tI)))+
  geom_point(aes(color=(w),size=(w))) 

+
  scale_y_continuous(lim=c(-5,6))


ggplot(tI,aes(a,w))+
  geom_point(aes(size=(tI)))+
  geom_abline(intercept=0, slope=1)
  

#stat_density2d(geom="point", aes(size = ..density..), contour = FALSE)
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE)