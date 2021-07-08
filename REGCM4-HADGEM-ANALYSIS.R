
require(Metrics)
require(rlist)
require(dplyr)


months<-function(df)
{
  
  df.m <-mutate(df,month = format(as.Date(dates), "%m")) %>% group_by(month)
  df.m <- summarise(df.m,pr.mo = mean(pr.d,na.rm = TRUE))
  df.m$chmo <- as.numeric(df.m$month)
  df.m<- transform(df.m, chmo = month.abb[chmo])
  df.m$month <- as.numeric(df.m$month)
  return(df.m)
}

jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  return(df.jja)
}

# Fast Analysis ======
#2A:: Antecedent Rainfall CCF===========
load('day_sc_bc_reg4hadgem.rdata')
day.sc<-x
load('day_val_bc_reg4hadgem.rdata')
day.ref<-x
r=1
ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  ref.mo<-months(ref.bc)
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  mid.mo<-months(mid.bc)
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  lat.mo<-months(lat.bc)
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]

#Antecedent Rainfall Climate Change Factors
ra.ccf<-ant.proj[ant.proj$Scn=='Ref',]
ma.ccf<-ant.proj[ant.proj$Scn=='Mid',]
la.ccf<-ant.proj[ant.proj$Scn=='Lat',]

#2B:: Extreme Event CCF===========
load('ext_sc_bc_reg4hadgem.rdata')
ext.sc<-x
load('ext_val_bc_reg4hadgem.rdata')
ext.ref<-x

#Extremes:: Climate Change Factors
ext.ccf<-data.frame(lat=0,lon=0,Scn=0)
y<-data.frame(c(2,5,c=seq(10,200,10)))
years<-data.frame(t(y))
colnames(years)<-y[,1]
ext.ccf<-cbind(ext.ccf,years)

re.tpr<-ext.ccf
ext.rlev<-data.frame(tyr=0,pr=0,Scn=0,lon=0,lat=0)


for (r in 1:length(ext.ref))
{
  temp<-data.frame(lat=0,lon=0,Scn=c('Mid','Lat'))
  temp$lat<-ext.ref[[r]][['lat']]
  temp$lon<-ext.ref[[r]][['lon']]
  
  bc.rl<-ext.ref[[r]][['ref.bc']][['return.levels']]
  bc<-bc.rl$bc
  msc.rl<-ext.sc[[r]][[5]][['return.levels']]
  msc<-msc.rl$bc
  lat.rl<-ext.sc[[r]][[6]][['return.levels']]
  lat<-lat.rl$bc
  years<-bc.rl$years
  rl<-data.frame(years,bc,msc,lat)
  
  tpr<-data.frame(lat=0,lon=0,Scn=c('Ref'))
  tpr$lat<-ext.ref[[r]][['lat']]
  tpr$lon<-ext.ref[[r]][['lon']]
  
  
  
  lev<-data.frame(years,bc)
  lev$Scn<-'ref'
  names(lev)<-c('tyr','pr','Scn')
  lev.msc<-data.frame(years,msc)
  lev.msc$Scn<-'mid'
  names(lev.msc)<-c('tyr','pr','Scn')
  lev.lat<-data.frame(years,lat)
  lev.lat$Scn<-'lat'
  names(lev.lat)<-c('tyr','pr','Scn')
  
  
  lev<-rbind(lev,lev.msc)
  lev<-rbind(lev,lev.lat)
  
  lev$lon<-ext.ref[[r]][['lon']]
  lev$lat<-ext.ref[[r]][['lat']]
  
  demp<-data.frame(t(years))
  names(demp)<-years
  
  temp<-cbind(temp,demp)
  tpr<-cbind(tpr,demp)
  #Mid Century CCF
  msc<-rl$msc/rl$bc
  #Late Century CCF
  lsc<-rl$lat/rl$bc
  
  
  #Reference BC Return Level Magnitudes
  
  
  tpr[1,4:25]<-bc
  temp[1,4:25]<-msc
  temp[2,4:25]<-lsc
  ext.ccf<-rbind(ext.ccf,temp)
  ext.rlev<-rbind(ext.rlev,lev)
  re.tpr<-rbind(re.tpr,tpr)
  
}  
ext.ccf<-ext.ccf[-1,]
ext.rlev<-ext.rlev[-1,]
re.tpr<-re.tpr[-1,]
#Extreme  Climate Change Factors
me.ccf<-ext.ccf[ext.ccf$Scn=='Mid',]
le.ccf<-ext.ccf[ext.ccf$Scn=='Lat',]







#2C:: Bilinear Interpolation Libraries====
require(rgeos)
require(raster)
require(ggplot2)
require(ggpubr)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))

#cols<-c(seq(0.5,2.5,0.3))=======
#2D:: Export T=100: Me, Le=======

cols<-c(seq(0.5,3,0.5))
require(RColorBrewer)
colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
colsdelta <- brewer.pal(n = length(cols), "Reds") 
colsbias <- brewer.pal(n = length(cols), "PiYG")
colssd <- brewer.pal(n = length(cols), "Blues")


df<-dplyr::select(me.ccf,c(1:2,15))
t.yr<-names(df[3])
names(df)[3]<-'val'
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
write.csv(dplyr::select(d2,c(4:5,3)),paste('mid_extreme_ccf_t_100.csv'))
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Extreme Precipication\nClimate Change Factor: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('me_ccf_inter_t_100.png'),dpi=300,width = 20, height = 20, units = "cm")



max(d2$value)
min(d2$value)


df<-dplyr::select(le.ccf,c(1:2,15))
t.yr<-names(df[3])
names(df)[3]<-'val'
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
write.csv(dplyr::select(d2,c(4:5,3)),paste('late_extreme_ccf_t_100.csv'))
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Extreme Precipitation\nClimate Change Factor: T=100 years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('le_ccf_inter_t_100.png'),dpi=300,width = 20, height = 20, units = "cm")



#2E:: Export Antecedent: Ma, La=========


max(d2$value)
min(d2$value)

cols<-c(seq(0.6,4.0,0.5))
require(RColorBrewer)
colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
colsdelta <- brewer.pal(n = length(cols), "Reds") 
colsbias <- brewer.pal(n = length(cols), "PiYG")
colssd <- brewer.pal(n = length(cols), "Blues")



# Mapping Antecedent CCFs
df<-dplyr::select(ma.ccf,c(1:2,9))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)

write.csv(dplyr::select(d2,c(4:5,3)),'mid_antecedent_ccf.csv')
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('ma_ccf_inter.png'),dpi=300,width = 20, height = 20, units = "cm")



df<-dplyr::select(la.ccf,c(1:2,9))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
write.csv(dplyr::select(d2,c(4:5,3)),'late_antecedent_ccf.csv')
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('la_ccf_inter.png'),dpi=300,width = 20, height = 20, units = "cm")


#======






#1::======== Validation and Analysis ========

#1::Testing Metrics >> Day Bias Correction
load('day_val_bc_reg4hadgem.rdata')
day.ref<-x

#Validation Metrics
require('Metrics')
require('stats')
metrics.ref<-data.frame(lon=0,lat=0,raw.mae=0,raw.rmse=0,raw.bias=0,ref.mae=0,ref.rmse=0,ref.bias=0,cv.mae=0,cv.rmse=0,cv.bias=0)
for (r in 1:length(day.ref))
{
  temp<-data.frame(lon=0,lat=0,raw.mae=0,raw.rmse=0,raw.bias=0,ref.mae=0,ref.rmse=0,ref.bias=0,cv.mae=0,cv.rmse=0,cv.bias=0)
  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  bc<-day.ref[[r]][['ref.bc']][[2]]
  bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  obs<-day.ref[[r]][['ref.bc']][[3]]
  train<-day.ref[[r]][['ref.bc']][[2]]
  
  #Raw RCM Performance 
  temp[,3]<-mae(obs$pr.d,train$pr.d)
  temp[,4]<-rmse(obs$pr.d,train$pr.d)
  temp[,5]<-bias(obs$pr.d,train$pr.d)
  #Bias-Correction Performance
  temp[,6]<-mae(obs$pr.d,bc$pr.d)
  temp[,7]<-rmse(obs$pr.d,bc$pr.d)
  temp[,8]<-bias(obs$pr.d,bc$pr.d)
  #Bias-Correction Performance
  cval<-day.ref[[r]][[6]]
  temp[,9]<-mean(cval$mae)
  temp[,10]<-mean(cval$rmse)
  temp[,11]<-mean(cval$bias)
  metrics.ref<-rbind(metrics.ref,temp)
}
metrics.ref<-metrics.ref[-1,]
day.met.ref<-metrics.ref

max(day.met.ref$cv.rmse)

#2::Testing Metrics >> Extreme Bias Correction
load('ext_val_bc_reg4hadgem.rdata')
ext.ref<-x

metrics.ref<-data.frame(lon=0,lat=0,raw.mae=0,raw.rmse=0,raw.bias=0,ref.mae=0,ref.rmse=0,ref.bias=0,cv.mae=0,cv.rmse=0,cv.bias=0)
for (r in 1:length(ext.ref))
{
  temp<-data.frame(lon=0,lat=0,raw.mae=0,raw.rmse=0,raw.bias=0,ref.mae=0,ref.rmse=0,ref.bias=0,cv.mae=0,cv.rmse=0,cv.bias=0)
  
  temp$lat<-ext.ref[[r]][['lat']]
  temp$lon<-ext.ref[[r]][['lon']]
  
  rlevels<-ext.ref[[r]][['ref.bc']][['return.levels']]
  
  bc<-rlevels$bc
  obs<-rlevels$obs
  train<-rlevels$proj
  
  #Raw RCM Performance 
  temp[,3]<-mae(obs ,train )
  temp[,4]<-rmse(obs ,train )
  temp[,5]<-bias(obs ,train )
  #Bias-Correction Performance
  temp[,6]<-mae(obs ,bc )
  temp[,7]<-rmse(obs ,bc )
  temp[,8]<-bias(obs ,bc )
  #Bias-Correction Performance
  cval<-ext.ref[[r]][[6]]
  temp[,9]<-mean(cval$mae)
  temp[,10]<-mean(cval$rmse)
  temp[,11]<-mean(cval$bias)
  metrics.ref<-rbind(metrics.ref,temp)
}
metrics.ref<-metrics.ref[-1,]
ext.met.ref<-metrics.ref
max(ext.met.ref$cv.mae)
# 
# 
boxplot(day.met.ref$cv.rmse,ext.met.ref$cv.rmse,las=1,main='Cross-Validation RMSE\n(1979-2018)',ylab='[mm/day]',names=c('Daily EQDM','Extremes PQDM'))

boxplot(day.met.ref$cv.mae,ext.met.ref$cv.mae,las=1,main='Cross-Validation MAE\n(1979-2018)',ylab='[mm/day]',names=c('Daily EQDM','Extremes PQDM'))
boxplot(day.met.ref$cv.bias,ext.met.ref$cv.bias,las=1,main='Cross-Validation Bias\n(1979-2018)',names=c('Daily EQDM','Extremes PQDM'))
# 



#1A::Analysis of Monthly Means and Accumulations
# r=1
# bc<-day.ref[[r]][['ref.bc']][[2]]
# bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
# obs<-day.ref[[r]][['ref.bc']][[3]]
# train<-day.ref[[r]][['ref.bc']][[2]]
# 
# 
# bc.mo<-months(ref.bc)
# val.mo<-months(val.obs)
# train.mo<-months(train.obs)
# 
# plot(bc.mo$month,bc.mo$pr.mo,type='l',col='red',ylim=c(0,15))
# lines(val.mo$month,val.mo$pr.mo,type='l',col='blue')
# lines(train.mo$month,train.mo$pr.mo,type='l')


#1B::Analysis of Reference Period Distributions
# plot(ecdf(bc$pr.d))
# lines(ecdf(obs$pr.d),col='blue')
# lines(ecdf(train$pr.d),col='red')
# 
# 
# plot(density(ref.bc$pr.d))
# lines(density(val.obs$pr.d),col='blue')
# lines(density(train.obs$pr.d),col='red')


#1C::Analysis on the Validation of Extreme Events



#2::======== Future Projections ========

#2A:: Antecedent Rainfall CCF===========
load('day_sc_bc_reg4hadgem.rdata')
day.sc<-x
load('day_val_bc_reg4hadgem.rdata')
day.ref<-x

ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  ref.mo<-months(ref.bc)
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  mid.mo<-months(mid.bc)
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  lat.mo<-months(lat.bc)
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]

#Antecedent Rainfall Climate Change Factors
ra.ccf<-ant.proj[ant.proj$Scn=='Ref',]
ma.ccf<-ant.proj[ant.proj$Scn=='Mid',]
la.ccf<-ant.proj[ant.proj$Scn=='Lat',]

#2B:: Extreme Event CCF===========
load('ext_sc_bc_reg4hadgem.rdata')
ext.sc<-x
load('ext_val_bc_reg4hadgem.rdata')
ext.ref<-x

#Extremes:: Climate Change Factors
ext.ccf<-data.frame(lat=0,lon=0,Scn=0)
y<-data.frame(c(2,5,c=seq(10,200,10)))
years<-data.frame(t(y))
colnames(years)<-y[,1]
ext.ccf<-cbind(ext.ccf,years)

re.tpr<-ext.ccf
ext.rlev<-data.frame(tyr=0,pr=0,Scn=0,lon=0,lat=0)


for (r in 1:length(ext.ref))
{
  temp<-data.frame(lat=0,lon=0,Scn=c('Mid','Lat'))
  temp$lat<-ext.ref[[r]][['lat']]
  temp$lon<-ext.ref[[r]][['lon']]
  
  bc.rl<-ext.ref[[r]][['ref.bc']][['return.levels']]
  bc<-bc.rl$bc
  msc.rl<-ext.sc[[r]][[5]][['return.levels']]
  msc<-msc.rl$bc
  lat.rl<-ext.sc[[r]][[6]][['return.levels']]
  lat<-lat.rl$bc
  years<-bc.rl$years
  rl<-data.frame(years,bc,msc,lat)
  
  tpr<-data.frame(lat=0,lon=0,Scn=c('Ref'))
  tpr$lat<-ext.ref[[r]][['lat']]
  tpr$lon<-ext.ref[[r]][['lon']]
  
  
  
  lev<-data.frame(years,bc)
  lev$Scn<-'ref'
  names(lev)<-c('tyr','pr','Scn')
  lev.msc<-data.frame(years,msc)
  lev.msc$Scn<-'mid'
  names(lev.msc)<-c('tyr','pr','Scn')
  lev.lat<-data.frame(years,lat)
  lev.lat$Scn<-'lat'
  names(lev.lat)<-c('tyr','pr','Scn')
  
  
  lev<-rbind(lev,lev.msc)
  lev<-rbind(lev,lev.lat)
  
  lev$lon<-ext.ref[[r]][['lon']]
  lev$lat<-ext.ref[[r]][['lat']]
  
  demp<-data.frame(t(years))
  names(demp)<-years
  
  temp<-cbind(temp,demp)
  tpr<-cbind(tpr,demp)
  #Mid Century CCF
  msc<-rl$msc/rl$bc
  #Late Century CCF
  lsc<-rl$lat/rl$bc
  
  
  #Reference BC Return Level Magnitudes
  
  
  tpr[1,4:25]<-bc
  temp[1,4:25]<-msc
  temp[2,4:25]<-lsc
  ext.ccf<-rbind(ext.ccf,temp)
  ext.rlev<-rbind(ext.rlev,lev)
  re.tpr<-rbind(re.tpr,tpr)
  
}  
ext.ccf<-ext.ccf[-1,]
ext.rlev<-ext.rlev[-1,]
re.tpr<-re.tpr[-1,]
#Extreme  Climate Change Factors
me.ccf<-ext.ccf[ext.ccf$Scn=='Mid',]
le.ccf<-ext.ccf[ext.ccf$Scn=='Lat',]

#****Exporting CSV Files==================


write.csv(dplyr::select(la.ccf,c(1:2,9)),'late_antecedent_ccf.csv')
write.csv(dplyr::select(ma.ccf,c(1:2,9)),'mid_antecedent_ccf.csv')
write.csv(dplyr::select(le.ccf,c(-3)),'late_extreme_ccf.csv')
write.csv(dplyr::select(me.ccf,c(-3)),'mid_extreme_ccf.csv')


#Random Mapping===========
require(ggplot2)
require(ggpubr)
ggplot()+
  geom_point(data=ext.rlev,aes(x=tyr,y=pr,color=Scn))+
  labs(title=paste("Return Period Curves\nGCM-RCM: HadGEM-RegCM4"),x="Return Period [years]",y="Precipitation [mm/day]",color="Scenario")+
  theme_pubr()


china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

#2C::Bilinear Interpolation================
require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))
cols<-c(seq(round_any(min(ext.rlev$pr),5)-5,round_any(max(ext.rlev$pr),5)+5,35))
require(RColorBrewer)
colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
colsdelta <- brewer.pal(n = length(cols), "Reds") 
colsbias <- brewer.pal(n = length(cols), "PiYG")
colssd <- brewer.pal(n = length(cols), "Blues")


for (i in 4:ncol(me.ccf))
{
#   df<-dplyr::select(re.tpr,c(1:2,all_of(i)))
#   t.yr<-names(df[3])
#   names(df)[3]<-'val'
#   
#   cols<-c(seq(round_any(min(df$val),5)-5,round_any(max(df$val),5)+5,5))
#   colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
#   d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
#   d2 <- reshape2::melt(d1$z, na.rm = TRUE)
#   d2$lon <- d1$x[d2$Var1]
#   d2$lat <- d1$y[d2$Var2]
#   d2<-na.omit(d2)
#   
#   
#   #Mapping
#   ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
#     geom_tile()+
#     scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='[mm/day]')+
#     geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
#     coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
#     labs(title=paste("Reference Scenario Extreme Precipitation\nReturn Level: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
#   ccf.map.2
#   ggsave(paste('re_tpr_inter_',t.yr,'.png'),dpi=300,width = 20, height = 20, units = "cm")
}
#3::======== Mapping CCFs ========!
#:::======== 
require(rgeos)
require(raster)
require(ggplot2)
require(ggpubr)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))
cols<-c(seq(0.5,round_any(max(le.ccf$'200'),0.1),0.3))


require(RColorBrewer)
colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
colsdelta <- brewer.pal(n = length(cols), "Reds") 
colsbias <- brewer.pal(n = length(cols), "PiYG")
colssd <- brewer.pal(n = length(cols), "Blues")
# MAPPING EXTREME PRECIPITATION CLIMATE CHANGE FACTORS ==============
# BILINEAR INTERPOLATION BETWEEN GCM TO CMFD GRID
# MAP EXPORTS VIA GGPLOT AND GGSAVE
for (i in 4:ncol(me.ccf))
{
  df<-dplyr::select(me.ccf,c(1:2,all_of(i)))
  t.yr<-names(df[3])
  names(df)[3]<-'val'
  
  d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
  d2 <- melt(d1$z, na.rm = TRUE)
  d2$lon <- d1$x[d2$Var1]
  d2$lat <- d1$y[d2$Var2]
  d2<-na.omit(d2)
  # write.csv(dplyr::select(d2,c(4:5,3)),paste('mid_extreme_ccf_t_',t.yr,'.csv'))
  # }
   ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
    geom_tile()+
    scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
    geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
    coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
    labs(title=paste("Mid-21st Century Scenario Extreme Precipication\nClimate Change Factor: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
  ccf.map.2
  ggsave(paste('me_ccf_inter_',t.yr,'.png'),dpi=300,width = 20, height = 20, units = "cm")
  
  p<-ggplot(df,aes(x = lon, y = lat))+
    geom_tile(aes(fill=val),size=10)+
    scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='CCF')+
    geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
    coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
    labs(title=paste("Mid-21st Century Scenario Extreme Precipication\nClimate Change Factor: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
  p
  #ggsave(paste('me_ccf_pts_',t.yr,'.png'),dpi=300,width = 20, height = 20, units = "cm")
}
for (i in 4:ncol(le.ccf))
{
  df<-dplyr::select(le.ccf,c(1:2,all_of(i)))
  t.yr<-names(df[3])
  names(df)[3]<-'val'
  
  d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
  d2 <- melt(d1$z, na.rm = TRUE)
  d2$lon <- d1$x[d2$Var1]
  d2$lat <- d1$y[d2$Var2]
  d2<-na.omit(d2)
# 
#   write.csv(dplyr::select(d2,c(4:5,3)),paste('late_extreme_ccf_t_',t.yr,'.csv'))
# }
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Extreme Precipitation\nClimate Change Factor: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('le_ccf_inter_',t.yr,'.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Extreme Precipitation\nClimate Change Factor: T=",t.yr,"years (June-August)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
#ggsave(paste('le_ccf_pts_',t.yr,'.png'),dpi=300,width = 20, height = 20, units = "cm")
}


# MAPPING ANTECEDENT PRECIPITATION CLIMATE CHANGE FACTORS ==============
# BILINEAR INTERPOLATION BETWEEN GCM TO CMFD GRID
# MAP EXPORTS VIA GGPLOT AND GGSAVE


# Mapping Antecedent CCFs
df<-dplyr::select(ma.ccf,c(1:2,9))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
write.csv(dplyr::select(d2,c(4:5,3)),'mid_antecedent_ccf.csv')
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('ma_ccf_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('ma_ccf_pts.png'),dpi=300,width = 20, height = 20, units = "cm")


df<-dplyr::select(la.ccf,c(1:2,9))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
write.csv(dplyr::select(d2,c(4:5,3)),'late_antecedent_ccf.csv')
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('la_ccf_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='CCF')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Antecedent Precipication\nClimate Change Factor: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('la_ccf_pts.png'),dpi=300,width = 20, height = 20, units = "cm")



#4::======== Mapping Antecedent Precipitation ========
require(rgeos)
require(raster)
require(ggplot2)
require(ggpubr)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))
cols<-c(seq(5,round_any(max(la.ccf$'MJJMean'),.5)+2,1))
require(RColorBrewer)
colstx <- rev(brewer.pal(n = length(cols), "Spectral")) 
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
colsdelta <- brewer.pal(n = length(cols), "Reds") 
colsbias <- brewer.pal(n = length(cols), "PiYG")
colssd <- brewer.pal(n = length(cols), "Blues")
#:::========


df<-dplyr::select(ra.ccf,c(1:2,'MJJMean'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)

#Mapping
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Reference Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('ra_ant_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colstx, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Reference Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('ra_ant_pts.png'),dpi=300,width = 20, height = 20, units = "cm")


df<-dplyr::select(ma.ccf,c(1:2,'MJJMean'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2) 
write.csv(dplyr::select(d2,c(4:5,3)),'mid_antecedent_ccf.csv')

ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('ma_ant_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Mid-21st Century Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('ma_ant_pts.png'),dpi=300,width = 20, height = 20, units = "cm")




df<-dplyr::select(la.ccf,c(1:2,'MJJMean'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2) 
write.csv(dplyr::select(d2,c(4:5,3)),'late_antecedent_ccf.csv')

ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('la_ant_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='mm/day')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Late-21st Century Scenario Antecedent Precipication\nSeasonal Mean Precipitation: May-July\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('la_ant_pts.png'),dpi=300,width = 20, height = 20, units = "cm")


#5::======== Mapping Extreme Cross-Validation ========
require(rgeos)
require(raster)
require(ggplot2)
require(ggpubr)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))

# #RMSE Maps
# cols<-c(seq(0,round_any(max(ext.met.ref$cv.rmse)+5,1),5))
# require(RColorBrewer)
# colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
# colsdelta <- brewer.pal(n = length(cols), "Reds") 
# colsbias <- brewer.pal(n = length(cols), "PiYG")
#:::========
#RMSE
#RMSE Maps
cols<-c(seq(0,round_any(max(ext.met.ref$cv.rmse)+5,1),5))
require(RColorBrewer)
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(ext.met.ref,c(1:2,'cv.rmse'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='RMSE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('rmse_cv_ext_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='RMSE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('rmse_cv_ext_pts.png'),dpi=300,width = 20, height = 20, units = "cm")

#MAE
cols<-c(seq(0,round_any(max(ext.met.ref$cv.mae)+5,1),5))
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(ext.met.ref,c(1:2,'cv.mae'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='MAE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('mae_cv_ext_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='MAE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('mae_cv_ext_pts.png'),dpi=300,width = 20, height = 20, units = "cm")


#bias
cols=NULL
cols<-c(seq(-20,10,5))
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(ext.met.ref,c(1:2,'cv.bias'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='Bias')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('bias_cv_ext_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='Bias')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Extreme Events Parametric Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('bias_cv_ext_pts.png'),dpi=300,width = 20, height = 20, units = "cm")
ccf.map.2

#6::======== Mapping Rainfall Cross-Validation ========
require(rgeos)
require(raster)
require(ggplot2)
require(ggpubr)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

require(akima)
require(RSAGA)
require(reshape2)
require(data.table)

load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))

#RMSE Maps
cols<-c(seq(0,round_any(max(day.met.ref$cv.rmse)+5,1),5))
require(RColorBrewer)
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(day.met.ref,c(1:2,'cv.rmse'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='RMSE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('rmse_cv_day_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='RMSE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('rmse_cv_day_pts.png'),dpi=300,width = 20, height = 20, units = "cm")

#MAE
# cols<-c(seq(0,round_any(max(day.met.ref$cv.mae)+5,1),5))
# colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(day.met.ref,c(1:2,'cv.mae'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='MAE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('mae_cv_day_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='MAE [mm/day]')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('mae_cv_day_pts.png'),dpi=300,width = 20, height = 20, units = "cm")


#bias
cols=NULL
cols<-c(seq(-20,10,5))
colsindex <- rev(brewer.pal(n =length(cols) , "RdYlBu")) 
df<-dplyr::select(day.met.ref,c(1:2,'cv.bias'))
names(df)[3]<-'val'
#Bilinear Interpolation
d1<-z.df<-interp(df$lon,df$lat,df$val,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
d2 <- reshape2::melt(d1$z, na.rm = TRUE)
d2$lon <- d1$x[d2$Var1]
d2$lat <- d1$y[d2$Var2]
d2<-na.omit(d2)
#Mapping
ccf.map.2 <-ggplot(d2,aes(x=lon,y=lat,fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name ='Bias')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
ccf.map.2
ggsave(paste('bias_cv_day_inter.png'),dpi=300,width = 20, height = 20, units = "cm")
p<-ggplot(df,aes(x = lon, y = lat))+
  geom_tile(aes(fill=val),size=10)+
  scale_fill_gradientn(colors=colsindex, guide= "colourbar",breaks=cols,limits=c(min(cols),max(cols)),name='Bias')+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=.8)+
  coord_sf(xlim = c(107.9,109), ylim = c(30.3,31.1)) + 
  labs(title=paste("Daily Precipitation Empirical Quantile Delta Mapping\nBias Correction Cross-Validation (1979-2018)\nGCM-RCM: HadGEM-RegCM4"),x="Longitude",y="Latitude",color="Legend",fill="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)  
p
ggsave(paste('bias_cv_day_pts.png'),dpi=300,width = 20, height = 20, units = "cm")
