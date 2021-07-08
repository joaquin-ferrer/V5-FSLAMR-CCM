#=====================================================================================
# ENSEMBLE GCM-RCM ASSEMBLY FOR CCFs
#=====================================================================================

require(Metrics)
require(rlist)
require(dplyr)
months<-function(df)
{
  
  df.m <-dplyr::mutate(df,month = format(as.Date(dates), "%m")) %>% group_by(month)
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


ccf.ext<-function(ext.sc,ext.ref)
{
  #Extremes:: Climate Change Factors (T=100)
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
  #ext.rlev<-ext.rlev[-1,]
  #re.tpr<-re.tpr[-1,]
  return(ext.ccf)
}

#1::Extract all CCFs from RAW RCM Resolution ========

#HADGEM-REMO
load('day_sc_bc_remohadgem.rdata')
day.sc<-x
load('day_val_bc_remohadgem.rdata')
day.ref<-x

load('ext_sc_bc_remohadgem.rdata')
ext.sc<-x
load('ext_val_bc_remohadgem.rdata')
ext.ref<-x


ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  
  ref.mo <-dplyr::mutate(ref.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  ref.mo <- summarise(ref.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  ref.mo$chmo <- as.numeric(ref.mo$month)
  ref.mo<- transform(ref.mo, chmo = month.abb[chmo])
  ref.mo$month <- as.numeric(ref.mo$month)
  
  
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  
  mid.mo <-dplyr::mutate(mid.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  mid.mo <- summarise(mid.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  mid.mo$chmo <- as.numeric(mid.mo$month)
  mid.mo<- transform(mid.mo, chmo = month.abb[chmo])
  mid.mo$month <- as.numeric(mid.mo$month)
  
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  
  lat.mo <-dplyr::mutate(lat.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  lat.mo <- summarise(lat.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  lat.mo$chmo <- as.numeric(lat.mo$month)
  lat.mo<- transform(lat.mo, chmo = month.abb[chmo])
  lat.mo$month <- as.numeric(lat.mo$month)
  
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]

hadrem.ant.ccf<-ant.proj
hadrem.ant.ccf$rcm<-'hadrem'
hadrem.ma.ccf<-hadrem.ant.ccf[hadrem.ant.ccf$Scn=='Mid',]
hadrem.la.ccf<-hadrem.ant.ccf[hadrem.ant.ccf$Scn=='Lat',]

hadrem.ext.ccf<-ccf.ext(ext.sc,ext.ref)
hadrem.ext.ccf$rcm<-"hadrem"
hadrem.me.ccf<-hadrem.ext.ccf[hadrem.ext.ccf$Scn=='Mid',]
hadrem.le.ccf<-hadrem.ext.ccf[hadrem.ext.ccf$Scn=='Lat',]



#HADGEM-REG
load('day_sc_bc_reg4hadgem.rdata')
day.sc<-x
load('day_val_bc_reg4hadgem.rdata')
day.ref<-x

load('ext_sc_bc_reg4hadgem.rdata')
ext.sc<-x
load('ext_val_bc_reg4hadgem.rdata')
ext.ref<-x


ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  
  ref.mo <-dplyr::mutate(ref.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  ref.mo <- summarise(ref.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  ref.mo$chmo <- as.numeric(ref.mo$month)
  ref.mo<- transform(ref.mo, chmo = month.abb[chmo])
  ref.mo$month <- as.numeric(ref.mo$month)
  
  
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  
  mid.mo <-dplyr::mutate(mid.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  mid.mo <- summarise(mid.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  mid.mo$chmo <- as.numeric(mid.mo$month)
  mid.mo<- transform(mid.mo, chmo = month.abb[chmo])
  mid.mo$month <- as.numeric(mid.mo$month)
  
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  
  lat.mo <-dplyr::mutate(lat.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  lat.mo <- summarise(lat.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  lat.mo$chmo <- as.numeric(lat.mo$month)
  lat.mo<- transform(lat.mo, chmo = month.abb[chmo])
  lat.mo$month <- as.numeric(lat.mo$month)
  
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]

hadreg.ant.ccf<-ant.proj
hadreg.ant.ccf$rcm<-'hadreg'
hadreg.ma.ccf<-hadreg.ant.ccf[hadreg.ant.ccf$Scn=='Mid',]
hadreg.la.ccf<-hadreg.ant.ccf[hadreg.ant.ccf$Scn=='Lat',]

hadreg.ext.ccf<-ccf.ext(ext.sc,ext.ref)
hadreg.ext.ccf$rcm<-"hadreg"
hadreg.me.ccf<-hadreg.ext.ccf[hadreg.ext.ccf$Scn=='Mid',]
hadreg.le.ccf<-hadreg.ext.ccf[hadreg.ext.ccf$Scn=='Lat',]


#MPI-ESM-REMO
load('day_sc_bc_remompiesm.rdata')
day.sc<-x
load('day_val_bc_remompiesm.rdata')
day.ref<-x

load('ext_sc_bc_remompiesm.rdata')
ext.sc<-x
load('ext_val_bc_remompiesm.rdata')
ext.ref<-x


ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  
  ref.mo <-dplyr::mutate(ref.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  ref.mo <- summarise(ref.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  ref.mo$chmo <- as.numeric(ref.mo$month)
  ref.mo<- transform(ref.mo, chmo = month.abb[chmo])
  ref.mo$month <- as.numeric(ref.mo$month)
  
  
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  
  mid.mo <-dplyr::mutate(mid.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  mid.mo <- summarise(mid.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  mid.mo$chmo <- as.numeric(mid.mo$month)
  mid.mo<- transform(mid.mo, chmo = month.abb[chmo])
  mid.mo$month <- as.numeric(mid.mo$month)
  
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  
  lat.mo <-dplyr::mutate(lat.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  lat.mo <- summarise(lat.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  lat.mo$chmo <- as.numeric(lat.mo$month)
  lat.mo<- transform(lat.mo, chmo = month.abb[chmo])
  lat.mo$month <- as.numeric(lat.mo$month)
  
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]

mpirem.ant.ccf<-ant.proj
mpirem.ant.ccf$rcm<-'mpirem'
mpirem.ma.ccf<-mpirem.ant.ccf[mpirem.ant.ccf$Scn=='Mid',]
mpirem.la.ccf<-mpirem.ant.ccf[mpirem.ant.ccf$Scn=='Lat',]

mpirem.ext.ccf<-ccf.ext(ext.sc,ext.ref)
mpirem.ext.ccf$rcm<-"mpirem"
mpirem.me.ccf<-mpirem.ext.ccf[mpirem.ext.ccf$Scn=='Mid',]
mpirem.le.ccf<-mpirem.ext.ccf[mpirem.ext.ccf$Scn=='Lat',]



#MPI-ESM-REG
load('day_sc_bc_reg4mpiesm.rdata')
day.sc<-x
load('day_val_bc_reg4mpiesm.rdata')
day.ref<-x

load('ext_sc_bc_reg4mpiesm.rdata')
ext.sc<-x
load('ext_val_bc_reg4mpiesm.rdata')
ext.ref<-x


ant.proj<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=0,MJJMean=0,CCF=0)
for (r in 1:length(day.ref))
{
  
  temp<-data.frame(lon=0,lat=0,May=0,June=0,Jul=0,Aug=0,Scn=c('Ref','Mid','Lat'),MJJMean=0,CCF=0)  
  temp$lat<-day.ref[[r]][['lat']]
  temp$lon<-day.ref[[r]][['lon']]
  
  #Reference Scenarios
  ref.bc<-day.ref[[r]][['ref.bc']][[2]]
  ref.bc$pr.d<-day.ref[[r]][['ref.bc']][[1]]
  
  ref.mo <-dplyr::mutate(ref.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  ref.mo <- summarise(ref.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  ref.mo$chmo <- as.numeric(ref.mo$month)
  ref.mo<- transform(ref.mo, chmo = month.abb[chmo])
  ref.mo$month <- as.numeric(ref.mo$month)
  
  
  temp[1,3:6]<-ref.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[1,8]<-sum(temp[1,3:5])/3
  #Mid Century Scenarios
  mid.bc<-day.sc[[r]][[5]][[1]]
  
  mid.mo <-dplyr::mutate(mid.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  mid.mo <- summarise(mid.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  mid.mo$chmo <- as.numeric(mid.mo$month)
  mid.mo<- transform(mid.mo, chmo = month.abb[chmo])
  mid.mo$month <- as.numeric(mid.mo$month)
  
  temp[2,3:6]<-mid.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[2,8]<-sum(temp[2,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[2,9]<-temp[2,8]/temp[1,8]
  #Late Century Scenarios
  lat.bc<-day.sc[[r]][[6]][[1]]
  
  lat.mo <-dplyr::mutate(lat.bc,month = format(as.Date(dates), "%m")) %>% group_by(month)
  lat.mo <- summarise(lat.mo,pr.mo = mean(pr.d,na.rm = TRUE))
  lat.mo$chmo <- as.numeric(lat.mo$month)
  lat.mo<- transform(lat.mo, chmo = month.abb[chmo])
  lat.mo$month <- as.numeric(lat.mo$month)
  
  temp[3,3:6]<-lat.mo[,2]
  #Antecedent Seasonal Mean (May-July)
  temp[3,8]<-sum(temp[3,3:5])/3
  #CCF - Relative to Reference Scenario
  temp[3,9]<-temp[3,8]/temp[1,8]
  
  ant.proj<-rbind(ant.proj,temp)
  
}  
ant.proj<-ant.proj[-1,]


mpireg.ant.ccf<-ant.proj
mpireg.ant.ccf$rcm<-'mpireg'
mpireg.ma.ccf<-mpireg.ant.ccf[mpireg.ant.ccf$Scn=='Mid',]
mpireg.la.ccf<-mpireg.ant.ccf[mpireg.ant.ccf$Scn=='Lat',]

mpireg.ext.ccf<-ccf.ext(ext.sc,ext.ref)
mpireg.ext.ccf$rcm<-"mpireg"
mpireg.me.ccf<-mpireg.ext.ccf[mpireg.ext.ccf$Scn=='Mid',]
mpireg.le.ccf<-mpireg.ext.ccf[mpireg.ext.ccf$Scn=='Lat',]


#2::Interpolate CCFS to CMFD Grids========


load('extent.jja.tr.cmfd.rdata')
ts.list<-x
cmfd.grid<-data.table::rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))

grid.to.cmfd<-function(df,cmfd.grid)
  #Interpolates Lon, Lat and CCF into a cmfd grid 
  #Bilinear AKIMA interpolation
{
  require(akima)
  
  d1<-interp(df$lon,df$lat,df$CCF,cmfd.grid$lon,cmfd.grid$lat,linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL)
  d2 <- reshape2::melt(d1$z, na.rm = TRUE)
  d2$lon <- d1$x[d2$Var1]
  d2$lat <- d1$y[d2$Var2]
  d2<-na.omit(d2)
  return(d2)
}

#Ensemble CCF Antecedent=========
#Ma

ma.comp<-rbind(hadreg.ma.ccf,hadrem.ma.ccf,mpireg.ma.ccf,mpirem.ma.ccf)
ma.comp<-dplyr::select(ma.comp,c('lon','lat','CCF','rcm'))
ma.ens<-grid.to.cmfd(ma.comp[ma.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(ma.comp[ma.comp$rcm=='hadreg',],cmfd.grid)
ma.ens<-rbind(ma.ens,t1)
t1<-grid.to.cmfd(ma.comp[ma.comp$rcm=='mpireg',],cmfd.grid)
ma.ens<-rbind(ma.ens,t1)
t1<-grid.to.cmfd(ma.comp[ma.comp$rcm=='mpirem',],cmfd.grid)
ma.ens<-rbind(ma.ens,t1)
names(ma.ens)<-c('i','j','CCF','lon','lat')
ens<-ma.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.ma<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))

write.csv(ens.ma,'ens_mid_ant.csv')
#La


la.comp<-rbind(hadreg.la.ccf,hadrem.la.ccf,mpireg.la.ccf,mpirem.la.ccf)

la.ens<-grid.to.cmfd(la.comp[la.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(la.comp[la.comp$rcm=='hadreg',],cmfd.grid)
la.ens<-rbind(la.ens,t1)
t1<-grid.to.cmfd(la.comp[la.comp$rcm=='mpireg',],cmfd.grid)
la.ens<-rbind(la.ens,t1)
t1<-grid.to.cmfd(la.comp[la.comp$rcm=='mpirem',],cmfd.grid)
la.ens<-rbind(la.ens,t1)
names(la.ens)<-c('i','j','CCF','lon','lat')
ens<-la.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens.la<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
ens<-ens[ens$lon>=107.8,]1
write.csv(ens.la,'ens_lat_ant.csv')

#Ensemble CCF T=100 Years========
#Me

me.comp<-rbind(hadreg.me.ccf,hadrem.me.ccf,mpireg.me.ccf,mpirem.me.ccf)
me.comp<-dplyr::select(me.comp,c(1:2,'100','rcm'))
names(me.comp)<-c('lat','lon','CCF','rcm')

me.ens<-grid.to.cmfd(me.comp[me.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='hadreg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpireg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpirem',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
names(me.ens)<-c('i','j','CCF','lon','lat')
ens<-me.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.me.100<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.me.100,'ens_mid_ext_t100.csv')


#Le

le.comp<-rbind(hadreg.le.ccf,hadrem.le.ccf,mpireg.le.ccf,mpirem.le.ccf)
le.comp<-dplyr::select(le.comp,c(1:2,'100','rcm'))
names(le.comp)<-c('lat','lon','CCF','rcm')

le.ens<-grid.to.cmfd(le.comp[le.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='hadreg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpireg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpirem',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
names(le.ens)<-c('i','j','CCF','lon','lat')
ens<-le.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.le.100<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.le.100,'ens_lat_ext_t100.csv')

#Ensemble CCF T=50 Years========
#Me

me.comp<-rbind(hadreg.me.ccf,hadrem.me.ccf,mpireg.me.ccf,mpirem.me.ccf)
me.comp<-dplyr::select(me.comp,c(1:2,'50','rcm'))
names(me.comp)<-c('lat','lon','CCF','rcm')

me.ens<-grid.to.cmfd(me.comp[me.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='hadreg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpireg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpirem',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
names(me.ens)<-c('i','j','CCF','lon','lat')
ens<-me.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.me.50<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.me.50,'ens_mid_ext_t50.csv')


#Le

le.comp<-rbind(hadreg.le.ccf,hadrem.le.ccf,mpireg.le.ccf,mpirem.le.ccf)
le.comp<-dplyr::select(le.comp,c(1:2,'50','rcm'))
names(le.comp)<-c('lat','lon','CCF','rcm')

le.ens<-grid.to.cmfd(le.comp[le.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='hadreg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpireg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpirem',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
names(le.ens)<-c('i','j','CCF','lon','lat')
ens<-le.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.le.50<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.le.50,'ens_lat_ext_t50.csv')

#Ensemble CCF T=20 Years========
#Me

me.comp<-rbind(hadreg.me.ccf,hadrem.me.ccf,mpireg.me.ccf,mpirem.me.ccf)
me.comp<-dplyr::select(me.comp,c(1:2,'20','rcm'))
names(me.comp)<-c('lat','lon','CCF','rcm')

me.ens<-grid.to.cmfd(me.comp[me.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='hadreg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpireg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpirem',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
names(me.ens)<-c('i','j','CCF','lon','lat')
ens<-me.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.me.20<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.me.20,'ens_mid_ext_t20.csv')


#Le

le.comp<-rbind(hadreg.le.ccf,hadrem.le.ccf,mpireg.le.ccf,mpirem.le.ccf)
le.comp<-dplyr::select(le.comp,c(1:2,'20','rcm'))
names(le.comp)<-c('lat','lon','CCF','rcm')

le.ens<-grid.to.cmfd(le.comp[le.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='hadreg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpireg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpirem',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
names(le.ens)<-c('i','j','CCF','lon','lat')
ens<-le.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.le.20<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.le.20,'ens_lat_ext_t20.csv')





#Ensemble CCF T=10 Years========
#Me

me.comp<-rbind(hadreg.me.ccf,hadrem.me.ccf,mpireg.me.ccf,mpirem.me.ccf)
me.comp<-dplyr::select(me.comp,c(1:2,'10','rcm'))
names(me.comp)<-c('lat','lon','CCF','rcm')

me.ens<-grid.to.cmfd(me.comp[me.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='hadreg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpireg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpirem',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
names(me.ens)<-c('i','j','CCF','lon','lat')
ens<-me.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.me.10<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.me.10,'ens_mid_ext_t10.csv')


#Le

le.comp<-rbind(hadreg.le.ccf,hadrem.le.ccf,mpireg.le.ccf,mpirem.le.ccf)
le.comp<-dplyr::select(le.comp,c(1:2,'10','rcm'))
names(le.comp)<-c('lat','lon','CCF','rcm')

le.ens<-grid.to.cmfd(le.comp[le.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='hadreg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpireg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpirem',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
names(le.ens)<-c('i','j','CCF','lon','lat')
ens<-le.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.le.10<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.le.10,'ens_lat_ext_t10.csv')





#Ensemble CCF T=5 Years========
#Me

me.comp<-rbind(hadreg.me.ccf,hadrem.me.ccf,mpireg.me.ccf,mpirem.me.ccf)
me.comp<-dplyr::select(me.comp,c(1:2,'5','rcm'))
names(me.comp)<-c('lat','lon','CCF','rcm')

me.ens<-grid.to.cmfd(me.comp[me.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='hadreg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpireg',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
t1<-grid.to.cmfd(me.comp[me.comp$rcm=='mpirem',],cmfd.grid)
me.ens<-rbind(me.ens,t1)
names(me.ens)<-c('i','j','CCF','lon','lat')
ens<-me.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.me.5<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.me.5,'ens_mid_ext_t5.csv')


#Le

le.comp<-rbind(hadreg.le.ccf,hadrem.le.ccf,mpireg.le.ccf,mpirem.le.ccf)
le.comp<-dplyr::select(le.comp,c(1:2,'5','rcm'))
names(le.comp)<-c('lat','lon','CCF','rcm')

le.ens<-grid.to.cmfd(le.comp[le.comp$rcm=='hadrem',],cmfd.grid)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='hadreg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpireg',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
t1<-grid.to.cmfd(le.comp[le.comp$rcm=='mpirem',],cmfd.grid)
le.ens<-rbind(le.ens,t1)
names(le.ens)<-c('i','j','CCF','lon','lat')
ens<-le.ens%>%group_by(lat,lon) %>%summarise(ccf=mean(CCF),sd=sd(CCF),.group='drop')
ens$cv<-ens$sd/ens$ccf
ens<-ens[ens$lon>=107.8,]
ens.le.5<-dplyr::select(ens,c('lon','lat','ccf','sd','cv'))
write.csv(ens.le.5,'ens_lat_ext_t5.csv')




#Summary of Antecedent and Extreme Statistics ---------

#Late 21st Century Extreme Events:
min(min(ens.le.5$ccf),min(ens.le.10$ccf),min(ens.le.20$ccf),min(ens.le.50$ccf),min(ens.le.100$ccf))
mean(mean(ens.le.5$ccf),mean(ens.le.10$ccf),mean(ens.le.20$ccf),mean(ens.le.50$ccf),mean(ens.le.100$ccf))
max(max(ens.le.5$ccf),max(ens.le.10$ccf),max(ens.le.20$ccf),max(ens.le.50$ccf),max(ens.le.100$ccf))

min(min(ens.le.5$sd),min(ens.le.10$sd),min(ens.le.20$sd),min(ens.le.50$sd),min(ens.le.100$sd))
mean(mean(ens.le.5$sd),mean(ens.le.10$sd),mean(ens.le.20$sd),mean(ens.le.50$sd),mean(ens.le.100$sd))
max(max(ens.le.5$sd),max(ens.le.10$sd),max(ens.le.20$sd),max(ens.le.50$sd),max(ens.le.100$sd))

#Mid 21st Century Extreme Events:
min(min(ens.me.5$ccf),min(ens.me.10$ccf),min(ens.me.20$ccf),min(ens.me.50$ccf),min(ens.me.100$ccf))
mean(mean(ens.me.5$ccf),mean(ens.me.10$ccf),mean(ens.me.20$ccf),mean(ens.me.50$ccf),mean(ens.me.100$ccf))
max(max(ens.me.5$ccf),max(ens.me.10$ccf),max(ens.me.20$ccf),max(ens.me.50$ccf),max(ens.me.100$ccf))

min(min(ens.me.5$sd),min(ens.me.10$sd),min(ens.me.20$sd),min(ens.me.50$sd),min(ens.me.100$sd))
mean(mean(ens.me.5$sd),mean(ens.me.10$sd),mean(ens.me.20$sd),mean(ens.me.50$sd),mean(ens.me.100$sd))
max(max(ens.me.5$sd),max(ens.me.10$sd),max(ens.me.20$sd),max(ens.me.50$sd),max(ens.me.100$sd))

min(min(ens.me.5$cv),min(ens.me.10$cv),min(ens.me.20$cv),min(ens.me.50$cv),min(ens.me.100$cv))
mean(mean(ens.me.5$cv),mean(ens.me.10$cv),mean(ens.me.20$cv),mean(ens.me.50$cv),mean(ens.me.100$cv))
max(max(ens.me.5$cv),max(ens.me.10$cv),max(ens.me.20$cv),max(ens.me.50$cv),max(ens.me.100$cv))


#Summary for Antecedent Events:
min(min(ens.la$ccf),min(ens.ma$ccf))
mean(mean(ens.la$ccf),mean(ens.ma$ccf))
max(max(ens.la$ccf),max(ens.ma$ccf))

min(min(ens.la$sd),min(ens.ma$sd))
mean(mean(ens.la$sd),mean(ens.ma$sd))
max(max(ens.la$sd),max(ens.ma$sd))


min(min(ens.la$cv),min(ens.ma$cv))
mean(mean(ens.la$cv),mean(ens.ma$cv))
max(max(ens.la$cv),max(ens.ma$cv))

#Mid 21st Century Antecedent Events:
min(ens.ma$ccf)
mean(ens.ma$ccf)
max(ens.ma$ccf)

min(ens.ma$sd)
mean(ens.ma$sd)
max(ens.ma$sd)

min(ens.ma$cv)
mean(ens.ma$cv)
max(ens.ma$cv)

#Late 21st Century Antecedent Events:
min(ens.la$ccf)
mean(ens.la$ccf)
max(ens.la$ccf)


min(ens.la$sd)
mean(ens.la$sd)
max(ens.la$sd)

min(ens.la$cv)
mean(ens.la$cv)
max(ens.la$cv)


#Careful Plotting for 2021/07/09 ======
#Conversion from CCF to Relative Change (%)
ens.le.100$ccf<-(ens.le.100$ccf-1)*100
ens.me.100$ccf<-(ens.me.100$ccf-1)*100
library(metR)
require(ggpubr)
require(ggplot2)

#Import Polygons for China, select region and clip
library(rgeos)
library(raster)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.4,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

boxplot(ens.le.100$cv,ens.me.100$cv)

cont.breaks<-c(seq(-20,80,10))

df.100<-ggplot(data=ens.le.100,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "RdYlBu", limits = c(-20,80), name='Climate Change Factor') +
  labs(title="Late 21st Century Projected Change in Precipitation\nEnsemble Mean for T=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.4,31))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_le_ccf_t100.png',dpi=300,width = 20, height = 20, units = "cm")


df.100<-ggplot(data=ens.me.100,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(-20,80), name='Climate Change Factor') +
  labs(title="Mid-21st Century Projected Change in Precipitation\nEnsemble Mean for T=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.4,31))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_me_ccf_t100.png',dpi=300,width = 20, height = 20, units = "cm")




#Le Plots=======
# CCF Plots 

library(metR)
require(ggpubr)
require(ggplot2)

#Import Polygons for China, select region and clip
library(rgeos)
library(raster)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.4,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)

cont.breaks<-c(seq(.5,3,.5))

df.5<-ggplot(data=ens.le.5,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF\nT=5-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_le_ccf_t5.png',dpi=300,width = 20, height = 20, units = "cm")

df.10<-ggplot(data=ens.le.10,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF\nT=10-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.10
ggsave('ens_le_ccf_t10.png',dpi=300,width = 20, height = 20, units = "cm")

df.50<-ggplot(data=ens.le.50,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF\nT=50-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.50
ggsave('ens_le_ccf_t50.png',dpi=300,width = 20, height = 20, units = "cm")

df.100<-ggplot(data=ens.le.100,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF\nT=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_le_ccf_t100.png',dpi=300,width = 20, height = 20, units = "cm")

# SD Plots 

cont.breaks<-c(seq(.1,1,.1))

df.5<-ggplot(data=ens.le.5,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF Standard Deviation\nT=5-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_le_sd_t5.png',dpi=300,width = 20, height = 20, units = "cm")

df.10<-ggplot(data=ens.le.10,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF Standard Deviation\nT=10-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.10
ggsave('ens_le_sd_t10.png',dpi=300,width = 20, height = 20, units = "cm")

df.50<-ggplot(data=ens.le.50,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCFStandard Deviation\nT=50-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.50
ggsave('ens_le_sd_t50.png',dpi=300,width = 20, height = 20, units = "cm")

df.100<-ggplot(data=ens.le.100,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble CCF Standard Deviation\nT=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_le_sd_t100.png',dpi=300,width = 20, height = 20, units = "cm")



#Me Plots=======
# CCF Plots 

library(metR)
require(ggpubr)
require(ggplot2)

#Import Polygons for China, select region and clip
library(rgeos)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)


cont.breaks<-c(seq(.5,3,.5))

df.5<-ggplot(data=ens.me.5,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF\nT=5-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_me_ccf_t5.png',dpi=300,width = 20, height = 20, units = "cm")

df.10<-ggplot(data=ens.me.10,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF\nT=10-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.10
ggsave('ens_me_ccf_t10.png',dpi=300,width = 20, height = 20, units = "cm")

df.50<-ggplot(data=ens.me.50,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF\nT=50-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.50
ggsave('ens_me_ccf_t50.png',dpi=300,width = 20, height = 20, units = "cm")

df.100<-ggplot(data=ens.me.100,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.5,3), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF\nT=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_me_ccf_t100.png',dpi=300,width = 20, height = 20, units = "cm")

# SD Plots 

cont.breaks<-c(seq(.1,1,.1))

df.5<-ggplot(data=ens.me.5,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF Standard Deviation\nT=5-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_me_sd_t5.png',dpi=300,width = 20, height = 20, units = "cm")

df.10<-ggplot(data=ens.me.10,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF Standard Deviation\nT=10-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.10
ggsave('ens_me_sd_t10.png',dpi=300,width = 20, height = 20, units = "cm")

df.50<-ggplot(data=ens.me.50,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCFStandard Deviation\nT=50-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.50
ggsave('ens_me_sd_t50.png',dpi=300,width = 20, height = 20, units = "cm")

df.100<-ggplot(data=ens.me.100,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble CCF Standard Deviation\nT=100-Year Extreme Events" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100
ggsave('ens_me_sd_t100.png',dpi=300,width = 20, height = 20, units = "cm")




#La Plots=======
# CCF Plots 

library(metR)
require(ggpubr)
require(ggplot2)

#Import Polygons for China, select region and clip
library(rgeos)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)


cont.breaks<-c(seq(.4,2,.2))

df.5<-ggplot(data=ens.la,aes(lon,lat))+
  geom_point()+
  coord_sf(xlim = c(107.85,109), ylim = c(30.3,31.05))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)
df.5  

df.5<-ggplot(data=ens.la,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.4,2), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble\nAntecedent Rainfall (May-July) CCF" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.905,109), ylim = c(30.4,31.05))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_la_ccf.png',dpi=300,width = 20, height = 20, units = "cm")

# SD Plots 


cont.breaks<-c(seq(.1,1,.1))
df.5<-ggplot(data=ens.la,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Blues", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Late 21st Century Scenario Ensemble Standard Deviation\nAntecedent Rainfall (May-July) CCF",x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.905,109), ylim = c(30.4,31.05))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5
ggsave('ens_la_sd.png',dpi=300,width = 20, height = 20, units = "cm")




#Ma Plots=======
# CCF Plots 

library(metR)
require(ggpubr)
require(ggplot2)

#Import Polygons for China, select region and clip
library(rgeos)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)


cont.breaks<-c(seq(.4,2,.2))

df.5<-ggplot(data=ens.ma,aes(lon,lat))+
  geom_contour_fill(aes(z = ccf))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.4,2), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble\nAntecedent Rainfall (May-July) CCF" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.905,109), ylim = c(30.4,31))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_ma_ccf.png',dpi=300,width = 20, height = 20, units = "cm")

# SD Plots 

cont.breaks<-c(seq(.1,1,.1))

df.5<-ggplot(data=ens.ma,aes(lon,lat))+
  geom_contour_fill(aes(z = sd))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(.1,1), name='Climate Change Factor') +
  labs(title="Mid 21st Century Scenario Ensemble Standard Deviation\nAntecedent Rainfall (May-July) CCF" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  
ggsave('ens_ma_sd.png',dpi=300,width = 20, height = 20, units = "cm")


#=======


