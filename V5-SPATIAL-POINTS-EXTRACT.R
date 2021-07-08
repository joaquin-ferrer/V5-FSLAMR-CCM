#0>>Global Setup--------
#0.Load libraries -----------------------------------------------------------------------------


#netcdy libraries
library(ncdf4)


#spatial libraries
library(sp)
library(geosphe3re)
library(rgdal)
library(raster)
library(RStoolbox)


#plotting libraries
library(ggplot2)


#spatial plotting libraries
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)


#list-data handling libraries
library(dplyr)
library(rlist)
library(tidyr)
library(purrr)
library(data.table)



#time series/zoo libraries
library(zoo)


#date-time libraries
library(lubridate)


#0.Plot Color Setup------------
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
r<-rf(100)

library(RColorBrewer)
colstx <- rev(brewer.pal(n = 9, "Spectral")) 
colsindex <- rev(brewer.pal(n = 9, "RdYlBu")) 
colsdelta <- brewer.pal(n = 9, "Reds") 
colsbias <- brewer.pal(n = 9, "PiYG")
colssd <- brewer.pal(n = 9, "Blues")
#
#1.Loading Raster Brick + Cropping---------- 

#1 - Aphrodite
#2 - CMFD
files <- list()
files = list.files('/Volumes/DAT/2_code/data',pattern='*.nc',full.names=TRUE)


#raster handling https://rpubs.com/markpayne/358146
cmfd.1 <- raster::brick(files[3])
# 
# aoi <- extent(107,110,30,32)
# cmfd.crop <- crop(cmfd.1,aoi)
# 
# pts <- rasterToPoints(cmfd.crop)
# pts.full <- data.frame(pts)

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(106,109, 30.0 ,35.0), "SpatialPolygons")
proj4string(CP) <- proj4string(wanzhou)
## Clip the map
library(rgeos)
wanzhou <- gIntersection(wanzhou, CP, byid=TRUE)
wanzhou_df<-fortify(wanzhou)


#crop tiles around area-of-interest using extent
# aoi.wan<-extent(wanzhou) #check aoi
aoi <- extent(107.5,109,30.35,31.1)

crop.wan<-crop(cmfd.1,aoi)
pts.wan <- rasterToPoints(crop.wan)

plot(pts)


pts <- data.frame(pts.wan)
#define coordinate data frame
df.coor <- data.frame(pts[,1],pts[,2])
names(df.coor) <-c("lon",'lat')

#create main list storing coordinate information 
main.list <- list()
main.list <- (split(df.coor, seq(nrow(df.coor))))   



pts <- (as.data.frame(t(pts)))
pts$dates<-rownames(pts)
pts$dates<-as.Date(pts$dates,format="X%Y.%m.%d")

pts<-pts[-1,]
pts<-pts[-2,]

for (i in 1:length(main.list))
{
  ts<-data.frame(pts$dates,pts[,i]*24) #CMFD comes in mm/hr := convert to mm/day
  names(ts) <-c('dates','pr')
  ts<-list(na.omit(ts))
  ts<-list.names(ts,"ts")
  main.list[[i]] <- append(main.list[[i]],ts)
  
}

list.save(main.list,'wan-extent-ts.rdata')


#antecedent ts dataset creation

pts <- data.frame(pts.wan)
#define coordinate data frame
df.coor <- data.frame(pts[,1],pts[,2])
names(df.coor) <-c("lon",'lat')

#create main list storing coordinate information 
main.list <- list()
main.list <- (split(df.coor, seq(nrow(df.coor))))   

pts <- (as.data.frame(t(pts)))
pts$dates<-rownames(pts)
pts$dates<-as.Date(pts$dates,format="X%Y.%m.%d")

pts<-pts[-1,]
pts<-pts[-2,]


for (i in 1:length(main.list))
{
  ts<-data.frame(pts$dates,pts[,i]*24) #CMFD comes in mm/hr := convert to mm/day
  names(ts) <-c('dates','pr')
  ts<-na.omit(ts)
  ts<- ts[as.numeric(strftime(ts$dates, "%m")) %in% 5:7,] 
  ant.30<-mean(ts$pr)
  ant.30<-list(ant.30)
  ant.30<-list.names(ant.30,"ant.30")
  ts<-list(ts)
  ts<-list.names(ts,"ts")
  main.list[[i]] <- append(main.list[[i]],ts)
  main.list[[i]] <- append(main.list[[i]],ant.30)
  
}



list.save(main.list,'wan-extent-jja-ant.30.rdata')


#2. Load Full Extent List ---------------

load('wan-extent-ts.rdata')
mn.inv <-x

#3 Extracting Annual Maximum Series at ALL Pixels ------------------


#3.2.1. Extracting AMS---------------
library(xts)
library(zoo)
tr.cmfd<-Map(c,lapply(mn.inv, '[', 'ts'), lapply(mn.inv, '[', 'lon'),lapply(mn.inv, '[', 'lat'))
#Expanded ams.jja function
# list<-tr.cmfd
# i=1
#with adf test statistic for diagnosing stationarity
library(tseries)

ams.jja <- function(list) 
{
  for(i in 1:length(list))
  {
    df<-list[[i]][['ts']]
    df<-na.omit(df)
    df<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,] 
    
    #ts2ts via xts
    ts<-xts(df$pr, order.by=as.POSIXct(df$dates))
    #apply yearly -> derive max
    
    yr.pr.max <- apply.monthly(ts, max)
    yr.pr.mean <- apply.monthly(ts, mean)
    
  
    #stationarity test
    adf<-adf.test(yr.pr.mean)
    print(adf$p.value)
    
    ts.yr<- fortify.zoo(yr.pr.max)
    
    ypr.max<-data.frame(ts.yr$yr.pr.max)
    names(ypr.max)[1] <- 'pr.d'
    ypr.max <- list(ypr.max)
    ypr.max <- list.names(ypr.max,"ypr.max")
    list[[i]] <- append(list[[i]],ypr.max)
  }
  return(list)
}

tr.cmfd<-ams.jja(tr.cmfd)
list.save(tr.cmfd,'extent.jja.tr.cmfd.rdata')

#3.2.1. Fitting Distributions (JJA)---------------

load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x

#3.2.1.1. Gumbel Fitting Functions---------------

library(fitdistrplus)

#Start:: Define Gumbel Functions
#a := beta


dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
#End :: Define Gumbel Functions

#3.2.1.2. Fitting to all CMFD Pixels--------------- 
for(i in 1:length(jja.tr.cmfd))
{
  df<-jja.tr.cmfd[[i]][['ypr.max']]
  df<-na.omit(df)
  fitgu<-fitdist(df$pr.d,"gumbel",start=list(a=10, b=10))
  fit.gum<-list(fitgu)
  fit.gum<<- list.names(fit.gum,"fit.gum")
  jja.tr.cmfd[[i]] <- append(jja.tr.cmfd[[i]],fit.gum)
}

list.save(jja.tr.cmfd,'extent.jja.tr.cmfd.rdata')

#3.2.1.3. Gumbel Goodness-of-fit Tests--------------

library(goftest)
load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x

for(i in 1:length(jja.tr.cmfd))
{
  fit<-jja.tr.cmfd[[i]][['fit.gum']]
  df<-jja.tr.cmfd[[i]][['ypr.max']]
  df<-na.omit(df)#flag for potential bug
  #Kolmogorov–Smirnov Goodness-of-fit test
  test.ks <- ks.test(df$pr.d,pgumbel,fit$estimate[[1]],fit$estimate[[2]], alternative = "two.sided")
  res.ks <- if(test.ks$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ks <- list(res.ks,test.ks)
  test.ks<-list.names(test.ks,'test.ks')
  names(test.ks) <- c("h.test","ks.stats")
  
  #Anderson–Darling test Goodness-of-fit test
  test.ad <- ad.test(df$pr.d,pgumbel,fit$estimate[[1]],fit$estimate[[2]])
  res.ad <- if(test.ad$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ad <- list(res.ad,test.ad)
  test.ad<-list.names(test.ad,'test.ad') 
  names(test.ad) <- c("h.test","ad.stats")
  
  #Cramer-von Mises Goodness-of-fit test
  test.cvm <- cvm.test(df$pr.d,pgumbel,fit$estimate[[1]],fit$estimate[[2]])
  res.cvm <- if(test.cvm$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.cvm <- list(res.cvm,test.cvm)
  test.cvm<-list.names(test.cvm,'test.cvm')
  names(test.cvm) <- c("h.test","cvm.stats")
  
  #gumbell test compilation and integration
  tests.gum<-list(test.ks,test.ad,test.cvm)
  names(tests.gum) <- c('test.ks','test.ad','test.cvm')
  tests.gum<-list(tests.gum)
  tests.gum<-list.names(tests.gum,'tests.gum')
  jja.tr.cmfd[[i]] <- append(jja.tr.cmfd[[i]],tests.gum)
  
}
list.save(jja.tr.cmfd,'extent.jja.tr.cmfd.rdata')



#3.2.2. Inspecting Gumbel Distributions--------------
load('extent.jja.tr.cmfd.rdata')



tg<-lapply(x, `[[`, 'tests.gum')
#checking cvm
cvm<-lapply(tg, `[`, 'test.cvm')
cv<-lapply(cvm, `[`, 'test.cvm')
cv<-lapply(cvm, `[[`, 'test.cvm')
h<-rbindlist(lapply(cv, `[`, 'h.test'))
#checking ad

cvm<-lapply(tg, `[`, 'test.ad')
cv<-lapply(cvm, `[`, 'test.ad')
cv<-lapply(cvm, `[[`, 'test.ad')
h<-rbindlist(lapply(cv, `[`, 'h.test'))

#checking ks
cvm<-lapply(tg, `[`, 'test.ks')
cv<-lapply(cvm, `[`, 'test.ks')
cv<-lapply(cvm, `[[`, 'test.ks')
h<-rbindlist(lapply(cv, `[`, 'h.test'))
jja.tr.cmfd<-x


fit<-lapply(x, `[[`, 'fit.gum')
plot(fit[[20]])
#checking cvm
par(mfrow=c(3,4))
for(i in 1:12)
{
  cdfcomp(list(fit[[i]]))
}

#3.2.3.  Mapping Return Periods-------------
load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x

#3.2.3.1 Gumbel Function for Return Periods-------------
#Rainfall function corresponding to recurrence period T:

gum.rt<-function(a,b,t.years)
{
  #a := beta
  #b := alpha
  #mu=beta+0.5772*alpha
  mu=a+(0.05772*b)
  Rt= mu-log(log(t.years/(t.years-1)))*b
  return(Rt)
}
a<-jja.tr.cmfd[[1]][['fit.gum']]



#3.2.3.2 Return Periods for all Rainfall-------------

load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x
# map.gum.tr<-Map(c, lapply(mn.inv, '[[', 'events'),lapply(mn.inv, '[[', 'lon_cmfd'),lapply(mn.inv, '[[', 'lat_cmfd'))
map.gum.tr<-Map(c,lapply(jja.tr.cmfd, '[', 'lon'),lapply(jja.tr.cmfd, '[', 'lat'))


for(i in 1:length(jja.tr.cmfd))
{
  #Retrieve Gumbel Fit Parameters
  fit<-jja.tr.cmfd[[i]][['fit.gum']]
  a<-fit$estimate[[1]]
  b<-fit$estimate[[2]]
  #Create df for return periods and rainfall
  T.yr<-(c(2,5,seq(10,200,10)))
  T.yr<-data.frame(T.yr)
  T.yr$R.t<-0
  for (j in 1: nrow(T.yr))
  {
    #Rainfall for each Return Period
    T.yr[j,2]<-gum.rt(a,b,T.yr[j,1])
  }
  tr.gum<-list(T.yr)
  tr.gum<-list.names(tr.gum,"tr.gum")
  jja.tr.cmfd[[i]] <- append(jja.tr.cmfd[[i]], tr.gum)
  map.gum.tr[[i]] <- append(map.gum.tr[[i]], tr.gum)
}  

list.save(jja.tr.cmfd,'extent.jja.tr.cmfd.rdata')
list.save(map.gum.tr,'extent.jja.map.gum.rdata')

#3.2.3.2 Mapping Return Periods-------------
load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x

load('extent.jja.map.gum.rdata')
jja.map.gum<-x



fit.list<-lapply(jja.tr.cmfd,`[`,'fit.gum')
library(fitdistrplus)

for (i in 1:length(fit.list))
{
plot(fit.list[[i]][['fit.gum']])
}


# stat_function(fun=dgumbel, args=list(loc=24, scale=11))


#Expand tr.gum dataframe for all coordinates
# list<-jja.map.gum
# i=1
# #Function to extract map data.frames with return period and coordinates
#t.pr -> return period of your mapped rainfall
#list -> seasonally-filtered distribution fitted - map list 
map.df.extract<-function(t.pr,list)
{
  for(i in 1:length(list))
  {
    df.map <- data.frame(list[[i]][['lon']],list[[i]][['lat']])
    colnames(df.map)<-c('lon','lat')
    df.tr<-list[[i]][['tr.gum']]#potential bug flag::one list for each distribution and season
    df.tr<-df.tr[df.tr$T.yr==t.pr,]
    df.map$pr=df.tr$R.t
    df.map<-list(df.map)
    df.map<-list.names(df.map,'df.map')
    list[[i]]<-append(list[[i]],df.map)
  }
  list<-lapply(list, `[`, 'df.map')
  return(list)
}


#20-Year Return Period 
list.map<-map.df.extract(20,jja.map.gum)
df.20.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.20.map,'t_20_pr.xyz')


list.map<-map.df.extract(2,jja.map.gum)
df.1.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.2.map,'t_2_pr.xyz')


list.map<-map.df.extract(50,jja.map.gum)
df.50.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.50.map,'t_50_pr.xyz')


list.map<-map.df.extract(2,jja.map.gum)
df.2.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.2.map,'t_2_pr.xyz')

list.map<-map.df.extract(5,jja.map.gum)
df.5.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.5.map,'t_5_pr.xyz')

list.map<-map.df.extract(10,jja.map.gum)
df.10.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.10.map,'t_10_pr.xyz')

list.map<-map.df.extract(100,jja.map.gum)
df.100.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
write.csv(df.100.map,'t_100_pr.xyz')

#Import Polygons for China, select region and clip
library(rgeos)
china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wanzhou))
wanzhou_df<-fortify(wanzhou)
# 
# Maps -inventory - cmfd - raingauge


load('rgs.rdata')
rgs <- x
r1<-rgs[[1]]
r1.rg<- data.frame(r1[["lat"]],r1[["lon"]])
names(r1.rg)[1]  <-'lat'
names(r1.rg)[2]  <-'lon'


load('mn7_days_cmfd.rdata')
mn.inv <-x
df.inv<-(lapply( mn.inv,`[[`,'events'))
df.inv.cor<-data.frame(rbindlist(lapply(df.inv,`[`,'lon_c')),rbindlist(lapply(df.inv,`[`,'lat_c')))



df1<-ggplot()+
  #geom_raster(data=df.5.map,aes(x=lon,y=lat,fill=pr))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c,color="Lanslide Events",fill="Lanslide Events"),size=1)+
  geom_point(data=r1.rg,aes(x=lon,y=lat,color="Rain Gauge",fill="Rain Gauge"),pch=22,size=4)+
  #geom_point(data=df.5.map,aes(x=lon,y=lat,color="CMFD Pixels",fill="CMFD Pixels"),pch=5,size=2)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  #scale_fill_gradientn(colours=colssd, breaks=c(50,70,90,110,130,150),  guide = "colourbar",name = "Rainfall [mm/day]")+
  #scale_fill_gradient2(low="white",mid="#1E90FF",high="#000080", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) +
  labs(title="Wanzhou County\nInventory Events and Rain Gauge Location" ,x="Longitude",y="Latitude",color="Legend",fill="Legend") +
  theme_pubr()
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
df1

#Raster maps--------
df.5<-ggplot()+
  geom_raster(data=df.5.map,aes(x=lon,y=lat,fill=pr))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  #geom_point(data=df.5.map,aes(x=lon,y=lat),pch=5)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  #scale_fill_gradientn(colours=colssd, breaks=c(50,70,90,110,130,150),  guide = "colourbar",name = "Rainfall [mm/day]")+
  #scale_fill_gradient2(low="white",mid="#1E90FF",high="#000080", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
  labs(title="Rainfall at T=5 years" ,x="Longitude",y="Latitude",color="Legend") 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1,legend.position = 'none')
df.5



df.20<-ggplot()+
  geom_raster(data=df.20.map,aes(x=lon,y=lat,fill=pr))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  #geom_point(data=df.20.map,aes(x=lon,y=lat),pch=5)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  #scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
  scale_fill_gradient2(low="white",mid="#1E90FF",high="#000080", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
  labs(title="Rainfall at T=20 years" ,x="Longitude",y="Latitude",color="Legend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
df.20




df.50<-ggplot()+
  geom_raster(data=df.50.map,aes(x=lon,y=lat,fill=pr))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  #geom_point(data=df.5.map,aes(x=lon,y=lat),pch=5)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  #scale_fill_gradientn(colours=colssd, breaks=c(50,70,90,110,130,150),  guide = "colourbar",name = "Rainfall [mm/day]")+
  scale_fill_gradient2(low="white",mid="#1E90FF",high="#000080", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
  labs(title="Rainfall at T=50 years" ,x="Longitude",y="Latitude",color="Legend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1,legend.position = 'none')
df.50



df.100<-ggplot()+
  geom_raster(data=df.100.map,aes(x=lon,y=lat,fill=pr))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  #geom_point(data=df.100.map,aes(x=lon,y=lat),pch=5)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  #scale_fill_gradientn(colours=colssd ,guide = "colourbar",name = "Rainfall [mm/day]")+
  scale_fill_gradient2(low="white",mid="#1E90FF",high="#000080", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  #scale_fill_gradient2(low="blue", mid="grey", high="red", midpoint=100,breaks=c(50,70,90,110,130,150), guide = "colourbar",name = "Rainfall [mm/day]",limits=c(floor(50), ceiling(150)))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
  labs(title="Rainfall at T=100 years" ,x="Longitude",y="Latitude",color="Legend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
df.100


library(gridExtra)
library(grid)
grid.arrange(df.5,df.20,df.50,df.100,ncol=2)

#Contour Maps=---------

library(metR)
require(ggpubr)

cont.breaks<-c(seq(50,150,10))

df.5<-ggplot(data=df.5.map,aes(lon,lat))+
  geom_contour_fill(aes(z = pr))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(50,150), name='Precipitation [mm/day]') +
  labs(title="Rainfall at T=5 years" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.5  

ggsave('p_t_5.png')

df.10<-ggplot(data=df.10.map,aes(lon,lat))+
  geom_contour_fill(aes(z = pr))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(50,150), name='Precipitation [mm/day]') +
  labs(title="Rainfall at T=10 years" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.10 
ggsave('p_t_10.png')
df.20<-ggplot(data=df.20.map,aes(lon,lat))+
  geom_contour_fill(aes(z = pr))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(50,150), name='Precipitation [mm/day]') +
  labs(title="Rainfall at T=20 years" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.20 
ggsave('p_t_20.png')
df.50<-ggplot(data=df.50.map,aes(lon,lat))+
  geom_contour_fill(aes(z = pr))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(50,150), name='Precipitation [mm/day]') +
  labs(title="Rainfall at T=50 years" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.50 
ggsave('p_t_50.png')
df.100<-ggplot(data=df.100.map,aes(lon,lat))+
  geom_contour_fill(aes(z = pr))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral", limits = c(50,150), name='Precipitation [mm/day]') +
  labs(title="Rainfall at T=100 years" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme_pubr()+theme(legend.position='bottom')+guides(fill = guide_colourbar(title.position='top',frame.colour ='white', ticks=FALSE,barwidth=15))
df.100 
ggsave('p_t_100.png')

library(gridExtra)
library(grid)
grid.arrange(df.5,df.10,ncol=1)
grid.arrange(df.20,df.50,ncol=1)


#4. Extracting Antecedent Rainfall For ALL Pixels --------------

load('wan-extent-jja-ant.30.rdata')
jja.ant<-x

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(106,109, 30.0 ,35.0), "SpatialPolygons")
proj4string(CP) <- proj4string(wanzhou)
## Clip the map
library(rgeos)
wanzhou <- gIntersection(wanzhou, CP, byid=TRUE)
wanzhou_df<-fortify(wanzhou)

ant.30.map<-data.frame(rbindlist(lapply(jja.ant,`[`,'lon')),rbindlist(lapply(jja.ant,`[`,'lat')),rbindlist(lapply(jja.ant,`[`,'ant.30')))
write.csv(ant.30.map,'ant_30.xyz')



df.ant<-ggplot()+
  geom_raster(data=ant.30.map,aes(x=lon,y=lat,fill=ant.30))+
  #geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
  geom_point(data=ant.30.map,aes(x=lon,y=lat),pch=5)+
  #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
  scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
  labs(title="Antecedent Rainfall (CMFD)\nJune-August" ,x="Longitude",y="Latitude",color="Legend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
df.ant

library(metR)

cont.breaks<-c(seq(5.7,6.5,0.1))
max(ant.30.map$ant.30)
df.ant<-ggplot(data=ant.30.map,aes(lon,lat))+
  geom_contour_fill(aes(z = ant.30))+
  scale_fill_fermenter(breaks = cont.breaks,palette = "Spectral",name='[mm/day]') +
  labs(title="Antecedent Rainfall\n30-Day Averaged Daily Precipitation (June-August)" ,x="Longitude",y="Latitude",color="Legend") +
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1))+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+       
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
df.ant

