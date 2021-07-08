#0>>Global Setup--------
#0.Load libraries -----------------------------------------------------------------------------


#netcdy libraries
library(ncdf4)


#spatial libraries
library(sp)
library(geosphere)
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

#0.Load Functions -----------------------------------------------------------------------------
ij.ind.fun <- function(lon_c,lat_c,lat,mesh_points,point_ID_matrix)
{
  city_point <- SpatialPoints(matrix(c(lon_c,lat_c), nrow=1, ncol=2), proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  closest_point <- which.min(distGeo(mesh_points, city_point))
  ij <- which( point_ID_matrix == closest_point, arr.ind=T)
  a <-ij[1,1]
  b<-ij[1,2]
  ab <- data.frame(a,b)
  return(ab)
}

pr <- function(i,ij)
{
  # Correction from seconds to days
  pr <- 86400*ncvar_get(outlist[[i]], varid = 'pr', start = c(as.integer(ij[1]),as.integer(ij[2]),1), count = c(1,1,-1))
  return(pr)
}

antecedent.extraction <- function(prts.re)
{
  for(i in 1:length(prts.re))
  {
    
    df <- data.frame(prts.re[[i]][['dates']], prts.re[[i]][['pr.d']])
    names(df)[1]  <-'dates'
    names(df)[2]  <-'pr.event'
    df$ant.3 <-0
    df$ant.5 <-0
    df$ant.8 <-0
    df$ant.9 <-0
    df$ant.10 <-0
    df$ant.15 <-0
    df$ant.20 <-0
    df$ant.30 <-0
    #df$ant.40 <-0
    dt=3
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.3[index] = df$ant.3[index] + df$pr.event[index-j]
        }
        df$ant.3[index] = df$ant.3[index]
      } 
    }
    
    dt=5
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.5[index] = df$ant.5[index] + df$pr.event[index-j]
        }
        df$ant.5[index] = df$ant.5[index]
      } 
    }
    
    dt=8
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.8[index] = df$ant.8[index] + df$pr.event[index-j]
        }
        df$ant.8[index] = df$ant.8[index]
      } 
    }
    
    
    dt=9
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.9[index] = df$ant.9[index] + df$pr.event[index-j]
        }
        df$ant.9[index] = df$ant.9[index]
      } 
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    dt=10
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.10[index] = df$ant.10[index] + df$pr.event[index-j]
        }
        df$ant.10[index] = df$ant.10[index]
      } 
    }
    dt=15
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > as.Date((df$date[1]-days(1))))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.15[index] = df$ant.15[index] + df$pr.event[index-j]
        }
        df$ant.15[index] = df$ant.15[index]
      } 
    }
    dt=20
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.20[index] = df$ant.20[index] + df$pr.event[index-j]
        }
        df$ant.20[index] = df$ant.20[index]
      } 
    }
    dt=30
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.30[index] = df$ant.30[index] + df$pr.event[index-j]
        }
        df$ant.30[index] = df$ant.30[index]
      } 
    }
    # dt=40
    # for(index in 1:nrow(df))
    # {
    #   
    #   if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
    #   {
    #     ind=df$date[index]
    #     for (j in 1:dt)
    #     {
    #       df$ant.40[index] = df$ant.40[index] + df$pr.event[index-j]
    #     }
    #     df$ant.40[index] = df$ant.40[index]
    #   } 
    # }
    # 
    
    df <-dplyr::select(df,-c(1,2))
    prts.re[[i]][['ant']] <-df
  }
  return(prts.re)
}


antecedent.extraction_long <- function(prts.re) #40,45,60,90
{
  for(i in 1:length(prts.re))
  {
    
    df <- data.frame(prts.re[[i]][['dates']], prts.re[[i]][['pr.d']])
    names(df)[1]  <-'dates'
    names(df)[2]  <-'pr.event'
    df$ant.40 <-0
    df$ant.50 <-0
    df$ant.60 <-0
    df$ant.90 <-0
    dt=40
    for(index in 1:nrow(df))
    {
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.40[index] = df$ant.40[index] + df$pr.event[index-j]
        }
        df$ant.40[index] = df$ant.40[index]
      } 
    }
    dt=50
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > as.Date((df$date[1]-days(1))))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.50[index] = df$ant.50[index] + df$pr.event[index-j]
        }
        df$ant.50[index] = df$ant.50[index]
      } 
    }
    dt=60
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.60[index] = df$ant.60[index] + df$pr.event[index-j]
        }
        df$ant.60[index] = df$ant.60[index]
      } 
    }
    dt=90
    for(index in 1:nrow(df))
    {
      
      if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
      {
        ind=df$date[index]
        for (j in 1:dt)
        {
          df$ant.90[index] = df$ant.90[index] + df$pr.event[index-j]
        }
        df$ant.90[index] = df$ant.90[index]
      } 
    }
    
    df <-dplyr::select(df,-c(1,2))
    prts.re[[i]][['ant']] <-df
  }
  return(prts.re)
}


#1>>Handling--------
#1.Load Inventory----------------------------------------------------------------------------
temp.dat <- read.csv("/Volumes/DAT/2_code/data/0_inventory.csv", header=TRUE) #FIX INVENTORY LAT/LONG LABELS WRONG
temp.dat$date <- as.Date(temp.dat$X ,format="%Y/%m/%d")


#main data.frame
main.df <- data.frame(temp.dat$date,temp.dat$lat,temp.dat$long,temp.dat$depth.m) %>%
  rename(date=temp.dat.date,lat_c=temp.dat.long,lon_c=temp.dat.lat, dep=temp.dat.depth.m)  

#spatial points 

coords_matrix <- cbind(as.vector(main.df[,2]),as.vector(main.df[,3]))
points_data <- SpatialPoints(coords_matrix, proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)

#1.Loading Raster Brick + Cropping---------- 

#1 - Aphrodite
#2 - CMFD
files <- list()
files = list.files('/Volumes/DAT/2_code/data',pattern='*.nc',full.names=TRUE)

#raster handling https://rpubs.com/markpayne/358146
cmfd.1 <- raster::brick(files[3])

#crop tiles around area-of-interest using extent
aoi <- extent(107,110,30,32)
cmfd.crop <- crop(cmfd.1,aoi)

pts <- rasterToPoints(cmfd.crop)
# pts.full <- data.frame(pts)

#1.Extract Time Series for Inventory -------------------

ints <- list()
ints <- raster::extract(cmfd.crop,points_data, buffer=.005,fun=max, df=TRUE) #buffer units degrees
ints <- (as.data.frame(t(ints)))
ints$dates<-rownames(ints)
ints$dates<-as.Date(ints$dates,format="X%Y.%m.%d")
ints <- data.frame(ints[-1,])
ints <- data.frame(ints) #ints -- time series data frame

int <- list()
int <- list(unique(as.list(ints)))

list.save(int, 'cmfd_grid_ts.rdata') #gridded time series based on 

# #Time Series Analysis of grid point TS
# require(ggplot2)
# require(reshape2)
# 
# grid.ts<- rbindlist(int)
# colnames(grid.ts)[32]<-'date'
# df.ts <- melt(grid.ts,  id.vars = 'date', variable.name = 'series')
# df.mean<-select(grid.ts,c(-32))
# a<-colMeans(df.mean)
# a
# b<-unique(a) as.data
# ncol(as.data.frame(b))==ncol(as.data.frame(a)) #use means


load('cmfd_grid_ts.rdata')
x <- data.frame(x)
x=x[,order(ncol(x):1)]
names(x)<-c('dates',1:9)
cmfd.grid.ts<-x




#2>>Processing &  Integration----------
#2.Engineering & Structuring Inventory List Data -------------------
 

pts <- data.frame(pts)
pts <- (as.data.frame(t(pts)))
pts$dates<-rownames(pts)
pts$dates<-as.Date(pts$dates,format="X%Y.%m.%d")
pts <- data.frame(pts[-1,])
pts <- data.frame(pts)

i=1
inv_list <- list()
inv_list <- (split(main.df, seq(nrow(main.df))))   
for (i in 1:length(inv_list)) 
{
  coords_matrix <- cbind(as.vector(inv_list[[i]][[2]]),as.vector(inv_list[[i]][[3]]))
  points_data <- SpatialPoints(coords_matrix, proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  ts <- raster::extract(cmfd.crop,points_data, buffer=.005,fun=max, df=TRUE) #matches to the nearest raster with a buffer in degrees
  ts <- as.data.frame(t(ts)) #transpose time series
  ts$dates<-rownames(ts)
  ts$dates<-as.Date(ts$dates,format="X%Y.%m.%d")
  ts <- data.frame(ts[-1,])
  names(ts)[1] <- 'pr.d'
  ts <- list(ts)
  inv_list[[i]] <- append(inv_list[[i]],ts)
  names(inv_list[[i]]) <- c('date_c', 'lon_c', 'lat_c','dep','ts')
}





list.save(inv_list,'cmfd_inv_list.rdata')

#start
load('cmfd_inv_list.rdata')
g1 <- list.group(x,ts)
names(g1) <- c(1:length(g1))
#modified index search May 5 2021 edit

for(i in 1:length(g1))
{
  a <- data.frame(g1[[i]][[1]][['lon_c']],g1[[i]][[1]][['lat_c']])
  names(a)<-c('x','y')
  sp::coordinates(a)<-1:2
  c.n<-raster::extract(cmfd.crop,a,cellnumbers=T)
  coor.cmfd <- raster::coordinates(cmfd.crop)[c.n[,1],]
  #matching mean from gridded ts to mean from inventory ts to derive CMFD coordinates
  ts <- list(g1[[i]][[1]][['ts']])
  cmfd<- list(coor.cmfd[1],coor.cmfd[2],g1[[i]][[1]][['ts']])
  names(cmfd) <- c('lon_cmfd','lat_cmfd','ts')
  for(j in 1:length(g1[[i]]))
  {
    g1[[i]][[j]][['ts']] <- NULL
  }
  g1[[i]]<-append(g1[[i]],cmfd)
}
list.save(g1,'grouped_cmfd_inv_list.rdata')
#end




#start
load('cmfd_inv_list.rdata')
g1 <- list.group(x,ts)
names(g1) <- c(1:length(g1))
g2 <- list()
for(i in 1:length(g1))
{
  a <- data.frame( map(g1[[i]],'date_c'))
  a <- data.frame(t(a))
  b <- data.frame(map(g1[[i]],'lon_c'))
  b <- data.frame(t(b))
  c <- data.frame(map(g1[[i]],'lat_c'))
  c <- data.frame(t(c))
  d <- data.frame(map(g1[[i]],'dep'))
  d <- data.frame(t(d))
  comp <- data.frame(a,b,c,d)
  colnames(comp) <- c('date_c','lon_c','lat_c','dep')
  comp$date_c<-as.Date(comp$date_c)
  events <- list(comp)
  g2[[i]] <- events
  g2[[i]]<-list.names(g2[[i]],'events')
}  

names(g2)<-c(1:length(g2))
load('grouped_cmfd_inv_list.rdata')
g1<-x

g3<-list()
a<-list()
for(i in 1:length(g1))
{
  a<-list(g1[[i]][['ts']],g1[[i]][['lon_cmfd']],g1[[i]][['lat_cmfd']])
  a[[1]][[1]] <- a[[1]][[1]]*24
  g3[[i]]<-append(g2[[i]],a)
  names(g3[[i]]) <- c('events','ts','lon_cmfd','lat_cmfd')
}
names(g3)<-c(1:length(g3))
#end
list.save(g3,'final_cmfd_inv_list.rdata')




#start - Pe to Events
load('final_cmfd_inv_list.rdata')
cmfd.inv <-x

for(i in 1:length(cmfd.inv))
{
  a <- data.frame(cmfd.inv[[i]][['events']])
  names(a)[1]  <-'dates'
  b <- data.frame(cmfd.inv[[i]][['ts']])
  names(b)[2]  <-'dates'
  names(b)[1]  <-'pr.d'
  c <- merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))
  c$dates <- as.Date(c$dates)
  c <- data.frame (c$dates,c$pr.d)
  colnames(c) <- c('dates','pr.d')
  pr.events <- list(c)
  pr.events <- list.names(pr.events,'pr.events')
  cmfd.inv[[i]]<-append(cmfd.inv[[i]],pr.events)
}
list.save(cmfd.inv,'final_cmfd_inv_list.rdata')
#end - Pe to Events



#start - Pa to ant
load('final_cmfd_inv_list.rdata')
cmfd.inv <-x

ptrs.re<-list()
for(i in 1:length(cmfd.inv))
{
  ptrs.re[[i]] <- cmfd.inv[[i]][['ts']]
}
names(ptrs.re)<-c(1:length(cmfd.inv))
ptrs.re<-antecedent.extraction(ptrs.re)

for(i in 1:length(cmfd.inv))
{
  a <- list(ptrs.re[[i]][[3]])
  a<-list.names(a,'ant.ts')
  cmfd.inv[[i]]<-append(cmfd.inv[[i]],a)
}
list.save(cmfd.inv,'final_cmfd_inv_list.rdata')
#end  - Pa to ant


#start - event antecedent
load('final_cmfd_inv_list.rdata')
cmfd.inv <-x
for(i in 1:length(cmfd.inv))
{
  a <- data.frame(cmfd.inv[[i]][['events']])
  names(a)[1]  <-'dates'
  b <- data.frame(cmfd.inv[[i]][['ts']][['dates']],cmfd.inv[[i]][['ant.ts']])
  names(b)[1]  <-'dates'
  c <- merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))
  c <- list(c[,4:12])
  cmfd.inv[[i]][['pr.events']]<-append(cmfd.inv[[i]][['pr.events']],c)
  cmfd.inv[[i]][['pr.events']] <- data.frame(cmfd.inv[[i]][['pr.events']])
}
list.save(cmfd.inv,'final_cmfd_inv_list.rdata')
#end - event antecedent





#3>>Record Date Uncertainty Incorporation----------
load('final_cmfd_inv_list.rdata')
cmfd.inv <-x
#3.Setting Uncertainty Parameters (Time Shift Window)---------
#Set Time Window for Uncertainty Incorporation
days.buffer.mn=7# Days prior.
days.buffer.pl=0 # Days after was not used, but can be an option
#3.Max Event Search Function---------
require(dplyr)
#Delay function optimized for maximum event search
delay.function <- function(cmfd.inv,days.buffer.mn,days.buffer.pl)
{
  for(i in 1:length(cmfd.inv))
  {
    events.v2 <- cmfd.inv[[i]][['pr.events']]
    daysdiff.v2 <- dplyr::select(events.v2, c(1:2,4:11))
    events.v2 <- dplyr::select(events.v2, c(1:2))
    events.v2[2]=0
    #names(events.v2)[3:10] <-names(cmfd.inv[[i]][['ant.ts']])
    daysdiff.v2[2:10] = 0
    
    colnames(events.v2)[1]<-'recdate'
    events.v2$dates<-events.v2$recdate
    
    
    for(j in 1:length(cmfd.inv[[i]][['pr.events']][['dates']]))
    {
      a <- cmfd.inv[[i]][['pr.events']][['dates']][j]
      b <- a-days.buffer.mn
      c <- a+days.buffer.pl
      
      ts <- cmfd.inv[[i]][['ts']]
      ts<-ts[ts$dates >= b & ts$dates <= c,]
      ts<-ts[ts$pr.d==max(ts$pr.d),]
      ts$date.diff<- ts$dates-a 
      
      #Update 1st loop list
      
      events.v2$pr.d[j]=max(ts$pr.d)
      events.v2$dates[j]=max(ts$dates)
      daysdiff.v2$pr.d[j] =max(ts$date.diff)
      }
    
  
      as.integer(daysdiff.v2$pr.d[j])
      dates <- data.frame(cmfd.inv[[i]][['ts']][['dates']])
      names(dates)[1] = 'dates'
      ant.ts <- cmfd.inv[[i]][['ant.ts']]
      ant.ts <- cbind(dates,ant.ts)
      
      #devnote::update for maximum pr.d dates only
      
      ma <- merge(transform(events.v2, dates = format(as.Date(dates), "%Y-%m-%d")), transform(ant.ts, dates = format(as.Date(dates), "%Y-%m-%d")))
      ma$recdate<-as.Date(ma$recdate)
      events.v2<-ma
      #--end of modification--
      
      #devnote::independent antecedent event search is also available
      #---function below line---
      # ant.ts<-ant.ts[ant.ts$dates >= b & ant.ts$dates <= c,]
      # 
      # for(k in 2:ncol(ant.ts))
      # {
      #   row <- ant.ts[ant.ts[,k]==max(ant.ts[,k]),]
      #   ts.ant <- data.frame(row[1],row[k])
      #   date.diff <- max(ts.ant$dates-a )
      #   #Update 1st loop list
      #   ts.ant <- max(row[,k])
      #   events.v2[j,k+1]=ts.ant
      #   daysdiff.v2[j,k+1]=max(date.diff )
      #   
      # }
      #---function above line---
    
    events.v2<-list(events.v2)
    daysdiff.v2<-list(daysdiff.v2)
    events.v2<-list.names(events.v2,'v2.events')
    daysdiff.v2<-list.names(daysdiff.v2,'v2.days')
    
    cmfd.inv[[i]]<-append(cmfd.inv[[i]],events.v2)
    cmfd.inv[[i]]<-append(cmfd.inv[[i]],daysdiff.v2)
    }
     return(cmfd.inv)
    }


#3. Creating the Time-Shifted List
inv<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
list.save(inv, 'mn7_days_cmfd.rdata')

#4>>Application of EasyBal Recharge-Rainfall Fraction Factors-----------

#Script for Processing Monthly Recharge-Rainfall Fraction Factors
#---commented below---
#Load EasyBal Results CSV
    # refrac <- read.csv('refrac.csv',header=T)
    # colnames(refrac) <- c('ym','p.tot','refrac') 
    # refrac$ym <- as.character(refrac$ym)
    # refrac$ym <- as.Date(paste0(as.character(refrac$ym), '01'),"%Y%m%d") #conversion trick day is artificial
    # 
    # 
    # mn7<- rbindlist(lapply(mn7.inv, `[[`, 'v2.events'))
    # mn7 <- dplyr::select(mn7,c(1,2,10)) #30-day antecedent rainfall analysis 
    # 
    # refrac <- dplyr::select(refrac,c(1,3))
    # colnames(refrac)<-c('dates','refrac')
    # ts.template<- rbindlist(lapply(mn7.inv, `[[`, 'ts'))
    # 
    # 
    # #Very Inefficient Loop but gets the job done.
    # ts.template$refrac<-0
    # for(i in ind:nrow(ts.template))
    # {
    #   k=0
    #   for(j in 1:nrow(refrac))
    #   {
    #     if(ts.template[i,2]>refrac[j,1])
    #     {
    #       k=j
    #     }
    #   }
    #   fac=refrac[k-1,2]
    #   if(k<=1)
    #   {
    #     fac=0
    #   }
    #   ts.template[i,3]=fac
    #   
    # }
    # 
    # ind<-i
    # 
    # #refrac.list <- list(dplyr::select(ts.template,c(2,3)))
    # #list.save(refrac.list,'refrac.rdata')
    # #---commented above---

#4.Creating Recharge Time Series --------------
load('mn7_days_cmfd.rdata')
mn.inv <-x

mn.ant30<- (lapply(mn.inv, `[[`, 'ant.ts'))
mn<- rbindlist(lapply(mn.inv, `[[`, 'v2.events'))
mn <- dplyr::select(mn,c(1:3,11)) #30-day antecedent rainfall analysis 


#4.Sorting Recharge for Events---------
#Start Create Pa.ts
load('refrac.rdata')
refrac <- data.frame(x)
refrac<-unique(refrac)
#refrac and 30-day ant
for(i in 1:length(mn.inv))
{
  a <- refrac
  names(a)[1]  <-'dates'
  b <- data.frame(mn.inv[[i]][['ts']][['dates']],mn.inv[[i]][['ant.ts']][['ant.30']])
  names(b)[1]  <-'dates'
  names(b)[2]  <-'ant.30'
  c <- merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))
  pa <- data.frame(mn.inv[[i]][['ts']][['dates']],c$refrac*c$ant.30,c$refrac*c$ant.30/30)
  names(pa)[1]  <-'dates'
  names(pa)[2]  <-'pa.30' 
  names(pa)[3]  <-'qa.30'
  pa<-list(pa)
  pa<-list.names(pa,'pa.ts')
  mn.inv[[i]]<-append(mn.inv[[i]],pa)
  
}
#list.save(mn.inv,'mn7_pa_inv.rdata')
#End Create Pa.ts
library(dplyr)
i=1
for(i in 1:length(mn.inv))
{
  a <- refrac
  names(a)[1]  <-'dates'
  b <- data.frame(dplyr::select(mn.inv[[i]][['v2.events']],c(1:3)),mn.inv[[i]][['v2.events']][['ant.30']])
  c <- merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))
  pa <- c
  pa$dates <- as.Date(pa$dates)
  names(pa)[5]<-'ant.30'
  pa$ant.30<-pa$ant.30/30
  pa$qa.30 <- pa$refrac*pa$ant.30
  pa<-list(pa)
  pa<-list.names(pa,'pa.events')
  mn.inv[[i]]<-append(mn.inv[[i]],pa)
  
}
list.save(mn.inv,'mn7_days_cmfd.rdata')
# #Start - Create Pa.events

#5>>Exporting to CSV
load('mn7_days_cmfd.rdata')
mn.inv <-x
list.events <- lapply(mn.inv, `[[`, 'pa.events')
cmfd.coor<-data.frame(rbindlist(lapply(mn.inv, `[`, 'lon_cmfd')),rbindlist(lapply(mn.inv, `[`, 'lat_cmfd')))
mn<- cbind(rbindlist(list.events),cmfd.coor)


coor<- rbindlist(lapply(mn.inv, `[[`, 'events'))
coor<-dplyr::select(coor,c(1:3))
colnames(coor)[1]<-'recdate'
coor<-arrange(coor,recdate)
mn<-arrange(mn,recdate)

export<-cbind(mn[,1],mn[,3],coor[,2:3],mn[,2],mn[,4:8])
colnames(export)[1]<-'actual_date'
colnames(export)[2]<-'recorded_date'
# 
# plot(export$lon_c,export$lat_c)
# points(export$lon_cmfd,export$lat_cmfd,col="blue")
write.csv(export,'cmfd-inventory-mn7-day-window.csv')

#=====

#5.1 Calculating Return Period each event date--------------

load('extent.jja.tr.cmfd.rdata')
jja.tr.cmfd<-x

load('extent.jja.map.gum.rdata')
jja.map.gum<-x


load('mn7_days_cmfd.rdata')
mn.inv <-x
list.events <- lapply(mn.inv, `[[`, 'pa.events')
cmfd.coor<-data.frame(rbindlist(lapply(mn.inv, `[`, 'lon_cmfd')),rbindlist(lapply(mn.inv, `[`, 'lat_cmfd')))
mn<- cbind(rbindlist(list.events),cmfd.coor)


coor<- rbindlist(lapply(mn.inv, `[[`, 'events'))
coor<-dplyr::select(coor,c(1:3))
colnames(coor)[1]<-'recdate'
coor<-arrange(coor,recdate)
mn<-arrange(mn,recdate)

export<-cbind(mn[,1],mn[,3],coor[,2:3],mn[,2],mn[,4:8])
colnames(export)[1]<-'actual_date'
colnames(export)[2]<-'recorded_date'













#5.2 Deriving return period from rainfall intensity---------



library(fitdistrplus)

#Start:: Define Gumbel Functions
#a := beta
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))


gumbel.tr<-function(a,b,r.level)
{
  #a := beta
  #b := alpha
  
  # mu=beta+0.5772*alpha
  Tr = 1/(1-(pgumbel(r.level,a,b)))
  # mu=a+(0.05772*b)
  #mu=a+(0.05772*b)
  
  
  # Tr=(exp(exp(((r.level-mu)/(-b))))-1)
  # Rt= mu-log(log(t.years/(t.years-1)))*b
  return(Tr)
}



fit.list<-lapply(jja.tr.cmfd,`[`,'fit.gum')
df.gum.fit<- lapply(fit.list, `[[`, 'fit.gum')
df.gum.fit<- (lapply(df.gum.fit, `[[`, 'estimate'))

df.est<-data.frame(rbindlist(lapply(jja.tr.cmfd,`[`,'lon')),rbindlist(lapply(jja.tr.cmfd,`[`,'lat')), data.frame(unlist(lapply(df.gum.fit, `[[`, 'a'))), data.frame(unlist(lapply(df.gum.fit, `[[`, 'b'))))
colnames(df.est)<-c('lon_cmfd','lat_cmfd','a','b')
# df.est.a<-  data.frame(unlist(lapply(df.gum.fit, `[[`, 'a')))
# df.est.b<- data.frame(unlist(lapply(df.gum.fit, `[[`, 'b')))
# df.est<-data.frame(df.est.a,df.est.b)

export.2<-export
export.2$a<-0
export.2$b<-0
export.2$mu<-0
export.2$Tr<-0

for(i in 1:nrow(export.2))
{
 
  for (j in 1:nrow(df.est))
  {
  
    if(  (as.character(round(df.est$lon_cmfd[j],4))==as.character(round(export.2$lon_cmfd[i],4))) &&  (as.character(round(df.est$lat_cmfd[j],4))==as.character(round(export.2$lat_cmfd[i],4))))
    {
    
      
      export.2$a[i]<-df.est$a[j]
      export.2$b[i]<-df.est$b[j]
     
      # a=df.est$a[j]
      # b=df.est$b[j]
    }
  }
  export.2$mu[i]<-export.2$a[i]+(0.05772*export.2$b[i])
  export.2$Tr[i]<-gumbel.tr(export.2$a[i],export.2$b[i],export.2$pr.d[i])
}
options(digits=4)


write.csv(export.2,'events_with_return_period_estimates.csv')

china3 <- readRDS('gadm36_CHN_3_sp.rds')
wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
CP <- as(extent(106,109, 30.0 ,35.0), "SpatialPolygons")
proj4string(CP) <- proj4string(wanzhou)
## Clip the map
library(rgeos)
wanzhou <- gIntersection(wanzhou, CP, byid=TRUE)
wanzhou_df<-fortify(wanzhou)



library(RColorBrewer)

#0.Plot Color Setup
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
r<-rf(100)

df<-export.2
df<- df[as.numeric(strftime(df$actual_date, "%m")) %in% 6:8,] 
#cols <- c("8" = "red", "4" = "blue", "6" = "darkgreen", "10" = "orange")

library(dplyr)
library(raster)
library(stringr)
library(lubridate)
library(tidyverse)


#Yearly events
df.y <-mutate(df,year = format(actual_date, "%Y")) 
df.m <-mutate(export.2,month = format(actual_date, "%m"))

hist(as.integer(df.y$year),xlim=c(1995,2005),breaks=10)
# plot(density(as.integer(df.y$year)))
 axis(1, at=1995:2005, labels=c(seq(1995,2005,1)))

#largest data sets
df.2005<- df.y[as.numeric(strftime(df.y$actual_date, "%Y")) %in% 2005,] 
df.1998<- df.y[as.numeric(strftime(df.y$actual_date, "%Y")) %in% 1998,]
df.2002<- df.y[as.numeric(strftime(df.y$actual_date, "%Y")) %in% 2002,] 

df.1999<- df.y[as.numeric(strftime(df.y$actual_date, "%Y")) %in% 1999,]

breaks<-c(1,1.5,2,5,10,20)
ggplot()+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="darkblue",size=.4)+
  geom_point(data=df.1998,aes(x=lon_c,y=lat_c,color=Tr,shape="Events"))+
  geom_point(data= df,aes(x=lon_cmfd,y=lat_cmfd,shape="Center of CMFD Pixel"),color='black')+
  scale_color_gradient(low='green', high = "red",breaks=breaks,guide = "colourbar")+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) +
  scale_x_continuous(breaks = seq(107.8,109, by = 0.1))+
  scale_y_continuous(breaks = seq(30.3,31.1, by = 0.1))+
  labs(title="Rainfall Return Period Esimates by Gumbel Fit\n for Inventory Events (June-August) " ,x="Longitude",y="Latitude",color="Event Rainfall Return Period", shape="") +
  scale_shape_manual(name = "", values = c(5, 19), labels = c("CMFD Pixel Centers", "Events")) +
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+#,aspect.ratio=1
  theme(panel.grid.major = element_line(colour = "grey"))





#June-August
df.refcase<-select(df,lon_cmfd,lat_cmfd,pr.d,Tr)
df.refcase<-select(df,lon_cmfd,lat_cmfd)
unique(df.refcase)

139/31
#---------------
ggplot()+
  geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="darkblue",size=.4)+
  geom_point(data= df,aes(x=lon_c,y=lat_c,color=Tr,shape="Events"),size=2.5)+
  geom_point(data= df,aes(x=lon_cmfd,y=lat_cmfd,shape="Center of CMFD Pixel"),color='black')+
  scale_color_gradient(low = "blue", high = "red")+
  coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) +
  scale_x_continuous(breaks = seq(107.8,109, by = 0.1))+
  scale_y_continuous(breaks = seq(30.3,31.1, by = 0.1))+
  labs(title="Rainfall Return Period Esimates by Gumbel Fit\nInventory Events (June-August) " ,x="Longitude",y="Latitude",color="Event Rainfall Return Period", shape="") +
  scale_shape_manual(name = "", values = c(5, 17), labels = c("CMFD Pixel Centers", "Events")) +
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+#,aspect.ratio=1
  theme(panel.grid.major = element_line(colour = "grey"))

  
# 3D Plotting
# colnames(wanzhou_df)[1]<-'lon_cmfd'
# colnames(wanzhou_df)[2]<-'lat_cmfd'
# wanzhou_df$Tr<-0
library(plotly)
# f<-plot_ly(wanzhou_df, x=~lon_cmfd, y=~lat_cmfd,type='scatter3d',mode="lines")
# f<-f %>% add_trace(export.2, z=~Tr, color=~Tr, type="scatter3d", colors = "BrBG",mode="markers")
# f
# 
# fig <- plot_ly(export.2, x=~lon_cmfd,y=~lat_cmfd, color=~Tr, type="scatter3d", colors = "BrBG",mode="markers")
# fig <- fig %>% add_trace(wanzhou_df, x=~lon_cmfd, y=~lat_cmfd, type="area")
# fig

#Plotly color sets Set1 Set2 Set3 Pastel1 Pastel2 Paired Dark2 Accent
# plot_ly(export.2, x=~lon_cmfd,y=~lat_cmfd,z=~Tr, color=~Tr, type="scatter3d", colors = "BrBG",mode="markers")
# plot_ly(export.2, x=~pr.d,y=~qa.30,z=~Tr, type="scatter3d", mode="markers")
par(mfrow=c(1,3))
boxplot(scale(export.2$Tr),main="Return Period")
boxplot(scale(export.2$pr.d), main="Event Rainfall")
boxplot(scale(export.2$qa.30),main="Recharge (30 days)")



par(mfrow=c(1,3))
boxplot( (export.2$Tr),main="Return Period")
boxplot( (export.2$pr.d), main="Event Rainfall")
boxplot( (export.2$qa.30),main="Recharge (30 days)")


#for all events
par(mfrow=c(3,1))
plot(density(export.2$Tr),main="Return Period")
plot(density(export.2$pr.d), main="Event Rainfall")
plot(density(export.2$qa.30), main="Recharge (30 days)")


#jja events
par(mfrow=c(2,1))
plot(density((df$Tr)),main="Return Period")
plot(density(scale(df$pr.d)), main="Event Rainfall\n(Scaled Values)")
par(mfrow=c(2,1))
plot(density(scale(df$ant.30)), main="Antecedent Rainfall\n(Scaled Values)")
plot(density(df$qa.30), main="Recharge (30 days)\n(Scaled Values)")


#jja events
par(mfrow=c(2,1))
hist((df$Tr),main="Return Period")
hist((df$pr.d), main="Event Rainfall")
par(mfrow=c(2,1))
hist(df$ant.30, main="Antecedent Rainfall")
hist(df$qa.30, main="Recharge (30 days)")



hist((export.2$Tr),main="Return Period")
hist((export.2$pr.d), main="Event Rainfall")
hist((export.2$qa.30), main="Recharge (30 days)")

# colnames(df.est)<-c('a','b')
# boxplot(df.est$a,df.est$b)#boxplot with a and b parameters


#End of Script----------
