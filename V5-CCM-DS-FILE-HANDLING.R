#0>>Load libraries -----------------------------------------------------------------------------


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








#Part 1:: CCM PRECIPITATION OUTPUT DATA HANDLING FOR DS===============
#Part 1A: CMFD Observations and Historical TS=========
# his.nc<-nc_open(files[4])
# fut.nc<-nc_open(files[5])
# 4-5 MOHC-HadGEM2 (His, RCP 8.5)
#============================================================
f=4 # f index := file this is  a trial for bigger funtion 
#============================================================


#1>>Load List of observation list files (CMFD) --------
load('wan-extent-ts.rdata')
obs.ts <-x
#1>>Open List of CCM Files -------------

files <- list()
files = list.files('/Volumes/DAT/2_code/data',pattern='*.nc',full.names=TRUE)



library(ncdf4)
nc<- nc_open(files[f], write=FALSE, readunlim=TRUE, verbose=TRUE, auto_GMT=TRUE, suppress_dimvals=FALSE)


#spatial ccm information - dimensions and grid

lat <- ncvar_get(nc,"lat")
lon <- ncvar_get(nc,"lon")


coords_matrix <- cbind(as.vector(lon),as.vector(lat))
mesh_points <- SpatialPoints(coords_matrix, proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
point_ID_matrix <- matrix(1:length(lat),nrow = nrow(lat),ncol = ncol(lat))

# plot(lon,lat,xlim=c(107.9,109),ylim=c(30.5,31),pch=22)
# points(mesh_points,col="red",pch=18)




ij.ind.fun <- 
#Finding corresponding RCM points to CMFD lon, CMFD lat 
function(lon_c,lat_c,lat,mesh_points,point_ID_matrix)
{
  city_point <- SpatialPoints(matrix(c(lon_c,lat_c), nrow=1, ncol=2), proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  closest_point <- which.min(distGeo(mesh_points, city_point))
  ij <- which( point_ID_matrix == closest_point, arr.ind=T)
  i <-ij[1,1]
  j<-ij[1,2]
  ab <- data.frame(i,j)
  return(ab)
}


for (k in 1:length(obs.ts))
{
  ab <- ij.ind.fun(obs.ts[[k]][['lon']],obs.ts[[k]][['lat']],lat,mesh_points,point_ID_matrix)
  ij<-list(ab)
  ij<-list.names(ij,'ij.ccm')
  obs.ts[[k]]<-append(obs.ts[[k]],ij)
}
lon<-rbindlist(lapply(obs.ts,`[`,'lon'))
lat<-rbindlist(lapply(obs.ts,`[`,'lat'))

cmfd<-cbind(lon,lat)
ij<-rbindlist(lapply(obs.ts,`[[`,'ij.ccm'))
df<-unique(data.frame(ij))

#extracting coordinates from ij matrix
for(r in 1:nrow(df))
{
  df[r,3] <- coords_matrix[point_ID_matrix[df[r,1],df[r,2]],1]
  df[r,4] <- coords_matrix[point_ID_matrix[df[r,1],df[r,2]],2]
}
names(df)[3]  <-'lon'
names(df)[4]  <-'lat'

main.list <- list()
main.list <- (split(df, seq(nrow(df))))   


#1>>Matching CMFD Observations to Each CCM Point -----------
#load CMFD data and extract precipitation
cmfd <- raster::brick(files[3])
aoi <- extent(107.5,109,30.35,31.1)
# cmfd.crop<-raster::crop(cmfd,aoi)
k=1
for (k in 1:length(main.list)) 
{
  coords_matrix <- cbind((main.list[[k]][['lon']]),(main.list[[k]][['lat']]))
  points_data <- SpatialPoints(coords_matrix, proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  ts <- raster::extract(cmfd,points_data, buffer=.005,fun=max, df=TRUE) #matches to the nearest raster with a buffer in degrees
  ts <- data.frame(t(ts)) #transpose time series
  ts$dates<-rownames(ts)
  ts$dates<-as.Date(ts$dates,format="X%Y.%m.%d")
  ts <- data.frame(ts[-1,])
  names(ts)[1] <- 'pr.d'
  ts$pr.d<-ts$pr.d*24 #bug fix conversion from mm/hr to mm/day
  ts <- list(ts)
  ts <- list.names(ts,'cmfd.ts')
  main.list[[k]] <- append(main.list[[k]],ts)
}




#precipitation extraction from CCM results


tpr <- ncvar_get(nc,"time")
dates <- as.Date(tpr, origin = '1949-12-01')

pr <- function(f,i,j)
{
  # Correction from seconds to days
  pr <- 86400*ncvar_get(nc, varid = 'pr', start = c(i,j,1), count = c(1,1,-1))
  return(pr)
}

#debugging:: main.list<-rcm.list

for(k in 1:length(main.list))
{
pr.d <- pr(f,main.list[[k]][['i']],main.list[[k]][['j']])
his.ts <- data.frame(dates,pr.d)
his.ts<-list(his.ts)
his.ts<-list.names(his.ts,'his.ts')
main.list[[k]]<-append(main.list[[k]],his.ts)
}


  
  
#Debug check ::rcm.list<-main.list
#Part 1B: Preparation of Future Projections TS==========
#=======================================================
f=f+1
#=======================================================
nc<- nc_open(files[f], write=FALSE, readunlim=TRUE, verbose=TRUE, auto_GMT=TRUE, suppress_dimvals=FALSE)
tpr <- ncvar_get(nc,"time")
dates <- as.Date(tpr, origin = '1949-12-01')

#debugging:: main.list<-rcm.list

for(k in 1:length(main.list))
{
  pr.d <- pr(f,main.list[[k]][['i']],main.list[[k]][['j']])
  fut.ts <- data.frame(dates,pr.d)
  fut.ts<-list(fut.ts)
  fut.ts<-list.names(fut.ts,'fut.ts')
  main.list[[k]]<-append(main.list[[k]],fut.ts)
}


list.save(main.list,'ts.mohc-hadgem.rdata')
# list.save(rcm.list,'ts.mohc-hadgem.rdata')



#Part 1C: Partitioning Training and Validation Lists==========

load('ts.mohc-hadgem.rdata')
rcm.list<-x
x<-NULL
ds.list<-Map(c,lapply(rcm.list, '[', 'i'), lapply(rcm.list, '[', 'j'),lapply(rcm.list, '[', 'lon'),lapply(rcm.list, '[', 'lat'))
#-- Training >> 1979-2005
for (k in 1:length(ds.list))
{
  a<-select(rcm.list[[k]][['cmfd.ts']], c('dates','pr.d'))
  names(a)[2]<-'cmfd.pr'
  b<-rcm.list[[k]][['his.ts']]
  names(b)[2]<-'his.pr'
  obs.ts<-merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))  
  obs.ts<-list(obs.ts)
  obs.ts<- list.names(obs.ts,'train.ts')
  ds.list[[k]]<-append(ds.list[[k]],obs.ts)
  }
#-- Validation >> 2005-2018
for (k in 1:length(ds.list))
{
  a<-dplyr::select(rcm.list[[k]][['cmfd.ts']], c('dates','pr.d'))
  names(a)[2]<-'cmfd.pr'
  b<-rcm.list[[k]][['fut.ts']]
  names(b)[2]<-'fut.pr'
  obs.ts<-merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))  
  obs.ts<-list(obs.ts)
  obs.ts<- list.names(obs.ts,'val.ts')
  ds.list[[k]]<-append(ds.list[[k]],obs.ts)
}


#-- Projections >> Future Time Series
main.list<-Map(c,ds.list,lapply(rcm.list, '[', 'fut.ts'))
list.save(main.list,'ds.mohc-hadgem.rdata')

#-- Anaylsis
load('ds.mohc-hadgem.rdata')
ds.list<-x

#End Script=================

