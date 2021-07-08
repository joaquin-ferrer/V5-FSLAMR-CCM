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




#0.Functions----------------------

jja.filter <-function(df)
{
  df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  return(df)
}

# #1.Under Construction - Clean Code RGTS Matching Script-----------------
# 
# 
# 
# 
# #Add antecedent rainfall to cmfd with pr.d and dates list prior
# antecedent.extraction <- function(prts.re)
# {
#   for(i in 1:length(prts.re))
#   {
#     
#     df <- data.frame(prts.re[[i]][['dates']], prts.re[[i]][['pr.d']])
#     names(df)[1]  <-'dates'
#     names(df)[2]  <-'pr.event'
#     df$ant.10 <-0
#     df$ant.15 <-0
#     df$ant.20 <-0
#     df$ant.30 <-0
#     df$ant.40 <-0
#     dt=10
#     for(index in 1:nrow(df))
#     {
#       if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
#       {
#         ind=df$date[index]
#         for (j in 1:dt)
#         {
#           df$ant.10[index] = df$ant.10[index] + df$pr.event[index-j]
#         }
#         df$ant.10[index] = df$ant.10[index]
#       } 
#     }
#     dt=15
#     for(index in 1:nrow(df))
#     {
#       
#       if ((df$date[index]-days(dt)) > as.Date((df$date[1]-days(1))))
#       {
#         ind=df$date[index]
#         for (j in 1:dt)
#         {
#           df$ant.15[index] = df$ant.15[index] + df$pr.event[index-j]
#         }
#         df$ant.15[index] = df$ant.15[index]
#       } 
#     }
#     dt=20
#     for(index in 1:nrow(df))
#     {
#       
#       if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
#       {
#         ind=df$date[index]
#         for (j in 1:dt)
#         {
#           df$ant.20[index] = df$ant.20[index] + df$pr.event[index-j]
#         }
#         df$ant.20[index] = df$ant.20[index]
#       } 
#     }
#     dt=30
#     for(index in 1:nrow(df))
#     {
#       
#       if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
#       {
#         ind=df$date[index]
#         for (j in 1:dt)
#         {
#           df$ant.30[index] = df$ant.30[index] + df$pr.event[index-j]
#         }
#         df$ant.30[index] = df$ant.30[index]
#       } 
#     }
#     dt=40
#     for(index in 1:nrow(df))
#     {
#       
#       if ((df$date[index]-days(dt)) > (df$date[1]-days(1)))
#       {
#         ind=df$date[index]
#         for (j in 1:dt)
#         {
#           df$ant.40[index] = df$ant.40[index] + df$pr.event[index-j]
#         }
#         df$ant.40[index] = df$ant.40[index]
#       } 
#     }
#     df <-select(df,-c(1,2))
#     prts.re[[i]][['ant']] <-df
#   }
#   return(prts.re)
# }
# 
# 
# cmfd.2<-antecedent.extraction(cmfd.1)
# list.save(cmfd.2,'cmfd_rgts.rdata')
# 
#1.Loading Variables---------

load('rgs.rdata')
rgs <- x
rg.pe<-data.frame(rgs$r1$dates,rgs$r1$pr.d)
names(rg.pg)[1]  <-'dates'
names(rg.pe)[2]  <-'pr'
rg.pa<-data.frame(rgs$r1$dates,rgs$r1$ant$ant.30)
names(rg.pa)[1]  <-'dates'
names(rg.pa)[2]  <-'ant.30'
rgs<-data.frame(rgs [[1]][['dates']],rgs [[1]][['pr.d']],rgs [[1]][['ant']][['ant.30']])
names(rgs)[1]  <-'dates'
names(rgs)[2]  <-'pr'
names(rgs)[3]  <-'ant.30'

# r2.rg<-data.frame(rgs$r2$dates,rgs$r2$pr.d,rgs$r2$ant$ant.30)
# names(r2.rg)[1]  <-'Date'
# names(r2.rg)[2]  <-'Pr'
# names(r2.rg)[3]  <-'ant.30'


load('cmfd_rgts.rdata')
cmfd.1 <- x
cmfd.pe<-data.frame(cmfd.1[[1]][['dates']],cmfd.1[[1]][['pr.d']])
names(cmfd.pe)[1]  <-'dates'
names(cmfd.pe)[2]  <-'pr'
cmfd.pa<-data.frame(cmfd.1[[1]][['dates']],cmfd.1[[1]][['ant']][['ant.30']])
names(cmfd.pe)[1]  <-'dates'
names(cmfd.pa)[2]  <-'ant.30'
cmfd<-data.frame(cmfd.1[[1]][['dates']],cmfd.1[[1]][['pr.d']],cmfd.1[[1]][['ant']][['ant.30']])
names(cmfd)[1]  <-'dates'
names(cmfd)[2]  <-'pr'
names(cmfd)[3]  <-'ant.30'

# #r2.cmfd<-data.frame(cmfd.1[[2]][['dates']],cmfd.1[[2]][['pr.d']])
# names(r2.cmfd)[1]  <-'Date'
# names(r2.cmfd)[2]  <-'Pr'



#2>>2D Histogram Analysis: Direct Data Analysis--------

#C2.olor Setup------------
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

#2.June-July-August Filters------------------
jja.rg.pa<-jja.filter(rg.pa)
jja.rg.pe<-jja.filter(rg.pe)

jja.cmfd.pa<-jja.filter(cmfd.pa)
jja.cmfd.pe<-jja.filter(cmfd.pe)


jja.rgs<-jja.filter(rgs)
jja.cmfd<-jja.filter(cmfd)



#2.1 (Annual) CMFD vs RG:: Figure 2D Density Plots Side by Side 2D Hist----------------


library(ggplot2)
library(gridExtra)
#Rain Gauge
df<-rgs[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30


rgs.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r ,name = "Percent", labels = scales::percent)+
  ggtitle("Rain Gauge Measurements") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
rgs.center

#CMFD
df<-cmfd[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30

cmfd.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r ,name = "Percent", labels = scales::percent)+
  ggtitle("CMFD at Rain Gauge Location") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
cmfd.center

grid.arrange(cmfd.center,rgs.center, ncol=2, nrow=1)#, widths=c(4, 1), heights=c(1, 4))




#2.2 (Annual) CMFD vs RG:: Figure 2D Density Plots Log Scale Transformations------------
library(ggplot2)
library(gridExtra)
#Rain Gauge
df<-rgs[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30


rgs.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r,trans='log',name = "Percent", labels = scales::percent)+
  ggtitle("Rain Gauge Measurements\nLog Scale Transformation") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
rgs.center

#CMFD
df<-cmfd[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30

cmfd.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r,trans='log',name = "Percent", labels = scales::percent)+
  ggtitle("CMFD at Rain Gauge Location\nLog Scale Transformation") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
cmfd.center

grid.arrange(cmfd.center,rgs.center, ncol=2, nrow=1)#, widths=c(4, 1), heights=c(1, 4))


#2.3 (Annual) CMFD vs RG::Figure 2scatter plot plus density plots-------------


library(ggpubr)
library(patchwork)


center<-ggplot()+
  geom_point(data=cmfd,aes(x=pr,y=ant.30/30,color='CMFD'),alpha=1,size=0.6)+
  geom_point(data=rgs,aes(x=pr,y=ant.30/30,color='Rain Gauge'),alpha=0.6,size=.6,shape=3)+
  theme_pubr()+
  ggtitle("CMFD vs Rain Gauge Measurements ") +
  labs(x="Event [mm]",y="Antecetend[mm/day]",color=" ") 
  

dens1 <- ggplot() + 
  geom_density(data=cmfd,aes(x=pr,fill='CMFD'),alpha = .6)+ 
  geom_density(data=rgs,aes(x=pr,fill='Rain Gauge'),alpha = 0.6)+ 
  theme_void()+ theme(legend.position = 'none')  +labs(fill=" ")

dens2 <- ggplot() + 
  geom_density(data=cmfd,aes(x=ant.30/30,fill='CMFD'),alpha=0.6)+ 
  geom_density(data=rgs,aes(x=ant.30/30,fill='Rain Gauge'),alpha = 0.6)+ 
  theme_void()+ theme(legend.position = "none") +labs(fill=" ") +coord_flip()

dens2
dens1



dens1 + plot_spacer() + center + dens2 + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))











#2.4 (JJA) CMFD vs RG:: Figure 2D Density Plots Side by Side 2D Hist----------------


library(ggplot2)
library(gridExtra)
#Rain Gauge
df<-jja.rgs[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30


rgs.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r ,name = "Percent", labels = scales::percent)+
  ggtitle("Rain Gauge Measurements") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
rgs.center

#CMFD
df<-jja.cmfd[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30

cmfd.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r ,name = "Percent", labels = scales::percent)+
  ggtitle("CMFD at Rain Gauge Location") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
cmfd.center

grid.arrange(cmfd.center,rgs.center, ncol=2, nrow=1)#, widths=c(4, 1), heights=c(1, 4))




#2.2 (JJA) CMFD vs RG:: Figure 2D Density Plots Log Scale Transformations------------
library(ggplot2)
library(gridExtra)
#Rain Gauge
df<-jja.rgs[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30


rgs.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r,trans='log',name = "Percent", labels = scales::percent)+
  ggtitle("Rain Gauge Measurements\nLog Scale Transformation") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
rgs.center

#CMFD
df<-jja.cmfd[,2:3]
colnames(df)<-c('x','y')
df$y<-df$y/30

cmfd.center <- ggplot(df, aes(x,y))+ 
  stat_bin2d(aes(fill = after_stat(density)),bins=c(40))+ 
  scale_fill_gradientn(colours=r,trans='log',name = "Percent", labels = scales::percent)+
  ggtitle("CMFD at Rain Gauge Location\nLog Scale Transformation") +
  labs(x="Event [mm]",y="Antecetend[mm/day]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
cmfd.center

grid.arrange(cmfd.center,rgs.center, ncol=2, nrow=1)#, widths=c(4, 1), heights=c(1, 4))


#2.3 (JJA) CMFD vs RG::Figure 2scatter plot plus density plots-------------








#3.CMFD vs RG:: Annual Maximum Extraction------------------

library(xts)
library(zoo)
# 
# <-xts(rgs[, -1], order.by=as.POSIXct(rgs$Date))
# 
# 
# 
# 
# 
# 
# 
# 
# ams <- apply.yearly(p_xts, max)
# ams <- fortify.zoo(p_xts)
# ams <- na.exclude((ams))
# ams$Index<-as.Date(ams$Index)
# class(ams$Index)
# jja.ams <- ams[as.numeric(strftime(ams$Index, "%m")) %in% 6:8,]
# head(jja.ams,n=2)
# max(jja.ams$p_xts)
# 
# 
# library(fitdistrplus)
# summary(fitdist(ams,'weibull'))





# #2>>2D Histogram Analysis: Direct Data Analysis--------
# 
# #ggplot approach
# 
# rg<-ggplot(rgs,aes(x=pr,y=ant.30/30))+
#   #geom_point()+
#   #stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
#   stat_density2d(aes(colour = ..level..))
# 
# 
# rg
# 
# 
# 
# #data prep
# df<-r1.cmfd[,2:3]
# colnames(df)<-c('x','y')
# df$y<-df$y/30
# 
# 
# 
# 
# # #gplot
# # library(gplots)
# # 
# # h1 <- hist(df$x, breaks=60, plot=F)
# # h2 <- hist(df$y, breaks=25, plot=F)
# # oldpar <- par()
# # par(mar=c(4,4,2,2))
# # layout(matrix(c(2,0,1,3),2,2,byrow=T),c(4,1), c(1,4))
# # #Log-scaled 2d Histogram
# # #h2d <- hist2d(df, nbins=c(60,25), xlab="Event [mm]",ylab="Antecedent [mm/day]" ,col=r, FUN=function(x) log(length(x)))
# # plot(rg)
# # par(mar=c(0,2,1,0))
# # barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
# # par(mar=c(2,0,0.5,1))
# # barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
# # 
# # library(MASS)
# # 
# # # plot.2d.his(df)
# # # 
# # # plot.2d.his<-function(df)
# # # {
# # 
# # 
# # top <- max(h1$counts, h2$counts)
# # #k <- kde2d(df$x, df$y, n=25)
# # k <- kde2d(df$x, df$y, n=10)
# # image(k,col=r)
# # 
# # # margins
# # oldpar <- par()
# # par(mar=c(3,3,1,1))
# # layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
# # image(k,col=r) #plot the image
# # par(mar=c(0,2,1,0))
# # barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
# # par(mar=c(2,0,0.5,1))
# # barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
# # 
# # 
# # #hexbin
# # library(hexbin)
# # # Create hexbin object and plot
# # h <- hexbin(df)
# # plot(h)
# # plot(h, colramp=rf)
# # 
# # 
# 
# # 
# # h2 <- hist2d(df, col=r)
# # # Coarser binsizing and add colouring
# # h2 <- hist2d(df, nbins=25, col=r)
# # # Scaling with log as before
# # h2d <- hist2d(df, col=r, FUN=function(x) log(length(x)))
# 
# #2>>2D Histogram Analysis: AMS (Max Yearly) Data Analysis--------
# 
# 
# 
# 
# library(xts)
# library(zoo)
# 
# r1.rg<-r1.rg[,-2]
# p_xts<-xts(r1.rg[, -1], order.by=as.POSIXct(r1.rg$Date))
# plot(p_xts)
# is.regular(p_xts)
# #p_xts<-xts(r1.cmfd[, -1], order.by=as.POSIXct(r1.cmfd$Date))
# ams <- apply.yearly(p_xts, max)
# plot(ams)
# 
# 
# ams <- fortify.zoo(ams)
# nrow(ams)
# 
# 
# ams.pr<-data.frame(ams[,2])
# names(ams.pr) <- 'pr'
# ams.pr <- na.exclude((ams.pr))
# hist(ams.pr$pr)
# 
# 
# #data prep
# df<-ams[,2:3]
# colnames(df)<-c('x','y')
# df$y<-df$y/30
# 
# 
# # 
# # df.rg<-df
# # df.cmfd<-df
# #Colors
# library(RColorBrewer)
# rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
# r <- rf(32)
# 
# 
# #hexbin
# library(hexbin)
# # Create hexbin object and plot
# h <- hexbin(df)
# plot(h)
# plot(h, colramp=rf)
# 
# 
# 
