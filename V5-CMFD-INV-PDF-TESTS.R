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
#1>> Variable Setup--------
load('mn7_days_cmfd.rdata')
mn.inv <-x
mn<- rbindlist(lapply(mn.inv, `[[`, 'pa.events') )
mn<-data.frame(mn$dates,mn$pr.d,mn$ant.30)
colnames(mn)<-c('dates','pr','ant.30')
jja.mn<- mn[as.numeric(strftime(mn$dates, "%m")) %in% 6:8,]


ggplot()+
  geom_point(data=mn,aes(x=pr,y=ant.30,color='ALL'),alpha=1,size=1.5,shape=19)+
  geom_point(data=jja.mn,aes(x=pr,y=ant.30,color='JJA'),alpha=1,size=3,shape=2)

# 
# list.qats <- lapply(mn.inv, `[[`, 'pa.ts')
# qa.30<- rbindlist(list.qats)


load('final_cmfd_inv_list.rdata')
cmfd.inv <-x
ts<-rbindlist(lapply(cmfd.inv, `[[`, 'ts'))
ant.ts<-rbindlist(lapply(cmfd.inv, `[[`, 'ant.ts'))
# data <- data.frame(ts,ant.ts$ant.30,qa.30$qa.30)
data <- data.frame(ts$dates,ts$pr.d,ant.ts$ant.30/30)
colnames(data)<-c('dates','pr','ant.30')
# data <- data[,-2]
# colnames(data)<-c('pr','ant.30','qa.30')
jja.data<- data[as.numeric(strftime(data$dates, "%m")) %in% 6:8,]

ts<-NULL
ant.ts=NULL
x=NULL



#2>>(Annual) Figure:: scatter plot plus density plots---------

library(ggpubr)
library(patchwork)




center<-ggplot()+
  geom_point(data=data,aes(x=pr,y=ant.30,color='All Event Pairs'),alpha=1,size=1)+
  geom_point(data=mn,aes(x=pr,y=ant.30,color='Triggering Events'),alpha=1,size=1.5,shape=19)+
  scale_color_manual(values=c("#56B4E9", "#E69F00"))+  theme_pubr()+
  labs(x="Event [mm]",y="Antecedent[mm/day]",color=" ") 
center
dens1 <- ggplot() + 
  geom_density(data=data,aes(x=pr,fill='All Event Pairs'),alpha = .6)+ 
  geom_density(data=mn,aes(x=pr,fill='Triggering Events'),alpha = 0.6)+ 
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+theme_void()+
  theme(legend.position = 'none') + labs(fill=" ")
dens1
dens2 <- ggplot() + 
  geom_density(data=data,aes(x=ant.30,fill='All Event Pairs'),alpha=0.6)+ 
  geom_density(data=mn,aes(x=ant.30,fill='Triggering Events'),alpha = 0.6)+ 
  theme_void()+scale_fill_manual(values=c("#56B4E9", "#E69F00"))+theme_void()+
  theme(legend.position = "none") + labs(fill=" ") + coord_flip()
   
center
              
dens1 + plot_spacer() + center + dens2 + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

#2>>(JJA) Figure:: scatter plot plus density plots-------------


library(ggpubr)
library(patchwork)



center<-ggplot()+
  geom_point(data=jja.data,aes(x=pr,y=ant.30,color='All Event Pairs'),alpha=1,size=1)+
  geom_point(data=jja.mn,aes(x=pr,y=ant.30,color='Triggering Events'),alpha=1,size=1.5,shape=19)+
  scale_color_manual(values=c("#56B4E9", "#E69F00"))+  theme_pubr()+
  labs(x="Event [mm]",y="Antecedent[mm/day]",color=" ") 

dens1 <- ggplot() + 
  geom_density(data=jja.data,aes(x=pr,fill='All Event Pairs'),alpha = .6)+ 
  geom_density(data=jja.mn,aes(x=pr,fill='Triggering Events'),alpha = 0.6)+ 
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+theme_void()+
  theme(legend.position = 'none') + labs(fill=" ")

dens2 <- ggplot() + 
  geom_density(data=jja.data,aes(x=ant.30,fill='All Event Pairs'),alpha=0.6)+ 
  geom_density(data=jja.mn,aes(x=ant.30,fill='Triggering Events'),alpha = 0.6)+ 
  theme_void()+scale_fill_manual(values=c("#56B4E9", "#E69F00"))+
  theme(legend.position = "none") + labs(fill=" ") + coord_flip()
  

center

dens1 + plot_spacer() + center + dens2 + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))



jja.dens1 <- ggplot() + 
  geom_density(data=jja.mn,aes(x=pr,fill='Triggering Events'),alpha = 0.6)+
  geom_density(data=jja.data,aes(x=pr,fill='All Event Pairs'),alpha = 0.6)+ 
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+theme_pubr()+
  theme(legend.position = 'top') + 
  labs(x="Event [mm]",y='Density',fill=" ")

jja.dens1 


jja.dens2 <- ggplot() + 
  geom_density(data=jja.mn,aes(x=ant.30,fill='Triggering Events'),alpha=0.6)+
  geom_density(data=jja.data,aes(x=ant.30,fill='All CMFD Data'),alpha=0.6)+ 
  geom_vline(data=jja.data,aes(xintercept=median(ant.30) ,color='Triggering Events'),linetype="dashed", size=0.5) +
  geom_vline(data=jja.mn,aes(xintercept=median(ant.30) ,color='All CMFD Data'),linetype="dashed", size=0.5) +
  theme_pubr()+ theme(legend.position = 'top') + guides(colour=FALSE)+
  labs(title="Distribution of Antecedent Rainfall within the Study Area",x="30-day Antecedent Rainfall [mm/day]",y='Density',fill=" ")
jja.dens2 

# scale_fill_manual(values=c("#56B4E9", "#E69F00"))
 

jja.dens2
geom_density

#2>>Measure Antecedent Skewness ()--------------------
library(e1071)  

skewness(jja.mn$ant.30)
skewness(jja.data$ant.30)

#2>>Wilcoxon (Rank-Sum Test)
boxplot(jja.mn$ant.30,jja.data$ant.30)

wilcox.test(jja.data$pr, jja.mn$pr, alternative = "two.sided")
wilcox.test(jja.data$pr, jja.mn$pr, mu=0,alt = "two.sided",conf.int=T,conf.level=0.95,paired=F,exact=F,correct=T)M


mean(jja.mn$ant.30)
mean(jja.data$ant.30)
#H0 the median of ant.30 in the inventory = the median of ant.30 in the entire data set for june,july and august.

rs.test<-wilcox.test(jja.mn$ant.30, jja.data$ant.30,alternative = "two.sided",conf.int=T,conf.level=0.95,paired=F,exact=T,correct=T)
rs.test$p.value > 0.05 #Accept H0: mu=0
#Reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
#wilcox.test(jja.data$ant.30,jja.mn$ant.30, mu=0,alt = "two.sided",

#3>>(JJA) Figure:: 2D Histogram Plots ----------
nrow(jja.data)*nrow(mn)




df<-jja.data[,2:3]
colnames(df)<-c('x','y')


jja.his.comp<- ggplot(df, aes(x,y))+ 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",bins =10)+ 
  scale_fill_gradientn(colours=r ,name = "Percent", labels = scales::percent)+
  #scale_fill_gradientn(colours=r ,name = "Prob",values = scales::rescale((c(0,0.0001,0.0005,0.001,0.05))))+
  #geom_point(data=jja.mn,aes(x=pr,y=ant.30,shape='Triggering Events'),alpha=.8,size=1.5,color='black')+
  ggtitle("2D Histogram of CMFD event pairs across the study area\nJune-August (1995-2005)") +
  labs(x="Event [mm]",y="Antecedent[mm/day]",shape=" ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)

jja.his.comp






jja.his.comp<- ggplot(df, aes(x,y))+ 
  #stat_density_2d(aes(fill = ..levels..), geom = "polygon", colour="white")+ 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE,bins=40) +
  geom_point(data=jja.mn,aes(x=pr,y=ant.30,shape='Triggering Events'),alpha=.8,size=1.5,color='black')+
  scale_fill_gradientn(colours=r ,name = "Percent",labels = scales::percent,values = scales::rescale((c(0,0.05))))+
  ggtitle("2D Density Plot of CMFD event pairs across the study area\nJune-August (1995-2005)") +
  labs(x="Event [mm]",y="Antecedent[mm/day]",shape=" ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
jja.his.comp

colnames(df)<-c('x','y')

  
jja.his.comp<- ggplot(data=jja.mn, aes(x=pr,y=ant.30))+ 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",bins=10)+ 
  scale_fill_gradientn(colours=r ,name = "")+
  
  ggtitle("2D Histogram of Reconstructed Triggering Rainfall Events\nJune-August (1995-2005)") +
  labs(x="Event [mm]",y="Antecedent[mm/day]",shape=" ") +
  


jja.his.comp

jja.his.2d<- ggplot(df, aes(x,y))+
        geom_de


jja.his.comp<- ggplot(df, aes(x,y))+ 
  geom_point(data=jja.mn,aes(x=pr,y=ant.30,shape='Triggering Events'),alpha=.8,size=1.5,color='black')+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE,bins=5) +
  
  scale_fill_gradientn(colours=r ,name = "Percent",labels = scales::percent)+
  geom_point(data=jja.mn,aes(x=pr,y=ant.30,shape='Triggering Events'),alpha=.8,size=1.5,color='black')+
  ggtitle("2D Density Plot from Triggering Events in the Inventory\nJune-August (1995-2005)") +
  labs(x="Event [mm]",y="Antecedent[mm/day]",shape=" ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
#theme_pubr()

jja.his.comp


#3>>AMS For JJA Event Rainfall----------
library(xts)
library(zoo)

#3.1. Extract Annual Maximum Series at RG Locations -------------
  
  #3.1.0. Load List Data ------------------------
  
  load('rgs.rdata')
  rgs <- x
  
  
  load('cmfd_rgts.rdata')
  cmfd.1 <- x
  
  #3.1.1. RG TS------------------------
  
  
  r1.rg<-data.frame(rgs$r1$dates,rgs$r1$pr.d)
  r2.rg<-data.frame(rgs$r2$dates,rgs$r2$pr.d)
  names(r1.rg)[1]  <-'dates'
  names(r2.rg)[1]  <-'dates'
  names(r1.rg)[2]  <-'Pr'
  names(r2.rg)[2]  <-'Pr'
  
  # #Analysis of annual rainfall
  # library(zoo)
  # library(xts)
  # df<-r1.rg
  # #ts2ts via xts
  # ts<-xts(df$Pr, order.by=as.POSIXct(df$dates))
  # #apply yearly -> derive max
  # sum <- apply.yearly(ts, sum)
  # ts.yr<- fortify.zoo(sum)
  # ts.yr<-ts.yr[-1,] #First year has incomplete data
  # plot(ts.yr,type='o')
  # mean(ts.yr$sum)
  # #end
  #3.1.2. CMFD TS------------------------
  
  
  r1.cmfd<-data.frame(cmfd.1[[1]][['dates']],cmfd.1[[1]][['pr.d']])
  r2.cmfd<-data.frame(cmfd.1[[2]][['dates']],cmfd.1[[2]][['pr.d']])
  names(r1.cmfd)[1]  <-'dates'
  names(r2.cmfd)[1]  <-'dates'
  names(r1.cmfd)[2]  <-'Pr'
  names(r2.cmfd)[2]  <-'Pr'
  
  #3.1.3. AMS Extraction Function:: Analysis ------------------------
  
  #June-July-August Analysis
  ams.jja <- function(df)
  {
    df<-na.omit(df)
    df<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,] 
    
    #ts2ts via xts
    ts<-xts(df$Pr, order.by=as.POSIXct(df$dates))
    #apply yearly -> derive max
    yr.pr.max <- apply.monthly(ts, max)
    ts.yr<- fortify.zoo(yr.pr.max)
    return(ts.yr)
    }
  ts.ams <- lapply(list(r1.cmfd,r1.rg,r2.cmfd,r2.rg),ams.jja) 
  # >>  1: cmfd, rg; 2: cmfd, rg
    
  #3.1.4. Exceedance Probabilities -----------------
  
  pdf.emp <-function (ams)
  {
    df <- data.frame(c(1,2,5,c=seq(10,200,10)))
    colnames(df)[1] <- 'pr'
    df$p <-0
    df$ne<-0
    df$cp <-0
    df$yr <-0 
    for (i in 1:nrow(df))
    {
      a<-nrow(ams[ams[,2]>=df[i,1],])
      b<-a/nrow(ams)
      df$p[i] = b
      df$ne[i]=1-df$p[i] 
      df$yr[i] =1/b
    }
    
    for (i in nrow(df):2)
    {
      df$cp[i] = df$ne[i]+df$cp[i-1]
    }
    return(df)
  }
  pdf.ams<-lapply(ts.ams,pdf.emp)
  
  plot(pdf.ams[[1]]$yr,pdf.ams[[1]]$pr,type='o')
  lines(pdf.ams[[2]]$yr,pdf.ams[[2]]$pr,type='o',col='red')
  
  plot(pdf.ams[[3]]$pr,pdf.ams[[3]]$cp,type='o')
  lines(pdf.ams[[4]]$pr,pdf.ams[[4]]$cp,type='o',col='red')
  abline(h=.95,col='blue')
  #3.1.5. Fitting Distributions ------------------------
  # library(evir)
  # #GEV
  # gev.fit <- function(df)
  # {
  # fit.gev<-gev(df$yr.pr.max)
  # return(fit.gev)
  # }
  # fit.gev <- lapply(ts.ams,gev.fit) 
  # 
  # #Gumbel
  # gumb.fit<-function(df)
  # {
  #   
  #   fit.gumb<-gumbel(df$yr.pr.max)
  #   return(fit.gumb)
  #   
  #   
  # }
  # fit.gumb <-lapply(ts.ams,gumb.fit) 
  # 
  # library(fitdistrplus)
  # #Weibull
  # weibull.fit<-function(df)
  # {
  #   
  #   fit.weibull<-fitdist(df$yr.pr.max, "weibull")
  #   return(fit.weibull)
  #   
  #   
  # }
  # fit.weibull <-lapply(ts.ams,gumb.fit) 
  
  
  library(gnFit)
  library(ismev)
  
  
  #Dist::CMFD,Comp::RG
  # #http://www-eio.upc.es/teaching/adtl/apunts/TrR_IDA-DF_1.html
  
  df<-ts.ams[[1]]
  # par=data.frame(t(fit.weibull[[1]]$par.ests))
  # 
  # sequence<-seq(0,1,by=0.02)
  # qualist<-quantile(df$yr.pr.max,sequence)
  # sequence;qualist
  # hist(df$yr.pr.max,freq=F,breaks=qualist,main="")
  # m=mean(df$yr.pr.max);std=sd(df$yr.pr.max);m;std
  # curve(dnorm(x,m,std),col="red",lwd=2,add=T)
  # curve(dexp(x,rate=1/m),col="green",lwd=2,add=T)
  # curve(dweibull(x, shape=2.415, scale =77.78))
  # lines(density(df$yr.pr.max),col="blue")
  
  
  #Start:: Define Gumbel Functions
  #a := beta
  #b := alpha
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  #End :: Define Gumbel Functions
  
  
  #Start: Pearson III
  library(e1071)
  m <- mean(df$yr.pr.max)
  v <- var(df$yr.pr.max)
  s <- sd(df$yr.pr.max)
  g <- e1071::skewness(df$yr.pr.max, type=1)
  
  # Correct the sample skew for bias using the recommendation of 
  # Bobee, B. and R. Robitaille (1977). "The use of the Pearson Type 3 and Log Pearson Type 3 distributions revisited." 
  # Water Resources Reseach 13(2): 427-443, as used by Kite
  
  n <- length(df$yr.pr.max)
  g <- g*(sqrt(n*(n-1))/(n-2))*(1+8.5/n)
  
  # We will use method of moment estimates as starting values for the MLE search
  
  my.shape <- (2/g)^2
  my.scale <- sqrt(v)/sqrt(my.shape)*sign(g) # modified as recommended by Carl Schwarz
  my.location <- m-sqrt(v * my.shape)
  
  my.param <- list(shape=my.shape, scale=my.scale, location=my.location)
  
  library(PearsonDS)
  dPIII<-function(x, shape, location, scale) PearsonDS::dpearsonIII(x, shape, location, scale, log=FALSE)
  pPIII<-function(q, shape, location, scale) PearsonDS::ppearsonIII(q, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
  qPIII<-function(p, shape, location, scale) PearsonDS::qpearsonIII(p, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
  
  fitpe<-fitdistrplus::fitdist(df$yr.pr.max, distr="PIII", method="mge", start=my.param)
  library(gsl)
  pearsonFitML(df$yr.pr.max)
  #End:: Pearson III
  
  #df<-ts.ams[[1]]
  library(fitdistrplus)
  fitwe<-fitdist(df$yr.pr.max,"weibull")
  fitln<-fitdist(df$yr.pr.max,"lnorm")
  fitga<-fitdist(df$yr.pr.max,"gamma")
  fitgu<-fitdist(df$yr.pr.max,"gumbel",start=list(a=10, b=10))
  library(actuar)
  fitig<-fitdist(df$yr.pr.max, "invgauss", start = list(mean = 5, shape = 1))
  
  
  #3.1.5. Goodness-of-fit Tests------------------------
  
  gofstat(fitln)
  plot(fitgu,demp = TRUE)
  gof.test<-gofstat(list(fitgu,fitpe))
  
  gof.test$cvmtest
  gof.test$kstest
  gof.test$adtest
  
  cdfcomp(list(fitga,fitgu,fitwe,fitpe,fitig,fitln), legendtext=c("Gamma","Gumbell","Weibull","Pearson III","Inverse Gaussian","Log Normal"))
  denscomp(list(fitga,fitgu,fitwe,fitpe,fitig,fitln), legendtext=c("Gamma","Gumbell","Weibull","Pearson III","Inverse Gaussian","Log Normal"))
  denscomp(list(fitgu,fitpe), legendtext=c("Gumbell","Pearson III"))
  ppcomp(list(fitga,fitgu,fitwe,fitpe,fitig), legendtext=c("Gamma","Gumbell","Weibull","Pearson III","Inverse Gaussian"))
  qqcomp(list(fitga,fitgu,fitwe,fitpe,fitig), legendtext=c("Gamma","Gumbell","Weibull","Pearson III","Inverse Gaussian"))
  
  #3.1.6. Diagnostic Tests -----------
  #Narrowing to Gumbell And Pearson
  #KS, AD, CVM Tests
  library(goftest)
  # library(reliaR)
  # ks.gumbel(df$yr.pr.max,fitgu$estimate[[1]],fitgu$estimate[[2]], alternative = "two.sided", plot = TRUE)
  
  ks.test(df$yr.pr.max,pgumbel,fitgu$estimate[[1]],fitgu$estimate[[2]], alternative = "two.sided")
  ks.test(df$yr.pr.max, pPIII, fitpe$estimate[[1]],fitpe$estimate[[2]], fitpe$estimate[[3]], alternative = "two.sided")
  
  
  cvm.test(df$yr.pr.max, pPIII, fitpe$estimate[[1]],fitpe$estimate[[2]], fitpe$estimate[[3]],estimated=TRUE)
  ad.test(df$yr.pr.max, pPIII, fitpe$estimate[[1]],fitpe$estimate[[2]], fitpe$estimate[[3]],estimated=TRUE)
  
  ks.test(df$yr.pr.max,pgumbel,fitgu$estimate[[1]],fitgu$estimate[[2]], alternative = "two.sided")
  cvm.test(df$yr.pr.max,pgumbel,fitgu$estimate[[1]],fitgu$estimate[[2]])
  ad.test(df$yr.pr.max,pgumbel,fitgu$estimate[[1]],fitgu$estimate[[2]])
  
  
  cdfcomp(list(fitgu,fitpe), legendtext=c("Gumbell","Pearson III"))
  denscomp(list(fitgu,fitpe), legendtext=c("Gumbell","Pearson III"))
  qqcomp(list(fitgu,fitpe), legendtext=c("Gumbell","Pearson III"))
  ppcomp(list(fitgu,fitpe), legendtext=c("Gumbell","Pearson III"))
  # #GEV
  # library('lmomco')
  # 
  # lmom <- lmoms(df$yr.pr.max,nmom=5)     #from lmomco package
  # para <- pargev(lmom, checklmom=TRUE)
  # dgev <- function(x,xi,alpha,kappa)
  #   {
  #   pdfgev(x,list(type="gev",para=c(xi,alpha,kappa),source="pargev"))
  #   }
  # pgev <- cdfgev
  # 
  # para[[2]]
  # 
  # fitgev <- fitdist(df$yr.pr.max, "gev", start=list(xi=para[[2]][[1]],alpha=para[[2]][[2]],kappa=para[[2]][[3]]))        
  # 
  # 
  # #log Pearson III
  # df$yr.pr.max<-log(df$yr.pr.max)
  
  
  
  # 
  # 
  # BIC Stats
  # exp((369.9703-370.1186)/2)
  # exp((369.9703-373.4932)/2)
  # exp((370.1186-373.4932)/2)
  
#3.2. Extracting Annual Maximum Series at ALL Inventory Locations ------------------
  #3.2.0. Load List Data-----------------
  load('mn7_days_cmfd.rdata')
  mn.inv <-x
  # mn<- rbindlist(lapply(mn.inv, `[[`, 'pa.events') )
  # mn<-data.frame(mn$dates,mn$pr.d,mn$ant.30)
  # colnames(mn)<-c('dates','pr','ant.30')
  
  #3.2.1. Extracting AMS---------------
  library(xts)
  library(zoo)
  tr.cmfd<-Map(c,lapply(mn.inv, '[', 'ts'), lapply(mn.inv, '[', 'lon_cmfd'),lapply(mn.inv, '[', 'lat_cmfd'))
  #Expanded ams.jja function
  list<-tr.cmfd
  i=1
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
    ts<-xts(df$pr.d, order.by=as.POSIXct(df$dates))
    #apply yearly -> derive max
    
    yr.pr.max <- apply.monthly(ts, max)
    yr.pr.mean <- apply.monthly(ts, mean)
    
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
  list.save(tr.cmfd,'jja.tr.cmfd.rdata')

  #3.2.1. Fitting Distributions (JJA)---------------
        
        load('jja.tr.cmfd.rdata')
        jja.tr.cmfd<-x
  
      #3.2.1.1. Gumbel Fitting Functions---------------
      
      library(fitdistrplus)

        #Start:: Define Gumbel Functions
        #a := beta
        #b := alpha
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
      
      list.save(jja.tr.cmfd,'jja.tr.cmfd.rdata')

      #3.2.1.3. Gumbel Goodness-of-fit Tests--------------
      
      library(goftest)
      load('jja.tr.cmfd.rdata')
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
      list.save(jja.tr.cmfd,'jja.tr.cmfd.rdata')
      
    
      
  #3.2.2. Inspecting Gumbel Distributions--------------
        load('jja.tr.cmfd.rdata')
        
        
      
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
        load('jja.tr.cmfd.rdata')
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

      #3.2.3.2 Return Periods for all Rainfall-------------
      
      load('mn7_days_cmfd.rdata')
      mn.inv <-x
      
      load('jja.tr.cmfd.rdata')
      jja.tr.cmfd<-x
      # map.gum.tr<-Map(c, lapply(mn.inv, '[[', 'events'),lapply(mn.inv, '[[', 'lon_cmfd'),lapply(mn.inv, '[[', 'lat_cmfd'))
      map.gum.tr<-Map(c,lapply(mn.inv, '[[', 'lon_cmfd'),lapply(mn.inv, '[[', 'lat_cmfd'))
      
      
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
      
      list.save(jja.tr.cmfd,'jja.tr.cmfd.rdata')
      list.save(map.gum.tr,'jja.map.gum.rdata')
      
      #3.2.3.2 Mapping Return Periods-------------
      load('jja.tr.cmfd.rdata')
      jja.tr.cmfd<-x
      
      load('jja.map.gum.rdata')
      jja.map.gum<-x
      
      #Expand tr.gum dataframe for all coordinates
      list<-jja.map.gum

      #Function to extract map data.frames with return period and coordinates
      #t.pr -> return period of your mapped rainfall
      #list -> seasonally-filtered distribution fitted - map list 
      map.df.extract<-function(t.pr,list)
      {
        for(i in 1:length(list))
        {
          df.map <- data.frame(list[[i]][['x']],list[[i]][['y']])
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
      
      
      #50-Year Return Period 
      list.map<-map.df.extract(50,jja.map.gum)
      df.50.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
      
      list.map<-map.df.extract(40,jja.map.gum)
      df.40.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
      
      
      
      
      #Import Polygons for China, select region and clip
      
      china3 <- readRDS('gadm36_CHN_3_sp.rds')
      wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
      CP <- as(extent(106,109, 30.0 ,35.0), "SpatialPolygons")
      proj4string(CP) <- proj4string(wanzhou)
      ## Clip the map
      library(rgeos)
      wanzhou <- gIntersection(wanzhou, CP, byid=TRUE)
      wanzhou_df<-fortify(wanzhou)
      
      
      
      
      
      
      load('mn7_days_cmfd.rdata')
      mn.inv <-x
      df.inv<-(lapply( mn.inv,`[[`,'events'))
      df.inv.cor<-data.frame(rbindlist(lapply(df.inv,`[`,'lon_c')),rbindlist(lapply(df.inv,`[`,'lat_c')))
      
      
      
      
      
      
      fortify(df.20.map)
      
      
      df.20<-ggplot(df.20.map)+
        geom_point(data= df.20.map,aes(x=lon,y=lat),pch=5)+
        geom_point(data= df.inv.cor,aes(x=lon_c,y=lat_c))+
      
      #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
      scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
        labs(title="Pixels-Inventory Points" ,x="Longitude",y="Latitude",color="Legend") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      df.20
      
      
      
      
      #convert from data frame to spatial points data frame
      library(sp)
      
      coordinates(df.20.map)<- 1:2
      coordinates(df.40.map)<- 1:2
      coordinates(df.50.map)<- 1:2
      
      #IDW interpolation Script
      library(gstat) # Use gstat's idw routine
      library(sp)    # Used for the spsample function
      
      #loading CMFD raster files    
      files <- list()
      files = list.files('/Volumes/DAT/2_code/data',pattern='*.nc',full.names=TRUE)
      
      a<-nc_open(files[3])
      #raster handling https://rpubs.com/markpayne/358146
      cmfd.1 <- raster::brick(files[3])
      
      #crop ciles around area-of-interest using extent
      aoi <- extent(107,110,30,32)
      cmfd.crop <- crop(cmfd.1[[1]],aoi)
     
      #Importing Grid fom CMFD
      
      grd.cmfd<-as.data.frame(cmfd.crop[[1]],xy=T)
      grd              <- grd.cmfd[,-3]
      names(grd)       <- c("X", "Y")
      coordinates(grd) <- c("X", "Y")
      gridded(grd)     <- TRUE  # Create SpatialPixel object
      fullgrid(grd)    <- TRUE  # Create SpatialGrid object
      
      plot(grd)
      
      # # Add P's projection information to the empty grid
      # 
      proj4string(grd) <- proj4string(wanzhou)
      
      idw.function<-function(grd,P)
      {
      proj4string(P) <- proj4string(grd) # Temp fix until new proj env is adopted
      #proj4string(grd) <- proj4string(wanzhou)
      # Interpolate the grid cells using a power value of 2 (idp=2.0)
      P.idw <- gstat::idw(pr ~ 1, P, newdata=grd, idp=2.0)
      # Convert to raster object then clip to Texas
      r       <- raster(P.idw)
      #r.m     <- mask(r, wanzhou)
      return(r)
      }
      
      names(ant.30.map)<-c("lon",'lat','pr')
      coordinates(ant.30.map)<- 1:2
      r.m<-idw.function(grd,ant.30.map)
      writeRaster(r.m.c, filename="antpr_uc.tif", format="GTiff", overwrite=TRUE,proj4string(wanzhou))
      
      
      r.m<-idw.function(grd,df.50.map)
      r.m.c<-crop(r.m,CP)
      proj4string(r.m.c) <- CRS(proj4string(wanzhou))
      #writeRaster(r.m.c, filename="rt50.tif", format="GTiff", overwrite=TRUE,proj4string(wanzhou))
      #writeRaster(r.m.c, filename="rt50_uc.tif", format="GTiff", overwrite=TRUE,proj4string(wanzhou))
      ras<-na.omit(fortify(r.m.c))
      
      # map.t.50<-ggplot(ras,aes(x=x,y=y,z=var1.pred))+geom_contour_filled()
      
      
      
      
      
      map.t.50<-ggplot(ras)+
        geom_raster(data=ras,aes(x=x,y=y,fill=var1.pred))+
       
        #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
        scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
        labs(title="50-Year Return Period Rainfall" ,x="Longitude",y="Latitude",color="Legend") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      map.t.50
      
      r.m<-idw.function(grd,df.20.map)
      r.m.c<-crop(r.m,CP)
      proj4string(r.m.c) <- CRS(proj4string(wanzhou))
      #writeRaster(r.m.c, filename="rt20.tif", format="GTiff", overwrite=TRUE)
      writeRaster(r.m.c, filename="rt20_uc.tif", format="GTiff", overwrite=TRUE)
      ras<-na.omit(fortify(r.m))
      
      map.t.20<-ggplot(ras)+
        geom_raster(data=ras,aes(x=x,y=y,fill=var1.pred))+
        #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
        scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.9,108.9), ylim = c(30.3,31.05)) + 
        labs(title="20-Year Return Period Rainfall" ,x="Longitude",y="Latitude",color="Legend") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      map.t.20
      
      r.m<-idw.function(grd,df.40.map)
      r.m.c<-crop(r.m,CP)
      proj4string(r.m.c) <- CRS(proj4string(wanzhou))
      #writeRaster(r.m.c, filename="rt40.tif", format="GTiff", overwrite=TRUE)
      writeRaster(r.m.c, filename="rt40_uc.tif", format="GTiff", overwrite=TRUE)
      ras<-na.omit(fortify(r.m))
      map.t.40<-ggplot(ras)+
        geom_raster(data=ras,aes(x=x,y=y,fill=var1.pred))+
        #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
        scale_fill_gradientn(colours=colssd,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
        labs(x="Longitude",y="Latitude",color="Legend") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      map.t.40 
      
      
      map.t.20
      map.t.40
      map.t.50
      
      
      t.r=200
      list.map<-map.df.extract(t.r,jja.map.gum)
      df.map<-rbindlist(lapply(list.map,`[[`,'df.map'))
      coordinates(df.map)<- 1:2
      r.m<-idw.function(wanzhou,df.map)
      ras<-na.omit(fortify(r.m))
      map.t0<-ggplot(ras)+
        geom_raster(data=ras,aes(x=x,y=y,fill=var1.pred))+
        scale_fill_gradientn(colours=r,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
        labs(x="Longitude",y="Latitude",color="Legend") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      map.t0
#3.3. Extracting Mean Antecedent Rainfall at ALL Inventory Locations----------
      
      load('mn7_days_cmfd.rdata')
      mn.inv <-x
      max(mn.inv[[1]][['ts']][['dates']])
      list<-mn.inv
      library(xts)
      library(zoo)

      ant.mean.jja <- function(list)
      {
        for(i in 1:length(list))
        {
          ts<-list[[i]][['ts']]
          ant.ts<-list[[i]][['ant.ts']]
          df <- data.frame(ts$dates,ant.ts$ant.30)
          colnames(df)<-c('dates','ant.30')
          df<-na.omit(df)
          df<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,] 
          
          #ts2ts via xts
          ts<-xts(df$ant.30, order.by=as.POSIXct(df$dates))
          #apply yearly -> derive max
          yr.pr.mean <- apply.monthly(ts, mean)
          ts.yr<- fortify.zoo(yr.pr.mean)
          
          ypr.mean<-data.frame(ts.yr$yr.pr.mean/30)
          names(ypr.mean)[1] <- 'ant.30'
          ypr.mean <- list(ypr.mean)
          ypr.mean <- list.names(ypr.mean,"ypr.mean")
          list[[i]] <- append(list[[i]],ypr.mean)
        }
        return(list)
      }
      jja.ant.mean<-ant.mean.jja(mn.inv)
      #jja.ant.mean<-(lapply(list, `[`, 'ypr.mean'))
      map.df.extract<-function(list)
      {
        for(i in 1:length(list))
        {
          df.map <- data.frame(list[[i]][['lon_cmfd']],list[[i]][['lat_cmfd']])
          colnames(df.map)<-c('lon','lat')
          pr.mean<-mean(list[[i]][['ypr.mean']][['ant.30']])#potential bug flag::one list for each distribution and season
          df.map$ant.30=pr.mean
          df.map<-list(df.map)
          df.map<-list.names(df.map,'ant.30.map')
          list[[i]]<-append(list[[i]],df.map)
        }
        list<-lapply(list, `[`, 'ant.30.map')
        return(list)
      }
      
      list.map<-map.df.extract(jja.ant.mean)
      ant.30.map<-rbindlist(lapply(list.map,`[[`,'ant.30.map'))
      write.csv(ant.30.map,file="ant_30.xyz")
      ras<-rasterFromXYZ(ant.30.map)
    
      #Import Polygons for China, select region and clip
      library(rgeos)
      china3 <- readRDS('gadm36_CHN_3_sp.rds')
      wanzhou<-china3[china3@data$NAME_3 == 'Wan',]
      CP <- as(extent(107.8,109, 30.38,31.05), "SpatialPolygons")
      proj4string(CP) <- CRS(proj4string(wanzhou))
      ## Clip the map
      wanzhou <- gIntersection(wanzhou, CP, byid=TRUE)
      library(sp)
      proj4string(ras) <- CRS(proj4string(wanzhou))
      wanzhou_df<-fortify(wanzhou)
      ras<-crop(ras,wanzhou)
      
      #writeRaster(ras, filename="ant.30.tif", format="GTiff", overwrite=TRUE)
  

      df.ant.map<-na.omit(as.data.frame(ras))
      map.ant<-ggplot()+
        geom_raster(data=ras,aes(x=x,y=y,fill=ant.30))+
        #geom_contour_filled(data=df.20.map,aes(x=lon,y=lat,z=pr),binwidth=.5)+
        #scale_fill_gradientn(colours=r,  guide = "colourbar",name = "Rainfall [mm/day]")+
        geom_polygon(data=wanzhou_df, aes(long,lat, group=group), fill=NA, color="black",size=1.5)+
        coord_sf(xlim = c(107.8,109), ylim = c(30.3,31.1)) + 
        labs(title="Average 30-day Antecedent Rainfall from June-August" ,x="Longitude",y="Latitude",fill="[mm/day]") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
      map.ant
      
#4.1 Calculating Return Period each event date--------------
      
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
      coor<-select(coor,c(1:3))
      colnames(coor)[1]<-'recdate'
      coor<-arrange(coor,recdate)
      mn<-arrange(mn,recdate)
      
      export<-cbind(mn[,1],mn[,3],coor[,2:3],mn[,2],mn[,4:8])
      colnames(export)[1]<-'actual_date'
      colnames(export)[2]<-'recorded_date'
      
      ƒ
      
      
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
      
      
#End of Script---------