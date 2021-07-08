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

#1.0 Load List Data-----------------
load('mn7_days_cmfd.rdata')
mn.inv <-x

load('refrac.rdata')
refrac <- data.frame(x)
refrac<-unique(refrac)

#2.0 Creating New Lists-----------------
#refrac and 30-day ant for 1994-2005
for(i in 1:length(mn.inv))
{
  a <- refrac
  names(a)[1]  <-'dates'
  b <- data.frame(mn.inv[[i]][['ts']][['dates']],mn.inv[[i]][['ant.ts']][['ant.30']])
  names(b)[1]  <-'dates'
  names(b)[2]  <-'ant.30'
  ant.30<-list(b)
  ant.30<-list.names(ant.30,'ant.30.ts')
  mn.inv[[i]]<-append(mn.inv[[i]],ant.30)
  b<- b[as.numeric(strftime(b$dates, "%Y")) %in% 1994:2005,] 
  c <- merge(transform(a, dates = format(as.Date(dates), "%Y-%m-%d")), transform(b, dates = format(as.Date(dates), "%Y-%m-%d")))
  qa <- data.frame(a$dates,c$refrac*c$ant.30/30)
  names(qa)[1]  <-'dates'
  names(qa)[2]  <-'qa.30' 
  qa<-list(qa)
  qa<-list.names(qa,'qa.ts')
  mn.inv[[i]]<-append(mn.inv[[i]],qa)
}



antqa.cmfd<-Map(c,lapply(mn.inv, '[', 'ant.30.ts'),lapply(mn.inv, '[', 'qa.ts'), lapply(mn.inv, '[', 'lon_cmfd'),lapply(mn.inv, '[', 'lat_cmfd'))
list.save(antqa.cmfd,'jja.antqa.cmfd.rdata')
#3.0 Fitting Distributions (JJA) LOG-NORM-----------------
load('jja.antqa.cmfd.rdata')
jja.antqa.cmfd<-x
library(fitdistrplus)
for(i in 1:length(jja.antqa.cmfd))
{
  df<-jja.antqa.cmfd[[i]][['ant.30.ts']]
  df<-df[df$ant.30!=0,]
  df<-na.omit(df)
  fit.lnorm<-fitdist(df$ant.30,"lnorm")
  fit.lnorm<-list(fit.lnorm)
  fit.lnorm<<- list.names(fit.lnorm,"fit.lnorm")
  jja.antqa.cmfd[[i]] <- append(jja.antqa.cmfd[[i]],fit.lnorm)
}


library(goftest)
library(stats)

dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
qlnorm(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
rlnorm(n, meanlog = 0, sdlog = 1)

for(i in 1:length(jja.antqa.cmfd))
{
  fit<-jja.antqa.cmfd[[i]][['fit.lnorm']]
  df<-jja.antqa.cmfd[[i]][['ant.30.ts']]
  df<-df[df$ant.30!=0,]
  df<-na.omit(df)#flag for potential bugs
  #Kolmogorov–Smirnov Goodness-of-fit test
  test.ks <- ks.test(df$ant.30,"plnorm",fit$estimate[[1]],fit$estimate[[2]], alternative = "two.sided")
  res.ks <- if(test.ks$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ks <- list(res.ks,test.ks)
  test.ks<-list.names(test.ks,'test.ks')
  names(test.ks) <- c("h.test","ks.stats")
  
  #Anderson–Darling test Goodness-of-fit test
  test.ad <- ad.test(df$ant.30,"plnorm",fit$estimate[[1]],fit$estimate[[2]])
  res.ad <- if(test.ad$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ad <- list(res.ad,test.ad)
  test.ad<-list.names(test.ad,'test.ad') 
  names(test.ad) <- c("h.test","ad.stats")
  
  #Cramer-von Mises Goodness-of-fit test
  test.cvm <- cvm.test(df$ant.30,"plnorm",fit$estimate[[1]],fit$estimate[[2]])
  res.cvm <- if(test.cvm$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.cvm <- list(res.cvm,test.cvm)
  test.cvm<-list.names(test.cvm,'test.cvm')
  names(test.cvm) <- c("h.test","cvm.stats")
  
  #test compilation and integration
  tests.lnorm<-list(test.ks,test.ad,test.cvm)
  names(tests.lnorm) <- c('test.ks','test.ad','test.cvm')
  tests.lnorm<-list(tests.lnorm)
  tests.lnorm<-list.names(tests.lnorm,'tests.lnorm')
  jja.antqa.cmfd[[i]] <- append(jja.antqa.cmfd[[i]],tests.lnorm)
  
}
x<-jja.antqa.cmfd

#3.0 Fitting Distributions (JJA) GAMMA-----------------
load('jja.antqa.cmfd.rdata')
jja.antqa.cmfd<-x
library(fitdistrplus)
for(i in 1:length(jja.antqa.cmfd))
{
  df<-jja.antqa.cmfd[[i]][['ant.30.ts']]
  df<-df[df$ant.30!=0,]
  df<-na.omit(df)
  fit.gamma<-fitdist(df$ant.30,"gamma")
  fit.gamma<-list(fit.gamma)
  fit.gamma<<- list.names(fit.gamma,"fit.gamma")
  jja.antqa.cmfd[[i]] <- append(jja.antqa.cmfd[[i]],fit.gamma)
}

library(goftest)
library(stats)
# 
# dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
# plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
# qlnorm(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
# rlnorm(n, meanlog = 0, sdlog = 1)

for(i in 1:length(jja.antqa.cmfd))
{
  fit<-jja.antqa.cmfd[[i]][['fit.gamma']]
  df<-jja.antqa.cmfd[[i]][['ant.30.ts']]
  df<-df[df$ant.30!=0,]
  df<-na.omit(df)#flag for potential bugs
  #Kolmogorov–Smirnov Goodness-of-fit test
  test.ks <- ks.test(df$ant.30,"pgamma",fit$estimate[[1]],fit$estimate[[2]], alternative = "two.sided")
  res.ks <- if(test.ks$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ks <- list(res.ks,test.ks)
  test.ks<-list.names(test.ks,'test.ks')
  names(test.ks) <- c("h.test","ks.stats")
  
  #Anderson–Darling test Goodness-of-fit test
  test.ad <- ad.test(df$ant.30,"pgamma",fit$estimate[[1]],fit$estimate[[2]])
  res.ad <- if(test.ad$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.ad <- list(res.ad,test.ad)
  test.ad<-list.names(test.ad,'test.ad') 
  names(test.ad) <- c("h.test","ad.stats")
  
  #Cramer-von Mises Goodness-of-fit test
  test.cvm <- cvm.test(df$ant.30,"pgamma",fit$estimate[[1]],fit$estimate[[2]])
  res.cvm <- if(test.cvm$p.value>=0.05) {'Accept'} else{'Reject'} # H-test at 0.05 confidence interval
  test.cvm <- list(res.cvm,test.cvm)
  test.cvm<-list.names(test.cvm,'test.cvm')
  names(test.cvm) <- c("h.test","cvm.stats")
  
  #test compilation and integration
  tests.gamma<-list(test.ks,test.ad,test.cvm)
  names(tests.gamma) <- c('test.ks','test.ad','test.cvm')
  tests.gamma<-list(tests.gamma)
  tests.gamma<-list.names(tests.gamma,'tests.gamma')
  jja.antqa.cmfd[[i]] <- append(jja.antqa.cmfd[[i]],tests.gamma)
  
}
x<-jja.antqa.cmfd
list.save(x,'jja.ant.gamma.cmfd.rdata')





#-------------------------#
#inspection
load('jja.ant.gamma.cmfd.rdata')


# tg<-lapply(x, `[[`, 'tests.gamma')
tg<-lapply(x, `[[`, 'tests.lnorm')
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

#2>>Measure Antecedent Skewness  ---------





ant.30<-antqa.cmfd[[1]][['ant.30.ts']]
#antecedent 
ant.30<- ant.30[as.numeric(strftime(ant.30$dates, "%m")) %in% 6:8,] 

library(e1071)  

skewness(ant.30$ant.30)
 
skewness(sqrt(ant.30$ant.30))
hist((ant.30$ant.30))
plot(density(ant.30$ant.30))

qa.30<-antqa.cmfd[[1]][['qa.ts']]

qa.30<- qa.30[as.numeric(strftime(qa.30$dates, "%m")) %in% 6:8,] 
hist((qa.30$qa.30))
plot(density(qa.30$qa.30))

skewness(qa.30$qa.30)



#3>>Exploring Distributions for ant.30---------

library(fitdistrplus)
fitwe<-fitdist(ant.30$ant.30,"weibull")
fitln<-fitdist(ant.30$ant.30,"lnorm")
fitga<-fitdist(ant.30$ant.30,"gauss")

fitn<-fitdist(ant.30$ant.30,"norm")


#Start: Pearson III
library(e1071)
m <- mean(ant.30$ant.30)
v <- var(ant.30$ant.30)
s <- sd(ant.30$ant.30)
g <- e1071::skewness(ant.30$ant.30, type=1)

# Correct the sample skew for bias using the recommendation of 
# Bobee, B. and R. Robitaille (1977). "The use of the Pearson Type 3 and Log Pearson Type 3 distributions revisited." 
# Water Resources Reseach 13(2): 427-443, as used by Kite

n <- length(ant.30$ant.30)
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

fitpe3<-fitdistrplus::fitdist(ant.30$ant.30, distr="PIII", method="mge", start=my.param)
library(gsl)
pearsonFitML(ant.30$ant.30)
#End:: Pearson III

plot(fitpe3) # Q-Q plot didn't work so well


#Start:: Define Gumbel Functions
#a := beta
#b := alpha
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
#End :: Define Gumbel Functions
fitgu<-fitdist(ant.30$ant.30,"gumbel",start=list(a=10, b=10))

plot(fitgu)


plot(fitln)


denscomp(list(fitwe,fitln,fitga), legendtext=c("Weibull","Log-Normal","Gamma"))
ppcomp(list(fitwe,fitln,fitga), legendtext=c("Weibull","Log-Normal","Gamma"))
qqcomp(list(fitwe,fitln,fitga), legendtext=c("Weibull","Log-Normal","Gamma"))




# 
# fitwe<-fitdist(sqrt(ant.30$ant.30),"weibull")
# fitln<-fitdist(sqrt(ant.30$ant.30),"lnorm")
# fitga<-fitdist(sqrt(ant.30$ant.30),"gamma")
# fitnorm<-fitdist(sqrt(ant.30$ant.30),"norm")
# 
# cdfcomp(list(fitwe,fitln,fitga,fitnorm), legendtext=c("Weibull","Log-Normal","Gamma","Normal"))
# qqcomp(list(fitwe,fitln,fitga,fitnorm), legendtext=c("Weibull","Log-Normal","Gamma","Normal"))





#4>>Log-Normal Fitting for ant.30 and GOF tests ---------

#fit
library(fitdistrplus)
fitln<-fitdist(ant.30$ant.30,"lnorm")
plot(fitln)
#gof tests
gofstat(fitln)


#5>>>Exploring Distributions for qa.30 and GOF tests ---------

qa.30[qa.30==0] <- NA
qa.30<-na.omit(qa.30)
library(fitdistrplus)
fitwe<-fitdist(qa.30$qa.30,"weibull")
fitln<-fitdist(qa.30$qa.30,"lnorm")
fitga<-fitdist(qa.30$qa.30,"gamma")



#Start:: Define Gumbel Functions
#a := beta
#b := alpha
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
#End :: Define Gumbel Functions
fitgu<-fitdist(qa.30$qa.30,"gumbel",start=list(a=10, b=10))




denscomp(list(fitwe,fitln,fitga,fitgu), legendtext=c("Weibull","Log-Normal","Gamma","Gumbel"))
ppcomp(list(fitwe,fitln,fitga,fitgu), legendtext=c("Weibull","Log-Normal","Gamma","Gumbel"))
qqcomp(list(fitwe,fitln,fitga,fitgu), legendtext=c("Weibull","Log-Normal","Gamma","Gumbel"))

plot(fitgu)

gofstat(fitgu)

