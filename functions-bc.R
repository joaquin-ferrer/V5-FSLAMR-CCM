#===================================================================#
# VERSION 1 EXTREME EDCDF FUNCTIONS
#===================================================================#




#===================================================================#
#Function 1: BILINEAR INTERPOLATION ----
#===================================================================#


#===================================================================#
#Function 2: AMS EXTRACTION ----
#===================================================================#

extract.ams<-function(df.pr,df.dates)
#INPUT DATAFRAME WITH PRECIPITATION AND DATAFRAME WITH DATES
{
  require(xts)
  require(zoo)
  ts<-xts(df.pr, order.by=as.POSIXct(df.dates))
  yr.pr.max <- apply.monthly(ts, max)
  ts.yr<- fortify.zoo(yr.pr.max)
  ypr.max<-data.frame(ts.yr$yr.pr.max)
  colnames(ypr.max)<-c('pr')
  return(ypr.max)
}
#===================================================================#
#Function 3: GUMBEL FIT DISTRIBUTION FUNCTION ----
#===================================================================#

#Define Gumbel Functions
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))


fit.gumb<-function(df.pr)
#INPUT DATAFRAME WITH MAXIMUM PR
{
  require(fitdistrplus)
  require(rlist)
  df<-as.matrix(na.omit(df.pr))
  #Method == Quantle Matching
  fitgu<-fitdist(as.vector(df),"gumbel",start=list(a=10, b=10),method='qme',probs=c(1/3,2/3))
  fit.gum<-list(fitgu)
  fit.gum<<- list.names(fit.gum,"fit.gum")
  return(fit.gum)
}

#===================================================================#
#CLASS 4:PARAMETRIC BIAS CORRECTION FUNCTIONS----
#===================================================================#


edcdf<-function(x.mod,par.obs,par.mod,par.proj)
#INPUT PARAMETERS ARE THE MODEL RESULTS
#GUMBEL PARAMETERS FOR OBSERVATION AND MODEL
#GUMBEL PARAMETER ESTIMATES:: ALPHA AND BETA
{
  a1<-pgumbel(x.mod,par.proj[1],par.proj[2])
  b1<-qgumbel(a1,par.obs[1],par.obs[2])
  b2<-qgumbel(a1,par.mod[1],par.mod[2])
  #c<-cbind(x.mod,a1,b1,b2)
  x.cor<-x.mod*b1/b2
  return(x.cor)

}


pqcdf<-function(x.mod,par.obs,par.mod,par.proj)
  #INPUT PARAMETERS ARE THE MODEL RESULTS
  #GUMBEL PARAMETERS FOR OBSERVATION AND MODEL
  #GUMBEL PARAMETER ESTIMATES:: ALPHA AND BETA
{
  a1<-pgumbel(x.mod,par.proj[1],par.proj[2])
  b1<-qgumbel(a1,par.obs[1],par.obs[2])
  #b2<-qgumbel(a1,par.mod[1],par.mod[2])
  #c<-cbind(x.mod,a1,b1,b2)
  #x.cor<-x.mod*b1BA
  x.cor<-b1
  return(x.cor)
  
}


pdqm<-function(x.mod,par.obs,par.mod,par.proj)
  #INPUT PARAMETERS ARE THE MODEL RESULTS
  #GUMBEL PARAMETERS FOR OBSERVATION AND MODEL
  #GUMBEL PARAMETER ESTIMATES:: ALPHA AND BETA
{
  a1<-pgumbel(x.mod,par.proj[1],par.proj[2])
  b1<-qgumbel(a1,par.obs[1],par.obs[2])
  b2<-qgumbel(a1,par.mod[1],par.mod[2])
  x.cor<-b1*x.mod/b2
  return(x.cor)
  
}



#===================================================================#
#CLASS 5:NON-ANRAMETRIC BIAS CORRECTION FUNCTIONS----
#===================================================================#
# 
# x.mod<-proj$pr.d
# obs.df<-obs$pr.d
# mod.df<-mod$pr.d

empqdm<-function(x.mod,obs.df,mod.df)
{
  #===================!
  ecdf.proj<-ecdf(x.mod)
  a1<-ecdf.proj(x.mod)
  b1<-quantile(obs.df,a1)
  b2<-quantile(mod.df,a1)
  x.cor<-b1*x.mod/b2
  # plot(x.cor)
  # # 
  # ana<-cbind(b1,b2,x.mod,x.cor)
  # 
  return(x.cor)
  
}

empqm<-function(x.mod,obs.df)
{
  #===================!

ecdf.mod<-ecdf(x.mod)
a1<-ecdf.mod(x.mod)
x.cor<-quantile(obs.df,a1)
return(x.cor)

}

