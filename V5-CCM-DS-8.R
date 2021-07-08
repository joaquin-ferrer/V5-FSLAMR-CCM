#===================================================================#
# VERSION 8 OF THE BIAS CORRECTION AND DOWNSCALING METHODOLOGY
#===================================================================#
require(ggplot2)
require(ggpubr)
require(rlist)
source('extreme-bc.R')
source('biasCorrection.R')
#!==========================!  
#Dat file handle
g=1
#!==========================!  
#RUN VALIDATION AGAIN
nval<-10
#!==========================!
# STEP 0: Data Handling ----------

dat.list <- list()
dat.list = list.files('/Volumes/DAT/2_code/V5-FSLAMR/bc-dat',pattern='*.rdata',full.names=TRUE)

# STEP 1A:: DAY BIAS CORRECTIONS: EMPIRICAL QDM--------
source('biasCorrection.R')
require('qmap')
require('dplyr')

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

#June-July-August Filter Functions
jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  return(df.jja)
}
mjja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 5:8,]
  return(df.jja)
}
mcs.jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  df.jja<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2021:2060,]
  return(df.jja)
}
lcs.jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  df.jja<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2061:2010,]
  return(df.jja)
}

months<-function(df)
#summarizes and returnes monthly means
#input:: df 'dates' and 'pr.d'
#output:: numeric months, monthly precipitation, character abbreviations of months
{
names(df)[1]<-'dates'
require('dplyr')
df.m <-mutate(df,month = format(as.Date(dates), "%m")) %>% group_by(month)
df.m <- summarise(df.m,pr.mo = mean(pr.d,na.rm = TRUE))
df.m$chmo <- as.numeric(df.m$month)
df.m<- transform(df.m, chmo = month.abb[chmo])
df.m$month <- as.numeric(df.m$month)
return(df.m)
}

dayref.bc.list<-Map(c,lapply(ts.list, '[', 'i'), lapply(ts.list, '[', 'j'),lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat'))
for (r in 1:length(dayref.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]

  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)

  obs.df<-dplyr::select(cmfd.df,c(2:1))

  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)
  
  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  for (j in 1:nval)
  {
    metrics.df<-data.frame(mae=0,rmse=0,bias=0,cor=0)
    #Training and Testing Data Split::
    test.fold<-createFolds(obs$pr.d, k = 3, list = TRUE, returnTrain = FALSE)
    
    #================K=1================!
    #Extract Test Data::
    obs.test<-obs[test.fold[[1]],]  
    mod.test<-mod[test.fold[[1]],]  
    #Extract Train Data::
    obs.train<-rbind(obs[test.fold[[2]],],obs[test.fold[[3]],])
    mod.train<-rbind(mod[test.fold[[2]],],mod[test.fold[[3]],])
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE
    rmse<-rmse(comp.df$bc,comp.df$pr.d)
    #Bias
    bias<-bias(comp.df$bc,comp.df$pr.d)
    #Correlation
    cor<-cor(comp.df$bc,comp.df$pr.d)
    temp<-data.frame(mae,rmse,bias,cor)
    metrics.df<-rbind(metrics.df,temp)
    #================K=1================!
    
    #================K=2================!
    #Extract Test Data::
    obs.test<-obs[test.fold[[2]],]  
    mod.test<-mod[test.fold[[2]],]  
    #Extract Train Data::
    obs.train<-rbind(obs[test.fold[[1]],],obs[test.fold[[3]],])
    mod.train<-rbind(mod[test.fold[[1]],],mod[test.fold[[3]],])
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE
    rmse<-rmse(comp.df$bc,comp.df$pr.d)
    #Bias
    bias<-bias(comp.df$bc,comp.df$pr.d)
    #Correlation
    cor<-cor(comp.df$bc,comp.df$pr.d)
    temp<-data.frame(mae,rmse,bias,cor)
    metrics.df<-rbind(metrics.df,temp)
    #================K=2================!
    
    
    #================K=3================!
    #Extract Test Data::
    obs.test<-obs[test.fold[[3]],]  
    mod.test<-mod[test.fold[[3]],]  
    #Extract Train Data::
    obs.train<-rbind(obs[test.fold[[1]],],obs[test.fold[[2]],])
    mod.train<-rbind(mod[test.fold[[1]],],mod[test.fold[[2]],])
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE
    rmse<-rmse(comp.df$bc,comp.df$pr.d)
    #Bias
    bias<-bias(comp.df$bc,comp.df$pr.d)
    #Correlation
    cor<-cor(comp.df$bc,comp.df$pr.d)
    temp<-data.frame(mae,rmse,bias,cor)
    metrics.df<-rbind(metrics.df,temp)
    #================K=3================!
    
    #================performance summary================!
    metrics.df<-metrics.df[-1,]
    mean.temp<-data.frame(mean(metrics.df$mae),mean(metrics.df$rmse),mean(metrics.df$bias),mean(metrics.df$cor))
    colnames(mean.temp)<-c('mae','rmse','bias','cor')
    mean.met<-rbind(mean.met,mean.temp)
    
  }
  #============k-fold cross-validation ===========!
  # df<-obs.df
  # df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  # obs<-na.omit(df)
  # #obs<-obs[obs$pr.d>=1,]
  # obs<-mjja(obs)
  # 
  # df<-mod.df
  # df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  # mod<-na.omit(df)
  # #mod<-mod[mod$pr.d>=1,]
  # mod<-mjja(mod)
  # 
  # df<-proj.df
  # df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  # proj<-na.omit(df)
  # #proj<-proj[proj$pr.d>=1,]
  # proj<-mjja(proj)
  # 
  # 
  # df<-obs.df
  # df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  # val<-na.omit(df)
  # #val<-val[val$pr.d>=1,]
  # val<-mjja(val)
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)

  bc<-empqm(mod$pr.d,obs$pr.d)
  rfit.list<-list(bc,mod,obs)
  names(rfit.list)<-c('bc','proj','obs')
  rfit.list<-list(rfit.list)
  rfit.list<-list.names(rfit.list,'ref.bc')
  dayref.bc.list[[r]]<-append(dayref.bc.list[[r]],rfit.list)
  metrics<-list(mean.met)
  metrics<-list.names(metrics,'crval.3fold')
  dayref.bc.list[[r]]<-append(dayref.bc.list[[r]],metrics)
   
} 
list.save(dayref.bc.list,'day_val_bc_hadgem.rdata')

# 
# bm<-mod
# bm$pr.d<-bc.list
# bc.mo<-months(bm)
# ob.mo<-months(obs.test)
# 
# plot(ob.mo$month,ob.mo$pr.mo,col='blue')
# points(bc.mo$month,bc.mo$pr.mo,type='l')
# BC for Mid 21st Century (2021-2060)
daysc.bc.list<-Map(c,lapply(ts.list, '[', 'i'), lapply(ts.list, '[', 'j'),lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat'))
for (r in 1:length(daysc.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]
  
  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)
  
  obs.df<-dplyr::select(cmfd.df,c(2:1))
  
  #day Day Bias Correction for Mid Century Scenario Projections
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)
  
  df<-proj.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2021:2060,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  #Empirical CDF QDM
  
  x.mod<-proj$pr.d
  x.cor<-empqdm(x.mod,obs$pr.d,mod$pr.d)
  bc.df<-proj
  bc.df$pr.d<-x.cor
  
  ref<-list(bc.df,proj)
  names(ref)<-c('day.bc','proj.mod')
  ref<-list(ref)
  ref<-list.names(ref,'day.midsc')
  daysc.bc.list[[r]]<-append(daysc.bc.list[[r]],ref)
  
} 
# BC for Late 21st Century (2061-2100)
for (r in 1:length(daysc.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]
  
  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)
  
  obs.df<-dplyr::select(cmfd.df,c(2:1))
  
  #day Day Bias Correction for Late Century Scenario Projections
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 1979:2005,]
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)
  
  df<-proj.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2061:2100,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  #Empirical CDF QDM
  
  x.mod<-proj$pr.d
  x.cor<-empqdm(x.mod,obs$pr.d,mod$pr.d)
  bc.df<-proj
  bc.df$pr.d<-x.cor
  
  ref<-list(bc.df,proj)
  names(ref)<-c('day.bc','proj.mod')
  ref<-list(ref)
  ref<-list.names(ref,'day.latsc')
  daysc.bc.list[[r]]<-append(daysc.bc.list[[r]],ref)
  
} 
list.save(daysc.bc.list,'day_sc_bc_hadgem.rdata')

  
# STEP 1B:: EXTREME VALUE CORRECTIONS::GUMBEL DQM ------------
#//Manually select RCM from dat.list
load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

#June-July-August Filter FunctionS
jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  return(df.jja)
}
mcs.jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  df.jja<- df.jja[as.numeric(strftime(df.jja$dates, "%Y")) %in% 2021:2060,]
  return(df.jja)
}
lcs.jja<-function(df)
{
  df.jja<- df[as.numeric(strftime(df$dates, "%m")) %in% 6:8,]
  df.jja<- df.jja[as.numeric(strftime(df.jja$dates, "%Y")) %in% 2061:2010,]
  return(df.jja)
}

#Gumbel Function for Return Periods
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

#Bias Correction Wrapper Function
extreme.bias.correct.wrapper<-function(obs.df,mod.df,proj.df)
  #bias correction performed by quantile correction
  #cmfd.df<-observations
  #mod.df<-historical rcm model results
  #val.df<-projection/validation rcm model results
{  
  source('extreme-bc.R')
  source('biasCorrection.R') 
  require(caret)
  
  proj.mom<-extract.ams(proj.df[,2],proj.df[,1])
  mod.mom <- extract.ams(mod.df[,2],mod.df[,1])
  obs.mom<- extract.ams(obs.df[,2],obs.df[,1])
  
  proj.fit<-fit.gumb(proj.mom)
  mod.fit<- fit.gumb(mod.mom)
  obs.fit<-fit.gumb(obs.mom)
  #Variable Preparation for Bias Correction Function
  #**non-essential, but a clear step**
  x.proj<-proj.mom$pr
  par.obs<-c(obs.fit[[1]][[1]])
  par.mod<-c(mod.fit[[1]][[1]])
  par.proj<-c(proj.fit[[1]][[1]])
  #Call Parametric CDF Bias Correction Function
  #===============================================
  bc.mom<-edcdf(x.proj,par.obs,par.mod,par.proj)
  #===============================================
  #Gumbel Fit for Bias Correction
  bc.fit<-fit.gumb(bc.mom)
  par.bc<-c(bc.fit[[1]][[1]])

  require(stats)
  require(goftest)
  #GOF Test for BC-MODEL OUTPUT
  bc.mom<-na.omit(bc.mom)
  ks<-stats::ks.test(bc.mom,pgumbel,bc.fit[[1]][[1]][[1]],bc.fit[[1]][[1]][[2]],alternative='two.sided')
  kstest<-if(ks$p.value>=0.05){"Accept"} else {"Reject"}
  ad<-goftest::ad.test(bc.mom,pgumbel,bc.fit[[1]][[1]][[1]],bc.fit[[1]][[1]][[2]])
  adtest<-if(ad$p.value>=0.05){"Accept"} else {"Reject"}
  cvm<-goftest::cvm.test(bc.mom,pgumbel,bc.fit[[1]][[1]][[1]],bc.fit[[1]][[1]][[2]])
  cvmtest<-if(cvm$p.value>=0.05){"Accept"} else {"Reject"}
  gof.bc<-list(kstest,adtest,cvmtest)
  names(gof.bc)<-c('test.ks','test.ad','test.cvm')
  #GOF TEST FOR PROJ
  ks<-stats::ks.test(proj.mom$pr,pgumbel,proj.fit[[1]][[1]][[1]],proj.fit[[1]][[1]][[2]],alternative='two.sided')
  kstest<-if(ks$p.value>=0.05){"Accept"} else {"Reject"}
  ad<-goftest::ad.test(proj.mom$pr,pgumbel,proj.fit[[1]][[1]][[1]],proj.fit[[1]][[1]][[2]])
  adtest<-if(ad$p.value>=0.05){"Accept"} else {"Reject"}
  cvm<-goftest::cvm.test(proj.mom$pr,pgumbel,proj.fit[[1]][[1]][[1]],proj.fit[[1]][[1]][[2]])
  cvmtest<-if(cvm$p.value>=0.05){"Accept"} else {"Reject"}
  gof.proj<-list(kstest,adtest,cvmtest)
  names(gof.proj)<-c('test.ks','test.ad','test.cvm')
  #GOF TEST FOR OBSERVATION DATASET
  ks<-stats::ks.test(obs.mom$pr,pgumbel,obs.fit[[1]][[1]][[1]],obs.fit[[1]][[1]][[2]],alternative='two.sided')
  kstest<-if(ks$p.value>=0.05){"Accept"} else {"Reject"}
  ad<-goftest::ad.test(obs.mom$pr,pgumbel,obs.fit[[1]][[1]][[1]],obs.fit[[1]][[1]][[2]])
  adtest<-if(ad$p.value>=0.05){"Accept"} else {"Reject"}
  cvm<-goftest::cvm.test(obs.mom$pr,pgumbel,obs.fit[[1]][[1]][[1]],obs.fit[[1]][[1]][[2]])
  cvmtest<-if(cvm$p.value>=0.05){"Accept"} else {"Reject"}
  gof.obs<-list(kstest,adtest,cvmtest)
  names(gof.obs)<-c('test.ks','test.ad','test.cvm')
  test.results<-list(gof.bc,gof.proj,gof.obs)
  names(test.results)<-c('gof.bc','gof.proj','gof.obs')
  
  #Gumbel-Derived Return Levels (Estimated)
  
  params<-data.frame(par.bc,par.obs,par.proj,par.mod)
  r.levels<-data.frame(c(2,5,c=seq(10,200,10)))
  for (i in 1:ncol(params))
  {
   r.levels[,i+1]<-gum.rt(params[1,i],params[2,i],r.levels[,1])
  }
  colnames(r.levels)<-c('years','bc','obs','proj','mod')
  prt.list<-list(r.levels)
  prt.list<-list.names(prt.list,'return.levels')
  #List Script::Gumbel-Derived Return Levels (Estimated)
  
  #List Script:: Monthly Maxima (MoM)
  #//bc.mom >> data.frame
  bc.mom<-data.frame(bc.mom)
  colnames(bc.mom)<-c('pr')
  mom.list<-list(bc.mom,obs.mom,proj.mom,mod.mom)
  names(mom.list)<-c('mom.bc','mom.obs','mom.proj','mom.mod')
  mom.list<-list(mom.list)
  mom.list<-list.names(mom.list,'mom.list')
  
  #List Script:: Gumbel Parameters 
  g.par<-list(par.bc,par.obs,par.proj,par.mod)
  names(g.par)<-c('par.bc','par.obs','par.proj','par.mod')
  g.par<-list(g.par)
  g.par<-list.names(g.par,'gum.par')
  fit.list<-list(bc.fit,obs.fit,proj.fit,mod.fit)
  names(fit.list)<-c('fit.bc','fit.obs','fit.proj','fit.mod')
  fit.list<-list(fit.list)
  fit.list<-list.names(fit.list,'fit.list')
  
  #List Script:: GOF Tests
  #Complicated listing script for test results // Can be simplified
  test.results<-list(test.results)
  test.results<-list.names(test.results,'test.results')
  test.list<-list(ks,ad,cvm)
  names(test.list)<-c('ks','ad','cvm')
  test.list<-list(test.list)
  test.list<-list.names(test.list,'test.list')
  goftests<-append(test.results,test.list)
  
  
  #これが多分大丈夫よ
  
  list<-append(mom.list,g.par)
  list<-append(list,fit.list)
  list<-append(list,goftests)
  list<-append(list,prt.list)
  return(list)
}

# Reference scenario (2013: 2005-2018):
ref.bc.list<-Map(c,lapply(ts.list, '[', 'i'), lapply(ts.list, '[', 'j'),lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat'))
for (r in 1:length(ref.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]

  #Apply JJA Filter
  cmfd.df<-jja(cmfd.df)
  train.df<-jja(train.df)
  proj.df<-jja(proj.df)
  val.df<-jja(val.df)

  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)

  #Train-Val Scheme Inputs::
  require(dplyr)
  # obs.df<-dplyr::select(train.df,c(1:2))
  obs.df<-dplyr::select(cmfd.df,c(2:1))
  #train.df<-dplyr::select(train.df,c(1,3))
  #val.df<-dplyr::select(val.df,c(1,3))

  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  
  for (j in 1:nval)
  {
  metrics.df<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  #Training and Testing Data Split::
  test.fold<-createFolds(obs.df$pr.d, k = 3, list = TRUE, returnTrain = FALSE)
  
  #================K=1================!
  #Extract Test Data::
  obs.test<-obs.df[test.fold[[1]],]  
  mod.test<-mod.df[test.fold[[1]],]  
  #Extract Train Data::
  obs.train<-rbind(obs.df[test.fold[[2]],],obs.df[test.fold[[3]],])
  mod.train<-rbind(mod.df[test.fold[[2]],],mod.df[test.fold[[3]],])
  #Perform Bias Correction::
  rfit.list<-extreme.bias.correct.wrapper(obs.train,mod.train,mod.test)
  #Extract Corrected Return Level data::
  rl.bc<-rfit.list[['return.levels']]
  rl.bc<-dplyr::select(rl.bc,c(1:2))
  #Perform Parametric EV Gumbel Model on Test Data::
  test.mom<- extract.ams(obs.test[,2],obs.test[,1])
  #Derive Test Return Levels::
  test.fit<-fit.gumb(test.mom)
  rl.test<-dplyr::select(rl.bc,c(1))
  rl.test$test<-gum.rt(test.fit[[1]][['estimate']][[1]],test.fit[[1]][['estimate']][[2]],rl.bc[,1])
  #Measure Error::
  comp.df<-cbind(rl.bc,rl.test)
  comp.df<-comp.df[,-3]
  #plot(comp.df$bc,comp.df$test)
  
  
  #MAE
  mae<-mae(comp.df$bc,comp.df$test)
  #RMSE
  rmse<-rmse(comp.df$bc,comp.df$test)
  #Bias
  bias<-bias(comp.df$bc,comp.df$test)
  #Correlation
  cor<-cor(comp.df$bc,comp.df$test)
  temp<-data.frame(mae,rmse,bias,cor)
  metrics.df<-rbind(metrics.df,temp)
  #================K=1================!
  
  #================K=2================!
  #Extract Test Data::
  obs.test<-obs.df[test.fold[[2]],]  
  mod.test<-mod.df[test.fold[[2]],]  
  #Extract Train Data::
  obs.train<-rbind(obs.df[test.fold[[1]],],obs.df[test.fold[[3]],])
  mod.train<-rbind(mod.df[test.fold[[1]],],mod.df[test.fold[[3]],])
  #Perform Bias Correction::
  rfit.list<-extreme.bias.correct.wrapper(obs.train,mod.train,mod.test)
  #Extract Corrected Return Level data::
  rl.bc<-rfit.list[['return.levels']]
  rl.bc<-dplyr::select(rl.bc,c(1:2))
  #Perform Parametric EV Gumbel Model on Test Data::
  test.mom<- extract.ams(obs.test[,2],obs.test[,1])
  #Derive Test Return Levels::
  test.fit<-fit.gumb(test.mom)
  rl.test<-dplyr::select(rl.bc,c(1))
  rl.test$test<-gum.rt(test.fit[[1]][['estimate']][[1]],test.fit[[1]][['estimate']][[2]],rl.bc[,1])
  #Measure Error::
  comp.df<-cbind(rl.bc,rl.test)
  comp.df<-comp.df[,-3]
  
  #MAE
  mae<-mae(comp.df$bc,comp.df$test)
  #RMSE
  rmse<-rmse(comp.df$bc,comp.df$test)
  #Bias
  bias<-bias(comp.df$bc,comp.df$test)
  #Correlation
  cor<-cor(comp.df$bc,comp.df$test)
  temp<-data.frame(mae,rmse,bias,cor)
  metrics.df<-rbind(metrics.df,temp)
  #================K=2================!
  
  
  #================K=3================!
  #Extract Test Data::
  obs.test<-obs.df[test.fold[[3]],]  
  mod.test<-mod.df[test.fold[[3]],]  
  #Extract Train Data::
  obs.train<-rbind(obs.df[test.fold[[1]],],obs.df[test.fold[[2]],])
  mod.train<-rbind(mod.df[test.fold[[1]],],mod.df[test.fold[[2]],])
  #Perform Bias Correction::
  rfit.list<-extreme.bias.correct.wrapper(obs.train,mod.train,mod.test)
  #Extract Corrected Return Level data::
  rl.bc<-rfit.list[['return.levels']]
  rl.bc<-dplyr::select(rl.bc,c(1:2))
  #Perform Parametric EV Gumbel Model on Test Data::
  test.mom<- extract.ams(obs.test[,2],obs.test[,1])
  #Derive Test Return Levels::
  test.fit<-fit.gumb(test.mom)
  rl.test<-dplyr::select(rl.bc,c(1))
  rl.test$test<-gum.rt(test.fit[[1]][['estimate']][[1]],test.fit[[1]][['estimate']][[2]],rl.bc[,1])
  #Measure Error::
  comp.df<-cbind(rl.bc,rl.test)
  comp.df<-comp.df[,-3]
  
  #MAE
  mae<-mae(comp.df$bc,comp.df$test)
  #RMSE
  rmse<-rmse(comp.df$bc,comp.df$test)
  #Bias
  bias<-bias(comp.df$bc,comp.df$test)
  #Correlation
  cor<-cor(comp.df$bc,comp.df$test)
  temp<-data.frame(mae,rmse,bias,cor)
  metrics.df<-rbind(metrics.df,temp)
  #================K=3================!
  
  #================performance summary================!
  metrics.df<-metrics.df[-1,]
  mean.temp<-data.frame(mean(metrics.df$mae),mean(metrics.df$rmse),mean(metrics.df$bias),mean(metrics.df$cor))
  colnames(mean.temp)<-c('mae','rmse','bias','cor')
  mean.met<-rbind(mean.met,mean.temp)
  
  
  }
  #============k-fold cross-validation ===========!
  metrics<-list(mean.met)
  metrics<-list.names(metrics,'crval.3fold')
  rfit.list<-extreme.bias.correct.wrapper(obs.df,mod.df,mod.df)
  rfit.list<-list(rfit.list)
  rfit.list<-list.names(rfit.list,'ref.bc')
  ref.bc.list[[r]]<-append(ref.bc.list[[r]],rfit.list)
  ref.bc.list[[r]]<-append(ref.bc.list[[r]],metrics)

}
list.save(ref.bc.list,'ext_val_bc_hadgem.rdata')
load('ext_val_bc_hadgem.rdata')

# BC for Mid 21st Century (2040: 2021-2060)
sc.bc.list<-Map(c,lapply(ts.list, '[', 'i'), lapply(ts.list, '[', 'j'),lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat'))
for (r in 1:length(sc.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]
  
  #Apply JJA Filter
  cmfd.df<-jja(cmfd.df)
  train.df<-jja(train.df)
  #Apply Mid-2st Century JJA Filter
  df<-proj.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2021:2060,]
  proj.df<-na.omit(df)
  proj.df<-jja(proj.df)
  
  val.df<-jja(val.df)
  
  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)
  
  
  #Fixing observation time series dataframe
  #Dataframe format dates[1] and pr.d [2]
  #Note... cmfd file is reversed
  require(dplyr)
  obs.df<-dplyr::select(cmfd.df,c(2:1))
  rfit.list<-extreme.bias.correct.wrapper(obs.df,mod.df,proj.df)
  rfit.list<-list(rfit.list)
  rfit.list<-list.names(rfit.list,'ext.midsc')
  sc.bc.list[[r]]<-append(sc.bc.list[[r]],rfit.list)
}

r=1
# BC for Late 21st Century (2061-2100)
for (r in 1:length(sc.bc.list))
{
  #Training Data Set::
  train.df<-ds.list[[r]][['train.ts']]
  #Validation Data Set::
  val.df<-ds.list[[r]][['val.ts']]
  #Projection Data Set::
  proj.df<-ts.list[[r]][['fut.ts']]
  #CMFD Data Set::
  cmfd.df<-ts.list[[r]][['cmfd.ts']]
  
  #Apply JJA Filter
  cmfd.df<-jja(cmfd.df)
  train.df<-jja(train.df)
  #Apply Mid-2st Century JJA Filter
  
  df<-proj.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2061:2100,]
  proj.df<-na.omit(df)
  proj.df<-jja(proj.df)
  
  val.df<-jja(val.df)
  
  #Combing the originally subdivided training data set and future model projection
  pr.mod<-c(train.df$his.pr,val.df$fut.pr)
  dates.mod<-c(train.df$dates,val.df$dates)
  mod.df<-data.frame(dates.mod,pr.mod)
  colnames(mod.df)<-c('dates','pr.d')
  mod.df$dates<-as.Date(mod.df$dates)
  
  
  #Fixing observation time series dataframe
  #Dataframe format dates[1] and pr.d [2]
  #Note... cmfd file is reversed
  require(dplyr)
  obs.df<-dplyr::select(cmfd.df,c(2:1))
  rfit.list<-extreme.bias.correct.wrapper(obs.df,mod.df,proj.df)
  rfit.list<-list(rfit.list)
  rfit.list<-list.names(rfit.list,'ext.latsc')
  sc.bc.list[[r]]<-append(sc.bc.list[[r]],rfit.list)
}
list.save(sc.bc.list,'ext_sc_bc_hadgem.rdata')
