#===================================================================
# ENSEMBLE MEMBER VALIDATION SCRIPT
#===================================================================!
require(ggplot2)
require(ggpubr)
require(rlist)
source('functions-bc.R')
source('biasCorrection.R')
require('qmap')
require('dplyr')

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
  source('functions-bc.R')
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

#!==========================!  
#RUN VALIDATION STANDARD
nval<-100
#!==========================!
# STEP 0: Data Handling ----------
dat.list <- list()
dat.list = list.files('/Volumes/DAT/2_code/V5-FSLAMR/bc-dat',pattern='*.rdata',full.names=TRUE)
dat.list


#=============HAD-REM=================
g=1
#!====================================!  


load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
list.save(dayref.bc.list,'day_val_bc_remohadgem.rdata')

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
list.save(ref.bc.list,'ext_val_bc_remohadgem.rdata')

#=============MPI-REM=================
g=3
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
list.save(dayref.bc.list,'day_val_bc_remompiesm.rdata')

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
list.save(ref.bc.list,'ext_val_bc_remompiesm.rdata')


#=============HAD-REG==================
g=5
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
  mod$pr.d[mod$pr.d<=0.001]<-0.001
  
  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  j=1
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
    
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE2
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
    mod.test$pr.d[mod.test$pr.d<=0.001]<-0.001
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
    mod.test$pr.d[mod.test$pr.d<=0.001]<-0.001
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  val<-na.omit(df)
  #val<-val[val$pr.d>=1,]
  val<-mjja(val)
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod$pr.d[mod$pr.d<=0.001]<-0.001
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
list.save(dayref.bc.list,'day_val_bc_reg4hadgem.rdata')
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
    test.fold<-createFolds(mod.df$pr.d, k = 3, list = TRUE, returnTrain = FALSE)
    
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
list.save(ref.bc.list,'ext_val_bc_reg4hadgem.rdata')


#=============MPI-REG==================
g=7
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
  #SPECIAL CASE CORRECTION ONLY FOR MPI-ESM-REGCM4 FROM 1980
  obs.df<-obs.df[as.numeric(strftime(obs.df$dates, "%Y")) %in% 1980:2018,]
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)
  mod$pr.d[mod$pr.d<=0.001]<-0.001
  
  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  j=1
  
  
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
    
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE2
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  val<-na.omit(df)
  #val<-val[val$pr.d>=1,]
  val<-mjja(val)
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod$pr.d[mod$pr.d<=0.001]<-0.001
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
list.save(dayref.bc.list,'day_val_bc_reg4mpiesm.rdata')
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
    test.fold<-createFolds(mod.df$pr.d, k = 3, list = TRUE, returnTrain = FALSE)
    
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
list.save(ref.bc.list,'ext_val_bc_reg4mpiesm.rdata')



#!==========================!  
#RUN DAY VALIDATION EXPERIMENT
nval<-1000
#!==========================!
# STEP 0: Data Handling ----------
dat.list <- list()
dat.list = list.files('/Volumes/DAT/2_code/V5-FSLAMR/bc-dat',pattern='*.rdata',full.names=TRUE)
dat.list

#=============HAD-REM=================
g=1
#!====================================!  


load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
list.save(dayref.bc.list,'day_val_bc_remohadgem.rdata')

#=============MPI-REM=================
g=3
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
list.save(dayref.bc.list,'day_val_bc_remompiesm.rdata')



#=============HAD-REG==================
g=5
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
  mod$pr.d[mod$pr.d<=0.001]<-0.001
  
  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  j=1
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
    
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE2
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
    mod.test$pr.d[mod.test$pr.d<=0.001]<-0.001
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
    mod.test$pr.d[mod.test$pr.d<=0.001]<-0.001
    #Perform Bias Correction::
    bc<-empqdm(mod.test$pr.d,obs.train$pr.d,mod.train$pr.d)
    #Measure Error::
    comp.df<-obs.test
    comp.df$bc<-bc
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  val<-na.omit(df)
  #val<-val[val$pr.d>=1,]
  val<-mjja(val)
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod$pr.d[mod$pr.d<=0.001]<-0.001
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
list.save(dayref.bc.list,'day_val_bc_reg4hadgem.rdata')


#=============MPI-REG==================
g=7
#!====================================!  

load(dat.list[g])
ds.list<-x
load(dat.list[g+1])
ts.list<-x

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
  #SPECIAL CASE CORRECTION ONLY FOR MPI-ESM-REGCM4 FROM 1980
  obs.df<-obs.df[as.numeric(strftime(obs.df$dates, "%Y")) %in% 1980:2018,]
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod<-mjja(mod)
  mod$pr.d[mod$pr.d<=0.001]<-0.001
  
  #============k-fold cross-validation ===========!
  require('Metrics')
  require('caret')
  mean.met<-data.frame(mae=0,rmse=0,bias=0,cor=0)
  j=1
  
  
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
    
    #MAE
    mae<-mae(comp.df$bc,comp.df$pr.d)
    #RMSE2
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
    comp.df<-comp.df[is.finite(comp.df$bc),]#!Correction step
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
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  proj<-na.omit(df)
  #proj<-proj[proj$pr.d>=1,]
  proj<-mjja(proj)
  
  
  df<-obs.df
  df<- df[as.numeric(strftime(df$dates, "%Y")) %in% 2006:2018,]
  val<-na.omit(df)
  #val<-val[val$pr.d>=1,]
  val<-mjja(val)
  
  #Day Bias Correction for Reference Scenario Projections
  df<-obs.df
  obs<-na.omit(df)
  #obs<-obs[obs$pr.d>=1,]
  obs<-mjja(obs)
  
  df<-mod.df
  mod<-na.omit(df)
  #mod<-mod[mod$pr.d>=1,]
  mod$pr.d[mod$pr.d<=0.001]<-0.001
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
list.save(dayref.bc.list,'day_val_bc_reg4mpiesm.rdata')




#=============Day Cross-Val Analysis==================
#CV Metrics:: BC for Day Correction
cval.list <- list()
cval.list = list.files('/Volumes/DAT/2_code/data_cv_1000',pattern='*.rdata',full.names=TRUE)
cval.list



cv.extraction<-function(day.ref)
{
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
return(metrics.ref)
}

load(cval.list[1])
hadreg.cv<-cv.extraction(x)
load(cval.list[2])
mpireg.cv<-cv.extraction(x)
load(cval.list[3])
hadrem.cv<-cv.extraction(x)
load(cval.list[4])
mpirem.cv<-cv.extraction(x)

boxplot(hadreg.cv$cv.rmse,hadrem.cv$cv.rmse,mpireg.cv$cv.rmse,mpirem.cv$cv.rmse,main='Daily Rainfall Bias Correction\nCross-Validation RMSE (1979-2018)',ylab='[mm/day]',names=c('HadGEM-RegCM4','MPI-RegCM4','HadGEM-REMO','MPI-REMO'))
boxplot(hadreg.cv$cv.mae,hadrem.cv$cv.mae,mpireg.cv$cv.mae,mpirem.cv$cv.mae,main='Daily Rainfall Bias Correction\nCross-Validation MAE (1979-2018)',ylab='[mm/day]',names=c('HadGEM-RegCM4','MPI-RegCM4','HadGEM-REMO','MPI-REMO'))

boxplot(hadreg.cv$cv.bias,hadrem.cv$cv.bias,mpireg.cv$cv.bias,mpirem.cv$cv.bias)


#=============Ext Cross-Val Analysis==================
#CV Metrics:: BC for Extreme Event Correction
cval.list <- list()
cval.list = list.files('/Volumes/DAT/2_code/data_cv_100',pattern='*.rdata',full.names=TRUE)
cval.list

cv.extraction<-function(ext.ref)
{
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
return(metrics.ref)
}

load(cval.list[5])
hadreg.cv<-cv.extraction(x)
load(cval.list[6])
mpireg.cv<-cv.extraction(x)
load(cval.list[7])
hadrem.cv<-cv.extraction(x)
load(cval.list[8])
mpirem.cv<-cv.extraction(x)

boxplot(hadreg.cv$cv.rmse,hadrem.cv$cv.rmse,mpireg.cv$cv.rmse,mpirem.cv$cv.rmse,main='Extreme Event Bias Correction\nCross-Validation RMSE (1979-2018)',ylab='[mm/day]',names=c('HadGEM-RegCM4','MPI-RegCM4','HadGEM-REMO','MPI-REMO'))
boxplot(hadreg.cv$cv.mae,hadrem.cv$cv.mae,mpireg.cv$cv.mae,mpirem.cv$cv.mae,main='Extreme Event Bias Correction\nCross-Validation MAE (1979-2018)',ylab='[mm/day]',names=c('HadGEM-RegCM4','MPI-RegCM4','HadGEM-REMO','MPI-REMO'))


boxplot(hadreg.cv$cv.bias,hadrem.cv$cv.bias,mpireg.cv$cv.bias,mpirem.cv$cv.bias)