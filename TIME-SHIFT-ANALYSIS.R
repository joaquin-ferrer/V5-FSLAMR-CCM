#3>>Record Date Uncertainty Incorporation----------
load('final_cmfd_inv_list.rdata')
cmfd.inv <-x

#3.Setting Uncertainty Parameters (Time Shift Window)---------
#Set Time Window for Uncertainty Incorporation
days.buffer.mn=7 # Days prior.
days.buffer.pl=0 # Days after was not used, but can be an option



#3.Max Event Search Function---------

#Delay function optimized for maximum event search
delay.function <- function(cmfd.inv,days.buffer.mn,days.buffer.pl)
{
  for(i in 1:length(cmfd.inv))
  {
    events.v2 <- cmfd.inv[[i]][['pr.events']]
    daysdiff.v2 <- select(events.v2, c(1:2,4:11))
    events.v2 <- select(events.v2, c(1:2))
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



days.buffer.mn=0 # Days prior.
mn.0<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.0<-rbindlist(lapply(mn.0, `[[`, 'v2.events'))
mn.0<-mn.0[as.numeric(strftime(mn.0$dates, "%m")) %in% 6:8,]


days.buffer.mn=1 # Days prior.
mn.1<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.1<-rbindlist(lapply(mn.1, `[[`, 'v2.events'))

mn.1<-mn.1[as.numeric(strftime(mn.1$dates, "%m")) %in% 6:8,]

days.buffer.mn=2 # Days prior.
mn.2<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.2<-rbindlist(lapply(mn.2, `[[`, 'v2.events'))

mn.2<-mn.2[as.numeric(strftime(mn.2$dates, "%m")) %in% 6:8,]

days.buffer.mn=3 # Days prior.
mn.3<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.3<-rbindlist(lapply(mn.3, `[[`, 'v2.events'))
mn.3<-mn.3[as.numeric(strftime(mn.3$dates, "%m")) %in% 6:8,]

days.buffer.mn=4 # Days prior.
mn.4<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.4<-rbindlist(lapply(mn.4, `[[`, 'v2.events'))

mn.4<-mn.4[as.numeric(strftime(mn.4$dates, "%m")) %in% 6:8,]


days.buffer.mn=5 # Days prior.
mn.5<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.5<-rbindlist(lapply(mn.5, `[[`, 'v2.events'))
mn.5<-mn.5[as.numeric(strftime(mn.5$dates, "%m")) %in% 6:8,]


days.buffer.mn=7 # Days prior.
mn.7<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.7<-rbindlist(lapply(mn.7, `[[`, 'v2.events'))
mn.7<-mn.7[as.numeric(strftime(mn.7$dates, "%m")) %in% 6:8,]


days.buffer.mn=11 # Days prior.
mn.11<-delay.function(cmfd.inv,days.buffer.mn,days.buffer.pl)
mn.11<-rbindlist(lapply(mn.11, `[[`, 'v2.events'))
mn.11<-mn.11[as.numeric(strftime(mn.11$dates, "%m")) %in% 6:8,]


a1<-ggplot() +
  geom_point(data=mn.0, aes(x=pr.d,y=ant.30/30,colour="Record Date",shape="Record Date"),size=4)+
  geom_point(data=mn.1, aes(x=pr.d,y=ant.30/30,colour="-1 Day"),size=2)+
  # geom_point(data=mn.2, aes(x=pr.d,y=ant.30/30,colour="-2 Days"),size=2)+    
  # geom_point(data=mn.5, aes(x=pr.d,y=ant.30/30,colour="-5 Days"),size=2)+
  ggtitle("CMFD Rainfall of Inventory Events \n 30-Day Antecedent Rainfall vs Event Rainfall \nJune-August (1995-2005) ") +
  labs(x="Pe [mm/day]",y="Pa [mm/day]",color="Legend",shape="Legend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)


a2<-ggplot() +
  geom_point(data=mn.0, aes(x=pr.d,y=ant.30/30,colour="Record Date",shape="Record Date"),size=3)+
  geom_point(data=mn.2, aes(x=pr.d,y=ant.30/30,colour="-2 Days",shape="-2 Days"),size=2.5)+ 
  geom_point(data=mn.7, aes(x=pr.d,y=ant.30/30,colour="-7 Days",shape="-7 Days"),size=2)+ 
  ggtitle("CMFD Rainfall of Inventory Events \nEvent Rainfall vs 30-Day Antecedent Rainfall\nJune-August (1995-2005) ") +
  labs(x="Pe [mm/day]",y="Pa [mm/day]",color="Legend",shape="Legend")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
a2

  geom_point(data=mn.1, aes(x=pr.d,y=ant.30/30,colour="-1 Day",shape="-1 Day"),size=3)+
  geom_point(data=mn.2, aes(x=pr.d,y=ant.30/30,colour="-2 Days",shape="-2 Days"),size=2)+    
  geom_point(data=mn.5, aes(x=pr.d,y=ant.30/30,colour="-5 Days",shape="-5 Days"),size=1)+
  geom_point(data=mn.7, aes(x=pr.d,y=ant.30/30,colour="-7 Days",shape="-7 Days"),size=4)+ 
  geom_point(data=mn.11, aes(x=pr.d,y=ant.30/30,colour="-11 Days",shape="-11 Days"),size=4)
  
  
  ggtitle("CMFD Rainfall of Inventory Events \nAverage 20-Day Antecedent Rainfall vs Event Rainfall \nJune-August (1995-2005) ") +
  labs(x="Pe [mm/day]",y="Pa [mm/day]",color="Legend",shape="Legend") +
  theme(aspect.ratio=1)



#mn2--------------
  
  #3.Setting Uncertainty Parameters (Time Shift Window)---------
  #Set Time Window for Uncertainty Incorporation
  days.buffer.mn=0 # Days prior.
  days.buffer.pl=0 # Days after was not used, but can be an option
  
  
  
  #3.Max Event Search Function---------
  
  #Delay function optimized for maximum event search
  delay.function <- function(cmfd.inv,days.buffer.mn,days.buffer.pl)
  {
    for(i in 1:length(cmfd.inv))
    {
      events.v2 <- cmfd.inv[[i]][['pr.events']]
      daysdiff.v2 <- select(events.v2, c(1:2,4:11))
      events.v2 <- select(events.v2, c(1:2))
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
  list.save(inv, 'mn0_days_cmfd.rdata')
  
  #4>>Application of EasyBal Recharge-Rainfall Fraction Factors-----------
  
  #4.Creating Recharge Time Series --------------
  load('mn0_days_cmfd.rdata')
  mn.inv <-x
  
  mn.ant30<- (lapply(mn.inv, `[[`, 'ant.ts'))
  mn<- rbindlist(lapply(mn.inv, `[[`, 'v2.events'))
  mn <- select(mn,c(1:3,11)) #30-day antecedent rainfall analysis 
  
  
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
  #list.save(mn.inv,'mn2_pa_inv.rdata')
  #End Create Pa.ts
  library(dplyr)
  i=1
  for(i in 1:length(mn.inv))
  {
    a <- refrac
    names(a)[1]  <-'dates'
    b <- data.frame(select(mn.inv[[i]][['v2.events']],c(1:3)),mn.inv[[i]][['v2.events']][['ant.30']])
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
  list.save(mn.inv,'mn0_days_cmfd.rdata')
  # #Start - Create Pa.events
  
  
  
  
  load('mn11_days_cmfd.rdata')
  mn11.r <-x
  
  
  load('mn0_days_cmfd.rdata')
  mn0.r <-x
  
  load('mn2_days_cmfd.rdata')
  mn2.r <-x
  
  load('mn7_days_cmfd.rdata')
  mn7.r <-x
  
  mn7.r<- rbindlist(lapply(mn7.r, `[[`, 'pa.events'))
  mn2.r<- rbindlist(lapply(mn2.r, `[[`, 'pa.events'))
  
  
  mn7.r<-mn7.r[as.numeric(strftime(mn7.r$dates, "%m")) %in% 6:8,]
  mn2.r<-mn2.r[as.numeric(strftime(mn2.r$dates, "%m")) %in% 6:8,]
  
  
  mn11.r<- rbindlist(lapply(mn11.r, `[[`, 'pa.events'))
  mn11.r<-mn11.r[as.numeric(strftime(mn11.r$dates, "%m")) %in% 6:8,]

  
  mn0.r<- rbindlist(lapply(mn0.r, `[[`, 'pa.events'))
  mn0.r<-mn0.r[as.numeric(strftime(mn0.r$dates, "%m")) %in% 6:8,]
  
  
  a3<-ggplot() +
    geom_point(data=mn0.r, aes(x=pr.d,y=qa.30,colour="Record Date",shape="Record Date"),size=3)+
    geom_point(data=mn2.r, aes(x=pr.d,y=qa.30,colour="-2 Days",shape="-2 Days"),size=2.5)+ 
    geom_point(data=mn7.r, aes(x=pr.d,y=qa.30,colour="-7 Days",shape="-7 Days"),size=2)+ 
    ggtitle("Event Rainfall vs Recharge Estimates\nJune-August (1995-2005) ") +
    labs(x="Pe [mm/day]",y="Pa [mm/day]",color="Legend",shape="Legend")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
  a3
  
  
  mn0.0.qa<-mn0.r[mn0.r$qa.30==0,]
  mn2.0.qa<-mn2.r[mn2.r$qa.30==0,]
  mn7.0.qa<-mn7.r[mn7.r$qa.30==0,]
  mn11.0.qa<-mn11.r[mn11.r$qa.30==0,]
  
  
  a4<-ggplot() +
    geom_point(data=mn2.0.qa, aes(x=pr.d,y=qa.30,colour="Estimated Recharge",shape="Estimated Recharge"),size=2.5)+
    geom_point(data=mn0.0.qa, aes(x=pr.d,y=ant.30,colour="-0 Days",shape="-0 Days"),size=3)+
    #geom_point(data=mn2.0.qa, aes(x=pr.d,y=ant.30,colour="-2 Days",shape="-2 Days"),size=2.5)+
    geom_point(data=mn7.0.qa, aes(x=pr.d,y=ant.30,colour="-7 Days",shape="-7 Days"),size=2)+ 
    #geom_point(data=mn11.0.qa, aes(x=pr.d,y=ant.30,colour="-11 Days",shape="-11 Days"),size=2)+ 
    ggtitle("Event Rainfall vs Recharge Estimates\nJune-August (1995-2005) ") +
    labs(x="Pe [mm/day]",y="Pa [mm/day]",color="Legend",shape="Legend")+
    geom_text(data=mn7.0.qa,aes(x=pr.d,y=ant.30,label=dates),nudge_x = 5, nudge_y = 0.01,size=2.5,check_overlap = T) +
    #geom_text(data=mn11.0.qa,aes(x=pr.d,y=ant.30,label=dates),nudge_x = 5, nudge_y = 0.01,size=2.5,check_overlap = T) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),aspect.ratio=1)
  a4
  
  
  