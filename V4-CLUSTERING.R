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

#clustering analysis
library(cluster)
library(factoextra)


#1>>Global Setup--------
load('mn7_days_cmfd.rdata')
mn.inv <-x
list.events <- lapply(mn.inv, `[[`, 'pa.events')
mn<- rbindlist(list.events )
jja.mn<- mn[as.numeric(strftime(mn$dates, "%m")) %in% 6:8,]


list.qats <- lapply(mn.inv, `[[`, 'pa.ts')
qa.30<- rbindlist(list.qats)


load('final_cmfd_inv_list.rdata')
cmfd.inv <-x
ts<-rbindlist(lapply(cmfd.inv, `[[`, 'ts'))
ant.ts<-rbindlist(lapply(cmfd.inv, `[[`, 'ant.ts'))
data <- data.frame(ts,ant.ts$ant.30,qa.30$qa.30)

data<- data[as.numeric(strftime(data$dates, "%m")) %in% 6:8,]
data <- data[,-2]
colnames(data)<-c('pr','ant.30','qa.30')


#2.3

#minus 7 days - centroid -rainfall
centroids <-  data.frame(jja.mn$pr.d,jja.mn$ant.30)
colnames(centroids)<-c('pr','ant.30')


#Plotting Rainfall vs Pe
plot(data$pr,data$ant.30/30,col='darkgrey',main="7-Day Adjusted CMFD Rainfall Combinations from June-August (1995-2005)\nComparison of Maximum Antecedent Event Selection Criteria",xlab="Pe [mm/day]", ylab="30-day Average Rainfall [mm/day]")
points(centroids$pr,centroids$ant.30,col='blue')
legend(95,3,legend=c('7-Day Maximum Event Criteria','CMFD Rainfall Combinations'),col=c("blue","darkgrey"), pch=1, cex=0.8)




#minus 7 days - centroid - recharge
centroids <-  data.frame(jja.mn$pr.d,jja.mn$ant.30)
colnames(centroids)<-c('pr','ant.30')


#Plotting recharge vs Pe
plot(data$pr,data$qa.30,col='darkgrey',main="7-Day Adjusted CMFD Rainfall Combinations from June-August (1995-2005)\nComparison of Maximum Antecedent Event Selection Criteria",xlab="Pe [mm/day]", ylab="30-day Average Rainfall [mm/day]")
points(centroids$pr,centroids$qa.30,col='blue')
legend(95,5,legend=c('7-Day Maximum Event Criteria','CMFD Rainfall Combinations'),col=c("blue","darkgrey"), pch=1, cex=0.8)


# 
# data <-scale(data)
# centroids <- scale(centroids)
set.seed(123)
kmeans_fit <- kmeans(centroids,centers=8,nstart = 25)
p3 <- fviz_cluster(kmeans_fit, geom = "point",  data = centroids, alpha = 0.4) + ggtitle("K-means Algorithm, k=8")               
p3
pam_fit <- pam(centroids,k=8)
p2 <- fviz_cluster(pam_fit, geom = "point",  data = centroids, alpha = 0.4) + ggtitle("k =5")               
p2





# 
# #Plotly Sucks
library(plotly)
# 
df<-centroids[,1:2]
df<-scale(df)
centroids$cluster<-factor(pam(df,8)$cluster)
 
# p3 <- plot_ly(centroids, x=~pr, y=~ant.30, color=~cluster) %>%add_markers(size=5)
# p3 <- p3 %>% add_markers(data=data, x =~pr, y = ~ant.30/30,color='blue')
# print(p3)
# 
#
   

#PAM Clustering and Plotting 
    df<-centroids[,1:2]
    df<-scale(df)
    centroids$cluster<-factor(pam(df,5)$cluster)
   
    library(plyr)
    df<-centroids
    find_hull <- function(df) df[chull(df$pr, df$ant.30), ]
    hulls <- ddply(df, 'cluster', find_hull)
   
    library(ggpubr)
    ggplot()+
    geom_point(data=data,aes(x=pr,y=ant.30/30),colour='grey',alpha=0.5)+
    
    geom_point(data=centroids,aes(x=pr,y=ant.30,color=cluster),size=1.5,alpha=1.2)+
    geom_polygon(data = hulls,aes(x=pr,y=ant.30,color=cluster), alpha = 0.5) +
    theme_pubr()



# centroids$cluster<-factor(pam(df,8)$cluster)
 # p3 <- plot_ly(data=data, x =~pr, y = ~ant.30/30) %>% add_markers(data=centroids, x=~pr, y=~ant.30, color=~cluster) 
 # print(p3)

