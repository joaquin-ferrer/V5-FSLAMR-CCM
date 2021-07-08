#Extracting RG RCM TS

require(rlist)
require(data.)
load('rgs.rdata')
ts.list<-x


rg<-data.table::rbindlist(Map(c,lapply(ts.list, '[', 'i'), lapply(ts.list, '[', 'j'),lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))
rg.coor<-rg[1,3:4]


cval.list <- list()
cval.list = list.files('/Volumes/DAT/2_code/data_cv_100',pattern='*.rdata',full.names=TRUE)
cval.list


load(cval.list[5])
ts.list<-x
rcm<-data.table::rbindlist(Map(c,lapply(ts.list, '[', 'lon'),lapply(ts.list, '[', 'lat')))

plot(rcm)
points(rg.coor)
#Use 10!


cval.list <- list()
cval.list = list.files('/Volumes/DAT/2_code/data_sc_day',pattern='*.rdata',full.names=TRUE)
cval.list

load(cval.list[1])
rcm<-x[[10]][['day.midsc']][['day.bc']]

write.csv(rcm,'hadgem_regcm4_mid.csv')
rcm<-x[[10]][['day.latsc']][['day.bc']]
write.csv(rcm,'hadgem_regcm4_lat.csv')

load(cval.list[2])

rcm<-x[[10]][['day.midsc']][['day.bc']]
write.csv(rcm,'mpiesm_regcm4_mid.csv')
rcm<-x[[10]][['day.latsc']][['day.bc']]
write.csv(rcm,'mpiesm_regcm4_lat.csv')

load(cval.list[3])

rcm<-x[[10]][['day.midsc']][['day.bc']]
write.csv(rcm,'hadgem_remo_mid.csv')
rcm<-x[[10]][['day.latsc']][['day.bc']]
write.csv(rcm,'hadgem_remo_lat.csv')

load(cval.list[4])

rcm<-x[[10]][['day.midsc']][['day.bc']]
write.csv(rcm,'mpiesm_remo_mid.csv')
rcm<-x[[10]][['day.latsc']][['day.bc']]
write.csv(rcm,'mpiesm_remo_lat.csv')
