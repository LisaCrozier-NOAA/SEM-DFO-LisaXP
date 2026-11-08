#Compare current model with previous sstarc model?
#https://connect.fisheries.noaa.gov/CMISST/

dat<-read.csv("Rcode_for_paper/metadata/sem_altprey_data_1998_2021_10212026.csv",na.strings = NA,row.names = 1);head(dat[,1:5])
write.csv(dat %>% select(year,capelin_WGoA_smoltyr),file="CMISST Rshiny/capelin.csv",row.names = FALSE)
write.csv(dat %>% select(year,X09_DFA_HakeAge5Plus_smoltyr),file="CMISST Rshiny/hake.csv",row.names = FALSE)
write.csv(dat %>% select(year,X07_cpue_IntSprJunHW_smoltyr),file="CMISST Rshiny/cpue.csv",row.names = FALSE)
