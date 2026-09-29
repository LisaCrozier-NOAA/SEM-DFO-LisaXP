#interpret alt prey code

print(all_sweep_df,n=Inf)
#add emb csl
path_smooth<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/datWide_1998_2021_indiv_indicators_noNA_MARSSsmoothed.csv"
path_raw<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/datWide_1995_2025_reprocessed_raw.csv"
path_dfa<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/datWide_1996_2025_MARSSextended_noNA_smoltyr_adultyr_labeled.csv"

emb_marss<-dat<-read.csv(path_smooth,row.names = NULL) %>%
  select(Year,AllSeaLionsEMB_2025,AllSeaLionsBonn_2025,ssl_est_wholerange) %>%
  rename(AllSeaLionsEMB_2025_marss=AllSeaLionsEMB_2025,
         AllSeaLionsBonn_2025_marss=AllSeaLionsBonn_2025,
         ssl_est_wholerange_marss=ssl_est_wholerange);head(dat)
emb_raw<-dat<-read.csv(path_raw,row.names = NULL) %>%
  select(Year,AllSeaLionsEMB_2025,AllSeaLionsBonn_2025)%>%
  rename(AllSeaLionsEMB_2025_raw=AllSeaLionsEMB_2025,
       AllSeaLionsBonn_2025_raw=AllSeaLionsBonn_2025);head(dat)
       
emb_dfa<-dat<-read.csv(path_dfa,row.names = NULL) %>%
  select(Year,X10.PredMammalNCC_DFA1_adultyr);head(dat)

emb_all<-left_join(emb_marss,emb_raw,join_by("Year"))
emb_all<-left_join(emb_all,emb_dfa,join_by("Year"))
n=ncol(emb_all)-1
matplot(emb_all$Year,scale(emb_all[,-1]),lty=1:n,col=1:n,type='b')
legend("topleft",legend=paste(1:n,names(emb_all[,-1]),sep="_"),bty='n',lty=1:n,col=1:n)

n=2
matplot(emb_all$Year,scale(emb_all[,c(2,5)]),lty=1:n,col=1:n,type='b')
legend("topleft",legend=paste(1:n,names(emb_all[,c(2,5)]),sep="_"),bty='n',lty=1:n,col=1:n)

matplot(emb_all$Year,scale(emb_all[,c(3,6)]),lty=1:n,col=1:n,type='b')
legend("topleft",legend=paste(1:n,names(emb_all[,c(3,6)]),sep="_"),bty='n',lty=1:n,col=1:n)

emb_new<-emb_all %>% select(Year,AllSeaLionsEMB_2025_raw,AllSeaLionsBonn_2025_raw) %>%
  rename(X10_AllSeaLionsEMB_2025_raw=AllSeaLionsEMB_2025_raw,
         X10_AllSeaLionsBonn_2025_raw=AllSeaLionsBonn_2025_raw)