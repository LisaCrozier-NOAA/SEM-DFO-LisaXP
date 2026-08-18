library(dplyr)
library(tidyr)

ak_dat<-read.csv("copilot/outputs_2/data_all_tested_columns_annual.csv");names(ak_dat)
    ak2<-ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom)


sem_master_data<-read.csv("outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
guild.dfas1<-read.csv("outputs_4/guild.dfa.NCC.AK.csv",row.names=1) %>%
  mutate( sarSR=sem_master_data$sarSR,
          sarUC=sem_master_data$sarUC) %>%
  inner_join(trend_df_1998_topbio_allwgoa,by="year") %>%
  inner_join(trend_df_1998_herring_egoa,by="year") %>%
  inner_join(ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom),by="year") 
names(guild.dfas1)

guild.dfas1 %>% rename(
  X13_wgoa_cap.pcod=wgoa_cap.pcod,
  X13_egoa_herring = egoa_herring,
  X13_mid_il_capelin = mid_il_capelin,
  X13_stka_herr_matbiom=stka_herr_matbiom)



tb_mycor <- guild.dfas1 %>%
  summarise(across(everything(), ~ round(cor(.x, X12_DFA_biomassEuphShelfSum, use = "pairwise.complete.obs"), 2))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "cor") %>%
  filter(abs(cor) > 0.5) %>%
  arrange(desc(cor))

(tb_mycor)
# variable                               cor
# 1 X12_DFA_biomassEuphShelfSum           1   
# 2 X13_pollockBiomassGoAage3plus_predAK  0.71
# 3 X04_PacificPompano_                   0.65
# 4 X04_marketsquid_GAM                   0.6 
# 5 X01_DFA_sumPreyOfPrey_planktonJun     0.57
# 6 X01_Red_necked_phal_WS                0.57
# 7 X15_salmonSharkGoA_predAK             0.55
# 8 X06_DFA_IGF_mu                        0.52
# 9 X13_gadid_WAI                        -0.56
# 10 X01_habCompInd                       -0.57
# 11 X05_DFA_abundSardine                 -0.59
# 12 X08_commonMurre_JSOES                -0.63
# 13 X02_DFA_SeaNettle_                   -0.7 

#positively correlated to X12_DFA_biomassEuphShelfSum
pos<- guild.dfas1 %>% select(year,
                             X12_DFA_biomassEuphShelfSum,
                             X15_salmonSharkGoA_predAK,
                             X10_DFA_ssl.est.wholerange_2yrLead,
                             X13_pollockBiomassGoAage3plus_predAK,
                             X12_copepodCom_WGoA,
                             X01_DFA_sumPreyOfPrey_planktonJun, 
                             X09_DFA_HakeAge5Plus)


X01_habCompInd, X02.krill, X05_DFA_abundSardine,
X15_DFA_sleeperSharkBSAI_predAK,
wgoa_cap.pcod,egoa_herring,
X07_DFA_cpue_IntSprJunHW,X16_SAR)

    dat=pos
n=(ncol(dat)-1)
matplot(dat[,1],scale(dat[,-1]),type="b",col=1:n)
legend("topleft",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)
mtext(side=3,"Positively corr w/ Players in DAG1A")
matlines(dat[,1],scale(dat[,"X12_DFA_biomassEuphShelfSum"]),lty=1,col=1,lwd=3)
matlines(dat[,1],scale(dat[,"X15_salmonSharkGoA_predAK"]),lty=1,col=2,lwd=3)
matlines(dat[,1],scale(dat[,"X10_DFA_ssl.est.wholerange_2yrLead"]),lty=1,col=3,lwd=3)

#shark
shark<- guild.dfas1 %>% select(year,
                             X12_DFA_biomassEuphShelfSum,
                             X15_salmonSharkGoA_predAK,
                    #         X10_DFA_ssl.est.wholerange_2yrLead,
                    #         sst_wgoa_coastwatch_junjulaug, 
                             pdo_djf)
ak2<-ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom)


dat=shark
n=(ncol(dat)-1)
matplot(dat[,1],scale(dat[,-1]),type="b",col=1:n)
legend("topleft",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)
mtext(side=3,"Positively corr w/ Players in DAG1A")
matlines(dat[,1],scale(dat[,"X12_DFA_biomassEuphShelfSum"]),lty=1,col=1,lwd=3)
matlines(dat[,1],scale(dat[,"X15_salmonSharkGoA_predAK"]),lty=1,col=2,lwd=3)
matlines(dat[,1],scale(dat[,"sst_wgoa_coastwatch_junjulaug"]),lty=1,col=3,lwd=3)
matlines(dat[,1],scale(dat[,"pdo_djf"]),lty=1,col=4,lwd=3)



#sharks================
        sem.dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
          clean_names()
          names(sem.dat)
        salmon_dat<-sem.dat %>% select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
          names(salmon_dat)
        
        shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
          clean_names() %>%
          select(year,contains("goa"))
          names(shark_dat)
        
        sst_dat<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
          clean_names() %>%
          select(c(1,19:32))
        names(sst_dat)
        
        dat=shark_dat
        n=(ncol(dat)-1)
        matplot(dat[,1],scale(dat[,-1]),type="b",col=1:n)
        legend("topleft",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)
        mtext(side=3,"Positively corr w/ Players in DAG1A")
        matlines(dat[,1],scale(dat[,"X12_DFA_biomassEuphShelfSum"]),lty=1,col=1,lwd=3)
        matlines(dat[,1],scale(dat[,"X15_salmonSharkGoA_predAK"]),lty=1,col=2,lwd=3)
        matlines(dat[,1],scale(dat[,"sst_wgoa_coastwatch_junjulaug"]),lty=1,col=3,lwd=3)
        matlines(dat[,1],scale(dat[,"pdo_djf"]),lty=1,col=4,lwd=3)
        
#csl_cr=================
        jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
        load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
        
        week_raw<-week_final_scaled %>% select(year,week,month,date,Spring_Achin_bonn_pass,csl_nonpup_total_emb,eulachon_ssb_4week_est)
        
        dat<-lre_dat_yearly %>% select(year,csl_nonpup_total_emb)
        n=(ncol(dat)-1)
        matplot(dat[,1],(dat[,-1]),type="b",col=1:n)
        legend("topleft",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)
        mtext(side=3,"Positively corr w/ Players in DAG1A")
        