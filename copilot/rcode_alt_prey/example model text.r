
dag1a_reduced_txt <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw 
'


sem_txt <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + issl_z + i_issl_sst + i_issl_herr + i_issl_cap
'

sem_txt_altprey <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw +  ssl_seak_pup_pred*x13_mid_il_capelin + ssl_seak_pup_pred*x13_stka_herr_matbiom
'

#sem_complete_data<- read.csv("copilot/outputs_9/sem_altprey_data_complete_1998_2021.csv", row.names = FALSE)
guildnames<-sort(tolower(names(sem_complete_data)[-1]));guildnames

guild.dfasAK <- read.csv("data_Lisa/guild.dfasAK.csv",row.names=1) %>% clean_names()
names(guild.dfasAK)
summary(guild.dfasAK)


lookup_pred_dfa_altprey <- read.csv("copilot/outputs_8/all_pred_dfa_altprey.csv",row.names=NULL) %>% clean_names()
names(lookup_pred_dfa_altprey)

lookup_pred_dfa_altprey <- read.csv("copilot/outputs_8/all_pred_dfa_altprey.csv",row.names=NULL) %>% clean_names()
names(lookup_pred_dfa_altprey)


summary(guild.dfasAK)summary(guild.dfasAK)summary(guild.dfasAK)

#DAG1a long top model reduced
DAG1a_long_txt <- '
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw +x12_dfa_biomasseuphshelfsum + x15_dfa_sleepersharks
'

#DAG1a long top model reduced
DAG1b_long_txt <- '
  x06_dfa_igf_mu ~ x05_dfa_abundsardine
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hake_age5plus
  x16_sar ~ x06_dfa_igf_mu + x07_dfa_cpue_intsprjunhw +x12_dfa_biomasseuphshelfsum + x15_pacificcodbiomass_predak
'

bespoke_altprey_txt <- '
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw +  x11_ssl_seak_pup_pred*x13_mid_il_capelin + x11_ssl_seak_pup_pred*x13_stka_herr_matbiom
'

#more recent???-------
source("functions/sem_output_summary_fxn.r")

sem_data_tolower<-sem_complete_data<- read.csv("copilot/outputs_9/sem_altprey_data_complete_1998_2021.csv", row.names = NULL)
names(sem_data_tolower)<-tolower(names(sem_complete_data))
sort(names(sem_data_tolower)[-1])

#DAG1a long top model reduced
DAG1a_long_txt <- '
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw +x12_dfa_biomasseuphshelfsum + x15_dfa_sleepersharks
'

#DAG1a long top model reduced
DAG1b_long_txt <- '
  x06_dfa_igf_mu ~ x05_dfa_abundsardine
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x06_dfa_igf_mu + x07_dfa_cpue_intsprjunhw +x12_dfa_biomasseuphshelfsum + x15_pacificcodbiomass_predak
'

bespoke_altprey_txt <- '
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw +  x11_ssl_seak_pup_pred*x13_mid_il_capelin + x11_ssl_seak_pup_pred*x13_stka_herr_matbiom
'

sem_text<-bespoke_altprey_txt
bespoke_altprey_fit<-fit

sem_text<-DAG1a_long_txt
sem_text<-DAG1b_long_txt

fit <-  sem(sem_text, data = sem_data_tolower, std.lv = TRUE, missing = "ML", warn = FALSE);sem_output_summary_fxn(fit)


DAG1a_long_fit<-fit
DAG1b_long_fit<-fit


sem_output_summary_fxn(fit)
sem_output_summary_fxn(DAG1a_long_fit)
sem_output_summary_fxn(DAG1b_long_fit)
sem_output_summary_fxn(bespoke_altprey_fit)



print(sweep_results,n=Inf)
sweep_results<-sweep_results %>% arrange(aic)
print(sweep_results %>% filter(model_type=="Full",pvalue>0.05),n=Inf)
print(sweep_results %>% filter(model_type=="Pruned",pvalue>0.05),n=Inf)
