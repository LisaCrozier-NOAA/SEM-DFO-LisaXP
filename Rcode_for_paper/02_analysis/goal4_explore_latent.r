
library(dplyr)
library(lavaan)
library(glue)
source("functions/sem_utils.R")
source("functions/sem_output_summary_fxn.r")
source("functions/get_top_modindices_fxn.r")



#Search more directly for the best parameters

#Data for modeling
sem_data_indiv<-read.csv("Rcode_for_paper/metadata/datWide_1998_2021_indiv_indicators_noNA_MARSSsmoothed.csv",row.names=NULL);names(sem_data_indiv)
sem_data_dfas<-read.csv("Rcode_for_paper/metadata/clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv",row.names=NULL) %>%
  filter(Year>=1998,Year<=2021);names(sem_data_dfas)
dat<-left_join(sem_data_dfas,sem_data_indiv,by="Year")
names(dat)

"C:\Users\Lisa.Crozier\Documents\Marine survival\SEM-DFO-LisaXP\Rcode_for_paper\metadata\clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv"
#get list of shortnames in hake DFA
crosswalk<-read.csv("Rcode_for_paper/metadata/master_name_crosswalk.csv",row.names=NULL);names(crosswalk)
crosswalk %>% filter(grepl("Hake|hake",LisaName))
crosswalk %>% filter(LisaName=="X09_DFA_HakeAge5Plus") %>% pull(rankedIndicators)

 # Latent_X09_DFA_HakeAge5Plus =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
 # 
 # cpue_IntSprJunHW_2025 ~ Latent_X09_DFA_HakeAge5Plus
 # 
 # grep("cpue",names(sem_data_indiv),value=T)
 # grep("SAR",names(sem_data_indiv),value=T)
 # cpue_IntSprJunHW_2025
 # SAR_2025
 
 
 
 
 
 single_model_text <- glue('

            # MEASUREMENT MODEL

                Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X15_DFA_sleeperSharks_adultyr
          
          ')
 
 fit_single <- sem(single_model_text, data = dat,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
 sem_output_summary_fxn(fit_single)
 get_top_mi(fit_single)
 summary(fit_single)
 
 
 # --- Model Fit Metrics ---
 #   pvalue     cfi     aic   rmsea    agfi 
 #    0.000   0.501 408.203   0.269   0.258 
 # 
 # --- R-Squared (Endogenous Variables) ---
 #   X07_DFA_cpue_IntSprJunHW                  X16_SAR 
 #                    0.405                    0.457 
 # 
 # --- Regression Path Ranking (By P-Value) ---
 #   lhs op                           rhs    est    se      z pvalue
 # X16_SAR  ~ X15_DFA_sleeperSharks_adultyr -0.153 0.066 -2.328  0.020
 # X07_DFA_cpue_IntSprJunHW  ~       Latent_09_PredFishNCC_b -0.470 0.140 -3.349  0.001
 # X16_SAR  ~      X07_DFA_cpue_IntSprJunHW  0.716 0.188  3.804  0.000
 
 dat$sharkxtemp_adult<-dat$X15_DFA_sleeperSharks_adultyr*dat$X21_sst_egoa_junjulaug_adultyr
 dat$sharkxtemp_smolt<-dat$X15_DFA_sleeperSharks_smoltyr*dat$X21_sst_egoa_junjulaug_smoltyr
 dat$hakextemp_smolt<-dat$X09_DFA_HakeAge5Plus*dat$X21_sst_egoa_junjulaug_smoltyr
 
 #after NS factors removed from latent path
 single_model_text <- glue('

            # MEASUREMENT MODEL

                Latent_09_PredFishNCC_b =~  HakeAge5Plus_2025 + darkblotchedRockfish_2025   + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
              X15_DFA_sleeperSharks_adultyr
         #                 sharkxtemp_adult
         #               X15_DFA_sleeperSharks_adultyr + X21_sst_egoa_junjulaug_adultyr + sharkxtemp_adult
         #               X15_DFA_sleeperSharks_smoltyr + X21_sst_egoa_junjulaug_smoltyr + sharkxtemp_smolt

          ')
 
 fit_single <- sem(single_model_text, data = dat,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
 summary(fit_single)
 sem_output_summary_fxn(fit_single)
 
 # w/ sharks (no temp)
 # pvalue     cfi     aic   rmsea    agfi 
 # 0.000   0.650 264.143   0.276   0.339 
 # --- R-Squared (Endogenous Variables) ---
 #   X07_DFA_cpue_IntSprJunHW                  X16_SAR 
 # 0.384                    0.459 
 # # 
 # sharkxtemp (w/o shark indiv)
 # pvalue     cfi     aic   rmsea    agfi 
 # 0.000   0.645 263.894   0.280   0.338 
 # X07_DFA_cpue_IntSprJunHW                  X16_SAR 
 # 0.384                    0.466 
 # 
 
 
 #BASE MODEL -- no AK
 single_model_text <- glue('

            # MEASUREMENT MODEL

                Latent_09_PredFishNCC_b =~  HakeAge5Plus_2025 + darkblotchedRockfish_2025   + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW 

          ')
 
 fit_single <- sem(single_model_text, data = dat,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
 summary(fit_single)
 sem_output_summary_fxn(fit_single)
 
 
 # pvalue     cfi     aic   rmsea    agfi 
 # 0.040   0.866 281.727   0.176   0.512 
 # 
 # --- R-Squared (Endogenous Variables) ---
 #   X07_DFA_cpue_IntSprJunHW                  X16_SAR 
 # 0.387                    0.300 
 # 
 # --- Regression Path Ranking (By P-Value) ---
 #   lhs op                      rhs    est    se      z pvalue
 # X07_DFA_cpue_IntSprJunHW  ~  Latent_09_PredFishNCC_b -0.459 0.141 -3.257  0.001
 # X16_SAR  ~ X07_DFA_cpue_IntSprJunHW  0.665 0.207  3.209  0.001
 
 
 #definitely better w/ sharks:
 # pvalue     cfi     aic   rmsea    agfi 
 # 0.001   0.731 278.840   0.233   0.404 
 # 
 # --- R-Squared (Endogenous Variables) ---
 #   X07_DFA_cpue_IntSprJunHW                  X16_SAR 
 # 0.387                    0.457 
 # 
 # --- Regression Path Ranking (By P-Value) ---
 #   lhs op                           rhs    est    se      z pvalue
 # X16_SAR  ~ X15_DFA_sleeperSharks_adultyr -0.153 0.066 -2.328  0.020
 # X07_DFA_cpue_IntSprJunHW  ~       Latent_09_PredFishNCC_b -0.459 0.141 -3.257  0.001
 # X16_SAR  ~      X07_DFA_cpue_IntSprJunHW  0.716 0.188  3.804  0.000