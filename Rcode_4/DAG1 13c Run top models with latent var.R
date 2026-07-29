
#CONCLUSIONS
#only 2 models still had latent variables when all NS removed
#or they didn't have any to start with
sem_output_summary_fxn(DAG1A_long_latent_reduced_final)
sem_output_summary_fxn(DAG1C_long_latent_reduced_final)


#This script 
#1. runs the top model from  Doug's output for 6 DAGs, 
#2. prunes links that are not significant
#3. saves output models for graphing

#Data for modeling
    sem_master_data<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1)
        names(sem_master_data)

    importance_topvar<-read.csv("LisaXP/outputs_4/importance_topvar_daic3.csv")
        print(importance_topvar %>% select(SEMnode, Lisaname.x, indNames,maxImport) ,n=Inf)
        
    sem_master_data<-read.csv("LisaXP/outputs_4/sem_master_data.csv",row.names = NULL)
        
#DAG1 top models
#find names in top models
topmodels<-read.csv("LisaXP/outputs_4/Doug.topmodels100.csv",row.names=1);

     x<-   topmodels %>%
            group_by(model_id) %>%
            slice(1) %>% 
            ungroup() 

           x %>% select(model_id,PreyAKindNames,PredAKindNames)
           x %>% select(model_id,PreyNCCindNames,PredNCCindNames)
           
           x4<-x %>% select(model_id,PreyNCCindNames,PredNCCindNames,PreyAKindNames,PredAKindNames);x4
           #    model_id    PreyNCCindNames                PredNCCindNames             PreyAKindNames           PredAKindNames
           # 1 DAG1A_long  05.ForageFishNCC_DFA1          09.PredFishNCC_b_DFA1       12.ZooPreyAK_DFA1          15.PredFishAK_b_DFA1       
           # 2 DAG1A_short 05.ForageFishNCC_DFA1          08.PredBirdNCC_b_0_smoothed 12.ZooPreyAK_0_smoothed    10.PredMammalNCC_1_smoothed
           # 3 DAG1B_long  05.ForageFishNCC_DFA1          09.PredFishNCC_b_DFA1       12.ZooPreyAK_DFA1          15.PredFishAK_3_smoothed   
           # 4 DAG1B_short 05.ForageFishNCC_DFA1          08.PredBirdNCC_b_0_smoothed 13.FishPreyAK_b_4_smoothed 15.PredFishAK_b_DFA1       
           # 5 DAG1C_long  05.ForageFishNCC_1_smoothed    09.PredFishNCC_b_DFA1       12.ZooPreyAK_1_smoothed    15.PredFishAK_b_2_smoothed 
           # 6 DAG1C_short 01.ZooPreyNCC_JSOES_1_smoothed 08.PredBirdNCC_b_0_smoothed 14.CompAK_0_smoothed       10.PredMammalNCC_1_smoothed
           
           # 1. Grab just the "Names" columns from your 4-model dataset and flatten them
           getnames <- x4 %>% 
             select(grep("Names", names(x4))) %>% 
             unlist()
           
           # 2. Filter the lookup table using the raw, unlisted names
           y <- var_lookup_NCC_AK %>% 
             filter(var %in% getnames)
           y
           
           #    SEMnode node_id                    Lisaname                         var                    guildname
           # 1 PreyNCC     X05        X05_DFA_abundSardine       05.ForageFishNCC_DFA1       X05.ForageFishNCC_DFA1
           # 2 PredNCC     X09        X09_DFA_HakeAge5Plus       09.PredFishNCC_b_DFA1       X09.PredFishNCC_b_DFA1
           # 3  PredAK     X10  X10_Harbor_seal_CR_2yrLead 10.PredMammalNCC_0_smoothed X10.PredMammalNCC_0_smoothed
           # 4  PredAK     X10    X10_Harbour_s_2yrLead_WS 10.PredMammalNCC_1_smoothed X10.PredMammalNCC_1_smoothed
           # 5  PreyAK     X12 X12_DFA_biomassEuphShelfSum           12.ZooPreyAK_DFA1           X12.ZooPreyAK_DFA1
           # 6  PreyAK     X12     X12_copepodBiomass_WGoA     12.ZooPreyAK_0_smoothed     X12.ZooPreyAK_0_smoothed 
 
 
#latent_syntax_blocks-------         
           cat(paste(latent_syntax_blocks, collapse = "\n\n"))        
           
           Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
           
           Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
           
           Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
           
           Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
           
           Latent_13_FishPreyAK_b =~ gadid_EAI + ammod_EAI + capelin_WGoA + estAbundHerringRecruits + WGOA_DFA_seabirds_2026
           
           Latent_13_FishPreyAK =~ rhinoAuk_EGoA + sitkaHerring_EGoA + WGOA_DFA_lowerTrophic_2026 + WGOA_DFA_midTrophic_2026 + pollockBiomassAIage1plus_predAK_2026
           
           Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026         
 
#DAG1C long  =========================
          model="DAG1C_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          y
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"));names(dat)
          
          

          
          single_model_text <- glue('

            # MEASUREMENT MODEL

               Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
                Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
                Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~ X05_herring_GAM  + Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        Latent_06_Cond1NCC +
                        X05_herring_GAM +
                        Latent_09_PredFishNCC_b +
                        X12_copepodCom_EGoA + 
                        Latent_15_PredFishAK_b
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          get_top_mi(fit_single)
          summary(fit_single)
          
          DAG1C_lofit_singleDAG1C_long_latent<-fit_single
 
          
    #after NS factors removed from structural path
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
 
            # MEASUREMENT MODEL

                Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
                Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X05_herring_GAM +
                         X12_copepodCom_EGoA + 
                        Latent_15_PredFishAK_b
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          summary(fit_reduced)
          
          DAG1C_long_latent_reduced<-fit_reduced
          
          #after NS factors removed from latent path
          summary(fit_single)
          single_model_reduced_latent <- glue('
 
            # MEASUREMENT MODEL

                Latent_09_PredFishNCC_b =~  HakeAge5Plus_2025 + darkblotchedRockfish_2025  +  jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 
                Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026

             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X05_herring_GAM +
                         X12_copepodCom_EGoA + 
                        Latent_15_PredFishAK_b
          
          ')
          
          fit_reduced_latent <- sem(single_model_reduced_latent, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced_latent)
          summary(fit_reduced_latent)
          
          DAG1C_long_latent_reduced_final<-fit_reduced_latent
          
          
          
##NO LATENT DAG1C short ------------------
          #predAK seals NS, removed
          
          model="DAG1C_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          #Doug top model          
          single_model_text <- glue('
             # Structural Model
              X07_DFA_cpue_IntSprJunHW ~ X01_habCompInd + X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_StomFull_May +
                        X01_habCompInd +
                        X08_commonMurre_JSOES +
                        X14_pinkSalmon + 
                        X10_Harbor_seal_CR_2yrLead
         
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1C_short_latent<-fit_single
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model
              X07_DFA_cpue_IntSprJunHW ~  X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_StomFull_May +
                        X01_habCompInd +
                        X08_commonMurre_JSOES +
                        X14_pinkSalmon + 
                        X10_Harbor_seal_CR_2yrLead
         
          ')
          
          fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1C_short_latent_reduced<-fit_reduced
          
#DAG1A long  =========================
          #none of the mammal predators on adults were significant, so I removed them
          model="DAG1A_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          y
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"));names(dat)
          
          # Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
          # 
          # Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
          # 
          # Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
          # 
          # Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
          # 
          # Latent_13_FishPreyAK_b =~ gadid_EAI + ammod_EAI + capelin_WGoA + estAbundHerringRecruits + WGOA_DFA_seabirds_2026
          # 
          # Latent_13_FishPreyAK =~ rhinoAuk_EGoA + sitkaHerring_EGoA + WGOA_DFA_lowerTrophic_2026 + WGOA_DFA_midTrophic_2026 + pollockBiomassAIage1plus_predAK_2026
          # 
          # Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026         
          
          
          single_model_text <- glue('
             #Measurement model
            
                Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
                Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
                Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
                Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
                Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026         
 
            
             # Structural Model
              Latent_06_Cond1NCC  ~    Latent_05_ForageFishNCC
              
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        Latent_12_ZooPreyAK  
                         + Latent_15_PredFishAK_b  
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single)
          
          #ERRORS-----------
          # Warning messages:
          #   1: lavaan->lav_mvnorm_missing_h1_estimate_moments():  
          #   Maximum number of iterations reached when computing the sample moments using EM; use the em.h1.iter.max= 
          #   argument to increase the number of iterations 
          # 2: lavaan->lav_mvnorm_missing_h1_estimate_moments():  
          #   The smallest eigenvalue of the EM estimated variance-covariance matrix (Sigma) is smaller than 1e-05; this 
          # may cause numerical instabilities; interpret the results with caution. 
          # 3: lavaan->lav_object_post_check():  
          #   some estimated ov variances are negative 
          # 
          DAG1A_long_latent<-fit_single
          
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             #Measurement model
            
                Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
                Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
                Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026         
 
            
             # Structural Model
              
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        Latent_12_ZooPreyAK  
                         + Latent_15_PredFishAK_b  
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          summary(fit_reduced)
          
          DAG1A_long_latent_reduced<-fit_reduced
          
          
          #after NS factors removed from latent path
          summary(fit_single)
          single_model_reduced_latent <- glue('
 
              #Measurement model
            
                Latent_09_PredFishNCC_b =~  HakeAge5Plus_2025 + darkblotchedRockfish_2025   + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025
                Latent_12_ZooPreyAK =~  leastAuk_WAI  + biomassAmphShelfSum + biomassEuphShelfSum 
         
            
             # Structural Model
              
              X07_DFA_cpue_IntSprJunHW ~  Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        Latent_12_ZooPreyAK  
                         + sleeperSharkGoA_predAK_2026  
          
          ')
          
          fit_reduced_latent <- sem(single_model_reduced_latent, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced_latent)
          summary(fit_reduced_latent)
          
          DAG1A_long_latent_reduced_final<-fit_reduced_latent
          #I removed the shark DFA, which solved the negative variance problem but it greatly reduced r2 on SAR
          
#NO LATENT LEFT DAG1A short ------------------
          #predAK seals NS, removed
          
          model="DAG1A_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          #Doug top model          
          single_model_text <- glue('
          
          #Measurement Model
                Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
                Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025

             # Structural Model
              Latent_06_Cond1NCC  ~    Latent_05_ForageFishNCC
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES + Latent_06_Cond1NCC
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_copepodBiomass_WGoA  
                         + X10_Harbor_seal_CR_2yrLead
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1A_short_latent<-fit_single
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
         #Measurement Model
 
             # Structural Model

              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_copepodBiomass_WGoA  
                         + X10_Harbor_seal_CR_2yrLead
           
          ')
          
          fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          summary(fit_reduced)
          
          DAG1A_short_latent_reduced<-fit_reduced
          #no latent variables left when removed cpue ~ growth
          
#NO LATENT LEFT DAG1B long ------------------
          #Measurement Model
          Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
          
          Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
          
          Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
          
          Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
          
          Latent_13_FishPreyAK_b =~ gadid_EAI + ammod_EAI + capelin_WGoA + estAbundHerringRecruits + WGOA_DFA_seabirds_2026
          
          Latent_13_FishPreyAK =~ rhinoAuk_EGoA + sitkaHerring_EGoA + WGOA_DFA_lowerTrophic_2026 + WGOA_DFA_midTrophic_2026 + pollockBiomassAIage1plus_predAK_2026
          
          Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026         
          
          
          model="DAG1B_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          #Doug top model          
          single_model_text <- glue('
          
              # Measurement Model
                Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
                
                Latent_09_PredFishNCC_b =~  HakeAge5Plus_2025 + darkblotchedRockfish_2025 +   jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 
                
                Latent_12_ZooPreyAK =~  leastAuk_WAI +  biomassAmphShelfSum + biomassEuphShelfSum 
                
  
             # Structural Model
              IGF_mu_2025  ~    Latent_05_ForageFishNCC
              
              X07_DFA_cpue_IntSprJunHW ~ Latent_09_PredFishNCC_b 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        IGF_mu_2025 + 
                        Latent_12_ZooPreyAK 
                        
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single)
          get_top_mi(fit_single)
          
          DAG1B_long_latent<-fit_single
          
          #I had to remove NS var in order to get it to converge without errors
          # #after NS factors removed
          # single_model_reduced <- glue('
          #    # Structural Model
          #     X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
          #     
          #     X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus 
          #     
          #     X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
          #               X06_DFA_IGF_mu + 
          #               X12_DFA_biomassEuphShelfSum  
          #                + X15_PacificCodBiomass_predAK
          # 
          # ')
          # 
          # fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          # sem_output_summary_fxn(fit_reduced)
          # #summary(fit_reduced,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          # get_top_mi(fit_reduced)
          # 
          # DAG1B_long_latent_reduced<-fit_reduced
          
          
          
#NO LATENT LEFT DAG1B short ------------------
          #predAK seals NS, removed
          
          model="DAG1B_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-sem_master_data %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          
          #Doug top model          
          single_model_text <- glue('
          
          # Measurement Model

              Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025

          # Structural Model
              IGF_mu_2025  ~    Latent_05_ForageFishNCC
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        IGF_mu_2025 + 
                        X13_pollockBiomassGoAage3plus_predAK  
                         + sleeperSharkGoA_predAK_2026
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1B_short_latent<-fit_single
          
          
          #after NS factors removed -- nothing left!
          # sem_output_summary_fxn(fit_single)
          # single_model_reduced <- glue('
          #    # Structural Model
          # 
          #     
          #  
          # ')
          # 
          # fit_reduced <- sem(single_model_reduced, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          # sem_output_summary_fxn(fit_reduced)
          # summary(fit_reduced)
          # 
          # DAG1B_short_latent_reduced<-fit_reduced
          
          
        
#RESULTS================
          sem_output_summary_fxn(DAG1C_long_latent_reduced_final)
          sem_output_summary_fxn(DAG1A_long_latent_reduced_final)

          
   save(
     DAG1A_long_latent_reduced_final,
     DAG1C_long_latent_reduced_final,

   file="LisaXP/outputs_4/DAG1.latents.rdata",
   envir = .GlobalEnv
     )   
   
   
   final_models_latent<-list(
     
     DAG1A_short_latent,
     DAG1A_long_latent,
     DAG1B_long_latent,
     DAG1B_short_latent,
     DAG1C_long_latent,
     DAG1C_short_latent,
     
     DAG1A_short_latent_reduced,
     DAG1A_long_latent_reduced,
     DAG1B_long_latent_reduced,
     DAG1B_short_latent_reduced,
     DAG1C_long_latent_reduced,
     DAG1C_short_latent_reduced,
     
     
     DAG1A_short_topmodel,
     DAG1A_long_topmodel,
     DAG1B_long_topmodel,
     DAG1B_short_topmodel,
     DAG1C_long_topmodel,
     DAG1C_short_topmodel,
     
     DAG1A_short_topmodel_reduced,
     DAG1A_long_topmodel_reduced,
     DAG1B_long_topmodel_reduced,
     DAG1B_short_topmodel_reduced,
     DAG1C_long_topmodel_reduced,
     DAG1C_short_topmodel_reduced)
     
#compare aic---------

   # 1. Give the list elements clean names so the table rows are labeled correctly
   names(final_models_latent) <- c(
     "DAG1A_short", "DAG1A_long", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short",
     "DAG1A_short_reduced", "DAG1A_long_reduced", "DAG1B_long_reduced", 
     "DAG1B_short_reduced", "DAG1C_long_reduced", "DAG1C_short_reduced"
   )
   
   # 2. Extract the metrics from each model fit object
   model_summary_table <- map_dfr(names(final_models_latent), function(mod_name) {
     fit <- final_models_latent[[mod_name]]
     
     # Extract standard fit measures safely
     measures <- fitMeasures(fit)
     
     tibble(
       Model = mod_name,
       N_Years = nobs(fit),                    # Number of observation years
       npar    = measures[["npar"]],            # Number of parameters
       df      = measures[["df"]],              # Degrees of freedom
       AIC     = round(measures[["aic"]], 2),   # Akaike Information Criterion
       AGFI    = round(measures[["agfi"]], 3)   # Adjusted Goodness of Fit Index
     )
   })
   
   # 3. View the completed table
   print(model_summary_table)   
   
   
   #    Model               N_Years  npar    df   AIC  AGFI
   # 1 DAG1A_short              19    12     9 104.  0.946
   # 2 DAG1A_long               24    12     9 156.  0.898
   # 3 DAG1B_long               24    12     9 147.  0.985
   # 4 DAG1B_short              19    12     9 104.  0.973
   # 5 DAG1C_long               24    12     3  74.8 0.981
   # 6 DAG1C_short              19    12     3  18.4 0.997
   # 7 DAG1A_short_reduced      19    12     9 105.  0.947
   # 8 DAG1A_long_reduced       24     8     3  96.2 0.485
   # 9 DAG1B_long_reduced       24    12     9 147.  0.985
   # 10 DAG1B_short_reduced     19    12     9 103.  0.985
   # 11 DAG1C_long_reduced      24    11     4  72.8 0.985
   # 12 DAG1C_short_reduced     19    11     4  17.1 0.998
   
   
#unacceptably bad fit:
   # 8 DAG1A_long_reduced       24     8     3  96.2 0.485
   
   
   model_summary_table %>%
     arrange(desc(N_Years),AIC,npar)
   
   #    Model               N_Years  npar    df   AIC  AGFI
   # 1 DAG1C_long_reduced       24    11     4  72.8 0.985
   # 2 DAG1C_long               24    12     3  74.8 0.981
   # 3 DAG1A_long_reduced       24     8     3  96.2 0.485
   # 4 DAG1B_long               24    12     9 147.  0.985
   # 5 DAG1B_long_reduced       24    12     9 147.  0.985
   # 6 DAG1A_long               24    12     9 156.  0.898
   # 7 DAG1C_short_reduced      19    11     4  17.1 0.998
   # 8 DAG1C_short              19    12     3  18.4 0.997
   # 9 DAG1B_short_reduced      19    12     9 103.  0.985
   # 10 DAG1A_short             19    12     9 104.  0.946
   # 11 DAG1B_short             19    12     9 104.  0.973
   # 12 DAG1A_short_reduced     19    12     9 105.  0.947   
   
   sem_output_summary_fxn(DAG1C_long_topmodel_reduced)
   sem_output_summary_fxn(DAG1C_long_topmodel)
   
   sem_output_summary_fxn(DAG1C_short_topmodel_reduced)
   sem_output_summary_fxn(DAG1C_long_topmodel_reduced)
  
   model_summary_table %>%
     arrange(desc(N_Years),npar,AIC)

   
   summary(DAG1C_long_topmodel) #df=3
   summary(DAG1B_long_topmodel) #df=9
   