
#This script 
#1. runs the top model from  Doug's output for 6 DAGs, 
#2. prunes links that are not significant
#3. saves output models for graphing

#Data for modeling
    guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1)
        names(guild.dfas1)

    importance_topvar<-read.csv("LisaXP/outputs_4/importance_topvar_daic3.csv")
        print(importance_topvar %>% select(SEMnode, Lisaname.x, indNames,maxImport) ,n=Inf)
        
        
#Top models ==========
#find names in top models
topmodels<-read.csv("LisaXP/outputs_4/Doug.topmodels100.csv",row.names=1);
        #remove HCI
        head(guilds) #habCompInd 01.ZooPreyNCC_JSOES_1

     x<-   topmodels %>% 
            filter(PreyNCCindNames!="01.ZooPreyNCC_JSOES_1_smoothed") %>%
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
 
 
         
         
         
 
#DAG1C long  =========================
          model="DAG1C_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          y
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"));names(dat)
          
          
          
          
          single_model_text <- glue('
             # Structural Model
               
              X07_DFA_cpue_IntSprJunHW ~ X05_herring_GAM  + X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_Lmu_IntSprJunH +
                        X05_herring_GAM +
                        X09_DFA_HakeAge5Plus +
                        X12_copepodCom_EGoA + 
                        X15_salmonSharkGoA_predAK
          
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          get_top_mi(fit_single)
          
          DAG1C_long_topmodel<-fit_single
 
          
    #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model

              X07_DFA_cpue_IntSprJunHW ~  X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X06_Lmu_IntSprJunH +
                        X05_herring_GAM +
                        X09_DFA_HakeAge5Plus +
                        X12_copepodCom_EGoA + 
                        X15_salmonSharkGoA_predAK
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1C_long_topmodel_reduced<-fit_reduced
#DAG1C short no HCI------------------
          #predAK seals NS, removed
          
          model="DAG1C_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          #Doug top model          
          single_model_text <- glue('
             # Structural Model
              X07_DFA_cpue_IntSprJunHW ~ X04_marketsquid_GAM + X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_StomFull_May +
                        X04_marketsquid_GAM +
                        X08_commonMurre_JSOES +
                        X13_ammod_WAI + 
                        X10_Harbor_seal_CR_2yrLead
         
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1C_short_noHCI_topmodel<-fit_single
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model
              X07_DFA_cpue_IntSprJunHW ~  X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_StomFull_May +
                        X04_marketsquid_GAM +
                        X08_commonMurre_JSOES +
                        X13_ammod_WAI + 
                        X10_Harbor_seal_CR_2yrLead
         
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1C_short_noHCI_topmodel_reduced<-fit_reduced
          
#DAG1C short ------------------
          #predAK seals NS, removed
          
          model="DAG1C_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"))
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
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1C_short_topmodel<-fit_single
          
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
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1C_short_topmodel_reduced<-fit_reduced
          
                    
#DAG1A long  =========================
          #none of the mammal predators on adults were significant, so I removed them
          model="DAG1A_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          y
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"));names(dat)
          
          
          
          
          single_model_text <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X06_DFA_IGF_mu + X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_DFA_biomassEuphShelfSum  
                         + X15_DFA_sleeperSharkBSAI_predAK  
          
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          get_top_mi(fit_single)
          
          DAG1A_long_topmodel<-fit_single
          
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model
            #  X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~  X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_DFA_biomassEuphShelfSum  
                         + X15_DFA_sleeperSharkBSAI_predAK  
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1A_long_topmodel_reduced<-fit_reduced
#DAG1A short ------------------
          #predAK seals NS, removed
          
          model="DAG1A_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          #Doug top model          
          single_model_text <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES + X06_DFA_IGF_mu
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_copepodBiomass_WGoA  
                         + X10_Harbor_seal_CR_2yrLead
          
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1A_short_topmodel<-fit_single
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES #+ X06_DFA_IGF_mu
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X12_copepodBiomass_WGoA  
                         + X10_Harbor_seal_CR_2yrLead
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1A_short_topmodel_reduced<-fit_reduced
          
          
#DAG1B long ------------------
          #predAK seals NS, removed
          
          model="DAG1B_long"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          #Doug top model          
          single_model_text <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_DFA_IGF_mu + 
                        X12_DFA_biomassEuphShelfSum  
                         + X15_PacificCodBiomass_predAK
          
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          # summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1B_long_topmodel<-fit_single
          
          #after NS factors removed
          single_model_reduced <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_DFA_IGF_mu + 
                        X12_DFA_biomassEuphShelfSum  
                         + X15_PacificCodBiomass_predAK
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          #summary(fit_reduced,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_reduced)
          
          DAG1B_long_topmodel_reduced<-fit_reduced
          
          
          
#DAG1B short ------------------
          #predAK seals NS, removed
          
          model="DAG1B_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          
          dat<-guild.dfas1 %>%  select( c(y$Lisaname,"X16_SAR"))
          names(dat)
          
          #Doug top model          
          single_model_text <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        X06_DFA_IGF_mu + 
                        X13_pollockBiomassGoAage3plus_predAK  
                         + X15_DFA_sleeperSharkBSAI_predAK
          
          ')
          
          fit_single <- sem(single_model_text, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          summary(fit_single,standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          get_top_mi(fit_single)
          
          DAG1B_short_topmodel<-fit_single
          
          
          #after NS factors removed
          sem_output_summary_fxn(fit_single)
          single_model_reduced <- glue('
             # Structural Model
              X06_DFA_IGF_mu  ~    X05_DFA_abundSardine
              
              X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
              
              X16_SAR ~ X06_DFA_IGF_mu + 
                        X13_pollockBiomassGoAage3plus_predAK  
                         + X15_DFA_sleeperSharkBSAI_predAK
          
          ')
          
          fit_reduced <- sem(single_model_reduced, data = guild.dfas1,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_reduced)
          get_top_mi(fit_reduced)
          
          DAG1B_short_topmodel_reduced<-fit_reduced
          
          
        
#RESULTS================
          sem_output_summary_fxn(DAG1C_short_topmodel)
          sem_output_summary_fxn(DAG1C_long_topmodel)

          sem_output_summary_fxn(DAG1C_short_topmodel_reduced)
          sem_output_summary_fxn(DAG1C_long_topmodel_reduced)
          
          # Extract just the number of parameters
          num_params <- fitMeasures(final_models_list[[1]], "npar")
          
          print(num_params)

   save(
     DAG1A_short_topmodel,
     DAG1A_long_topmodel,
     DAG1B_long_topmodel,
     DAG1B_short_topmodel,
     DAG1C_long_topmodel,
     DAG1C_short_topmodel,
     DAG1C_short_noHCI_topmodel,
     
     DAG1A_short_topmodel_reduced,
     DAG1A_long_topmodel_reduced,
     DAG1B_long_topmodel_reduced,
     DAG1B_short_topmodel_reduced,
     DAG1C_long_topmodel_reduced,
     DAG1C_short_topmodel_reduced,
     DAG1C_short_noHCI_topmodel_reduced,
     
   file="LisaXP/outputs_4/DAG1.topmodels_noHCI.rdata",
   envir = .GlobalEnv
     )   
   
   
   final_models<-list(
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
   myaic<-AIC(
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
     
   
   myaic %>% arrange(df,)
   
   fitmeasures(DAG1A_short_topmodel)
   library(tidyverse)
   library(lavaan)
   
   # 1. Give the list elements clean names so the table rows are labeled correctly
   names(final_models) <- c(
     "DAG1A_short", "DAG1A_long", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short",
     "DAG1A_short_reduced", "DAG1A_long_reduced", "DAG1B_long_reduced", 
     "DAG1B_short_reduced", "DAG1C_long_reduced", "DAG1C_short_reduced"
   )
   
   # 2. Extract the metrics from each model fit object
   model_summary_table <- map_dfr(names(final_models), function(mod_name) {
     fit <- final_models[[mod_name]]
     
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
   