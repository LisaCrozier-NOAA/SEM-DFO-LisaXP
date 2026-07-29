

trend_goa_bio_1998<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/mybio_DFAtrends.csv")



#This script 
#1. runs the top model from  Doug's output for 6 DAGs, 
#2. prunes links that are not significant
#3. saves output models for graphing

#Data for modeling
    guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1)
        names(guild.dfas1)

    importance_topvar<-read.csv("LisaXP/outputs_4/importance_topvar_daic3.csv")
        print(importance_topvar %>% select(SEMnode, Lisaname.x, indNames,maxImport) ,n=Inf)
        
    topmodels<-read.csv("LisaXP/outputs_4/Doug.topmodels100.csv",row.names=1);
        
        
    #add columns for SR and UC separately------------
        sar.dfa<-clusData %>% filter (SEMlatent=="SAR") %>% select(date,finalVal)
        sar.raw<-rawData %>% filter (SEMlatent=="SAR") %>% select(indicator,date,finalVal)
        sar.dat<- pivot_wider(sar.raw,names_from=indicator,values_from=finalVal) %>%
          left_join(sar.dfa,by="date") %>%
          rename(DFA=finalVal,sar.SRHW=SAR_2025) %>%
          drop_na(DFA)
        head(sar.dat)

        plot(DFA~date, data=sar.dat,type='l')
            points(scale(sar.SRHW) ~date, data=sar.dat,col=2,pch=16)
            points(scale(sar.UCsprCh_2025) ~date, data=sar.dat,col=3,pch=16)

        guild.dfas1$sarSR<-sar.dat$sar.SRHW
        guild.dfas1$sarUC<-sar.dat$sar.UCsprCh_2025
        
#DAG1C long  =========================
          model="DAG1C_long"  

          
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

          DAG1C_long_topmodel_reduced<-fit_reduced
##DAG1C short ------------------
          #predAK seals NS, removed
          
          model="DAG1C_short"  

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

          DAG1C_short_topmodel_reduced<-fit_reduced
          
#DAG1A long  =========================
          #none of the mammal predators on adults were significant, so I removed them
          model="DAG1A_long"  

          
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

          DAG1A_long_topmodel_reduced<-fit_reduced
#DAG1A short ------------------
          #predAK seals NS, removed
          
          model="DAG1A_short"  

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

          DAG1A_short_topmodel_reduced<-fit_reduced
          
          
#DAG1B long ------------------
          #predAK seals NS, removed
          
          model="DAG1B_long"  

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

          DAG1B_long_topmodel_reduced<-fit_reduced
          
          
          
#DAG1B short ------------------
          #predAK seals NS, removed
          
          model="DAG1B_short"  
          getnames<-x %>% filter(model_id==model) %>% select(grep("Names",names(x)))
          y<-var_lookup_NCC_AK %>% filter(var %in% unlist(getnames))
          y

          
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

          DAG1B_short_topmodel_reduced<-fit_reduced
          
          
        
#RESULTS================
          sem_output_summary_fxn(DAG1C_short_topmodel)
          sem_output_summary_fxn(DAG1C_long_topmodel)

          sem_output_summary_fxn(DAG1C_short_topmodel_reduced)
          sem_output_summary_fxn(DAG1C_long_topmodel_reduced)
          

   save(
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
     DAG1C_short_topmodel_reduced,
     
   file="LisaXP/outputs_4/DAG1.topmodels.rdata",
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
   
   model_summary_table %>%
     arrange(desc(N_Years),npar,AIC)

   
   
   
