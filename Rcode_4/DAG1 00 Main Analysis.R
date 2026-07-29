#DAG1 Analysis

#libraries----------
library(dplyr)
library(tidyr)
library(tibble) #for deframe
library(tidyverse)
library(tidytext) # For reorder_within
library(ggplot2)
library(extrafont)
library(ggcorrplot)
library(ggdendro)
library(reshape2) # For melting the matrix

library(readxl)
library(lubridate)
library(stringr) #for str_remove

library(lavaan)
library(glue)
library(MARSS)

source("LisaXP/functions/sem_utils.R")
source("LisaXP/functions/sem_output_summary_fxn.r")
source("LisaXP/functions/get_top_modindices_fxn.r")
source("LisaXP/functions/find_mds_neighbors_fxn.r")



#main outputs from this script----
    #Helper files
        var_lookup_NCC_AK<-read.csv("LisaXP/outputs_4/var_lookup_NCC_AK.SAR.csv",row.names=1);names(var_lookup_NCC_AK)
        guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1);head(guild.dfas1)
        sem_master_data<-read.csv("LisaXP/outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
        guild.dfas1$sarSR<-sem_master_data$sarSR


    #Result files
        # topmodels<-read.csv("LisaXP/outputs/Doug.topmodels.csv",row.names=1);topmodels
        # importance_summary<-read.csv("LisaXP/outputs/Doug0409SEMd.importance_summary.csv",row.names=1);importance_summary
        # performance_comparison<-read.csv("LisaXP/outputs/aic_DAG1A_v_DAG1B.csv",row.names=1);performance_comparison
        # importance_topvar<-read.csv("LisaXP/outputs/importance_topvar.csv",row.names=1)
        #     importance_topvar<-importance_topvar %>% select(SEMnode, plot_names, percent_occurrence, DAG1A_long, DAG1A_short,  DAG1B_long,  DAG1B_short) %>%
        #           mutate(across(DAG1A_long:DAG1B_short, /(x) round(x, 2)))
        #           write.csv(importance_topvar,"LisaXP/outputs/importance_topvar.daic2.csv")
#==================================
#=====================================
                  path0<- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/"
                  path<- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
                  

#1. Raw Data examine -----
rawData<-readRDS(paste0(path0,"allData.rds")) ;head(rawData)
guilds<-read.csv(paste0(path,"DFA/guildsDFA.csv"),na.strings = "") ;head(guilds) #note: use "guild" to have a known column name for Doug's names: shiftLisa_step3_26jun26 -- i.e., it just pulls the right column from the guilds.csv and renames it "guild"
  table(guilds$guild)
guild.dfas<-read.csv(paste0(path,"DFA/newGuilds.csv"))%>%
  select(2, last_col(offset = 2):last_col()) %>%
  rename(DFAguild=newGuild)
  head(guild.dfas)  #note: use "newGuild" to easily pull only the DFAs that combined multiple taxa for the main guild DFA: excludes smoothed single time series
  table(guild.dfas$DFAguild)
  
#reduce raw data to time series in the final analysis
  dat<-left_join(rawData,guilds,by="shortName");head(dat)
  dat<-left_join(dat,guild.dfas,by="shortName");head(dat)
  dat<-dat %>%
    drop_na(guild)
  head(dat)
  guild.data<-dat
  head(guild.data)
  
  write.csv(guild.data,file="LisaXP/outputs_4/guild.data.csv")

clusData<-readRDS(paste0(path,"DFA/clusDataDFA.rds")) 
head(clusData)



sar.dfa<-clusData %>% filter (SEMlatent=="SAR") %>% select(date,finalVal)
sar.raw<-rawData %>% filter (SEMlatent=="SAR") %>% select(indicator,date,finalVal)
sar.dat<- pivot_wider(sar.raw,names_from=indicator,values_from=finalVal) %>%
          left_join(sar.dfa,by="date") %>%
          rename(DFA=finalVal,sar.SRHW=SAR_2025) %>%
          drop_na(DFA)
head(sar.dat)
#sar$sarDFA<-sar.dfa$finalVal
plot(DFA~date, data=sar.dat,type='l')
    points(scale(sar.SRHW) ~date, data=sar.dat,col=2,pch=16)
    points(scale(sar.UCsprCh_2025) ~date, data=sar.dat,col=3,pch=16)

plot(sar.SRHW~date, data=sar.dat,type='l')

#Set up lookup table

            model="DAG1A_short"
            
            var0<-read.csv(paste0(path,model,"/varImportance.csv")) %>% 
                    mutate(ind1 = indicators %>% 
                             str_remove(",.*") %>% 
                             str_remove("_202[56]")) %>%
                    mutate(prefix = str_replace(var, "^([^\\.]*)\\..*", "\\1"))  %>%           
                    mutate(prefix2=case_when(
                      str_detect(var,"DFA") &  !str_detect(prefix,"DFA") ~ paste0(prefix,"_DFA"),
                      .default =  prefix
                    ),
                    Lisaname=case_when(
                      str_detect(indicators,"Tspin") ~"X02.krill",
                      .default =paste0("X",prefix2,"_",ind1)
                    )) %>%
                    mutate(node_id=case_when(
                      Lisaname=="X02.krill"   ~"X02",
                      .default= paste0("X",prefix)
                    )) %>%
                    mutate(guildname=paste0("X",var)) %>%
                    select(SEMnode,node_id,Lisaname,var,guildname) %>%
                    arrange(Lisaname)
            
            tail(var0)
           
            var_lookup_NCC_AK<-rbind(var0,c("SAR","X16","X16_SAR","16.SAR","X16.SAR"))
    
    #check direction of DFAs-- all good
            
    # #Add column for which variable to reverse  ---------          
    #         var_lookup_NCC_AK<-var_lookup_NCC_AK %>%
    #           mutate(rev=case_when(
    #             str_detect(guildname,"DFA") & 
    #               str_detect(guildname, 
    #                        "ZooPreyNCC_JSOES|Tspin|ForageFishNCC_DFA|PredFishNCC_b|X12.ZooPreyAK_DFA1|X11.PredMammalSmolt|X13.FishPreyAK_DFA1|X13.FishPreyAK_b_DFA1") ~ 1,
    #             .default=0
    #           ))
    #         
            write.csv(var_lookup_NCC_AK,"LisaXP/outputs_4/var_lookup_NCC_AK.SAR.csv")
            
            nrow(var_lookup_NCC_AK)
            head(var_lookup_NCC_AK)
            
      #      var_lookup_NCC_AK<-read.csv("LisaXP/outputs/var_lookup_NCC_AK.SAR.csv",row.names=1);names(var_lookup_NCC_AK)
            

#Enter data file & Reverse sign of some DFAs---------
            guild.dfas0 <- read.csv(paste0(path, model, "/completeness.csv")) %>%
              select(1, any_of(sort(var_lookup_NCC_AK$guildname))) %>%
              mutate(year = year(as.Date(date))) %>% 
              select(year, everything(), -date)
            
            #Rename the columns in guild.dfas0
            translator <- setNames(var_lookup_NCC_AK$Lisaname, var_lookup_NCC_AK$guildname)
            guild.dfas1<-guild.dfas0
            names(guild.dfas1) <- names(guild.dfas1) %>%
              str_trim() %>% # Remove hidden spaces from column names
              { ifelse(. %in% names(translator), translator[.], .) }
            
            # Check results
            names(guild.dfas1)
            
            guild.dfas1<-guild.dfas1 %>%
              mutate(X16_SAR=sar.dat$DFA)

            write.csv(guild.dfas1,"LisaXP/outputs_4/guild.dfa.NCC.AK.csv")
            
       #  # 1. Identify which columns need to be reversed from your lookup table
      # cols_to_reverse <- var_lookup_NCC_AK %>%
      #   filter(rev == 1) %>%
      #   pull(guildname) %>%
      #   intersect(names(guild.dfas0)) # Safety check: only columns that exist
      # 
      # # 2. Create the reversed dataset
      # guild.dfa.rev0 <- guild.dfas0 %>%
      #   mutate(across(all_of(cols_to_reverse), ~ .x * -1))
      # 
      # # 3. Verify a couple of columns flipped
      # # Compare the first few rows of one of the flipped columns
      # head(guild.dfas0[[cols_to_reverse[1]]])
      # head(guild.dfa.rev0[[cols_to_reverse[1]]])
      # 
      # 
      # #4. Rename the columns in guild.dfas.rev0
      # translator <- setNames(var_lookup_NCC_AK$Lisaname, var_lookup_NCC_AK$guildname)
      # guild.dfa.rev1<-guild.dfa.rev0
      # names(guild.dfa.rev1) <- names(guild.dfa.rev1) %>%
      #   str_trim() %>% # Remove hidden spaces from column names
      #   { ifelse(. %in% names(translator), translator[.], .) }
      # 
      # guild.dfa.rev1<-guild.dfa.rev1 %>% 
      #   mutate(X16.SAR=sar$sar.sc)
      # 
      # # Check results
      # names(guild.dfa.rev1)
      
#      write.csv(guild.dfa.rev1,"LisaXP/guild.dfa.NCC.AK.rev1.csv")
      
#      guild.dfa.rev1<-read.csv("LisaXP/guild.dfa.NCC.AK.rev1.csv",row.names=1)
      
#==============
#2. DAG1 Merge aic and importance files from all models ----

# Define the models and folders
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")

# --- Function to read and label files ---
load_sem_data <- function(name, file_type) {
  file_path <- paste0(path,name, "/", ifelse(file_type == "results", "SEMresultsByClus.csv", "varImportance.csv"))
  
  read_csv(file_path) %>%
    mutate(
      model_id = name,
      dag_type = str_extract(name, "DAG1[ABC]"), # Extract DAG1A or DAG1B
      length   = str_extract(name, "long|short") # Extract long or short
    )
}

# --- Load everything into two main tables ----
all_performance <- map_df(model_names, ~load_sem_data(.x, "results"))
all_importance  <- map_df(model_names, ~load_sem_data(.x, "importance"))

importance_summary <- all_importance %>%
      # 1. Reshape and sort the initial data
      select(model_id, SEMnode, var, importance) %>%
      pivot_wider(names_from = model_id, values_from = importance) %>%
      arrange(SEMnode, var) %>%
      
      # 2. Add the custom SAR row, setting the 1s directly inside the new row
      # (Using your 'model_names' vector ensures every model column gets a 1 safely)
      add_row(
        SEMnode = "SAR", 
        var = "16.SAR", 
        !!!setNames(rep(1, length(model_names)), model_names)
      ) %>%
      
      # 3. Join the metadata lookup
      left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = "var") %>%
      
      # 4. Clean, round, and calculate max importances
      mutate(across(DAG1A_long:DAG1C_short, as.numeric)) %>% 
      mutate(across(DAG1A_long:DAG1C_short, ~ round(., 3))) %>%
      rowwise() %>%
      mutate(maxImport = round(max(c_across(DAG1A_long:DAG1C_short), na.rm = TRUE), 2)) %>%
      ungroup() %>%
      
      # 5. Final positioning and sorting
      relocate(SEMnode, Lisaname, var, maxImport, DAG1A_long:DAG1C_short) %>%
      arrange(desc(maxImport))
importance_summary


          # importance_tmp <- all_importance %>%
          #           select(model_id, SEMnode, var, importance) %>%
          #           pivot_wider(names_from = model_id, values_from = importance) %>%
          #           arrange(SEMnode, var)%>%
          #               # Add the new row directly at the end of the pipe
          #           add_row(SEMnode = "SAR", var = "16.SAR") %>%
          #               # Fill all the missing model columns in this new row with 1
          #           mutate(across(everything(), ~ replace_na(., 1)))
          # 
          # 
          # importance_summary <- left_join(importance_tmp, var_lookup_NCC_AK %>% select(var, Lisaname), by = "var") %>%
          #   mutate(across(DAG1A_long:DAG1C_short, as.numeric)) %>% 
          #   mutate(across(DAG1A_long:DAG1C_short, ~ round(., 3))) %>%
          #   rowwise() %>%
          #   mutate(maxImport = round(max(c_across(DAG1A_long:DAG1C_short), na.rm = TRUE),2)) %>%
          #   ungroup() %>%
          #   relocate(SEMnode,Lisaname,var,maxImport,DAG1A_long:DAG1C_short) %>%
          #   arrange(desc(maxImport))

importance_crop<-importance_summary %>% filter(maxImport>0.2)

importance_crop %>% arrange(SEMnode,desc(DAG1A_long))
importance_crop %>% arrange(desc(DAG1C_long))

importance_crop %>% arrange(desc(DAG1C_long))

        print(importance_summary)#,n=Inf)
        
        print(importance_summary %>% filter(maxImport>0.2))
        
        # SEMnode   Lisaname                    var                         maxImport DAG1A_long DAG1A_short DAG1B_long DAG1B_short
        # <chr>     <chr>                       <chr>                           <dbl>      <dbl>       <dbl>      <dbl>       <dbl>
        #   1 Abundance X07_DFA_cpue_IntSprJunHW    07.Cond2NCC_DFA1                 1         1         1            1         1      
        # 2 Growth    X06_DFA_IGF_mu              06.Cond1NCC_DFA1                 1         1.000     0.909        1.000     0.906  
        # 3 PredAK    X10_Harbor_seal_CR_2yrLead  10.PredMammalNCC_0_smoothed      0.99      0.304     0.993        0.211     0.992  
        # 4 PredAK    X10_Harbour_s_2yrLead_WS    10.PredMammalNCC_1_smoothed      0.51      0.513     0.00242      0.495     0.00322
        # 5 PredNCC   X08_commonMurre_JSOES       08.PredBirdNCC_b_0_smoothed      0.48     NA         0.483       NA         0.443  
        # 6 PredNCC   X08_Loons_8_WS              08.PredBirdNCC_b_2_smoothed      0.22      0.222     0.145        0.189     0.104  
        # 7 PredNCC   X09_DFA_HakeAge5Plus        09.PredFishNCC_b_DFA1            0.79      0.742     0.326        0.792     0.436  
        # 8 PreyAK    X12_copepodBiomass_WGoA     12.ZooPreyAK_0_smoothed          0.36      0.262     0.365        0.172     0.292  
        # 9 PreyAK    X12_DFA_biomassEuphShelfSum 12.ZooPreyAK_DFA1                0.46      0.137     0.0119       0.457     0.0993 
        # 10 PreyNCC   X05_DFA_abundSardine        05.ForageFishNCC_DFA1            0.97      0.973     0.837        0.973     0.835  
        # 11 SAR       X16.SAR                     16.SAR                           1         1         1            1         1      

my_performance<-all_performance %>%  select(modNum, model_id,length, dag_type, AIC, CFI, AGFI, p_Chi2,
                                            PreyNCCindNames,       PredNCCindNames, GrowthIndNames, AbundanceIndNames, PreyAKindNames, PredAKindNames) 

topmodel1<-my_performance %>%
        group_by(model_id) %>%
        slice(1) %>% 
        ungroup() %>%
        arrange(length)
topmodel1

print(topmodel1)

        topmodels<-my_performance %>%
          group_by(model_id) %>%
          slice(1:100) %>% 
          ungroup() %>%
          arrange(length)
        
        print(topmodels)

write.csv(topmodels,"LisaXP/outputs_4/Doug.topmodels100.csv")


write.csv(importance_summary,paste("LisaXP/outputs_4/Doug","06_29_SEM_AKPred","importance_summary.csv",sep="_"))


#3. RESULTS -------
#AIC DAG1A vs DAG1B-------------
performance_comparison <- all_performance %>%
  # 1. Ensure we only have ONE row per unique model/length combination
  group_by(model_id) %>%
  slice(1) %>% 
  ungroup() %>%
  
  # 2. Select only what we need
  select(length, dag_type, AIC, CFI, AGFI, p_Chi2) %>%
  
  # 3. Pivot side-by-side
  pivot_wider(names_from = dag_type, 
              values_from = c(AIC, CFI, AGFI, p_Chi2)) %>%
  
  # 4. ONLY convert the metric columns to numeric (avoiding 'length')
  mutate(across(c(contains("AIC"), contains("CFI"), contains("AGFI"), contains("p_Chi2")), 
                ~as.numeric(as.character(unlist(.))))) %>%  
  # 5. Now calculate the delta
  mutate(delta_AIC = abs(AIC_DAG1A - AIC_DAG1B)) %>%
  select(length, contains("AIC"), contains("CFI"), everything())

print(performance_comparison)


# length AIC_DAG1A AIC_DAG1B delta_AIC CFI_DAG1A CFI_DAG1B AGFI_DAG1A AGFI_DAG1B p_Chi2_DAG1A p_Chi2_DAG1B
# 1 long       133.      131.       2.25     0.949     0.969      0.946      0.952       0.160         0.245
# 2 short       64.9      61.6      3.30     0.835     0.925      0.950      0.989       0.0160        0.125
write.csv(performance_comparison,"LisaXP/outputs_4/aic_DAG1A_v_DAG1B_DAG1C.csv")


#Most important variables------------
#note that I did check that a small number of models did substantially better than the rest to justify an aic cutoff (see fig in our project notes doc), ---------
#but I can't find the code so we will have to remake that or discuss whether to use a cutoff this way-------

        # 1. Calculate Delta AIC and Rank
        top_aic_data <- topmodels %>%
          mutate(AIC = as.numeric(as.character(AIC))) %>%
          group_by(model_id) %>%
          # Ensure they are sorted by AIC before calculating rank/delta
          arrange(AIC) %>%
          mutate(
            delta_AIC = AIC - min(AIC, na.rm = TRUE),
            rank = row_number()
          ) %>%
          ungroup()



        #topmodel index importance score--------

        # 1. Filter for your "High Support" models (dAIC < 3)
        top_tier_models <- top_aic_data %>%
          filter(delta_AIC <= 3)
        
        # 2. Reshape and calculate frequencies
        indicator_frequency_table <- top_tier_models %>%
          # Select the metadata and the indicator columns
          select(model_id, length, dag_type, 
                 PreyNCCindNames, PredNCCindNames, GrowthIndNames, AbundanceIndNames) %>%
          # Pivot from wide to long
          pivot_longer(
            cols = c(PreyNCCindNames, PredNCCindNames, GrowthIndNames, AbundanceIndNames),
            names_to = "node",
            values_to = "indNames"
          ) %>%
          # Group by Node and Indicator within each Model Type
          group_by(model_id, node, indNames) %>%
          summarise(count = n(), .groups = "drop_last") %>%
          # Calculate percentage based on the total number of models in that model_id
          mutate(
            total_models_in_tier = sum(count),
            percent_occurrence = round((count / total_models_in_tier) * 100, 1)
          ) %>%
          # Clean up column names and sort
          select(model_id, node, indNames, percent_occurrence, count) %>%
          arrange(model_id, node, desc(percent_occurrence))
        
        indicator_frequency_table %>%
          arrange(desc(percent_occurrence),node)
        # 3. View the table
        print(indicator_frequency_table)
        
        
        # 1. Create the 'long' data once to use for both calculations
        top_tier_long <- top_tier_models %>%
          select(model_id, length, dag_type, 
                 PreyNCCindNames, PredNCCindNames,PreyAKindNames,PredAKindNames, GrowthIndNames, AbundanceIndNames) %>%
          pivot_longer(
            cols = c(PreyNCCindNames, PredNCCindNames, PreyAKindNames,PredAKindNames,GrowthIndNames, AbundanceIndNames),
            names_to = "node",
            values_to = "indNames"
          )
        
        # 2. Calculate the OVERALL frequency (across all 4 model_ids combined)
        overall_freq <- top_tier_long %>%
          group_by(node, indNames) %>%
          summarise(overall_count = n(), .groups = "drop") %>%
          mutate(
            # Denominator is total number of top-tier models found across the study
            overall_percent = round((overall_count / nrow(top_tier_models)) * 100, 1)
          )
        
        # 3. Calculate the PER-MODEL_ID frequency and join with the Overall column
        indicator_frequency_table <- top_tier_long %>%
          group_by(model_id, node, indNames) %>%
          summarise(count = n(), .groups = "drop_last") %>%
          mutate(
            total_models_in_model_id = sum(count),
            percent_occurrence = round((count / total_models_in_model_id) * 100, 1)
          ) %>%
          ungroup() %>%
          # --- JOIN STEP ---
          left_join(overall_freq, by = c("node", "indNames")) %>%
          # Organize columns for readability
          select(node, indNames, overall_percent, percent_occurrence, model_id, count, overall_count) %>%
          # Sort by Node and the most frequent Overall winners
          arrange(node, desc(overall_percent), indNames, model_id)
        
        # View the results
        print(indicator_frequency_table)
        
        
        #unified indicator summary table-------
        # 1. Create the simplified summary
        indicator_summary_once <- top_tier_long %>%
          # Group only by node and the indicator string
          group_by(node, indNames) %>%
          summarise(
            total_count = n(), 
            .groups = "drop"
          ) %>%
          mutate(
            # Percentage is based on the total number of models that made the dAIC < 3 cut
            percent_occurrence = round((total_count / nrow(top_tier_models)) * 100, 1)
          ) %>%
          # 2. Sort to put the most "important" indicators at the top
          arrange(node, desc(percent_occurrence))
        
        # 3. View the clean table
        print(indicator_summary_once)
        
        

        #add Lisanames to summary table--------
        # 1. Clean up the lookup table to avoid the 'node' conflict
        lookup_clean <- var_lookup_NCC_AK %>%
          select(var, node_id, Lisaname) # Keep only necessary columns
        
        # 2. Join with your summary table
        # We match 'var' from the lookup with 'indNames' from your summary
        indicator_summary_lisa <- indicator_summary_once %>%
          left_join(lookup_clean, by = c("indNames" = "var")) %>%
          # 3. Organize the table so biology is front and center
          #select(node, indNames, Lisaname, percent_occurrence, total_count) %>%
          arrange(node_id, desc(percent_occurrence))
        
        indicator_summary_lisa <-indicator_summary_lisa %>%
          mutate(order=case_when(
            node=="AbundanceIndNames" ~ 0,
            node=="GrowthIndNames" ~ 0.5,
            node=="PreyNCCindNames" ~ 1,
            node=="PredNCCindNames" ~ 2,
            node=="PreyAKindNames" ~ 3,
            node=="PredAKindNames" ~ 4
          )) %>%
          arrange(order)
        
        indicator_summary_lisa <-indicator_summary_lisa %>%
          arrange(order,desc(percent_occurrence))
        
        
        # #3. add plotting names
        # indicator_summary_lisa$plot_names<-c("Salmon_Abund_JSOES","Salmon_IGF_JSOES",
        #                                      "Plankton_JSOES","Herring_WCVI",
        #                                      "Murre_JSOES","Hake_Mackerel",
        #   "Prey_WCVI","Copepods_WGOA","Greenlings_AK_EAI","Low-Mid_trophic_DFA_WGOA","FishPreyAK_b_hexagram_EAI","Seabirds_DFA_WGOA",
        #   "Pink_Salmon_All","Pink_Salmon_Asia",
        #   "Harbor_seal_CR_adult","Harbor_seal_WS_adult","Northern_Fur_seal_WS","SSL_wholerange_adult"
        # )
        
        
#FINAL summary of cumulative importance---------        
        # 4. View the result
        head(indicator_summary_lisa)
        
        write.csv(indicator_summary_lisa,"LisaXP/outputs_4/DAG1_summarytable_dAIC3.csv")
        #write.csv(indicator_summary_lisa,"LisaXP/DAG_final/DAG1_summarytable_dAIC3.csv")
        
        #merge w/ actual importance
       importance_topvar<- inner_join(indicator_summary_lisa,importance_summary %>% select(-Lisaname),join_by(indNames==var)) %>%
         filter(order>0)
       importance_topvar<-importance_topvar %>%
       mutate(Trophic=case_when(
         str_detect(SEMnode,"Pred") ~ "Predator",
         str_detect(SEMnode,"Prey") ~ "Prey",
         .default = SEMnode
       )) %>%
         mutate(Region=case_when(
           str_detect(SEMnode,"AK") ~ "BC/AK",
           str_detect(SEMnode,"NCC") ~ "NCC",
           .default = "NCC"
         ))
       head(importance_topvar)
         
       write.csv(importance_topvar,"LisaXP/outputs_4/importance_topvar_daic3.csv")
        
  print(importance_topvar,n=Inf)
            