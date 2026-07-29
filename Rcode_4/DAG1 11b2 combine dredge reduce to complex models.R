#Apples to apples comparison between dredge and DAG1s

library(tidyverse)

      #trusted indicators
      ecological_synthesis<-read.csv("LisaXP/outputs_4/Ecological_Synthesis_Summary_Table.csv")
      
      #combine 2 dredge products
      dredge_results_noX06X07<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names = NULL) %>%
        mutate(adjaic2=Raw_AIC + 2*Missing_Yrs) %>%
        filter(!grepl("X01_habCompInd",Predictors)) %>%
        arrange(adjaic2)%>%
        mutate(
          Growth  = if_else(grepl("Growth", Guild_Sources), 1, 0),
          Abundance=if_else(grepl("Abundance",  Guild_Sources), 1, 0),
          PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
          PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
          PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
          PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
        ) %>%
        rowwise() %>%
        mutate(Nguilds=sum(c_across(Growth:PredAK))) %>%
        ungroup() # <-- CRITICAL FIX: Strip rowwise grouping here
      
      dredge_results_X06X07<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results_with_X06X07.csv",row.names = NULL)%>%
        mutate(adjaic2=Raw_AIC + 2*Missing_Yrs) %>%
        filter(!grepl("X01_habCompInd",Predictors)) %>%
        arrange(adjaic2) %>%
        mutate(
          Growth  = if_else(grepl("Growth", Guild_Sources), 1, 0),
          Abundance=if_else(grepl("Abundance",  Guild_Sources), 1, 0),
          PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
          PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
          PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
          PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
        ) %>%
        rowwise() %>%
        mutate(Nguilds=sum(c_across(PreyNCC:PredAK))) %>%
        ungroup() # <-- CRITICAL FIX: Strip rowwise grouping here
      
      dredge_results<-rbind(dredge_results_X06X07,dredge_results_noX06X07)
      
      
      dredge_results <-dredge_results %>%
        mutate(
          Growth  = if_else(grepl("Growth", Guild_Sources), 1, 0),
          Abundance=if_else(grepl("Abundance",  Guild_Sources), 1, 0),
          PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
          PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
          PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
          PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
        ) %>%
        rowwise() %>%
        mutate(Nguilds=sum(c_across(Growth:PredAK))) %>%
        ungroup() # <-- CRITICAL FIX: Strip rowwise grouping here
      
       test0<- dredge_results %>% 
        head() %>% 
        relocate(Growth,Abundance, .before=PreyNCC) %>%
        relocate(Nguilds,.after = PredAK) %>%
         rowwise() %>%
         mutate(Nguilds=sum(c_across(Growth:PredAK))) %>%
         ungroup() # <-- CRITICAL FIX: Strip rowwise grouping here
       test0
     
       dredge_results<- dredge_results %>% 
         relocate(Growth,Abundance, .before=PreyNCC) %>%
         relocate(Nguilds,.after = PredAK) %>%
         rowwise() %>%
         mutate(Nguilds=sum(c_across(Growth:PredAK))) %>%
         ungroup() # <-- CRITICAL FIX: Strip rowwise grouping here
       
       dredge_results %>% 
         head()
       
      head(dredge_results)
      
      complex_models_only <- dredge_results %>%
        filter(Num_Params > 2) %>%
        filter(Nguilds > 2) %>%
        arrange(adjaic2)

      nrow(complex_models_only) #132494
      
      p4my0 <- complex_models_only %>%
        filter(Num_Params == 4) %>%
        filter(Missing_Yrs == 0) %>%
        arrange(adjaic2) %>%
        head(100)
      
      nrow(p4my0)
      
      p4my0 %>% head()

      #separate groups of models by Num_Param and Missing_Yrs
      p4my0_popularity <- p4my0 %>%
        head(100) %>%
        separate_rows(Predictors, sep = " \\+ ") %>%
        group_by(Predictors) %>%
        summarise(Times_Selected_In_Top_100_Dredged = n(), .groups = "drop") %>%
        arrange(desc(Times_Selected_In_Top_100_Dredged)) %>%
        filter(Times_Selected_In_Top_100_Dredged>=10)
      
      print(p4my0_popularity,n=Inf)
      
      p4my0[,1] %>% head()
      
      # Predictors                  Times_Selected_In_Top_100_Dredged
      # 1 X09_DFA_HakeAge5Plus                                       90
      # 2 X06_Lmu_IntSprJunH                                         69
      # 3 X04_Sablefish__JSOES                                       27
      # 4 X12_DFA_biomassEuphShelfSum                                20
      # 5 X13_ammod_WAI                                              20
      # 6 X08_DFA_DC_corm_3_WS                                       18
      # 7 X15_sablefishBiomass_predAK                                17
      # 8 X15_salmonSharkGoA_predAK                                  14
      # 9 X12_copepodCom_EGoA                                        13
      # 10 X13_gadid_WAI                                              12
      # 11 X08_Large_gulls_7_WS                                       11
      # 12 X11_DFA_Harbour_p_WS                                       10      
      
      #I think what's happening is that hake is standing in for cpue, and either sablefish or the next 2 preyAK items are standing in for growth.
      #so I think you do need something that is representing both in all of the best models.

#these should be sort of comparable to the DAG1B long      
      Formula                                                                                                
      1 X16_SAR ~ X04_Sablefish__JSOES + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus + X13_ammod_WAI           
      2 X16_SAR ~ X06_Lmu_IntSprJunH + X05_anchovy_GAM + X12_DFA_biomassEuphShelfSum + X12_copepodCom_EGoA     
      3 X16_SAR ~ X06_Lmu_IntSprJunH + X04_Sablefish__JSOES + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus      
      4 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum + X12_copepodCom_EGoA
      5 X16_SAR ~ X05_anchovy_GAM + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI         
      6 X16_SAR ~ X06_Lmu_IntSprJunH + X04_Sablefish__JSOES + X09_DFA_ChinAbundSnakeFall + X09_DFA_HakeAge5Plus
      
      
# --- 1. Get your trusted indicators ---
trusted_indicators <- ecological_synthesis %>%
  distinct(Indicator, Node, Explanation) %>%
  group_by(Indicator) %>%
  slice(1) %>%
  ungroup()

# --- 2. Calculate selection frequency in your NEW combined dredge ---
combined_dredge_popularity <- complex_models_only %>%
  # Zero in on the top 100 models from the combined stack
  head(100) %>%
  separate_rows(Predictors, sep = " \\+ ") %>%
  group_by(Predictors) %>%
  summarise(Times_Selected_In_Top_100_Dredged = n(), .groups = "drop")

print(combined_dredge_popularity,n=Inf)

# --- 3. Join them together ---
comparison_table <- trusted_indicators %>%
  left_join(combined_dredge_popularity, by = c("Indicator" = "Predictors")) %>%
  mutate(
    Times_Selected_In_Top_100_Dredged = replace_na(Times_Selected_In_Top_100_Dredged, 0)
  ) %>%
  # Arrange by trophic Node, then by how popular they were in the dredge
  arrange(Node, desc(Times_Selected_In_Top_100_Dredged))

# --- 4. Save & View the Results ---
print(as_tibble(comparison_table), n = Inf)
write_csv(comparison_table, "LisaXP/outputs_4/MainRegression_vs_CombinedDredge_Comparison.csv")