#186286
#130778
dredge_results<-rbind(dredge_results_X06X07,dredge_results_noX06X07) %>%
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

dredge_results_noHCI<- dredge_results %>%
  filter(!grepl("X01_habCompInd",Predictors))

write_csv(dredge_results_noHCI, "LisaXP/outputs_4/NCC_AK_dredge_results_combined_noHCI.csv")


#this script reexamines dregde results to
#1. remove HCI
#2. answer specific questions for Eve

#Input files
    dredge_results<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names = NULL);head(dredge_results)
    adjusted_sorting<-read.csv(file="LisaXP/outputs_4/adjusted_aic_dredge_results.csv",row.names=FALSE)
    head(adjusted_sorting)
    head(dredge_results)
    head(guilds)
#no HCI----------    
    dredge_results_noHCI<- dredge_results %>%
      filter(!grepl("X01_habCompInd",Predictors))
    
    adjusted_sorting<- dredge_results_noHCI %>%
      # Apply your +1 penalty per missing year to create a custom sorting index
      mutate(Adjusted_AIC = Raw_AIC + (1 * Missing_Yrs))

#3 & 4 param models only---------
# 1. Isolate the complex models using your preferred flat penalty (e.g., Mult_0.4)
complex_models_only <- dredge_results_noHCI %>%
  filter(Num_Params %in% c(3, 4)) %>%
  mutate(
    # Applying the stable flat penalty at your sensitivity sweet spot
    Adjusted_AIC = Raw_AIC + (1.0  * Missing_Yrs),
    Adj_AIC2 = Raw_AIC + (2.0 *  Missing_Yrs)
  ) %>%
  # 2. Add 4 clean binary columns to check for specific guild presence
  mutate(
    PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
    PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
    PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
    PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
  ) %>%
  rowwise() %>%
  mutate(Nguilds=sum(c_across(PreyNCC:PredAK))) %>%
  
  arrange(Adj_AIC2)

    
    print(complex_models_only)
    complex_models_only[1,1]
    
    
    #Subsets
    #top5 models by complete_Yrs
   top1_adj <- complex_models_only %>%
      arrange(Adjusted_AIC) %>%
    group_by(Complete_Yrs) %>%
      slice(1) %>% 
      ungroup() 
      
    print(top1_adj) #same as raw aic, bc we have grouped by missing years and the adjustment add 1/yr, so it doesn't reorder within those groups
    
    # 1 X16_SAR ~ X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X09_canaryRockfish
    # 2 X16_SAR ~ X01_NHLlogSum_win_05 + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X13_ammod_WAI                         
    # 3 X16_SAR ~ X01_DFA_sumPreyOfPrey_planktonJun + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X10_Harbour_s_2yrLead_WS 
    # 4 X16_SAR ~ X12_DFA_biomassEuphShelfSum + X13_ammod_WAI + X13_hexagram_EAI + X13_hexagram_WAI               
    #
    
    #     Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources                      Adjusted_AIC PreyNCC PredNCC PreyAK PredAK Nguilds
    # 1          4    20.1           19           5 PreyNCC, PreyNCC, PredNCC, PredNCC         25.1       1       1      0      0       2
    # 2          4    29.4           21           3 PreyNCC, PreyNCC, PredNCC, PreyAK          32.4       1       1      1      0       3
    # 3          4    40.9           23           1 PreyNCC, PreyNCC, PredNCC, PredAK          41.9       1       1      0      1       3
    # 4          4    30.1           24           0 PreyAK, PreyAK, PreyAK, PreyAK             30.1       0       0      1      0       1
    
    #re 1 missing yrs, you would have to add 10 aic/yr for the 23 yr time series to beat the 24 yr time series, so that is not realistic. we can eliminate the JSOES sum of prey one
    #re 3 missing yrs, it would take a 0 penalty over 3 missing years to beat the 24 yr time series, so that is not realistic. we can eliminate the X01_NHLlogSum_win_05
    #re 5 missing years, adding aic=2/yr would have made that worse than 0 missing years, but 1/yr keeps it better. So those seem to be better fits for real
    
    top1_200 <- 
      complex_models_only %>%
      filter(Adjusted_AIC<35) %>%
      arrange(Missing_Yrs)
    top1_200

    stem(top1_200$Missing_Yrs) #-- 5 good long models, a couple of 3missing, then all the 5 missing years
    stem(top1_200$Num_Params)
    
    top1_100 <- 
      complex_models_only %>%
      mutate(
        Adj_AIC2 = Raw_AIC + (2.0  * Missing_Yrs)
        )%>%
        filter(Adj_AIC2<35) %>%
      arrange(Missing_Yrs)
    top1_100 #actually only 32 meet this criterion, the 3missing yrs dropped out, but all the original 5 long models and the same top 5 missing yrs in there
    
    
    print(top1_100,n=Inf)
    
    top1_100[,1]
    # 1 X16_SAR ~ X12_DFA_biomassEuphShelfSum + X13_ammod_WAI + X13_hexagram_EAI + X13_hexagram_WAI                                
    # 2 X16_SAR ~ X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X12_copepodCom_WGoA + X13_ammod_WAI                         
    # 3 X16_SAR ~ X04_Sablefish__JSOES + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus + X13_ammod_WAI                               
    # 4 X16_SAR ~ X09_DFA_HakeAge5Plus + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI                        
    # 5 X16_SAR ~ X05_anchovy_GAM + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI                             
    # 6 X16_SAR ~ X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X09_canaryRockfish      
    # 7 X16_SAR ~ X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X09_chilipepper         
    # 8 X16_SAR ~ X12_DFA_biomassEuphShelfSum + X12_copepodCom_EGoA + X15_DFA_sleeperSharkBSAI_predAK + X15_spinyDogfishBSAI_predAK
    # 9 X16_SAR ~ X12_DFA_biomassEuphShelfSum + X13_hexagram_EAI + X15_DFA_sleeperSharkBSAI_predAK + X15_spinyDogfishBSAI_predAK   
    # 10 X16_SAR ~ X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X09_canaryRockfish + X15_DFA_sleeperSharkBSAI_predAK               

    
    #so what the hell is going on in these models? 
    #we could require hake as the predNCC. that would cover JSOES, and they are in a lot of the top models
    
    nrow(complex_models_only) #129371
#hake only-----------        
    hake_models<-complex_models_only %>% 
      filter(grepl("X09_DFA_HakeAge5Plus",Predictors)) %>%
      arrange(Adj_AIC2)
    
    nrow(hake_models) #8545
    
    top100_hake_models<- hake_models[1:100,]
    
    t(names(top100_hake_models))
    
    apply(top100_hake_models[,c(3,5,6,4,8:14)],2,sum)/100
    # Num_Params Complete_Yrs  Missing_Yrs      Raw_AIC Adjusted_AIC     Adj_AIC2      PreyNCC      PredNCC       PreyAK       PredAK      Nguilds 
    # 3.98000     19.16000      4.84000     27.00835     31.84835     30.88035      0.80000      1.00000      0.31000      0.34000      2.45000 
    #all were required to have a predNCC; 
    # almost all of them got a preyNCC
    #there were none w/ all 4 guilds
    #about 1/3 had preyAK and 1/3 had predAK the remaining 1/3 had only NCC
    
    top100_hake_models$Predictors
    
    
    hist(top100_hake_models$Nguilds)
    
    top100_hake_models$Predictors
    
    
    #summary stats-----
    #Step 2: Summarizing Guild Representation in Top Models---------
    guild_representation_noHCI <- top100_hake_models %>%
      group_by(Num_Params) %>%
      # Isolate the top 20 best models for each parameter tier
      slice_min(order_by = Adjusted_AIC, n = 20) %>%
      
      # Calculate what percentage of these top models contain each guild
      summarize(
        Total_Top_Models = n(),
        Pct_Prey_NCC     = round(mean(PreyNCC) * 100, 1),
        Pct_Pred_NCC     = round(mean(PredNCC) * 100, 1),
        Pct_Prey_AK      = round(mean(PreyAK) * 100, 1),
        Pct_Pred_AK      = round(mean(PredAK) * 100, 1),
        .groups = "drop"
      )
    
    print(guild_representation_noHCI)
    
    #     Num_Params Total_Top_Models Pct_Prey_NCC Pct_Pred_NCC Pct_Prey_AK Pct_Pred_AK
    # 1          3                2          100          100           0           0
    # 2          4               20           90          100          25          40
        
    
    #Step 3: Tracking "Winning" Predictor Consistency Across Tiers------
    # Unnest the predictors column to score individual variable popularity
    variable_popularity_summary <- top100_hake_models %>%
      # group_by(Num_Params) %>%
      # slice_min(order_by = Adjusted_AIC, n = 100) %>%
      # ungroup() %>%
      # Separate "var1 + var2" strings into individual rows
      separate_rows(Predictors, sep = " \\+ ") %>%
      
      # Count occurrences
      group_by(Predictors) %>%
      summarize(Times_Selected_In_Top100 = n(), .groups = "drop") %>%
      
      # Sort to see the heaviest hitters float to the top of each parameter tier
      arrange( desc(Times_Selected_In_Top100)) 
    
    # View the top winning variables per tier
    print(variable_popularity_summary, n = 30)

    # Predictors                                Times_Selected_In_Top100
    # 1 X09_DFA_HakeAge5Plus                                           100
    # 2 X04_Sablefish__JSOES                                            78
    # 3 X09_chilipepper                                                 51
    # 4 X09_canaryRockfish                                              45
    # 5 X15_sablefishBiomass_predAK                                     11
    # 6 X12_DFA_biomassEuphShelfSum                                      7
    # 7 X13_ammod_WAI                                                    7
    # 8 X15_DFA_sleeperSharkBSAI_predAK                                  7
    # 9 X15_halibutBiomassAge8plus_2yrLead_predAK                        5
    # 10 X10_Californian_s_l_2yrLead_WS                                   4
    # 11 X11_DFA_Harbour_p_WS                                             4
    # 12 X12_copepodCom_EGoA                                              4
    # 13 X15_salmonSharkGoA_predAK                                        4
    # 14 X01_NHLlogSum_win_05                                             3
    # 15 X01_NHLlogSum_win_15                                             3
    # 16 X01_NHLlogSum_win_25                                             3
    # 17 X02_DFA_SeaNettle_                                               3
    # 18 X04_PacificPompano_                                              3
    # 19 X08_DFA_DC_corm_3_WS                                             3
    # 20 X11_Harbor_seal_CR                                               3
    # 21 X13_gadid_WAI                                                    3
    # 22 X15_ArrowtoothFlounderBiomass_predAK                             3
    # 23 X01_DFA_sumPreyOfPrey_planktonJun                                2
    # 24 X01_Red_phal_WS                                                  2
    # 25 X01_Sabine_gull_WS                                               2
    # 26 X02.krill                                                        2
    # 27 X03_DFA_comMurreDietHerrSard_Yaquina_NCC                         2
    # 28 X03_HakeAge1                                                     2
    # 29 X04_marketsquid_GAM                                              2
    # 30 X05_DFA_abundSardine                                             2
    
    
    #I am really suspicious of those straggler rockfish
#Remover chilipepper and canary rockfish-----------
   
    hake_norockfish_models<-complex_models_only %>% 
      filter(grepl("X09_DFA_HakeAge5Plus",Predictors)) %>%
      filter(!grepl("X09_chilipepper|X09_canaryRockfish",Predictors)) %>%
      arrange(Adj_AIC2)
    
    hake_norockfish_models
    hake_norockfish_models[,1]
    # 1 X16_SAR ~ X04_Sablefish__JSOES + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus + X13_ammod_WAI                  
    # 2 X16_SAR ~ X01_NHLlogSum_win_05 + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X13_ammod_WAI                  
    # 3 X16_SAR ~ X01_NHLlogSum_win_25 + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X13_ammod_WAI                  
    # 4 X16_SAR ~ X01_NHLlogSum_win_15 + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X13_ammod_WAI                  
    # 5 X16_SAR ~ X09_DFA_HakeAge5Plus + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI           
    # 6 X16_SAR ~ X08_Large_gulls_7_WS + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X13_gadid_WAI                         
    # 7 X16_SAR ~ X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI + X13_hexagram_EAI               
    # 8 X16_SAR ~ X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X13_hexagram_WAI                      
    # 9 X16_SAR ~ X08_Large_gulls_7_WS + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X15_sablefishBiomass_predAK           
    # 10 X16_SAR ~ X01_NHLlogSum_win_15 + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X15_DFA_sleeperSharkBSAI_predAK
    
    #   Formula              Predictors Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources Adjusted_AIC Adj_AIC2 PreyNCC PredNCC PreyAK PredAK Nguilds
    # 1 X16_SAR ~ X04_Sable… X04_Sable…          4    30.3           24           0 PreyNCC, Pre…         30.3     30.3       1       1      1      0       3
    # 2 X16_SAR ~ X01_NHLlo… X01_NHLlo…          4    29.4           21           3 PreyNCC, Pre…         32.4     31.8       1       1      1      0       3
    # 3 X16_SAR ~ X01_NHLlo… X01_NHLlo…          4    30.0           21           3 PreyNCC, Pre…         33.0     32.4       1       1      1      0       3
    # 4 X16_SAR ~ X01_NHLlo… X01_NHLlo…          4    30.0           21           3 PreyNCC, Pre…         33.0     32.4       1       1      1      0       3
    # 5 X16_SAR ~ X09_DFA_H… X09_DFA_H…          4    33.3           24           0 PredNCC, Pre…         33.3     33.3       0       1      1      0       2
    # 6 X16_SAR ~ X08_Large… X08_Large…          4    35.9           24           0 PredNCC, Pre…         35.9     35.9       0       1      1      0       2
    # 7 X16_SAR ~ X09_DFA_H… X09_DFA_H…          4    35.9           24           0 PredNCC, Pre…         35.9     35.9       0       1      1      0       2
    # 8 X16_SAR ~ X04_Sable… X04_Sable…          4    36.2           24           0 PreyNCC, Pre…         36.2     36.2       1       1      1      0       3
    # 9 X16_SAR ~ X08_Large… X08_Large…          4    36.4           24           0 PredNCC, Pre…         36.4     36.4       0       1      1      1       3
    # 10 X16_SAR ~ X01_NHLlo… X01_NHLlo…         4    34.1           21           3 PreyNCC, Pre…         37.1     36.5       1       1      0      1       3
    
    guild_representation_noHCI_norockfish <- hake_norockfish_models %>%
      group_by(Num_Params) %>%
      # Isolate the top 20 best models for each parameter tier
      slice_min(order_by = Adj_AIC2, n = 100) %>%
      
      # Calculate what percentage of these top models contain each guild
      summarize(
        Total_Top_Models = n(),
        Pct_Prey_NCC     = round(mean(PreyNCC) * 100, 1),
        Pct_Pred_NCC     = round(mean(PredNCC) * 100, 1),
        Pct_Prey_AK      = round(mean(PreyAK) * 100, 1),
        Pct_Pred_AK      = round(mean(PredAK) * 100, 1),
        .groups = "drop"
      )
    
    print(guild_representation_noHCI_norockfish)
    
    #     Num_Params Total_Top_Models Pct_Prey_NCC Pct_Pred_NCC Pct_Prey_AK Pct_Pred_AK
    # 1          3              100           72          100          49          39
    # 2          4              100           79          100          64          39
         
    
    #Step 3: Tracking "Winning" Predictor Consistency Across Tiers------
    # Unnest the predictors column to score individual variable popularity
    variable_popularity_summary_norockfish <- hake_norockfish_models %>%
          group_by(Num_Params) %>%
          slice_min(order_by = Adjusted_AIC, n = 100) %>%
          ungroup() %>%
#      Separate "var1 + var2" strings into individual rows
      separate_rows(Predictors, sep = " \\+ ") %>%
      
      # Count occurrences
      group_by(Predictors) %>%
      summarize(Times_Selected_In_Top100 = n(), .groups = "drop") %>%
      
      # Sort to see the heaviest hitters float to the top of each parameter tier
      arrange( desc(Times_Selected_In_Top100)) 
    
    # View the top winning variables per tier
    print(variable_popularity_summary_norockfish, n = Inf)
    
    # Predictors                                Times_Selected_In_Top100
    # 1 X09_DFA_HakeAge5Plus                                           200
    # 2 X04_Sablefish__JSOES                                            93
    # 3 X13_ammod_WAI                                                   71
    # 4 X01_NHLlogSum_win_25                                            27
    # 5 X01_NHLlogSum_win_15                                            25
    # 6 X10_Harbour_s_2yrLead_WS                                        23
    # 7 X13_gadid_WAI                                                   21
    # 8 X15_sablefishBiomass_predAK                                     20
    # 9 X12_DFA_biomassEuphShelfSum                                     19
    # 10 X08_DFA_DC_corm_3_WS                                            18
    # 11 X01_NHLlogSum_win_05                                            15
    # 12 X15_DFA_sleeperSharkBSAI_predAK                                 12
    # 13 X13_hexagram_EAI                                                10
    # 14 X13_hexagram_WAI                                                10
    # 15 X04_PacificPompano_                                              9
    # 16 X08_Large_gulls_7_WS                                             9
    # 17 X11_DFA_Harbour_p_WS                                             8
    # 18 X02_DFA_SeaNettle_                                               6
    # 19 X03_DFA_comMurreDietHerrSard_Yaquina_NCC                         6
    # 20 X10_Harbor_seal_CR_2yrLead                                       6
    # 21 X12_copepodCom_EGoA                                              6
    # 22 X15_halibutBiomassAge8plus_2yrLead_predAK                        6
    # 23 X09_DFA_ChinAbundSnakeFall                                       5
    # 24 X13_DFA_WGOA_DFA_seabirds                                        5
    # 25 X14_pinkSalmon                                                   5
    # 26 X14_pinkSalmonAsia                                               5
    # 27 X15_sharkCatchGoA_predAK                                         5
    # 28 X01_Red_phal_WS                                                  4
    # 29 X01_Sabine_gull_WS                                               4
    # 30 X02.krill                                                        4
    # 31 X10_Northern_f_s_2yrLead_WS                                      4
    # 32 X11_Harbor_seal_CR                                               4
    # 33 X15_salmonSharkGoA_predAK                                        4
    # 34 X01_DFA_sumPreyOfPrey_planktonJun                                3
    # 35 X01_Red_necked_phal_WS                                           3
    # 36 X10_DFA_ssl.est.wholerange_2yrLead                               3
    # 37 X13_pollockBiomassGoAage3plus_predAK                             3
    # 38 X14_pinkSalmonNorthAmerica                                       3
    # 39 X15_spinyDogfishGoA_predAK                                       3
    # 40 X03_HakeAge1                                                     2
    # 41 X04_marketsquid_GAM                                              2
    # 42 X05_DFA_abundSardine                                             2
    # 43 X05_herring_GAM                                                  2
    # 44 X05_herring_NCC                                                  2
    # 45 X15_ArrowtoothFlounderBiomass_predAK                             2
    # 46 X10_Californian_s_l_2yrLead_WS                                   1
    # >     