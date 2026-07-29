library(tidyverse)





# 1. Generate complex model list----------
dredge_results_noHCI<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results_combined_noHCI.csv",row.names = NULL);head(dredge_results)
head(dredge_results_noHCI)

complex_models_only <- dredge_results_noHCI %>%
  filter(Num_Params %in% c(3, 4)) %>%
  # mutate(
  #   # Applying the stable flat penalty at your sensitivity sweet spot
  #   Adjusted_AIC = Raw_AIC + (1.0  * Missing_Yrs),
  #   adjaic2 = Raw_AIC + (2.0 *  Missing_Yrs)
  # ) %>%
  # # 2. Add 4 clean binary columns to check for specific guild presence
  # mutate(
  #   PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
  #   PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
  #   PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
  #   PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
  # ) %>%
  # rowwise() %>%
  # mutate(Nguilds=sum(c_across(PreyNCC:PredAK))) %>%
  # ungroup() %>% # <-- CRITICAL FIX: Strip rowwise grouping here! 
  arrange(adjaic2)


guild_representation <- complex_models_only %>%
  group_by(Num_Params) %>%
  # Isolate the top 20 best models for each parameter tier
  slice_min(order_by = adjaic2, n = 100) %>%
  
  # Calculate what percentage of these top models contain each guild
  summarize(
    Total_Top_Models = n(),
    Pct_Prey_NCC     = round(mean(PreyNCC) * 100, 1),
    Pct_Pred_NCC     = round(mean(PredNCC) * 100, 1),
    Pct_Prey_AK      = round(mean(PreyAK) * 100, 1),
    Pct_Pred_AK      = round(mean(PredAK) * 100, 1),
    .groups = "drop"
  )

#print(guild_representation)-----------
# Num_Params Total_Top_Models Pct_Prey_NCC Pct_Pred_NCC Pct_Prey_AK Pct_Pred_AK
# 1          3              100           41           78          58          41
# 2          4              100           63           74          64          46


head(complex_models_only)

complex_models_only %>% filter(PreyAK==1) %>% select(Formula)
complex_models_only %>% filter(PredAK==1) %>% select(Formula)


# ==============================================================================
# SECTION 1: MASTER PROXY DEFINITION (CALCULATED ONCE FOR ALL SUBSETS)
# ==============================================================================
cat("Calculating Master Correlation Matrix and Proxy Vectors...\n")

# Calculate pairwise correlation matrix across ALL variables in guild.dfas1
cor_matrix_master <- cor(guild.dfas1, use = "pairwise.complete.obs")

cor_threshold <- 0.50

# Clean atomic vectors (using unname() to remove structural metadata)
growth_proxies    <- unname(names(which(abs(cor_matrix_master["X06_Lmu_IntSprJunH", ]) >= cor_threshold)))
abundance_proxies <- unname(names(which(abs(cor_matrix_master["X07_DFA_cpue_IntSprJunHW", ]) >= cor_threshold)))

cat(sprintf("Master Growth Proxies (|r| >= %.2f): %d variables\n", cor_threshold, length(growth_proxies)))
cat(sprintf("Master Abundance Proxies (|r| >= %.2f): %d variables\n\n", cor_threshold, length(abundance_proxies)))


# ==============================================================================
# SECTION 2: SUBSET EVALUATION PIPELINE
# ==============================================================================

# --- Helper Function to Test Any Subset ---
test_proxy_hypothesis <- function(model_data, model_set_label) {
  
  cat("=========================================================\n")
  cat(sprintf(" RUNNING HYPOTHESIS TEST FOR MODEL SET: %s\n", model_set_label))
  cat("=========================================================\n")
  
  tested_df <- model_data %>%
    mutate(
      # Split predictors cleanly
      Pred_List = str_split(Predictors, " \\+ "),
      
      # Robust scalar boolean evaluations per row
      Has_Growth_Proxy    = map_lgl(Pred_List, function(x) any(x %in% growth_proxies)),
      Has_Abundance_Proxy = map_lgl(Pred_List, function(x) any(x %in% abundance_proxies)),
      
      Meets_Both = Has_Growth_Proxy & Has_Abundance_Proxy
    )
  
  percent_supported <- mean(tested_df$Meets_Both) * 100
  
  print(summary(tested_df[, c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")]))
  
  cat(sprintf("\nRESULTS (%s): %.1f%% of Top 100 models contain BOTH proxies.\n", 
              model_set_label, percent_supported))
  cat("=========================================================\n\n")
  
  # Remove list column before saving to CSV
  output_df <- tested_df %>% select(-Pred_List)
  
  file_out <- sprintf("LisaXP/outputs_4/dredge_%s_chose_proxies_for_growth_abund.csv", model_set_label)
  write_csv(output_df, file_out)
  
  return(tested_df)
}


# ------------------------------------------------------------------------------
# MODEL SET 1: p4my0 (4 Params, 0 Missing Years)
# ------------------------------------------------------------------------------
          p4my0_dat <- complex_models_only %>%
            filter(Num_Params == 4, Missing_Yrs == 0) %>%
            arrange(adjaic2) %>%
            head(100)
          
          hypothesis_testing_p4my0 <- test_proxy_hypothesis(p4my0_dat, "p4my0")
          
          # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
          # Mode :logical    Mode :logical       Mode :logical  
          # FALSE:2          FALSE:20            FALSE:21       
          # TRUE :98         TRUE :80            TRUE :79       
          # 
          # RESULTS (p4my0): 79.0% of Top 100 models contain BOTH proxies.
          # 
# ------------------------------------------------------------------------------
# MODEL SET 2: p4my5 (4 Params, 5 Missing Years)
# ------------------------------------------------------------------------------
          p4my5_dat <- complex_models_only %>%
            filter(Num_Params == 4, Missing_Yrs == 5) %>%
            arrange(adjaic2) %>%
            head(100)
          
          hypothesis_testing_p4my5 <- test_proxy_hypothesis(p4my5_dat, "p4my5")
          
          # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
          # Mode :logical    Mode :logical       Mode :logical  
          # FALSE:26         FALSE:39            FALSE:65       
          # TRUE :74         TRUE :61            TRUE :35       
          # 
          # RESULTS (p4my5): 35.0% of Top 100 models contain BOTH proxies.
# ------------------------------------------------------------------------------
# MODEL SET 3: p3my0 (3 Params, 0 Missing Years)
# ------------------------------------------------------------------------------
          p3my0_dat <- complex_models_only %>%
            filter(Num_Params == 3, Missing_Yrs == 0) %>%
            arrange(adjaic2) %>%
            head(100)
          
          hypothesis_testing_p3my0 <- test_proxy_hypothesis(p3my0_dat, "p3my0")
          
          # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
          # Mode :logical    Mode :logical       Mode :logical  
          # FALSE:29         FALSE:16            FALSE:38       
          # TRUE :71         TRUE :84            TRUE :62       
          # 
          # RESULTS (p3my0): 62.0% of Top 100 models contain BOTH proxies.
          
          
          
# ------------------------------------------------------------------------------
# MODEL SET 4: p3my5 (3 Params, 5 Missing Years)
# ------------------------------------------------------------------------------
          p3my5_dat <- complex_models_only %>%
            filter(Num_Params == 3, Missing_Yrs == 5) %>%
            arrange(adjaic2) %>%
            head(100)
          
          hypothesis_testing_p3my5 <- test_proxy_hypothesis(p3my5_dat, "p3my5")
          
          # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
          # Mode :logical    Mode :logical       Mode :logical  
          # FALSE:54         FALSE:55            FALSE:79       
          # TRUE :46         TRUE :45            TRUE :21       
          # 
          # RESULTS (p3my5): 21.0% of Top 100 models contain BOTH proxies.
          # 
          

#Top models========
          p4my0_dat[1:10,1]
          p3my0_dat[1:10,1]
#exclude all stragglers?------
         complex_models_DFAs_only <- complex_models_only %>%
            # 1. Split the "var1 + var2" string into a list vector of individual variable names
            mutate(Pred_List = str_split(Predictors, " \\+ ")) %>%
            
            # 2. Filter rows where EVERY non-X06/X07 predictor contains "DFA"
            filter(
              map_lgl(Pred_List, function(preds) {
                # Keep only predictors that do NOT start with X06 or X07
                env_preds <- preds[!grepl("^X06|^X07", preds)]
                
                # If all non-mediator predictors are empty (e.g. model only had X06/X07), decide if TRUE/FALSE
                if (length(env_preds) == 0) return(FALSE) 
                
                # Check if ALL remaining environmental predictors contain "DFA"
                all(grepl("DFA", env_preds))
              })
            ) %>%
            
            # 3. Clean up the temporary list-column
            select(-Pred_List)          
 
          complex_models_DFAs_only[1:10,1]
          
          # 1 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_HakeAge5Plus + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum                    
          # 2 X16_SAR ~ X06_Lmu_IntSprJunH + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum                    
          # 3 X16_SAR ~ X06_Lmu_IntSprJunH + X02_DFA_SeaNettle_ + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus                             
          # 4 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum                                           
          # 5 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X13_DFA_WGOA_DFA_seabirds         
          # 6 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum + X13_DFA_WGOA_DFA_seabirds               
          # 7 X16_SAR ~ X06_Lmu_IntSprJunH + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus                                                  
          # 8 X16_SAR ~ X06_Lmu_IntSprJunH + X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum
          # 9 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_ChinAbundSnakeFall + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum              
          # 10 X16_SAR ~ X06_Lmu_IntSprJunH + X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus      
          # 
          
          complex_models_DFAs_only %>% head()
          
  guild_representation <- complex_models_DFAs_only %>%
            group_by(Num_Params) %>%
            # Isolate the top 20 best models for each parameter tier
            slice_min(order_by = adjaic2, n = 100) %>%
            
            # Calculate what percentage of these top models contain each guild
            summarize(
              Total_Top_Models = n(),
              Pct_Prey_NCC     = round(mean(PreyNCC) * 100, 1),
              Pct_Pred_NCC     = round(mean(PredNCC) * 100, 1),
              Pct_Prey_AK      = round(mean(PreyAK) * 100, 1),
              Pct_Pred_AK      = round(mean(PredAK) * 100, 1),
              .groups = "drop"
            )
          
          print(guild_representation)
          guild_representation_DFA<-guild_representation
          #   Num_Params Total_Top_Models Pct_Prey_NCC Pct_Pred_NCC Pct_Prey_AK Pct_Pred_AK
          # 1          3              100           34           79          57          21
          # 2          4              100           38           86          84          23  
          
          variable_popularity <- complex_models_DFAs_only %>%
            group_by(Num_Params,Missing_Yrs) %>%
            slice_min(order_by = adjaic2, n = 100) %>%
            ungroup() %>%
            #      Separate "var1 + var2" strings into individual rows
            separate_rows(Predictors, sep = " \\+ ") %>%
            
            # Count occurrences
            group_by(Predictors) %>%
            summarize(Times_Selected_In_Top100 = n(), .groups = "drop") %>%
            
            # Sort to see the heaviest hitters float to the top of each parameter tier
            arrange( desc(Times_Selected_In_Top100)) 
          
          # View the top winning variables per tier
          print(variable_popularity, n = Inf)
          
          Predictors                               Times_Selected_In_Top100
          1 X09_DFA_HakeAge5Plus                                          159
          2 X12_DFA_biomassEuphShelfSum                                   148
          3 X06_StomFull_May                                              144
          4 X13_DFA_WGOA_DFA_seabirds                                     135
          5 X03_DFA_comMurreDietHerrSard_Yaquina_NCC                      127
          6 X07_DFA_cpue_IntSprJunHW                                      118
          7 X11_DFA_Harbour_p_WS                                          104
          8 X08_DFA_DC_corm_3_WS                                           84
          9 X06_Lmu_IntSprJunH                                             77
          10 X09_DFA_ChinAbundSnakeFall                                     51
          11 X06_Lmu_IntSprMayW                                             49
          12 X15_DFA_sleeperSharkBSAI_predAK                                38
          13 X10_DFA_ssl.est.wholerange_2yrLead                             36
          14 X06_DFA_IGF_mu                                                 29
          15 X01_DFA_sumPreyOfPrey_planktonJun                              28
          16 X02_DFA_SeaNettle_                                             27
          17 X13_DFA_WGOA_DFA_midTrophic                                    24
          18 X05_DFA_abundSardine                                           22
          
        #separate tables-
          library(tidyverse)
          
          variable_popularity_grid <- complex_models_DFAs_only %>%
            # 1. Filter specifically for 3 & 4 params AND 0 & 5 missing years
            filter(Num_Params %in% c(3, 4), Missing_Yrs %in% c(0, 5)) %>%
            
            # 2. Grab the top 100 models by adjaic2 WITHIN each of the 4 combinations
            group_by(Num_Params, Missing_Yrs) %>%
            slice_min(order_by = adjaic2, n = 100, with_ties = FALSE) %>%
            ungroup() %>%
            
            # 3. Unnest individual predictors from the formula string
            separate_rows(Predictors, sep = " \\+ ") %>%
            
            # 4. Count occurrences WITHIN each of the 4 parameter/missing year tiers
            group_by(Num_Params, Missing_Yrs, Predictors) %>%
            summarize(Times_Selected = n(), .groups = "drop") %>%
            
            # 5. Create a clean Label for the 4 combinations (e.g., "P3_MY0", "P4_MY5")
            mutate(Group = paste0("P", Num_Params, "_MY", Missing_Yrs)) %>%
            select(-Num_Params, -Missing_Yrs) %>%
            
            # 6. Pivot so each of the 4 combinations gets its own column
            pivot_wider(
              names_from = Group,
              values_from = Times_Selected,
              values_fill = 0
            ) %>%
            
            # 7. Add a Total count across all 4 groups and sort by heaviest hitters
            rowwise() %>%
            mutate(Total_Selection_Count = sum(c_across(-Predictors))) %>% ungroup() %>%
            arrange(desc(Total_Selection_Count))
          
          # View the formatted table
          print(variable_popularity_grid, n = Inf)
          
          complex_models_DFAs_only %>%
            filter(Num_Params %in% c(3, 4), Missing_Yrs %in% c(0, 5)) %>%
            count(Num_Params, Missing_Yrs)
          
          # Save the 4-combination table
          write_csv(variable_popularity_grid, "LisaXP/outputs_4/DFA_Only_Variable_Popularity_0MY_3or4Param.csv")
          
          
          #   Predictors                               P3_MY0 P4_MY0 Total_Selection_Count
          # 1 X12_DFA_biomassEuphShelfSum                  35     56                    91
          # 2 X13_DFA_WGOA_DFA_seabirds                    35     52                    87
          # 3 X09_DFA_HakeAge5Plus                         40     45                    85
          # 4 X07_DFA_cpue_IntSprJunHW                     35     32                    67
          # 5 X11_DFA_Harbour_p_WS                         25     42                    67
          # 6 X03_DFA_comMurreDietHerrSard_Yaquina_NCC     22     41                    63
          # 7 X08_DFA_DC_corm_3_WS                         23     33                    56
          # 8 X06_Lmu_IntSprJunH                           20     33                    53
          # 9 X09_DFA_ChinAbundSnakeFall                   13     18                    31
          # 10 X10_DFA_ssl.est.wholerange_2yrLead           11     11                    22
          # 11 X15_DFA_sleeperSharkBSAI_predAK              11     11                    22
          # 12 X02_DFA_SeaNettle_                            9      8                    17
          # 13 X06_DFA_IGF_mu                                9      7                    16
          # 14 X13_DFA_WGOA_DFA_midTrophic                   6      8                    14
          # 15 X05_DFA_abundSardine                          6      3                     9
          
          
          
          library(tidyverse)
          
          library(tidyverse)
#Sort by guilds-------          
          # --- 1. Dynamically Create the Predictor-to-Guild Map from your Active Model List ---
          predictor_guild_map <- complex_models_only %>%
            # Keep only the formula and the guild metadata
            select(Predictors, Guild_Sources) %>%
            distinct() %>%
            
            # Separate predictors and guild sources into parallel rows
            separate_rows(Predictors, sep = " \\+ ") %>%
            separate_rows(Guild_Sources, sep = ", ") %>%
            
            # Group by Predictor to get its distinct primary Guild
            group_by(Predictors) %>%
            summarize(Guild = first(Guild_Sources), .groups = "drop")
          
          
          # --- 2. Build and Order the 4-Combination Table by Guild ---
          variable_popularity_grid_by_guild <- complex_models_DFAs_only %>%
            filter(Num_Params %in% c(3, 4), Missing_Yrs %in% c(0, 5)) %>%
            
            # Grab top 100 models by adjaic2
            group_by(Num_Params, Missing_Yrs) %>%
            slice_min(order_by = adjaic2, n = 100, with_ties = FALSE) %>%
            ungroup() %>%
            
            # Unnest individual predictors
            separate_rows(Predictors, sep = " \\+ ") %>%
            
            # Count occurrences
            group_by(Num_Params, Missing_Yrs, Predictors) %>%
            summarize(Times_Selected = n(), .groups = "drop") %>%
            
            # Ensure all 4 grid combinations exist
            complete(
              nesting(Predictors), 
              Num_Params = c(3, 4), 
              Missing_Yrs = c(0, 5), 
              fill = list(Times_Selected = 0)
            ) %>%
            
            # Format column labels
            mutate(Group = paste0("P", Num_Params, "_MY", Missing_Yrs)) %>%
            select(-Num_Params, -Missing_Yrs) %>%
            
            # Pivot out to 4 distinct columns
            pivot_wider(
              names_from = Group,
              values_from = Times_Selected,
              values_fill = 0
            ) %>%
            
            # Calculate total selection count safely
            rowwise() %>%
            mutate(Total_Selection_Count = sum(c_across(-Predictors))) %>%
            ungroup() %>%
            
            # --- 3. Join the Dynamic Guild Map & Order ---
            left_join(predictor_guild_map, by = "Predictors") %>%
            
            # Set custom factor levels for ecological sorting
            mutate(Guild = factor(Guild, levels = c("Growth", "Abundance", "PreyNCC", "PredNCC", "PreyAK", "PredAK"))) %>%
            
            # Arrange first by Guild, then by heaviest hitters within each Guild
            arrange(Guild, desc(Total_Selection_Count)) %>%
            
            # Relocate Guild column right next to Predictors
            relocate(Guild, .after = Predictors)
          
          # View formatted output
          print(variable_popularity_grid_by_guild, n = Inf)
          
          # Save the guild-ordered table
          write_csv(variable_popularity_grid_by_guild, "LisaXP/outputs_4/DFA_Only_Variable_Popularity_By_Guild.csv")
          
          
          # Predictors                               Guild     P3_MY0 P3_MY5 P4_MY0 P4_MY5 Total_Selection_Count
          # <chr>                                    <fct>      <int>  <int>  <int>  <int>                 <int>
          #   1 X06_Lmu_IntSprJunH                       Growth        20      0     33      0                    53
          # 2 X09_DFA_ChinAbundSnakeFall               Growth        13      0     18      0                    31
          # 3 X06_DFA_IGF_mu                           Growth         9      0      7      0                    16
          # 4 X07_DFA_cpue_IntSprJunHW                 Abundance     35      0     32      0                    67
          # 5 X13_DFA_WGOA_DFA_seabirds                PreyNCC       35      0     52      0                    87
          # 6 X09_DFA_HakeAge5Plus                     PreyNCC       40      0     45      0                    85
          # 7 X03_DFA_comMurreDietHerrSard_Yaquina_NCC PreyNCC       22      0     41      0                    63
          # 8 X08_DFA_DC_corm_3_WS                     PreyNCC       23      0     33      0                    56
          # 9 X02_DFA_SeaNettle_                       PreyNCC        9      0      8      0                    17
          # 10 X05_DFA_abundSardine                     PreyNCC        6      0      3      0                     9
          # 11 X11_DFA_Harbour_p_WS                     PredNCC       25      0     42      0                    67
          # 12 X10_DFA_ssl.est.wholerange_2yrLead       PredNCC       11      0     11      0                    22
          # 13 X13_DFA_WGOA_DFA_midTrophic              PredNCC        6      0      8      0                    14
          # 14 X12_DFA_biomassEuphShelfSum              PreyAK        35      0     56      0                    91
          # 15 X15_DFA_sleeperSharkBSAI_predAK          PreyAK        11      0     11      0                    22
          # > 