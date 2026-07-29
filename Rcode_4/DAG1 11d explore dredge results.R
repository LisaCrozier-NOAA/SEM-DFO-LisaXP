

#Block 3: 
#1. playing with adjusting aic to adequately penalize shorter models, that always rose to the top,  
#2. exploring variable importance ("popularity")
#3. comparing guild performance


dredge_results<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names = NULL)

adjusted_sorting <- dredge_results %>%
  # Apply your +2 penalty per missing year to create a custom sorting index
  mutate(Adjusted_AIC = Raw_AIC + (2 * Missing_Yrs)) %>%
  
  # Example Post-Hoc Filtering using your new metadata column:
  # Let's say you only want models that have at least one Alaska Predator
  filter(grepl("PredAK", Guild_Sources)) %>% 
  
  # Sort by your new adjusted metric
  arrange(Adjusted_AIC)


# 1. Identify the top 10 complete models for each parameter count (1 to 4)
top_24yr_models <- dredge_results %>%
  filter(Complete_Yrs == 24) %>%
  group_by(Num_Params) %>%
  slice_min(order_by = Raw_AIC, n = 10) %>%
  ungroup() %>%
  mutate(Model_ID = row_number()) 

p1<-top_24yr_models %>% filter(Num_Params==1)
summary()
print(p1,n=10)

print(top_19yr_models ,n=Inf)
print(top_19yr_models %>% filter(Num_Params==1),n=Inf)

library(tidyverse)


#========================
#Adjust aic-------
# Define your exact empirical grand mean penalties
penalty_lookup <- c("1" = 3.0, "2" = 2.5, "3" = 1.5, "4" = 1.0)

adjusted_sorting <- dredge_results %>%
  # 1. Apply the parameter-specific penalty dynamically
  mutate(
    Penalty_Rate = penalty_lookup[as.character(Num_Params)],
    Adjusted_AIC = Raw_AIC + (Penalty_Rate * Missing_Yrs)
  ) %>%
  
  # 2. Add 4 clean binary columns to check for specific guild presence
  mutate(
    PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
    PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
    PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
    PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
  ) %>%
  
  # 3. Sort by your new adjusted metric
  arrange(Adjusted_AIC)


adjusted_sorting
write.csv(adjusted_sorting,file="LisaXP/outputs_4/adjusted_aic_dredge_results.csv",row.names=FALSE)

#Step 2: Summarizing Guild Representation in Top Models---------
guild_representation_summary <- adjusted_sorting %>%
  # Group by parameter count to keep comparisons fair
  group_by(Num_Params) %>%
  # Isolate the top 20 best models for each parameter tier
  slice_min(order_by = Adjusted_AIC, n = 20) %>%
  
  # Calculate what percentage of these top models contain each guild
  summarize(
    Total_Top_Models = n(),
    Pct_Prey_NCC     = round(mean(Has_PreyNCC) * 100, 1),
    Pct_Pred_NCC     = round(mean(Has_PredNCC) * 100, 1),
    Pct_Prey_AK      = round(mean(Has_PreyAK) * 100, 1),
    Pct_Pred_AK      = round(mean(Has_PredAK) * 100, 1),
    .groups = "drop"
  )

print(guild_representation_summary)

#It looks like Prey_NCC play superbly when combined with Pred_NCC and Pred_AK in 3 & 4 parameter models, but they do worse or = than any other node on their own (1 & 2 param models)


#   Num_Params Total_Top_Models   Pct_Prey_NCC Pct_Pred_NCC Pct_Prey_AK Pct_Pred_AK
# 1          1               20           20           25          35          20
# 2          2               20           40           60          45          50
# 3          3               20           95           60          15          80
# 4          4               20           80           70          35          65


#Step 3: Tracking "Winning" Predictor Consistency Across Tiers------
# Unnest the predictors column to score individual variable popularity
variable_popularity_summary <- adjusted_sorting %>%
  group_by(Num_Params) %>%
  slice_min(order_by = Adjusted_AIC, n = 20) %>%
  ungroup() %>%
  # Separate "var1 + var2" strings into individual rows
  separate_rows(Predictors, sep = " \\+ ") %>%
  
  # Count occurrences
  group_by(Num_Params, Predictors) %>%
  summarize(Times_Selected_In_Top20 = n(), .groups = "drop") %>%
  
  # Sort to see the heaviest hitters float to the top of each parameter tier
  arrange( desc(Times_Selected_In_Top20)) %>%
  #arrange(Num_Params, desc(Times_Selected_In_Top20))
  filter(Times_Selected_In_Top20>2)

# View the top winning variables per tier
print(variable_popularity_summary, n = 30)


# Num_Params Predictors                                Times_Selected_In_Top20
# 1          3 X01_habCompInd                                                 17
# 2          4 X01_habCompInd                                                 13
# 3          4 X04_Sablefish__JSOES                                           12
# 4          4 X09_canaryRockfish                                              9
# 5          3 X09_canaryRockfish                                              7
# 6          4 X15_DFA_sleeperSharkBSAI_predAK                                 7
# 7          2 X01_habCompInd                                                  6
# 8          2 X09_DFA_HakeAge5Plus                                            6
# 9          2 X12_copepodBiomass_WGoA                                         6
# 10          4 X03_DFA_comMurreDietHerrSard_Yaquina_NCC                        6
# 11          3 X15_spinyDogfishBSAI_predAK                                     5
# 12          3 X15_ArrowtoothFlounderBiomass_predAK                            4
# 13          3 X15_halibutBiomassAge8plus_2yrLead_predAK                       4
# 14          3 X15_salmonSharkBSAI_predAK                                      4
# 15          4 X09_chilipepper                                                 4
# 16          4 X12_DFA_biomassEuphShelfSum                                     4
# 17          4 X15_salmonSharkBSAI_predAK                                      4
# 18          4 X15_spinyDogfishBSAI_predAK                                     4
# 19          2 X09_canaryRockfish                                              3
# 20          2 X15_salmonSharkBSAI_predAK                                      3
# 21          3 X04_Sablefish__JSOES                                            3
# 22          3 X09_chilipepper                                                 3
# 23          4 X09_DFA_HakeAge5Plus                                            3
# 24          4 X12_copepodCom_EGoA                                             3

#aic Penalty Sensitivity Analysis-------
#Because I don't believe those winners, I am running a sensitivity analysis on the aic penality
library(tidyverse)

# Define your baseline empirical penalty vector
base_penalty <- c("1" = 3.0, "2" = 2.5, "3" = 1.5, "4" = 1.0)
base_penalty <- c("1" = 2.0, "2" = 2, "3" = 2, "4" = 2)

multiplier_sweep <- seq(0, 2, by = 0.2)

sensitivity_results <- list()

for (mult in multiplier_sweep) {
  
  # Recalculate adjusted AIC using the current scaled penalty
  temp_adjusted <- dredge_results %>%
    mutate(
      Current_Penalty_Rate = base_penalty[as.character(Num_Params)] * mult,
      Adjusted_AIC = Raw_AIC + (Current_Penalty_Rate * Missing_Yrs)
    ) %>%
    group_by(Num_Params) %>%
    slice_min(order_by = Adjusted_AIC, n = 20) %>%
    ungroup() %>%
    separate_rows(Predictors, sep = " \\+ ")
  
  # Count occurrences across your expanded candidate list
  counts <- temp_adjusted %>%
    group_by(Predictors) %>%
    summarize(Times_Selected = n()) %>%
    filter(grepl(paste0("habCompInd|Sablefish|X05_DFA_abundSardine|",
                        "X15_DFA_sleeperSharkBSAI_predAK|X09_DFA_HakeAge5Plus|",
                        "X08_commonMurre_JSOES|X10_Harbor_seal_CR_2yrLead|",
                        "X10_Harbour_s_2yrLead_WS"), Predictors)) %>%
    mutate(Penalty_Multiplier = mult)
  
  sensitivity_results[[length(sensitivity_results) + 1]] <- counts
}

df_sensitivity <- bind_rows(sensitivity_results)

# --- CRITICAL FIX: Collapse any multi-tier duplicates before pivoting ---
pivot_sensitivity <- df_sensitivity %>%
  group_by(Predictors, Penalty_Multiplier) %>%
  summarize(Times_Selected = sum(Times_Selected), .groups = "drop") %>%
  
  # Now pivot_wider will execute cleanly with values_fill = 0
  pivot_wider(
    names_from = Penalty_Multiplier, 
    values_from = Times_Selected,
    values_fill = 0,
    names_prefix = "Mult_"
  )


#

print(base_penalty)
print(pivot_sensitivity)

#Original aic adj from 1-3

# Predictors                      Mult_0 Mult_0.2 Mult_0.4 Mult_0.6 Mult_0.8 Mult_1 Mult_1.2 Mult_1.4 Mult_1.6 Mult_1.8 Mult_2
# 1 X01_habCompInd                      35       36       37       37       37     37       39       39       34       31     30
# 2 X04_Sablefish__JSOES                15       15       16       16       16     16       16       16       16       18     19
# 3 X05_DFA_abundSardine                 1        1        1        1        1      1        1        1        1        1      1
# 4 X08_commonMurre_JSOES                6        6        6        5        3      2        1        0        0        0      0
# 5 X09_DFA_HakeAge5Plus                 8        8        7        9       10     12       13       16       18       21     22
# 6 X10_Harbor_seal_CR_2yrLead           4        4        4        3        3      2        1        1        1        1      1
# 7 X10_Harbour_s_2yrLead_WS             1        2        2        3        3      3        3        3        3        4      4
# 8 X15_DFA_sleeperSharkBSAI_predAK     11       11       11       11       10     10       10       11        9        9      9


#Std aic adj all=2
print(base_penalty)
1 2 3 4 
2 2 2 2 



> print(pivot_sensitivity)
# A tibble: 8 × 12
#   Predictors                      Mult_0 Mult_0.2 Mult_0.4 Mult_0.6 Mult_0.8 Mult_1 Mult_1.2 Mult_1.4 Mult_1.6 Mult_1.8 Mult_2
# 1 X01_habCompInd                      35       36       37       37       35     35       36       35       35       33     28
# 2 X04_Sablefish__JSOES                15       15       16       16       16     16       14       13       12        8      9
# 3 X05_DFA_abundSardine                 1        1        1        1        1      1        1        1        0        0      0
# 4 X08_commonMurre_JSOES                6        6        6        6        5      3        3        2        2        0      0
# 5 X09_DFA_HakeAge5Plus                 8        8        7        8       10     11       12       16       19       20     23
# 6 X10_Harbor_seal_CR_2yrLead           4        4        4        3        3      2        2        1        1        1      1
# 7 X10_Harbour_s_2yrLead_WS             1        2        2        2        3      3        3        4        5        6      8
# 8 X15_DFA_sleeperSharkBSAI_predAK     11       11       11       11        8      7        5        3        2        3      3


#3 & 4 param models heavily favored, ignore 1 & 2 param models---------
# 1. Isolate the complex models using your preferred flat penalty (e.g., Mult_0.4)
complex_models_only <- dredge_results %>%
  filter(Num_Params %in% c(3, 4)) %>%
  mutate(
    # Applying the stable flat penalty at your sensitivity sweet spot
    Adjusted_AIC = Raw_AIC + (2.0  * Missing_Yrs)
  #  Adjusted_AIC = Raw_AIC + (2.0 * 0.4 * Missing_Yrs)
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
           
  arrange(Adjusted_AIC)

head(complex_models_only)

#require all 4 guilds------
  complex_models_only %>% filter(Nguilds==4)
        # Formula                                                                             Predictors Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources Adjusted_AIC PreyNCC PredNCC PreyAK PredAK Nguilds
        # 1 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X12_copepodCom_EGoA + X15_DFA_slee… X01_habCo…          4    22.8           19           5 PreyNCC, Pre…         32.8       1       1      1      1       4
        # 2 X16_SAR ~ X01_habCompInd + X09_chilipepper + X13_gadid_WAI + X15_DFA_sleeperSharkB… X01_habCo…          4    23.5           19           5 PreyNCC, Pre…         33.5       1       1      1      1       4
        # 3 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X12_copepodCom_EGoA + X15_halibutB… X01_habCo…          4    23.9           19           5 PreyNCC, Pre…         33.9       1       1      1      1       4
        # 4 X16_SAR ~ X01_habCompInd + X09_chilipepper + X13_gadid_WAI + X15_halibutBiomassAge… X01_habCo…          4    24.5           19           5 PreyNCC, Pre…         34.5       1       1      1      1       4
        # 5 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X13_hexagram_EAI + X15_spinyDogfis… X01_habCo…          4    24.6           19           5 PreyNCC, Pre…         34.6       1       1      1      1       4
        # 6 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X12_copepodCom_EGoA + X15_spinyDog… X01_habCo…          4    24.9           19           5 PreyNCC, Pre…         34.9       1       1      1      1       4
        # 7 X16_SAR ~ X01_habCompInd + X09_chilipepper + X12_copepodCom_EGoA + X15_halibutBiom… X01_habCo…          4    25.0           19           5 PreyNCC, Pre…         35.0       1       1      1      1       4
        # 8 X16_SAR ~ X01_habCompInd + X08_DFA_DC_corm_3_WS + X13_ammod_WAI + X15_sablefishRec… X01_habCo…          4    35.1           24           0 PreyNCC, Pre…         35.1       1       1      1      1       4
        # 9 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X13_hexagram_EAI + X15_halibutBiom… X01_habCo…          4    25.1           19           5 PreyNCC, Pre…         35.1       1       1      1      1       4
        # 10 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X13_hexagram_WAI + X15_halibutBiom… X01_habCo…          4    25.2           19           5 PreyNCC, Pre…         35.2       1       1      1      1       4

complex_models_only %>% filter(Nguilds==3, Num_Params==3)
        #   Formula                                                                             Predictors Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources Adjusted_AIC PreyNCC PredNCC PreyAK PredAK Nguilds
        # 1 X16_SAR ~ X01_habCompInd + X08_DFA_DC_corm_3_WS + X13_ammod_WAI                     X01_habCo…          3    33.6           24           0 PreyNCC, Pre…         33.6       1       1      1      0       3
        # 2 X16_SAR ~ X01_habCompInd + X09_chilipepper + X15_halibutBiomassAge8plus_2yrLead_pr… X01_habCo…          3    24.1           19           5 PreyNCC, Pre…         34.1       1       1      0      1       3
        # 3 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X15_halibutBiomassAge8plus_2yrLead… X01_habCo…          3    24.5           19           5 PreyNCC, Pre…         34.5       1       1      0      1       3
        # 4 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X15_ArrowtoothFlounderBiomass_pred… X01_habCo…          3    25.4           19           5 PreyNCC, Pre…         35.4       1       1      0      1       3
        # 5 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X15_DFA_sleeperSharkBSAI_predAK     X01_habCo…          3    26.6           19           5 PreyNCC, Pre…         36.6       1       1      0      1       3
        # 6 X16_SAR ~ X01_habCompInd + X12_copepodCom_EGoA + X15_salmonSharkBSAI_predAK         X01_habCo…          3    26.8           19           5 PreyNCC, Pre…         36.8       1       0      1      1       3
        # 7 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X15_spinyDogfishGoA_predAK          X01_habCo…          3    27.0           19           5 PreyNCC, Pre…         37.0       1       1      0      1       3
        # 8 X16_SAR ~ X01_habCompInd + X09_chilipepper + X15_ArrowtoothFlounderBiomass_predAK   X01_habCo…          3    27.1           19           5 PreyNCC, Pre…         37.1       1       1      0      1       3
        # 9 X16_SAR ~ X01_habCompInd + X09_canaryRockfish + X10_DFA_ssl.est.wholerange_2yrLead  X01_habCo…          3    27.3           19           5 PreyNCC, Pre…         37.3       1       1      0      1       3
        # 10 X16_SAR ~ X01_habCompInd + X11_DFA_Harbour_p_WS + X15_spinyDogfishBSAI_predAK       X01_habCo…          3    27.4           19           5 PreyNCC, Pre…         37.4       1       1      0      1       3
#

#I think maybe I should exclude HCI
x<-complex_models_only %>% filter(!grepl("X01_habCompInd|X09_chilipepper", Predictors), Nguilds==4);x[,1:6]
#   Formula                                                                                                                               Predictors                Num_Params Raw_AIC Complete_Yrs Missing_Yrs
# 1 X16_SAR ~ X03_DFA_comMurreDietHerrSard_Yaquina_NCC + X11_DFA_Harbour_p_WS + X12_DFA_biomassEuphShelfSum + X15_spinyDogfishBSAI_predAK X03_DFA_comMurreDietHerr…          4    27.9           19           5
# 2 X16_SAR ~ X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X10_Harbour_s_2yrLead_WS + X13_hexagram_WAI                                   X04_Sablefish__JSOES + X…          4    38.8           24           0
# 3 X16_SAR ~ X05_anchovy_GAM + X09_canaryRockfish + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK                        X05_anchovy_GAM + X09_ca…          4    29.8           19           5
# 4 X16_SAR ~ X01_Red_phal_WS + X09_canaryRockfish + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK                        X01_Red_phal_WS + X09_ca…          4    30.4           19           5
# 5 X16_SAR ~ X01_NHLlogSum_win_05 + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X15_sablefishBiomass_predAK                                   X01_NHLlogSum_win_05 + X…          4    34.5           21           3
# 6 X16_SAR ~ X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X10_Harbour_s_2yrLead_WS + X13_hexagram_EAI                                   X04_Sablefish__JSOES + X…          4    40.6           24           0
# 7 X16_SAR ~ X04_PacificPompano_ + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X15_sablefishBiomass_predAK                                    X04_PacificPompano_ + X0…          4    40.7           24           0
# 8 X16_SAR ~ X04_PacificPompano_ + X08_commonMurre_JSOES + X13_gadid_WAI + X15_DFA_sleeperSharkBSAI_predAK                               X04_PacificPompano_ + X0…          4    31.0           19           5
# 9 X16_SAR ~ X01_NHLlogSum_win_25 + X09_DFA_HakeAge5Plus + X13_ammod_WAI + X15_sablefishBiomass_predAK                                   X01_NHLlogSum_win_25 + X…          4    35.1           21           3
# 10 X16_SAR ~ X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X10_Harbour_s_2yrLead_WS + X12_copepodCom_EGoA                                X04_Sablefish__JSOES + X…          4    41.2           24           0
#

# 2. See the ultimate Top 20 combinations left standing
top_complex_combinations <- complex_models_only %>%
  head(20) %>%
  select(Num_Params, Predictors, Raw_AIC, Missing_Yrs, Adjusted_AIC, Guild_Sources)%>%
  
  
  # 3. Sort by your new adjusted metric
  arrange(Adjusted_AIC)

print(top_complex_combinations)

#Adjusted_AIC = Raw_AIC + (2.0 * 0.4 * Missing_Yrs)
#favors short models only -- all 20 top models are missing 5 years
#so I'm boosting the penalty to 2. That pulled in 3 longer models
summary(top_complex_combinations$Missing_Yrs)

#a theme here is that 
