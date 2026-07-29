



#model="DAG1D_with_Mediators"
print(Sys.time())
library(tidyverse)
library(lavaan)

# --- 1. Define Variables & Metadata ---
response_var <- "X16_SAR"

# Pull the core mediator variables (X06 and X07)
x06_x07_meta <- var_lookup_NCC_AK %>% 
  filter(node_id %in% c("X06", "X07")) %>% 
  select(Lisaname, Role = SEMnode)

x06_x07_vars <- x06_x07_meta$Lisaname  # These MUST be included (at least one)

# Pull the other candidate environmental/trophic variables
other_candidate_meta <- var_lookup_NCC_AK %>% 
  filter(!node_id %in% c("X06", "X07", "X16")) %>% 
  select(Lisaname, Role = SEMnode)

other_candidates <- other_candidate_meta$Lisaname

# Master metadata map for role tracking
candidate_meta <- bind_rows(x06_x07_meta, other_candidate_meta)
all_candidates <- unique(c(x06_x07_vars, other_candidates))

# Calculate absolute correlation matrix across ALL active variables
cor_matrix <- abs(cor(guild.dfas1[, all_candidates], use = "pairwise.complete.obs"))

# --- 2. Generate Valid Combinations (1 to 4 predictors, requiring X06/X07) ---
model_list <- list()

for (m in 1:4) {
  # Determine how many X06/X07 variables (k) we want to force into this model size
  # k must be at least 1, and cannot exceed the total model size (m) or the number of available mediator vars
  max_k <- min(m, length(x06_x07_vars))
  
  for (k in 1:max_k) {
    # Step A: Get combinations of the X06/X07 mediators of size k
    x_combos <- combn(x06_x07_vars, k, simplify = FALSE)
    
    # Step B: Determine how many other environmental predictors are needed to reach size m
    other_size <- m - k
    
    if (other_size == 0) {
      # The model consists ONLY of X06 and/or X07 variables
      for (x_combo in x_combos) {
        if (k > 1) {
          # Collinearity check within the forced mediators
          sub_cor <- cor_matrix[x_combo, x_combo]
          if (max(sub_cor[lower.tri(sub_cor)]) > 0.5) next
        }
        model_list[[length(model_list) + 1]] <- x_combo
      }
    } else {
      # The model has k mediators and (m - k) other environmental variables
      other_combos <- combn(other_candidates, other_size, simplify = FALSE)
      
      for (x_combo in x_combos) {
        for (oth_combo in other_combos) {
          combo <- c(x_combo, oth_combo)
          
          # Collinearity check across the entire combined set
          sub_cor <- cor_matrix[combo, combo]
          if (max(sub_cor[lower.tri(sub_cor)]) > 0.5) next
          
          model_list[[length(model_list) + 1]] <- combo
        }
      }
    }
  }
}

# Quick validation check in console
cat("Generated", length(model_list), "valid model formulas to test.\n")

# # --- 3. Run the Custom Search and Record Metadata ---
# dredge_results <- map_dfr(model_list, function(predictors) {
#   
#   # Build formula text
#   formula_text <- paste(response_var, "~", paste(predictors, collapse = " + "))
#   
#   # Track missingness dynamically for the variables in this specific model
#   sub_data <- guild.dfas1[, c(response_var, predictors), drop = FALSE]
#   complete_yrs <- sum(complete.cases(sub_data))
#   missing_yrs <- 24 - complete_yrs  # Tracking the exact gap
#   
#   # Pull the metadata roles of the predictors in this model
#   roles_included <- candidate_meta %>% 
#     filter(Lisaname %in% predictors) %>% 
#     pull(Role) %>% 
#     paste(collapse = ", ")
#   
#   # Fit using FIML to maximize data usage
#   fit <- tryCatch({
#     sem(formula_text, data = guild.dfas1, missing = "ML")
#   }, error = function(e) NULL)
#   
#   if (is.null(fit)) return(tibble())
#   
#   measures <- fitMeasures(fit)
#   
#   tibble(
#     Formula = formula_text,
#     Predictors = paste(predictors, collapse = " + "),
#     Num_Params = length(predictors),
#     Raw_AIC = measures[["aic"]],
#     Complete_Yrs = complete_yrs,
#     Missing_Yrs = missing_yrs,
#     Guild_Sources = roles_included
#   )
# })

# --- 3. Run the Custom Search with a 5-Minute Progress Reporter ---

total_models <- length(model_list)
model_counter <- 0
last_report_time <- Sys.time()

cat("Starting model fitting at:", format(last_report_time, "%H:%M:%S"), "\n")
cat("Total combinations to fit:", total_models, "\n\n")

dredge_results_X06X07 <- map_dfr(model_list, function(predictors) {
  
  # Increment global counter
  model_counter <<- model_counter + 1
  
  # Check if 5 minutes have passed since the last report
  current_time <- Sys.time()
  elapsed_mins <- as.numeric(difftime(current_time, last_report_time, units = "mins"))
  
  # Print progress on the very first model, the last model, or every 5 minutes
  if (model_counter == 1 || model_counter == total_models || elapsed_mins >= 5) {
    pct_complete <- (model_counter / total_models) * 100
    cat(sprintf("[%s] %d/%d models completed (%.1f%%)\n", 
                format(current_time, "%Y-%m-%d %H:%M:%S"), 
                model_counter, 
                total_models, 
                pct_complete))
    # Reset the clock tracker
    last_report_time <<- current_time
  }
  
  # --- Standard Model Fitting ---
  formula_text <- paste(response_var, "~", paste(predictors, collapse = " + "))
  
  # Track missingness dynamically
  sub_data <- guild.dfas1[, c(response_var, predictors), drop = FALSE]
  complete_yrs <- sum(complete.cases(sub_data))
  missing_yrs <- 24 - complete_yrs  
  
  # Pull the metadata roles
  roles_included <- candidate_meta %>% 
    filter(Lisaname %in% predictors) %>% 
    pull(Role) %>% 
    paste(collapse = ", ")
  
  # Fit using FIML
  fit <- tryCatch({
    sem(formula_text, data = guild.dfas1, missing = "ML")
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(tibble())
  
  measures <- fitMeasures(fit)
  
  tibble(
    Formula = formula_text,
    Predictors = paste(predictors, collapse = " + "),
    Num_Params = length(predictors),
    Raw_AIC = measures[["aic"]],
    Complete_Yrs = complete_yrs,
    Missing_Yrs = missing_yrs,
    Guild_Sources = roles_included
  )
})
print(Sys.time())

#combine with previous results?--------
        topdredge<-read.csv(file = "LisaXP/outputs_4/Top_100_dredgemodels_r2.csv", row.names = NULL)
        
        top5<-topdredge[1:5,] %>% select(Formula,Complete_Yrs,Num_Params,R2,All_Coefficients_Sig,Raw_AIC,Adjusted_AIC,Adj_AIC2,P_Values,Coefficients)
        top5[,1]
        top5
        #4/5 are all signif, including the top model
        #the one that wasn't all signif had only 21 yrs (the rest were 24 yrs)
        
        sig<-detailed_results %>%
          filter(All_Coefficients_Sig=="Yes")


# --- 4. Save and inspect results ---
        dredge_results_X06X07<-dredge_results_X06X07 %>%
          mutate(adjaic2=Raw_AIC + 2*Missing_Yrs) %>%
          arrange(adjaic2)
        head(dredge_results_X06X07)
write.csv(as.data.frame(dredge_results_X06X07), 
          file = "LisaXP/outputs_4/NCC_AK_dredge_results_with_X06X07.csv", 
          row.names = FALSE)

#remove hci
        dredge_results_X06X07_noHCI <- dredge_results_X06X07 %>% filter(!grepl("X01_habCompInd",Predictors))
        dredge_results_X06X07_noHCI %>% head()

        # Formula                                                                                                                         Predictors Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources adjaic2
        # 1 X16_SAR ~ X06_Lmu_IntSprJunH + X05_anchovy_GAM + X12_DFA_biomassEuphShelfSum + X12_copepodCom_EGoA                              X06_Lmu_I…          4    31.6           24           0 Growth, Prey…    31.6
        # 2 X16_SAR ~ X06_Lmu_IntSprJunH + X04_Sablefish__JSOES + X08_DFA_DC_corm_3_WS + X09_DFA_HakeAge5Plus                               X06_Lmu_I…          4    32.7           24           0 Growth, Prey…    32.7
        # 3 X16_SAR ~ X06_Lmu_IntSprJunH + X09_DFA_HakeAge5Plus + X12_DFA_biomassEuphShelfSum + X12_copepodCom_EGoA                         X06_Lmu_I…          4    33.7           24           0 Growth, Pred…    33.7
        # 4 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK + X15_spinyDogfishBSAI_pred… X07_DFA_c…          4    23.8           19           5 Abundance, P…    33.8
        # 5 X16_SAR ~ X06_Lmu_IntSprJunH + X04_Sablefish__JSOES + X09_DFA_ChinAbundSnakeFall + X09_DFA_HakeAge5Plus                         X06_Lmu_I…          4    34.2           24           0 Growth, Prey…    34.2
        # 6 X16_SAR ~ X06_DFA_IGF_mu + X04_Sablefish__JSOES + X09_DFA_HakeAge5Plus + X09_canaryRockfish                                     X06_DFA_I…          4    24.2           19           5 Growth, Prey…    34.2
        
        #ineresting that growth did better than cpue
        
        dredge_results_X06X07_cpue <- dredge_results_X06X07 %>% filter(grepl("X07_DFA_cpue_IntSprJunHW",Predictors))
        dredge_results_X06X07_cpue %>% head()
        
        # Formula                                                                                                                         Predictors Num_Params Raw_AIC Complete_Yrs Missing_Yrs Guild_Sources adjaic2
        # 1 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK + X15_spinyDogfishBSAI_pred… X07_DFA_c…          4    23.8           19           5 Abundance, P…    33.8
        # 2 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X09_chilipepper + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK            X07_DFA_c…          4    24.5           19           5 Abundance, P…    34.5
        # 3 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_ArrowtoothFlounderBiomass_predAK + X15_PacificCodBiomas… X07_DFA_c…          4    36.9           24           0 Abundance, P…    36.9
        # 4 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X04_Sablefish__JSOES + X09_canaryRockfish + X10_Harbor_seal_CR_2yrLead                     X07_DFA_c…          4    28.0           19           5 Abundance, P…    38.0
        # 5 X16_SAR ~ X06_Lmu_IntSprJunH + X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_spinyDogfishBSAI_predAK             X06_Lmu_I…          4    28.1           19           5 Growth, Abun…    38.1
        # 6 X16_SAR ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X13_ammod_WAI + X13_hexagram_EAI                             X07_DFA_c…          4    38.5           24           0 Abundance, P…    38.5
        
        
                