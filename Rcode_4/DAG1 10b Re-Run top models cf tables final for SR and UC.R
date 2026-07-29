library(tidyverse)
library(lavaan)
library(glue)

# --- 1. Define base formulas with {SAR_VAR} placeholder ---
dag_defs <- list(
  DAG1C_long          = 'X07_DFA_cpue_IntSprJunHW ~ X05_herring_GAM + X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_Lmu_IntSprJunH + X05_herring_GAM + X09_DFA_HakeAge5Plus + X12_copepodCom_EGoA + X15_salmonSharkGoA_predAK',
  DAG1C_long_reduced  = 'X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X06_Lmu_IntSprJunH + X05_herring_GAM + X09_DFA_HakeAge5Plus + X12_copepodCom_EGoA + X15_salmonSharkGoA_predAK',
  DAG1C_short         = 'X07_DFA_cpue_IntSprJunHW ~ X01_habCompInd + X08_commonMurre_JSOES\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_StomFull_May + X01_habCompInd + X08_commonMurre_JSOES + X14_pinkSalmon + X10_Harbor_seal_CR_2yrLead',
  DAG1C_short_reduced = 'X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_StomFull_May + X01_habCompInd + X08_commonMurre_JSOES + X14_pinkSalmon + X10_Harbor_seal_CR_2yrLead',
  DAG1A_long          = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X06_DFA_IGF_mu + X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK',
  DAG1A_long_reduced  = 'X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK',
  DAG1A_short         = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES + X06_DFA_IGF_mu\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_copepodBiomass_WGoA + X10_Harbor_seal_CR_2yrLead',
  DAG1A_short_reduced = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_copepodBiomass_WGoA + X10_Harbor_seal_CR_2yrLead',
  DAG1B_long          = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_DFA_IGF_mu + X12_DFA_biomassEuphShelfSum + X15_PacificCodBiomass_predAK',
  DAG1B_long_reduced  = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_DFA_IGF_mu + X12_DFA_biomassEuphShelfSum + X15_PacificCodBiomass_predAK',
  DAG1B_short         = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES\n{SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X06_DFA_IGF_mu + X13_pollockBiomassGoAage3plus_predAK + X15_DFA_sleeperSharkBSAI_predAK',
  DAG1B_short_reduced = 'X06_DFA_IGF_mu ~ X05_DFA_abundSardine\nX07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES\n{SAR_VAR} ~ X06_DFA_IGF_mu + X13_pollockBiomassGoAage3plus_predAK + X15_DFA_sleeperSharkBSAI_predAK'
)

target_responses <- c("X16_SAR", "sarSR", "sarUC")

# Create an empty environment to store the raw lavaan objects cleanly
fitted_models <- list()
summary_rows <- list()

# --- 2. Loop, Fit, and Store ---
for (dag in names(dag_defs)) {
  fitted_models[[dag]] <- list()
  
  for (resp in target_responses) {
    model_syntax <- glue(dag_defs[[dag]], SAR_VAR = resp)
    
    fit <- tryCatch({
      sem(model_syntax, data = guild.dfas1, std.lv = TRUE, missing = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      # 1. Save the actual model object to our nested list
      fitted_models[[dag]][[resp]] <- fit
      
      # 2. Collect metadata for the quick summary table
      measures <- fitMeasures(fit)
      pe <- standardizedSolution(fit)
      
      coefs_to_sar <- pe %>% 
        filter(op == "~", lhs == resp) %>% 
        mutate(coef_label = paste0("beta_", rhs)) %>% 
        select(coef_label, est.std) %>% 
        pivot_wider(names_from = coef_label, values_from = est.std)
      
      row_meta <- tibble(DAG = dag, Response = resp, N_Years = nobs(fit), AIC = round(measures[["aic"]], 1), AGFI = round(measures[["agfi"]], 3))
      summary_rows[[length(summary_rows) + 1]] <- bind_cols(row_meta, coefs_to_sar)
    }
  }
}

# Combine the table summary rows
master_results <- bind_rows(summary_rows)
print(master_results, width = Inf)

# --- 3. Save the actual objects to disk ---
save(fitted_models, file = "LisaXP/outputs_4/DAG1.SRvUC_responses.rdata")


# Load your models back in
load("LisaXP/outputs_4/DAG1.swapped_responses.rdata")

# Inspect the summary of DAG1C_long using the raw sarSR column
summary(fitted_models$DAG1C_long$sarSR, standardized = TRUE, fit.measures = TRUE)

# Compare parameter estimates side-by-side for a specific DAG structure
parameterEstimates(fitted_models$DAG1C_long$X16_SAR) # Original DFA
parameterEstimates(fitted_models$DAG1C_long$sarSR)   # Complete raw
parameterEstimates(fitted_models$DAG1C_long$sarUC)   # Missing data raw


library(tidyverse)
library(lavaan)
library(glue)

# (Assuming 'dag_defs', 'target_responses', and 'guild.dfas1' are loaded from previous turns)

summary_rows <- list()

for (dag in names(dag_defs)) {
  for (resp in target_responses) {
    
    # 1. Inject the specific response variable into the syntax
    model_syntax <- glue(dag_defs[[dag]], SAR_VAR = resp)
    
    # 2. Extract every variable explicitly mentioned in this specific model structure
    # This prevents counting missing data in columns you aren't even using!
    model_vars <- lavaanify(model_syntax) %>% 
      filter(op == "~") %>% 
      select(lhs, rhs) %>% 
      pivot_longer(cols = everything()) %>% 
      distinct() %>% 
      pull(value)
    
    # Extract the subset of data used for this specific model configuration
    sub_data <- guild.dfas1[, model_vars, drop = FALSE]
    
    # Calculate exact missingness metrics for this specific DAG + Response combo
    complete_years <- sum(complete.cases(sub_data))
    total_missing_cells <- sum(is.na(sub_data))
    
    # 3. Fit the model via FIML
    fit <- tryCatch({
      sem(model_syntax, data = guild.dfas1, std.lv = TRUE, missing = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      measures <- fitMeasures(fit)
      pe <- standardizedSolution(fit)
      
      # Extract standardized paths leading directly to the response variable
      coefs_to_sar <- pe %>% 
        filter(op == "~", lhs == resp) %>% 
        mutate(coef_label = paste0("beta_", rhs)) %>% 
        select(coef_label, est.std) %>% 
        pivot_wider(names_from = coef_label, values_from = est.std)
      
      # Build the row metadata
      row_meta <- tibble(
        DAG = dag,
        Response = resp,
        FIML_N = nobs(fit),                # Total rows entering FIML engine (24)
        Complete_Yrs = complete_years,     # True years with zero missing data in this model
        Missing_Cells = total_missing_cells, # Total empty cells across all model variables
        AIC = round(measures[["aic"]], 1),
        AGFI = round(measures[["agfi"]], 3)
      )
      
      summary_rows[[length(summary_rows) + 1]] <- bind_cols(row_meta, coefs_to_sar)
    }
  }
}

# --- 4. Format and calculate localized dAIC ---
final_comparison_table <- bind_rows(summary_rows) %>%
  group_by(Response) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  arrange(Response, dAIC) %>%
  mutate(dAIC = round(dAIC, 1)) %>%
  select(Response, DAG, dAIC, AIC, AGFI, Complete_Yrs, Missing_Cells, FIML_N, everything()) %>%
  ungroup()

# View the clean comparison matrix
print(final_comparison_table, width = Inf, n = Inf)


#compare dAIC only===========
# Assuming 'final_comparison_table' has already been generated by the previous script

pivot_daic_table <- final_comparison_table %>%
  # 1. Keep only the core identifying information we want to pivot
  select(DAG, Response, dAIC) %>%
  
  # 2. Reshape the table from long format to wide format
  pivot_wider(
    names_from = Response, 
    values_from = dAIC,
    names_prefix = "dAIC_" # Adds a clear prefix to your column names
  ) %>%
  
  # 3. Optional: Sort the rows by the performance on your original baseline model (X16_SAR)
  arrange(dAIC_X16_SAR)

# Print the clean pivot table
print(pivot_daic_table)

#results:
# ok, so the good news is that the rank order of the models is almost perfectly identical, 
# with a single swap for UC, 
# where DAG1A_long_reduced was a little better than DAG1C_long_reduced although it was much worse for SR. 



#     DAG                 dAIC_X16_SAR dAIC_sarSR dAIC_sarUC
# 1 DAG1C_short_reduced          0          0          0  
# 2 DAG1C_short                  1.4        1.3        1.3
# 3 DAG1C_long_reduced          55.7       34.8       36  
# 4 DAG1C_long                  57.7       36.7       37.9
# 5 DAG1A_long_reduced          79.1       46.2       33.6
# 6 DAG1B_short_reduced         85.6       55         58  
# 7 DAG1A_short                 86.7       58.9       60.4
# 8 DAG1B_short                 87.3       56.4       60  
# 9 DAG1A_short_reduced         87.4       59         61.1
# 10 DAG1B_long                 129.       102.        91.3
# 11 DAG1B_long_reduced         129.       102.        91.3
# 12 DAG1A_long                 138.       106.        93  



#check usefulness of cpue============
# Assuming 'final_comparison_table' has already been generated by the previous script
names(final_comparison_table)
pivot_cpue_table <- final_comparison_table %>%
  # 1. Keep only the core identifying information we want to pivot
  select(DAG, Response, beta_X07_DFA_cpue_IntSprJunHW) %>%
  
  # 2. Reshape the table from long format to wide format
  pivot_wider(
    names_from = Response, 
    values_from = beta_X07_DFA_cpue_IntSprJunHW,
    names_prefix = "cpue_" # Adds a clear prefix to your column names
  ) %>%
  
  # 3. Optional: Sort the rows by the performance on your original baseline model (X16_SAR)
  arrange(cpue_X16_SAR)

# Print the clean pivot table
print(pivot_cpue_table)


#     DAG                 cpue_X16_SAR cpue_sarSR cpue_sarUC
# 1 DAG1C_long                 0.117      0.321     0.0899
# 2 DAG1B_short                0.199      0.371     0.166 
# 3 DAG1A_short                0.318      0.519     0.192 
# 4 DAG1A_short_reduced        0.320      0.513     0.202 
# 5 DAG1B_long                 0.320      0.457     0.212 
# 6 DAG1B_long_reduced         0.320      0.457     0.212 
# 7 DAG1A_long_reduced         0.472      0.543     0.269 
# 8 DAG1A_long                 0.473      0.544     0.270 
# 9 DAG1C_short                0.680      0.754     0.734 
# 10 DAG1C_short_reduced       0.709      0.781     0.763 
# 11 DAG1C_long_reduced        NA         NA        NA     
# 12 DAG1B_short_reduced       NA         NA        NA     