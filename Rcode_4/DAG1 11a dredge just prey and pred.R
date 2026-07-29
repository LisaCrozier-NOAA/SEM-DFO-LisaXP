
#model="DAG1D"

library(tidyverse)
library(lavaan)

# --- 1. Define Variables & Metadata ---
response_var <- "X16_SAR"

# Pull candidate names and keep their metadata tags
candidate_meta <- var_lookup_NCC_AK %>% 
  filter(!node_id %in% c("X06", "X07", "X16")) %>% 
  select(Lisaname, Role = SEMnode) # e.g., PreyNCC, PredAK, etc.

candidates <- candidate_meta$Lisaname

# Calculate absolute correlation matrix
cor_matrix <- abs(cor(guild.dfas1[, candidates], use = "pairwise.complete.obs"))

# --- 2. Generate Valid All-Subsets Combinations (up to 4 predictors) ---
model_list <- list()

for (m in 1:4) {
  combos <- combn(candidates, m, simplify = FALSE)
  
  for (combo in combos) {
    # Check for collinearity within this specific combination
    if (m > 1) {
      sub_cor <- cor_matrix[combo, combo]
      if (max(sub_cor[lower.tri(sub_cor)]) > 0.5) next # Skip if pairs > 0.5
    }
    
    model_list[[length(model_list) + 1]] <- combo
  }
}

# --- 3. Run the Custom Search and Record Metadata ---
dredge_results <- map_dfr(model_list, function(predictors) {
  
  # Build formula text
  formula_text <- paste(response_var, "~", paste(predictors, collapse = " + "))
  
  # Track missingness dynamically for the variables in this specific model
  sub_data <- guild.dfas1[, c(response_var, predictors), drop = FALSE]
  complete_yrs <- sum(complete.cases(sub_data))
  missing_yrs <- 24 - complete_yrs  # Tracking the exact gap
  
  # Pull the metadata roles of the predictors in this model
  roles_included <- candidate_meta %>% 
    filter(Lisaname %in% predictors) %>% 
    pull(Role) %>% 
    paste(collapse = ", ")
  
  # Fit using FIML to maximize data usage
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



head(dredge_results)
write.csv(as.data.frame(dredge_results),file="LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names=FALSE)
dredge_results<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names = NULL)


