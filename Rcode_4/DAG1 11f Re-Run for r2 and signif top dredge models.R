#collect r2 and sign of coefficients from top 100 models--------

library(tidyverse)
library(lavaan)

# 1. Ensure your response variable is defined
response_var <- "X16_SAR"

# 2. Apply your preferred penalty (e.g., flat +2 at the 0.4 multiplier) to identify the true top 100
top_100_models <- dredge_results %>%
  mutate(Adjusted_AIC = Raw_AIC + (1.0  * Missing_Yrs)) %>%
  mutate(Adj_AIC2 = Raw_AIC + (2.0  * Missing_Yrs)) %>%
  arrange(Adj_AIC2)%>% 
  filter(grepl("X09_DFA_HakeAge5Plus",Predictors)) %>%
  filter(!grepl("X09_chilipepper|X09_canaryRockfish",Predictors))  %>%
  head(100)

# 3. Refit ONLY these 100 models and extract the detailed statistics
detailed_results <- top_100_models %>%
  split(1:nrow(.)) %>%  # Split row-by-row to map cleanly
  map_dfr(function(row) {
    
    formula_text <- row$Formula
    predictors <- unlist(strsplit(row$Predictors, " \\+ "))
    
    # Fit the model
    fit <- tryCatch({
      sem(formula_text, data = guild.dfas1, missing = "ML")
    }, error = function(e) NULL)
    
    # If a model fails, return the original row with blank stats
    if (is.null(fit)) {
      return(bind_cols(row, tibble(
        R2 = NA, Coefficients = NA, Std_Errors = NA, P_Values = NA, All_Coefficients_Sig = NA
      )))
    }
    
    # Extract R-squared for X16_SAR
    r2_val <- inspect(fit, "rsquare")[[response_var]]
    
    # Extract coefficients table
    est_table <- parameterEstimates(fit) %>%
      filter(lhs == response_var & op == "~" & rhs %in% predictors)
    
    # Ensure the statistics match the exact order of your Predictors string
    est_table <- est_table %>%
      slice(match(predictors, rhs))
    
    coefs_str  <- paste(round(est_table$est, 4), collapse = ", ")
    se_str     <- paste(round(est_table$se, 4), collapse = ", ")
    p_vals_str <- paste(round(est_table$pvalue, 4), collapse = ", ")
    
    # Check if all terms are significant (p < 0.05)
    all_significant <- if_else(all(est_table$pvalue < 0.05, na.rm = TRUE), "Yes", "No")
    
    # Bind the new metrics right onto the existing metadata columns
    bind_cols(row, tibble(
      R2 = r2_val,
      Coefficients = coefs_str,
      Std_Errors = se_str,
      P_Values = p_vals_str,
      All_Coefficients_Sig = all_significant
    ))
  })

# 4. Save the top 100 detailed summary
write.csv(detailed_results, file = "LisaXP/outputs_4/Top_100_dredgemodels_r2.csv", row.names = FALSE)

sig<-detailed_results %>%
  filter(All_Coefficients_Sig=="Yes")

nrow(sig) # 35
hist(sig$R2) #0.72-0.84
head(sig)
sig[,1]


