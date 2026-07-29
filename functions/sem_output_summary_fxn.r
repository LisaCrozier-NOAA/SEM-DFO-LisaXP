sem_output_summary_fxn <- function(model_fit) {
  
  # 1. Extract R-Squared for available endogenous variables
  r2_all <- lavInspect(model_fit, "rsquare")
  target_vars <- c("X16_SAR", "X06_DFA_IGF_mu  ", "X07_DFA_cpue_IntSprJunHW")
#  target_vars <- c("X15.SAR", "X06.Cond1NCC_DFA1_IGF_mu", "X07.Cond2NCC_DFA1_cpue_IntSprJunHW")
  existing_r2 <- r2_all[names(r2_all) %in% target_vars]
  
  # 2. Extract, Round, and Rank Regressions
  regression_ranking <- parameterEstimates(model_fit) %>%
    filter(op == "~") %>%
    select(lhs, op, rhs, est, se, z, pvalue) %>%
    # Round the numeric columns to 3 decimal places
    mutate(across(c(est, se, z, pvalue), ~ round(.x, 3))) %>%
    arrange(desc(pvalue))
  
  # 3. Extract Fit Measures
  all_metrics <- fitMeasures(model_fit)
  metrics_to_show <- c("aic", "pvalue", "cfi", "agfi", "rmsea")
  metrics <- round(all_metrics[names(all_metrics) %in% metrics_to_show], 3)
  
  # --- Printed Output ---
  cat("\n==============================================\n")
  cat("           MODEL SUMMARY REPORT               ")
  cat("\n==============================================\n")
  
  cat("\n--- Model Fit Metrics ---\n")
  print(metrics)
  
  cat("\n--- R-Squared (Endogenous Variables) ---\n")
  if(length(existing_r2) > 0) {
    print(round(existing_r2, 3))
  } else {
    cat("No target R2 variables found in this model.\n")
  }

    cat("\n--- Regression Path Ranking (By P-Value) ---\n")
    
    
  # Printing as data frame to ensure columns align in console
  print(as.data.frame(regression_ranking), row.names = FALSE)
  
  cat("\n==============================================\n")
  
  return(invisible(list(
    r2 = existing_r2,
    paths = regression_ranking,
    fit = metrics
  )))
}
