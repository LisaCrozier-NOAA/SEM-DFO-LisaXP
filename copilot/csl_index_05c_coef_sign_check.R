#coefficient sign directions against theoretical expectations, adds directional pass/fail columns to the all_models table, builds a category-level summary, and updates the var_importance table with sign-consistency metrics.Directional Expectations Rules AppliedClimate (x21_): Negative ($<0$)Prey (x12_, x13_): Positive ($>0$)Competitor (x14_): Negative ($<0$)Predator (x10_, x15_) & Risk Indices (ishark, issl, icsl): Negative ($<0$
#original files saved with all_models
#all_models_sub <- all_models 
  
# ==============================================================================
# Script Extension: Directional Expectation Verification & Variable Importance
# ==============================================================================

# Helper function to classify variable group and expected sign
get_expected_sign <- function(var_name) {
  if (is.na(var_name) || var_name == "") return(NA_character_)
  
  case_when(
    str_detect(var_name, "ishark|issl|icsl") ~ "Negative",
    str_detect(var_name, "^x10_|^x15_")       ~ "Negative",
    str_detect(var_name, "^x21_")              ~ "Negative",
    str_detect(var_name, "^x12_|^x13_")        ~ "Positive",
    str_detect(var_name, "^x14_")              ~ "Negative",
    TRUE                                       ~ "Unknown"
  )
}

# Helper to check if estimated coefficient meets expectation
check_sign <- function(b_est, exp_sign) {
  if (is.na(b_est) || is.na(exp_sign) || exp_sign == "Unknown") return(NA)
  if (exp_sign == "Positive") return(b_est > 0)
  if (exp_sign == "Negative") return(b_est < 0)
  return(NA)
}

# -----------------------------------------------------------------------------
# 1. Annotate `all_models` with Expected Signs and Pass/Fail Booleans
# -----------------------------------------------------------------------------

all_models_sub <- all_models_sub %>%
  mutate(
    # Expected signs for terms 1, 2, and 3
    exp_sign_1 = sapply(term1, get_expected_sign),
    exp_sign_2 = sapply(term2, get_expected_sign),
    exp_sign_3 = sapply(term3, get_expected_sign),
    
    # Check if estimate matches expectation (TRUE/FALSE)
    pass_sign_1 = mapply(check_sign, b1, exp_sign_1),
    pass_sign_2 = mapply(check_sign, b2, exp_sign_2),
    pass_sign_3 = mapply(check_sign, b3, exp_sign_3),
    
    # Model-level count and proportion of signs meeting expectations
    n_valid_terms = (!is.na(pass_sign_1)) + (!is.na(pass_sign_2)) + (!is.na(pass_sign_3)),
    n_passed_terms = replace_na(pass_sign_1, FALSE) + replace_na(pass_sign_2, FALSE) + replace_na(pass_sign_3, FALSE),
    prop_model_signs_met = n_passed_terms / n_valid_terms,
    all_signs_met = (n_passed_terms == n_valid_terms)
  )

# -----------------------------------------------------------------------------
# 2. Reshape into Term-Level Dataset for Summary Statistics
# -----------------------------------------------------------------------------

terms_evaluated <- bind_rows(
  all_models_sub %>% select(model_label, aic_weight, variable = term1, estimate = b1, exp_sign = exp_sign_1, sign_met = pass_sign_1),
  all_models_sub %>% select(model_label, aic_weight, variable = term2, estimate = b2, exp_sign = exp_sign_2, sign_met = pass_sign_2),
  all_models_sub %>% select(model_label, aic_weight, variable = term3, estimate = b3, exp_sign = exp_sign_3, sign_met = pass_sign_3)
) %>%
  filter(!is.na(variable)) %>%
  mutate(
    var_group = case_when(
      str_detect(variable, "ishark|issl|icsl") ~ "Index",
      str_detect(variable, "^x10_|^x15_")       ~ "Predator",
      str_detect(variable, "^x21_")              ~ "Climate",
      str_detect(variable, "^x12_|^x13_")        ~ "Prey",
      str_detect(variable, "^x14_")              ~ "Competitor",
      TRUE                                       ~ "Other"
    )
  )

# -----------------------------------------------------------------------------
# 3. Create Summary Table by Functional Category
# -----------------------------------------------------------------------------

category_expectation_summary <- terms_evaluated %>%
  group_by(var_group, exp_sign) %>%
  summarise(
    total_evaluations     = n(),
    n_expectations_met    = sum(sign_met, na.rm = TRUE),
    raw_prop_met          = mean(sign_met, na.rm = TRUE),
    aic_weighted_prop_met = sum(aic_weight[sign_met == TRUE], na.rm = TRUE) / sum(aic_weight, na.rm = TRUE),
    .groups               = "drop"
  ) %>%
  arrange(desc(aic_weighted_prop_met))

# -----------------------------------------------------------------------------
# 4. Update Variable Importance Table with Proportion of Time Expectations Met
# -----------------------------------------------------------------------------

var_sign_stats <- terms_evaluated %>%
  group_by(variable) %>%
  summarise(
    expected_sign         = first(exp_sign),
    times_in_models       = n(),
    n_sign_met            = sum(sign_met, na.rm = TRUE),
    prop_sign_met_unwtd   = mean(sign_met, na.rm = TRUE),
    prop_sign_met_aic_wtd = sum(aic_weight[sign_met == TRUE], na.rm = TRUE) / sum(aic_weight, na.rm = TRUE),
    .groups               = "drop"
  )

var_importance <- var_importance %>%
  left_join(var_sign_stats, by = "variable") %>%
  relocate(expected_sign, prop_sign_met_aic_wtd, prop_sign_met_unwtd, .after = importance)%>%
  arrange(desc(importance))

# -----------------------------------------------------------------------------
# 5. Export Updated Tables
# -----------------------------------------------------------------------------

# write.csv(all_models, file.path(out_dir, "sem_global_models_with_sign_checks.csv"), row.names = FALSE)
# write.csv(category_expectation_summary, file.path(out_dir, "sem_category_expectation_summary.csv"), row.names = FALSE)
# write.csv(var_importance, file.path(out_dir, "sem_global_variable_importance_with_sign_props.csv"), row.names = FALSE)

# Print Summary Tables to Console
cat("\n========================================================\n")
cat("   SUMMARY: EXPECTATION SUCCESS RATE BY CATEGORY        \n")
cat("========================================================\n")
print(category_expectation_summary)

cat("\n========================================================\n")
cat("   TOP 15 VARIABLE IMPORTANCE WITH SIGN SUCCESS RATES   \n")
cat("========================================================\n")
print(
  var_importance %>% 
    select(variable, var_group, expected_sign, importance, prop_sign_met_aic_wtd, prop_sign_met_unwtd) %>%
 #   filter(importance>0.2) %>%
    slice_head(n = 15) %>%
    arrange(desc(importance))
)
