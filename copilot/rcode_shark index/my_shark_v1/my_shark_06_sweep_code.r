suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})


#allow T into the model. First put enso into guild.dfas1 dataset

# -----------------------------------------------------------------------------
# 1. Helper Function: Standardized Fitting with Strict Year Tracking
# -----------------------------------------------------------------------------
fit_sem_two_topdown <- function(pred1_col, pred2_col, df_data) {
  
  # Clean subset for the target model
  df_model <- df_data %>%
    filter(year >= 1998, year <= 2021) %>%
    mutate(
      Pred1 = .data[[pred1_col]],
      Pred2 = .data[[pred2_col]]
    ) %>%
    filter(
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(x16_sar),
      !is.na(Pred1),
      !is.na(Pred2)
    )
  
  n_years <- nrow(df_model)
  
  # Skip if sample size is truncated due to NAs
  if (n_years < 12) return(NULL)
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + Pred1 + Pred2
  '
  
  fit <- tryCatch({
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "listwise", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea", "npar"))
  pe <- parameterEstimates(fit)
  
  get_est <- function(lhs_v, rhs_v, field) {
    val <- pe %>% filter(lhs == lhs_v, op == "~", rhs == rhs_v) %>% pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  tibble(
    term1_col     = pred1_col,
    term2_col     = pred2_col,
    n_years       = n_years,
    aic           = fm[["aic"]],
    bic           = fm[["bic"]],
    cfi           = fm[["cfi"]],
    rmsea         = fm[["rmsea"]],
    
    b_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    
    b_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_term1_sar   = get_est("x16_sar", "Pred1", "est"),
    p_term1_sar   = get_est("x16_sar", "Pred1", "pvalue"),
    
    b_term2_sar   = get_est("x16_sar", "Pred2", "est"),
    p_term2_sar   = get_est("x16_sar", "Pred2", "pvalue")
  )
}

# -----------------------------------------------------------------------------
# 2. Run Iterative Model Sweeps
# -----------------------------------------------------------------------------
# Sweep 1: 2 Predators
pred_pairs <- combn(ak_pred_cands, 2, simplify = FALSE)
res_pred_pairs <- purrr::map_dfr(pred_pairs, ~ fit_sem_two_topdown(.x[1], .x[2], guild.dfas1)) %>%
  mutate(model_type = "2 Predators")

# Sweep 2: 1 Predator + 1 Prey
prey_pred_grid <- expand.grid(pred_col = ak_pred_cands, prey_col = ak_prey_cands, stringsAsFactors = FALSE)
res_pred_prey <- purrr::map2_dfr(
  prey_pred_grid$pred_col, 
  prey_pred_grid$prey_col, 
  function(p1, p2) fit_sem_two_topdown(p1, p2, guild.dfas1)
) %>%
  mutate(model_type = "1 Predator + 1 Prey")
# -----------------------------------------------------------------------------
# 3. Filter for Complete Cases (Max N) & Rank Fairly
# -----------------------------------------------------------------------------
all_models <- bind_rows(res_pred_pairs, res_pred_prey)

# Determine maximum complete years in the dataset (usually 24 for 1998-2021)
max_n <- max(all_models$n_years, na.rm = TRUE)

cat("Maximum years available across models: N =", max_n, "\n")

# Fit baseline model on the exact same full N dataset
df_base <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(x16_sar))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw
', data = df_base, std.lv = TRUE, missing = "listwise", warn = FALSE)

base_aic <- fitMeasures(fit_base, "aic")

# Rank ONLY models that have the full sample size (n_years == max_n)
fair_ranked_models <- all_models %>%
  filter(n_years == max_n) %>% # Eliminates bias from missing years!
  mutate(
    delta_aic_vs_base = aic - base_aic,
    contains_new_shark = (term1_col == "x15_shark_enso_roll2" | term2_col == "x15_shark_enso_roll2")
  ) %>%
  arrange(aic)

# Print Top 15 Fairly-Compared Models
cat("\n=== TOP 15 MODELS EVALUATED ON COMPLETE DATASET (N =", max_n, "YEARS) ===\n")
print(
  fair_ranked_models %>% 
    select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
           b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
    head(15)
)

print(
  fair_ranked_models %>%
    filter(model_type=="2 Predators") %>% 
    select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
           b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
    head(15)
)
