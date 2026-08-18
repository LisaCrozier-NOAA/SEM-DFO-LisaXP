#Base model for sar_uc:
# fit_base <- sem('
#   x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
#   sar_uc ~ x09_dfa_hake_age5plus
# sem_text <- '
#     x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
#     sar_uc ~ x09_dfa_hake_age5plus + Pred1 + Pred2
#   '

#Results:
  # The new shark index (x15_shark_enso_roll2) appears in 8 out of the top 20 fair models.
  # === BEST MODEL SUMMARY ===
  # model_type          term1_col            term2_col        n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
  # 1 Predator + 1 Prey x15_shark_enso_roll2 x13_egoa_herring      18  79.9             -12.6      -0.542     0.00208       0.463    0.000113 TRUE              
  
#===========================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
out_dir <- "copilot/outputs_4"

# Load shark predictions output
all_sims <- read.csv(file.path(out_dir, "shark_predictions.csv"))

top_shark_index <- all_sims %>%
  filter(Scenario == "enso_dj | z_roll2") %>%
  select(year, x15_shark_enso_roll2 = I_Shark)

# Load, clean, and scale guild dataset
guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% 
  clean_names() %>% 
  select(-any_of("x10_harbor_seal_cr_2yr_lead")) %>%
  left_join(top_shark_index, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

# -----------------------------------------------------------------------------
# 2. Identify Candidate Variable Sets
# -----------------------------------------------------------------------------
# Temperature (x21_) and Prey candidates (x12_, x13_, x14_)
ak_prey_cands <- names(guild.dfas1)[grepl("^x21_|^x12_|^x13_|^x14_", names(guild.dfas1))]

# Predator candidates (x10_, x15_ including your new x15_shark_enso_roll2)
ak_pred_cands <- names(guild.dfas1)[grepl("^x10_|^x15_", names(guild.dfas1))]

message("Found ", length(ak_prey_cands), " Prey/Temp candidates and ", 
        length(ak_pred_cands), " Predator candidates (including new shark index).")

# -----------------------------------------------------------------------------
# 3. Helper Function with Strict Year Tracking & Listwise Fitting
# -----------------------------------------------------------------------------
fit_sem_two_topdown <- function(pred1_col, pred2_col, df_data) {
  
  df_model <- df_data %>%
    filter(year >= 1998, year <= 2021) %>%
    mutate(
      Pred1 = .data[[pred1_col]],
      Pred2 = .data[[pred2_col]]
    ) %>%
    filter(
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(sar_uc),
      !is.na(Pred1),
      !is.na(Pred2)
    )
  
  n_years <- nrow(df_model)
  if (n_years < 12) return(NULL)
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    sar_uc ~ x09_dfa_hake_age5plus + Pred1 + Pred2
  '
  
  fit <- tryCatch({
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "listwise", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea"))
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
    
    b_hake_sar   = get_est("sar_uc", "x09_dfa_hake_age5plus", "est"),
    p_hake_sar   = get_est("sar_uc", "x09_dfa_hake_age5plus", "pvalue"),
    # b_cpue_sar    = get_est("sar_uc", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    # p_cpue_sar    = get_est("sar_uc", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_term1_sar   = get_est("sar_uc", "Pred1", "est"),
    p_term1_sar   = get_est("sar_uc", "Pred1", "pvalue"),
    
    b_term2_sar   = get_est("sar_uc", "Pred2", "est"),
    p_term2_sar   = get_est("sar_uc", "Pred2", "pvalue")
  )
}

# -----------------------------------------------------------------------------
# 4A. SWEEP 1: Two Predator Combinations
# -----------------------------------------------------------------------------
pred_pairs <- combn(ak_pred_cands, 2, simplify = FALSE)

res_pred_pairs <- purrr::map_dfr(pred_pairs, function(pair) {
  fit_sem_two_topdown(pair[1], pair[2], guild.dfas1)
}) %>%
  mutate(model_type = "2 Predators")

# -----------------------------------------------------------------------------
# 4B. SWEEP 2: One Predator + One Prey/Temp Combinations
# -----------------------------------------------------------------------------
prey_pred_grid <- expand.grid(
  pred_col = ak_pred_cands,
  prey_col = ak_prey_cands,
  stringsAsFactors = FALSE
)

res_pred_prey <- purrr::map2_dfr(
  prey_pred_grid$pred_col, 
  prey_pred_grid$prey_col, 
  function(p1, p2) {
    fit_sem_two_topdown(p1, p2, guild.dfas1)
  }
) %>%
  mutate(model_type = "1 Predator + 1 Prey")

# -----------------------------------------------------------------------------
# 5. Combine, Filter for Complete N, and Rank Fairly
# -----------------------------------------------------------------------------
all_competing_models <- bind_rows(res_pred_pairs, res_pred_prey)

# Identify maximum complete years across runs
max_n <- max(all_competing_models$n_years, na.rm = TRUE)
cat("Filtering models to evaluate strictly on complete data: N =", max_n, "years.\n")

# Baseline model evaluated on complete full dataset
df_base <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(sar_uc))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  sar_uc ~ x09_dfa_hake_age5plus
', data = df_base, std.lv = TRUE, missing = "listwise", warn = FALSE)

base_aic <- fitMeasures(fit_base, "aic")

# Rank ONLY models that have full sample size (drops truncated time-series bias)
fair_ranked_models <- all_competing_models %>%
  filter(n_years == max_n) %>%
  mutate(
    delta_aic_vs_base = aic - base_aic,
    contains_new_shark = (term1_col == "x15_shark_enso_roll2" | term2_col == "x15_shark_enso_roll2")
  ) %>%
  arrange(aic)

 #write.csv(fair_ranked_models,file.path(out_dir,"fair_ranked_models_x16_sar.csv"))
 write.csv(fair_ranked_models,file.path(out_dir,"fair_ranked_models_sar_uc.csv"))
 
# -----------------------------------------------------------------------------
# 6. Display Model Selection Summary
# -----------------------------------------------------------------------------
cat("\n=== TOP 15 FAIRLY-EVALUATED SEM MODELS (N =", max_n, "YEARS) ===\n")
topmod <- fair_ranked_models %>% 
  select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
         b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
  head(15)

top2pred <- fair_ranked_models %>% 
  filter(model_type=="2 Predators") %>% 
  select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
         b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
  head(15)

topmod[,1:5]
top2pred[,1:5]

print(topmod) %>% head()
print(top2pred) %>% head()

# Frequency of New Shark Index in Top Models
shark_in_top20 <- sum(fair_ranked_models$contains_new_shark[1:min(20, nrow(fair_ranked_models))])
cat("\nThe new shark index (x15_shark_enso_roll2) appears in", shark_in_top20, "out of the top 20 fair models.\n")

# Inspect #1 Model
cat("\n=== BEST MODEL SUMMARY ===\n")
print(topmod[1, ])


#Results=======================
# print(topmod) %>% head()
# # A tibble: 15 × 11
# model_type          term1_col                               term2_col                                n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>               <chr>                                   <chr>                                      <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 1 Predator + 1 Prey x15_shark_enso_roll2                    x13_egoa_herring                              18  79.9            -12.6       -0.542  0.00208          0.463  0.000113   TRUE              
# 2 1 Predator + 1 Prey x15_shark_enso_roll2                    x13_stka_herr_matbiom                         18  80.5            -12.0       -0.556  0.00184          0.449  0.000202   TRUE              
# 3 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x13_stka_herr_matbiom                         18  81.3            -11.2       -0.483  0.00334          0.418  0.000811   FALSE             
# 4 1 Predator + 1 Prey x15_pacific_cod_biomass_pred_ak         x13_hexagram_eai                              18  81.4            -11.0       -0.592  0.000328        -0.757  0.00000283 FALSE             
# 5 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x13_ammod_wai                                 18  81.8            -10.6       -1.66   0.00000408       1.99   0.000953   FALSE             
# 6 2 Predators         x15_sablefish_biomass_pred_ak           x15_shark_enso_roll2                          18  81.9            -10.6        0.450  0.000731        -0.678  0.000265   TRUE              
# 7 2 Predators         x15_salmon_shark_go_a_pred_ak           x15_shark_enso_roll2                          18  82.1            -10.4       -0.439  0.000811        -0.923  0.00000860 TRUE              
# 8 1 Predator + 1 Prey x15_shark_enso_roll2                    x13_pollock_biomass_go_aage3plus_pred_ak      18  82.5            -10.0       -1.51   0.00000587      -0.840  0.00116    TRUE              
# 9 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x13_egoa_herring                              18  82.8             -9.66      -0.414  0.0199           0.404  0.00290    FALSE             
# 10 1 Predator + 1 Prey x15_arrowtooth_flounder_biomass_pred_ak x13_stka_herr_matbiom                         18  83.0             -9.46      -0.413  0.0141           0.391  0.00362    FALSE             
# 11 2 Predators         x15_spiny_dogfish_go_a_pred_ak          x15_shark_catch_go_a_pred_ak                  18  83.2             -9.28       4.53   0.0000273       -4.25   0.0000761  FALSE             
# 12 1 Predator + 1 Prey x10_californian_s_l_2yr_lead_ws         x13_stka_herr_matbiom                         18  83.2             -9.24       0.331  0.0169           0.481  0.000209   FALSE             
# 13 1 Predator + 1 Prey x15_sablefish_recruitment_pred_ak       x13_hexagram_eai                              18  83.5             -9.02       0.501  0.00200         -0.540  0.000245   FALSE             
# 14 1 Predator + 1 Prey x15_shark_enso_roll2                    x12_copepod_biomass_w_go_a                    18  84.1             -8.39      -0.525  0.00850         -0.523  0.00440    TRUE              
# 15 1 Predator + 1 Prey x15_shark_enso_roll2                    x21_sst_wgoa_coastwatch_junjulaug             18  84.3             -8.15      -0.852  0.0000671       -0.705  0.00507    TRUE              
# # A tibble: 6 × 11
# model_type          term1_col                          term2_col             n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>               <chr>                              <chr>                   <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 1 Predator + 1 Prey x15_shark_enso_roll2               x13_egoa_herring           18  79.9             -12.6      -0.542  0.00208          0.463  0.000113   TRUE              
# 2 1 Predator + 1 Prey x15_shark_enso_roll2               x13_stka_herr_matbiom      18  80.5             -12.0      -0.556  0.00184          0.449  0.000202   TRUE              
# 3 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak x13_stka_herr_matbiom      18  81.3             -11.2      -0.483  0.00334          0.418  0.000811   FALSE             
# 4 1 Predator + 1 Prey x15_pacific_cod_biomass_pred_ak    x13_hexagram_eai           18  81.4             -11.0      -0.592  0.000328        -0.757  0.00000283 FALSE             
# 5 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak x13_ammod_wai              18  81.8             -10.6      -1.66   0.00000408       1.99   0.000953   FALSE             
# 6 2 Predators         x15_sablefish_biomass_pred_ak      x15_shark_enso_roll2       18  81.9             -10.6       0.450  0.000731        -0.678  0.000265   TRUE              
# > print(top2pred) %>% head()
# # A tibble: 15 × 11
# model_type  term1_col                               term2_col                                     n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>       <chr>                                   <chr>                                           <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 2 Predators x15_sablefish_biomass_pred_ak           x15_shark_enso_roll2                               18  81.9            -10.6        0.450   0.000731       -0.678  0.000265   TRUE              
# 2 2 Predators x15_salmon_shark_go_a_pred_ak           x15_shark_enso_roll2                               18  82.1            -10.4       -0.439   0.000811       -0.923  0.00000860 TRUE              
# 3 2 Predators x15_spiny_dogfish_go_a_pred_ak          x15_shark_catch_go_a_pred_ak                       18  83.2             -9.28       4.53    0.0000273      -4.25   0.0000761  FALSE             
# 4 2 Predators x15_arrowtooth_flounder_biomass_pred_ak x15_shark_catch_go_a_pred_ak                       18  84.9             -7.61      -1.18    0.000133       -0.686  0.0150     FALSE             
# 5 2 Predators x15_sablefish_biomass_pred_ak           x15_dfa_sleeper_shark_bsai_pred_ak                 18  85.0             -7.48       0.349   0.0163         -0.526  0.00365    FALSE             
# 6 2 Predators x15_salmon_shark_go_a_pred_ak           x15_dfa_sleeper_shark_bsai_pred_ak                 18  86.4             -6.10      -0.271   0.0448         -0.656  0.000577   FALSE             
# 7 2 Predators x15_arrowtooth_flounder_biomass_pred_ak x15_spiny_dogfish_go_a_pred_ak                     18  86.5             -5.96      -1.21    0.00148        -0.683  0.0510     FALSE             
# 8 2 Predators x10_harbour_s_2yr_lead_ws               x15_dfa_sleeper_shark_bsai_pred_ak                 18  86.5             -5.95       0.484   0.0480         -1.00   0.000348   FALSE             
# 9 2 Predators x10_harbour_s_2yr_lead_ws               x15_halibut_biomass_age8plus_2yr_lead_pred_ak      18  87.0             -5.50       0.778   0.0159         -1.74   0.000633   FALSE             
# 10 2 Predators x10_californian_s_l_2yr_lead_ws         x15_sablefish_biomass_pred_ak                      18  87.1             -5.37       0.358   0.0206          0.423  0.00581    FALSE             
# 11 2 Predators x15_arrowtooth_flounder_biomass_pred_ak x15_sablefish_biomass_pred_ak                      18  87.5             -4.98      -0.434   0.0270          0.270  0.104      FALSE             
# 12 2 Predators x15_arrowtooth_flounder_biomass_pred_ak x15_salmon_shark_go_a_pred_ak                      18  87.6             -4.85      -0.574   0.00184        -0.218  0.113      FALSE             
# 13 2 Predators x10_dfa_ssl_est_wholerange_2yr_lead     x15_sablefish_biomass_pred_ak                      18  87.8             -4.64       0.411   0.0352          0.375  0.0166     FALSE             
# 14 2 Predators x10_harbour_s_2yr_lead_ws               x15_arrowtooth_flounder_biomass_pred_ak            18  88.0             -4.50       0.339   0.147          -0.805  0.00158    FALSE             
# 15 2 Predators x15_shark_catch_go_a_pred_ak            x15_dfa_sleeper_shark_bsai_pred_ak                 18  88.1             -4.36      -0.358   0.156          -0.892  0.00212    FALSE             
# # A tibble: 6 × 11
# model_type  term1_col                               term2_col                          n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>       <chr>                                   <chr>                                <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 2 Predators x15_sablefish_biomass_pred_ak           x15_shark_enso_roll2                    18  81.9            -10.6        0.450   0.000731       -0.678  0.000265   TRUE              
# 2 2 Predators x15_salmon_shark_go_a_pred_ak           x15_shark_enso_roll2                    18  82.1            -10.4       -0.439   0.000811       -0.923  0.00000860 TRUE              
# 3 2 Predators x15_spiny_dogfish_go_a_pred_ak          x15_shark_catch_go_a_pred_ak            18  83.2             -9.28       4.53    0.0000273      -4.25   0.0000761  FALSE             
# 4 2 Predators x15_arrowtooth_flounder_biomass_pred_ak x15_shark_catch_go_a_pred_ak            18  84.9             -7.61      -1.18    0.000133       -0.686  0.0150     FALSE             
# 5 2 Predators x15_sablefish_biomass_pred_ak           x15_dfa_sleeper_shark_bsai_pred_ak      18  85.0             -7.48       0.349   0.0163         -0.526  0.00365    FALSE             
# 6 2 Predators x15_salmon_shark_go_a_pred_ak           x15_dfa_sleeper_shark_bsai_pred_ak      18  86.4             -6.10      -0.271   0.0448         -0.656  0.000577   FALSE             
# > 
#   > # Frequency of New Shark Index in Top Models
#   > shark_in_top20 <- sum(fair_ranked_models$contains_new_shark[1:min(20, nrow(fair_ranked_models))])
# > cat("\nThe new shark index (x15_shark_enso_roll2) appears in", shark_in_top20, "out of the top 20 fair models.\n")
# 
# The new shark index (x15_shark_enso_roll2) appears in 8 out of the top 20 fair models.
# > 
#   > # Inspect #1 Model
#   > cat("\n=== BEST MODEL SUMMARY ===\n")
# 
# === BEST MODEL SUMMARY ===
#   > print(topmod[1, ])
# # A tibble: 1 × 11
# model_type          term1_col            term2_col        n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>               <chr>                <chr>              <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 1 Predator + 1 Prey x15_shark_enso_roll2 x13_egoa_herring      18  79.9             -12.6      -0.542     0.00208       0.463    0.000113 TRUE              
# > 
