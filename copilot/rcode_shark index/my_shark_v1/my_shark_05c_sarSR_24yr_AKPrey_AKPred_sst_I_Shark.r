#Base model for sar_sr:
# fit_base <- sem('
#   x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
#   sar_sr ~ x09_dfa_hake_age5plus
# sem_text <- '
#     x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
#     sar_sr ~ x09_dfa_hake_age5plus + Pred1 + Pred2
#   '

#Results:
# The new shark index (x15_shark_enso_roll2) appears in 4 out of the top 20 fair models.
# 
# === BEST MODEL SUMMARY ===
# model_type          term1_col                     term2_col     n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak x13_ammod_wai      24  108.             -10.2       0.304      0.0179      -0.417    0.000253 FALSE             

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
      !is.na(sar_sr),
      !is.na(Pred1),
      !is.na(Pred2)
    )
  
  n_years <- nrow(df_model)
  if (n_years < 12) return(NULL)
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    sar_sr ~ x09_dfa_hake_age5plus + Pred1 + Pred2
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
    
    b_hake_sar   = get_est("sar_sr", "x09_dfa_hake_age5plus", "est"),
    p_hake_sar   = get_est("sar_sr", "x09_dfa_hake_age5plus", "pvalue"),
    # b_cpue_sar    = get_est("sar_sr", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    # p_cpue_sar    = get_est("sar_sr", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_term1_sar   = get_est("sar_sr", "Pred1", "est"),
    p_term1_sar   = get_est("sar_sr", "Pred1", "pvalue"),
    
    b_term2_sar   = get_est("sar_sr", "Pred2", "est"),
    p_term2_sar   = get_est("sar_sr", "Pred2", "pvalue")
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
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(sar_sr))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  sar_sr ~ x09_dfa_hake_age5plus
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
 write.csv(fair_ranked_models,file.path(out_dir,"fair_ranked_models_sar_sr.csv"))
 
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
# model_type          term1_col                                     term2_col                      n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# 1 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak                 x13_ammod_wai                       24  108.            -10.2       0.304      0.0179       -0.417   0.000253  FALSE             
# 2 1 Predator + 1 Prey x15_salmon_shark_go_a_pred_ak                 x13_ammod_wai                       24  110.             -9.09     -0.244      0.0389       -0.476   0.0000606 FALSE             
# 3 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak            x13_ammod_wai                       24  110.             -8.82      0.383      0.0328       -0.734   0.0000483 FALSE             
# 4 1 Predator + 1 Prey x15_shark_enso_roll2                          x13_gadid_wai                       24  111.             -7.56     -0.368      0.00327       0.413   0.00129   TRUE              
# 5 1 Predator + 1 Prey x10_dfa_ssl_est_wholerange_2yr_lead           x13_ammod_wai                       24  111.             -7.43     -0.287      0.116        -0.624   0.000268  FALSE             
# 6 1 Predator + 1 Prey x15_shark_enso_roll2                          x12_dfa_biomass_euph_shelf_sum      24  112.             -6.96     -0.414      0.00169      -0.447   0.00204   TRUE              
# 7 1 Predator + 1 Prey x10_californian_s_l_2yr_lead_ws               x13_ammod_wai                       24  112.             -6.93     -0.199      0.170        -0.542   0.000214  FALSE             
# 8 1 Predator + 1 Prey x15_pacific_cod_biomass_pred_ak               x13_ammod_wai                       24  112.             -6.79     -0.186      0.194        -0.423   0.000534  FALSE             
# 9 1 Predator + 1 Prey x10_northern_f_s_2yr_lead_ws                  x13_ammod_wai                       24  112.             -6.59     -0.186      0.222        -0.520   0.000285  FALSE             
# 10 1 Predator + 1 Prey x15_halibut_biomass_age8plus_2yr_lead_pred_ak x13_ammod_wai                       24  112.             -6.20      0.171      0.295        -0.527   0.000641  FALSE             
# 11 1 Predator + 1 Prey x10_harbour_s_2yr_lead_ws                     x13_gadid_wai                       24  113.             -5.77     -0.314      0.0117        0.292   0.0224    FALSE             
# 12 1 Predator + 1 Prey x15_shark_enso_roll2                          x21_pdo_djf                         24  113.             -5.76     -0.349      0.00682      -0.431   0.00500   TRUE              
# 13 1 Predator + 1 Prey x15_shark_enso_roll2                          x13_ammod_wai                       24  113.             -5.72      0.143      0.432        -0.535   0.00324   TRUE              
# 14 1 Predator + 1 Prey x15_sablefish_recruitment_pred_ak             x13_ammod_wai                       24  113.             -5.71      0.114      0.455        -0.401   0.00193   FALSE             
# 15 1 Predator + 1 Prey x15_arrowtooth_flounder_biomass_pred_ak       x13_ammod_wai                       24  113.             -5.43      0.0959     0.602        -0.482   0.00305   FALSE             
# # A tibble: 6 × 11
# model_type          term1_col                           term2_col                      n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>               <chr>                               <chr>                            <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak       x13_ammod_wai                       24  108.            -10.2        0.304     0.0179       -0.417   0.000253  FALSE             
# 2 1 Predator + 1 Prey x15_salmon_shark_go_a_pred_ak       x13_ammod_wai                       24  110.             -9.09      -0.244     0.0389       -0.476   0.0000606 FALSE             
# 3 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak  x13_ammod_wai                       24  110.             -8.82       0.383     0.0328       -0.734   0.0000483 FALSE             
# 4 1 Predator + 1 Prey x15_shark_enso_roll2                x13_gadid_wai                       24  111.             -7.56      -0.368     0.00327       0.413   0.00129   TRUE              
# 5 1 Predator + 1 Prey x10_dfa_ssl_est_wholerange_2yr_lead x13_ammod_wai                       24  111.             -7.43      -0.287     0.116        -0.624   0.000268  FALSE             
# 6 1 Predator + 1 Prey x15_shark_enso_roll2                x12_dfa_biomass_euph_shelf_sum      24  112.             -6.96      -0.414     0.00169      -0.447   0.00204   TRUE              
# > print(top2pred) %>% head()
# # A tibble: 15 × 11
# model_type  term1_col                       term2_col                                     n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>       <chr>                           <chr>                                           <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 2 Predators x10_harbour_s_2yr_lead_ws       x15_sablefish_biomass_pred_ak                      24  114.             -4.27     -0.292       0.0246      0.272       0.0639 FALSE             
# 2 2 Predators x15_sablefish_biomass_pred_ak   x15_shark_enso_roll2                               24  115.             -3.48      0.329       0.0256     -0.263       0.0438 TRUE              
# 3 2 Predators x10_harbour_s_2yr_lead_ws       x15_salmon_shark_go_a_pred_ak                      24  116.             -2.76     -0.342       0.0101     -0.175       0.186  FALSE             
# 4 2 Predators x15_salmon_shark_go_a_pred_ak   x15_shark_enso_roll2                               24  116.             -2.36     -0.272       0.0555     -0.351       0.0138 TRUE              
# 5 2 Predators x10_harbour_s_2yr_lead_ws       x15_pacific_cod_biomass_pred_ak                    24  117.             -1.84     -0.310       0.0233     -0.142       0.376  FALSE             
# 6 2 Predators x15_sablefish_biomass_pred_ak   x15_shark_catch_go_a_pred_ak                       24  117.             -1.70      0.279       0.0737      0.216       0.154  FALSE             
# 7 2 Predators x10_harbour_s_2yr_lead_ws       x15_shark_enso_roll2                               24  117.             -1.40     -0.271       0.100      -0.0982      0.552  TRUE              
# 8 2 Predators x10_harbour_s_2yr_lead_ws       x15_sablefish_recruitment_pred_ak                  24  117.             -1.40     -0.295       0.0461      0.102       0.560  FALSE             
# 9 2 Predators x10_harbour_s_2yr_lead_ws       x10_northern_f_s_2yr_lead_ws                       24  117.             -1.35     -0.370       0.0174     -0.0900      0.588  FALSE             
# 10 2 Predators x10_californian_s_l_2yr_lead_ws x10_harbour_s_2yr_lead_ws                          24  117.             -1.26     -0.0695      0.653      -0.363       0.0195 FALSE             
# 11 2 Predators x15_sablefish_biomass_pred_ak   x15_spiny_dogfish_go_a_pred_ak                     24  117.             -1.20      0.270       0.0921      0.192       0.221  FALSE             
# 12 2 Predators x10_harbour_s_2yr_lead_ws       x15_halibut_biomass_age8plus_2yr_lead_pred_ak      24  118.             -1.15     -0.313       0.0337     -0.0451      0.773  FALSE             
# 13 2 Predators x10_harbour_s_2yr_lead_ws       x15_shark_catch_go_a_pred_ak                       24  118.             -1.11     -0.370       0.104      -0.0553      0.825  FALSE             
# 14 2 Predators x10_harbour_s_2yr_lead_ws       x15_arrowtooth_flounder_biomass_pred_ak            24  118.             -1.07     -0.317       0.0762     -0.0208      0.918  FALSE             
# 15 2 Predators x10_harbour_s_2yr_lead_ws       x10_dfa_ssl_est_wholerange_2yr_lead                24  118.             -1.07     -0.323       0.0427      0.0120      0.944  FALSE             
# # A tibble: 6 × 11
# model_type  term1_col                     term2_col                       n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# <chr>       <chr>                         <chr>                             <int> <dbl>             <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <lgl>             
#   1 2 Predators x10_harbour_s_2yr_lead_ws     x15_sablefish_biomass_pred_ak        24  114.             -4.27      -0.292      0.0246       0.272      0.0639 FALSE             
# 2 2 Predators x15_sablefish_biomass_pred_ak x15_shark_enso_roll2                 24  115.             -3.48       0.329      0.0256      -0.263      0.0438 TRUE              
# 3 2 Predators x10_harbour_s_2yr_lead_ws     x15_salmon_shark_go_a_pred_ak        24  116.             -2.76      -0.342      0.0101      -0.175      0.186  FALSE             
# 4 2 Predators x15_salmon_shark_go_a_pred_ak x15_shark_enso_roll2                 24  116.             -2.36      -0.272      0.0555      -0.351      0.0138 TRUE              
# 5 2 Predators x10_harbour_s_2yr_lead_ws     x15_pacific_cod_biomass_pred_ak      24  117.             -1.84      -0.310      0.0233      -0.142      0.376  FALSE             
# 6 2 Predators x15_sablefish_biomass_pred_ak x15_shark_catch_go_a_pred_ak         24  117.             -1.70       0.279      0.0737       0.216      0.154  FALSE             
# > 
# 
# The new shark index (x15_shark_enso_roll2) appears in 4 out of the top 20 fair models.
# 
# === BEST MODEL SUMMARY ===
# model_type          term1_col                     term2_col     n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak x13_ammod_wai      24  108.             -10.2       0.304      0.0179      -0.417    0.000253 FALSE             
