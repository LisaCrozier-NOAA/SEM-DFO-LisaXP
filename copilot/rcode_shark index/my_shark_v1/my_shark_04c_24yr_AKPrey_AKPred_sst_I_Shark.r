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
      !is.na(x16_sar),
      !is.na(Pred1),
      !is.na(Pred2)
    )
  
  n_years <- nrow(df_model)
  if (n_years < 12) return(NULL)
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + Pred1 + Pred2
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
    
    b_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_term1_sar   = get_est("x16_sar", "Pred1", "est"),
    p_term1_sar   = get_est("x16_sar", "Pred1", "pvalue"),
    
    b_term2_sar   = get_est("x16_sar", "Pred2", "est"),
    p_term2_sar   = get_est("x16_sar", "Pred2", "pvalue")
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
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(x16_sar))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw
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

# -----------------------------------------------------------------------------
# 6. Display Model Selection Summary
# -----------------------------------------------------------------------------
cat("\n=== TOP 15 FAIRLY-EVALUATED SEM MODELS (N =", max_n, "YEARS) ===\n")
topmod <- fair_ranked_models %>% 
  select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
         b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
  head(15)

print(topmod)

# Frequency of New Shark Index in Top Models
shark_in_top20 <- sum(fair_ranked_models$contains_new_shark[1:min(20, nrow(fair_ranked_models))])
cat("\nThe new shark index (x15_shark_enso_roll2) appears in", shark_in_top20, "out of the top 20 fair models.\n")

# Inspect #1 Model
cat("\n=== BEST MODEL SUMMARY ===\n")
print(topmod[1, ])

topmod[,1:5]
# A tibble: 15 × 5
# model_type          term1_col                               term2_col                         n_years   aic
# 1 1 Predator + 1 Prey x15_shark_enso_roll2                    x21_pdo_djf                            24  106.
# 2 1 Predator + 1 Prey x15_shark_enso_roll2                    x21_sst_wgoa_coastwatch_junjulaug      24  107.
# 3 1 Predator + 1 Prey x15_shark_enso_roll2                    x12_dfa_biomass_euph_shelf_sum         24  109.
# 4 1 Predator + 1 Prey x15_arrowtooth_flounder_biomass_pred_ak x21_sst_wgoa_coastwatch_junjulaug      24  110.
# 5 1 Predator + 1 Prey x10_harbour_s_2yr_lead_ws               x21_sst_wgoa_coastwatch_junjulaug      24  110.
# 6 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x21_sst_wgoa_coastwatch_junjulaug      24  110.
# 7 1 Predator + 1 Prey x10_harbour_s_2yr_lead_ws               x21_pdo_djf                            24  111.
# 8 1 Predator + 1 Prey x10_dfa_ssl_est_wholerange_2yr_lead     x21_sst_wgoa_coastwatch_junjulaug      24  112.
# 9 1 Predator + 1 Prey x15_spiny_dogfish_go_a_pred_ak          x21_sst_wgoa_coastwatch_junjulaug      24  112.
# 10 1 Predator + 1 Prey x15_shark_catch_go_a_pred_ak            x21_sst_wgoa_coastwatch_junjulaug      24  113.
# 11 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x12_dfa_biomass_euph_shelf_sum         24  113.
# 12 1 Predator + 1 Prey x10_californian_s_l_2yr_lead_ws         x21_sst_wgoa_coastwatch_junjulaug      24  113.
# 13 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak           x21_sst_wgoa_coastwatch_junjulaug      24  113.
# 14 1 Predator + 1 Prey x15_shark_enso_roll2                    x12_copepod_biomass_w_go_a             24  114.
# 15 1 Predator + 1 Prey x15_sablefish_recruitment_pred_ak       x21_sst_wgoa_coastwatch_junjulaug      24  114.
# > 

head(topmod)
#    model_type          term1_col                               term2_col                         n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# 1 1 Predator + 1 Prey x15_shark_enso_roll2                    x21_pdo_djf                            24  106.             -16.5      -0.500   0.0000153      -0.576  0.00000267 TRUE              
# 2 1 Predator + 1 Prey x15_shark_enso_roll2                    x21_sst_wgoa_coastwatch_junjulaug      24  107.             -15.2      -0.383   0.000848       -0.596  0.00000499 TRUE              
# 3 1 Predator + 1 Prey x15_shark_enso_roll2                    x12_dfa_biomass_euph_shelf_sum         24  109.             -13.5      -0.563   0.0000122      -0.535  0.0000584  TRUE              
# 4 1 Predator + 1 Prey x15_arrowtooth_flounder_biomass_pred_ak x21_sst_wgoa_coastwatch_junjulaug      24  110.             -12.1      -0.339   0.00947        -0.663  0.00000324 FALSE             
# 5 1 Predator + 1 Prey x10_harbour_s_2yr_lead_ws               x21_sst_wgoa_coastwatch_junjulaug      24  110.             -12.0      -0.326   0.00999        -0.506  0.000443   FALSE             
# 6 1 Predator + 1 Prey x15_dfa_sleeper_shark_bsai_pred_ak      x21_sst_wgoa_coastwatch_junjulaug      24  110.             -12.0      -0.317   0.00960        -0.551  0.0000765  FALSE             
# 13 1 Predator + 1 Prey x15_sablefish_biomass_pred_ak           x21_sst_wgoa_coastwatch_junjulaug      24  113.             -9.50       0.255   0.0599         -0.618  0.0000290  FALSE             
# 14 1 Predator + 1 Prey x15_shark_enso_roll2                    x12_copepod_biomass_w_go_a             24  114.             -8.49      -0.274   0.0425         -0.419  0.00321    TRUE              
# 15 1 Predator + 1 Prey x15_sablefish_recruitment_pred_ak       x21_sst_wgoa_coastwatch_junjulaug      24  114.             -8.40       0.219   0.128          -0.645  0.0000304  FALSE             


top2pred <- fair_ranked_models %>% 
  filter(model_type=="2 Predators") %>% 
  select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
         b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
  head(15)
top2pred[1:5,1:5]

# model_type  term1_col                                     term2_col                          n_years   aic
# 1 2 Predators x15_salmon_shark_go_a_pred_ak                 x15_shark_enso_roll2                    24  117.
# 2 2 Predators x10_dfa_ssl_est_wholerange_2yr_lead           x15_dfa_sleeper_shark_bsai_pred_ak      24  117.
# 3 2 Predators x15_halibut_biomass_age8plus_2yr_lead_pred_ak x15_dfa_sleeper_shark_bsai_pred_ak      24  117.
# 4 2 Predators x10_harbour_s_2yr_lead_ws                     x15_shark_catch_go_a_pred_ak            24  118.
# 5 2 Predators x10_harbour_s_2yr_lead_ws                     x15_salmon_shark_go_a_pred_ak           24  118.

top2pred[1:5,]
# model_type  term1_col                      term2_col n_years   aic delta_aic_vs_base b_term1_sar p_term1_sar b_term2_sar p_term2_sar contains_new_shark
# 1 2 Predators x15_salmon_shark_go_a_pred_ak  x15_shar…      24  117.             -5.65      -0.333     0.0256       -0.487     0.00110 TRUE              
# 2 2 Predators x10_dfa_ssl_est_wholerange_2y… x15_dfa_…      24  117.             -5.62      -0.736     0.0287       -0.998     0.00234 FALSE             
# 3 2 Predators x15_halibut_biomass_age8plus_… x15_dfa_…      24  117.             -5.01       0.550     0.0353       -0.798     0.00179 FALSE             
# 4 2 Predators x10_harbour_s_2yr_lead_ws      x15_shar…      24  118.             -4.64      -0.657     0.00138      -0.330     0.120   FALSE             
# 5 2 Predators x10_harbour_s_2yr_lead_ws      x15_salm…      24  118.             -4.32      -0.429     0.00294      -0.195     0.175   FALSE             
