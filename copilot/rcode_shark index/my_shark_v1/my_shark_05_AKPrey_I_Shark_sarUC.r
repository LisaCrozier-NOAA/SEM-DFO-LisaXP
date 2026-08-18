suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(glue)
})

# -----------------------------------------------------------------------------
# 0. Organize database
# -----------------------------------------------------------------------------

# ak_dat<-read.csv("copilot/outputs_2/data_all_tested_columns_annual.csv");names(ak_dat)
# ak2<-ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom)


# sem_master_data<-read.csv("outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
# guild.dfas1<-read.csv("outputs_4/guild.dfa.NCC.AK.csv",row.names=1) %>%
#   mutate( sarSR=sem_master_data$sarSR,
#           sarUC=sem_master_data$sarUC) %>%
#   inner_join(trend_df_1998_topbio_allwgoa,by="year") %>%
#   inner_join(trend_df_1998_herring_egoa,by="year") %>%
#   inner_join(ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom),by="year") 
# names(guild.dfas1)
# 
# guild.dfasAK<-guild.dfas1 %>% rename(
#   X13_wgoa_cap.pcod=wgoa_cap.pcod,
#   X13_egoa_herring = egoa_herring,
#   X13_mid_il_capelin = mid_il_capelin,
#   X13_stka_herr_matbiom=stka_herr_matbiom)
# 
# write.csv(guild.dfasAK,file.path("data_Lisa/guild.dfasAK.csv"))




# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
out_dir <- "copilot/outputs_4"

all_sims<-read.csv(file.path(out_dir, "shark_predictions.csv"))

top_shark_index <- all_sims %>%
  filter(Scenario == "enso_dj | z_roll2") %>%
  select(year, x15_shark_enso_roll2 = I_Shark)

guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% clean_names()

# Merge the top shark index into guild dataset by year
guild.dfas1 <- guild.dfas1 %>%
  left_join(top_shark_index, by = "year")

# -----------------------------------------------------------------------------
# 2. SCALE ENTIRE DATASET (Except year)
# -----------------------------------------------------------------------------
guild.dfas1 <- guild.dfas1 %>%
  mutate(across(-year, ~ as.vector(scale(.))))

# -----------------------------------------------------------------------------
# 3. Identify AK Prey Candidates (x12_, x13_, x14_)
# -----------------------------------------------------------------------------
ak_prey_cands <- names(guild.dfas1)[grepl("^x12_|^x13_|^x14_", names(guild.dfas1))]
ak_pred_cands <- names(guild.dfas1)[grepl("^x15_|^x13_|^x14_", names(guild.dfas1))]

if (length(ak_prey_cands) == 0) {
  stop("No prey candidates found matching prefixes x12_, x13_, or x14_ in guild.dfas1.")
}

message("Found ", length(ak_prey_cands), " candidate AK Prey columns to evaluate.")

# -----------------------------------------------------------------------------
# 4. Model Selection Loop Across AK Prey Candidates
# -----------------------------------------------------------------------------
prey_results <- list()

for (i in seq_along(ak_prey_cands)) {
  prey_col <- ak_prey_cands[i]
  
  df_model <- guild.dfas1 %>%
    filter(year >= 1998, year <= 2021) %>%
    mutate(AKPrey = .data[[prey_col]]) %>%
    filter(
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(sar_uc),
#      !is.na(x16_sar),
      !is.na(x15_shark_enso_roll2),
      !is.na(AKPrey)
    )
  
  if (nrow(df_model) < 12) next
  
  # SEM Specification
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    sar_uc ~ x07_dfa_cpue_int_spr_jun_hw + AKPrey + x15_shark_enso_roll2
  '
  
  fit <- tryCatch({
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit) || !lavInspect(fit, "converged")) next
  
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea"))
  pe <- parameterEstimates(fit)
  
  get_est <- function(lhs_v, rhs_v, field) {
    val <- pe %>% filter(lhs == lhs_v, op == "~", rhs == rhs_v) %>% pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  prey_results[[i]] <- tibble(
    prey_variable = prey_col,
    n             = nrow(df_model),
    aic           = fm[["aic"]],
    bic           = fm[["bic"]],
    cfi           = fm[["cfi"]],
    rmsea         = fm[["rmsea"]],
    
    b_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    
    b_cpue_sar    = get_est("sar_uc", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar    = get_est("sar_uc", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_prey_sar    = get_est("sar_uc", "AKPrey", "est"),
    p_prey_sar    = get_est("sar_uc", "AKPrey", "pvalue"),
    
    b_shark_sar   = get_est("sar_uc", "x15_shark_enso_roll2", "est"),
    p_shark_sar   = get_est("sar_uc", "x15_shark_enso_roll2", "pvalue")
  )
}

res_prey <- bind_rows(prey_results)

# -----------------------------------------------------------------------------
# 5. Baseline Comparison & Ranking
# -----------------------------------------------------------------------------
df_base <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), 
         !is.na(sar_uc), !is.na(x15_shark_enso_roll2))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  sar_uc ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2
', data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)

base_aic <- fitMeasures(fit_base, "aic")

ranked_prey <- res_prey %>%
  mutate(
    delta_aic_vs_base = aic - base_aic,
    prey_sign         = case_when(
      b_prey_sar > 0  ~ "Positive",
      b_prey_sar < 0  ~ "Negative",
      TRUE            ~ "Zero"
    ),
    sig_prey          = !is.na(p_prey_sar) & p_prey_sar < 0.05
  ) %>%
  arrange(aic)

# Print top 10 models
print(ranked_prey %>% 
        select(prey_variable, aic, delta_aic_vs_base, b_prey_sar, p_prey_sar, b_shark_sar, p_shark_sar) %>% 
        head(10))


print(t(ranked_prey[1:2,]))
#ranked_prey_X16_sar<-ranked_prey
ranked_prey_sarUC<-ranked_prey

write.csv(ranked_prey, file.path(out_dir, "shark_SEM_sarUC_modelcomparison.csv"))


#for sar_UC response variable:

# prey_variable                              aic delta_aic_vs_base b_prey_sar p_prey_sar b_shark_sar p_shark_sar
# 1 x12_copepod_biomass_w_go_a                86.6            -8.53      -0.593   0.000155      -0.510  0.0115    
# 2 x13_pollock_biomass_go_aage3plus_pred_ak  87.7            -7.40      -0.963   0.000442      -1.64   0.00000488
# 3 x12_dfa_biomass_euph_shelf_sum            90.0            -5.06      -0.487   0.00328       -0.878  0.000237  
# 4 x13_stka_herr_matbiom                     91.8            -3.33       0.384   0.0128        -0.556  0.0169    
# 5 x13_mid_il_capelin                        91.9            -3.24       0.482   0.0136        -0.548  0.0188    
# 6 x13_egoa_herring                          92.7            -2.37       0.349   0.0262        -0.550  0.0215    
# 7 x13_wgoa_cap_pcod                         93.0            -2.08       0.396   0.0323        -0.530  0.0288    
# 8 x13_hexagram_eai                          93.5            -1.64      -0.384   0.0447        -0.414  0.113     
# 9 x12_copepod_com_e_go_a                    95.4             0.304      0.256   0.182         -0.536  0.0407    
# 10 x14_pink_salmon                           95.8             0.669      0.240   0.240         -0.536  0.0437    

# prey_variable     "x12_copepod_biomass_w_go_a" "x13_pollock_biomass_go_aage3plus_pred_ak"
# n                 "18"                         "18"                                      
# aic               "86.57066"                   "87.69928"                                
# bic               "93.69363"                   "94.82225"                                
# cfi               "1.0000000"                  "0.9800978"                               
# rmsea             "0.00000000"                 "0.09336129"                              
# b_hake_cpue       "-0.5913743"                 "-0.5913743"                              
# p_hake_cpue       "0.0001764453"               "0.0001764453"                            
# b_cpue_sar        "0.2414923"                  "0.3382555"                               
# p_cpue_sar        "0.13550242"                 "0.03380262"                              
# b_prey_sar        "-0.5933115"                 "-0.9628024"                              
# p_prey_sar        "0.0001551168"               "0.0004420633"                            
# b_shark_sar       "-0.5103214"                 "-1.6382456"                              
# p_shark_sar       "1.152894e-02"               "4.880882e-06"                            
# delta_aic_vs_base "-8.529565"                  "-7.400945"                               
# prey_sign         "Negative"                   "Negative"                                
# sig_prey          "TRUE"                       "TRUE"                         