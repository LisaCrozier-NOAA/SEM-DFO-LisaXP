suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(glue)
})

# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
out_dir <- "copilot/outputs_4"

ak_dat<-read.csv("copilot/outputs_2/data_all_tested_columns_annual.csv");names(ak_dat)
ak2<-ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,mid_il_capelin,stka_herr_matbiom)


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

# Assuming 'all_sims' contains the scenario output from the top model (enso_dj | z_roll2)
top_shark_index <- all_sims %>%
  filter(Scenario == "enso_dj | z_roll2") %>%
  select(year, x15_shark_enso_roll2 = I_Shark)

# Load and clean guild dataset
guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% clean_names()

# Merge the top shark index into guild dataset by year
guild.dfas1 <- guild.dfas1 %>%
  left_join(top_shark_index, by = "year")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(glue)
})

# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
top_shark_index <- all_sims %>%
  filter(Scenario == "enso_dj | z_roll2") %>%
  select(year, x15_shark_enso_roll2 = I_Shark)

guild.dfas1 <- read.csv("data_Lisa/guild_dfas1.csv", row.names = NULL) %>% clean_names()

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
      !is.na(x16_sar),
      !is.na(x15_shark_enso_roll2),
      !is.na(AKPrey)
    )
  
  if (nrow(df_model) < 12) next
  
  # SEM Specification
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + AKPrey + x15_shark_enso_roll2
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
    
    b_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_prey_sar    = get_est("x16_sar", "AKPrey", "est"),
    p_prey_sar    = get_est("x16_sar", "AKPrey", "pvalue"),
    
    b_shark_sar   = get_est("x16_sar", "x15_shark_enso_roll2", "est"),
    p_shark_sar   = get_est("x16_sar", "x15_shark_enso_roll2", "pvalue")
  )
}

res_prey <- bind_rows(prey_results)

# -----------------------------------------------------------------------------
# 5. Baseline Comparison & Ranking
# -----------------------------------------------------------------------------
df_base <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), 
         !is.na(x16_sar), !is.na(x15_shark_enso_roll2))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2
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
ranked_prey_X16_sar<-ranked_prey


write.csv(ranked_prey, file.path(out_dir, "shark_SEM_modelcomparison.csv"))
          
print(ranked_prey_X16_sar %>% 
        select(prey_variable, aic, delta_aic_vs_base, b_prey_sar, p_prey_sar, b_shark_sar, p_shark_sar) %>% 
        head(10))


#   prey_variable                              aic delta_aic_vs_base b_prey_sar p_prey_sar b_shark_sar p_shark_sar
# 1 x12_dfa_biomass_euph_shelf_sum            113.         -10.4         -0.535  0.0000547      -0.563   0.0000120
# 2 x12_copepod_biomass_w_go_a                118.          -5.37        -0.419  0.00331        -0.274   0.0451   
# 3 x13_pollock_biomass_go_aage3plus_pred_ak  119.          -4.24        -0.536  0.00758        -0.751   0.000139 
# 4 x13_stka_herr_matbiom                     121.          -1.99         0.303  0.0372         -0.323   0.0251   
# 5 x13_egoa_herring                          121.          -1.82         0.298  0.0418         -0.319   0.0277   
# 6 x13_wgoa_cap_pcod                         122.          -1.42         0.309  0.0552         -0.300   0.0439   
# 7 x13_mid_il_capelin                        122.          -1.35         0.311  0.0579         -0.312   0.0342   
# 8 x13_gadid_wai                             122.          -1.00         0.269  0.0738         -0.441   0.00343  
# 9 x14_pink_salmon_north_america             123.           0.00811     -0.218  0.150          -0.320   0.0359   
# 10 x13_dfa_wgoa_dfa_seabirds                 123.           0.0145      -0.237  0.150          -0.276   0.0884   
# >