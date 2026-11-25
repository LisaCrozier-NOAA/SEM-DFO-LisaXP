

#RESULTS
# Model                                     npar    df Chisq p_Chisq   AIC   BIC   CFI  SRMR Growth_R2 CPUE_R2 SAR_R2 Notes                                                                        
# 1 M2fish_preyAK_WGOAmidtrophic                18    17  65.1       0  277.  298. 0.644 0.23      0.44    0.749  0.63  2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)     
# 2 M2fish_preyAK_pinkNA                        18    17  67.4       0  278.  299. 0.629 0.229     0.455   0.749  0.565 2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth            
# 3 M2fish_preyAK_pinkAsia                      18    17  82.2       0  279.  300. 0.565 0.241     0.438   0.749  0.555 2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)     
# 4 M2fish_preyAK_WGOAmid_preyNCC               30    42 155.        0  386.  422. 0.699 0.201     0.905   0.749  0.565 fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth        
# 5 M2fish_preyNCC                              29    34 107.        0  390.  424. 0.779 0.181     0.898   0.749  0.486 fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth        
# 6 M1mammal_noPreyAK                           25    35 126.        0  398.  428. 0.663 0.365     0.451   0.749  0.823 5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE                    
# 7 M2fish_preyAK_latentEGOA                    28    32 101.        0  478.  511. 0.598 0.252     0.453   0.749  0.71  2-fish/shark predAK + EGOA preyAK, no growth                                 
# 8 M1mammals_preyNCC                           38    64 178.        0  506.  550. 0.758 0.159     0.898   0.749  0.701 5-mammal predAK (var floor 0.01), latent NCC, growth, no preyAK, fixed.x=TRUE
# 9 M2fish_preyAK_latentEGOA_preyNCC            41    61 173.        0  583.  631. 0.718 0.173     0.899   0.749  0.576 fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth        
# 10 M2fish_preyAK_WGOAmid_preyNCC_predNCC       44    86 309.        0  639.  691. 0.575 0.258     0.901   0.749  0.563 fish/shark predAK + WGOAmid observed + NCC prey+ NCC pred, no growth         
# 11 M2fish_preyNCC_predNCC                      43    74 255.        0  642.  693. 0.621 0.251     0.896   0.749  0.488 fish/shark + NCC prey + NCC pred, no growth, no AK prey                      
# 12 M1mammals_preyNCC_predNCC                   52   100 251.        0  781.  842. 0.713 0.135     0.895   0.674  0.642 Mammals,latent NCCpred                                                       
# 13 M2fish_preyAK_latentEGOA_preyNCC_predNCC    56    96 215.        0  840.  906. 0.729 0.142     0.889   0.829  0.667 fish/shark predAK + latent EGOA preyAK + NCC prey + NCC pred, no growth      

save.image(file.path(data_out_dir,"07_latent_model_fits.Rdata"))

# ==============================================================================
# Script: 02_analysis/07f_fit_latent_models_NCCpred.R
# Step 1: get raw data from allData.rds
# Step 2: fit SEM
# ==============================================================================


library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)

# ==============================================================================
# 0. DIRECTORY & DATA SETUP
# ==============================================================================
proj_dir      <- file.path(getwd(), "Rcode_for_paper")
master_out    <- file.path(proj_dir, "Outputs")
data_out_dir  <- file.path(master_out, "data")
table_out_dir <- file.path(master_out, "tables")
fig_out_dir   <- file.path(master_out, "figures")
plot_out_dir  <- file.path(master_out, "plots")

data_in_path  <- file.path(data_out_dir, "sem_master_peff_zscored.csv")
data_out_path  <- file.path(data_out_dir,"sem_master_peff_zscored_nccpred_latent.csv")
all_data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/allData.rds"
cw_path       <- file.path(proj_dir, "metadata/crosswalk_names.csv")
fn_path       <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/functions/sem_output_summary_fxn.r"

if (file.exists(fn_path)) source(fn_path)

all_data   <- readRDS(all_data_path)
cw_df      <- read.csv(cw_path, stringsAsFactors = FALSE)
sem_master <- read.csv(data_in_path, stringsAsFactors = FALSE)

# ==============================================================================
# 1. PARSE & EXTRACT CONSTITUENTS (<= 2 MISSING YEARS, 1998-2021)
# ==============================================================================
parse_cw_tokens <- function(node_key, cw) {
  row_match <- cw %>%
    filter(LisaName == node_key | DFAname == node_key | dfa_cols == node_key)
  raw_tokens <- trimws(unlist(strsplit(row_match$rankedIndicators[1], "\\s*-\\s*|\\s*,\\s*")))
  raw_tokens <- raw_tokens[raw_tokens != ""]
  unique(gsub("_[0-9]{4}$", "", raw_tokens))
}

raw_x05_tokens <- parse_cw_tokens("X05_DFA_abundSardine", cw_df)
raw_x09_tokens <- parse_cw_tokens("X09_DFA_HakeAge5Plus", cw_df)
raw_x12_tokens <- parse_cw_tokens("X12_DFA_EGOA_DFA_lowerTrophic", cw_df)
raw_x13_tokens <- parse_cw_tokens("X13_DFA_WGOA_DFA_midTrophic", cw_df)
raw_x15_tokens <- parse_cw_tokens("X15_DFA_halibutBiomassAge8plus_predAK_adultyr", cw_df)

extract_constituents <- function(all_df, token_list, node_prefix, max_na = 2) {
  # Ensure a complete 1998-2021 year grid so missing years count properly as NAs
  year_grid <- tibble(year = 1998:2021)
  
  df_wide <- all_df %>%
    mutate(
      year            = year(ymd(date)),
      shortName_clean = gsub("_[0-9]{4}$", "", shortName)
    ) %>%
    filter(
      year >= 1998, year <= 2021,
      yearType == "smolt",
      shortName_clean %in% token_list
    ) %>%
    select(year, shortName = shortName_clean, finalVal) %>%
    distinct(year, shortName, .keep_all = TRUE) %>%
    pivot_wider(names_from = shortName, values_from = finalVal)
  
  df_wide <- year_grid %>%
    left_join(df_wide, by = "year") %>%
    arrange(year)
  
  cand_cols <- setdiff(names(df_wide), "year")
  col_nas   <- sapply(cand_cols, function(col) sum(is.na(df_wide[[col]])))
  
  kept_cols    <- names(col_nas[col_nas <= max_na])
  dropped_cols <- col_nas[col_nas > max_na]
  
  if (length(dropped_cols) > 0) {
    message(node_prefix, " dropped (> ", max_na, " NAs): ",
            paste0(names(dropped_cols), " (", dropped_cols, " NAs)", collapse = ", "))
  }
  
  df_wide %>%
    select(year, all_of(kept_cols)) %>%
    rename_with(~ paste0(node_prefix, .x), -year)
}

df_x05_comp <- extract_constituents(all_data, raw_x05_tokens, "X05_latent_", max_na = 2)
df_x09_comp <- extract_constituents(all_data, raw_x09_tokens, "X09_latent_", max_na = 2)
df_x12_comp <- extract_constituents(all_data, raw_x12_tokens, "X12_latent_", max_na = 2)
df_x13_comp <- extract_constituents(all_data, raw_x13_tokens, "X13_latent_", max_na = 2)
df_x15_comp <- extract_constituents(all_data, raw_x15_tokens, "X15_latent_", max_na = 2)

dfa_constituents <- list(df_x05_comp, df_x09_comp, df_x12_comp, df_x13_comp, df_x15_comp) %>%
  reduce(full_join, by = "year") %>%
  arrange(year)

# Z-scale all latent constituent columns (leaving year unscaled)
dfa_constituents_z <- dfa_constituents %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# Join onto sem_master_peff_zscored.csv
sem_data_merged <- sem_master %>%
  left_join(dfa_constituents_z, by = "year")

# Check which X09_latent_ columns survived the <= 2 missing years filter
print(names(df_x09_comp))






# ==============================================================================
# STEP 1: EXTRACT & Z-SCALE WGOA MID-TROPHIC (X13) CONSTITUENTS
# ==============================================================================
parse_cw_tokens <- function(node_key, cw) {
  row_match <- cw %>%
    filter(LisaName == node_key | DFAname == node_key | dfa_cols == node_key)
  raw_tokens <- trimws(unlist(strsplit(row_match$rankedIndicators[1], "\\s*-\\s*|\\s*,\\s*")))
  raw_tokens <- raw_tokens[raw_tokens != ""]
  unique(gsub("_[0-9]{4}$", "", raw_tokens))
}

extract_constituents <- function(all_df, token_list, node_prefix, max_na = 2) {
  year_grid <- tibble(year = 1998:2021)
  
  df_wide <- all_df %>%
    mutate(
      year            = year(ymd(date)),
      shortName_clean = gsub("_[0-9]{4}$", "", shortName)
    ) %>%
    filter(
      year >= 1998, year <= 2021,
      yearType == "smolt",
      shortName_clean %in% token_list
    ) %>%
    select(year, shortName = shortName_clean, finalVal) %>%
    distinct(year, shortName, .keep_all = TRUE) %>%
    pivot_wider(names_from = shortName, values_from = finalVal)
  
  df_wide <- year_grid %>%
    left_join(df_wide, by = "year") %>%
    arrange(year)
  
  cand_cols    <- setdiff(names(df_wide), "year")
  col_nas      <- sapply(cand_cols, function(col) sum(is.na(df_wide[[col]])))
  kept_cols    <- names(col_nas[col_nas <= max_na])
  dropped_cols <- col_nas[col_nas > max_na]
  
  if (length(dropped_cols) > 0) {
    message(node_prefix, " dropped (> ", max_na, " NAs): ",
            paste0(names(dropped_cols), " (", dropped_cols, " NAs)", collapse = ", "))
  }
  
  df_wide %>%
    select(year, all_of(kept_cols)) %>%
    rename_with(~ paste0(node_prefix, .x), -year)
}

# Parse tokens for X13_DFA_WGOA_DFA_midTrophic and extract <= 2 NA series
raw_x13_tokens <- parse_cw_tokens("X13_DFA_WGOA_DFA_midTrophic", cw_df)
print(raw_x13_tokens)

df_x13_comp <- extract_constituents(all_data, raw_x13_tokens, "X13_latent_", max_na = 2)

# Z-scale and merge onto sem_data_merged (dropping old X13_latent_ cols first if re-running)
df_x13_z <- df_x13_comp %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

sem_data_merged <- sem_data_merged %>%
  select(-starts_with("X13_latent_")) %>%
  left_join(df_x13_z, by = "year")

# Save updated master file so X13_latent_* columns persist
write_csv(sem_data_merged, data_out_path)

# Print exact column names to paste into preyAK =~ and preyAK_ind
x13_latent_cols <- grep("^X13_latent_", names(sem_data_merged), value = TRUE)
cat("\nSurviving X13_latent_ columns:\n")
cat(paste(x13_latent_cols, collapse = " +\n"), "\n")

# ==============================================================================
# EXTRACT & Z-SCALE EGOA LOWER TROPHIC (X12) CONSTITUENTS
# ==============================================================================
raw_x12_tokens <- parse_cw_tokens("X12_DFA_EGOA_DFA_lowerTrophic", cw_df)
print(raw_x12_tokens)

df_x12_comp <- extract_constituents(all_data, raw_x12_tokens, "X12_latent_", max_na = 2)

# Z-scale and merge onto sem_data_merged (dropping old X12_latent_ cols first if re-running)
df_x12_z <- df_x12_comp %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

sem_data_merged <- sem_data_merged %>%
  select(-starts_with("X12_latent_")) %>%
  left_join(df_x12_z, by = "year")

# Persist updated dataset
write_csv(sem_data_merged, data_out_path)

# Print exact column names to paste into preyAK =~ and build_master_layout(preyAK_ind = ...)
x12_latent_cols <- grep("^X12_latent_", names(sem_data_merged), value = TRUE)
cat("\nSurviving X12_latent_ columns:\n")
cat(paste(x12_latent_cols, collapse = " +\n"), "\n")

# ==============================================================================
# EXTRACT & Z-SCALE SARDINE / NCC PREY (X05) CONSTITUENTS
# ==============================================================================
raw_x05_tokens <- parse_cw_tokens("X05_DFA_abundSardine", cw_df)
print(raw_x05_tokens)

df_x05_comp <- extract_constituents(all_data, raw_x05_tokens, "X05_latent_", max_na = 2)

# Z-scale and merge onto sem_data_merged
df_x05_z <- df_x05_comp %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

sem_data_merged <- sem_data_merged %>%
  select(-starts_with("X05_latent_")) %>%
  left_join(df_x05_z, by = "year")

# Persist updated dataset
write_csv(sem_data_merged, data_out_path)

# Print exact surviving column names to swap into preyNCC =~
x05_latent_cols <- grep("^X05_latent_", names(sem_data_merged), value = TRUE)
cat("\nSurviving X05_latent_ columns:\n")
cat(paste(x05_latent_cols, collapse = " +\n"), "\n")


#SAVE NEW DATA FILE==============
write.csv(sem_data_merged, file.path(data_out_dir,"sem_master_peff_zscored_nccpred_latent.csv"),row.names = FALSE)

# ==============================================================================
# 2. SEM MODEL DEFINITION (predNCC with X09_latent_ constituents)
# ==============================================================================
# NOTE: If any of the 8 X09_latent_ indicators below had > 2 missing years and
# were dropped above, remove that line from predNCC =~ before fitting.
################################################
#FIT M1 MODELS W/ NCC PRED LATENT------------
##################################################
model_M1mammals_preyNCC_predNCC <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyNCC =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
             X02_DFA_SeaNettle_smoltyr +
             X04_marketsquid_GAM_smoltyr +
             X05_DFA_abundSardine_smoltyr

  predNCC =~ X09_latent_HakeAge5Plus +
             X09_latent_pacificmackerel_GAM +
             X09_latent_jackmackerel_GAM +
             X09_latent_hakeFisheryWA  

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M1mammals_preyNCC_predNCC <- sem(
  model_M1mammals_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_M1mammals_preyNCC_predNCC)
summary(fit_M1mammals_preyNCC_predNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
# Log kept model--------
log_sem_model(fit_M1mammals_preyNCC_predNCC,   "M1mammals_preyNCC_predNCC",   "Mammals,latent NCCpred")
print(model_comparison_table)
write.csv(model_comparison_table,file.path(data_out_dir,"latent_model_comparison_table.csv"),row.names=FALSE)

##############################
#EXPAND SARDINE DFA INTO LATENT CONSITUENTS------------
model_M1mammals_latentsardine_predNCC <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyNCC =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
             X02_DFA_SeaNettle_smoltyr +
             X04_marketsquid_GAM_smoltyr +
             X05_latent_northAnch_NCC +
              X05_latent_sardine_NCC +
              X05_latent_herring_GAM +
              X05_latent_sardine_GAM 

  predNCC =~ X09_latent_HakeAge5Plus +
             X09_latent_pacificmackerel_GAM +
             X09_latent_jackmackerel_GAM +
             X09_latent_hakeFisheryWA  

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M1mammals_latentsardine_predNCC <- sem(
  model_M1mammals_latentsardine_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)
fit<-fit_M1mammals_latentsardine_predNCC

sem_output_summary_fxn(fit)
parameterEstimates(fit, standardized = TRUE) %>%
  as_tibble() %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.lv, std.all) %>%
  print(n = Inf)
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

################################################
#FIT M2 MODELS W/ NCC PRED LATENT------------
##################################################
# ------------------------------------------------------------------------------
# M2 var 1: group 1: 2-Fish/Shark predAK + Multi-Indicator EGOA preyAK (k3_group == 1)----------
# ------------------------------------------------------------------------------

#replaced X12_DFA_EGOA_DFA_lowerTrophic_smoltyr with :
              # X12_latent_crestedAuk_WAI +
              # X12_latent_leastAuk_WAI +
              #   X12_latent_zoop_EGoA +
              #   X12_latent_biomassAmphShelfSum +
              #   X12_latent_biomassEuphShelfSum +
              #   X12_latent_EGOA_DFA_lowerTrophic 
model_M2fish_preyAK_latentEGOA_preyNCC_predNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
                          X12_biomassMysidShelfSum_smoltyr +
                          X12_latent_crestedAuk_WAI +
                          X12_latent_leastAuk_WAI +
                          X12_latent_zoop_EGoA +
                          X12_latent_biomassAmphShelfSum +
                          X12_latent_biomassEuphShelfSum +
                          X12_latent_EGOA_DFA_lowerTrophic 
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  

            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ predNCC
                X16_DFA_sar_UCsprCh_smoltyr  ~X07_cpue_IntSprJunHW_smoltyr + 
                                               preyAK + 
                                               predAK 
              "

fit_M2fish_preyAK_latentEGOA_preyNCC_predNCC <- sem(
  model_M2fish_preyAK_latentEGOA_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

fit<-fit_M2fish_preyAK_latentEGOA_preyNCC_predNCC

sem_output_summary_fxn(fit)
parameterEstimates(fit, standardized = TRUE) %>%
  as_tibble() %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.lv, std.all) %>%
  print(n = Inf)

summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

#now akprey is not significant


# Log kept model--------
log_sem_model(fit_M2fish_preyAK_latentEGOA_preyNCC_predNCC,   "M2fish_preyAK_latentEGOA_preyNCC_predNCC",   "fish/shark predAK + latent EGOA preyAK + NCC prey + NCC pred, no growth")
#M2 var 2----------------
#sig: X13_DFA_WGOA_DFA_midTrophic_smoltyr, X14_DFA_pinkSalmonAsia_smoltyr,X14_pinkSalmonNorthAmerica_smoltyr

model_M2fish_preyAK_WGOAmid_preyNCC_predNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr

 
                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  
                           
                 preyAK =~X13_latent_EGOA_DFA_midTrophic +
                          X13_latent_WGOA_DFA_midTrophic +
                          X13_latent_HerringSSB_BCcoast_wcvi_cc_prd           
                           
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ predNCC
                X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr +
                                                X07_cpue_IntSprJunHW_smoltyr + 
                                               preyAK + 
                                               predAK
              "

fit_M2fish_preyAK_WGOAmid_preyNCC_predNCC <- sem(
  model_M2fish_preyAK_WGOAmid_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

fit<-fit_M2fish_preyAK_WGOAmid_preyNCC_predNCC

sem_output_summary_fxn(fit)
# Print full, untruncated parameter table (filtered to just latent loadings if desired)
parameterEstimates(fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  select(lhs, rhs, est, se, z, pvalue, std.lv, std.all)

summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

#hum, now AKprey is not significant
# --- Regression Path Ranking (By P-Value) ---
#   lhs op                          rhs    est    se      z pvalue
# X16_DFA_sar_UCsprCh_smoltyr  ~                       preyAK -0.238 0.223 -1.070  0.285
# X16_DFA_sar_UCsprCh_smoltyr  ~       X06_DFA_IGF_mu_smoltyr -0.509 0.319 -1.597  0.110
# X16_DFA_sar_UCsprCh_smoltyr  ~                       predAK -0.689 0.292 -2.360  0.018
# X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr  0.502 0.205  2.449  0.014
# X06_DFA_IGF_mu_smoltyr  ~                      preyNCC  0.927 0.149  6.239  0.000
# X07_cpue_IntSprJunHW_smoltyr  ~                      predNCC -0.904 0.156 -5.799  0.000



# Log kept model--------
log_sem_model(fit_M2fish_preyAK_WGOAmid_preyNCC_predNCC,   "M2fish_preyAK_WGOAmid_preyNCC_predNCC",   "fish/shark predAK + WGOAmid observed + NCC prey+ NCC pred, no growth")
print(model_comparison_table)

#M2 var 3----------------
#sig: X13_DFA_WGOA_DFA_midTrophic_smoltyr, X14_DFA_pinkSalmonAsia_smoltyr,X14_pinkSalmonNorthAmerica_smoltyr

model_M2fish_preyAK_pinkAsia_preyNCC_predNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr


                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  
 
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X14_DFA_pinkSalmonAsia_smoltyr + 
                                               predAK +
                                               X06_DFA_IGF_mu_smoltyr
              "

fit_M2fish_preyAK_pinkAsia_preyNCC_predNCC <- sem(
  model_M2fish_preyAK_pinkAsia_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_M2fish_preyAK_pinkAsia_preyNCC_predNCC)
summary(fit_M2fish_preyAK_pinkAsia_preyNCC_predNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


#pink NS if NCC prey are in the model -- very highly correlted. But because growth isn't in there, it makes more sense to keep pink

#M2 var 4----------------
model_M2fish_preyAK_pinkNA_preyNCC_predNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr


                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  
  
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X14_pinkSalmonNorthAmerica_smoltyr + 
                                               predAK +
                                               X06_DFA_IGF_mu_smoltyr
              "

fit_M2fish_preyAK_pinkNA_preyNCC_predNCC <- sem(
  model_M2fish_preyAK_pinkNA_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_M2fish_preyAK_pinkNA_preyNCC_predNCC)
summary(fit_M2fish_preyAK_pinkNA_preyNCC_predNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


#pink Asia and NA NS if NCC prey are in the model -- very highly correlated. But because growth isn't in there, it makes more sense to keep pink

#M3 Pink no AKpred -- no NCC prey or growth----------------
model_M3Pink_pinkNA_preyNCC_predNCC <- "
               preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr

                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  
  
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ predNCC
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X14_pinkSalmonNorthAmerica_smoltyr+
                                               X06_DFA_IGF_mu_smoltyr
              "

fit_model_M3Pink_pinkNA_preyNCC_predNCC <- sem(
  model_M3Pink_pinkNA_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)
fit<-fit_model_M3Pink_pinkNA_preyNCC_predNCC

sem_output_summary_fxn(fit)
summary(fit_M2fish_preyAK_pinkNA_preyNCC_predNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


#pink Asia and NA NS if NCC prey are in the model -- very highly correlated. But because growth isn't in there, it makes more sense to keep pink

#M2 reduced to NCC prey (no pink) ----------------
model_M2fish_preyNCC_predNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
 
                 predNCC =~ X09_latent_HakeAge5Plus +
                           X09_latent_pacificmackerel_GAM +
                           X09_latent_jackmackerel_GAM +
                           X09_latent_hakeFisheryWA  
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ predNCC
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               predAK +
                                               X06_DFA_IGF_mu_smoltyr
              "

fit_M2fish_preyNCC_predNCC <- sem(
  model_M2fish_preyNCC_predNCC,
  data = sem_data_merged, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_M2fish_preyNCC_predNCC)
summary(fit_M2fish_preyNCC_predNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
# Log kept model--------
log_sem_model(fit_M2fish_preyNCC_predNCC,   "M2fish_preyNCC_predNCC",   "fish/shark + NCC prey + NCC pred, no growth, no AK prey")
print(model_comparison_table %>% arrange(AIC))          

write.csv(model_comparison_table,file.path(data_out_dir,"latent_model_comparison_table.csv"),row.names=FALSE)


#SAVE RDATA FILE FOR MODEL FITS-------------
# Define the 13 model names from your summary table
model_names <- c(
  "1B",
  "M2fish_preyAK_WGOAmidtrophic",
  "M2fish_preyAK_pinkNA",
  "M2fish_preyAK_pinkAsia",
  "M2fish_preyAK_WGOAmid_preyNCC",
  "M2fish_preyNCC",
  "M1mammal_noPreyAK",
  "M2fish_preyAK_latentEGOA",
  "M1mammals_preyNCC",
  "M2fish_preyAK_WGOAmid_preyNCC_predNCC",
  "M2fish_preyNCC_predNCC",
  "M1mammals_preyNCC_predNCC",
  "M2fish_preyAK_latentEGOA_preyNCC_predNCC"
)

# Prepend "fit_" to each name
fit_names <- paste0("fit_", model_names)

# Keep only the objects currently present in your environment (warns if any are missing)
present_fits <- fit_names[sapply(fit_names, exists)]
missing_fits <- setdiff(fit_names, present_fits)

if (length(missing_fits) > 0) {
  warning("Not found in environment: ", paste(missing_fits, collapse = ", "))
}
#Not found in environment: fit_M2fish_preyAK_latentEGOA, fit_M2fish_preyAK_latentEGOA_preyNCC, fit_M2fish_preyAK_latentEGOA_preyNCC_predNCC 

# Save to your RData path (or swap in your preferred output filename)
save(
  list = present_fits,
  file = file.path(data_out_dir, "07f_latentmodel_fits.RData")
)
