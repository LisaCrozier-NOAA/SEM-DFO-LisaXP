

# ==============================================================================
# Script: 02_analysis/08_plot_final_3_latent_models.R
# Fit & Plot Final 3 Alternative SEM Models with Locked Coordinates
# ==============================================================================
library(tidyverse)
library(lavaan)
library(tidySEM)
library(ggplot2)

proj_dir      <- file.path(getwd(), "Rcode_for_paper")
master_out    <- file.path(proj_dir, "Outputs")
data_out_dir  <- file.path(master_out, "data")
table_out_dir <- file.path(master_out, "tables")
plot_out_dir  <- file.path(master_out, "plots")

source("functions/Lisa_sem_graph_v4_noCI.fxn.r")

data_in_path    <- file.path(data_out_dir, "sem_master_peff_zscored_nccpred_latent.csv")
sem_data_merged <- read.csv(data_in_path, stringsAsFactors = FALSE)

# ==============================================================================
# 1. SHARED MASTER LABELS & INDICATOR LISTS
# ==============================================================================
# Verify the exact 4 X05_latent_* column names in your dataset
x05_latent_cols <- grep("^X05_latent_", names(sem_data_merged), value = TRUE)

preyNCC_indicators <- c(
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr",
  "X02_DFA_SeaNettle_smoltyr",
  "X04_marketsquid_GAM_smoltyr",
  x05_latent_cols
)

predNCC_indicators <- c(
  "X09_latent_HakeAge5Plus",
  "X09_latent_pacificmackerel_GAM",
  "X09_latent_jackmackerel_GAM",
  "X09_latent_hakeFisheryWA"
)

mammal_indicators <- c(
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"
)

fish_indicators <- c(
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z",
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z"
)

# Single unified label dictionary across all three models
labels_master <- c(
  # Latent & Core Structural Nodes
  "preyNCC"                            = "NCC Prey",
  "predNCC"                            = "NCC Fish\nPredators",
  "X14_pinkSalmonNorthAmerica_smoltyr" = "Pink Salmon\n(N. America)",
  "X06_DFA_IGF_mu_smoltyr"             = "Growth",
  "X07_cpue_IntSprJunHW_smoltyr"       = "CPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"        = "SAR",
  
  # preyNCC Indicators (update the 4 X05_latent_* keys to match your exact names)
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea\nNettle",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_latent_northAnch_NCC"                  = "Northern\nAnchovy",
  "X05_latent_sardine_NCC"                    = "Sardine\n(NCC)",
  "X05_latent_herring_GAM"                    = "Herring\n(GAM)",
  "X05_latent_sardine_GAM"                    = "Sardine\n(GAM)",
  
  # predNCC Indicators
  "X09_latent_HakeAge5Plus"        = "Hake\nAge 5+",
  "X09_latent_pacificmackerel_GAM" = "Pacific\nMackerel",
  "X09_latent_jackmackerel_GAM"    = "Jack\nMackerel",
  "X09_latent_hakeFisheryWA"       = "WA Hake\nFishery",
  
  # Mammal predAK Indicators
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"   = "Steller SL\n(WS)",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = "Sea Lions\n(EMB)",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = "NR Killer\nWhales",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = "Harbor Seal\n(WCVI)",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = "SEAK \nSSL",
  
  # Fish/Shark predAK Indicators
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"    = "Sablefish\n(3+)",
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z" = "Sleeper\nShark"
)

preyNCC_rhs <- paste(preyNCC_indicators, collapse = " + ")
predNCC_rhs <- paste(predNCC_indicators, collapse = " + ")

# ==============================================================================
# 2. MODEL 1: MAMMAL AK PREDATORS + DECOMPOSED NCC PREY & PREDATORS
# ==============================================================================
model_M1 <- paste0("
  predAK  =~ ", paste(mammal_indicators, collapse = " + "), "
  preyNCC =~ ", preyNCC_rhs, "
  predNCC =~ ", predNCC_rhs, "

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + X07_cpue_IntSprJunHW_smoltyr + predAK
")

fit_M1 <- sem(model_M1, data = sem_data_merged, std.lv = TRUE, missing = "ML",
              fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000))

layout_M1 <- build_master_layout_13(
  preyNCC_node = "preyNCC", predNCC_node = "predNCC", preyAK_node = "", predAK_node = "predAK",
  preyNCC_ind = preyNCC_indicators, predNCC_ind = predNCC_indicators, predAK_ind = mammal_indicators
)

gg_M1 <- Lisa_sem_graph_v4_noCI.fxn(
  fit               = fit_M1,
  layout            = layout_M1,
  node_labels       = c(labels_master, "predAK" = "Mammal\nPredators"),
  title             = "M1: Mammal Predators + Decomposed NCC Prey & Predators",
  save              = TRUE,
  savename          = file.path(plot_out_dir, "M1_mammals_preyNCC7_predNCC.png"),
  export_table_path = file.path(table_out_dir, "M1_mammals_preyNCC7_predNCC_params.csv")
)

# ==============================================================================
# 3. MODEL 2: FISH/SHARK AK PREDATORS + DECOMPOSED NCC PREY & PREDATORS
# ==============================================================================
model_M2 <- paste0("
  predAK  =~ ", paste(fish_indicators, collapse = " + "), "
  preyNCC =~ ", preyNCC_rhs, "
  predNCC =~ ", predNCC_rhs, "

  DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + X07_cpue_IntSprJunHW_smoltyr + predAK
")

fit_M2 <- sem(model_M2, data = sem_data_merged, std.lv = TRUE, missing = "ML",
              fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000))

layout_M2 <- build_master_layout_13(
  preyNCC_node = "preyNCC", predNCC_node = "predNCC", preyAK_node = "", predAK_node = "predAK",
  preyNCC_ind = preyNCC_indicators, predNCC_ind = predNCC_indicators, predAK_ind = fish_indicators
)

gg_M2 <- Lisa_sem_graph_v4_noCI.fxn(
  fit               = fit_M2,
  layout            = layout_M2,
  node_labels       = c(labels_master, "predAK" = "AK Fish\nPredators"),
  title             = "M2: AK Fish/Shark Predators + Decomposed NCC Prey & Predators",
  save              = TRUE,
  savename          = file.path(plot_out_dir, "M2_fish_preyNCC7_predNCC.png"),
  export_table_path = file.path(table_out_dir, "M2_fish_preyNCC7_predNCC_params.csv")
)

# ==============================================================================
# 4. MODEL 3: NORTH AMERICAN PINK SALMON + DECOMPOSED NCC PREY & PREDATORS
# ==============================================================================
model_M3 <- paste0("
  preyNCC =~ ", preyNCC_rhs, "
  predNCC =~ ", predNCC_rhs, "

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                X14_pinkSalmonNorthAmerica_smoltyr +
                                X06_DFA_IGF_mu_smoltyr
  # Explicit covariances between Pink Salmon and the two NCC latent factors
  X14_pinkSalmonNorthAmerica_smoltyr ~~ preyNCC + predNCC
")

fit_M3 <- sem(model_M3, data = sem_data_merged, std.lv = TRUE, missing = "ML",
              fixed.x = FALSE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000))
sem_output_summary_fxn(fit_M3)

layout_M3 <- build_master_layout_13(
  preyNCC_node = "preyNCC", predNCC_node = "predNCC",
  preyAK_node  = "X14_pinkSalmonNorthAmerica_smoltyr", predAK_node = "",
  growth_node  = "X06_DFA_IGF_mu_smoltyr", # Growth omitted in M3 so Pink Salmon effect is unconfounded
  preyNCC_ind  = preyNCC_indicators, predNCC_ind = predNCC_indicators
)

gg_M3 <- Lisa_sem_graph_v4_noCI.fxn(
  fit               = fit_M3,
  layout            = layout_M3,
   node_labels       = labels_master,
  title             = "M3: North American Pink Salmon + Decomposed NCC Prey & Predators",
  save              = TRUE,
  savename          = file.path(plot_out_dir, "M3_pinkNA_preyNCC7_predNCC.png"),
  export_table_path = file.path(table_out_dir, "M3_pinkNA_preyNCC7_predNCC_params.csv")
)
