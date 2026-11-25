# ==============================================================================
# Script: 02_analysis/06j_fit_latent_models_AKprey.R
# Step 10: Test Multi-Indicator preyAK (k3_group == 1) in Both Model Branches

#store history of how eliminated  AKprey
#1. cluster analysis showed 3 groups: 
# EGOA lowertrophic+mysid+krill
#        mammals -- prey_latent NA, predAK ~~ preyAK; krill loading= 0.496 within latent, removing krill did not help latent
              # Covariances:
              # Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
              # predAK ~~                                                             
              #   preyAK           -0.748    0.121   -6.206    0.000   -0.748   -0.748

#       sleeper/sablefish -- all 3 load high on latent, growth highly NS, remove growth

# WGOA lowertrophic + pinkNA
# WGOA midtrohpic+ deca + pinkAsia

# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# Model M4a: 5-Mammal predAK + Multi-Indicator EGOA preyAK (k3_group == 1)
# ------------------------------------------------------------------------------
model_M4a_mammals_preyAK5 <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr +
            X12_egoa_krill_smoltyr 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M4a_mammals_preyAK5 <- sem(
  model_M4a_mammals_preyAK5,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M4a_mammals_preyAK5 (5 Mammals + k3_group==1 preyAK)\n")
cat("====================================================================\n")
summary(fit_M4a_mammals_preyAK5, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# Model M4b: 2-Fish/Shark predAK + Multi-Indicator EGOA preyAK (k3_group == 1)
# ------------------------------------------------------------------------------
model_M4b_fish2_preyAK5 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr +
            X12_egoa_krill_smoltyr 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M4b_fish2_preyAK5 <- sem(
  model_M4b_fish2_preyAK5,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M4b_fish2_preyAK5 (2 Fish/Sharks + k3_group==1 preyAK)\n")
cat("====================================================================\n")
summary(fit_M4b_fish2_preyAK5, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


#Trimmed======================
# ==============================================================================
# Step 11: Trim Low-Loading (< 0.50) & adultyr Indicators from preyAK
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# Model M5a: 5-Mammal predAK + Trimmed 2-Indicator EGOA preyAK (smoltyr only)
# ------------------------------------------------------------------------------
model_M5a_mammals_preyAK2 <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M5a_mammals_preyAK2 <- sem(
  model_M5a_mammals_preyAK2,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M5a_mammals_preyAK2 (5 Mammals + Trimmed 2-Ind preyAK)\n")
cat("====================================================================\n")
summary(fit_M5a_mammals_preyAK2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# Model M5b: 2-Fish/Shark predAK + Trimmed 2-Indicator EGOA preyAK (with Growth)
# ------------------------------------------------------------------------------
model_M5b_fish2_preyAK2 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr  +
            X12_egoa_krill_smoltyr 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M5b_fish2_preyAK2 <- sem(
  model_M5b_fish2_preyAK2,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M5b_fish2_preyAK2 (2 Fish/Sharks + Trimmed 2-Ind preyAK)\n")
cat("====================================================================\n")
summary(fit_M5b_fish2_preyAK2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Log both trimmed models
log_sem_model(fit_M5a_mammals_preyAK2, "fit_M5a_mammals_preyAK2", "5-mammal predAK + trimmed 2-ind EGOA preyAK (Std.all >= 0.5)")
log_sem_model(fit_M5b_fish2_preyAK2,   "fit_M5b_fish2_preyAK2",   "2-fish/shark predAK + trimmed 2-ind EGOA preyAK (Std.all >= 0.5)")

# 1. Fit & log clean Branch B (2-fish/shark predAK + 2-ind EGOA preyAK, no Growth->SAR)
model_M5c_fish2_preyAK2_noGrowth <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr  +
            X12_egoa_krill_smoltyr 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M5c_fish2_preyAK2_noGrowth <- sem(
  model_M5c_fish2_preyAK2_noGrowth,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

summary(fit_M5c_fish2_preyAK2_noGrowth, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(
  fit_M5c_fish2_preyAK2_noGrowth,
  "fit_M5c_fish2_preyAK2_noGrowth",
  "Branch B: 2-fish/shark predAK + 2-ind EGOA preyAK + CPUE (no Growth->SAR)"
)

# 2. Check a priori prefixes for preyNCC & predNCC (smoltyr only, 0 missing years)
other_blocks <- sem_data %>%
  select(matches("^X0[1-589]_|^X1[01]_"), -matches("_adultyr$")) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing")

cat("\nA priori candidates for preyNCC / predNCC (smoltyr only):\n")
print(as.data.frame(other_blocks), row.names = FALSE)


#next set of akprey clusters
# ==============================================================================
# Step 12: Evaluate Remaining smoltyr preyAK Clusters (k3_group == 2 & 3)-------------
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# 1. Test k3_group == 3 (WGOA Forage: LowerTrophic, NA Pinks, Capelin, Pollock)
# ------------------------------------------------------------------------------
# 1a. With 5-Mammal predAK
model_M6a_mammals_preyAK_k3g3 <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X12_WGOA_DFA_lowerTrophic_smoltyr +
            X14_pinkSalmonNorthAmerica_smoltyr +
            X13_mid_il_capelin_smoltyr +
            X13_pollockGOA_age1and2_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M6a_mammals_preyAK_k3g3 <- sem(
  model_M6a_mammals_preyAK_k3g3,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M6a_mammals_preyAK_k3g3 (5 Mammals + k3_group == 3 WGOA)\n")
cat("====================================================================\n")
summary(fit_M6a_mammals_preyAK_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# 1b. With 2-Fish/Shark predAK
model_M6b_fish2_preyAK_k3g3 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_WGOA_DFA_lowerTrophic_smoltyr +
            X14_pinkSalmonNorthAmerica_smoltyr +
            X13_mid_il_capelin_smoltyr +
            X13_pollockGOA_age1and2_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M6b_fish2_preyAK_k3g3 <- sem(
  model_M6b_fish2_preyAK_k3g3,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M6b_fish2_preyAK_k3g3 (2 Fish/Sharks + k3_group == 3 WGOA)\n")
cat("====================================================================\n")
summary(fit_M6b_fish2_preyAK_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# 2. Test k3_group == 2 (Mid-Trophic / Seabirds / Asian Pink Salmon Competitors)
# ------------------------------------------------------------------------------
# 2a. With 5-Mammal predAK
model_M7a_mammals_preyAK_k3g2 <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X12_biomassDecaShelfSum_smoltyr +
            X13_DFA_WGOA_DFA_midTrophic_smoltyr +
            X13_EGOA_DFA_seabirds_smoltyr +
            X14_DFA_pinkSalmonAsia_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M7a_mammals_preyAK_k3g2 <- sem(
  model_M7a_mammals_preyAK_k3g2,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M7a_mammals_preyAK_k3g2 (5 Mammals + k3_group == 2 Comp/Mid)\n")
cat("====================================================================\n")
summary(fit_M7a_mammals_preyAK_k3g2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# 2b. With 2-Fish/Shark predAK
model_M7b_fish2_preyAK_k3g2 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_biomassDecaShelfSum_smoltyr +
            X13_DFA_WGOA_DFA_midTrophic_smoltyr +
            X13_EGOA_DFA_seabirds_smoltyr +
            X14_DFA_pinkSalmonAsia_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M7b_fish2_preyAK_k3g2 <- sem(
  model_M7b_fish2_preyAK_k3g2,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M7b_fish2_preyAK_k3g2 (2 Fish/Sharks + k3_group == 2 Comp/Mid)\n")
cat("====================================================================\n")
summary(fit_M7b_fish2_preyAK_k3g2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Log all four models into the comparison table
log_sem_model(fit_M6a_mammals_preyAK_k3g3, "fit_M6a_mammals_preyAK_k3g3", "5-mammal predAK + k3_group==3 (WGOA forage)")
log_sem_model(fit_M6b_fish2_preyAK_k3g3,   "fit_M6b_fish2_preyAK_k3g3",   "2-fish/shark predAK + k3_group==3 (WGOA forage)")
log_sem_model(fit_M7a_mammals_preyAK_k3g2, "fit_M7a_mammals_preyAK_k3g2", "5-mammal predAK + k3_group==2 (MidTrophic/PinkAsia)")
log_sem_model(fit_M7b_fish2_preyAK_k3g2,   "fit_M7b_fish2_preyAK_k3g2",   "2-fish/shark predAK + k3_group==2 (MidTrophic/PinkAsia)")


#trimmed k3_group--------------------
# ==============================================================================
# Step 13: Clean Trimmed Test of k3_group == 2 (MidTrophic / Seabirds / PinkAsia)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# 1. Trimmed k3_group == 2 (decapods dropped, p = 0.113) + 5-Mammal predAK
model_M7a_trimmed <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X13_DFA_WGOA_DFA_midTrophic_smoltyr +
            X14_DFA_pinkSalmonAsia_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M7a_trimmed <- sem(
  model_M7a_trimmed,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M7a_trimmed (5 Mammals + Trimmed 3-Ind k3_group == 2)\n")
cat("====================================================================\n")
summary(fit_M7a_trimmed, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(fit_M7a_trimmed, "fit_M7a_trimmed", "5-mammal predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")

# 2. Trimmed k3_group == 2 + 2-Fish/Shark predAK (using NLMINB for clean convergence)
model_M7b_trimmed <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X13_DFA_WGOA_DFA_midTrophic_smoltyr +
            X13_EGOA_DFA_seabirds_smoltyr +
            X14_DFA_pinkSalmonAsia_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M7b_trimmed <- sem(
  model_M7b_trimmed,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M7b_trimmed (2 Fish/Sharks + Trimmed 3-Ind k3_group == 2)\n")
cat("====================================================================\n")
summary(fit_M7b_trimmed, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(fit_M7b_trimmed, "fit_M7b_trimmed", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")


#cluster 3======================
# ==============================================================================
# Step 13b: Test Cleaned smoltyr-Only Cluster 3 (k3_group == 3: WGOA Forage)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# 1. Clean smoltyr-only k3_group == 3 + 5-Mammal predAK
model_M6a_smoltyr_k3g3 <- "
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  preyAK =~ X14_pinkSalmonNorthAmerica_smoltyr +
            X13_pollockGOA_age1and2_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M6a_smoltyr_k3g3 <- sem(
  model_M6a_smoltyr_k3g3,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M6a_smoltyr_k3g3 (5 Mammals + smoltyr-only k3_group == 3)\n")
cat("====================================================================\n")
summary(fit_M6a_smoltyr_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(fit_M6a_smoltyr_k3g3, "fit_M6a_smoltyr_k3g3", "5-mammal predAK + smoltyr-only k3_group==3 (WGOA forage)")

# 2. Clean smoltyr-only k3_group == 3 + 2-Fish/Shark predAK
model_M6b_smoltyr_k3g3 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ #X14_pinkSalmonNorthAmerica_smoltyr +
            X13_mid_il_capelin_smoltyr +
            X13_pollockGOA_age1and2_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M6b_smoltyr_k3g3 <- sem(
  model_M6b_smoltyr_k3g3,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M6b_smoltyr_k3g3 (2 Fish/Sharks + smoltyr-only k3_group == 3)\n")
cat("====================================================================\n")
summary(fit_M6b_smoltyr_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(fit_M6b_smoltyr_k3g3, "fit_M6b_smoltyr_k3g3", "2-fish/shark predAK + smoltyr-only k3_group==3 (WGOA forage)")

print(as.data.frame(model_comparison_table), row.names = FALSE)
