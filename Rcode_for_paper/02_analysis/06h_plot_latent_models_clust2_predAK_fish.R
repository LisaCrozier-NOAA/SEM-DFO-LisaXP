
# ==============================================================================
# Script: 02_analysis/06h_plot_latent_models_clust2_predAK_fish.R
#       see 06a_cluster_diagrams_semnodes_v2.R for cluster diagram
# ==============================================================================
#are growth and akprey correlated? not really. it is coincidental with sleepers
      # print(round(cor(cor_vars, use = "complete.obs"), 3))

      #           Growth preyAK sablefish sleeper harbor_WS
      # Growth     1.000 -0.370    -0.503  -0.783    -0.520
      # preyAK    -0.370  1.000     0.123   0.029    -0.071
      # sablefish -0.503  0.123     1.000   0.787     0.333
      # sleeper   -0.783  0.029     0.787   1.000     0.588
      # harbor_WS -0.520 -0.071     0.333   0.588     1.000
# ==============================================================================
# Step 7: Test Alternative predAK Cluster (k2_group == 2: Offshore Fish / Sharks)
# ==============================================================================

# ------------------------------------------------------------------------------
# Model M2a: Full k2_group == 2 Cluster (5 indicators from dendrogram Cluster 2)
# ------------------------------------------------------------------------------
model_M2a_predAK_k2grp2 <- "
  predAK =~ DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z +
            sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_halibutBiomassAge8plus_predAK_adultyr__peff_TempPrey1_z +
            harbour_seal_DFO_Coastwide_smoltyr__peff_TempPrey1_z +
            Harbour_s_WS_adultyr__peff_Base_z

  # Central Structural Backbone (same structure as M1 for direct comparison)
  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M2a_predAK_k2grp2 <- sem(
  model_M2a_predAK_k2grp2,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M2a_predAK_k2grp2 (Full k2_group == 2 Cluster)\n")
cat("====================================================================\n")
summary(fit_M2a_predAK_k2grp2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# Model M2b: Pure AK Offshore Fish/Shark Sub-cluster (3 indicators, no harbor seals)
# ------------------------------------------------------------------------------
model_M2b_predAK_fish3 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z +
            DFA_halibutBiomassAge8plus_predAK_adultyr__peff_TempPrey1_z

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M2b_predAK_fish3 <- sem(
  model_M2b_predAK_fish3,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M2b_predAK_fish3 (3 AK Offshore Fish/Shark Indicators)\n")
cat("====================================================================\n")
summary(fit_M2b_predAK_fish3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# Model M2c: Pruned AK Offshore Fish/Shark Sub-cluster (2 indicators, no halibut)
# ------------------------------------------------------------------------------
model_M2c_predAK_fish3 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                 predAK
"

fit_M2c_predAK_fish3 <- sem(
  model_M2c_predAK_fish3,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

cat("\n====================================================================\n")
cat("SUMMARY: fit_M2b_predAK_fish3 (3 AK Offshore Fish/Shark Indicators)\n")
cat("====================================================================\n")
summary(fit_M2c_predAK_fish3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# Log both models into the persistent comparison table
# ------------------------------------------------------------------------------
log_sem_model(
  fit_M2a_predAK_k2grp2,
  model_name = "fit_M2a_predAK_k2grp2",
  notes      = "Dendrogram k2_group==2 (5 indicators), no preyAK, fixed.x=TRUE"
)

log_sem_model(
  fit_M2b_predAK_fish3,
  model_name = "fit_M2b_predAK_fish3",
  notes      = "3 AK offshore fish/shark indicators, no preyAK, fixed.x=TRUE"
)

#final models----------
# 1. Exact bivariate correlations among Growth, preyAK, and the fish/shark predators
cor_vars <- sem_data %>%
  select(
    Growth     = X06_DFA_IGF_mu_smoltyr,
    preyAK     = X12_DFA_EGOA_DFA_lowerTrophic_smoltyr,
    sablefish  = sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z,
    sleeper    = DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z,
    harbor_WS  = Harbour_s_WS_adultyr__peff_Base_z
  )

cat("Bivariate Correlation Matrix (Growth, preyAK, and Fish/Shark/WS Predators):\n")
print(round(cor(cor_vars, use = "complete.obs"), 3))
#           Growth preyAK sablefish sleeper harbor_WS
# Growth     1.000 -0.370    -0.503  -0.783    -0.520
# preyAK    -0.370  1.000     0.123   0.029    -0.071
# sablefish -0.503  0.123     1.000   0.787     0.333
# sleeper   -0.783  0.029     0.787   1.000     0.588
# harbor_WS -0.520 -0.071     0.333   0.588     1.000

# 2. Fit & Log M2a_trimmed (sleepers + sablefish + WS harbor seals, Growth, no preyAK)
model_M2a_trimmed <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z +
            Harbour_s_WS_adultyr__peff_Base_z

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M2a_trimmed <- sem(
  model_M2a_trimmed, data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

# 3. Fit & Log M2c (sleepers + sablefish + preyAK, no Growth->SAR)
model_M2c_fish2_preyAK <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                 predAK
"

fit_M2c_fish2_preyAK <- sem(
  model_M2c_fish2_preyAK, data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

log_sem_model(fit_M2a_trimmed, "fit_M2a_trimmed", "3-ind predAK (sablefish, sleeper, WS seal), no preyAK")
log_sem_model(fit_M2c_fish2_preyAK, "fit_M2c_fish2_preyAK", "2-ind predAK (sablefish, sleeper) + preyAK, no Growth->SAR")
