ncc_candidates <- sem_data %>%
  select(matches("^X0[1-589]_|^X1[01]_"), -matches("_adultyr$")) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing")

print(as.data.frame(ncc_candidates), row.names = FALSE)

#Biological interpretation of fit_M8a_mammals_preyNCC4 the +0.802 covariance: Basin-scale North Pacific climate/regime forcing simultaneously creates favorable feeding conditions in the Northern California Current (preyNCC $\uparrow$ $\rightarrow$ Growth $\uparrow$) while coastwide marine mammal predator recovery/abundance (predAK $\uparrow$) exerts an offsetting top-down drag on survival (SAR $\downarrow$). Including that exogenous covariance explicitly captures the competing bottom-up vs. top-down forces across the coast!
  
# ==============================================================================
# ==============================================================================
# Script: 02_analysis/06k_fit_latent_models_NCCprey.R
# Step 14: Fit 4-Indicator preyNCC (X01, X02, X04, X05) in Both Leading Branches
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# Model M8a: Branch A (5-Mammal predAK) + 4-Indicator preyNCC
# ------------------------------------------------------------------------------
model_M8a_mammals_preyNCC4 <- "
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

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M8a_mammals_preyNCC4 <- sem(
  model_M8a_mammals_preyNCC4,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M8a_mammals_preyNCC4 (5 Mammals + 4-Ind preyNCC)\n")
cat("====================================================================\n")
summary(fit_M8a_mammals_preyNCC4, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
log_sem_model(fit_M8a_mammals_preyNCC4, "fit_M8a_mammals_preyNCC4", "Branch A (5-mammal predAK) + 4-ind preyNCC (X01, X02, X04, X05)")


# ------------------------------------------------------------------------------
# Model M8b: Branch B (2-Fish/Shark predAK + EGOA preyAK) + 4-Indicator preyNCC
# ------------------------------------------------------------------------------
# ==============================================================================
# Step 14b: Re-run Branch B Using the Winning 2-Indicator EGOA preyAK Factor
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

model_M8b_fish2_preyAK2_preyNCC4 <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr

  preyNCC =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
             X02_DFA_SeaNettle_smoltyr +
             X04_marketsquid_GAM_smoltyr +
             X05_DFA_abundSardine_smoltyr

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M8b_fish2_preyAK2_preyNCC4 <- sem(
  model_M8b_fish2_preyAK2_preyNCC4,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
)

cat("====================================================================\n")
cat("SUMMARY: fit_M8b_fish2_preyAK2_preyNCC4 (2-Ind preyAK + 4-Ind preyNCC)\n")
cat("====================================================================\n")
summary(fit_M8b_fish2_preyAK2_preyNCC4, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

log_sem_model(
  fit_M8b_fish2_preyAK2_preyNCC4,
  "fit_M8b_fish2_preyAK2_preyNCC4",
  "Branch B: 2-fish predAK + 2-ind EGOA preyAK + 4-ind preyNCC (no Growth->SAR)"
)