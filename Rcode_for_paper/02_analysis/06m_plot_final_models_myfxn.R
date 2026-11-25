# ==============================================================================
# Step 15 (Final): Plot Both Competing SEM Models via Lisa_sem_graph.v3
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(tidySEM)
  library(ggplot2)
})

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Shared Pretty Node Labels for Both Models
# ------------------------------------------------------------------------------
labels_final <- c(
  # preyNCC indicators (4)
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea Nettles",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_DFA_abundSardine_smoltyr"              = "Sardine\nAbund.",
  # predNCC single indicator
  "X09_DFA_HakeAge5Plus_smoltyr"              = "NCC Pred\n(Hake 5+)",
  # Mammal predAK indicators (5)
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"    = "Steller SL\n(WS)",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = "Sea Lions\n(EMB)",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = "Killer Whales\n(NR BC)",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = "Harbor Seal\n(WCVI)",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = "Steller SL\n(SEAK Pup)",
  # Fish/Shark predAK indicators (2)
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"   = "Sablefish\n(AK 3+)",
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z" = "Sleeper Shark\n(BSAI)",
  # preyAK indicators (2)
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr" = "EGOA Lower\nTrophic",
  "X12_biomassMysidShelfSum_smoltyr"      = "Shelf Mysids\n(EGOA)",
  # Latent factors
  "preyNCC" = "NCC Prey",
  "preyAK"  = "AK Prey\n(EGOA)",
  "predAK"  = "AK Pred",
  # Core endogenous nodes
  "X06_DFA_IGF_mu_smoltyr"       = "Early Marine\nGrowth",
  "X07_cpue_IntSprJunHW_smoltyr" = "Smolt\nCPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"  = "Chinook\nSAR"
)

# ------------------------------------------------------------------------------
# 2. BRANCH A: fit_M8a_mammals_preyNCC4 (Mammals + Early Growth + CPUE)
# ------------------------------------------------------------------------------
layout_M8a <- matrix(c(
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr", ""       , ""                      , ""                            , ""                           , ""                           , "X09_DFA_HakeAge5Plus_smoltyr",
  "X02_DFA_SeaNettle_smoltyr"                , "preyNCC", "X06_DFA_IGF_mu_smoltyr", ""                            , "X07_cpue_IntSprJunHW_smoltyr", ""                           , ""                            ,
  "X04_marketsquid_GAM_smoltyr"              , ""       , ""                      , "X16_DFA_sar_UCsprCh_smoltyr" , ""                            , ""                           , ""                            ,
  "X05_DFA_abundSardine_smoltyr"             , ""       , ""                      , ""                            , ""                            , ""                           , ""                            ,
  ""                                         , ""       , ""                      , ""                            , ""                            , "predAK"                     , "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z",
  ""                                         , ""       , ""                      , ""                            , ""                            , ""                           , "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z",
  ""                                         , ""       , ""                      , ""                            , ""                            , ""                           , "Killer_whales_NR_BC_adultyr__peff_Temp_z",
  ""                                         , ""       , ""                      , ""                            , ""                            , ""                           , "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z",
  ""                                         , ""       , ""                      , ""                            , ""                            , ""                           , "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"
), nrow = 9, ncol = 7, byrow = TRUE)

labels_M8a <- labels_final
labels_M8a["predAK"] <- "AK Pred\n(Mammals)"

plot_M8a <- Lisa_sem_graph.v3(
  fit               = fit_M8a_mammals_preyNCC4,
  title             = "Branch A: Coastal Marine Mammal Predation & Early Growth SEM",
  layout            = layout_M8a,
  node_labels       = labels_M8a,
  show_coefficients = TRUE,
  show_r2           = TRUE,
  r2_nodes          = c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"),
  width             = 15,
  height            = 10,
  save              = TRUE,
  savename          = file.path(fig_dir, "SEM_BranchA_Mammals.png"),
  display           = TRUE
)

# ------------------------------------------------------------------------------
# 3. BRANCH B: fit_M8b_fish2_preyAK2_preyNCC4 (Fish/Sharks + EGOA preyAK)
# ------------------------------------------------------------------------------
layout_M8b <- matrix(c(
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr", ""       , ""                      , ""                            , ""                            , ""      , "X09_DFA_HakeAge5Plus_smoltyr",
  "X02_DFA_SeaNettle_smoltyr"                , "preyNCC", "X06_DFA_IGF_mu_smoltyr", ""                            , "X07_cpue_IntSprJunHW_smoltyr", ""      , ""                            ,
  "X04_marketsquid_GAM_smoltyr"              , ""       , ""                      , "X16_DFA_sar_UCsprCh_smoltyr" , ""                            , ""      , ""                            ,
  "X05_DFA_abundSardine_smoltyr"             , ""       , ""                      , ""                            , ""                            , ""      , ""                            ,
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"    , "preyAK" , ""                      , ""                            , ""                            , "predAK", "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z",
  "X12_biomassMysidShelfSum_smoltyr"         , ""       , ""                      , ""                            , ""                            , ""      , "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z"
), nrow = 6, ncol = 7, byrow = TRUE)

labels_M8b <- labels_final
labels_M8b["predAK"] <- "AK Pred\n(Fish/Shark)"

plot_M8b <- Lisa_sem_graph.v3(
  fit               = fit_M8b_fish2_preyAK2_preyNCC4,
  title             = "Branch B: Subarctic Fish/Shark Predation & EGOA Prey SEM",
  layout            = layout_M8b,
  node_labels       = labels_M8b,
  show_coefficients = TRUE,
  show_r2           = TRUE,
  r2_nodes          = c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"),
  width             = 15,
  height            = 9,
  save              = TRUE,
  savename          = file.path(fig_dir, "SEM_BranchB_FishPreyAK.png"),
  display           = TRUE
)

cat("Saved both figures to:", fig_dir, "\n")