# ==============================================================================
# Script: 02_analysis/06j2_plot_full_latent_models_semplot_stillwrong.R
# Step 15 (Fixed): Exact Name-Mapped Path Diagrams Saved to figures/
# ==============================================================================
suppressPackageStartupMessages({
  library(lavaan)
  library(semPlot)
})


# Dictionary mapping exact lavaan variable names -> clean plot labels
label_dict <- c(
  
  # Mammal predAK indicators
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"    = "Steller SL\n(WS)",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = "Sea Lions\n(EMB)",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = "Killer Whales\n(NR BC)",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = "Harbor Seal\n(WCVI)",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = "Steller SL\n(SEAK Pup)",
  
  # Fish/Shark predAK indicators
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"   = "Sablefish\n(AK 3+)",
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z" = "Sleeper Shark\n(BSAI)",
  
  # preyAK indicators
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr" = "EGOA Lower\nTrophic",
  "X12_biomassMysidShelfSum_smoltyr"      = "Shelf Mysids\n(EGOA)",
  
  # preyNCC indicators
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea Nettles",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_DFA_abundSardine_smoltyr"              = "Sardine\nAbund.",
  
  # Exogenous & Endogenous structural nodes
  "X09_DFA_HakeAge5Plus_smoltyr" = "Hake\nAge 5+",
  "X06_DFA_IGF_mu_smoltyr"       = "Early Growth\n(IGF-1)",
  "X07_cpue_IntSprJunHW_smoltyr" = "Juvenile\nCPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"  = "Chinook\nSAR",
  
  # Latents
  "predAK"  = "predAK",
  "preyAK"  = "preyAK\n(EGOA)",
  "preyNCC" = "preyNCC"
)

# Helper to build exact (x, y) layout & edge colors from semPlotModel
plot_clean_sem <- function(fit, title_text, predAK_label, xy_map) {
  pm <- semPlotModel(fit)
  node_names <- pm@Vars$name
  
  # 1. Exact label lookup by variable name
  node_labels <- unname(label_dict[node_names])
  node_labels[node_names == "predAK"] <- predAK_label
  
  # 2. Exact (x, y) coordinate matrix
  layout_mat <- t(sapply(node_names, function(v) xy_map[[v]]))
  
  # 3. Color edges by sign of standardized estimate
  stds <- standardizedSolution(fit)
  edge_df <- pm@Pars[pm@Pars$edge %in% c("~", "=~", "<->") & pm@Pars$lhs != pm@Pars$rhs, ]
  edge_cols <- ifelse(edge_df$std >= 0, "#1b7837", "#d95f02")
  
  semPaths(
    pm,
    what           = "std",
    whatLabels     = "std",
    style          = "ram",
    layout         = layout_mat,
    nodeLabels     = node_labels,
    sizeMan        = 7.5,
    sizeMan2       = 4.5,
    sizeLat        = 10.5,
    sizeLat2       = 6.5,
    label.cex      = 0.85,
    edge.label.cex = 0.85,
    edge.color     = edge_cols,
    residuals      = FALSE,
    intercepts     = FALSE,
    exoCov         = TRUE,
    curve          = 1.6,
    nCharNodes     = 0,
    mar            = c(3, 4, 4, 4)
  )
  title(title_text, line = 2.0, cex.main = 1.15)
}

# ------------------------------------------------------------------------------
# 1. Branch A Layout & Plot (Mammals + Growth + CPUE -> SAR)
# ------------------------------------------------------------------------------
xy_M8a <- list(
  
  # Column 1 (x = -1.0): Indicators & Hake
  "X09_DFA_HakeAge5Plus_smoltyr"                    = c(-1.0,  0.95),
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr"       = c(-1.0,  0.65),
  "X02_DFA_SeaNettle_smoltyr"                       = c(-1.0,  0.45),
  "X04_marketsquid_GAM_smoltyr"                     = c(-1.0,  0.25),
  "X05_DFA_abundSardine_smoltyr"                    = c(-1.0,  0.05),
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"    = c(-1.0, -0.25),
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = c(-1.0, -0.43),
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = c(-1.0, -0.61),
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = c(-1.0, -0.79),
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = c(-1.0, -0.97),
  
  # Column 2 (x = -0.25): Latent Factors
  "preyNCC"                                         = c(-0.25,  0.35),
  "predAK"                                          = c(-0.25, -0.60),
  
  # Column 3 (x = 0.40): Intermediate Endogenous Nodes
  "X07_cpue_IntSprJunHW_smoltyr"                    = c( 0.40,  0.85),
  "X06_DFA_IGF_mu_smoltyr"                          = c( 0.40,  0.35),
  
  # Column 4 (x = 1.0): Final Response
  "X16_DFA_sar_UCsprCh_smoltyr"                     = c( 1.00,  0.00)
)

png("figures/SEM_BranchA_Mammals.png", width = 3200, height = 2100, res = 300)
plot_clean_sem(fit_M8a_mammals_preyNCC4, "Branch A: Coastal Marine Mammal Predation & Early Growth SEM", "predAK\n(Mammals)", xy_M8a)
dev.off()
plot_clean_sem(fit_M8a_mammals_preyNCC4, "Branch A: Coastal Marine Mammal Predation & Early Growth SEM", "predAK\n(Mammals)", xy_M8a)

# ------------------------------------------------------------------------------
# 2. Branch B Layout & Plot (Fish/Sharks + EGOA preyAK + CPUE -> SAR)
# ------------------------------------------------------------------------------
xy_M8b <- list(
  # Column 1 (x = -1.0): Indicators & Hake
  "X09_DFA_HakeAge5Plus_smoltyr"                          = c(-1.0,  0.95),
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr"             = c(-1.0,  0.68),
  "X02_DFA_SeaNettle_smoltyr"                             = c(-1.0,  0.50),
  "X04_marketsquid_GAM_smoltyr"                           = c(-1.0,  0.32),
  "X05_DFA_abundSardine_smoltyr"                          = c(-1.0,  0.14),
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"                 = c(-1.0, -0.15),
  "X12_biomassMysidShelfSum_smoltyr"                      = c(-1.0, -0.35),
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"    = c(-1.0, -0.65),
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z" = c(-1.0, -0.85),
  # Column 2 (x = -0.25): Latent Factors
  "preyNCC"                                               = c(-0.25,  0.40),
  "preyAK"                                                = c(-0.25, -0.25),
  "predAK"                                                = c(-0.25, -0.75),
  # Column 3 (x = 0.40): Intermediate Endogenous Nodes
  "X07_cpue_IntSprJunHW_smoltyr"                          = c( 0.40,  0.85),
  "X06_DFA_IGF_mu_smoltyr"                                = c( 0.40,  0.40),
  # Column 4 (x = 1.0): Final Response
  "X16_DFA_sar_UCsprCh_smoltyr"                           = c( 1.00,  0.00)
)

png("figures/SEM_BranchB_FishPreyAK.png", width = 3200, height = 2100, res = 300)
plot_clean_sem(fit_M8b_fish2_preyAK2_preyNCC4, "Branch B: Subarctic Fish/Shark Predation & EGOA Prey SEM", "predAK\n(Fish/Shark)", xy_M8b)
dev.off()
plot_clean_sem(fit_M8b_fish2_preyAK2_preyNCC4, "Branch B: Subarctic Fish/Shark Predation & EGOA Prey SEM", "predAK\n(Fish/Shark)", xy_M8b)

cat("Saved cleanly to:", file.path(getwd(), "figures"), "\n")