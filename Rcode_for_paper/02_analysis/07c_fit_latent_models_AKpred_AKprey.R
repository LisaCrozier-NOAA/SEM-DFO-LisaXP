# ==============================================================================
# Script: 02_analysis/07c_fit_latent_models_AKpred_AKprey.R
# Step 10: Test Multi-Indicator predAK and then preyAK in Both Model Branches

#Take homes: 
#w/ mammals, it doesn't matter what you use for AKprey, they are never significant
#w/ fish, EGOA made a latent variable, but the others had convergence problems and were not coherent latents.
# entered 1 at a time, I kept the preyAK that had p<0.1 (WGOAmid, pink asia and pink NA)

# Log all kept models------
log_sem_model(fit_M1mammal_noPreyAK, "M1mammal_noPreyAK","5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE")
log_sem_model(fit_M2fish_preyAK_latentEGOA,   "M2fish_preyAK_latentEGOA",   "2-fish/shark predAK + EGOA preyAK, no growth")
log_sem_model(fit_M2fish_preyAK_pinkAsia, "M2fish_preyAK_pinkAsia", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
log_sem_model(fit_M2fish_preyAK_pinkNA,   "M2fish_preyAK_pinkNA",   "2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth")
log_sem_model(fit_M2fish_preyAK_WGOAmidtrophic, "M2fish_preyAK_WGOAmidtrophic", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")


suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# HELPER TO LOG & COMPARE MODELS ACROSS ITERATIONS
# Initialize model storage bank & comparison log (if not already in workspace)
if (!exists("sem_model_bank")) sem_model_bank <- list()
if (!exists("model_comparison_table")) model_comparison_table <- tibble()

log_sem_model <- function(fit_obj, model_name, notes = "") {
  sem_model_bank[[model_name]] <<- fit_obj
  r2 <- inspect(fit_obj, "r2")
  
  new_row <- tibble(
    Model     = model_name,
    npar      = as.integer(fitMeasures(fit_obj, "npar")),
    df        = as.integer(fitMeasures(fit_obj, "df")),
    Chisq     = round(as.numeric(fitMeasures(fit_obj, "chisq")), 2),
    p_Chisq   = round(as.numeric(fitMeasures(fit_obj, "pvalue")), 3),
    AIC       = round(as.numeric(fitMeasures(fit_obj, "aic")), 2),
    BIC       = round(as.numeric(fitMeasures(fit_obj, "bic")), 2),
    CFI       = round(as.numeric(fitMeasures(fit_obj, "cfi")), 3),
    SRMR      = round(as.numeric(fitMeasures(fit_obj, "srmr")), 3),
    Growth_R2 = round(as.numeric(r2["X06_DFA_IGF_mu_smoltyr"]), 3),
    CPUE_R2   = round(as.numeric(r2["X07_cpue_IntSprJunHW_smoltyr"]), 3),
    SAR_R2    = round(as.numeric(r2["X16_DFA_sar_UCsprCh_smoltyr"]), 3),
    Notes     = notes
  )
  
  if (nrow(model_comparison_table) > 0) {
    model_comparison_table <<- model_comparison_table %>%
      filter(Model != model_name) %>%
      bind_rows(new_row)
  } else {
    model_comparison_table <<- new_row
  }
  
  print(as.data.frame(model_comparison_table), row.names = FALSE)
}

# ==============================================================================
# Model M1 predAK models: 5 mammal cluster M1_predAK_mammal5_noPreyAK---------
M1mammal_noPreyAK <- "
              # Latent Measurement Model (predAK)
              predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
                        DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
                        Killer_whales_NR_BC_adultyr__peff_Temp_z +
                        harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
                        ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z
            
              # Fix boundary residual variance to a small positive floor (eliminates Heywood warning)
              DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z
            
              # Central Structural Backbone (without preyAK)
              X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
              X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
              X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                             X07_cpue_IntSprJunHW_smoltyr + 
                                             predAK
            "

fit_M1mammal_noPreyAK <- sem(
  M1mammal_noPreyAK,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_M1mammal_noPreyAK)
summary(fit_M1mammal_noPreyAK)
#log kept--------      
log_sem_model(fit_M1mammal_noPreyAK, "M1mammal_noPreyAK","5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE")

print(model_comparison_table)
# log_sem_model(fit_M2fish_preyAK_WGOAmidtrophic, "M2fish_preyAK_WGOAmidtrophic", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
# log_sem_model(fit_M2fish_preyAK_pinkAsia, "M2fish_preyAK_pinkAsia", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
# log_sem_model(fit_M2fish_preyAK_pinkNA,   "M2fish_preyAK_pinkNA",   "2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth")
# log_sem_model(fit_M2fish_preyAK_latentEGOA,   "M2fish_preyAK_latentEGOA",   "2-fish/shark predAK + EGOA preyAK, no growth")


# Model M2: AK offshore fish M2_predAK_fish2_noGrowth------
# Model M2: Pruned AK Offshore Fish/Shark Sub-cluster (2 indicators, no halibut (NS), no harbor seals (for ideological reasons))

M2_predAK_fish2_noGrowth <- "
                    predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                              DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z 
                  
                    X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
                    X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                    X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                                   X07_cpue_IntSprJunHW_smoltyr + 
                                                   X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                                   predAK
"

fit_M2_predAK_fish2 <- sem(
  M2_predAK_fish2_noGrowth,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_M2_predAK_fish2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
#not logged because this was replaced with a latent EGOA term   

# ------------------------------------------------------------------------------
# Model M2 akprey: group 1: 2-Fish/Shark predAK + Multi-Indicator EGOA preyAK (k3_group == 1)----------
# ------------------------------------------------------------------------------
model_M2fish_preyAK_EGOA <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            X12_biomassMysidShelfSum_smoltyr +
            X12_egoa_krill_smoltyr 

  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                 preyAK + 
                                 predAK
"

fit_M2fish_preyAK_latentEGOA <- sem(
  model_M2fish_preyAK_EGOA,
  data = sem_data, std.lv = TRUE, missing = "ML",
  fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
)

summary(fit_M2fish_preyAK_latentEGOA, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
sem_output_summary_fxn(fit_M2fish_preyAK_latentEGOA)



# Log kept model--------
log_sem_model(fit_M2fish_preyAK_latentEGOA,   "M2fish_preyAK_latentEGOA",   "2-fish/shark predAK + EGOA preyAK, no growth")


# ==============================================================================
# Model M2 akprey: group 2: 2-Fish/Shark predAK + WGOA Forage (LowerTrophic, NA Pinks, Capelin, Pollock) (k3_group == 2)----------
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Test k3_group == 3 (WGOA Forage: LowerTrophic, NA Pinks, Capelin, Pollock)
# ------------------------------------------------------------------------------
# 1a. Show 5-Mammal preyAK NS
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
sem_output_summary_fxn(fit_M6a_mammals_preyAK_k3g3)

# 1b.  2-Fish/Shark predAK ---- this preyAK group does not converge no matter how many you remove
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
          summary(fit_M6b_fish2_preyAK_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# 1b_trimmed. Final 2-Fish/Shark predAK ---- trim out pollock, then pink, then add them back individually
          #NS: X12_WGOA_DFA_lowerTrophic_smoltyr, X13_mid_il_capelin_smoltyr, X13_pollockGOA_age1and2_smoltyr at p=0.1
          #Sig close: X14_pinkSalmonNorthAmerica_smoltyr
            # preyAK =~ X12_WGOA_DFA_lowerTrophic_smoltyr +
            #           X14_pinkSalmonNorthAmerica_smoltyr +
            #           X13_mid_il_capelin_smoltyr +
            #           X13_pollockGOA_age1and2_smoltyr
        model_M2fish_preyAK_pinkNA <- "
          predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                    DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
        
        
          X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
          X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
          X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                         X14_pinkSalmonNorthAmerica_smoltyr + 
                                         predAK
        "
        
        fit_M2fish_preyAK_pinkNA <- sem(
          model_M2fish_preyAK_pinkNA,
          data = sem_data, std.lv = TRUE, missing = "ML",
          fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
        )
      #  summary(fit_M6b_fish2_preyAK_k3g3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
        sem_output_summary_fxn(fit_M2fish_preyAK_pinkNA)
        
        # Log kept model------
        log_sem_model(fit_M2fish_preyAK_pinkNA,   "M2fish_preyAK_pinkNA",   "2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth")
        #fit_M2fish_preyAK_latentEGOA
# ==============================================================================
# Model M2 akprey: group 3: 2-Fish/Shark predAK + (Mid-Trophic / Seabirds / Asian Pink Salmon Competitors) (k3_group == 3)----------
# ==============================================================================
# 2a. show that preyAK still NS w/ Mammal predAK
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
          
          sem_output_summary_fxn(fit_M7a_mammals_preyAK_k3g2)
          
          summary(fit_M7a_mammals_preyAK_k3g2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# 2b. With 2-Fish/Shark predAK -- doesn't converge, rm growth, rm deca
        model_M7b_fish2_preyAK_k3g2 <- "
          predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                    DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
        
          preyAK =~ X12_biomassDecaShelfSum_smoltyr +
                    X13_DFA_WGOA_DFA_midTrophic_smoltyr +
                    X13_EGOA_DFA_seabirds_smoltyr +
                    X14_DFA_pinkSalmonAsia_smoltyr
        
          X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
          X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
          X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                         preyAK + 
                                         predAK
        "
        
        fit_M7b_fish2_preyAK_k3g2 <- sem(
          model_M7b_fish2_preyAK_k3g2,
          data = sem_data, std.lv = TRUE, missing = "ML",
          fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
        )
        
        sem_output_summary_fxn(fit_M7b_fish2_preyAK_k3g2)
        summary(fit_M7b_fish2_preyAK_k3g2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)




# 2. Trimmed k3_group == NS preyAK, not stored , keep 1 at a time
  # preyAK =~ X13_DFA_WGOA_DFA_midTrophic_smoltyr +
  #           X13_EGOA_DFA_seabirds_smoltyr +
  #           X14_DFA_pinkSalmonAsia_smoltyr
        #NS: X13_EGOA_DFA_seabirds_smoltyr
        #sig: X13_DFA_WGOA_DFA_midTrophic_smoltyr, X14_DFA_pinkSalmonAsia_smoltyr
        
        M2fish_preyAK_pinkAsia <- "
        predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                  DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
      
      
      
        X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
        X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
        X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                       X14_DFA_pinkSalmonAsia_smoltyr +
                                       predAK
      "
      
      fit_M2fish_preyAK_pinkAsia <- sem(
        M2fish_preyAK_pinkAsia,
        data = sem_data, std.lv = TRUE, missing = "ML",
        fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
      )
      
      sem_output_summary_fxn(fit_M2fish_preyAK_pinkAsia)
      
      summary(fit_M2fish_preyAK_pinkAsia, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

      # Log kept model------
      log_sem_model(fit_M2fish_preyAK_pinkAsia, "M2fish_preyAK_pinkAsia", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
      log_sem_model(fit_M2fish_preyAK_pinkNA,   "M2fish_preyAK_pinkNA",   "2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth")
      
      M2fish_preyAK_WGOAmidtrophic <- "
        predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                  DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
      
      
      
        X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
        X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
        X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                       X13_DFA_WGOA_DFA_midTrophic_smoltyr +
                                       predAK
      "
      
      fit_M2fish_preyAK_WGOAmidtrophic <- sem(
        M2fish_preyAK_WGOAmidtrophic,
        data = sem_data, std.lv = TRUE, missing = "ML",
        fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
      )
      
      sem_output_summary_fxn(fit_M2fish_preyAK_WGOAmidtrophic)
      
      summary(M2fish_preyAK_WGOAmidtrophic, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
      
      # Log kept model------
      log_sem_model(fit_M2fish_preyAK_WGOAmidtrophic, "M2fish_preyAK_WGOAmidtrophic", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")


