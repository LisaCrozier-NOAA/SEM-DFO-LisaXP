

# ==============================================================================
# Script: 02_analysis/07d_fit_latent_models_NCCprey.R
# Step 1: Fit current best latent models 
# ==============================================================================

#Add NCCprey latent to all current models
# Log all kept models------
# log_sem_model(fit_M1mammal_noPreyAK, "M1mammal_noPreyAK","5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE")
# log_sem_model(fit_M2fish_preyAK_EGOA,   "M2fish_preyAK_latentEGOA",   "2-fish/shark predAK + EGOA preyAK, no growth")
# log_sem_model(fit_M2fish_preyAK_WGOAmidtrophic, "M2fish_preyAK_WGOAmidtrophic", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
# log_sem_model(fit_M2fish_preyAK_pinkAsia, "M2fish_preyAK_pinkAsia", "2-fish/shark predAK + trimmed k3_group==2 (MidTrophic/Seabirds/PinkAsia)")
# log_sem_model(fit_M2fish_preyAK_pinkNA,   "M2fish_preyAK_pinkNA",   "2-fish/shark predAK + pink NorthAmer, no preyAK latent, no growth")

# ------------------------------------------------------------------------------
# Model M1a: Branch A (5-Mammal predAK) + 4-Indicator preyNCC------
# ------------------------------------------------------------------------------
          model_M1mammals_preyNCC <- "
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
          
          fit_M1mammals_preyNCC <- sem(
            model_M1mammals_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "NLMINB", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M1mammals_preyNCC)
          summary(fit_M1mammals_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          
          #log kept------
          log_sem_model(
            fit_M1mammals_preyNCC,
            model_name = "M1mammals_preyNCC",
            notes      = "5-mammal predAK (var floor 0.01), latent NCC, growth, no preyAK, fixed.x=TRUE"
          )

          # ------------------------------------------------------------------------------
          # M2 var 1: group 1: 2-Fish/Shark predAK + Multi-Indicator EGOA preyAK (k3_group == 1)----------
          # ------------------------------------------------------------------------------
          model_M2fish_preyAK_EGOA_preyNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyAK =~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
                          X12_biomassMysidShelfSum_smoltyr +
                          X12_egoa_krill_smoltyr 
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               preyAK + 
                                               predAK
              "
          
          fit_M2fish_preyAK_EGOA_preyNCC <- sem(
            model_M2fish_preyAK_EGOA_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M2fish_preyAK_EGOA_preyNCC)
          summary(fit_M2fish_preyAK_EGOA_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          
          
          
          # Log kept model--------
          log_sem_model(fit_M2fish_preyAK_EGOA_preyNCC,   "M2fish_preyAK_latentEGOA_preyNCC",   "fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth")
  #M2 var 2----------------
          #sig: X13_DFA_WGOA_DFA_midTrophic_smoltyr, X14_DFA_pinkSalmonAsia_smoltyr,X14_pinkSalmonNorthAmerica_smoltyr
          
          model_M2fish_preyAK_WGOAmid_preyNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X13_DFA_WGOA_DFA_midTrophic_smoltyr + 
                                               predAK
              "
          
          fit_M2fish_preyAK_WGOAmid_preyNCC <- sem(
            model_M2fish_preyAK_WGOAmid_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M2fish_preyAK_WGOAmid_preyNCC)
          summary(fit_M2fish_preyAK_WGOAmid_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          
          
          
          # Log kept model--------
          log_sem_model(fit_M2fish_preyAK_EGOA_preyNCC,   "M2fish_preyAK_latentEGOA_preyNCC",   "fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth")
          log_sem_model(fit_M2fish_preyAK_WGOAmid_preyNCC,   "M2fish_preyAK_WGOAmid_preyNCC",   "fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth")
          print(model_comparison_table)
          
          #M2 var 3----------------
          #sig: X13_DFA_WGOA_DFA_midTrophic_smoltyr, X14_DFA_pinkSalmonAsia_smoltyr,X14_pinkSalmonNorthAmerica_smoltyr
          
          model_M2fish_preyAK_pinkAsia_preyNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X14_DFA_pinkSalmonAsia_smoltyr + 
                                               predAK
              "
          
          fit_M2fish_preyAK_pinkAsia_preyNCC <- sem(
            model_M2fish_preyAK_pinkAsia_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M2fish_preyAK_pinkAsia_preyNCC)
          summary(fit_M2fish_preyAK_pinkAsia_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          
          
          #pink NS if NCC prey are in the model -- very highly correlted. But because growth isn't in there, it makes more sense to keep pink
          
          #M2 var 4----------------
          model_M2fish_preyAK_pinkNA_preyNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               X14_pinkSalmonNorthAmerica_smoltyr + 
                                               predAK
              "
          
          fit_M2fish_preyAK_pinkNA_preyNCC <- sem(
            model_M2fish_preyAK_pinkNA_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M2fish_preyAK_pinkNA_preyNCC)
          summary(fit_M2fish_preyAK_pinkNA_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          
          
          #pink Asia and NA NS if NCC prey are in the model -- very highly correlated. But because growth isn't in there, it makes more sense to keep pink
          
  #M2 reduced to NCC prey (no pink) ----------------
          model_M2fish_preyNCC <- "
                predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
                          DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                preyNCC=~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
                          X02_DFA_SeaNettle_smoltyr +
                          X04_marketsquid_GAM_smoltyr +
                          X05_DFA_abundSardine_smoltyr
 
             DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z
              
                X06_DFA_IGF_mu_smoltyr       ~ preyNCC
                X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
                X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                               predAK
              "
          
          fit_M2fish_preyNCC <- sem(
            model_M2fish_preyNCC,
            data = sem_data, std.lv = TRUE, missing = "ML",
            fixed.x = TRUE, bounds = "pos.var", optim.method = "BFGS", control = list(iter.max = 2000)
          )
          
          sem_output_summary_fxn(fit_M2fish_preyNCC)
          summary(fit_M2fish_preyNCC, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
          # Log kept model--------
          log_sem_model(fit_M2fish_preyNCC,   "M2fish_preyNCC",   "fish/shark predAK (fixed pos var) + EGOA preyAK + NCC prey, no growth")
print(model_comparison_table)          
