

# Script: 02_analysis/06c_build_sem_latent_models.R




#Selected models:
#Results for NCC PREY-------
        # Latent       Indicator    Est Std_Loading    SE       z pval
        # Latent_sardine      NCCsardine  1.000       0.937 0.000      NA   NA
        # Latent_sardine  JSOES_plankton -1.048      -0.992 0.086 -12.238    0
        # Latent_sardine marketsquid_GAM -1.041      -0.976 0.094 -11.108    0
        # Latent_sardine       SeaNettle  0.676       0.633 0.179   3.783    0
        
        # ==============================================
        #   MODEL SUMMARY REPORT               
        # ==============================================
        #   
        #   --- Model Fit Metrics ---
        #   pvalue     cfi     aic   rmsea 
        # 0.005   0.934 187.809   0.315 
        # 
        # --- R-Squared (Endogenous Variables) ---
        #   No target R2 variables found in this model.
        # 
        # --- Regression Path Ranking (By P-Value) ---
        #   lhs op            rhs    est    se      z pvalue
        # Growth  ~ Latent_sardine -1.007 0.106 -9.505      0
        # 
        # ==============================================

# Results for AKprey -- just pink
        # 
        # ==============================================
        #   MODEL SUMMARY REPORT               
        # ==============================================
        #   
        #   --- Model Fit Metrics ---
        #   pvalue    cfi    aic  rmsea 
        # NA  1.000 59.387  0.000 
        # 
        # --- R-Squared (Endogenous Variables) ---
        #   No target R2 variables found in this model.
        # 
        # --- Regression Path Ranking (By P-Value) ---
        #   lhs op                                rhs    est    se      z pvalue
        # SAR  ~ X14_pinkSalmonNorthAmerica_smoltyr -0.315 0.148 -2.128  0.033
        # SAR  ~                               CPUE  0.647 0.148  4.379  0.000
        # 
        #==============================================
#Alternative AKprey----
        #   Latent Variables:
        #   Estimate  Std.Err  z-value  P(>|z|)
        # my_latent =~                                        
        #   X12_bmssDcShS_    1.000                           
        # X13_EGOA_DFA__    0.708    0.382    1.856    0.063
        # X12_DFA_EGOA_D   -0.919    0.424   -2.169    0.030
        # X12_bmssMysSS_   -0.899    0.471   -1.907    0.057
        # 
        # Regressions:
        #   Estimate  Std.Err  z-value  P(>|z|)
        # SAR ~                                               
        #   CPUE              0.418    0.186    2.248    0.025
        # Growth            0.643    0.167    3.859    0.000
        # my_latent        -1.373    0.503   -2.729    0.006
        
        # ==============================================
        #   MODEL SUMMARY REPORT               
        # ==============================================
        #   
        #   --- Model Fit Metrics ---
        #   pvalue     cfi     aic   rmsea 
        # 0.000   0.475 330.479   0.337 
        # 
        # --- R-Squared (Endogenous Variables) ---
        #   No target R2 variables found in this model.
        # 
        # --- Regression Path Ranking (By P-Value) ---
        #   lhs op       rhs    est    se      z pvalue
        # SAR  ~      CPUE  0.418 0.186  2.248  0.025
        # SAR  ~ my_latent -1.373 0.503 -2.729  0.006
        # SAR  ~    Growth  0.643 0.167  3.859  0.000
        # 
        # ==============================================
        #   
#AK Pred---------
# sem_latent <- "
#           my_latent =~ Peff_sablefishAK_biomass3plus_smoltyr + 
#                           Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
#                             Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
#                             Peff_Harbour_s_WS_adultyr +
#                             Peff_DFA_halibutBiomassAge8plus_predAK_adultyr
#           SAR  ~  CPUE + Growth + my_latent
#         
#         "

#pruned growth, now everything signif but some var are neg
# Latent Variables:
#   Estimate  Std.Err  z-value  P(>|z|)
# my_latent =~                                        
#   Pff_sblfsAK_3_    1.000                           
# P_DFA_SBSAI_AK    0.936    0.235    3.987    0.000
# P_DFA_SBSAI_AK    1.271    0.213    5.965    0.000
# Pff_Hrbr_s_WS_    0.747    0.248    3.008    0.003
# Pf_DFA_BA8_AK_    0.582    0.263    2.208    0.027
# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue     cfi     aic   rmsea 
# 0.000   0.725 361.050   0.280 
# 
# --- R-Squared (Endogenous Variables) ---
#   No target R2 variables found in this model.
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op       rhs    est    se      z pvalue
# SAR  ~ my_latent -0.526 0.201 -2.613  0.009
# SAR  ~      CPUE  0.549 0.140  3.909  0.000
# 
# ==============================================

  
  
  
  
  
  
#----------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

raw_out_path    <- file.path(data_dir, "sem_modeling_data_screened_raw_1998_2021.csv")
scaled_out_path <- file.path(data_dir, "sem_modeling_data_screened_zscored_1998_2021.csv")

write_csv(df_master, raw_out_path)
write_csv(df_master_scaled, scaled_out_path)

# Models--------
# NCC PREY-------

sardine_latent <- "
  Latent_sardine =~  
                   X01_DFA_sumPreyOfPrey_planktonJun_smoltyr + 
                   X02_DFA_SeaNettle_smoltyr +
                   X04_marketsquid_GAM_smoltyr + 
                   X05_DFA_abundSardine_smoltyr
  Growth ~   Latent_sardine

"


fit_sardine <- sem(
  sardine_latent,
  data         = df_master_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

sem_output_summary_fxn(fit_sardine)
summary(fit_sardine)

parameterEstimates(fit_sardine, standardized = TRUE) %>%
  filter(op == "=~") %>%
  select(Latent = lhs, Indicator = rhs, Est = est, 
         Std_Loading = std.all, SE = se, z = z, pval = pvalue) %>%
  print(row.names = FALSE)


#NCC PRED not nec (only hake)-------
#AK PREY-------
names(df_master_scaled)

        sem_latent <- "
          my_latent =~ X14_pinkSalmonNorthAmerica_smoltyr + 
                          X13_pollockGOA_age1and2_smoltyr + 
                            X12_WGOA_DFA_lowerTrophic_smoltyr 
          SAR  ~  CPUE + Growth + my_latent
        
        "
        
        sem_latent <- "
          my_latent =~ X14_pinkSalmonNorthAmerica_smoltyr + 
                          X13_pollockGOA_age1and2_smoltyr 
          SAR  ~  CPUE + Growth + my_latent
        
        "
        sem_latent <- "
         SAR  ~  CPUE + Growth + X14_pinkSalmonNorthAmerica_smoltyr
        
        "
        sem_latent <- "
         SAR  ~  CPUE + X14_pinkSalmonNorthAmerica_smoltyr
        
        "
       #full latent NS, prune WGOA prune pollock (leaves pink observed), prune growth    -----------
        
        
        #Start 2n cluster-------
        sem_latent <- "
          my_latent =~ X12_biomassDecaShelfSum_smoltyr + 
                          X13_DFA_WGOA_DFA_midTrophic_smoltyr + 
                            X13_EGOA_DFA_seabirds_smoltyr +
                            X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
                            X12_biomassMysidShelfSum_smoltyr
          SAR  ~  CPUE + Growth + my_latent
        
        "
        sem_latent <- "
          my_latent =~ X12_biomassDecaShelfSum_smoltyr + 
                            X13_EGOA_DFA_seabirds_smoltyr +
                            X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
                            X12_biomassMysidShelfSum_smoltyr
          SAR  ~  CPUE + Growth + my_latent
        
        "
        
        
        fit <- sem(
          sem_latent,
          data         = df_master_scaled,
          missing      = "ML",
          bounds       = "pos.var",
          optim.method = "BFGS",
          control      = list(iter.max = 2000)
        )
        
        sem_output_summary_fxn(fit)
        summary(fit)
        
        #prune(midtroph, leave the rest)

        
        
        
#AK Predators=============
        names(df_master_scaled)
        sem_latent <- "
          my_latent =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                          Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                            Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                            Peff_Harbour_s_WS_adultyr +
                            Peff_DFA_halibutBiomassAge8plus_predAK_adultyr
          SAR  ~  CPUE + Growth + my_latent
        
        "
        
        sem_latent <- "
          my_latent =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                          Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                            Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                            Peff_Harbour_s_WS_adultyr +
                            Peff_DFA_halibutBiomassAge8plus_predAK_adultyr
          SAR  ~  CPUE +  my_latent 
        
        "
        
        fit <- sem(
          sem_latent,
          data         = df_master_scaled,
          missing      = "ML",
          bounds       = "pos.var",
          optim.method = "BFGS",
          control      = list(iter.max = 2000)
        )
        
        sem_output_summary_fxn(fit)
        summary(fit)
        plot(fit)

        
        #pruned growth, now everything signif but some var are neg
#start 2nd cluster akpred-------
        sem_latent <- "
          my_latent =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                          Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                            Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                            Peff_Harbour_s_WS_adultyr +
                            Peff_DFA_halibutBiomassAge8plus_predAK_adultyr +
                            Peff_DFA_AllSeaLionsEMB_adultyr
          SAR  ~  CPUE +  my_latent 
        
        "
        
        fit <- sem(
          sem_latent,
          data         = df_master_scaled,
          missing      = "ML",
          bounds       = "pos.var",
          optim.method = "BFGS",
          control      = list(iter.max = 2000)
        )
        
        sem_output_summary_fxn(fit)
        summary(fit)
        
        #pruned fur seals, but left everything else in the model
        
        
        
#final confirmatory model-------------
        names(df_master_scaled)
        
        sem_latent <- "
            Latent_sardine =~  
                             X01_DFA_sumPreyOfPrey_planktonJun_smoltyr + 
                             X02_DFA_SeaNettle_smoltyr +
                             X04_marketsquid_GAM_smoltyr + 
                             X05_DFA_abundSardine_smoltyr
            Growth ~   Latent_sardine
            
            CPUE ~ Hake_DFA
            

  
            Latent_Peff =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                          Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                            Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                            Peff_Harbour_s_WS_adultyr +
                            Peff_DFA_halibutBiomassAge8plus_predAK_adultyr +
                            Peff_DFA_AllSeaLionsEMB_adultyr
          
          SAR  ~  Growth + CPUE +  Latent_Peff 

"
        
        #adding pink and akprey to the rest of the model doesn't work
#          SAR  ~  Growth +  Latent_AKprey+ Latent_Peff 
            # Latent_AKprey =~ X12_biomassDecaShelfSum_smoltyr + 
            #                 X13_EGOA_DFA_seabirds_smoltyr +
            #                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr +
            #                 X12_biomassMysidShelfSum_smoltyr
        
        fit <- sem(
          sem_latent,
          data         = df_master_scaled,
          missing      = "ML",
          bounds       = "pos.var",
          optim.method = "BFGS",
          control      = list(iter.max = 2000)
        )
        
        sem_output_summary_fxn(fit)
        summary(fit)
        
        
        
        