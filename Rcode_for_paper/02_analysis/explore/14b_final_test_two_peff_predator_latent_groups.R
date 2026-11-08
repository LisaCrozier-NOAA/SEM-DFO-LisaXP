

#14b_final_test_two_peff_predator_latent_groups.R
#sharks or pink alone worked best:

# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue    cfi    aic  rmsea 
# NA  1.000 56.005  0.000 
# 
# --- R-Squared (Endogenous Variables) ---
#   No target R2 variables found in this model.
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op                                 rhs    est    se      z pvalue
# SAR  ~ DFA_sleeperSharkBSAI_predAK_smoltyr -0.414 0.139 -2.974  0.003
# SAR  ~                                CPUE  0.549 0.139  3.948  0.000
# 
# ==============================================

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
#   lhs op  rhs    est    se      z pvalue
# SAR  ~ Pink -0.315 0.148 -2.128  0.033
# SAR  ~ CPUE  0.647 0.148  4.379  0.000
# 
# ==============================================
# 

# ------------------------------------------------------------------------------
# 1b. SPECIFY DAG 1B WITH ONE BIPOLAR LATENT VARIABLES -- mammals------
# ------------------------------------------------------------------------------
# Model 1: Group 1 Alone (Mammal Latent)
mammal_latent <- "
  Latent_Mammal =~ DFA_Steller_s_l_WS_adultyr + 
                   harbour_seal_DFO_WCVI_smoltyr + 
                   Killer_whales_NR_BC_adultyr + 
                   ssl_seak_pup_pred_smoltyr

  SAR ~ CPUE + Pink + Latent_Mammal
"

fit_mammal <- sem(
  mammal_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_mammal)
sem_output_summary_fxn(fit_mammal)

#all latent var p=0, but Latent_Mammal  p=0.45 -- NS!!!
#prune 1 at a time?
  Latent_Mammal =~ DFA_Steller_s_l_WS_adultyr + 
                   harbour_seal_DFO_WCVI_smoltyr + 
                   Killer_whales_NR_BC_adultyr + 
                   ssl_seak_pup_pred_smoltyr
mammal_latent <- "
  SAR ~ CPUE + Pink 
"

fit_mammal <- sem(
  mammal_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_mammal)
sem_output_summary_fxn(fit_mammal)

# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue    cfi    aic  rmsea 
# NA   1.00  59.78   0.00 
# 
# --- R-Squared (Endogenous Variables) ---
#   No target R2 variables found in this model.
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op                       rhs    est    se      z pvalue
# SAR  ~ ssl_seak_pup_pred_smoltyr -0.184 0.143 -1.289  0.197
# SAR  ~                      Pink -0.300 0.143 -2.089  0.037
# SAR  ~                      CPUE  0.639 0.143  4.469  0.000

==============================================
# prune DFA_Steller_s_l_WS_adultyr, harbour_seal_DFO_WCVI_smoltyr, Killer_whales_NR_BC_adultyr
#ssl_seak_pup_pred_smoltyr was NS, but closer than any of the others

#Fish group========
#prune pink: + Pink
fish_latent <- "
  Latent_Demersal =~ a*DFA_sleeperSharkBSAI_predAK_smoltyr + 
                     a*sablefishAK_biomass3plus_smoltyr

  SAR ~ CPUE  + Latent_Demersal
"
fit_fish <- sem(
  fish_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_fish)
sem_output_summary_fxn(fit_fish)

#the fish latent was better than pink!
#try w/o latent, just observed

fish_latent <- "
  SAR ~ CPUE  + DFA_sleeperSharkBSAI_predAK_smoltyr 
"
fit_fish_obs <- sem(
  fish_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_fish_obs)
sem_output_summary_fxn(fit_fish_obs)
#prune sablefish-------

# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue    cfi    aic  rmsea 
# NA  1.000 56.005  0.000 
# 
# --- R-Squared (Endogenous Variables) ---
#   No target R2 variables found in this model.
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op                                 rhs    est    se      z pvalue
# SAR  ~ DFA_sleeperSharkBSAI_predAK_smoltyr -0.414 0.139 -2.974  0.003
# SAR  ~                                CPUE  0.549 0.139  3.948  0.000
# 
# ==============================================