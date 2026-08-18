#Simplify SSL index


# MODEL COMPARISON SUMMARY:
# 1. Empirically Tuned I_SSL Model: AIC = 103.46, SAR R2 = 0.739
# 2. Simplified Fixed-Weight Model:  AIC = 102.61, SAR R2 = 0.748
# 3. Baseline PDO Model:            AIC = 105.74, SAR R2 = 0.706

#when I used herr_avg, it didn't quite beat pdo. but still damn good
#2. Simplified Fixed-Weight Model:  AIC = 105.97, SAR R2 = 0.706

#see next script, but here are those results too:
#I could have used x13_wgoa_cap_pcod instead of capelin_avg, but not both of the DFAs
#   ssl_variable      herring_variable      capelin_variable     aic delta_aic_vs_pdo r2_sar b_issl_sar  p_issl_sar
# 1 ssl_seak_pup_pred x13_stka_herr_matbiom x13_wgoa_cap_pcod   103.            -2.72  0.759     -0.531 0.000000342
# 2 ssl_seak_pup_pred x13_stka_herr_matbiom x13_mid_il_capelin  104.            -1.52  0.733     -0.536 0.000000821
# 3 ssl_east_pup_pred x13_stka_herr_matbiom x13_wgoa_cap_pcod   107.             1.43  0.716     -0.482 0.0000143  
# 4 ssl_east_pup_pred x13_stka_herr_matbiom x13_mid_il_capelin  108.             2.41  0.693     -0.477 0.0000337  
# 5 ssl_west_np_pred  x13_stka_herr_matbiom x13_wgoa_cap_pcod   110.             3.92  0.677     -0.449 0.000122   
# 6 ssl_seak_pup_pred x13_egoa_herring      x13_wgoa_cap_pcod   110.             4.13  0.668     -0.450 0.000160   


suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# -----------------------------------------------------------------------------
# 1. Build Simplified Fixed-Weight Index (I_SSL_simple)
# -----------------------------------------------------------------------------
df_top_ssl <- df_top_ssl %>%
  mutate(
    # Standardize underlying raw variables
    ssl    = as.vector(scale(ssl_seak_pup_pred)),
    sst    = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
  #  f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    f_herr = as.vector(scale(herr_avg)),
    f_cap  = as.vector(scale(capelin_avg)),
    
    # Calculate Simplified Hazard Index (No empirical regression tuning!)
    # Negated so higher = higher predation hazard
    I_SSL_simple = -1 * ( 1.0 * ssl - 0.3 * (ssl * sst) + 2.0 * (ssl * f_herr) + 1.0 * (ssl * f_cap) ),
    
    # Standardize the resulting composite index
    I_SSL_simple = as.vector(scale(I_SSL_simple))
  )

# -----------------------------------------------------------------------------
# 2. Fit SEM with Fixed-Weight Index
# -----------------------------------------------------------------------------
fixed_model_text <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2 + I_SSL_simple
'

fit_fixed <- sem(fixed_model_text, data = df_top_ssl, missing = "listwise")

# -----------------------------------------------------------------------------
# 3. Model Report & Comparison
# -----------------------------------------------------------------------------
cat("\n=== FIXED-WEIGHT I_SSL SEM REPORT ===\n")
summary(fit_fixed, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# Compare AIC against empirical I_SSL model (103.464) and PDO model (105.735)
aic_fixed <- fitMeasures(fit_fixed, "aic")
r2_sar_fixed <- inspect(fit_fixed, "rsquare")["x16_sar"]

cat("\n------------------------------------------------\n")
cat("MODEL COMPARISON SUMMARY:\n")
cat(sprintf("1. Empirically Tuned I_SSL Model: AIC = 103.46, SAR R2 = 0.739\n"))
cat(sprintf("2. Simplified Fixed-Weight Model:  AIC = %.2f, SAR R2 = %.3f\n", aic_fixed, r2_sar_fixed))
cat(sprintf("3. Baseline PDO Model:            AIC = 105.74, SAR R2 = 0.706\n"))
cat("------------------------------------------------\n")
