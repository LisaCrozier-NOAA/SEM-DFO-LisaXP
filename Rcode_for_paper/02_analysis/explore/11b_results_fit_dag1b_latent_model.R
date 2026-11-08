sem_output_summary_fxn(fit_dag1b_wgoa)

==============================================
  MODEL SUMMARY REPORT               
==============================================
  
  --- Model Fit Metrics ---
  pvalue     cfi     aic   rmsea 
0.000   0.573 242.640   0.431 

--- R-Squared (Endogenous Variables) ---
  No target R2 variables found in this model.

--- Regression Path Ranking (By P-Value) ---
  lhs op         rhs    est    se      z pvalue
SAR  ~ Latent_WGOA  0.801 0.398  2.013  0.044
SAR  ~ P_ocean_eff -0.374 0.151 -2.479  0.013
SAR  ~        CPUE  0.764 0.177  4.310  0.000

==============================================
  > sem_output_summary_fxn(fit_dag1b_egoa)

==============================================
  MODEL SUMMARY REPORT               
==============================================
  
  --- Model Fit Metrics ---
  pvalue     cfi     aic   rmsea 
0.066   0.859 188.944   0.224 

--- R-Squared (Endogenous Variables) ---
  No target R2 variables found in this model.

--- Regression Path Ranking (By P-Value) ---
  lhs op         rhs    est    se      z pvalue
SAR  ~        CPUE  0.299 0.154  1.941  0.052
SAR  ~ Latent_EGOA  0.635 0.288  2.203  0.028
SAR  ~ P_ocean_eff -0.532 0.138 -3.857  0.000

==============================================
  > sem_output_summary_fxn(fit_dag1b_pink)

==============================================
  MODEL SUMMARY REPORT               
==============================================
  
  --- Model Fit Metrics ---
  pvalue    cfi    aic  rmsea 
NA  1.000 54.457  0.000 

--- R-Squared (Endogenous Variables) ---
  No target R2 variables found in this model.

--- Regression Path Ranking (By P-Value) ---
  lhs op         rhs    est    se      z pvalue
SAR  ~        Pink -0.288 0.128 -2.246  0.025
SAR  ~ P_ocean_eff -0.388 0.137 -2.834  0.005
SAR  ~        CPUE  0.506 0.137  3.681  0.000

==============================================