#take home -- go back to ssl and try to get herring to work


# Extract residuals from the top model fit
# 1 1 Predator + 1 Prey                    x15_shark_enso_roll2                       x21_pdo_djf      24 105.7353         -16.53250  -0.5002842 1.533218e-05  -0.5755494 2.666867e-06               TRUE
source("functions/sem_output_summary_fxn.r")

#get I_ssl, created in my_ssl_03_create_I_ssl.r
df_top_ssl<-read.csv(file.path(out_dir,"ssl_index_dat.csv"))
#recall:
# df_top_ssl <- ssl.dat %>%
#   mutate(
#     # Core predictors
#     ssl = as.vector(scale(ssl_seak_pup_pred)),
#     sst = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
#     f_cap = as.vector(scale(capelin_avg)),
#     f_herr = as.vector(scale(x13_stka_herr_matbiom)),
#     
#     # Linear interaction products
#     ssl_x_sst  = ssl * sst,
#     ssl_x_cap  = ssl * f_cap,
#     ssl_x_herr = ssl * f_herr
#   )

# x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + 
#   x15_shark_enso_roll2 + 
#   ssl + 
#   ssl_x_sst + 
#   ssl_x_herr + 
#   ssl_x_cap

# Regressions:
#                               Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# ssl                            0.352    0.151    2.324    0.020    0.352    0.346
# ssl_x_sst                     -0.103    0.142   -0.725    0.469   -0.103   -0.087
# ssl_x_herr                     0.736    0.168    4.392    0.000    0.736    0.549
# ssl_x_cap                      0.351    0.138    2.548    0.011    0.351    0.305


ishark_issl_model_text <- '
            x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
            x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2 + I_SSL
                      '


ishark_issl_fit <- sem(ishark_issl_model_text, data = df_top_ssl, missing = "listwise")
sem_output_summary_fxn(ishark_issl_fit)


# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue     cfi     aic      rmsea    agfi 
#   0.215     0.963   103.464   0.143   0.700 
# 
# --- R-Squared (Endogenous Variables) ---
#   x07_dfa_cpue_int_spr_jun_hw                     x16_sar 
# 0.376                       0.739 
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op                         rhs    est    se      z pvalue
# x16_sar  ~        x15_shark_enso_roll2 -0.308 0.107 -2.888  0.004
# x07_dfa_cpue_int_spr_jun_hw  ~       x09_dfa_hake_age5plus -0.613 0.161 -3.802  0.000
# x16_sar  ~ x07_dfa_cpue_int_spr_jun_hw  0.474 0.110  4.330  0.000
# x16_sar  ~                       I_SSL -1.231 0.242 -5.085  0.000
# 
# ==============================================
#   > 



df_top<-guild.dfas1
single_model_text <- '
            x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
            x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2 + x21_pdo_djf
                      '
                               
          
          top_fit <- sem(single_model_text, data = df_top, missing = "listwise")
          sem_output_summary_fxn(top_fit)
 
          
          # ==============================================
          #   MODEL SUMMARY REPORT               
          # ==============================================
          #   
          #   --- Model Fit Metrics ---
          #   pvalue     cfi     aic   rmsea    agfi 
          # 0.441   1.000 105.735   0.000   0.788 
          # 
          # --- R-Squared (Endogenous Variables) ---
          #   x07_dfa_cpue_int_spr_jun_hw                     x16_sar 
          # 0.376                       0.706 
          # 
          # --- Regression Path Ranking (By P-Value) ---
          #   lhs op                         rhs    est    se      z pvalue
          # x16_sar  ~ x07_dfa_cpue_int_spr_jun_hw  0.398 0.118  3.364  0.001
          # x07_dfa_cpue_int_spr_jun_hw  ~       x09_dfa_hake_age5plus -0.613 0.161 -3.802  0.000
          # x16_sar  ~        x15_shark_enso_roll2 -0.500 0.116 -4.324  0.000
          # x16_sar  ~                 x21_pdo_djf -0.576 0.123 -4.695  0.000
          # 
          # ==============================================
          #   >           
          
                   
          # Residuals for SAR_Adult
          # Extract parameter estimates table
          pe <- parameterEstimates(top_fit)
          
          # Extract specific structural coefficients
          b_cpue  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x07_dfa_cpue_int_spr_jun_hw"]
          b_prey  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x21_pdo_djf"] # or whatever prey term is in top_fit
          b_shark <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x15_shark_enso_roll2"]
          
          # Intercept for SAR (if model includes intercepts, otherwise default is 0 for standardized variables)
          b0_sar  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~1"]
          if (length(b0_sar) == 0) b0_sar <- 0
          
          # Compute fitted values
          df_top <- df_top %>%
            mutate(
              sar_hat       = b0_sar + (b_cpue * x07_dfa_cpue_int_spr_jun_hw) + (b_prey * x21_pdo_djf) + (b_shark * x15_shark_enso_roll2),
              sar_residuals = x16_sar - sar_hat
            )
          
          # Check summary: should have mean ~ 0 and non-zero variation
          summary(df_top$sar_residuals)
          
          
          # Plot residual time-series
          x16_ishark_pdo_plot<-
            ggplot(df_top, aes(x = year, y = sar_residuals)) +
            geom_line() + geom_point() +
            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
            labs(title = "Unexplained SAR Residuals from top model: x16 ~ ishark + pdo", y = "Residual (Observed - Predicted)")
          
          
#w/o pdo---------------------------------
          df_sharkonly<-guild.dfas1 %>% select(year,x07_dfa_cpue_int_spr_jun_hw,
                                               x09_dfa_hake_age5plus,
                                               x21_pdo_djf,
                                               x15_shark_enso_roll2,
                                               x16_sar)
          single_model_text <- '
            x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
            x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2 
                      '
          
          
          top_fit <- sem(single_model_text, data = df_sharkonly, missing = "listwise")
          sem_output_summary_fxn(top_fit)
          
          # Residuals for SAR_Adult
          # Extract parameter estimates table
          pe <- parameterEstimates(top_fit)
          
          # Extract specific structural coefficients
          b_cpue  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x07_dfa_cpue_int_spr_jun_hw"]
          b_prey  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x21_pdo_djf"] # or whatever prey term is in top_fit
          b_shark <- pe$est[pe$lhs == "x16_sar" & pe$op == "~" & pe$rhs == "x15_shark_enso_roll2"]
          
          # Intercept for SAR (if model includes intercepts, otherwise default is 0 for standardized variables)
          b0_sar  <- pe$est[pe$lhs == "x16_sar" & pe$op == "~1"]
          if (length(b0_sar) == 0) b0_sar <- 0
          
          # Compute fitted values
          df_sharkonly <- df_sharkonly %>%
            mutate(
              sar_hat       = b0_sar + (b_cpue * x07_dfa_cpue_int_spr_jun_hw) + (b_shark * x15_shark_enso_roll2),
              sar_residuals = x16_sar - sar_hat
            )
          
          # Check summary: should have mean ~ 0 and non-zero variation
          summary(df_sharkonly$sar_residuals)
          
          write.csv(df_sharkonly,file.path(out_dir,"resid_x16_ishark_cpue.csv"))
          
          # Plot residual time-series
          x16_ishark_plot<-
            ggplot(df_sharkonly, aes(x = year, y = sar_residuals)) +
            geom_line() + geom_point() +
            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
            labs(title = "Unexplained SAR Residuals from shark-only model: x16 ~ ishark", y = "Residual (Observed - Predicted)")
          
          print(x16_ishark_plot)
          
          
  #cor plot------------
          # Calculate correlation between true residuals and all dataset variables
          residual_cors <- guild.dfas1 %>%
            filter(year >= 1998, year <= 2021) %>%
            pivot_longer(cols = -year, names_to = "variable", values_to = "value") %>%
            group_by(variable) %>%
            filter(sum(!is.na(value)) >= 12) %>% # Ensure sufficient data
            summarise(
              r_val = cor(value, df_sharkonly$sar_residuals, use = "complete.obs"),
              p_val = cor.test(value, df_sharkonly$sar_residuals)$p.value,
              .groups = "drop"
            ) %>%
            arrange(desc(abs(r_val)))
          
          # View the top 15 candidate variables that explain the leftover SAR variation
          print(head(residual_cors, 15))

        #take home -- go back to ssl and try to get herring to work
          
          # variable                           r_val     p_val
          # 1 x16_sar                            0.751 0.0000231
          # 2 sar_uc                             0.751 0.000332 
          # 3 sar_sr                             0.650 0.000589 
          # 4 x21_pdo_djf                       -0.619 0.00127  
          # 5 x12_dfa_biomass_euph_shelf_sum    -0.567 0.00384  
          # 6 x21_sst_wgoa_coastwatch_junjulaug -0.560 0.00439  
          # 7 x02_krill                          0.541 0.00636  
          # 8 x01_hab_comp_ind                   0.520 0.00915  
          # 9 x04_pacific_pompano               -0.491 0.0148   
          # 10 x12_copepod_biomass_w_go_a        -0.474 0.0192   
          # 11 x09_dfa_hake_age5plus             -0.400 0.0529   
          # 12 x15_salmon_shark_go_a_pred_ak     -0.388 0.0608   
          # 13 x13_stka_herr_matbiom              0.380 0.0668   
          # 14 x03_hake_age1                      0.377 0.0691   
          # 15 x13_egoa_herring                   0.371 0.0740             
                    