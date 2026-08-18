

#results -- plausible sem
df_top_ssl<-read.csv(file.path(out_dir,"ssl_index_dat.csv"))


#setup the data----------
path_ssl_dat <- "copilot/outputs_2/ssl.dat.csv"   # <- update if needed
out_dir      <- "copilot/outputs_4"

ssl.dat <- read_csv(path_ssl_dat, show_col_types = FALSE) %>% clean_names() %>%
  select(-salmon_resid)

df_sharkonly <- read.csv(file.path("copilot/outputs_4/resid_x16_ishark_cpue.csv"))%>%
        rename(sar_ishark_resid=sar_residuals)

guild.dfasAK <- read.csv(file.path("data_Lisa/guild.dfasAK.csv")) %>% clean_names()


ssl.dat<-left_join(ssl.dat, guild.dfasAK %>% 
                     select(year, x13_stka_herr_matbiom,x13_egoa_herring),
                   join_by(year)
          ) %>%
        left_join(df_sharkonly %>% 
                    select(year,x09_dfa_hake_age5plus,sar_ishark_resid,x15_shark_enso_roll2),
                  join_by(year)) 
names(ssl.dat)                            



# Standardize inputs and construct product terms
df_top_ssl <- ssl.dat %>%
  mutate(
    # Core predictors
    ssl = as.vector(scale(ssl_seak_pup_pred)),
    sst = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_cap = as.vector(scale(capelin_avg)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    
    # Linear interaction products
    ssl_x_sst  = ssl * sst,
    ssl_x_cap  = ssl * f_cap,
    ssl_x_herr = ssl * f_herr
  )

#fit SEM w/ linear interactions to define I-ssl-------

# Updated SEM model syntax incorporating SSL interactions
single_model_text <- '
  # Structural paths
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + 
            x15_shark_enso_roll2 + 
            ssl + 
            ssl_x_sst + 
            ssl_x_herr + 
            ssl_x_cap
'

# Fit model with complete cases
fit_single <- sem(
  single_model_text, 
  data = df_top_ssl, 
  std.lv = TRUE, 
  missing = "listwise", 
  warn = FALSE
)

# Output summary
summary(fit_single, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

#coef in I_ssl--------
# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# x07_dfa_cpue_int_spr_jun_hw ~                                                      
#   x09_df_hk_g5pl                -0.613    0.161   -3.802    0.000   -0.613   -0.613
# x16_sar ~                                                                          
#   x07_df_cp_n___                 0.477    0.114    4.172    0.000    0.477    0.469
# x15_shrk_ns_r2                -0.341    0.138   -2.469    0.014   -0.341   -0.335
# ssl                            0.352    0.151    2.324    0.020    0.352    0.346
# ssl_x_sst                     -0.103    0.142   -0.725    0.469   -0.103   -0.087
# ssl_x_herr                     0.736    0.168    4.392    0.000    0.736    0.549
# ssl_x_cap                      0.351    0.138    2.548    0.011    0.351    0.305


#create index_ssl--------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# -----------------------------------------------------------------------------
# 1. Calculate Integrated Sea Lion Index (I_SSL)
# -----------------------------------------------------------------------------
# Pull standardized coefficients from SEM fit
pe <- parameterEstimates(fit_single, standardized = TRUE)

get_std <- function(rhs_var) {
  pe %>% 
    filter(lhs == "x16_sar", op == "~", rhs == rhs_var) %>% 
    pull(std.all)
}

b_ssl      <- get_std("ssl")
b_ssl_sst  <- get_std("ssl_x_sst")
b_ssl_herr <- get_std("ssl_x_herr")
b_ssl_cap  <- get_std("ssl_x_cap")

# Construct I_SSL as the net sea lion response term
df_top_ssl <- df_top_ssl %>%
  mutate(
    # Net effective SSL impact
    # Negate so higher values reflect higher net predation hazard
    I_SSL = -1 * (ssl * (b_ssl + b_ssl_sst * sst + b_ssl_herr * f_herr + b_ssl_cap * f_cap))
  )


write.csv(df_top_ssl,file.path(out_dir,"ssl_index_dat.csv"))
# -----------------------------------------------------------------------------
# 2. Fit Parsimonious SEM with Integrated Indices
# -----------------------------------------------------------------------------
integrated_sem_text <- '
  # Structural paths
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + 
            x15_shark_enso_roll2 + 
            I_SSL
'

fit_integrated <- sem(
  integrated_sem_text, 
  data = df_top_ssl, 
  missing = "listwise", 
  warn = FALSE
)

# Print Summary of Simplified Model
cat("\n=== PARSIMONIOUS SEM SUMMARY (With I_Shark and I_SSL) ===\n")
summary(fit_integrated, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# -----------------------------------------------------------------------------
# 3. Time-Series Plot of I_Shark vs I_SSL vs SAR
# -----------------------------------------------------------------------------
p_indices <- df_top_ssl %>%
  select(year, Observed_SAR = x16_sar, I_Shark = x15_shark_enso_roll2, I_SSL) %>%
  pivot_longer(-year, names_to = "Index", values_to = "Value") %>%
  ggplot(aes(x = year, y = Value, color = Index, linetype = Index)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("Observed_SAR" = "black", "I_Shark" = "firebrick", "I_SSL" = "dodgerblue3")) +
  scale_linetype_manual(values = c("Observed_SAR" = "solid", "I_Shark" = "dashed", "I_SSL" = "dashed")) +
  labs(
    title = "Comparison of Integrated Predator/Buffering Indices over Time",
    x = "Year",
    y = "Standardized Anomaly (z-score)",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(p_indices)



#SEM summary-------------
# lavaan 0.6-20 ended normally after 1 iteration
# 
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of model parameters                         6
# 
# Number of observations                            24
# 
# Model Test User Model:
#   
#   Test statistic                                 4.470
# Degrees of freedom                                 3
# P-value (Chi-square)                           0.215
# 
# Model Test Baseline Model:
#   
#   Test statistic                                47.181
# Degrees of freedom                                 7
# P-value                                        0.000
# 
# User Model versus Baseline Model:
#   
#   Comparative Fit Index (CFI)                    0.963
# Tucker-Lewis Index (TLI)                       0.915
# 
# Loglikelihood and Information Criteria:
#   
#   Loglikelihood user model (H0)                -45.732
# Loglikelihood unrestricted model (H1)        -43.497
# 
# Akaike (AIC)                                 103.464
# Bayesian (BIC)                               110.533
# Sample-size adjusted Bayesian (SABIC)         91.945
# 
# Root Mean Square Error of Approximation:
#   
#   RMSEA                                          0.143
# 90 Percent confidence interval - lower         0.000
# 90 Percent confidence interval - upper         0.398
# P-value H_0: RMSEA <= 0.050                    0.239
# P-value H_0: RMSEA >= 0.080                    0.724
# 
# Standardized Root Mean Square Residual:
#   
#   SRMR                                           0.048
# 
# Parameter Estimates:
#   
#   Standard errors                             Standard
# Information                                 Expected
# Information saturated (h1) model          Structured
# 
# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# x07_dfa_cpue_int_spr_jun_hw ~                                                      
#   x09_df_hk_g5pl                -0.613    0.161   -3.802    0.000   -0.613   -0.613
# x16_sar ~                                                                          
#   x07_df_cp_n___                 0.474    0.110    4.330    0.000    0.474    0.466
# x15_shrk_ns_r2                -0.308    0.107   -2.888    0.004   -0.308   -0.303
# I_SSL                         -1.231    0.242   -5.085    0.000   -1.231   -0.548
# 
# Variances:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# .x07_df_cp_n___    0.598    0.173    3.464    0.001    0.598    0.624
# .x16_sar           0.259    0.075    3.464    0.001    0.259    0.261
# 
# R-Square:
#   Estimate
# x07_df_cp_n___    0.376
# x16_sar           0.739
