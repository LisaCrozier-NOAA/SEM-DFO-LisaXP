#RESULTS
guild.dfasAK <- read.csv(file.path("data_Lisa/guild.dfasAK.csv")) %>% clean_names()

#Define index----------
df_top_ssl <- df_top_ssl %>%
  mutate(
    # Standardize underlying raw variables
    ssl    = as.vector(scale(ssl_seak_pup_pred)),
    sst    = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    f_cap  = as.vector(scale(x13_mid_il_capelin)),
    
    # Calculate Simplified Hazard Index (No empirical regression tuning!)
    # Negated so higher = higher predation hazard
    I_SSL_simple = -1 * ( 1.0 * ssl - 0.3 * (ssl * sst) + 2.0 * (ssl * f_herr) + 1.0 * (ssl * f_cap) ),
    
    # Standardize the resulting composite index
    I_SSL_simple = as.vector(scale(I_SSL_simple)),
    issl = as.vector(scale(I_SSL_simple))
  )

guild.dfasAK_ishark_issl<-full_join(guild.dfasAK,
                                    df_top_ssl %>% select(year,ishark,issl),
                                    by="year")
write.csv(guild.dfasAK_ishark_issl,   file.path(out_dir, "guild.dfasAK_ishark_issl.csv"), row.names = FALSE)

#==================
out_dir    <- "copilot/outputs_ssl"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#setup the data ----

guild.dfasAK <- read.csv(file.path("data_Lisa/guild.dfasAK.csv")) %>% clean_names()

ssl.dat <- read.csv("data_Lisa/ssl.dat.csv") %>% clean_names() %>%
  select(year,ssl_seak_pup_pred,sst_wgoa_coastwatch_junjulaug) %>%
  left_join(guild.dfasAK %>%
              select(year, x13_stka_herr_matbiom,x13_egoa_herring,x13_mid_il_capelin),
              join_by("year"))

#Shark function------
ishark_fxn <- function(N_shark, 
                       T_mt,
                       Tref_mt=mean(T_mt, na.rm = TRUE), 
                       T_ot, 
                       Q10 = 2, 
                       overlap_slope = 1, 
                       transform_shark = "log1p", 
                       scale_final = TRUE) {
  
  # 1. Input checks
  if (length(N_shark) != length(T_mt) || length(N_shark) != length(T_ot)) {
    stop("Input vectors (N_shark, T_mt, T_ot) must be of equal length.")
  }
  
  # 2. Transform Abundance (Keep non-negative)
  if (transform_shark == "log1p") {
    minv <- suppressWarnings(min(N_shark, na.rm = TRUE))
    N_pos <- log1p(N_shark - minv)
  } else if (transform_shark == "raw") {
    N_pos <- N_shark
  } else {
    stop("transform_shark must be either 'log1p' or 'raw'.")
  }
  
  # 3. Calculate Metabolic Multiplier (Mt) using raw T_mt differential
  M_t     <- Q10^((T_mt - Tref_mt) / 10)
  
  # 4. Calculate Spatial Overlap Multiplier (Ot) using Z-scaled T_ot for Logistic curve
  # Safe Z-scaling function
  scale_vector <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  
  T_ot_z <- scale_vector(T_ot)
  O_t    <- plogis(overlap_slope * T_ot_z) * 2
  
  # 5. Compute Integrated Predation Index in Positive Space
  I_Shark_raw <- N_pos * M_t * O_t
  
  # 6. Standardize or Return Raw Index
  if (scale_final) {
    I_Shark_out <- scale_vector(I_Shark_raw)
  } else {
    I_Shark_out <- I_Shark_raw
  }
  
  return(I_Shark_out)
}


shark_raw<-read.csv(file.path("data_Lisa/shark_wide.csv")) %>% 
  clean_names() %>%
  select(year,goa_pacific_sleeper_shark) %>%
  full_join(guild.dfasAK %>%
            select(year,x21_sst_wgoa_coastwatch_junjulaug),
            by="year") %>%
  mutate(ishark=ishark_fxn(
    N_shark=goa_pacific_sleeper_shark,
    T_mt=x21_sst_wgoa_coastwatch_junjulaug,
    T_ot=x21_sst_wgoa_coastwatch_junjulaug))
head(shark_raw)

fit_dat<-guild.dfasAK %>% full_join(shark_raw %>% select(year,ishark),by="year") %>%
  select(year,x07_dfa_cpue_int_spr_jun_hw,x09_dfa_hake_age5plus,x16_sar,ishark)


# Standardize inputs and construct product terms
df_top_ssl <- ssl.dat %>%
  mutate(
    # Core predictors
    ssl = as.vector(scale(ssl_seak_pup_pred)),
    sst = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_cap = as.vector(scale(x13_mid_il_capelin)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    
    # Linear interaction products
    ssl_x_sst  = ssl * sst,
    ssl_x_cap  = ssl * f_cap,
    ssl_x_herr = ssl * f_herr
  )%>%
  full_join(fit_dat,by="year")

#fit SEM w/ linear interactions to define I-ssl-------

#  SEM model syntax incorporating SSL interactions
single_model_text <- '
  # Structural paths
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + 
            ishark + 
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

# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# x07_dfa_cpue_int_spr_jun_hw ~                                                      
#   x09_df_hk_g5pl                -0.372    0.098   -3.802    0.000   -0.372   -0.613
# x16_sar ~                                                                          
#   x07_df_cp_n___                 0.373    0.120    3.111    0.002    0.373    0.302
# ishark                        -0.420    0.104   -4.026    0.000   -0.420   -0.469
# ssl                            0.377    0.107    3.508    0.000    0.377    0.421
# ssl_x_sst                      0.113    0.101    1.118    0.264    0.113    0.109
# ssl_x_herr                     0.596    0.125    4.761    0.000    0.596    0.505
# ssl_x_cap                      0.195    0.138    1.417    0.156    0.195    0.168


#Define index----------
df_top_ssl <- df_top_ssl %>%
  mutate(
    # Standardize underlying raw variables
    ssl    = as.vector(scale(ssl_seak_pup_pred)),
    sst    = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    f_cap  = as.vector(scale(x13_mid_il_capelin)),
    
    # Calculate Simplified Hazard Index (No empirical regression tuning!)
    # Negated so higher = higher predation hazard
    I_SSL_simple = -1 * ( 1.0 * ssl - 0.3 * (ssl * sst) + 2.0 * (ssl * f_herr) + 1.0 * (ssl * f_cap) ),
    
    # Standardize the resulting composite index
    I_SSL_simple = as.vector(scale(I_SSL_simple))
  )

#double check fit------

single_model_text <- '
  # Structural paths
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + 
            ishark + 
            I_SSL_simple
'

# Fit model with complete cases
fit_single2 <- sem(
  single_model_text,
  data = df_top_ssl,
  std.lv = TRUE,
  missing = "listwise",
  warn = FALSE
)



# Output summary
summary(fit_single2, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

#Results
# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# x07_dfa_cpue_int_spr_jun_hw ~                                                      
#   x09_df_hk_g5pl                -0.372    0.098   -3.802    0.000   -0.372   -0.613
# x16_sar ~                                                                          
#   x07_df_cp_n___                 0.361    0.135    2.673    0.008    0.361    0.290
# ishark                        -0.395    0.108   -3.651    0.000   -0.395   -0.438
# I_SSL_simple                  -0.376    0.099   -3.811    0.000   -0.376   -0.417
