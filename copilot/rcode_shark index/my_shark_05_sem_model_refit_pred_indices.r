#Run optimizer on both indices:
# beta_sst   beta_herr    beta_cap    alpha_ot 
# 5.0000000 -41.1277691 -20.0000000   0.3723296 
#1 ssl: 5 sst: -40 herr: -20 cap
#0.25: 1: -8: -4



library(lavaan)
library(tidyverse)
library(janitor)


#refit both ssl and shark indices
    source("functions/shark_index_fxn.r")
    source("functions/ssl_index_fxn.r")


#load data---------
salmon_dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
  clean_names()  %>% 
  select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
names(salmon_dat)


shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
  clean_names() %>%
  select(year,goa_pacific_sleeper_shark,goa_salmon_shark)  
head(shark_dat)


sst<-read.csv(file.path("data_Lisa/ak_yr.csv"),row.names=NULL) %>%
  select(year,contains("sst"),contains("176to226m"),contains("pdo"),contains("enso"),contains("195to205m"),contains("246to255")) %>%
  filter(year>1997);names(sst)

#too many NAs in the other deep water temps, swln_temp_spr_176to226m had no NAs 

summary(sst)

sst0<-sst %>%
  select(year,contains("coastwatch"),contains("swln_temp_spr_176to226m"),contains("pdo"),contains("enso"))
summary(sst0)

ssl_dat<-read.csv(file.path("copilot/outputs_2/ssl.dat.csv"))%>% 
  clean_names() %>%
  select(!contains("_avg"))
names(ssl_dat)

goa_prey<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
  clean_names() %>%
  select(1:12)
names(goa_prey)

names(shark_dat)
names(sst0)

pred_index_dat <- left_join(sst0 %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf),
                            shark_dat %>% select(year,goa_pacific_sleeper_shark),by="year") %>%
  left_join(ssl_dat %>% select(year,ssl_seak_pup_pred),by="year") %>%
  left_join(goa_prey %>% select(year,egoa_bio_stka_herr_matbiom,wgoa_bio_mid_il_capelin),by="year") %>%
  left_join(salmon_dat, by="year")
names(pred_index_dat)

#SEM---------
fit_model <- sem('
            x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
            x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + i_ssl + i_shark
          ', data = pred_index_dat, std.lv = TRUE, missing = "ML", warn = FALSE)

#optimization code-----
library(lavaan)
library(tidyverse)

# 1. Helper Function: Safe Scale with Variance Threshold
safe_scale <- function(x) {
  # Check if vector is flat / near zero variance
  v <- var(x, na.rm = TRUE)
  if (is.na(v) || v < 1e-6) {
    return(NULL) # Explicitly signal invalid variance
  }
  return(as.vector(scale(x)))
}

# 2. Index Calculation Functions
calc_i_ssl <- function(p_ssl, ssl, sst, herr, cap) {
  # p_ssl = c(beta_sst, beta_herr, beta_cap)
  ssl * (1.0 + p_ssl[1]*sst + p_ssl[2]*herr + p_ssl[3]*cap)
}

calc_i_shark <- function(p_shark, shark_s, t_raw, t_ref, t_ot, Q10 = 2) {
  # p_shark = c(alpha_ot)
  Mt <- Q10^((t_raw - t_ref) / 10)
  Ot_scaled <- 1.0 + p_shark[1] * t_ot
  shark_s * Mt * Ot_scaled
}

# 3. Objective Function with Crash-Protection
sem_objective_fxn <- function(params, df_data) {
  p_ssl   <- params[1:3]
  p_shark <- params[4]
  
  # Standardize raw components
  ssl_s   <- safe_scale(df_data$ssl_seak_pup_pred)
  sst_s   <- safe_scale(df_data$sst_wgoa_coastwatch_junjulaug)
  herr_s  <- safe_scale(df_data$egoa_bio_stka_herr_matbiom)
  cap_s   <- safe_scale(df_data$wgoa_bio_mid_il_capelin)
  shark_s <- safe_scale(df_data$goa_pacific_sleeper_shark)
  
  t_raw   <- df_data$sst_wgoa_coastwatch_junjulaug
  t_ref   <- mean(t_raw, na.rm = TRUE)
  t_ot    <- df_data$pdo_djf
  
  # Calculate raw indices
  i_ssl_raw   <- calc_i_ssl(p_ssl, ssl_s, sst_s, herr_s, cap_s)
  i_shark_raw <- calc_i_shark(p_shark, shark_s, t_raw, t_ref, t_ot)
  
  # Standardize compound indices
  i_ssl_scaled   <- safe_scale(i_ssl_raw)
  i_shark_scaled <- safe_scale(i_shark_raw)
  
  # Heavy penalty if either index collapsed to zero variance
  if (is.null(i_ssl_scaled) || is.null(i_shark_scaled)) {
    return(1e8)
  }
  
  df_temp <- df_data %>%
    mutate(
      i_ssl   = i_ssl_scaled,
      i_shark = i_shark_scaled
    )
  
  # Fit SEM model
  fit <- try(
    sem('
      x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
      x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + i_ssl + i_shark
    ', data = df_temp, std.lv = TRUE, missing = "ML", warn = FALSE),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error") || !fit@Fit@converged) {
    return(1e8)
  }
  
  return(-2 * logLik(fit)) 
}

# 4. Run Optimization with L-BFGS-B Bounded Search
init_params <- c(
  beta_sst  =  0.3, 
  beta_herr = -2.0, 
  beta_cap  = -1.0, 
  alpha_ot  =  0.5
)

opt_results <- optim(
  par     = init_params,
  fn      = sem_objective_fxn,
  df_data = pred_index_dat,
  method  = "L-BFGS-B",
  lower   = c(-10, -50, -20, -10), # Keep search space reasonable
  upper   = c( 5,  10,  10,  5),
  control = list(maxit = 500)
)

# 5. Extract Optimal Parameters
p_opt <- opt_results$par
names(p_opt) <- c("beta_sst", "beta_herr", "beta_cap", "alpha_ot")

# 6. Reconstruct Fitted Dataset
ssl_s   <- safe_scale(pred_index_dat$ssl_seak_pup_pred)
sst_s   <- safe_scale(pred_index_dat$sst_wgoa_coastwatch_junjulaug)
herr_s  <- safe_scale(pred_index_dat$egoa_bio_stka_herr_matbiom)
cap_s   <- safe_scale(pred_index_dat$wgoa_bio_mid_il_capelin)
shark_s <- safe_scale(pred_index_dat$goa_pacific_sleeper_shark)

t_raw   <- pred_index_dat$sst_wgoa_coastwatch_junjulaug
t_ref   <- mean(t_raw, na.rm = TRUE)
t_ot    <- pred_index_dat$pdo_djf

pred_index_dat_fitted <- pred_index_dat %>%
  mutate(
    i_ssl   = safe_scale(calc_i_ssl(p_opt[1:3], ssl_s, sst_s, herr_s, cap_s)),
    i_shark = safe_scale(calc_i_shark(p_opt[4], shark_s, t_raw, t_ref, t_ot))
  )

# 7. Fit Final SEM Model & Print Summary
final_sem_fit <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + i_ssl + i_shark
', data = pred_index_dat_fitted, std.lv = TRUE, missing = "ML", warn = FALSE)

cat("--- Optimized Parameters ---\n")
print(p_opt)

cat("\n--- SEM Model Summary ---\n")
summary(final_sem_fit, fit.measures = TRUE, standardized = TRUE)

# beta_sst   beta_herr    beta_cap    alpha_ot 
# 5.0000000 -41.1277691 -20.0000000   0.3723296 
#1 ssl: 5 sst: -40 herr: -20 cap
#0.25: 1: -8: -4
