#Here is the complete script to sweep across all combinations of Sea Lion, Herring, and Capelin candidate series.
#It constructs the standardized I_ssl index using fixed integer weights ($\text{Hazard} = -1 \times [\text{SSL} - 0.3(\text{SSL} \cdot \text{SST}) + 2.0(\text{SSL} \cdot \text{Herring}) + 1.0(\text{SSL} \cdot \text{Capelin})]$) for every candidate trio, 
#fits the SEM, and compares each model's AIC directly against the PDO benchmark ($\text{AIC} = 105.74$).


save.image(file="outputs_4/myshark_myssl.rdata")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

# =============================================================================
# 1. IDENTIFY ALL CANDIDATE VARIABLES
# =============================================================================

# A. Sea Lion Candidates from ssl.dat
ssl_cands <- c(
  "ssl_model_eric",
  "ssl_west_pup_pred",
  "ssl_seak_pup_pred",
  "ssl_west_np_pred",
  "ssl_east_pup_pred",
  "ssl_east_np_pred",
  "ssl_cent_pup_pred",
  "ssl_cent_np_pred"
)
ssl_cands <- ssl_cands[ssl_cands %in% names(ssl.dat)]

# B. Herring Candidates from guild.dfasAK
#herr_cands <- names(guild.dfas1)[grepl("herr", names(guild.dfas1), ignore.case = TRUE)]
herr_cands <- names(ssl.dat)[grepl("herr", names(ssl.dat), ignore.case = TRUE)][1]
#herr_cands <- names(ak_yr)[grepl("herr", names(ak_yr), ignore.case = TRUE)]

# C. Capelin Candidates from guild.dfasAK
#cap_cands  <- names(guild.dfas1)[grepl("capelin|wgoa_cap", names(guild.dfas1), ignore.case = TRUE)]
cap_cands  <- names(ak_yr)[grepl("capelin|wgoa_cap", names(ak_yr), ignore.case = TRUE)][2]

# D. Shared Environmental & Structural Predictors
sst_col <- "sst_wgoa_coastwatch_junjulaug"

cat(sprintf("Found %d SSL, %d Herring, and %d Capelin candidates.\n", 
            length(ssl_cands), length(herr_cands), length(cap_cands)))
cat(sprintf("Total permutations to test: %d models.\n", 
            length(ssl_cands) * length(herr_cands) * length(cap_cands)))

# =============================================================================
# 2. MERGE CANDIDATE DATASETS & CLEAN NAMES
# =============================================================================

# Merge ssl.dat columns with main guild.dfas1 by year
# (Assuming year is the primary key)
combined_df <- guild.dfas1 %>%
  left_join(
    ssl.dat %>% select(year, all_of(c(ssl_cands, sst_col,herr_cands))), 
    by = "year"
  )%>%
left_join(
  ak_yr %>% select(year, all_of(cap_cands)), 
  by = "year"
)

# Define PDO Baseline AIC Benchmark
pdo_benchmark_aic <- 105.735

# =============================================================================
# 3. HELPER FUNCTION TO FIT FIXED-WEIGHT SEM FOR A SINGLE TRIO
# =============================================================================
fit_fixed_issl_combination <- function(ssl_v, herr_v, cap_v, data_in) {
  
  # Select and clean variables for 1998-2021
  d_sub <- data_in %>%
    filter(year >= 1998, year <= 2021) %>%
    select(
      year,
      x16_sar,
      x07_dfa_cpue_int_spr_jun_hw,
      x09_dfa_hake_age5plus,
      x15_shark_enso_roll2,
      ssl_raw  = all_of(ssl_v),
      sst_raw  = all_of(sst_col),
      herr_raw = all_of(herr_v),
      cap_raw  = all_of(cap_v)
    ) %>%
    filter(if_all(everything(), ~ !is.na(.x) & is.finite(.x)))
  
  n_obs <- nrow(d_sub)
  if (n_obs < 24) return(NULL) # Strict year-tracking: Only complete N = 24 series
  
  # Calculate Fixed Integer Weight Hazard Index
  # Hazard = -1 * [ 1.0*SSL - 0.3*(SSL*SST) + 2.0*(SSL*Herr) + 1.0*(SSL*Cap) ]
  d_mod <- d_sub %>%
    mutate(
      ssl_z  = as.vector(scale(ssl_raw)),
      sst_z  = as.vector(scale(sst_raw)),
      herr_z = as.vector(scale(herr_raw)),
      cap_z  = as.vector(scale(cap_raw)),
      
      # Uncalibrated Fixed-Weight Hazard Sum
      I_SSL_raw = -1 * ( 1.0 * ssl_z - 0.3 * (ssl_z * sst_z) + 2.0 * (ssl_z * herr_z) + 1.0 * (ssl_z * cap_z) ),
      
      # Standardize final index for SEM stability
      I_SSL_simple = as.vector(scale(I_SSL_raw))
    )
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_shark_enso_roll2 + I_SSL_simple
  '
  
  fit <- tryCatch({
    sem(sem_text, data = d_mod, missing = "listwise", std.lv = TRUE, warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "cfi", "rmsea", "srmr"))
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  get_param <- function(lhs_v, rhs_v, field) {
    val <- pe %>% filter(lhs == lhs_v, op == "~", rhs == rhs_v) %>% pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  r2_sar <- inspect(fit, "rsquare")["x16_sar"]
  
  tibble(
    ssl_variable     = ssl_v,
    herring_variable = herr_v,
    capelin_variable = cap_v,
    n_years          = n_obs,
    
    aic              = fm[["aic"]],
    delta_aic_vs_pdo = fm[["aic"]] - pdo_benchmark_aic,
    beats_pdo        = fm[["aic"]] <= pdo_benchmark_aic,
    
    cfi              = fm[["cfi"]],
    rmsea            = fm[["rmsea"]],
    srmr             = fm[["srmr"]],
    r2_sar           = r2_sar,
    
    # Path Estimates
    b_issl_sar       = get_param("x16_sar", "I_SSL_simple", "std.all"),
    p_issl_sar       = get_param("x16_sar", "I_SSL_simple", "pvalue"),
    
    b_shark_sar      = get_param("x16_sar", "x15_shark_enso_roll2", "std.all"),
    p_shark_sar      = get_param("x16_sar", "x15_shark_enso_roll2", "pvalue")
  )
}

# =============================================================================
# 4. RUN ALL PERMUTATIONS VIA EXPAND.GRID
# =============================================================================
combo_grid <- expand.grid(
  ssl_v  = ssl_cands,
  herr_v = herr_cands,
  cap_v  = cap_cands,
  stringsAsFactors = FALSE
)

message("Running fixed-weight SEM sweep across ", nrow(combo_grid), " permutations...")

sweep_results2 <- purrr::pmap_dfr(
  list(combo_grid$ssl_v, combo_grid$herr_v, combo_grid$cap_v),
  function(s, h, c) fit_fixed_issl_combination(s, h, c, combined_df)
) %>%
  arrange(aic)

# =============================================================================
# 5. DISPLAY RESULTS & PDO BENCHMARK SUMMARY
# =============================================================================

total_models <- nrow(sweep_results)
beating_pdo  <- sum(sweep_results$beats_pdo, na.rm = TRUE)
pct_beating  <- (beating_pdo / total_models) * 100

cat("\n=================================================================\n")
cat(sprintf("SWEEP RESULTS SUMMARY (N = %d Valid Permutations)\n", total_models))
cat("=================================================================\n")
cat(sprintf("PDO Model Benchmark AIC:  %.3f\n", pdo_benchmark_aic))
cat(sprintf("Models Beating PDO Bar:   %d out of %d (%.1f%%)\n", beating_pdo, total_models, pct_beating))
cat(sprintf("Best Overall Model AIC:   %.3f (Delta AIC vs PDO: %.3f)\n", 
            min(sweep_results$aic), min(sweep_results$delta_aic_vs_pdo)))
cat("=================================================================\n\n")

# Top 15 Models
cat("=== TOP 15 FIXED-WEIGHT SEM COMBINATIONS ===\n")
print(
  sweep_results %>%
    select(ssl_variable, herring_variable, capelin_variable, aic, delta_aic_vs_pdo, r2_sar, b_issl_sar, p_issl_sar) %>%
    head(15)
)

# Export Full Results to CSV
write_csv(sweep_results, file.path(out_dir, "ssl_index_fixed_weight_full_permutation_sweep.csv"))


#I got lucky:
#   ssl_variable      herring_variable      capelin_variable     aic delta_aic_vs_pdo r2_sar b_issl_sar  p_issl_sar
# 1 ssl_seak_pup_pred x13_stka_herr_matbiom x13_wgoa_cap_pcod   103.            -2.72  0.759     -0.531 0.000000342
# 2 ssl_seak_pup_pred x13_stka_herr_matbiom x13_mid_il_capelin  104.            -1.52  0.733     -0.536 0.000000821
# 3 ssl_east_pup_pred x13_stka_herr_matbiom x13_wgoa_cap_pcod   107.             1.43  0.716     -0.482 0.0000143  
# 4 ssl_east_pup_pred x13_stka_herr_matbiom x13_mid_il_capelin  108.             2.41  0.693     -0.477 0.0000337  
# 5 ssl_west_np_pred  x13_stka_herr_matbiom x13_wgoa_cap_pcod   110.             3.92  0.677     -0.449 0.000122   
# 6 ssl_seak_pup_pred x13_egoa_herring      x13_wgoa_cap_pcod   110.             4.13  0.668     -0.450 0.000160   
# 7 ssl_west_np_pred  x13_stka_herr_matbiom x13_mid_il_capelin  110.             4.41  0.662     -0.447 0.000185   
# 8 ssl_west_pup_pred x13_stka_herr_matbiom x13_wgoa_cap_pcod   110.             4.47  0.666     -0.444 0.000194   
# 9 ssl_west_pup_pred x13_stka_herr_matbiom x13_mid_il_capelin  111.             5.16  0.648     -0.438 0.000346   
# 10 ssl_seak_pup_pred x13_egoa_herring      x13_mid_il_capelin  112.             5.84  0.626     -0.434 0.000582   
# 11 ssl_east_pup_pred x13_egoa_herring      x13_wgoa_cap_pcod   113.             6.92  0.629     -0.402 0.00137    
# 12 ssl_east_pup_pred x13_egoa_herring      x13_mid_il_capelin  114.             7.86  0.603     -0.387 0.00288    
# 13 ssl_east_np_pred  x13_stka_herr_matbiom x13_wgoa_cap_pcod   114.             8.01  0.621     -0.377 0.00494    
# 14 ssl_east_np_pred  x13_stka_herr_matbiom x13_mid_il_capelin  114.             8.17  0.605     -0.384 0.00448    
# 15 ssl_west_np_pred  x13_egoa_herring      x13_wgoa_cap_pcod   114.             8.17  0.605     -0.377 0.00352    
