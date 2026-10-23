# ==============================================================================
# Script: 10d_run_ak_adult_temp_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Execute AK & Adult-Stage Temperature/HCI Gating Sweep on SAR:
#          CPUE ~ Hake_smoltyr
#          SAR  ~ CPUE + Peff
#          Sweeps lower bounds a in [0.0, 0.5].
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION & DATA LOADING
# ------------------------------------------------------------------------------
years_included <- "19yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

# Load raw unscaled data
raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)

# Clean empty/unnamed columns
raw_df <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Invert HCI for NCC predators so High = Warm / Compressed
if ("X01_habCompInd_smoltyr" %in% names(raw_df)) {
  raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr
}
if ("X01_habCompInd_adultyr" %in% names(raw_df)) {
  raw_df$X01_habCompInd_adultyr <- -1 * raw_df$X01_habCompInd_adultyr
}

# Filter complete cases
if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PREDATORS & MATCH SPATIAL ENVIRONMENTAL COLUMNS
# ------------------------------------------------------------------------------
all_cols <- names(filtered_df)

# Blacklist
blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

# Extract X10 (adultyr ONLY), X15, and X11_ssl_seak_pup_pred
candidate_preds <- all_cols[grepl("^(X10_|X15_)", all_cols) | all_cols == "X11_ssl_seak_pup_pred"]
candidate_preds <- candidate_preds[!grepl("^X10_.*smoltyr$", candidate_preds)]
candidate_preds <- candidate_preds[!sapply(candidate_preds, function(x) any(sapply(blacklist, grepl, x)))]

# Build matching table for predators and their corresponding environmental gate
pred_env_lookup <- tibble(pred_col = candidate_preds) %>%
  mutate(
    stage_suffix = ifelse(grepl("adultyr$", pred_col), "adultyr", "smoltyr"),
    region = case_when(
      grepl("^X15_", pred_col) | pred_col == "X11_ssl_seak_pup_pred" ~ "AK",
      grepl("^X10_", pred_col) ~ "NCC_adult",
      TRUE ~ "AK"
    ),
    spatial_col = case_when(
      region == "AK" & stage_suffix == "smoltyr" ~ "X21_sst_egoa_junjulaug_smoltyr",
      region == "AK" & stage_suffix == "adultyr" ~ "X21_sst_egoa_junjulaug_adultyr",
      region == "NCC_adult" ~ paste0("X01_habCompInd_", stage_suffix)
    )
  ) %>%
  filter(spatial_col %in% all_cols)

cat(sprintf("Identified %d AK/Adult stage predators for SAR sweep.\n", nrow(pred_env_lookup)))

# ------------------------------------------------------------------------------
# STEP 2: MIN-MAX SCALING FUNCTION
# ------------------------------------------------------------------------------
scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Standard 0-1 scaling for base model variables
minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# ------------------------------------------------------------------------------
# STEP 3: AK / ADULT TEMPERATURE SWEEP FUNCTION
# ------------------------------------------------------------------------------
test_ak_adult_temp_gating <- function(pred_col, spatial_col, df_raw, df_base) {
  
  P_norm  <- df_base[[pred_col]]
  raw_env <- df_raw[[spatial_col]]
  
  # --- Raw Baseline Model (CPUE ~ Hake; SAR ~ CPUE + Raw Predator) ---
  base_model_str <- paste0(
    target_cpue, " ~ ", hake_col, "\n",
    target_sar,  " ~ ", target_cpue, " + ", pred_col
  )
  
  base_fit <- tryCatch(
    sem(base_model_str, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), 
    error = function(e) NULL
  )
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- tryCatch(parameterEstimates(base_fit, standardized = TRUE), error = function(e) NULL)
  if (is.null(base_pe) || !"rhs" %in% names(base_pe)) return(NULL)
  
  base_pe_sub <- base_pe %>% filter(op == "~" & lhs == target_sar & rhs == pred_col)
  if (nrow(base_pe_sub) == 0) return(NULL)
  
  base_est <- as.numeric(base_pe_sub$est[1])
  base_std <- as.numeric(base_pe_sub$std.all[1])
  base_p   <- as.numeric(base_pe_sub$pvalue[1])
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_sar]])
  
  # Grid range for dynamic scaling bound (a)
  a_bounds <- seq(0.0, 0.5, by = 0.05)
  
  best_abs     <- NULL
  best_abs_aic <- Inf
  
  best_neg     <- NULL
  best_neg_aic <- Inf
  
  # --- Grid Search over Scaling Bounds (a) ---
  for (a in a_bounds) {
    # Scale environmental gate to [a, 1.0]
    env_gate <- scale_custom_bounds(raw_env, a = a, b = 1.0)
    
    # Direct Multiplier: High SST/HCI (1.0) = Full Overlap; Cold (a) = Attenuated Overlap
    P_eff <- P_norm * env_gate
    
    if (is.na(sd(P_eff)) || sd(P_eff) < 1e-6) next
    
    df_temp <- df_base
    df_temp$P_ak_eff <- P_eff
    
    gated_model_str <- paste0(
      target_cpue, " ~ ", hake_col, "\n",
      target_sar,  " ~ ", target_cpue, " + P_ak_eff"
    )
    
    fit <- tryCatch(
      sem(gated_model_str, data = df_temp, std.lv = TRUE, missing = "ML", warn = FALSE), 
      error = function(e) NULL
    )
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      pe <- tryCatch(parameterEstimates(fit, standardized = TRUE), error = function(e) NULL)
      
      if (!is.null(pe) && "rhs" %in% names(pe)) {
        pe_sub <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_ak_eff")
        
        if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
          curr_aic <- as.numeric(fitMeasures(fit, "aic"))
          curr_est <- as.numeric(pe_sub$est[1])
          curr_std <- as.numeric(pe_sub$std.all[1])
          curr_p   <- as.numeric(pe_sub$pvalue[1])
          curr_r2  <- as.numeric(inspect(fit, "r2")[[target_sar]])
          
          res_obj <- list(
            minmax_a  = a,
            aic       = curr_aic,
            delta_aic = curr_aic - base_aic,
            beta_peff = curr_est,
            beta_std  = curr_std,
            pvalue    = curr_p,
            r2        = curr_r2
          )
          
          # 1. Track Absolute Best AIC Model
          if (curr_aic < best_abs_aic) {
            best_abs_aic <- curr_aic
            best_abs     <- res_obj
          }
          
          # 2. Track Best AIC Model with Negative Beta
          if (curr_est < 0 && curr_aic < best_neg_aic) {
            best_neg_aic <- curr_aic
            best_neg     <- res_obj
          }
        }
      }
    }
  }
  
  # --- Safe Unpacking of Best ABS Model ---
  if (!is.null(best_abs)) {
    abs_a     <- best_abs$minmax_a
    abs_beta  <- best_abs$beta_peff
    abs_std   <- best_abs$beta_std
    abs_p     <- best_abs$pvalue
    abs_aic   <- best_abs$aic
    abs_daic  <- best_abs$delta_aic
    abs_r2    <- best_abs$r2
  } else {
    abs_a <- abs_beta <- abs_std <- abs_p <- abs_aic <- abs_daic <- abs_r2 <- NA_real_
  }
  
  # --- Safe Unpacking of Best NEG Model ---
  if (!is.null(best_neg)) {
    neg_a     <- best_neg$minmax_a
    neg_beta  <- best_neg$beta_peff
    neg_p     <- best_neg$pvalue
    neg_aic   <- best_neg$aic
    neg_daic  <- best_neg$delta_aic
    neg_r2    <- best_neg$r2
  } else {
    neg_a <- neg_beta <- neg_p <- neg_aic <- neg_daic <- neg_r2 <- NA_real_
  }
  
  # Diagnostic Evaluators
  abs_aic_improved <- if (!is.na(abs_daic)) (abs_daic < 0) else FALSE
  abs_is_negative  <- if (!is.na(abs_beta)) (abs_beta < 0) else FALSE
  abs_sign_flipped <- if (!is.na(abs_beta)) (base_est > 0 && abs_beta < 0) else FALSE
  abs_more_neg     <- if (!is.na(abs_beta)) (abs_beta < base_est) else FALSE
  
  return(tibble(
    predator            = pred_col,
    spatial_col         = spatial_col,
    base_beta           = base_est,
    base_beta_std       = base_std,
    base_pvalue         = base_p,
    base_aic            = base_aic,
    base_r2             = base_r2,
    
    # Absolute Best AIC Model Results
    abs_minmax_a        = abs_a,
    abs_beta_peff       = abs_beta,
    abs_beta_std        = abs_std,
    abs_pvalue          = abs_p,
    abs_final_aic       = abs_aic,
    abs_delta_aic       = abs_daic,
    abs_final_r2        = abs_r2,
    
    # 4 Diagnostic Evaluators
    Q1_aic_improved     = abs_aic_improved,
    Q2_beta_is_negative = abs_is_negative,
    Q3_sign_flipped_pos_to_neg = abs_sign_flipped,
    Q4_impact_increased = abs_more_neg,
    
    # Best Negative Beta Model Results
    neg_minmax_a        = neg_a,
    neg_beta_peff       = neg_beta,
    neg_pvalue          = neg_p,
    neg_final_aic       = neg_aic,
    neg_delta_aic       = neg_daic,
    neg_final_r2        = neg_r2
  ))
}

# ------------------------------------------------------------------------------
# STEP 4: EXECUTE SWEEP WITH MONITORING
# ------------------------------------------------------------------------------
cat(sprintf("Executing AK/Adult Temperature Gating Sweep on %s dataset...\n", years_included))

sweep_results <- list()

for (i in 1:nrow(pred_env_lookup)) {
  pred_col    <- pred_env_lookup$pred_col[i]
  spatial_col <- pred_env_lookup$spatial_col[i]
  
  cat(sprintf("\n[%d/%d] Testing AK/Adult predator: %s | Spatial: %s | Time: %s\n", 
              i, nrow(pred_env_lookup), pred_col, spatial_col, format(Sys.time(), "%H:%M:%S")))
  
  if (pred_col %in% names(minmax_base) && spatial_col %in% names(filtered_df)) {
    res_df <- test_ak_adult_temp_gating(pred_col = pred_col, spatial_col = spatial_col, df_raw = filtered_df, df_base = minmax_base)
    
    if (!is.null(res_df)) {
      sweep_results[[length(sweep_results) + 1]] <- res_df
      
      print(
        res_df %>% select(predator, base_beta, abs_minmax_a, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg)
      )
    }
  }
}

ak_sweep_table <- bind_rows(sweep_results) %>% arrange(abs_final_aic)

# Print Final Summary Table
cat("\n==============================================================================\n")
cat(sprintf("   AK & ADULT TEMPERATURE GATING SWEEP RESULTS (SAR ~ Peff) - %s DATASET    \n", years_included))
cat("==============================================================================\n\n")

print(
  ak_sweep_table %>% 
    select(predator,base_aic,abs_final_aic,base_r2,abs_final_r2, base_beta, abs_minmax_a, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg),
  n = 50
)

# Export Results
out_csv <- file.path(data_out_dir, paste0("Goal3_AK_Adult_Temperature_Gating_Sweep_", years_included, ".csv"))
write_csv(ak_sweep_table, out_csv)

cat(sprintf("\nSweep results successfully saved to: %s\n", out_csv))



#results-------------
# AK & ADULT TEMPERATURE GATING SWEEP RESULTS (SAR ~ Peff) - 24yr DATASET    
# predator                                          base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_minmax_a abs_beta_peff abs_delta_aic Q1_aic_improved Q3_sign_flipped_pos_to_neg
# <chr>                                                <dbl>         <dbl>   <dbl>        <dbl>     <dbl>        <dbl>         <dbl>         <dbl> <lgl>           <lgl>                     
#   1 X15_DFA_sleeperSharks_smoltyr                      -1.05          -5.97    0.409        0.516   -0.253          0          -0.657        -4.91   TRUE            FALSE                     
# 2 X10_Harbour_s_2yrLead_WS_adultyr                   -3.87          -4.64    0.443        0.464   -0.546          0.45       -0.642        -0.766  TRUE            FALSE                     
# 3 X15_DFA_sleeperSharks_adultyr                      -2.38          -4.49    0.453        0.497   -0.285          0.1        -0.507        -2.12   TRUE            FALSE                     
# 4 X15_ArrowtoothFlounderBiomass_predAK_smoltyr        0.299         -4.07    0.369        0.454   -0.202          0          -0.467        -4.37   TRUE            FALSE                     
# 5 X10_Harbour_s_2yrLead_WS_smoltyr                    1.66          -2.79    0.345        0.464   -0.226          0          -0.970        -4.45   TRUE            FALSE                     
# 6 X10_Northern_f_s_2yrLead_WS_smoltyr                 2.49          -1.94    0.301        0.423    0.0256         0          -0.488        -4.43   TRUE            TRUE                      
# 7 X15_ArrowtoothFlounderBiomass_predAK_adultyr        0.992         -1.74    0.346        0.413   -0.187          0          -0.417        -2.73   TRUE            FALSE                     
# 8 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr   2.33          -0.290   0.308        0.362   -0.0802         0          -0.585        -2.62   TRUE            FALSE                     
# 9 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr   1.96           0.274   0.321        0.366   -0.138          0          -0.426        -1.69   TRUE            FALSE                     
# 10 X10_Californian_s_l_2yrLead_WS_smoltyr             -5.23           0.419   0.515        0.372    0.485          0.5         0.339         5.65   FALSE           FALSE                     
# 11 X15_salmonSharkGoA_predAK_adultyr                   1.60           0.974   0.328        0.344   -0.222          0          -0.420        -0.623  TRUE            FALSE                     
# 12 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr          0.883          1.18    0.361        0.339    0.238          0          -0.304         0.296  FALSE           TRUE                      
# 13 X10_AllSeaLionsEMB_2025_raw_sc_smoltyr              1.78           1.35    0.322        0.336   -0.201          0          -0.252        -0.428  TRUE            FALSE                     
# 14 X10_AllSeaLionsEMB_smoltyr                          1.78           1.35    0.322        0.336   -0.201          0          -0.252        -0.428  TRUE            FALSE                     
# 15 X15_sablefishBiomass_predAK_adultyr                 1.39           1.43    0.333        0.333    0.228          0.5         0.262         0.0352 FALSE           FALSE                     
# 16 X10_Northern_f_s_2yrLead_WS_adultyr                 0.211          1.55    0.367        0.328    0.279          0.5         0.242         1.34   FALSE           FALSE                     
# 17 X15_PacificCodBiomass_predAK_smoltyr                0.0557         1.58    0.377        0.341    0.270          0.5         0.234         1.52   FALSE           FALSE                     
# 18 X15_spinyDogfishGoA_predAK_smoltyr                  2.50           1.69    0.301        0.322    0.0226         0          -0.221        -0.803  TRUE            TRUE                      
# 19 X15_PacificCodBiomass_predAK_adultyr                2.46           1.99    0.302        0.314    0.0389         0          -0.148        -0.474  TRUE            TRUE                      
# 20 X15_sablefishRecruitment_predAK_smoltyr             2.36           2.18    0.306        0.310   -0.0674         0          -0.132        -0.181  TRUE            FALSE                     
# 21 X10_Californian_s_l_2yrLead_WS_adultyr              0.507          2.20    0.363        0.311    0.260          0.5         0.107         1.69   FALSE           FALSE                     
# 22 X15_sablefishBiomass_predAK_smoltyr                 1.99           2.20    0.316        0.310    0.131          0.5         0.125         0.213  FALSE           FALSE                     
# 23 X15_salmonSharkGoA_predAK_smoltyr                   2.46           2.21    0.298        0.315    0.0565         0          -0.197        -0.241  TRUE            TRUE                      
# 24 X11_ssl_seak_pup_pred                               1.00           2.21    0.354        0.315    0.183          0.5         0.105         1.21   FALSE           FALSE                     
# 25 X15_sablefishRecruitment_predAK_adultyr             2.32           2.26    0.305        0.308    0.0774         0           0.136        -0.0530 TRUE            FALSE                     
# 26 X10_DFA_ssl.est.wholerange_2yrLead_adultyr          1.94           2.34    0.319        0.307    0.132          0          -0.0753        0.396  FALSE           TRUE                      
# 27 X10_AllSeaLionsEMB_adultyr                          2.49           2.41    0.301        0.303    0.0301         0          -0.0729       -0.0822 TRUE            TRUE                      
# 28 X10_AllSeaLionsEMB_2025_raw_sc_adultyr              2.49           2.41    0.301        0.303    0.0301         0          -0.0729       -0.0822 TRUE            TRUE                      
# > 