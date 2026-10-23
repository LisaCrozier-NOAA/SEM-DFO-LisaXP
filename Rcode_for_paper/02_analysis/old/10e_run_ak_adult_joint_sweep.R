# ==============================================================================
# Script: 10e_run_ak_adult_joint_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Execute Combined SST/HCI + Alternative Prey Attenuation Sweep on SAR:
#          CPUE ~ Hake_smoltyr
#          SAR  ~ CPUE + Peff
#          Fixed lower bound a = 0.2. Sweeps 1 and 2 alternative prey pools.
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
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

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

guild_pred_lookup <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PREDATORS & MATCH SPATIAL & PREY COLUMNS
# ------------------------------------------------------------------------------
all_cols <- names(filtered_df)

# Blacklist
blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

# Extract candidate predators
candidate_preds <- all_cols[grepl("^(X10_|X15_)", all_cols) | all_cols == "X11_ssl_seak_pup_pred"]
candidate_preds <- candidate_preds[!grepl("^X10_.*smoltyr$", candidate_preds)]
candidate_preds <- candidate_preds[!sapply(candidate_preds, function(x) any(sapply(blacklist, grepl, x)))]

# Standard min-max scaling function
scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# ------------------------------------------------------------------------------
# STEP 2: COMBINED SST + ALT PREY SWEEP FUNCTION (FIXED a = 0.2)
# ------------------------------------------------------------------------------
test_ak_adult_joint_gating <- function(pred_col, df_raw, df_base, max_combo_size = 2) {
  
  P_norm <- df_base[[pred_col]]
  
  row_data <- guild_pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  stage_suffix <- ifelse(grepl("adultyr$", pred_col), "adultyr", "smoltyr")
  region <- ifelse(grepl("^X15_", pred_col) | pred_col == "X11_ssl_seak_pup_pred", "AK", "NCC_adult")
  
  spatial_col <- case_when(
    region == "AK" & stage_suffix == "smoltyr" ~ "X21_sst_egoa_junjulaug_smoltyr",
    region == "AK" & stage_suffix == "adultyr" ~ "X21_sst_egoa_junjulaug_adultyr",
    region == "NCC_adult" ~ paste0("X01_habCompInd_", stage_suffix)
  )
  if (!spatial_col %in% names(df_raw)) return(NULL)
  
  raw_env <- df_raw[[spatial_col]]
  
  # Environmental Gate with FIXED a = 0.2
  env_gate <- scale_custom_bounds(raw_env, a = 0.2, b = 1.0)
  P_temp   <- P_norm * env_gate
  
  # Extract candidate prey (excluding spatial col and predator itself!)
  all_cands <- unlist(strsplit(row_data$all_candidates_str, ", "))
  prey_candidates <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
  prey_candidates <- intersect(prey_candidates, names(df_base))
  
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
  
  best_abs     <- NULL
  best_abs_aic <- Inf
  
  best_neg     <- NULL
  best_neg_aic <- Inf
  
  w_grid <- seq(0.1, 1.0, by = 0.1)
  
  # --- 1. Test Temperature-Only Variant (Fixed a = 0.2, No Alt Prey) ---
  df_base$P_eff_temp <- P_temp
  temp_model_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_temp")
  
  temp_fit <- tryCatch(sem(temp_model_str, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  if (!is.null(temp_fit) && lavInspect(temp_fit, "converged")) {
    pe_t <- tryCatch(parameterEstimates(temp_fit, standardized = TRUE), error = function(e) NULL)
    if (!is.null(pe_t) && "rhs" %in% names(pe_t)) {
      pe_t_sub <- pe_t %>% filter(op == "~" & lhs == target_sar & rhs == "P_eff_temp")
      if (nrow(pe_t_sub) > 0 && !is.na(pe_t_sub$pvalue[1])) {
        t_aic <- as.numeric(fitMeasures(temp_fit, "aic"))
        t_est <- as.numeric(pe_t_sub$est[1])
        t_std <- as.numeric(pe_t_sub$std.all[1])
        t_p   <- as.numeric(pe_t_sub$pvalue[1])
        t_r2  <- as.numeric(inspect(temp_fit, "r2")[[target_sar]])
        
        res_t <- list(
          variant_type = "SST Only",
          prey_combo   = "None",
          weight_w     = NA_real_,
          aic          = t_aic,
          delta_aic    = t_aic - base_aic,
          beta_peff    = t_est,
          beta_std     = t_std,
          pvalue       = t_p,
          r2           = t_r2
        )
        if (t_aic < best_abs_aic) { best_abs_aic <- t_aic; best_abs <- res_t }
        if (t_est < 0 && t_aic < best_neg_aic) { best_neg_aic <- t_aic; best_neg <- res_t }
      }
    }
  }
  
  # --- 2. Test Combined SST + Alternative Prey Variants ---
  if (length(prey_candidates) > 0) {
    max_k <- min(length(prey_candidates), max_combo_size)
    
    for (k in 1:max_k) {
      combos <- combn(prey_candidates, k, simplify = FALSE)
      
      for (cb in combos) {
        # Build additive forage pool
        if (length(cb) == 1) {
          forage_pool <- df_base[[cb[1]]]
        } else {
          raw_sum <- rowSums(df_base[, cb, drop = FALSE], na.rm = TRUE)
          forage_pool <- scale_custom_bounds(raw_sum, 0.0, 1.0)
        }
        
        for (w in w_grid) {
          # Combined Gating: P_eff = P_temp * (1 - w * forage_pool)
          P_joint <- P_temp * (1 - w * forage_pool)
          
          if (is.na(sd(P_joint)) || sd(P_joint) < 1e-6) next
          
          df_base$P_joint_eff <- P_joint
          joint_model_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_joint_eff")
          
          fit <- tryCatch(sem(joint_model_str, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
          
          if (!is.null(fit) && lavInspect(fit, "converged")) {
            pe <- tryCatch(parameterEstimates(fit, standardized = TRUE), error = function(e) NULL)
            
            if (!is.null(pe) && "rhs" %in% names(pe)) {
              pe_sub <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_joint_eff")
              
              if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
                curr_aic <- as.numeric(fitMeasures(fit, "aic"))
                curr_est <- as.numeric(pe_sub$est[1])
                curr_std <- as.numeric(pe_sub$std.all[1])
                curr_p   <- as.numeric(pe_sub$pvalue[1])
                curr_r2  <- as.numeric(inspect(fit, "r2")[[target_sar]])
                
                res_obj <- list(
                  variant_type = paste0("SST + ", k, " AltPrey"),
                  prey_combo   = paste(cb, collapse = " + "),
                  weight_w     = w,
                  aic          = curr_aic,
                  delta_aic    = curr_aic - base_aic,
                  beta_peff    = curr_est,
                  beta_std     = curr_std,
                  pvalue       = curr_p,
                  r2           = curr_r2
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
      }
    }
  }
  
  # --- Safe Unpacking of Best ABS Model ---
  if (!is.null(best_abs)) {
    abs_variant <- best_abs$variant_type
    abs_prey    <- best_abs$prey_combo
    abs_w       <- best_abs$weight_w
    abs_beta    <- best_abs$beta_peff
    abs_std     <- best_abs$beta_std
    abs_p       <- best_abs$pvalue
    abs_aic     <- best_abs$aic
    abs_daic    <- best_abs$delta_aic
    abs_r2      <- best_abs$r2
  } else {
    abs_variant <- abs_prey <- "None"
    abs_w <- abs_beta <- abs_std <- abs_p <- abs_aic <- abs_daic <- abs_r2 <- NA_real_
  }
  
  # --- Safe Unpacking of Best NEG Model ---
  if (!is.null(best_neg)) {
    neg_variant <- best_neg$variant_type
    neg_prey    <- best_neg$prey_combo
    neg_w       <- best_neg$weight_w
    neg_beta    <- best_neg$beta_peff
    neg_p       <- best_neg$pvalue
    neg_aic     <- best_neg$aic
    neg_daic    <- best_neg$delta_aic
    neg_r2      <- best_neg$r2
  } else {
    neg_variant <- neg_prey <- "None"
    neg_w <- neg_beta <- neg_p <- neg_aic <- neg_daic <- neg_r2 <- NA_real_
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
    abs_variant_type    = abs_variant,
    abs_prey            = abs_prey,
    abs_weight_w        = abs_w,
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
    neg_variant_type    = neg_variant,
    neg_prey            = neg_prey,
    neg_weight_w        = neg_w,
    neg_beta_peff       = neg_beta,
    neg_pvalue          = neg_p,
    neg_final_aic       = neg_aic,
    neg_delta_aic       = neg_daic,
    neg_final_r2        = neg_r2
  ))
}

# ------------------------------------------------------------------------------
# STEP 3: EXECUTE COMBINED SWEEP WITH MONITORING
# ------------------------------------------------------------------------------
cat(sprintf("Executing AK/Adult Joint SST + AltPrey Sweep on %s dataset...\n", years_included))

sweep_results <- list()

for (i in seq_along(candidate_preds)) {
  pred_col <- candidate_preds[i]
  cat(sprintf("\n[%d/%d] Testing joint sweep for: %s | Time: %s\n", 
              i, length(candidate_preds), pred_col, format(Sys.time(), "%H:%M:%S")))
  
  if (pred_col %in% names(minmax_base)) {
    res_df <- test_ak_adult_joint_gating(pred_col = pred_col, df_raw = filtered_df, df_base = minmax_base, max_combo_size = 2)
    
    if (!is.null(res_df)) {
      sweep_results[[length(sweep_results) + 1]] <- res_df
      
      print(
        res_df %>% select(predator, base_beta, abs_variant_type, abs_prey, abs_beta_peff, abs_delta_aic, Q1_aic_improved)
      )
    }
  }
}

ak_joint_sweep_table <- bind_rows(sweep_results) %>% arrange(abs_final_aic)

# Print Final Summary Table
cat("\n==============================================================================\n")
cat(sprintf("   AK & ADULT JOINT SST + ALTPREY SWEEP RESULTS (SAR ~ Peff) - %s DATASET    \n", years_included))
cat("==============================================================================\n\n")

print(
  ak_joint_sweep_table %>% 
    select(predator,base_aic,abs_final_aic,base_r2,abs_final_r2, base_beta, abs_variant_type, abs_prey, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved),
  n = 50
)

# Export Results
out_csv <- file.path(data_out_dir, paste0("Goal3_AK_Adult_Joint_SST_AltPrey_Sweep_", years_included, ".csv"))
write_csv(ak_joint_sweep_table, out_csv)

cat(sprintf("\nJoint sweep results successfully saved to: %s\n", out_csv))


#RESULTS 24YR-----------
# predator                                          base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_variant_type abs_prey             abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved
# 1 X10_Harbour_s_2yrLead_WS_adultyr                   -3.87         -9.39     0.443        0.570   -0.546  SST + 2 AltPrey  X12_DFA_biomassEuph…          0.9        -2.32         -5.52  TRUE           
# 2 X15_ArrowtoothFlounderBiomass_predAK_adultyr        0.992        -6.35     0.346        0.555   -0.187  SST + 1 AltPrey  X13_capelin_WGoA_ad…          1          -0.561        -7.34  TRUE           
# 3 X15_ArrowtoothFlounderBiomass_predAK_smoltyr        0.299        -5.93     0.369        0.503   -0.202  SST + 1 AltPrey  X13_capelin_WGoA_sm…          1          -0.496        -6.23  TRUE           
# 4 X15_sablefishBiomass_predAK_smoltyr                 1.99         -5.63     0.316        0.504    0.131  SST + 2 AltPrey  X13_pollock_age1plu…          1          -2.35         -7.62  TRUE           
# 5 X15_DFA_sleeperSharks_adultyr                      -2.38         -5.29     0.453        0.529   -0.285  SST + 2 AltPrey  X13_pollock_age1plu…          1          -0.541        -2.92  TRUE           
# 6 X15_DFA_sleeperSharks_smoltyr                      -1.05         -4.82     0.409        0.495   -0.253  SST + 1 AltPrey  X13_capelin_WGoA_sm…          1          -0.586        -3.77  TRUE           
# 7 X11_ssl_seak_pup_pred                               1.00         -4.08     0.354        0.458    0.183  SST + 2 AltPrey  X13_sitkaHerring_EG…          1          -0.544        -5.08  TRUE           
# 8 X15_salmonSharkGoA_predAK_smoltyr                   2.46         -3.90     0.298        0.456    0.0565 SST + 1 AltPrey  X13_capelin_WGoA_sm…          1          -1.63         -6.36  TRUE           
# 9 X10_Californian_s_l_2yrLead_WS_adultyr              0.507        -2.36     0.363        0.438    0.260  SST + 2 AltPrey  X04_marketsquid_GAM…          1           1.15         -2.87  TRUE           
# 10 X15_spinyDogfishGoA_predAK_smoltyr                  2.50         -1.59     0.301        0.384    0.0226 SST + 2 AltPrey  X13_pollock_age1plu…          1          -1.26         -4.08  TRUE           
# 11 X15_PacificCodBiomass_predAK_smoltyr                0.0557       -1.19     0.377        0.372    0.270  SST + 2 AltPrey  X13_sitkaHerring_EG…          1          -0.706        -1.25  TRUE           
# 12 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr   1.96         -0.911    0.321        0.416   -0.138  SST + 1 AltPrey  X13_capelin_WGoA_ad…          1          -0.493        -2.87  TRUE           
# 13 X10_Northern_f_s_2yrLead_WS_adultyr                 0.211        -0.865    0.367        0.393    0.279  SST + 1 AltPrey  X09_DFA_HakeAge5Plu…          1           1.17         -1.08  TRUE           
# 14 X10_DFA_ssl.est.wholerange_2yrLead_adultyr          1.94         -0.681    0.319        0.390    0.132  SST + 2 AltPrey  X04_marketsquid_GAM…          1           1.04         -2.62  TRUE           
# 15 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr   2.33         -0.591    0.308        0.382   -0.0802 SST + 1 AltPrey  X13_capelin_WGoA_sm…          1          -0.531        -2.92  TRUE           
# 16 X15_sablefishBiomass_predAK_adultyr                 1.39         -0.513    0.333        0.404    0.228  SST + 2 AltPrey  X13_pollock_age1plu…          1          -1.52         -1.90  TRUE           
# 17 X15_salmonSharkGoA_predAK_adultyr                   1.60          0.0163   0.328        0.399   -0.222  SST + 2 AltPrey  X13_pollock_age1plu…          0.9        -1.36         -1.58  TRUE           
# 18 X15_PacificCodBiomass_predAK_adultyr                2.46          0.213    0.302        0.394    0.0389 SST + 2 AltPrey  X13_pollock_age1plu…          1          -0.432        -2.25  TRUE           
# 19 X10_AllSeaLionsEMB_adultyr                          2.49          0.261    0.301        0.362    0.0301 SST + 2 AltPrey  X04_marketsquid_GAM…          1           1.75         -2.23  TRUE           
# 20 X10_AllSeaLionsEMB_2025_raw_sc_adultyr              2.49          0.261    0.301        0.362    0.0301 SST + 2 AltPrey  X04_marketsquid_GAM…          1           1.75         -2.23  TRUE           
# 21 X15_sablefishRecruitment_predAK_smoltyr             2.36          1.12     0.306        0.339   -0.0674 SST + 2 AltPrey  X13_sitkaHerring_EG…          1          -0.278        -1.24  TRUE           
# 22 X15_sablefishRecruitment_predAK_adultyr             2.32          2.00     0.305        0.324    0.0774 SST + 2 AltPrey  X13_pollock_age1plu…          1          -0.646        -0.319 TRUE           
# > 