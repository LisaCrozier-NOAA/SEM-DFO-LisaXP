# ==============================================================================
# Script: 10b_run_ncc_cpue_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Execute NCC Predator Joint Attenuation Sweep targeting CPUE (CPUE ~ Peff).
#          Tracks BOTH absolute minimum AIC AND lowest AIC with negative beta.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SWITCH & LOAD PREPARED DATA
# ------------------------------------------------------------------------------
years_included <- "19yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

minmax_complete_data <- read.csv(file.path(data_out_dir, paste0("minmax_complete_data_", years_included, ".csv")), stringsAsFactors = FALSE)
guild_pred_lookup    <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)

target_cpue <- "X07_DFA_cpue_IntSprJunHW"

# Filter lookup table strictly for NCC Predators
ncc_lookup <- guild_pred_lookup %>% filter(region == "NCC")

# ------------------------------------------------------------------------------
# STEP 1: DEFINE MINMAX JOINT GATING SWEEP FUNCTION FOR NCC (CPUE TARGET)
# ------------------------------------------------------------------------------
find_best_joint_altprey_temp_minmax_nccpred <- function(pred_col, df_model, max_combo_size = 2) {
  P_norm <- df_model[[pred_col]]
  row_data <- ncc_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  candidates <- unlist(strsplit(row_data$all_candidates_str, ", "))
  candidates <- intersect(candidates[!is.na(candidates) & candidates != ""], names(df_model))
  
  # --- Raw Baseline Model (CPUE ~ Raw Predator) ---
  base_fit <- tryCatch(
    sem(paste0(target_cpue, " ~ ", pred_col), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), 
    error = function(e) NULL
  )
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- parameterEstimates(base_fit) %>% filter(op == "~" & rhs == pred_col)
  base_est <- as.numeric(base_pe$est[1])
  base_p   <- as.numeric(base_pe$pvalue[1])
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_cpue]])
  
  if (length(candidates) == 0) return(NULL)
  
  Anorm_mat <- as.matrix(df_model[, candidates, drop = FALSE])
  w_grid    <- seq(0.1, 1.0, by = 0.1)
  
  # Trackers for Absolute Best AIC vs. Best AIC with Negative Beta
  best_abs <- NULL
  best_abs_aic <- Inf
  
  best_neg <- NULL
  best_neg_aic <- Inf
  
  # --- Grid Search over Gating Combinations ---
  for (k in 1:min(length(candidates), max_combo_size)) {
    combos <- combn(candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      w_combo_grid <- expand.grid(rep(list(w_grid), k))
      
      for (g_idx in 1:nrow(w_combo_grid)) {
        w_vals <- as.numeric(w_combo_grid[g_idx, ])
        
        w1_gate <- 1 - w_vals[1] * Anorm_mat[, cb[1]]
        P_joint <- P_norm * w1_gate
        
        if (k >= 2) {
          w2_gate <- 1 - w_vals[2] * Anorm_mat[, cb[2]]
          P_joint <- P_joint * w2_gate
        }
        
        df_model$P_joint_eff <- P_joint
        
        fit <- tryCatch(
          sem(paste0(target_cpue, " ~ P_joint_eff"), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), 
          error = function(e) NULL
        )
        
        if (!is.null(fit) && lavInspect(fit, "converged")) {
          pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_joint_eff")
          
          if (nrow(pe) > 0 && !is.na(pe$pvalue[1])) {
            curr_aic <- as.numeric(fitMeasures(fit, "aic"))
            curr_est <- as.numeric(pe$est[1])
            curr_p   <- as.numeric(pe$pvalue[1])
            curr_r2  <- as.numeric(inspect(fit, "r2")[[target_cpue]])
            
            res_obj <- list(
              selected_prey = paste(cb, collapse = " + "),
              weights       = paste(round(w_vals, 2), collapse = "; "),
              aic           = curr_aic,
              delta_aic     = curr_aic - base_aic,
              beta_peff     = curr_est,
              pvalue        = curr_p,
              final_r2      = curr_r2
            )
            
            # 1. Update Absolute Best Model
            if (curr_aic < best_abs_aic) {
              best_abs_aic <- curr_aic
              best_abs     <- res_obj
            }
            
            # 2. Update Best Negative Beta Model
            if (curr_est < 0 && curr_aic < best_neg_aic) {
              best_neg_aic <- curr_aic
              best_neg     <- res_obj
            }
          }
        }
      }
    }
  }
  
  # Return combined record
  return(tibble(
    predator         = pred_col,
    base_beta        = base_est,
    base_pvalue      = base_p,
    base_aic         = base_aic,
    base_r2          = base_r2,
    
    # Absolute Lowest AIC Model Results
    abs_prey         = ifelse(!is.null(best_abs), best_abs$selected_prey, NA_character_),
    abs_weights      = ifelse(!is.null(best_abs), best_abs$weights, NA_character_),
    abs_beta_peff    = ifelse(!is.null(best_abs), best_abs$beta_peff, NA_real_),
    abs_pvalue       = ifelse(!is.null(best_abs), best_abs$pvalue, NA_real_),
    abs_final_aic    = ifelse(!is.null(best_abs), best_abs$aic, NA_real_),
    abs_delta_aic    = ifelse(!is.null(best_abs), best_abs$delta_aic, NA_real_),
    abs_final_r2     = ifelse(!is.null(best_abs), best_abs$final_r2, NA_real_),
    
    # Best Negative Beta Model Results
    neg_prey         = ifelse(!is.null(best_neg), best_neg$selected_prey, NA_character_),
    neg_weights      = ifelse(!is.null(best_neg), best_neg$weights, NA_character_),
    neg_beta_peff    = ifelse(!is.null(best_neg), best_neg$beta_peff, NA_real_),
    neg_pvalue       = ifelse(!is.null(best_neg), best_neg$pvalue, NA_real_),
    neg_final_aic    = ifelse(!is.null(best_neg), best_neg$aic, NA_real_),
    neg_delta_aic    = ifelse(!is.null(best_neg), best_neg$delta_aic, NA_real_),
    neg_final_r2     = ifelse(!is.null(best_neg), best_neg$final_r2, NA_real_)
  ))
}

# ------------------------------------------------------------------------------
# STEP 2: EXECUTE ATTENUATION SWEEP ACROSS NCC PREDATORS
# ------------------------------------------------------------------------------
cat(sprintf("Executing NCC Attenuation Sweep (CPUE ~ Peff) on %s dataset...\n", years_included))

ncc_predators       <- ncc_lookup$pred_data_col
joint_sweep_results <- list()

for (i in seq_along(ncc_predators)) {
  print(Sys.time())
  pred_col <- ncc_predators[i]
  cat(sprintf("[%d/%d] Testing NCC predator: %s ...\n", i, length(ncc_predators), pred_col))
  
  if (pred_col %in% names(minmax_complete_data)) {
    res_df <- find_best_joint_altprey_temp_minmax_nccpred(pred_col = pred_col, df_model = minmax_complete_data, max_combo_size = 2)
    if (!is.null(res_df)) {
      joint_sweep_results[[length(joint_sweep_results) + 1]] <- res_df
    }
  }
  print(joint_sweep_results[i])
}

ncc_sweep_table <- bind_rows(joint_sweep_results) %>% arrange(abs_final_aic)

# Print Summary
cat("\n==============================================================================\n")
cat(sprintf("   NCC DUAL SWEEP RESULTS (CPUE ~ Peff) - %s DATASET                         \n", years_included))
cat("==============================================================================\n\n")
print(ncc_sweep_table %>% select(predator, base_beta, abs_prey, abs_beta_peff, abs_delta_aic, neg_prey, neg_beta_peff, neg_delta_aic), n = 50)

# Export Results
out_csv <- file.path(data_out_dir, paste0("Goal3_NCC_CPUE_Sweep_Results_minmax_combo_size_2_", years_included, ".csv"))
write_csv(ncc_sweep_table, out_csv)

cat(sprintf("\nNCC CPUE sweep results successfully written to: %s\n", out_csv))