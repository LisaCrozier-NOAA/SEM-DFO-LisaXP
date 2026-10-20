
#Results: negligible improvement, not worth the effort
#most improved 24 yr model was loons (they passed hake for lowest aic) but the sign was still positive so it doesn't solve my problem.
#hake was a tiny bit improved (r2 from 3.73 to 4.1) if you add krill and marketsquid. but probably not worth it. daic=

#19yrs: big improvement for birds and fur seals but all the betas are still positive. so it doesn't solve my mechanism problem. not a single sign was flipped
#hake only marginally improved (again) daic = -0.5
# ==============================================================================
# Script: 10c_run_ncc_altprey_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Execute NCC Alternative Prey Attenuation Sweep targeting CPUE (CPUE ~ Peff).
#          Enforces predator self-exclusion and tests additive forage pools.
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

# Load raw dataset and lookup table
raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

guild_pred_lookup <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)
target_cpue       <- "X07_DFA_cpue_IntSprJunHW"

# Filter lookup table strictly for NCC Predators
ncc_lookup <- guild_pred_lookup %>% filter(region == "NCC")

# Standard 0-1 min-max scaling function
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# ------------------------------------------------------------------------------
# STEP 1: ALTERNATIVE PREY SWEEP FUNCTION (WITH SELF-EXCLUSION)
# ------------------------------------------------------------------------------
test_altprey_attenuation_ncc <- function(pred_col, df_base, max_combo_size = 2) {
  
  if (is.na(pred_col) || pred_col == "" || !pred_col %in% names(df_base)) return(NULL)
  
  P_norm <- df_base[[pred_col]]
  
  row_data <- ncc_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  # Extract candidate prey variables (Excluding spatial overlap AND the predator itself!)
  all_cands   <- unlist(strsplit(row_data$all_candidates_str, ", "))
  spatial_col <- trimws(as.character(row_data$spatial_overlap_col))
  
  prey_candidates <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
  prey_candidates <- intersect(prey_candidates, names(df_base))
  
  # --- Raw Baseline Model (CPUE ~ Raw Predator) ---
  base_fit <- tryCatch(
    sem(paste0(target_cpue, " ~ ", pred_col), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), 
    error = function(e) NULL
  )
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- tryCatch(parameterEstimates(base_fit, standardized = TRUE), error = function(e) NULL)
  if (is.null(base_pe) || !"rhs" %in% names(base_pe)) return(NULL)
  
  base_pe_sub <- base_pe %>% filter(op == "~" & rhs == pred_col)
  if (nrow(base_pe_sub) == 0) return(NULL)
  
  base_est <- as.numeric(base_pe_sub$est[1])
  base_std <- as.numeric(base_pe_sub$std.all[1])
  base_p   <- as.numeric(base_pe_sub$pvalue[1])
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_cpue]])
  
  if (length(prey_candidates) == 0) {
    return(tibble(
      predator = pred_col, base_beta = base_est, base_beta_std = base_std,
      base_pvalue = base_p, base_aic = base_aic, base_r2 = base_r2,
      abs_prey = "None", abs_weight_w = NA_real_, abs_beta_peff = base_est,
      abs_beta_std = base_std, abs_pvalue = base_p, abs_final_aic = base_aic,
      abs_delta_aic = 0, abs_final_r2 = base_r2,
      Q1_aic_improved = FALSE, Q2_beta_is_negative = (base_est < 0),
      Q3_sign_flipped_pos_to_neg = FALSE, Q4_impact_increased = FALSE,
      neg_prey = ifelse(base_est < 0, "None", NA_character_), neg_weight_w = NA_real_,
      neg_beta_peff = ifelse(base_est < 0, base_est, NA_real_),
      neg_pvalue = ifelse(base_est < 0, base_p, NA_real_),
      neg_final_aic = ifelse(base_est < 0, base_aic, NA_real_),
      neg_delta_aic = ifelse(base_est < 0, 0, NA_real_),
      neg_final_r2 = ifelse(base_est < 0, base_r2, NA_real_)
    ))
  }
  
  w_grid <- seq(0.1, 1.0, by = 0.1)
  
  best_abs     <- NULL
  best_abs_aic <- Inf
  
  best_neg     <- NULL
  best_neg_aic <- Inf
  
  # --- Loop over Candidate Forage Combinations ---
  max_k <- min(length(prey_candidates), max_combo_size)
  
  for (k in 1:max_k) {
    combos <- combn(prey_candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      # Build additive forage pool: sum raw candidate columns and re-scale to [0, 1]
      if (length(cb) == 1) {
        forage_pool <- df_base[[cb[1]]]
      } else {
        raw_sum <- rowSums(df_base[, cb, drop = FALSE], na.rm = TRUE)
        forage_pool <- min_max_scale(raw_sum)
      }
      
      for (w in w_grid) {
        # Attenuation Gating Equation: P_eff = P_norm * (1 - w * forage_pool)
        gate_factor <- 1 - w * forage_pool
        P_eff       <- P_norm * gate_factor
        
        if (is.na(sd(P_eff)) || sd(P_eff) < 1e-6) next
        
        df_base$P_alt_eff <- P_eff
        
        fit <- tryCatch(
          sem(paste0(target_cpue, " ~ P_alt_eff"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), 
          error = function(e) NULL
        )
        
        if (!is.null(fit) && lavInspect(fit, "converged")) {
          pe <- tryCatch(parameterEstimates(fit, standardized = TRUE), error = function(e) NULL)
          
          if (!is.null(pe) && "rhs" %in% names(pe)) {
            pe_sub <- pe %>% filter(op == "~" & rhs == "P_alt_eff")
            
            if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
              curr_aic <- as.numeric(fitMeasures(fit, "aic"))
              curr_est <- as.numeric(pe_sub$est[1])
              curr_std <- as.numeric(pe_sub$std.all[1])
              curr_p   <- as.numeric(pe_sub$pvalue[1])
              curr_r2  <- as.numeric(inspect(fit, "r2")[[target_cpue]])
              
              res_obj <- list(
                prey_combo = paste(cb, collapse = " + "),
                weight_w   = w,
                aic        = curr_aic,
                delta_aic  = curr_aic - base_aic,
                beta_peff  = curr_est,
                beta_std   = curr_std,
                pvalue     = curr_p,
                r2         = curr_r2
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
  
  # --- Safe Unpacking of Best ABS Model ---
  if (!is.null(best_abs)) {
    abs_prey  <- best_abs$prey_combo
    abs_w     <- best_abs$weight_w
    abs_beta  <- best_abs$beta_peff
    abs_std   <- best_abs$beta_std
    abs_p     <- best_abs$pvalue
    abs_aic   <- best_abs$aic
    abs_daic  <- best_abs$delta_aic
    abs_r2    <- best_abs$r2
  } else {
    abs_prey <- "None"
    abs_w <- abs_beta <- abs_std <- abs_p <- abs_aic <- abs_daic <- abs_r2 <- NA_real_
  }
  
  # --- Safe Unpacking of Best NEG Model ---
  if (!is.null(best_neg)) {
    neg_prey  <- best_neg$prey_combo
    neg_w     <- best_neg$weight_w
    neg_beta  <- best_neg$beta_peff
    neg_p     <- best_neg$pvalue
    neg_aic   <- best_neg$aic
    neg_daic  <- best_neg$delta_aic
    neg_r2    <- best_neg$r2
  } else {
    neg_prey <- "None"
    neg_w <- neg_beta <- neg_p <- neg_aic <- neg_daic <- neg_r2 <- NA_real_
  }
  
  # Calculate 4 Diagnostic Evaluators
  abs_aic_improved <- if (!is.na(abs_daic)) (abs_daic < 0) else FALSE
  abs_is_negative  <- if (!is.na(abs_beta)) (abs_beta < 0) else FALSE
  abs_sign_flipped <- if (!is.na(abs_beta)) (base_est > 0 && abs_beta < 0) else FALSE
  abs_more_neg     <- if (!is.na(abs_beta)) (abs_beta < base_est) else FALSE
  
  return(tibble(
    predator            = pred_col,
    base_beta           = base_est,
    base_beta_std       = base_std,
    base_pvalue         = base_p,
    base_aic            = base_aic,
    base_r2             = base_r2,
    
    # Absolute Best AIC Model Results
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
# STEP 2: EXECUTE ATTENUATION SWEEP WITH MONITORING
# ------------------------------------------------------------------------------
cat(sprintf("Executing NCC Alternative Prey Attenuation Sweep on %s dataset...\n", years_included))

ncc_predators  <- ncc_lookup$pred_data_col
sweep_results  <- list()

for (i in seq_along(ncc_predators)) {
  pred_col <- ncc_predators[i]
  cat(sprintf("[%d/%d] Testing alt prey attenuation for: %s | Time: %s\n", 
              i, length(ncc_predators), pred_col, format(Sys.time(), "%H:%M:%S")))
  
  if (!is.na(pred_col) && pred_col %in% names(minmax_base)) {
    res_df <- test_altprey_attenuation_ncc(pred_col = pred_col, df_base = minmax_base, max_combo_size = 2)
    if (!is.null(res_df)) {
      sweep_results[[length(sweep_results) + 1]] <- res_df
      
      print(
        res_df %>% select(predator, base_beta, abs_prey, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg)
      )
    }
  }
}

ncc_altprey_table <- bind_rows(sweep_results) %>% arrange(abs_final_aic)

# Print Summary Table
cat("\n==============================================================================\n")
cat(sprintf("   NCC ALTERNATIVE PREY SWEEP RESULTS (CPUE ~ Peff) - %s DATASET              \n", years_included))
cat("==============================================================================\n\n")

print(
  ncc_altprey_table %>% 
    select(predator,base_aic,abs_final_aic,base_r2,abs_final_r2, base_beta, abs_prey, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg),
  n = 50
)

# Export Results
out_csv <- file.path(data_out_dir, paste0("Goal3_NCC_AltPrey_Sweep_Results_", years_included, ".csv"))
write_csv(ncc_altprey_table, out_csv)

cat(sprintf("\nSweep results successfully saved to: %s\n", out_csv))