
#Results: this script run for both24yrs and 19yrs showed that 
#there is no reason to add hci or krill to general script. 
#top pred were murre, hake and loons, still positive

# ==============================================================================
# Script: 10b_run_ncc_hci_sweep_10212026.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Execute NCC Temperature/HCI Gating Sweep (CPUE ~ Peff) with dynamic
#          min-max lower bounds [a, 1.0] and weight grid search [0.1, 1.0].
#          Answers 4 core diagnostic questions on AIC, sign flips, and impact.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION & DATA LOADING
# ------------------------------------------------------------------------------
years_included <- "24yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

# Load raw unscaled filtered data (prepared prior to 0-1 minmax)
raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)

# Clean empty/unnamed columns
raw_df <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# # Invert HCI so that High Value = High Temperature / High Compression
# if ("X01_habCompInd_smoltyr" %in% names(raw_df)) {
#   raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr
# }
# if ("X01_habCompInd_adultyr" %in% names(raw_df)) {
#   raw_df$X01_habCompInd_adultyr <- -1 * raw_df$X01_habCompInd_adultyr
# }

# Filter complete cases by year switch
if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

guild_pred_lookup <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)
target_cpue       <- "X07_cpue_IntSprJunHW_smoltyr"

# Filter lookup table strictly for NCC Predators
ncc_lookup <- guild_pred_lookup %>% filter(region == "NCC", stage_suffix == "smoltyr", !grepl("DFO",pred_data_col));ncc_lookup

# ------------------------------------------------------------------------------
# STEP 1: DYNAMIC MIN-MAX SCALING FUNCTION
# ------------------------------------------------------------------------------
scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Standard 0-1 minmax for raw predictors and targets
minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# ------------------------------------------------------------------------------
# STEP 2: TEMPERATURE/HCI GATING SWEEP FUNCTION
# ------------------------------------------------------------------------------
test_temperature_gating_ncc <- function(pred_col, df_raw, df_base) {
  
  if (is.na(pred_col) || pred_col == "" || !pred_col %in% names(df_base)) return(NULL)
  
  P_norm <- df_base[[pred_col]]
  
  row_data <- ncc_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  spatial_col <- trimws(as.character(row_data$spatial_overlap_col))
  if (is.na(spatial_col) || spatial_col == "" || !spatial_col %in% names(df_raw)) return(NULL)
  
  raw_hci <- df_raw[[spatial_col]]
  
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
  
  # Grid ranges for dynamic scaling bound (a) and weight (w)
  a_bounds <- seq(0.0, 0.5, by = 0.1) # Min-Max lower range bounds
  w_grid   <- seq(0.1, 1.0, by = 0.1) # Attenuation weight grid
  
  best_abs     <- NULL
  best_abs_aic <- Inf
  
  best_neg     <- NULL
  best_neg_aic <- Inf
  
  # --- Double Grid Search (Lower Bound 'a' x Weight 'w') ---
  for (a in a_bounds) {
    hci_scaled <- scale_custom_bounds(raw_hci, a = a, b = 1.0)
    
    for (w in w_grid) {
      gate_factor <- 1 - w * (1 - hci_scaled)
      P_eff       <- P_norm * gate_factor
      
      if (is.na(sd(P_eff)) || sd(P_eff) < 1e-6) next
      
      df_base$P_temp_eff <- P_eff
      
      fit <- tryCatch(
        sem(paste0(target_cpue, " ~ P_temp_eff"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), 
        error = function(e) NULL
      )
      
      if (!is.null(fit) && lavInspect(fit, "converged")) {
        pe <- tryCatch(
          parameterEstimates(fit, standardized = TRUE),
          error = function(e) NULL
        )
        
        if (!is.null(pe) && "rhs" %in% names(pe)) {
          pe_sub <- pe %>% filter(op == "~" & rhs == "P_temp_eff")
          
          if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
            curr_aic <- as.numeric(fitMeasures(fit, "aic"))
            curr_est <- as.numeric(pe_sub$est[1])
            curr_std <- as.numeric(pe_sub$std.all[1])
            curr_p   <- as.numeric(pe_sub$pvalue[1])
            curr_r2  <- as.numeric(inspect(fit, "r2")[[target_cpue]])
            
            res_obj <- list(
              minmax_a  = a,
              weight_w  = w,
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
  }
  
  # --- Safe Unpacking of Best ABS Model ---
  if (!is.null(best_abs)) {
    abs_a     <- best_abs$minmax_a
    abs_w     <- best_abs$weight_w
    abs_beta  <- best_abs$beta_peff
    abs_std   <- best_abs$beta_std
    abs_p     <- best_abs$pvalue
    abs_aic   <- best_abs$aic
    abs_daic  <- best_abs$delta_aic
    abs_r2    <- best_abs$r2
  } else {
    abs_a <- abs_w <- abs_beta <- abs_std <- abs_p <- abs_aic <- abs_daic <- abs_r2 <- NA_real_
  }
  
  # --- Safe Unpacking of Best NEG Model ---
  if (!is.null(best_neg)) {
    neg_a     <- best_neg$minmax_a
    neg_w     <- best_neg$weight_w
    neg_beta  <- best_neg$beta_peff
    neg_p     <- best_neg$pvalue
    neg_aic   <- best_neg$aic
    neg_daic  <- best_neg$delta_aic
    neg_r2    <- best_neg$r2
  } else {
    neg_a <- neg_w <- neg_beta <- neg_p <- neg_aic <- neg_daic <- neg_r2 <- NA_real_
  }
  
  # Calculate 4 Diagnostic Evaluators
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
    abs_weight_w        = abs_w,
    abs_beta_peff       = abs_beta,
    abs_beta_std        = abs_std,
    abs_pvalue          = abs_p,
    abs_final_aic       = abs_aic,
    abs_delta_aic       = abs_daic,
    abs_final_r2        = abs_r2,
    
    # Diagnostic Evaluators
    Q1_aic_improved     = abs_aic_improved,
    Q2_beta_is_negative = abs_is_negative,
    Q3_sign_flipped_pos_to_neg = abs_sign_flipped,
    Q4_impact_increased = abs_more_neg,
    
    # Best Negative Beta Model Results
    neg_minmax_a        = neg_a,
    neg_weight_w        = neg_w,
    neg_beta_peff       = neg_beta,
    neg_pvalue          = neg_p,
    neg_final_aic       = neg_aic,
    neg_delta_aic       = neg_daic,
    neg_final_r2        = neg_r2
  ))
}

# ------------------------------------------------------------------------------
# STEP 3: EXECUTE SWEEP WITH MONITORING
# ------------------------------------------------------------------------------
cat(sprintf("Executing NCC Temperature Gating Sweep on %s dataset...\n", years_included))

ncc_predators  <- ncc_lookup$pred_data_col
sweep_results  <- list()

for (i in seq_along(ncc_predators)) {
  pred_col <- ncc_predators[i]
  cat(sprintf("\n[%d/%d] Testing temperature gating for: %s | Time: %s\n", 
              i, length(ncc_predators), pred_col, format(Sys.time(), "%H:%M:%S")))
  
  if (!is.na(pred_col) && pred_col %in% names(minmax_base)) {
    res_df <- test_temperature_gating_ncc(pred_col = pred_col, df_raw = filtered_df, df_base = minmax_base)
    if (!is.null(res_df)) {
      sweep_results[[length(sweep_results) + 1]] <- res_df
      
      # Monitoring output
      print(
        res_df %>% select(predator, base_beta, abs_minmax_a, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg)
      )
    }
  }
}

ncc_sweep_table <- bind_rows(sweep_results) %>% arrange(abs_final_aic)

# Print Final Summary Table
cat("\n==============================================================================\n")
cat(sprintf("   NCC TEMPERATURE GATING SWEEP RESULTS (CPUE ~ Peff) - %s DATASET           \n", years_included))
cat("==============================================================================\n\n")

print(
  ncc_sweep_table %>% 
    select(predator,base_aic,abs_final_aic,  base_beta, abs_minmax_a, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg),
  n = 50
)

# Export Results
out_csv <- file.path(data_out_dir, paste0("Goal3_NCC_Temperature_Gating_Sweep_", years_included, ".csv"))
write_csv(ncc_sweep_table, out_csv)

cat(sprintf("\nSweep results successfully saved to: %s\n", out_csv))