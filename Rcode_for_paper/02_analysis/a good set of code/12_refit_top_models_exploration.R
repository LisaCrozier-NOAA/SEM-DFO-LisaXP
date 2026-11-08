# ==============================================================================
# Script: 12_refit_top_models_exploration.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Refit user-selected top Min-Max gated models, calculate detailed 
#          fit statistics (R^2 SAR, R^2 CPUE, CFI, p-value), and export summaries.
# Output: Rcode_for_paper/Routput_for_paper/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & RAW DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_file <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
lookup_file   <- file.path(data_out_dir, "guild_predator_lookup_table.csv")
sem_orig      <- read.csv("data_Lisa/sem_master_data.csv", stringsAsFactors = FALSE)

raw_sem_data  <- read.csv(sem_data_file, stringsAsFactors = FALSE) %>% filter(complete.cases(.))

# Align Sitka Herring columns
raw_sem_data$X13_sitkaHerring_EGoA_smoltyr <- sem_orig$sitkaHerring_EGoA
raw_sem_data$X13_sitkaHerring_EGoA_adultyr <- sem_orig$sitkaHerring_EGoA[c(3:24, 24, 24)]

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"

# Load Guild Lookup Table
guild_pred_lookup <- read_csv(lookup_file, show_col_types = FALSE)

# Scale all numeric indicators EXCEPT targets and Year into [0, 1]
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_complete_data <- raw_sem_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TOP PREDATOR MODELS TO EXPLORE
# ------------------------------------------------------------------------------
# Option A: Specify exact predator columns of interest manually
selected_predators <- c(
  "X15_sablefishBiomass_predAK_adultyr",
  "X10_Harbour_s_2yrLead_WS_adultyr",
  "X15_ArrowtoothFlounderBiomass_predAK_adultyr",
  "X15_DFA_sleeperSharks_adultyr",
  "X11_ssl_seak_pup_pred",
  "X15_spinyDogfishGoA_predAK_smoltyr"
)

# Option B (Optional): Automatically load top N from checkpoint file if available
checkpoint_path <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_CHECKPOINT.csv")
if (file.exists(checkpoint_path)) {
  chk <- read_csv(checkpoint_path, show_col_types = FALSE)
  # Take top 10 best models by AIC
  top_chk <- chk %>% arrange(aic) %>% head(10) %>% pull(predator)
  selected_predators <- unique(c(selected_predators, top_chk))
}

cat(sprintf("Exploring %d selected top-performing predator models...\n", length(selected_predators)))

# ------------------------------------------------------------------------------
# STEP 2: DETAILED REFITTING & R^2 EXTRACTION FUNCTION
# ------------------------------------------------------------------------------
refit_and_evaluate_model <- function(pred_col, df_model) {
  row_data <- guild_pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  P_norm <- df_model[[pred_col]]
  candidates <- unlist(strsplit(row_data$all_candidates_str, ", "))
  candidates <- intersect(candidates[!is.na(candidates) & candidates != ""], names(df_model))
  
  # Base Model Fit
  base_syntax <- paste0(target_sar, " ~ ", target_cpue, " + ", pred_col)
  base_fit    <- tryCatch(sem(base_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_sar]])
  
  best_res <- list(
    predator = pred_col, selected_prey = "None (Unattenuated)", weights = "None",
    raw_aic = base_aic, base_aic = base_aic, delta_aic = 0.0,
    beta = as.numeric(parameterEstimates(base_fit) %>% filter(op == "~" & rhs == pred_col) %>% pull(est)),
    pval = as.numeric(parameterEstimates(base_fit) %>% filter(op == "~" & rhs == pred_col) %>% pull(pvalue)),
    sar_r2 = base_r2, cfi = as.numeric(fitMeasures(base_fit, "cfi")), global_p = as.numeric(fitMeasures(base_fit, "pvalue"))
  )
  
  if (length(candidates) == 0) return(best_res)
  
  Anorm_mat <- as.matrix(df_model[, candidates, drop = FALSE])
  w_grid    <- seq(0.1, 1.0, by = 0.1)
  
  # Sweep 1-prey and 2-prey joint gates
  for (k in 1:min(length(candidates), 2)) {
    combos <- combn(candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      w_combo_grid <- expand.grid(rep(list(w_grid), k))
      
      for (g_idx in 1:nrow(w_combo_grid)) {
        w_vals <- as.numeric(w_combo_grid[g_idx, ])
        w1_gate  <- 1 - w_vals[1] * Anorm_mat[, cb[1]]
        P_stage1 <- P_norm * w1_gate
        
        P_joint <- if (k > 1) P_norm * (w1_gate * (1 - w_vals[2] * Anorm_mat[, cb[2]])) else P_stage1
        
        df_model$P_joint_eff <- P_joint
        
        fit <- tryCatch(
          sem(paste0(target_sar, " ~ ", target_cpue, " + P_joint_eff"), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
          error = function(e) NULL
        )
        
        if (!is.null(fit) && lavInspect(fit, "converged")) {
          pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_joint_eff")
          
          if (nrow(pe) > 0 && !is.na(pe$pvalue[1]) && pe$est[1] < 0) {
            curr_aic <- as.numeric(fitMeasures(fit, "aic"))
            
            if (curr_aic < best_res$raw_aic) {
              r2_vals <- inspect(fit, "r2")
              
              best_res <- list(
                predator      = pred_col,
                selected_prey = paste(cb, collapse = " + "),
                weights       = paste(round(w_vals, 2), collapse = "; "),
                raw_aic       = curr_aic,
                base_aic      = base_aic,
                delta_aic     = curr_aic - base_aic,
                beta          = as.numeric(pe$est[1]),
                pval          = as.numeric(pe$pvalue[1]),
                sar_r2        = as.numeric(r2_vals[[target_sar]]),
                cfi           = as.numeric(fitMeasures(fit, "cfi")),
                global_p      = as.numeric(fitMeasures(fit, "pvalue"))
              )
            }
          }
        }
      }
    }
  }
  return(best_res)
}

# ------------------------------------------------------------------------------
# STEP 3: RUN EXPLORATION AND PRINT SUMMARY
# ------------------------------------------------------------------------------
detailed_results <- list()

for (pred in selected_predators) {
  res <- refit_and_evaluate_model(pred, minmax_complete_data)
  if (!is.null(res)) {
    detailed_results[[length(detailed_results) + 1]] <- as_tibble(res)
  }
}

exploration_table <- bind_rows(detailed_results) %>%
  arrange(raw_aic)

cat("\n==============================================================================\n")
cat("      DETAILED FIT EXPLORATION FOR SELECTED TOP MODELS (INCLUDES R^2 SAR)     \n")
cat("==============================================================================\n\n")

print(
  exploration_table %>% 
    select(predator, selected_prey, weights, beta, pval, sar_r2, cfi, raw_aic, delta_aic),
  n = 30
)

write_csv(exploration_table, file.path(tbl_out_dir, "Selected_Top_Models_Detailed_Fit_with_R2.csv"))