# ==============================================================================
# Script: 05_ak_predators_logistic_switching_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Logistic switching/attenuation grid sweep across all Alaska predators 
#          (Guilds X10, X11, X15) with multiple prey items and temperature.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS & DATA SETUP
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_file <- file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv")
lookup_file   <- file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv")

if (!file.exists(sem_data_file)) stop("Missing input file: ", sem_data_file)
if (!file.exists(lookup_file))   stop("Missing input file: ", lookup_file)

sem_data_noNA     <- read.csv(sem_data_file, stringsAsFactors = FALSE)
sem_complete_data <- sem_data_noNA %>% filter(complete.cases(.))

cat(sprintf("Loaded dataset with %d rows and %d columns.\n", nrow(sem_complete_data), ncol(sem_complete_data)))

# Apply +5.0 domain shift across numeric columns (excluding Year)
sem_complete_data <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))

pred_lookup <- read.csv(lookup_file, stringsAsFactors = FALSE)

# Targets & Drivers
target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
temp_ak     <- "X21_sst_egoa_junjulaug_smoltyr"

# Filter to Alaska Predators (X10, X11, X15)
ak_candidates <- pred_lookup %>% 
  filter(grepl("^(X10_|X11_|X15_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(sem_complete_data)) %>%
  distinct(pred_data_col, .keep_all = TRUE)

cat(sprintf("Found %d unique Alaska predators for evaluation.\n", nrow(ak_candidates)))

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

logistic_gate <- function(a_norm, gamma) {
  1 / (1 + exp(-gamma * (a_norm - 0.5)))
}

get_valid_candidates <- function(row_data, df_model) {
  preys <- c(row_data$altprey1_data_col, row_data$altprey2_data_col, row_data$altprey3_data_col)
  preys <- unique(preys[!is.na(preys) & preys != ""])
  all_candidates <- unique(c(preys, temp_ak))
  intersect(all_candidates, names(df_model))
}

# ------------------------------------------------------------------------------
# STEP 1: LOGISTIC SWITCHING SINGLE-PREY/TEMP EVALUATOR
# ------------------------------------------------------------------------------
evaluate_logistic_switching <- function(pred_col, candidate_col, df_model) {
  
  P <- df_model[[pred_col]]
  A <- df_model[[candidate_col]]
  A_norm <- min_max_scale(A)
  
  # Grid bounds: w in [0.0, 1.0], gamma in [1, 10]
  w_grid     <- seq(0.1, 1.0, by = 0.1)
  gamma_grid <- c(1, 2, 4, 6, 8, 10)
  
  best_fit <- NULL
  best_aic <- Inf
  best_params <- list(w = 0, gamma = 0, est = NA_real_, pval = NA_real_)
  
  # Fit baseline un-moderated model first
  base_syntax <- paste0(target_sar, " ~ ", target_cpue, " + ", pred_col)
  base_sem    <- tryCatch(sem(base_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (!is.null(base_sem) && lavInspect(base_sem, "converged")) {
    best_aic <- fitMeasures(base_sem, "aic")
    pe <- parameterEstimates(base_sem) %>% filter(op == "~" & rhs == pred_col)
    if (nrow(pe) > 0) {
      best_params <- list(w = 0, gamma = 0, est = pe$est[1], pval = pe$pvalue[1])
    }
  }
  
  # Grid search over w and gamma
  for (w_val in w_grid) {
    for (g_val in gamma_grid) {
      
      # Construct Effective Predator Impact P_eff
      weight_mod <- 1 - (w_val * logistic_gate(A_norm, g_val))
      df_model$P_eff <- P * weight_mod
      
      sem_syntax <- paste0(target_sar, " ~ ", target_cpue, " + P_eff")
      fit <- tryCatch(sem(sem_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
      
      if (!is.null(fit) && lavInspect(fit, "converged")) {
        pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_eff")
        
        if (nrow(pe) > 0 && !is.na(pe$pvalue[1]) && pe$est[1] < 0) { # Must be negative predation
          current_aic <- fitMeasures(fit, "aic")
          
          if (current_aic < best_aic) {
            best_aic <- current_aic
            best_fit <- fit
            best_params <- list(w = w_val, gamma = g_val, est = pe$est[1], pval = pe$pvalue[1])
          }
        }
      }
    }
  }
  
  return(list(aic = best_aic, params = best_params))
}

# ------------------------------------------------------------------------------
# STEP 2: MULTI-CANDIDATE LOGISTIC SWEEP PER ALASKA PREDATOR
# ------------------------------------------------------------------------------
cat("Executing Logistic Switching Grid Sweep across Alaska Predators...\n")

sweep_results <- list()

for (i in seq_len(nrow(ak_candidates))) {
  ak_row   <- ak_candidates[i, ]
  pred_col <- ak_row$pred_data_col
  candidates <- get_valid_candidates(ak_row, sem_complete_data)
  
  cat(sprintf("\nProcessing Predator [%d/%d]: %s (Candidates: %s)\n", 
              i, nrow(ak_candidates), pred_col, paste(candidates, collapse = ", ")))
  
  # 1. Baseline Model (No Moderator)
  base_res <- evaluate_logistic_switching(pred_col, candidates[1], sem_complete_data)
  
  sweep_results[[length(sweep_results) + 1]] <- tibble(
    predator      = pred_col,
    moderator     = "None (Unattenuated)",
    best_w        = 0.0,
    best_gamma    = 0.0,
    pred_est      = base_res$params$est,
    pvalue        = base_res$params$pval,
    aic           = base_res$aic,
    delta_aic     = 0.0
  )
  
  # 2. Sweep over all alternate prey & temperature candidates
  for (cand_col in candidates) {
    res <- evaluate_logistic_switching(pred_col, cand_col, sem_complete_data)
    
    sweep_results[[length(sweep_results) + 1]] <- tibble(
      predator      = pred_col,
      moderator     = cand_col,
      best_w        = res$params$w,
      best_gamma    = res$params$gamma,
      pred_est      = res$params$est,
      pvalue        = res$params$pval,
      aic           = res$aic,
      delta_aic     = NA_real_
    )
  }
}

# ------------------------------------------------------------------------------
# STEP 3: CONSOLIDATE & EXPORT RESULTS
# ------------------------------------------------------------------------------
all_sweep_df <- map_dfr(sweep_results, ~ .x) %>%
  group_by(predator) %>%
  mutate(
    base_aic  = min(aic[moderator == "None (Unattenuated)"], na.rm = TRUE),
    delta_aic = aic - base_aic
  ) %>%
  ungroup() %>%
  arrange(predator, aic)

# Top Moderated Models (where attenuation improved AIC: delta_aic < 0)
improved_models <- all_sweep_df %>%
  filter(moderator != "None (Unattenuated)" & delta_aic < -2.0 & pvalue < 0.05) %>%
  arrange(delta_aic)

cat("\n======================================================\n")
cat("       LOGISTIC SWITCHING SWEEP RESULTS SUMMARY        \n")
cat("======================================================\n")
cat(sprintf("Total Predator-Moderator Configurations Tested: %d\n", nrow(all_sweep_df)))
cat(sprintf("Models with Substantial AIC Improvement (delta_aic < -2.0): %d\n", nrow(improved_models)))
cat("======================================================\n\n")

if (nrow(improved_models) > 0) {
  cat("--- TOP MODERATED ALASKA PREDATOR MODELS ---\n")
  print(improved_models, n = 20)
} else {
  cat("No moderated configurations achieved delta_aic < -2.0 over baseline unattenuated predation.\n")
}

# Export Output Data
out_all_file <- file.path(data_out_dir, "Goal3_AK_Predators_Logistic_Switching_All.csv")
out_top_file <- file.path(data_out_dir, "Goal3_AK_Predators_Logistic_Switching_Improved.csv")

write_csv(all_sweep_df, out_all_file)
write_csv(improved_models, out_top_file)

improved_models$moderator


cat("\nFull results exported to:", data_out_dir, "\n")