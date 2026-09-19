# ==============================================================================
# Script: 04_goal3_runmodels_altprey.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Step-down Alternative Prey interaction selection (+5.0 shift)
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA (+5.0 DOMAIN SHIFT)
# ------------------------------------------------------------------------------
proj_dir     <- paste0(getwd(), "/Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_noNA  <- read.csv(file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_short <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_complete_1998_2021.csv"), stringsAsFactors = FALSE)
names(sem_data_short)

# Choose dataset
sem_data <- sem_data_short

# Apply +5.0 domain shift
sem_complete_data <- sem_data %>%
  filter(complete.cases(.)) %>%
  mutate(across(where(is.numeric) & !matches("year", ignore.case = TRUE), ~ .x + 5.0))

pred_lookup <- read.csv(file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv"), stringsAsFactors = FALSE)

target_ncc <- "X07_DFA_cpue_IntSprJunHW"
target_ak  <- "X16_SAR"

# FULL CANDIDATE SELECTION (Including X15 Alaska Predators)
ncc_candidates <- pred_lookup %>% filter(region == "NCC" & grepl("^(X08_|X09_|X10_|X11_)", pred_data_col))
ak_candidates  <- pred_lookup %>% filter(region == "AK"  & grepl("^(X11_|X15_)", pred_data_col))

is_sig <- function(fit, term_name, alpha = 0.05) {
  pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == term_name)
  if (nrow(pe) == 0) return(FALSE)
  return(!is.na(pe$pvalue[1]) && pe$pvalue[1] < alpha)
}

is_sig_negative <- function(fit, term_name, alpha = 0.05) {
  pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == term_name)
  if (nrow(pe) == 0) return(FALSE)
  return(!is.na(pe$pvalue[1]) && pe$pvalue[1] < alpha && pe$est[1] < 0)
}

get_safe_r2 <- function(fit, target_var) {
  r2_vals <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
  if (!is.null(r2_vals) && target_var %in% names(r2_vals)) {
    return(as.numeric(r2_vals[[target_var]]))
  }
  return(NA_real_)
}

# ------------------------------------------------------------------------------
# STEP 1: STEP-DOWN ALTPREY SELECTION LOGIC
# ------------------------------------------------------------------------------
run_altprey_stepdown <- function(ncc_row, ak_row, df_model, pair_id) {
  
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  ncc_prey <- ncc_row$altprey1_data_col
  ak_prey  <- ak_row$altprey1_data_col
  
  req_cols <- c(target_ncc, target_ak, ncc_pred, ak_pred)
  if (!all(req_cols %in% names(df_model))) return(NULL)
  
  has_ncc_prey <- !is.na(ncc_prey) && ncc_prey %in% names(df_model)
  has_ak_prey  <- !is.na(ak_prey)  && ak_prey  %in% names(df_model)
  
  if (has_ncc_prey) df_model$ncc_alt_int <- df_model[[ncc_pred]] * df_model[[ncc_prey]]
  if (has_ak_prey)  df_model$ak_alt_int  <- df_model[[ak_pred]]  * df_model[[ak_prey]]
  
  # ============================================================================
  # REGION 1: NCC STEP-DOWN
  # ============================================================================
  if (has_ncc_prey) {
    # 1A. Full Interaction
    ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred, " + ", ncc_prey, " + ncc_alt_int")
    fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ncc_stage <- "AltPrey Interaction"
    
    if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig(fit_ncc, "ncc_alt_int")) {
      # 1B. Additive Model (Predator + Prey)
      ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred, " + ", ncc_prey)
      fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
      ncc_stage <- "Additive (Pred + Prey)"
      
      if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig_negative(fit_ncc, ncc_pred)) {
        # 1C. Univariate Predator Effect
        ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred)
        fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
        ncc_stage <- "Predator Main Effect (Negative)"
        
        if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig_negative(fit_ncc, ncc_pred)) {
          ncc_formula <- paste0(target_ncc, " ~ 1")
          ncc_stage <- "Null (Dropped)"
        }
      }
    }
  } else {
    ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred)
    fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ncc_stage <- "Predator Main Effect (Negative)"
    
    if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig_negative(fit_ncc, ncc_pred)) {
      ncc_formula <- paste0(target_ncc, " ~ 1")
      ncc_stage <- "Null (Dropped)"
    }
  }
  
  # ============================================================================
  # REGION 2: AK STEP-DOWN
  # ============================================================================
  if (has_ak_prey) {
    # 2A. Full Interaction
    ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred, " + ", ak_prey, " + ak_alt_int")
    fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ak_stage <- "AltPrey Interaction"
    
    if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig(fit_ak, "ak_alt_int")) {
      # 2B. Additive Model
      ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred, " + ", ak_prey)
      fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
      ak_stage <- "Additive (Pred + Prey)"
      
      if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig_negative(fit_ak, ak_pred)) {
        # 2C. Univariate Predator Effect
        ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred)
        fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
        ak_stage <- "Predator Main Effect (Negative)"
        
        if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig_negative(fit_ak, ak_pred)) {
          ak_formula <- paste0(target_ak, " ~ ", target_ncc)
          ak_stage <- "Stage Link Only"
        }
      }
    }
  } else {
    ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred)
    fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ak_stage <- "Predator Main Effect (Negative)"
    
    if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig_negative(fit_ak, ak_pred)) {
      ak_formula <- paste0(target_ak, " ~ ", target_ncc)
      ak_stage <- "Stage Link Only"
    }
  }
  
  # ============================================================================
  # COMBINED FIT & EVALUATION
  # ============================================================================
  final_syntax <- paste0(ncc_formula, "\n", ak_formula)
  final_fit    <- tryCatch(sem(final_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (is.null(final_fit) || !lavInspect(final_fit, "converged")) return(NULL)
  
  fm <- fitMeasures(final_fit, c("aic", "bic", "cfi", "pvalue"))
  
  summary_row <- tibble(
    pair_id          = pair_id,
    ncc_predator     = ncc_pred,
    ak_predator      = ak_pred,
    ncc_altprey      = if (has_ncc_prey) ncc_prey else NA_character_,
    ak_altprey       = if (has_ak_prey)  ak_prey  else NA_character_,
    ncc_final_stage  = ncc_stage,
    ak_final_stage   = ak_stage,
    aic              = fm[["aic"]],
    bic              = fm[["bic"]],
    cfi              = fm[["cfi"]],
    global_p         = fm[["pvalue"]],
    cpue_r2          = get_safe_r2(final_fit, target_ncc),
    sar_r2           = get_safe_r2(final_fit, target_ak)
  )
  
  params_df <- parameterEstimates(final_fit) %>%
    filter(op == "~") %>%
    select(lhs, term = rhs, est, se, z, pvalue) %>%
    mutate(
      pair_id      = pair_id,
      ncc_predator = ncc_pred,
      ak_predator  = ak_pred,
      aic          = fm[["aic"]],
      cfi          = fm[["cfi"]]
    )
  
  list(summary = summary_row, params = params_df)
}

# ------------------------------------------------------------------------------
# STEP 2: EXECUTE ALTPREY SWEEP
# ------------------------------------------------------------------------------
cat("Running Corrected Step-Down AltPrey Interaction Sweep...\n")

comparison_grid <- expand_grid(
  ncc_idx = seq_len(nrow(ncc_candidates)),
  ak_idx  = seq_len(nrow(ak_candidates))
) %>% mutate(pair_id = row_number())

all_results <- purrr::map(
  seq_len(nrow(comparison_grid)),
  function(idx) {
    p_id  <- comparison_grid$pair_id[idx]
    n_row <- ncc_candidates[comparison_grid$ncc_idx[idx], ]
    a_row <- ak_candidates[comparison_grid$ak_idx[idx], ]
    run_altprey_stepdown(n_row, a_row, sem_complete_data, p_id)
  }
)

summary_results <- map_dfr(compact(all_results), ~ .x$summary)
all_parameters  <- map_dfr(compact(all_results), ~ .x$params)

# Retain well-fitting models (global_p >= 0.05)
summary_results_sig <- summary_results %>% 
  filter(global_p >= 0.05) %>% 
  arrange(aic)

print(summary_results_sig[,1:5])
all_parameters %>% filter(pair_id==165)
all_parameters %>% filter(pair_id==112)

# Save Outputs
write_csv(summary_results_sig, file.path(data_out_dir, "Goal3_AltPrey_Stepdown_Summaries.csv"))
write_csv(all_parameters,      file.path(data_out_dir, "Goal3_AltPrey_Stepdown_Parameters.csv"))

cat("\nAltPrey Step-Down Selection Complete!")
cat("\n- Total evaluated pairs:", nrow(comparison_grid))
cat("\n- Retained well-fitting models (global_p >= 0.05):", nrow(summary_results_sig), "\n")