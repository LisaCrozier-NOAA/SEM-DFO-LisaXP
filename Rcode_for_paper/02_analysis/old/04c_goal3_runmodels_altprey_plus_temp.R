# ==============================================================================
# Script: 04c_goal3_runmodels_altprey_plus_temp.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Unified Step-Down Interaction Selection including both AltPrey and 
#          Temperature/HCI interactions without main effects of prey/temp.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA
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

# Apply +5.0 domain shift across numeric columns (excluding Year)
sem_complete_data <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))

pred_lookup <- read.csv(lookup_file, stringsAsFactors = FALSE)

# Core Targets and Temperature Drivers
target_ncc <- "X07_DFA_cpue_IntSprJunHW"
target_ak  <- "X16_SAR"

temp_ncc   <- "X01_habCompInd_smoltyr"
temp_ak    <- "X21_sst_egoa_junjulaug_smoltyr"

# Candidate Selection Rerouting
ncc_candidates <- pred_lookup %>% 
  filter(grepl("^(X08_|X09_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(sem_complete_data))

ak_candidates <- pred_lookup %>% 
  filter(grepl("^(X10_|X11_|X15_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(sem_complete_data))

cat(sprintf("Candidates Evaluated -> NCC (Pre-CPUE): %d, AK (Post-CPUE): %d\n", 
            nrow(ncc_candidates), nrow(ak_candidates)))

# Helpers
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

get_valid_preys_and_temp <- function(row_data, temp_col, df_model) {
  preys <- c(row_data$altprey1_data_col, row_data$altprey2_data_col, row_data$altprey3_data_col)
  preys <- unique(preys[!is.na(preys) & preys != ""])
  
  # Append region temperature/HCI as a potential 4th interaction candidate
  all_candidates <- unique(c(preys, temp_col))
  intersect(all_candidates, names(df_model))
}

# ------------------------------------------------------------------------------
# STEP 1: STEP-DOWN INTERACTION LOGIC (PREY + TEMP)
# ------------------------------------------------------------------------------
stepdown_region_interactions <- function(target_var, pred_var, candidate_vars, df_model, base_covars = c()) {
  
  null_covs <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else "1"
  null_return <- list(
    formula = paste0(target_var, " ~ ", null_covs), 
    stage   = "Null / Stage Link Only", 
    fit     = NULL, 
    data    = df_model
  )
  
  if (length(candidate_vars) == 0) {
    cov_str <- paste(c(base_covars, pred_var), collapse = " + ")
    formula_str <- paste0(target_var, " ~ ", cov_str)
    fit <- tryCatch(sem(formula_str, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    
    if (!is.null(fit) && lavInspect(fit, "converged") && is_sig_negative(fit, pred_var)) {
      return(list(formula = formula_str, stage = "Predator Main Effect (Negative)", fit = fit, data = df_model))
    } else {
      return(null_return)
    }
  }
  
  # Create named interaction variables in model data frame
  int_names <- c()
  for (i in seq_along(candidate_vars)) {
    c_col <- candidate_vars[i]
    clean_c_col <- gsub("^X\\d+_", "", c_col)
    i_name <- paste0("int_", pred_var, "_x_", clean_c_col)
    
    df_model[[i_name]] <- df_model[[pred_var]] * df_model[[c_col]]
    int_names <- c(int_names, i_name)
  }
  
  active_ints <- int_names
  
  # Iterative Culling Loop
  while (length(active_ints) > 0) {
    rhs_terms <- c(base_covars, pred_var, active_ints)
    formula_str <- paste0(target_var, " ~ ", paste(rhs_terms, collapse = " + "))
    
    fit <- tryCatch(sem(formula_str, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    
    if (is.null(fit) || !lavInspect(fit, "converged")) {
      active_ints <- active_ints[-length(active_ints)]
      next
    }
    
    pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs %in% active_ints)
    if (nrow(pe) == 0) break
    
    non_sig <- pe %>% filter(is.na(pvalue) | pvalue >= 0.05)
    pred_is_negative <- is_sig_negative(fit, pred_var)
    
    # Check if all remaining interactions are significant AND predator main effect is negative
    if (nrow(non_sig) == 0 && pred_is_negative) {
      stage_desc <- paste0("Significant Interaction (", length(active_ints), " term", if(length(active_ints)>1) "s" else "", ")")
      return(list(formula = formula_str, stage = stage_desc, fit = fit, data = df_model))
    }
    
    # Drop interactions if predator main effect fails directionality
    if (!pred_is_negative) {
      break
    }
    
    # Remove least significant interaction term (highest p-value)
    worst_term <- non_sig %>% arrange(desc(pvalue)) %>% slice(1) %>% pull(rhs)
    active_ints <- setdiff(active_ints, worst_term)
  }
  
  # Fallback 1: Predator Main Effect Only
  cov_str <- paste(c(base_covars, pred_var), collapse = " + ")
  formula_str <- paste0(target_var, " ~ ", cov_str)
  fit <- tryCatch(sem(formula_str, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (!is.null(fit) && lavInspect(fit, "converged") && is_sig_negative(fit, pred_var)) {
    return(list(formula = formula_str, stage = "Predator Main Effect (Negative)", fit = fit, data = df_model))
  }
  
  # Fallback 2: Null / Stage Link
  return(null_return)
}

# ------------------------------------------------------------------------------
# STEP 2: PAIRWISE MODEL RUNNER
# ------------------------------------------------------------------------------
run_altprey_temp_pair <- function(ncc_row, ak_row, df_model, pair_id) {
  
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  
  ncc_candidates_list <- get_valid_preys_and_temp(ncc_row, temp_ncc, df_model)
  ak_candidates_list  <- get_valid_preys_and_temp(ak_row, temp_ak, df_model)
  
  # Region 1 (NCC / Pre-CPUE)
  ncc_res <- stepdown_region_interactions(target_ncc, ncc_pred, ncc_candidates_list, df_model, base_covars = c())
  if (!is.null(ncc_res$data)) df_model <- ncc_res$data
  
  # Region 2 (AK / Post-CPUE)
  ak_res  <- stepdown_region_interactions(target_ak, ak_pred, ak_candidates_list, df_model, base_covars = c(target_ncc))
  
  # Combined SEM Fit
  final_syntax <- paste0(ncc_res$formula, "\n", ak_res$formula)
  final_fit    <- tryCatch(sem(final_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (is.null(final_fit) || !lavInspect(final_fit, "converged")) return(NULL)
  
  fm <- fitMeasures(final_fit, c("aic", "bic", "cfi", "pvalue"))
  
  summary_row <- tibble(
    pair_id          = pair_id,
    ncc_predator     = ncc_pred,
    ak_predator      = ak_pred,
    ncc_candidates   = paste(ncc_candidates_list, collapse = "; "),
    ak_candidates    = paste(ak_candidates_list, collapse = "; "),
    ncc_final_stage  = ncc_res$stage,
    ak_final_stage   = ak_res$stage,
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
# STEP 3: EXECUTE FACTORIAL SWEEP
# ------------------------------------------------------------------------------
cat("Running Combined AltPrey + Temperature Interaction Sweep...\n")

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
    run_altprey_temp_pair(n_row, a_row, sem_complete_data, p_id)
  }
)

summary_results <- map_dfr(compact(all_results), ~ .x$summary)
all_parameters  <- map_dfr(compact(all_results), ~ .x$params)

# Filter to well-fitting global models (global_p >= 0.05)
summary_results_sig <- summary_results %>% 
  filter(global_p >= 0.05) %>% 
  arrange(aic)

cat("\n--- SIGNIFICANT & WELL-FITTING MODELS (Global p >= 0.05) ---\n")
print(summary_results_sig[, 1:min(6, ncol(summary_results_sig))])

# Export Outputs
write_csv(summary_results,     file.path(data_out_dir, "Goal3_AltPreyTemp_AllSummaries.csv"))
write_csv(summary_results_sig, file.path(data_out_dir, "Goal3_AltPreyTemp_Summaries.csv"))
write_csv(all_parameters,      file.path(data_out_dir, "Goal3_AltPreyTemp_Parameters.csv"))

cat("\nCombined Selection Sweep Complete!")
cat("\n- Total evaluated pairs:", nrow(comparison_grid))
cat("\n- Retained well-fitting models (global_p >= 0.05):", nrow(summary_results_sig), "\n")