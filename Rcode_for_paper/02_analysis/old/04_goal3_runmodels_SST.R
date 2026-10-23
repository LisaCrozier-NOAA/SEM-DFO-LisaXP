# ==============================================================================
# Script: 04_goal3_runmodels_SST.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Step-down thermal interaction selection with biological direction checks
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

cat(sprintf("Loaded dataset with %d rows and %d columns.\n", nrow(sem_complete_data), ncol(sem_complete_data)))

# Apply +5.0 domain shift across numeric columns (excluding Year)
sem_complete_data <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))

pred_lookup <- read.csv(lookup_file, stringsAsFactors = FALSE)

# --- EXACT TARGET VARIABLES & TEMPERATURE DRIVERS FROM MATRIX ---
target_ncc <- "X07_DFA_cpue_IntSprJunHW"
target_ak  <- "X16_SAR"

temp_ncc   <- "X01_habCompInd_smoltyr"
temp_ak    <- "X21_sst_egoa_junjulaug_smoltyr"

# Verify target variables exist
req_targets <- c(target_ncc, target_ak, temp_ncc, temp_ak)
missing_targets <- setdiff(req_targets, names(sem_complete_data))
if (length(missing_targets) > 0) {
  stop("Missing core target or temperature columns: ", paste(missing_targets, collapse = ", "))
}

# --- CANDIDATE SELECTION (Filtering only candidates present in sem_data) ---
ncc_candidates <- pred_lookup %>% 
  filter(
    region == "NCC" & (
      grepl("^(X08_|X09_|X11_)", pred_data_col) | 
        (grepl("^X10_", pred_data_col) & grepl("smoltyr", pred_data_col, ignore.case = TRUE))
    )
  ) %>%
  filter(pred_data_col %in% names(sem_complete_data))

ak_candidates <- pred_lookup %>% 
  filter(
    region == "AK" & (
      grepl("^X11_", pred_data_col) | 
        grepl("smoltyr", pred_data_col, ignore.case = TRUE)
    )
  ) %>%
  filter(pred_data_col %in% names(sem_complete_data))

cat(sprintf("Evaluated Candidates Present in Data -> NCC: %d, AK: %d\n", nrow(ncc_candidates), nrow(ak_candidates)))

# Helpers
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
# STEP 1: THERMAL INTERACTION & STEP-DOWN LOGIC
# ------------------------------------------------------------------------------
run_thermal_stepdown <- function(ncc_pred, ak_pred, df_model, pair_id) {
  
  req_cols <- c(target_ncc, target_ak, ncc_pred, ak_pred, temp_ncc, temp_ak)
  if (!all(req_cols %in% names(df_model))) return(NULL)
  
  # Construct interaction terms
  df_model <- df_model %>%
    mutate(
      ncc_sst_int = .data[[ncc_pred]] * .data[[temp_ncc]],
      ak_sst_int  = .data[[ak_pred]]  * .data[[temp_ak]]
    )
  
  # ============================================================================
  # REGION 1: NCC STEP-DOWN (Target: CPUE, Temp: temp_ncc)
  # ============================================================================
  # Step 1: Test Full Interaction (Pred + Temp + Pred*Temp)
  ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred, " + ", temp_ncc, " + ncc_sst_int")
  fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
  ncc_stage <- "Thermal Interaction"
  
  # If interaction is NOT significant, skip additive entirely and drop to Predator Main Effect
  if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig(fit_ncc, "ncc_sst_int")) {
    
    # Step 2: Test Predator Main Effect Only
    ncc_formula <- paste0(target_ncc, " ~ ", ncc_pred)
    fit_ncc <- tryCatch(sem(ncc_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ncc_stage <- "Predator Main Effect (Negative)"
    
    # Step 3: Must be significant AND negative (predation effect)
    if (is.null(fit_ncc) || !lavInspect(fit_ncc, "converged") || !is_sig_negative(fit_ncc, ncc_pred)) {
      ncc_formula <- paste0(target_ncc, " ~ 1")
      ncc_stage <- "Null (Dropped)"
    }
  }
  
  # ============================================================================
  # REGION 2: AK STEP-DOWN (Target: SAR, Temp: temp_ak)
  # ============================================================================
  # Step 1: Test Full Interaction (Stage Link + Pred + Temp + Pred*Temp)
  ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred, " + ", temp_ak, " + ak_sst_int")
  fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
  ak_stage <- "Thermal Interaction"
  
  # If interaction is NOT significant, skip additive entirely and drop to Predator Main Effect
  if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig(fit_ak, "ak_sst_int")) {
    
    # Step 2: Test Predator Main Effect Only
    ak_formula <- paste0(target_ak, " ~ ", target_ncc, " + ", ak_pred)
    fit_ak <- tryCatch(sem(ak_formula, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    ak_stage <- "Predator Main Effect (Negative)"
    
    # Step 3: Must be significant AND negative (predation effect)
    if (is.null(fit_ak) || !lavInspect(fit_ak, "converged") || !is_sig_negative(fit_ak, ak_pred)) {
      ak_formula <- paste0(target_ak, " ~ ", target_ncc)
      ak_stage <- "Stage Link Only"
    }
  }
  
  # ============================================================================
  # FINAL COMBINED MODEL FIT
  # ============================================================================
  final_syntax <- paste0(ncc_formula, "\n", ak_formula)
  final_fit    <- tryCatch(sem(final_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (is.null(final_fit) || !lavInspect(final_fit, "converged")) return(NULL)
  
  fm <- fitMeasures(final_fit, c("aic", "bic", "cfi", "pvalue"))
  
  summary_row <- tibble(
    pair_id          = pair_id,
    ncc_predator     = ncc_pred,
    ak_predator      = ak_pred,
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
# STEP 2: EXECUTE FACTORIAL SWEEP
# ------------------------------------------------------------------------------
cat("Running Thermal Step-Down Sweep...\n")

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
    run_thermal_stepdown(n_row$pred_data_col, a_row$pred_data_col, sem_complete_data, p_id)
  }
)

summary_results <- map_dfr(compact(all_results), ~ .x$summary) %>% arrange(aic)
all_parameters  <- map_dfr(compact(all_results), ~ .x$params)

summary_results_sig <- summary_results %>% filter(global_p >= 0.05)

cat("\n--- SIGNIFICANT MODELS (Global p >= 0.05) ---\n")
print(summary_results_sig[, 1:min(5, ncol(summary_results_sig))])

write_csv(summary_results, file.path(data_out_dir, "Goal3_Thermal_Stepdown_Summaries.csv"))
write_csv(all_parameters,  file.path(data_out_dir, "Goal3_Thermal_Stepdown_Parameters.csv"))

cat("\nStep-Down Thermal Model Selection Complete!")
cat("\n- Total evaluated pairs:", nrow(comparison_grid))
cat("\n- Outputs exported to:", data_out_dir, "\n")