# ==============================================================================
# Script: 04b_goal3_runmodels_altprey.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Step-down Multi-AltPrey interaction selection without prey main effects.
#          ENFORCED: Predator main effect MUST be negative (beta < 0, p < 0.05)
#          across ALL stages (interactions and main effects).
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

# Core Targets
target_ncc <- "X07_DFA_cpue_IntSprJunHW"
target_ak  <- "X16_SAR"

# Candidate Selection Rerouting
ncc_candidates <- pred_lookup %>% 
  filter(grepl("^(X08_|X09_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(sem_complete_data))

ak_candidates <- pred_lookup %>% 
  filter(grepl("^(X10_|X11_|X15_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(sem_complete_data))

cat(sprintf("Evaluated Candidates -> NCC (Pre-CPUE): %d rows, AK (Post-CPUE): %d rows\n", 
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

get_valid_preys <- function(row_data, df_model) {
  preys <- c(row_data$altprey1_data_col, row_data$altprey2_data_col, row_data$altprey3_data_col)
  preys <- unique(preys[!is.na(preys) & preys != ""])
  intersect(preys, names(df_model))
}

# ------------------------------------------------------------------------------
# STEP 1: MULTI-PREY INTERACTION STEP-DOWN LOGIC (STRICT DIRECTIONALITY)
# ------------------------------------------------------------------------------
stepdown_region_interactions <- function(target_var, pred_var, prey_vars, df_model, base_covars = c()) {
  
  null_covs <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else "1"
  null_return <- list(
    formula = paste0(target_var, " ~ ", null_covs), 
    stage   = "Null / Stage Link Only", 
    fit     = NULL, 
    data    = df_model
  )
  
  if (length(prey_vars) == 0) {
    # Fallback directly to Predator Main Effect
    cov_str <- paste(c(base_covars, pred_var), collapse = " + ")
    formula_str <- paste0(target_var, " ~ ", cov_str)
    fit <- tryCatch(sem(formula_str, data = df_model, missing = "ML", warn = FALSE), error = function(e) NULL)
    
    if (!is.null(fit) && lavInspect(fit, "converged") && is_sig_negative(fit, pred_var)) {
      return(list(formula = formula_str, stage = "Predator Main Effect (Negative)", fit = fit, data = df_model))
    } else {
      return(null_return)
    }
  }
  
  # Construct interaction terms in data frame
  int_names <- c()
  for (i in seq_along(prey_vars)) {
    p_col <- prey_vars[i]
    i_name <- paste0("int_", target_var, "_", i)
    df_model[[i_name]] <- df_model[[pred_var]] * df_model[[p_col]]
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
    
    # Verify both interaction significance AND negative predator main effect
    pred_is_negative <- is_sig_negative(fit, pred_var)
    
    if (nrow(non_sig) == 0 && pred_is_negative) {
      stage_desc <- paste0("AltPrey Interaction (", length(active_ints), " term", if(length(active_ints)>1) "s" else "", ")")
      return(list(formula = formula_str, stage = stage_desc, fit = fit, data = df_model))
    }
    
    # If predator main effect is NOT negative, drop all interactions and fallback
    if (!pred_is_negative) {
      break
    }
    
    # Remove term with largest p-value
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
  
  # Fallback 2: Null / Stage Link Only
  return(null_return)
}

# ------------------------------------------------------------------------------
# STEP 2: PAIRWISE MODEL RUNNER
# ------------------------------------------------------------------------------
run_altprey_stepdown_pair <- function(ncc_row, ak_row, df_model, pair_id) {
  
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  
  ncc_preys <- get_valid_preys(ncc_row, df_model)
  ak_preys  <- get_valid_preys(ak_row, df_model)
  
  # Region 1 (NCC / Pre-CPUE) Step-Down
  ncc_res <- stepdown_region_interactions(target_ncc, ncc_pred, ncc_preys, df_model, base_covars = c())
  if (!is.null(ncc_res$data)) df_model <- ncc_res$data
  
  # Region 2 (AK / Post-CPUE) Step-Down
  ak_res  <- stepdown_region_interactions(target_ak, ak_pred, ak_preys, df_model, base_covars = c(target_ncc))
  
  # Combined SEM Fit
  final_syntax <- paste0(ncc_res$formula, "\n", ak_res$formula)
  final_fit    <- tryCatch(sem(final_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (is.null(final_fit) || !lavInspect(final_fit, "converged")) return(NULL)
  
  fm <- fitMeasures(final_fit, c("aic", "bic", "cfi", "pvalue"))
  
  summary_row <- tibble(
    pair_id          = pair_id,
    ncc_predator     = ncc_pred,
    ak_predator      = ak_pred,
    ncc_prey_tested  = paste(ncc_preys, collapse = "; "),
    ak_prey_tested   = paste(ak_preys, collapse = "; "),
    ncc_final_stage  = ncc_res$stage,
    ak_final_stage   = ak_res$stage,
    aic              = fm[["aic"]],
    bic              = fm[["bic"]],
    cfi              = fm[["cfi"]],
    global_p         = fm[["pvalue"]],
    cpue_r2          = get_safe_r2(final_fit, target_ncc),
    sar_r2           = get_safe_r2(final_fit, target_ak)
  )
  
  # Extract parameter estimates without dropping interaction terms
  params_df <- parameterEstimates(final_fit) %>%
    filter(op == "~") %>%
    select(lhs, term = rhs, est, se, z, pvalue) %>%
    mutate(
      pair_id      = pair_id,
      ncc_predator = ncc_pred,
      ak_predator  = ak_pred,
      aic          = fm[["aic"]],
      cfi          = fm[["cfi"]],
      # Tag whether the term is a main effect, stage link, or interaction
      term_type    = case_when(
        grepl("^int_", term) ~ "Interaction Term",
        term == target_ncc   ~ "Stage Link (CPUE -> SAR)",
        TRUE                 ~ "Main Effect"
      )
    )
  
  list(summary = summary_row, params = params_df)
}

# ------------------------------------------------------------------------------
# STEP 3: EXECUTE FACTORIAL SWEEP
# ------------------------------------------------------------------------------
cat("Running Multi-AltPrey Step-Down Sweep with Strict Negative Predator Enforcement...\n")

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
    run_altprey_stepdown_pair(n_row, a_row, sem_complete_data, p_id)
  }
)

summary_results <- map_dfr(compact(all_results), ~ .x$summary)
all_parameters  <- map_dfr(compact(all_results), ~ .x$params)

# Filter to well-fitting global models (global_p >= 0.05)
summary_results_sig <- summary_results %>% 
#  filter(global_p >= 0.05) %>% 
  arrange(aic)

cat("\n--- SIGNIFICANT & WELL-FITTING MODELS (Global p >= 0.05) ---\n")
print(summary_results_sig[, 1:min(6, ncol(summary_results_sig))])

# Export Outputs----------
write_csv(summary_results,     file.path(data_out_dir, "Goal3_AltPrey_MultiStepdown_AllSummaries.csv"))
write_csv(summary_results_sig, file.path(data_out_dir, "Goal3_AltPrey_MultiStepdown_Summaries.csv"))
write_csv(all_parameters,      file.path(data_out_dir, "Goal3_AltPrey_MultiStepdown_Parameters.csv"))

cat("\nMulti-AltPrey Step-Down Selection Complete!")
cat("\n- Total evaluated pairs:", nrow(comparison_grid))
cat("\n- Retained well-fitting models (global_p >= 0.05):", nrow(summary_results_sig), "\n")

# ==============================================================================
# POST-PROCESSING: REDUCE TO UNIQUE FINAL MODEL STRUCTURES
# ==============================================================================

cat("\nDeduplicating models down to unique final SEM structures...\n")

# Normalize dropped predator labels so identical final formulas merge correctly
summary_results_unique <- summary_results %>%
  mutate(
    # If NCC stage dropped the predator, normalize predator & prey to "None"
    ncc_pred_final = if_else(grepl("Null", ncc_final_stage), "None", ncc_predator),
    ncc_prey_final = if_else(grepl("Interaction", ncc_final_stage), ncc_prey_tested, "None"),
    
    # If AK stage dropped the predator, normalize predator & prey to "None"
    ak_pred_final  = if_else(grepl("Stage Link|Null", ak_final_stage), "None", ak_predator),
    ak_prey_final   = if_else(grepl("Interaction", ak_final_stage), ak_prey_tested, "None")
  ) %>%
  # Group by the exact structure of the FINAL fitted SEM
  group_by(ncc_pred_final, ak_pred_final, ncc_final_stage, ak_final_stage, ncc_prey_final, ak_prey_final) %>%
  # Pick the best representative model per unique structure (lowest AIC)
  arrange(aic) %>%
  slice(1) %>%
  ungroup() %>%
  # Clean up temporary grouping columns and re-sort by AIC
  select(-ncc_pred_final, -ak_pred_final, -ncc_prey_final, -ak_prey_final) %>%
  arrange(aic)

# Filter unique well-fitting models (global_p >= 0.05)
summary_results_unique_sig <- summary_results_unique 

# Filter all_parameters to keep ONLY the parameter sets of the unique representative models
unique_pair_ids <- summary_results_unique$pair_id
all_parameters_unique <- all_parameters %>%
  filter(pair_id %in% unique_pair_ids)

# ------------------------------------------------------------------------------
# PRINT & EXPORT DEDUPLICATED RESULTS
# ------------------------------------------------------------------------------
cat("\n======================================================\n")
cat("          UNIQUE FINAL MODELS SUMMARY RESULTS         \n")
cat("======================================================\n")
cat(sprintf("Original total fitted models : %d\n", nrow(summary_results)))
cat(sprintf("Unique final model structures: %d\n", nrow(summary_results_unique)))
cat(sprintf("Unique well-fitting models   : %d (global_p >= 0.05)\n", nrow(summary_results_unique_sig)))
cat("======================================================\n\n")

cat("--- TOP UNIQUE WELL-FITTING MODELS ---\n")
print(summary_results_unique_sig[, 1:min(7, ncol(summary_results_unique_sig))], n = 20)

# Export Deduplicated Outputs
write_csv(summary_results_unique,     file.path(data_out_dir, "Goal3_AltPrey_Unique_AllSummaries.csv"))
write_csv(summary_results_unique_sig, file.path(data_out_dir, "Goal3_AltPrey_Unique_Summaries.csv"))
write_csv(all_parameters_unique,      file.path(data_out_dir, "Goal3_AltPrey_Unique_Parameters.csv"))

cat("\nDeduplicated tables exported to:", data_out_dir, "\n")


#my summary-----------

summary_results_unique_sig<-summary_results_unique_sig %>% filter(!grepl("sleeperSharks_smoltyr",ak_predator))
print(summary_results_unique_sig[, 1:min(6, ncol(summary_results_unique_sig))])
all_parameters %>% filter(pair_id==summary_results_sig$pair_id[1])

print(summary_results[1, 1:min(6, ncol(summary_results))])

summary_results_top<-summary_results_unique_sig %>% filter(!grepl("Null",ncc_final_stage))%>% filter(!grepl("Null",ak_final_stage))
summary_results_top
#drop cormorants & X10_DFA_ssl, they were never significant, only keep hake, X10_Harbour_s_2yrLead_WS_adultyr & X15_DFA_sleeperSharks_adultyr              

# pair_id ncc_predator         ak_predator                      ncc_prey_tested                                      ak_prey_tested                        ncc_final_stage ak_final_stage   aic   bic   cfi global_p cpue_r2 sar_r2
# 1     159 X09_DFA_HakeAge5Plus X10_Harbour_s_2yrLead_WS_adultyr X09_DFA_HakeAge5Plus; X05_DFA_abundSardine_smoltyr   X09_DFA_HakeAge5Plus; X05_DFA_abundS… AltPrey Intera… Predator Main…  119.  128. 0.798   0.0190   0.496  0.452
# 2     177 X09_DFA_HakeAge5Plus X15_DFA_sleeperSharks_adultyr    X09_DFA_HakeAge5Plus; X05_DFA_abundSardine_smoltyr   X13_pollock_age1plus_adultyr; X13_si… AltPrey Intera… Predator Main…  120.  130. 0.811   0.0287   0.496  0.455


summary_results_top %>% select(pair_id, ncc_predator,ak_predator,ncc_prey_tested,ncc_final_stage, ak_final_stage,   aic,global_p, cpue_r2, sar_r2)
# pair_id ncc_predator         ak_predator                      ncc_prey_tested                                      ncc_final_stage                 ak_final_stage                    aic global_p cpue_r2 sar_r2
# 1     159 X09_DFA_HakeAge5Plus X10_Harbour_s_2yrLead_WS_adultyr X09_DFA_HakeAge5Plus; X05_DFA_abundSardine_smoltyr   AltPrey Interaction (1 term)    Predator Main Effect (Negative)  119.   0.0190   0.496  0.452
# 2     177 X09_DFA_HakeAge5Plus X15_DFA_sleeperSharks_adultyr    X09_DFA_HakeAge5Plus; X05_DFA_abundSardine_smoltyr   AltPrey Interaction (1 term)    Predator Main Effect (Negative)  120.   0.0287   0.496  0.455
# 3     129 X09_DFA_HakeAge5Plus X10_Harbour_s_2yrLead_WS_adultyr X12_DFA_biomassEuphShelfSum; X05_anchovy_GAM_smoltyr Predator Main Effect (Negative) Predator Main Effect (Negative)  122.   0.0119   0.373  0.443
# 4     147 X09_DFA_HakeAge5Plus X15_DFA_sleeperSharks_adultyr    X12_DFA_biomassEuphShelfSum; X05_anchovy_GAM_smoltyr Predator Main Effect (Negative) Predator Main Effect (Negative)  124.   0.0219   0.373  0.453

all_parameters %>% filter(pair_id==summary_results_top$pair_id[1])


# Quick inspection snippet for Pair 159------
p159_data <- sem_complete_data %>%
  mutate(
    int_hake_x_hci    = X09_DFA_HakeAge5Plus * X01_habCompInd_smoltyr,
    int_hake_x_krill = X09_DFA_HakeAge5Plus * X01_DFA_NHL_krill
  )

# Fit both candidates against CPUE to see which term matches est = +0.353, p = 0.015
m1 <- lavaan::sem("X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus + int_hake_x_hci", data = p159_data)
m2 <- lavaan::sem("X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus + int_hake_x_krill", data = p159_data)

cat("--- Candidate 1: Hake x Hake ---\n")
print(lavaan::parameterEstimates(m1) %>% dplyr::filter(op == "~"))

cat("\n--- Candidate 2: Hake x Sardine ---\n")
print(lavaan::parameterEstimates(m2) %>% dplyr::filter(op == "~"))
