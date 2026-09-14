



# ==============================================================================
# Script: 01_goal3_parameter_extraction_and_comparison.R
# Directory: Rcode_for_paper/04_goal3_interaction_analysis/
# Purpose: Goal 3 - Complete Parameter Extraction & Non-Linear Model Benchmarking
# Output: Rcode_for_paper/Routput_for_paper/ (data and tables)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(gt)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DIRECTORIES------
# ------------------------------------------------------------------------------
proj_dir     <- getwd()
master_out   <- file.path(proj_dir, "Rcode_for_paper", "Routput_for_paper")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_path <- "copilot/outputs_9/sem_altprey_data_complete_1998_2021.csv"
ak_yr_path    <- "C:/Users/Lisa.Crozier C drive from work/Marine survival/SEM-DFO-LisaXP/data_Lisa/ak_yr.csv"
if (!file.exists(ak_yr_path)) ak_yr_path <- "data_Lisa/ak_yr.csv"

lookup_path   <- "copilot/outputs_9/all_pred_dfa_altprey.csv"

cat("Preparing standardized comparison dataset for Goal 3...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD & PREPARE UNIFIED DATASET (+5.0 SHIFT)------
# ------------------------------------------------------------------------------
sem_raw <- read.csv(sem_data_path, row.names = NULL)
names(sem_raw) <- tolower(names(sem_raw))

ak_yr_df <- read.csv(ak_yr_path) %>%
  clean_names() %>%
  select(year, sst_egoa = sst_egoa_coastwatch_junjulaug)

sem_joined <- sem_raw %>% left_join(ak_yr_df, by = "year")

# Apply +5.0 domain shift across numeric columns (excluding year)
sem_complete_data <- sem_joined %>%
  mutate(across(where(is.numeric) & !matches("year"), ~ .x + 5.0))

pred_lookup <- read.csv(lookup_path, row.names = NULL) %>%
  clean_names() %>%
  mutate(across(everything(), tolower))

target_ncc <- "x07_dfa_cpue_intsprjunhw"
target_ak  <- "x16_sar"

ncc_candidates <- pred_lookup %>% filter(region == "ncc" & grepl("^(x08_|x09_|x11_)", pred_data_col))
ak_candidates  <- pred_lookup %>% filter(region == "ak" | grepl("^x10_", pred_data_col))

get_safe_r2 <- function(fit, target_var) {
  r2_vals <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
  if (!is.null(r2_vals) && target_var %in% names(r2_vals)) {
    return(as.numeric(r2_vals[[target_var]]))
  }
  return(NA_real_)
}

# ------------------------------------------------------------------------------
# STEP 2: FIT & EXTRACT ALL PARAMETERS ACROSS COMPARISON STRUCTURES------
# ------------------------------------------------------------------------------
evaluate_comparison_trio <- function(ncc_row, ak_row, df_model, pair_id) {
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  
  ncc_prey <- ncc_row$altprey1_data_col
  ak_prey  <- ak_row$altprey1_data_col
  
  # Thermal Interaction Setup
  df_model$ncc_sst_int <- df_model[[ncc_pred]] * df_model$sst_egoa
  df_model$ak_sst_int  <- df_model[[ak_pred]]  * df_model$sst_egoa
  
  # AltPrey Interaction Setup
  if (!is.na(ncc_prey) && ncc_prey %in% names(df_model)) {
    df_model$ncc_alt_int <- df_model[[ncc_pred]] * df_model[[ncc_prey]]
  } else {
    df_model$ncc_alt_int <- NA
  }
  
  if (!is.na(ak_prey) && ak_prey %in% names(df_model)) {
    df_model$ak_alt_int <- df_model[[ak_pred]] * df_model[[ak_prey]]
  } else {
    df_model$ak_alt_int <- NA
  }
  
  # ----------------------------------------------------------------------------
  # SYNTAX DEFINITIONS----------
  # ----------------------------------------------------------------------------
  # 1. Baseline DAG 1A: CPUE ~ PredNCC | SAR ~ CPUE + PreyAK + PredAK-------
  syn_base <- paste0(
    target_ncc, " ~ ", ncc_pred, "\n",
    target_ak,  " ~ ", target_ncc, " + ", ak_pred, 
    if (!is.na(ak_prey) && ak_prey %in% names(df_model)) paste0(" + ", ak_prey) else ""
  )
  
  # 2. Thermal Interaction DAG 1D---------
  syn_thermal <- paste0(
    target_ncc, " ~ ", ncc_pred, " + ncc_sst_int\n",
    target_ak,  " ~ ", target_ncc, " + ", ak_pred, " + ak_sst_int"
  )
  
  # 3. AltPrey Interaction DAG 1D---------
  syn_altprey <- if (!is.na(df_model$ncc_alt_int[1]) && !is.na(df_model$ak_alt_int[1])) {
    paste0(
      target_ncc, " ~ ", ncc_pred, " + ncc_alt_int\n",
      target_ak,  " ~ ", target_ncc, " + ", ak_pred, " + ak_alt_int"
    )
  } else NULL
  
  # Helper to fit, summarize metrics, and preserve ALL parameters
  fit_and_extract_all <- function(syn, model_type) {
    if (is.null(syn)) return(NULL)
    fit <- tryCatch(
      sem(syn, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
    
    fm <- fitMeasures(fit, c("aic", "bic", "cfi", "pvalue"))
    
    summary_row <- tibble(
      pair_id      = pair_id,
      model_type   = model_type,
      ncc_predator = ncc_pred,
      ak_predator  = ak_pred,
      aic          = fm[["aic"]],
      bic          = fm[["bic"]],
      cfi          = fm[["cfi"]],
      global_p     = fm[["pvalue"]],
      cpue_r2      = get_safe_r2(fit, target_ncc),
      sar_r2       = get_safe_r2(fit, target_ak)
    )
    
    params_df <- parameterEstimates(fit) %>%
      filter(op == "~") %>%
      select(lhs, term = rhs, est, se, z, pvalue) %>%
      mutate(
        pair_id      = pair_id,
        model_type   = model_type,
        ncc_predator = ncc_pred,
        ak_predator  = ak_pred,
        aic          = fm[["aic"]],
        cfi          = fm[["cfi"]],
        term_type    = case_when(
          term %in% c(ncc_pred, ak_pred) ~ "Main Predator Effect",
          term %in% c("ncc_sst_int", "ak_sst_int") ~ "Thermal Interaction (Pred*SST)",
          term %in% c("ncc_alt_int", "ak_alt_int") ~ "AltPrey Interaction (Pred*Prey)",
          term == target_ncc ~ "Stage Link (CPUE -> SAR)",
          TRUE ~ "Prey Effect"
        )
      )
    
    list(summary = summary_row, params = params_df)
  }
  
  res_base    <- fit_and_extract_all(syn_base, "Baseline DAG 1A (Main Effects)")
  res_thermal <- fit_and_extract_all(syn_thermal, "Thermal Moderation (Pred*SST)")
  res_altprey <- fit_and_extract_all(syn_altprey, "AltPrey Moderation (Pred*Prey)")
  
  list(
    summaries = bind_rows(res_base$summary, res_thermal$summary, res_altprey$summary),
    parameters = bind_rows(res_base$params, res_thermal$params, res_altprey$params)
  )
}

# ------------------------------------------------------------------------------
# STEP 3: RUN FACTORIAL SWEEP------------
# ------------------------------------------------------------------------------
cat("Running Goal 3 Factorial Sweep across all model variants...\n")

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
    evaluate_comparison_trio(n_row, a_row, sem_complete_data, p_id)
  }
)

# Extract and bind complete tables
summary_results <- map_dfr(all_results, ~ .x$summaries) %>%
  group_by(pair_id) %>%
  mutate(delta_aic_in_pair = round(aic - min(aic, na.rm = TRUE), 2)) %>%
  ungroup() %>%
  arrange(pair_id, aic)

all_parameters_df <- map_dfr(all_results, ~ .x$parameters)

# ------------------------------------------------------------------------------
# STEP 4: EXPORT FULL PARAMETERS & MODEL SUMMARIES---------
# ------------------------------------------------------------------------------
write_csv(summary_results,   file.path(data_out_dir, "Goal3_Model_Summaries_Comparison.csv"))
write_csv(all_parameters_df, file.path(data_out_dir, "Goal3_All_Parameter_Estimates.csv"))

cat("Saved full parameter estimates to: Goal3_All_Parameter_Estimates.csv\n")

# # ------------------------------------------------------------------------------
# # STEP 5: RENDER MANUSCRIPT TABLE 5 (GT)
# # ------------------------------------------------------------------------------
# top_comparisons <- summary_results %>%
#   group_by(model_type) %>%
#   slice_min(order_by = aic, n = 5) %>%
#   ungroup() %>%
#   arrange(aic) %>%
#   mutate(delta_aic_overall = round(aic - min(aic), 2))
# 
# table5_gt <- top_comparisons %>%
#   select(model_type, ncc_predator, ak_predator, aic, delta_aic_overall, cfi, cpue_r2, sar_r2) %>%
#   gt(groupname_col = "model_type") %>%
#   tab_header(
#     title = md("**Table 5. Goal 3 Model Evaluation: Baseline DAG 1A vs. Non-Linear Interaction Structures**"),
#     subtitle = "Standardized evaluation on complete 1998–2021 observations"
#   ) %>%
#   cols_label(
#     ncc_predator      = md("**NCC Predator**"),
#     ak_predator       = md("**AK Predator**"),
#     aic               = md("**AIC**"),
#     delta_aic_overall = md("**ΔAIC**"),
#     cfi               = md("**CFI**"),
#     cpue_r2           = md("**R² (CPUE)**"),
#     sar_r2            = md("**R² (SAR)**")
#   ) %>%
#   cols_align(align = "center", columns = c(aic, delta_aic_overall, cfi, cpue_r2, sar_r2)) %>%
#   fmt_number(columns = c(aic, delta_aic_overall, cfi, cpue_r2, sar_r2), decimals = 2) %>%
#   tab_options(
#     table.font.size = px(9.5),
#     heading.title.font.size = px(11.5),
#     column_labels.font.weight = "bold",
#     row_group.font.weight = "bold"
#   )
# 
# # Save Outputs
# gtsave(table5_gt, file.path(tbl_out_dir, "Table5_Model_Comparison_Baseline_vs_Interactions.html"))
# gtsave(
#   data     = table5_gt,
#   filename = file.path(tbl_out_dir, "Table5_Model_Comparison_Baseline_vs_Interactions.png"),
#   vwidth   = 1100,
#   vheight  = 600
# )
# 
# message("Goal 3 Pipeline Complete! Table 5 and full parameters exported to: ", master_out)
# 
# 
# #check parameters==============
# library(tidyverse)
# 
# # Load Goal 3 Parameters
# params <- read.csv("Rcode_for_paper/Routput_for_paper/data/Goal3_All_Parameter_Estimates.csv")
# 
# # Filter for the top Thermal Moderation model (pair_id = 1 or matching top predators)
# top_thermal_params <- params %>%
#   filter(
#     model_type == "Thermal Moderation (Pred*SST)",
#     ncc_predator == "x11_dfa_harbour_p_ws",
#     ak_predator == "x10_harbour_s_2yrlead_ws"
#   ) %>%
#   select(lhs, term, term_type, est, se, z, pvalue)
# 
# print(top_thermal_params)
# 
# params %>%
#   filter(pvalue<0.05)
