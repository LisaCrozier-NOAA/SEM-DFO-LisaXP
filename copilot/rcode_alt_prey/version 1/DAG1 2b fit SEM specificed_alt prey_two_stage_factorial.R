suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_9"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Lookup Table (Lower-case Standard)
# -----------------------------------------------------------------------------
sem_complete_data <- read.csv(file.path(out_dir, "sem_altprey_data_complete_1998_2021.csv"), row.names = NULL)
names(sem_complete_data) <- tolower(names(sem_complete_data))

pred_lookup <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv"), row.names = NULL) %>%
  clean_names() %>%
  mutate(
    pred_data_col     = tolower(pred_data_col),
    altprey1_data_col = tolower(altprey1_data_col),
    altprey2_data_col = tolower(altprey2_data_col),
    altprey3_data_col = tolower(altprey3_data_col)
  )

target_ncc <- "x07_dfa_cpue_intsprjunhw"
target_ak  <- "x16_sar"

# Separate NCC Stage vs. AK Stage Candidate Lists
ncc_candidates <- pred_lookup %>% filter(region == "NCC" & grepl("^(x08_|x09_|x11_)", pred_data_col))
ak_candidates  <- pred_lookup %>% filter(region == "AK" | grepl("^x10_", pred_data_col))

# -----------------------------------------------------------------------------
# 2. Factorial Two-Stage SEM Fitting & Backward Pruning Function
# -----------------------------------------------------------------------------
fit_two_stage_pair <- function(ncc_row, ak_row, df_raw) {
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  
  df_model <- df_raw
  
  # --- Build NCC Interaction Products (ncc_int_1, ncc_int_2, ...) ---
  ncc_preys <- c(ncc_row$altprey1_data_col, ncc_row$altprey2_data_col, ncc_row$altprey3_data_col)
  ncc_preys <- ncc_preys[!is.na(ncc_preys) & ncc_preys != "" & ncc_preys != "na"]
  ncc_ints  <- c()
  
  for (i in seq_along(ncc_preys)) {
    pc <- ncc_preys[i]
    if (pc %in% names(df_model)) {
      in_name <- paste0("ncc_int_", i)
      df_model[[in_name]] <- df_model[[ncc_pred]] * df_model[[pc]]
      ncc_ints <- c(ncc_ints, in_name)
    }
  }
  
  # --- Build AK Interaction Products (ak_int_1, ak_int_2, ...) ---
  ak_preys <- c(ak_row$altprey1_data_col, ak_row$altprey2_data_col, ak_row$altprey3_data_col)
  ak_preys <- ak_preys[!is.na(ak_preys) & ak_preys != "" & ak_preys != "na"]
  ak_ints  <- c()
  
  for (j in seq_along(ak_preys)) {
    pc <- ak_preys[j]
    # Generate _2yrlead on the fly if missing
    if (!pc %in% names(df_model)) {
      bp <- str_remove(pc, "_2yrlead$")
      if (bp %in% names(df_model)) df_model[[pc]] <- dplyr::lead(df_model[[bp]], 2)
    }
    if (pc %in% names(df_model)) {
      in_name <- paste0("ak_int_", j)
      df_model[[in_name]] <- df_model[[ak_pred]] * df_model[[pc]]
      ak_ints <- c(ak_ints, in_name)
    }
  }
  
  # Helper to construct two-stage lavaan syntax string
  build_twostage_syntax <- function(active_ncc_ints, active_ak_ints) {
    rhs_ncc <- paste(c(ncc_pred, active_ncc_ints), collapse = " + ")
    rhs_ak  <- paste(c(target_ncc, ak_pred, active_ak_ints), collapse = " + ")
    
    paste0(
      target_ncc, " ~ ", rhs_ncc, "\n",
      target_ak,  " ~ ", rhs_ak
    )
  }
  
  # --- STEP 1: Fit Full Two-Stage Model ---
  curr_ncc_ints <- ncc_ints
  curr_ak_ints  <- ak_ints
  full_syntax   <- build_twostage_syntax(curr_ncc_ints, curr_ak_ints)
  
  fit_full <- tryCatch(
    sem(full_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit_full) || !lavInspect(fit_full, "converged")) return(NULL)
  
  fm_full <- fitMeasures(fit_full, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2_full_ncc <- inspect(fit_full, "r2")[[target_ncc]]
  r2_full_ak  <- inspect(fit_full, "r2")[[target_ak]]
  
  full_res <- tibble(
    ncc_predator = ncc_pred,
    ak_predator  = ak_pred,
    model_type   = "Full",
    active_ncc_ints = paste(curr_ncc_ints, collapse = ", "),
    active_ak_ints  = paste(curr_ak_ints, collapse = ", "),
    aic          = fm_full[["aic"]],
    bic          = fm_full[["bic"]],
    cfi          = fm_full[["cfi"]],
    pvalue       = fm_full[["pvalue"]],
    cpue_r2      = as.numeric(r2_full_ncc),
    sar_r2       = as.numeric(r2_full_ak)
  )
  
  # --- STEP 2: Stepwise Backward P-Value Pruning Across Both Stages ---
  fit_curr <- fit_full
  while ((length(curr_ncc_ints) + length(curr_ak_ints)) > 0) {
    pe <- parameterEstimates(fit_curr) %>%
      filter((lhs == target_ncc & rhs %in% curr_ncc_ints) | (lhs == target_ak & rhs %in% curr_ak_ints))
    
    if (nrow(pe) == 0) break
    
    # Identify least significant interaction across either stage
    max_p_row <- pe %>% filter(pvalue == max(pvalue, na.rm = TRUE)) %>% slice(1)
    
    if (is.na(max_p_row$pvalue) || max_p_row$pvalue <= 0.05) break
    
    worst_term <- max_p_row$rhs
    
    # Drop worst term from its respective active list
    if (worst_term %in% curr_ncc_ints) curr_ncc_ints <- setdiff(curr_ncc_ints, worst_term)
    if (worst_term %in% curr_ak_ints)  curr_ak_ints  <- setdiff(curr_ak_ints, worst_term)
    
    pruned_syntax <- build_twostage_syntax(curr_ncc_ints, curr_ak_ints)
    fit_curr <- tryCatch(
      sem(pruned_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(fit_curr) || !lavInspect(fit_curr, "converged")) break
  }
  
  fm_pruned <- fitMeasures(fit_curr, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2_pruned_ncc <- inspect(fit_curr, "r2")[[target_ncc]]
  r2_pruned_ak  <- inspect(fit_curr, "r2")[[target_ak]]
  
  pruned_res <- tibble(
    ncc_predator = ncc_pred,
    ak_predator  = ak_pred,
    model_type   = "Pruned",
    active_ncc_ints = if (length(curr_ncc_ints) > 0) paste(curr_ncc_ints, collapse = ", ") else "None",
    active_ak_ints  = if (length(curr_ak_ints) > 0)  paste(curr_ak_ints, collapse = ", ")  else "None",
    aic          = fm_pruned[["aic"]],
    bic          = fm_pruned[["bic"]],
    cfi          = fm_pruned[["cfi"]],
    pvalue       = fm_pruned[["pvalue"]],
    cpue_r2      = as.numeric(r2_pruned_ncc),
    sar_r2       = as.numeric(r2_pruned_ak)
  )
  
  bind_rows(full_res, pruned_res)
}

# -----------------------------------------------------------------------------
# 3. Run Factorial Sweep (All NCC Preds x All AK Preds)
# -----------------------------------------------------------------------------
cat("Executing Factorial Two-Stage SEM Sweep...\n")

two_stage_grid <- expand_grid(
  ncc_idx = seq_len(nrow(ncc_candidates)),
  ak_idx  = seq_len(nrow(ak_candidates))
)

all_two_stage_models <- two_stage_grid %>%
  mutate(
    res = purrr::map2(ncc_idx, ak_idx, ~ fit_two_stage_pair(ncc_candidates[.x, ], ak_candidates[.y, ], sem_complete_data))
  ) %>%
  pull(res) %>%
  bind_rows()

# Save complete results
write.csv(all_two_stage_models, file.path(out_dir, "two_stage_factorial_sem_sweep.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Summary Displays (Filtered for Goodness-of-Fit pvalue >= 0.05)
# -----------------------------------------------------------------------------
cat("\n====================================================================\n")
cat(" TOP VALID TWO-STAGE PRUNED MODELS (pvalue >= 0.05 & CFI >= 0.90)   \n")
cat("====================================================================\n")

top_valid_twostage <- all_two_stage_models %>%
  filter(model_type == "Pruned", pvalue >= 0.05, cfi >= 0.90) %>%
  arrange(aic) %>%
  mutate(delta_aic = aic - min(aic, na.rm = TRUE))

print(
  top_valid_twostage %>%
    select(ncc_predator, ak_predator, active_ncc_ints, active_ak_ints, aic, delta_aic, sar_r2, cfi, pvalue),
  n = Inf
)