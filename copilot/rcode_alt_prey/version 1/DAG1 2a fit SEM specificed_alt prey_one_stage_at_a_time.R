suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_9"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Lookup Table (Lower-case standard)
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

# -----------------------------------------------------------------------------
# 2. Automated SEM Fitting with Corrected Region & Prefix Routing
# -----------------------------------------------------------------------------
run_predator_sem_sweep <- function(row) {
  pred_col <- row$pred_data_col
  region   <- row$region
  
  # ROUTING FIX: NCC smolt predators -> x07; All AK predators + NCC Guild 10 -> x16_sar
  is_ncc_smolt_pred <- (region == "NCC" & grepl("^(x08_|x09_|x11_)", pred_col))
  target_stage      <- if (is_ncc_smolt_pred) target_ncc else target_ak
  
  prey_cols <- c(row$altprey1_data_col, row$altprey2_data_col, row$altprey3_data_col)
  prey_cols <- prey_cols[!is.na(prey_cols) & prey_cols != "" & prey_cols != "na"]
  
  df_model <- sem_complete_data
  int_cols <- c()
  
  for (i in seq_along(prey_cols)) {
    prey_c <- prey_cols[i]
    if (!prey_c %in% names(df_model)) {
      base_prey <- str_remove(prey_c, "_2yrlead$")
      if (base_prey %in% names(df_model)) {
        df_model[[prey_c]] <- dplyr::lead(df_model[[base_prey]], 2)
      }
    }
    
    if (prey_c %in% names(df_model)) {
      int_name <- paste0("int_", i)
      df_model[[int_name]] <- df_model[[pred_col]] * df_model[[prey_c]]
      int_cols <- c(int_cols, int_name)
    }
  }
  
  if (length(int_cols) == 0) return(NULL)
  
  # Construct lavaan equation based on stage
  build_sem_text <- function(active_ints) {
    if (target_stage == target_ncc) {
      # NCC Smolt Predators (NCC x08_, x09_, x11_) predict x07_dfa_cpue_intsprjunhw
      rhs_terms <- paste(c(pred_col, active_ints), collapse = " + ")
      paste0(
        target_ncc, " ~ ", rhs_terms, "\n",
        target_ak, " ~ ", target_ncc
      )
    } else {
      # AK Predators & Guild 10 Mammals predict x16_sar
      rhs_terms <- paste(c(target_ncc, pred_col, active_ints), collapse = " + ")
      paste0(
        target_ncc, " ~ x09_dfa_hakeage5plus\n",
        target_ak, " ~ ", rhs_terms
      )
    }
  }
  
  # --- STEP 1: Fit Full Model ---
  current_ints <- int_cols
  full_sem_text <- build_sem_text(current_ints)
  
  fit_full <- tryCatch(
    sem(full_sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit_full) || !lavInspect(fit_full, "converged")) return(NULL)
  
  fm_full <- fitMeasures(fit_full, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2_full <- inspect(fit_full, "r2")[[target_ak]]
  
  full_res <- tibble(
    predator = pred_col,
    region = region,
    stage = target_stage,
    model_type = "Full",
    active_interactions = paste(current_ints, collapse = ", "),
    aic = fm_full[["aic"]],
    bic = fm_full[["bic"]],
    cfi = fm_full[["cfi"]],
    pvalue = fm_full[["pvalue"]],
    sar_r2 = as.numeric(r2_full)
  )
  
  # --- STEP 2: Stepwise Backward P-Value Pruning ---
  fit_curr <- fit_full
  while (length(current_ints) > 0) {
    pe <- parameterEstimates(fit_curr) %>%
      filter(lhs == target_stage, op == "~", rhs %in% current_ints)
    
    if (nrow(pe) == 0) break
    
    max_p_row <- pe %>% filter(pvalue == max(pvalue, na.rm = TRUE)) %>% slice(1)
    
    if (is.na(max_p_row$pvalue) || max_p_row$pvalue <= 0.05) break
    
    worst_int <- max_p_row$rhs
    current_ints <- setdiff(current_ints, worst_int)
    
    pruned_sem_text <- build_sem_text(current_ints)
    fit_curr <- tryCatch(
      sem(pruned_sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(fit_curr) || !lavInspect(fit_curr, "converged")) break
  }
  
  fm_pruned <- fitMeasures(fit_curr, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2_pruned <- inspect(fit_curr, "r2")[[target_ak]]
  
  pruned_res <- tibble(
    predator = pred_col,
    region = region,
    stage = target_stage,
    model_type = "Pruned",
    active_interactions = if (length(current_ints) > 0) paste(current_ints, collapse = ", ") else "None (Main Effect Only)",
    aic = fm_pruned[["aic"]],
    bic = fm_pruned[["bic"]],
    cfi = fm_pruned[["cfi"]],
    pvalue = fm_pruned[["pvalue"]],
    sar_r2 = as.numeric(r2_pruned)
  )
  
  bind_rows(full_res, pruned_res)
}

# -----------------------------------------------------------------------------
# 3. Execute & Output
# -----------------------------------------------------------------------------
sweep_results <- pred_lookup %>%
  split(seq(nrow(.))) %>%
  map_dfr(run_predator_sem_sweep)

write.csv(sweep_results, file.path(out_dir, "sem_altprey_sweep_results_corrected_routing.csv"), row.names = FALSE)

cat("\n====================================================================\n")
cat(" CORRECTED SWEEP - FULL MODELS (pvalue >= 0.05)                      \n")
cat("====================================================================\n")

print(
  sweep_results %>%
    filter(model_type == "Full", pvalue >= 0.05) %>%
    arrange(aic) %>%
    select(predator, region, stage, active_interactions, aic, sar_r2, cfi, pvalue),
  n = Inf
)

cat("\n====================================================================\n")
cat(" CORRECTED SWEEP - PRUNED MODELS (pvalue >= 0.05)                    \n")
cat("====================================================================\n")

print(
  sweep_results %>%
    filter(model_type == "Pruned", pvalue >= 0.05) %>%
    arrange(aic) %>%
    select(predator, region, stage, active_interactions, aic, sar_r2, cfi, pvalue),
  n = Inf
)