suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_9"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data, Exclude Year, and Shift Numeric Variables (+5.0)
# -----------------------------------------------------------------------------
sem_raw <- read.csv(file.path(out_dir, "sem_altprey_data_complete_1998_2021.csv"), row.names = NULL)
names(sem_raw) <- tolower(names(sem_raw))

# Shift all numeric columns EXCEPT year into positive domain
sem_complete_data <- sem_raw %>%
  mutate(across(where(is.numeric) & !matches("year"), ~ .x + 5.0))

pred_lookup <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv"), row.names = NULL) %>%
  clean_names() %>%
  mutate(across(everything(), tolower))

target_ncc <- "x07_dfa_cpue_intsprjunhw"
target_ak  <- "x16_sar"

ncc_candidates <- pred_lookup %>% filter(region == "ncc" & grepl("^(x08_|x09_|x11_)", pred_data_col))
ak_candidates  <- pred_lookup %>% filter(region == "ak" | grepl("^x10_", pred_data_col))

# Vectorized helper to map interaction codes to real prey column names
get_prey_col <- function(pred_row, int_index) {
  cols <- c(pred_row$altprey1_data_col, pred_row$altprey2_data_col, pred_row$altprey3_data_col)
  cols <- cols[!is.na(cols) & cols != "" & cols != "na"]
  cols[int_index]
}

# Safe R2 extraction helper
get_safe_r2 <- function(fit, target_var) {
  r2_vals <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
  if (!is.null(r2_vals) && target_var %in% names(r2_vals)) {
    return(as.numeric(r2_vals[[target_var]]))
  }
  return(NA_real_)
}

# -----------------------------------------------------------------------------
# 2. Combined Fitting & Parameter Storage Function
# -----------------------------------------------------------------------------
fit_and_extract_twostage_pair <- function(ncc_row, ak_row, df_positive, pair_id) {
  ncc_pred <- ncc_row$pred_data_col
  ak_pred  <- ak_row$pred_data_col
  
  df_model <- df_positive
  
  # Build NCC interaction columns
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
  
  # Build AK interaction columns
  ak_preys <- c(ak_row$altprey1_data_col, ak_row$altprey2_data_col, ak_row$altprey3_data_col)
  ak_preys <- ak_preys[!is.na(ak_preys) & ak_preys != "" & ak_preys != "na"]
  ak_ints  <- c()
  for (j in seq_along(ak_preys)) {
    pc <- ak_preys[j]
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
  
  build_syntax <- function(active_ncc, active_ak) {
    rhs_ncc <- if (length(active_ncc) > 0) paste(active_ncc, collapse = " + ") else "1"
    rhs_ak  <- paste(c(target_ncc, active_ak), collapse = " + ")
    paste0(target_ncc, " ~ ", rhs_ncc, "\n", target_ak, " ~ ", rhs_ak)
  }
  
  curr_ncc <- c(ncc_pred, ncc_ints)
  curr_ak  <- c(ak_pred, ak_ints)
  
  fit_curr <- tryCatch(
    sem(build_syntax(curr_ncc, curr_ak), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit_curr) || !lavInspect(fit_curr, "converged")) return(NULL)
  
  # Stepwise Backward Elimination across both stages
  while ((length(curr_ncc) + length(curr_ak)) > 0) {
    pe <- parameterEstimates(fit_curr) %>%
      filter((lhs == target_ncc & rhs %in% curr_ncc) | (lhs == target_ak & rhs %in% curr_ak))
    
    if (nrow(pe) == 0) break
    max_p_row <- pe %>% filter(pvalue == max(pvalue, na.rm = TRUE)) %>% slice(1)
    if (is.na(max_p_row$pvalue) || max_p_row$pvalue <= 0.05) break
    
    worst_term <- max_p_row$rhs
    if (worst_term %in% curr_ncc) curr_ncc <- setdiff(curr_ncc, worst_term)
    if (worst_term %in% curr_ak)  curr_ak  <- setdiff(curr_ak, worst_term)
    
    fit_curr <- tryCatch(
      sem(build_syntax(curr_ncc, curr_ak), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(fit_curr) || !lavInspect(fit_curr, "converged")) break
  }
  
  # Extract fit metrics safely
  fm <- fitMeasures(fit_curr, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2_cpue <- get_safe_r2(fit_curr, target_ncc)
  r2_sar  <- get_safe_r2(fit_curr, target_ak)
  
  summary_row <- tibble(
    pair_id          = pair_id,
    ncc_predator     = ncc_pred,
    ak_predator      = ak_pred,
    active_ncc_terms = if (length(curr_ncc) > 0) paste(curr_ncc, collapse = ", ") else "None",
    active_ak_terms  = if (length(curr_ak) > 0)  paste(curr_ak, collapse = ", ")  else "None",
    aic              = fm[["aic"]],
    bic              = fm[["bic"]],
    cfi              = fm[["cfi"]],
    pvalue           = fm[["pvalue"]],
    cpue_r2          = r2_cpue,
    sar_r2           = r2_sar
  )
  
  # Extract and annotate parameter estimates
  params_df <- parameterEstimates(fit_curr) %>%
    filter(op == "~") %>%
    select(lhs, term = rhs, est, se, z, pvalue) %>%
    mutate(
      pair_id      = pair_id,
      ncc_predator = ncc_pred,
      ak_predator  = ak_pred,
      aic          = fm[["aic"]],
      cfi          = fm[["cfi"]],
      global_p     = fm[["pvalue"]],
      resolved_species = case_when(
        term == ncc_pred ~ ncc_pred,
        term == ak_pred ~ ak_pred,
        term == target_ncc ~ "Early CPUE Link",
        grepl("^ncc_int_", term) ~ get_prey_col(ncc_row, as.numeric(str_remove(term, "ncc_int_"))),
        grepl("^ak_int_", term)  ~ get_prey_col(ak_row, as.numeric(str_remove(term, "ak_int_"))),
        TRUE ~ term
      ),
      term_type = case_when(
        term %in% c(ncc_pred, ak_pred) ~ "Main Predator Effect",
        term == target_ncc ~ "Stage Link (CPUE -> SAR)",
        TRUE ~ "Alternate Prey Interaction"
      ),
      effect_direction = case_when(
        term_type == "Main Predator Effect" & est < 0 ~ "Top-Down Predation (-)",
        term_type == "Main Predator Effect" & est > 0 ~ "Co-variance / Bottom-Up (+)",
        term_type == "Alternate Prey Interaction" & est > 0 ~ "Buffering / Prey-Switching (+)",
        term_type == "Alternate Prey Interaction" & est < 0 ~ "Apparent Competition (-)",
        TRUE ~ "Positive Link"
      )
    )
  
  list(summary = summary_row, params = params_df)
}

# -----------------------------------------------------------------------------
# 3. Run Factorial Sweep
# -----------------------------------------------------------------------------
cat("Running Complete Two-Stage Factorial Sweep with Direct Parameter Storage...\n")

two_stage_grid <- expand_grid(
  ncc_idx = seq_len(nrow(ncc_candidates)),
  ak_idx  = seq_len(nrow(ak_candidates))
) %>% mutate(pair_id = row_number())

all_results <- purrr::map(
  seq_len(nrow(two_stage_grid)),
  function(idx) {
    p_id <- two_stage_grid$pair_id[idx]
    n_row <- ncc_candidates[two_stage_grid$ncc_idx[idx], ]
    a_row <- ak_candidates[two_stage_grid$ak_idx[idx], ]
    fit_and_extract_twostage_pair(n_row, a_row, sem_complete_data, p_id)
  }
)

# Bind output tables
summary_results <- map_dfr(all_results, ~ .x$summary) %>% as_tibble()
all_params_df   <- map_dfr(all_results, ~ .x$params) %>% as_tibble()

# Save master files
write.csv(summary_results, file.path(out_dir, "two_stage_model_summaries_positive.csv"), row.names = FALSE)
write.csv(all_params_df, file.path(out_dir, "two_stage_all_parameters_positive.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Display Outputs
# -----------------------------------------------------------------------------
cat("\n====================================================================\n")
cat(" TOP VALID TWO-STAGE MODELS SUMMARY (pvalue >= 0.05 & CFI >= 0.90)   \n")
cat("====================================================================\n")

top_models_summary <- summary_results %>%
  filter(pvalue >= 0.05, cfi >= 0.90) %>%
  arrange(aic) %>%
  mutate(delta_aic = aic - min(aic, na.rm = TRUE))

print(top_models_summary, n = Inf)

cat("\n====================================================================\n")
cat(" PARAMETER ESTIMATES FOR TOP VALID MODELS                          \n")
cat("====================================================================\n")

top_params <- all_params_df %>%
  filter(pair_id %in% top_models_summary$pair_id) %>%
  left_join(top_models_summary %>% select(pair_id, delta_aic), by = "pair_id") %>%
  arrange(delta_aic, pair_id)

print(
  top_params %>%
    select(pair_id, delta_aic, stage = lhs, term, resolved_species, term_type, est, pvalue, effect_direction),
  n = Inf
)