# -----------------------------------------------------------------------------
# Extract Detailed Interaction Coefficients & Check Signs (Bug-Free)
# -----------------------------------------------------------------------------
extract_twostage_coefficients <- function(top_models_df, raw_data, ncc_lookup, ak_lookup) {
  coef_records <- list()
  
  for (i in seq_len(nrow(top_models_df))) {
    row <- top_models_df[i, ]
    ncc_pred <- row$ncc_predator
    ak_pred  <- row$ak_predator
    
    df_m <- raw_data
    
    # Re-build NCC interaction terms
    ncc_row <- ncc_lookup %>% filter(pred_data_col == ncc_pred) %>% slice(1)
    ncc_preys <- c(ncc_row$altprey1_data_col, ncc_row$altprey2_data_col, ncc_row$altprey3_data_col)
    ncc_preys <- ncc_preys[!is.na(ncc_preys) & ncc_preys != "" & ncc_preys != "na"]
    ncc_ints <- c()
    for (k in seq_along(ncc_preys)) {
      pc <- ncc_preys[k]
      if (pc %in% names(df_m)) {
        in_name <- paste0("ncc_int_", k)
        df_m[[in_name]] <- df_m[[ncc_pred]] * df_m[[pc]]
        ncc_ints <- c(ncc_ints, in_name)
      }
    }
    
    # Re-build AK interaction terms
    ak_row <- ak_lookup %>% filter(pred_data_col == ak_pred) %>% slice(1)
    ak_preys <- c(ak_row$altprey1_data_col, ak_row$altprey2_data_col, ak_row$altprey3_data_col)
    ak_preys <- ak_preys[!is.na(ak_preys) & ak_preys != "" & ak_preys != "na"]
    ak_ints <- c()
    for (k in seq_along(ak_preys)) {
      pc <- ak_preys[k]
      if (!pc %in% names(df_m)) {
        bp <- str_remove(pc, "_2yrlead$")
        if (bp %in% names(df_m)) df_m[[pc]] <- dplyr::lead(df_m[[bp]], 2)
      }
      if (pc %in% names(df_m)) {
        in_name <- paste0("ak_int_", k)
        df_m[[in_name]] <- df_m[[ak_pred]] * df_m[[pc]]
        ak_ints <- c(ak_ints, in_name)
      }
    }
    
    # Active interaction terms in pruned model
    active_ncc <- unlist(strsplit(row$active_ncc_ints, ",\\s*"))
    active_ncc <- active_ncc[active_ncc %in% ncc_ints]
    
    active_ak <- unlist(strsplit(row$active_ak_ints, ",\\s*"))
    active_ak <- active_ak[active_ak %in% ak_ints]
    
    all_active <- c(active_ncc, active_ak)
    
    # Skip fitting coefficient inspection if no active interaction terms exist
    if (length(all_active) == 0) {
      coef_records[[i]] <- tibble(
        model_rank   = i,
        ncc_predator = ncc_pred,
        ak_predator  = ak_pred,
        lhs          = "None",
        term         = "No Interactions",
        est          = NA_real_,
        se           = NA_real_,
        z            = NA_real_,
        pvalue       = NA_real_,
        aic          = row$aic,
        sign_check   = "No Interaction Terms"
      )
      next
    }
    
    # Build syntax
    rhs_ncc <- paste(c(ncc_pred, active_ncc), collapse = " + ")
    rhs_ak  <- paste(c("x07_dfa_cpue_intsprjunhw", ak_pred, active_ak), collapse = " + ")
    
    sem_syntax <- paste0(
      "x07_dfa_cpue_intsprjunhw ~ ", rhs_ncc, "\n",
      "x16_sar ~ ", rhs_ak
    )
    
    fit <- tryCatch(
      sem(sem_syntax, data = df_m, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      pe <- parameterEstimates(fit) %>%
        filter(op == "~", rhs %in% all_active) %>%
        select(lhs, term = rhs, est, se, z, pvalue) %>%
        mutate(
          model_rank   = i,
          ncc_predator = ncc_pred,
          ak_predator  = ak_pred,
          aic          = row$aic,
          sign_check   = ifelse(est < 0, "Buffering (Negative)", "Positive (Check Overfit)")
        )
      
      coef_records[[i]] <- pe
    }
  }
  
  bind_rows(coef_records)
}

# -----------------------------------------------------------------------------
# Execute & Print
# -----------------------------------------------------------------------------
top_valid_models <- all_two_stage_models %>%
  filter(model_type == "Pruned", pvalue >= 0.05, cfi >= 0.90) %>%
  arrange(aic) %>%
  slice_head(n = 20)

interaction_coef_table <- extract_twostage_coefficients(
  top_valid_models, 
  sem_complete_data, 
  ncc_candidates, 
  ak_candidates
) %>%
  as_tibble()

write.csv(interaction_coef_table, file.path(out_dir, "two_stage_interaction_coefficients.csv"), row.names = FALSE)

cat("\n====================================================================\n")
cat(" INTERACTION COEFFICIENT SIGNS FOR TOP VALID TWO-STAGE MODELS       \n")
cat("====================================================================\n")

print(
  interaction_coef_table %>%
    select(model_rank, ncc_predator, ak_predator, lhs, term, est, pvalue, sign_check),
  n = Inf
)
