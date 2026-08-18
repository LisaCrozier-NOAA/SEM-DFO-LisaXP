fit_shark_variant_sem_fxn <- function(data_base,col_shark, shark_temp_col, Q10=2, overlap_form, overlap_slope, shark_transform) {
  
  # 1. Shark Abundance Transformation
  if (!col_shark %in% names(data_base)) return(NULL)
  shark_vec <- data_base[[col_shark]]
  
  if (shark_transform == "z") {
    Shark_use <- safe_scale(shark_vec)
  } else if (shark_transform == "log1p_z") {
    minv <- suppressWarnings(min(shark_vec, na.rm = TRUE))
    if (!is.finite(minv)) return(NULL)
    Shark_use <- safe_scale(log1p(shark_vec - minv))
  } else if (shark_transform == "z_roll2") {
    z0   <- safe_scale(shark_vec)
    zlag <- dplyr::lag(z0, 1)
    Shark_use <- rowMeans(cbind(z0, zlag), na.rm = TRUE)
  } else {
    return(NULL)
  }
  
  # 2. Shark Temperature Driver & Metabolic Response
  if (!shark_temp_col %in% names(data_base)) return(NULL)
  T_shark_raw <- data_base[[shark_temp_col]]
  T_shark     <- safe_scale(T_shark_raw)
  
  Tref_shark  <- mean(T_shark_raw, na.rm = TRUE)
  M_t         <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  # 3. Spatial Overlap Form
  O_t <- if (overlap_form == "constant") {
    rep(1, nrow(data_base))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark) * 2
  } else {
    return(NULL)
  }
  
  # 4. Construct Integrated Index & Prep Model Frame
  fit_dat <- data_base %>%
    mutate(I_Shark = Shark_use * M_t * O_t) %>%
    filter(year >= 1998, year <= 2021) %>%
    rename(
      x07_dfa_cpue_int_spr_jun_hw = !!sym(col_salmon_juv),
      x16_sar                  = !!sym(col_salmon_adult)
    ) %>%
    filter(
      !is.na(x16_sar),
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  # 5. Lavaan SEM Model Syntax
  single_model_text <- '
    # Structural Model
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus 
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_Shark
  '
  
  # 6. Fit SEM with Error Handling
  fit_sem <- tryCatch({
    lavaan::sem(single_model_text, data = fit_dat, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit_sem) || !lavaan::lavInspect(fit_sem, "converged")) return(NULL)
  
  # 7. Extract Model Metrics and Parameters
  fit_measures <- lavaan::fitMeasures(fit_sem, c("aic", "bic", "cfi", "rmsea"))
  pe           <- lavaan::parameterEstimates(fit_sem)
  
  get_path <- function(lhs_val, rhs_val, field) {
    val <- pe %>% 
      filter(lhs == lhs_val, op == "~", rhs == rhs_val) %>% 
      pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  tibble(
    shark_col        = col_shark,
    shark_temp_col   = shark_temp_col,
    Q10              = Q10,
    overlap_form     = overlap_form,
    overlap_slope    = overlap_slope,
    shark_transform  = shark_transform,
    n                = nrow(fit_dat),
    
    aic              = fit_measures[["aic"]],
    bic              = fit_measures[["bic"]],
    cfi              = fit_measures[["cfi"]],
    rmsea            = fit_measures[["rmsea"]],
    
    # Structural Path Coefficients & P-Values
    b_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    b_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    b_shark_sar      = get_path("x16_sar", "I_Shark", "est"),
    p_shark_sar      = get_path("x16_sar", "I_Shark", "pvalue")
  )
}
