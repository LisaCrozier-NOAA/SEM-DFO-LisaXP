# ==============================================================================
# STEP 2: EXTEND ALT PREY SERIES THROUGH 2023 (LOCKED 1998-2021 SCALING)
# ==============================================================================

# Load saved MARSS model fits from 1998-2021 baseline run
marss_fits <- readRDS(file.path(outputDir, "marss_fits.rds"))

extended_results <- list()
years_full       <- raw_data$Year  # 1998 to 2023

# ------------------------------------------------------------------------------
# BLOCK A: SINGLE-SERIES STRAGGLERS (Squid, Anchovy, Pollock)--------
# ------------------------------------------------------------------------------
straggler_targets <- tribble(
  ~LisaName,              ~raw_col,
  "X04_marketsquid_GAM",  "marketsquid_GAM_2025",
  "X05_anchovy_GAM",      "anchovy_GAM_2025",
  "X13_pollock_age1plus", "pollockBiomassAIage1plus_predAK_2026"
)

for (i in seq_len(nrow(straggler_targets))) {
  l_name  <- straggler_targets$LisaName[i]
  r_col   <- straggler_targets$raw_col[i]
  
  cat("Processing Straggler Extension for:", l_name, "-> Raw Col:", r_col, "\n")
  
  if (r_col %in% names(raw_data)) {
    raw_vec  <- raw_data[[r_col]]
    
    # 1. Standardize using locked 1998-2021 baseline mean/SD
    z_scaled <- scale_with_baseline_lock(raw_vec, years_full)
    
    # 2. Kalman Smooth via MARSSkf (fit = FALSE)
    fit_single <- MARSS(matrix(z_scaled, nrow = 1), fit = FALSE)
    fit_single$par <- fit_single$start
    kfList <- MARSSkf(fit_single)
    
    # 3. Store extended smoothed output
    extended_results[[l_name]] <- tibble(
      Year = years_full,
      LisaName = l_name,
      value = as.numeric(kfList$xtT)
    )
  } else {
    stop("CRITICAL ERROR: Straggler raw column missing from raw_data: ", r_col)
  }
}

# ------------------------------------------------------------------------------
# BLOCK B: MULTI-INDICATOR DFAs (Sardine DFA, ChinAbund DFA, Hake DFA, Euph DFA)--------
# ------------------------------------------------------------------------------
dfa_targets <- tribble(
  ~LisaName,                     ~DFAname,
  "X05_DFA_abundSardine",        "05.ForageFishNCC_DFA1",
  "X09_DFA_ChinAbundSnakeFall",  "09.PredFishNCC_DFA1",
  "X09_DFA_HakeAge5Plus",        "09.PredFishNCC_b_DFA1",
  "X12_DFA_biomassEuphShelfSum", "12.ZooPreyAK_DFA1"
)

for (i in seq_len(nrow(dfa_targets))) {
  l_name <- dfa_targets$LisaName[i]
  d_name <- dfa_targets$DFAname[i]
  
  cat("Processing DFA Projection for:", l_name, "(Fit key:", d_name, ")\n")
  
  fit_obj <- marss_fits[[d_name]]
  if (is.null(fit_obj)) {
    stop("CRITICAL ERROR: MARSS model fit not found in marss_fits.rds for key: ", d_name)
  }
  
  fitted_rows <- rownames(fit_obj$model$data)
  
  # 1. Build matrix for 26 years (1998–2023)
  mat_data <- matrix(NA_real_, nrow = length(fitted_rows), ncol = length(years_full))
  rownames(mat_data) <- fitted_rows
  colnames(mat_data) <- years_full
  
  for (r_idx in seq_along(fitted_rows)) {
    r_name <- fitted_rows[r_idx]
    
    target_col <- case_when(
      r_name %in% names(raw_data) ~ r_name,
      r_name == "abundSardine" && "sardine_GAM_2025" %in% names(raw_data) ~ "sardine_GAM_2025",
      r_name == "abundSardine" && "sardine_NCC" %in% names(raw_data) ~ "sardine_NCC",
      r_name == "abundHerring" && "herring_GAM_2025" %in% names(raw_data) ~ "herring_GAM_2025",
      r_name == "abundHerring" && "herring_NCC" %in% names(raw_data) ~ "herring_NCC",
      TRUE ~ NA_character_
    )
    
    if (!is.na(target_col) && target_col %in% names(raw_data)) {
      mat_data[r_idx, ] <- scale_with_baseline_lock(raw_data[[target_col]], years_full)
    }
  }
  
  # 2. Build model shell over 26-year matrix
  marss_model_26 <- MARSS(mat_data, model = list(m = 1), form = "dfa", fit = FALSE, silent = TRUE)
  
  # Assign locked parameter estimates directly to the model object
  marss_model_26$par <- fit_obj$par
  
  # 3. Run Kalman Smoother across extended 26-year window
  kf_ext <- MARSSkf(marss_model_26)
  ext_trend <- as.numeric(kf_ext$xtT[1, ])
  
  # Preserve factor sign convention
  Z_matrix <- fit_obj$par$Z
  maxIndex <- which.max(abs(Z_matrix[, 1]))
  if (Z_matrix[maxIndex, 1] < 0) {
    ext_trend <- -ext_trend
  }
  
  extended_results[[l_name]] <- tibble(
    Year = years_full,
    LisaName = l_name,
    value = ext_trend
  )
}

# Combine all 7 extended series (Block A + Block B)
all_extended_altprey <- bind_rows(extended_results)
cat("\nStep 2 Complete! Extended all 7 target series through 2023 cleanly.\n")
