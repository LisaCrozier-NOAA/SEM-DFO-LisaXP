# ==============================================================================
# Script: 03b_run_sst_altprey_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Unified Full 4-Way Competition Sweep (Raw vs Temp vs AltPrey vs Joint)
#          Evaluates NCC (CPUE target) and AK/Adult (SAR target) across 24yr and 19yr.
#          Fixed parameters: a = 0.2 (both Temp and AltPrey), w = 1.0.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & UTILITIES
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
if (!file.exists(raw_csv_path)) stop("Raw data missing at: ", raw_csv_path)

raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Hardcoded exclusion list from 10e
blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

# Scaling helper with custom bounds [a, b]
scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

target_cpue <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar  <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: GENERALIZED SWEEP ENGINE FUNCTION
# ------------------------------------------------------------------------------
run_single_predator_competition <- function(pred_col, region_type, df_raw, df_base, lookup_row, max_combo_size = 2) {
  
  if (is.na(pred_col) || pred_col == "" || !pred_col %in% names(df_base)) return(NULL)
  
  P_norm <- df_base[[pred_col]]
  
  # Thermal driver from lookup table
  spatial_col <- trimws(as.character(lookup_row$spatial_overlap_col[1]))
  has_temp <- (!is.na(spatial_col) && spatial_col != "" && spatial_col %in% names(df_raw))
  
  if (has_temp) {
    raw_env  <- df_raw[[spatial_col]]
    env_gate <- scale_custom_bounds(raw_env, a = 0.2, b = 1.0)
    P_temp   <- P_norm * env_gate
  } else {
    P_temp <- NULL
  }
  
  # Forage candidates (excluding spatial overlap and predator itself)
  all_cands <- unlist(strsplit(lookup_row$all_candidates_str[1], ", "))
  prey_candidates <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
  prey_candidates <- intersect(prey_candidates, names(df_base))
  
  # Set up SEM model string depending on region topology
  build_sem_string <- function(p_var) {
    if (region_type == "NCC") {
      paste0(target_cpue, " ~ ", p_var)
    } else {
      paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + ", p_var)
    }
  }
  
  eval_target <- if (region_type == "NCC") target_cpue else target_sar
  
  # --- Raw Baseline Model ---
  base_fit <- tryCatch(
    sem(build_sem_string(pred_col), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- tryCatch(parameterEstimates(base_fit, standardized = TRUE), error = function(e) NULL)
  if (is.null(base_pe) || !"rhs" %in% names(base_pe)) return(NULL)
  
  base_pe_sub <- base_pe %>% filter(op == "~" & lhs == eval_target & rhs == pred_col)
  if (nrow(base_pe_sub) == 0) return(NULL)
  
  base_est <- as.numeric(base_pe_sub$est[1])
  base_std <- as.numeric(base_pe_sub$std.all[1])
  base_p   <- as.numeric(base_pe_sub$pvalue[1])
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[eval_target]])
  
  best_abs <- list(variant_type = "Raw Baseline", prey_combo = "None",
                   aic = base_aic, delta_aic = 0, beta_peff = base_est,
                   beta_std = base_std, pvalue = base_p, r2 = base_r2)
  best_abs_aic <- base_aic
  
  best_neg <- if (base_est < 0) best_abs else NULL
  best_neg_aic <- if (base_est < 0) base_aic else Inf
  
  # --- VARIANT 1: SST Only (w = 1.0, a = 0.2) ---
  if (has_temp && !is.null(P_temp) && !is.na(sd(P_temp)) && sd(P_temp) >= 1e-6) {
    df_base$P_eff_test <- P_temp
    fit_temp <- tryCatch(sem(build_sem_string("P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
    if (!is.null(fit_temp) && lavInspect(fit_temp, "converged")) {
      pe <- tryCatch(parameterEstimates(fit_temp, standardized = TRUE), error = function(e) NULL)
      if (!is.null(pe) && "rhs" %in% names(pe)) {
        pe_sub <- pe %>% filter(op == "~" & lhs == eval_target & rhs == "P_eff_test")
        if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
          curr_aic <- as.numeric(fitMeasures(fit_temp, "aic"))
          curr_est <- as.numeric(pe_sub$est[1])
          res_obj  <- list(variant_type = "Temp Only", prey_combo = "None", aic = curr_aic, delta_aic = curr_aic - base_aic,
                           beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]), pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_temp, "r2")[[eval_target]]))
          if (curr_aic < best_abs_aic) { best_abs_aic <- curr_aic; best_abs <- res_obj }
          if (curr_est < 0 && curr_aic < best_neg_aic) { best_neg_aic <- curr_aic; best_neg <- res_obj }
        }
      }
    }
  }
  
  # --- VARIANTS 2 & 3: AltPrey Only AND Joint (Temp + AltPrey) ---
  # Attenuation gate uses bounds [0.2, 1.0] and fixed w = 1.0
  if (length(prey_candidates) > 0) {
    max_k <- min(length(prey_candidates), max_combo_size)
    
    for (k in 1:max_k) {
      combos <- combn(prey_candidates, k, simplify = FALSE)
      
      for (cb in combos) {
        if (length(cb) == 1) {
          raw_forage <- df_base[[cb[1]]]
        } else {
          raw_forage <- rowSums(df_base[, cb, drop = FALSE], na.rm = TRUE)
        }
        
        # Scale alternative prey to [0.2, 1.0]
        forage_scaled <- scale_custom_bounds(raw_forage, a = 0.2, b = 1.0)
        gate_alt <- 1 - 1.0 * forage_scaled # Fixed w = 1.0
        
        # A) AltPrey Only
        P_alt_only <- P_norm * gate_alt
        if (!is.na(sd(P_alt_only)) && sd(P_alt_only) >= 1e-6) {
          df_base$P_eff_test <- P_alt_only
          fit_alt <- tryCatch(sem(build_sem_string("P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
          if (!is.null(fit_alt) && lavInspect(fit_alt, "converged")) {
            pe <- tryCatch(parameterEstimates(fit_alt, standardized = TRUE), error = function(e) NULL)
            if (!is.null(pe) && "rhs" %in% names(pe)) {
              pe_sub <- pe %>% filter(op == "~" & lhs == eval_target & rhs == "P_eff_test")
              if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
                curr_aic <- as.numeric(fitMeasures(fit_alt, "aic"))
                curr_est <- as.numeric(pe_sub$est[1])
                res_obj  <- list(variant_type = paste0("AltPrey Only (", k, " Prey)"), prey_combo = paste(cb, collapse = " + "),
                                 aic = curr_aic, delta_aic = curr_aic - base_aic, beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]),
                                 pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_alt, "r2")[[eval_target]]))
                if (curr_aic < best_abs_aic) { best_abs_aic <- curr_aic; best_abs <- res_obj }
                if (curr_est < 0 && curr_aic < best_neg_aic) { best_neg_aic <- curr_aic; best_neg <- res_obj }
              }
            }
          }
        }
        
        # B) Joint Temp + AltPrey
        if (has_temp && !is.null(P_temp)) {
          P_joint <- P_temp * gate_alt
          if (!is.na(sd(P_joint)) && sd(P_joint) >= 1e-6) {
            df_base$P_eff_test <- P_joint
            fit_jnt <- tryCatch(sem(build_sem_string("P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
            if (!is.null(fit_jnt) && lavInspect(fit_jnt, "converged")) {
              pe <- tryCatch(parameterEstimates(fit_jnt, standardized = TRUE), error = function(e) NULL)
              if (!is.null(pe) && "rhs" %in% names(pe)) {
                pe_sub <- pe %>% filter(op == "~" & lhs == eval_target & rhs == "P_eff_test")
                if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
                  curr_aic <- as.numeric(fitMeasures(fit_jnt, "aic"))
                  curr_est <- as.numeric(pe_sub$est[1])
                  res_obj  <- list(variant_type = paste0("Joint Temp + ", k, " AltPrey"), prey_combo = paste(cb, collapse = " + "),
                                   aic = curr_aic, delta_aic = curr_aic - base_aic, beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]),
                                   pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_jnt, "r2")[[eval_target]]))
                  if (curr_aic < best_abs_aic) { best_abs_aic <- curr_aic; best_abs <- res_obj }
                  if (curr_est < 0 && curr_aic < best_neg_aic) { best_neg_aic <- curr_aic; best_neg <- res_obj }
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(tibble(
    predator            = pred_col,
    region_type         = region_type,
    spatial_col         = ifelse(has_temp, spatial_col, "None"),
    base_beta           = base_est,
    base_beta_std       = base_std,
    base_pvalue         = base_p,
    base_aic            = base_aic,
    base_r2             = base_r2,
    
    # Absolute Best Model
    abs_variant_type    = best_abs$variant_type,
    abs_prey            = best_abs$prey_combo,
    abs_beta_peff       = best_abs$beta_peff,
    abs_beta_std        = best_abs$beta_std,
    abs_pvalue          = best_abs$pvalue,
    abs_final_aic       = best_abs$aic,
    abs_delta_aic       = best_abs$delta_aic,
    abs_final_r2        = best_abs$r2,
    
    # Diagnostics
    Q1_aic_improved     = (best_abs$delta_aic < 0),
    Q2_beta_is_negative = (best_abs$beta_peff < 0),
    Q3_sign_flipped     = (base_est > 0 && best_abs$beta_peff < 0),
    Q4_impact_increased = (best_abs$beta_peff < base_est),
    
    # Best Negative Beta Model
    neg_variant_type    = if (!is.null(best_neg)) best_neg$variant_type else "None",
    neg_prey            = if (!is.null(best_neg)) best_neg$prey_combo else "None",
    neg_beta_peff       = if (!is.null(best_neg)) best_neg$beta_peff else NA_real_,
    neg_pvalue          = if (!is.null(best_neg)) best_neg$pvalue else NA_real_,
    neg_final_aic       = if (!is.null(best_neg)) best_neg$aic else NA_real_,
    neg_delta_aic       = if (!is.null(best_neg)) best_neg$delta_aic else NA_real_,
    neg_final_r2        = if (!is.null(best_neg)) best_neg$r2 else NA_real_
  ))
}

# ------------------------------------------------------------------------------
# STEP 2: MASTER ORCHESTRATION PIPELINE
# ------------------------------------------------------------------------------
# Execution sequence ordered by priority
execution_queue <- list(
  list(years = "24yr", region = "NCC"),
  list(years = "24yr", region = "AK"),
  list(years = "19yr", region = "NCC"),
  list(years = "19yr", region = "AK")
)

for (cfg in execution_queue) {
  cur_years  <- cfg$years
  cur_region <- cfg$region
  
  cat("\n==============================================================================\n")
  cat(sprintf("STARTING SWEEP: %s Dataset | %s Region | %s\n", cur_years, cur_region, format(Sys.time(), "%H:%M:%S")))
  cat("==============================================================================\n")
  
  # 1. Filter Data by Year Threshold
  req_count <- if (cur_years == "24yr") 24 else 19
  valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) >= req_count]
  
  df_filtered <- raw_df_master %>% 
    select(all_of(valid_cols)) %>% 
    filter(complete.cases(.))
  
  # Standard 0-1 min-max scaling for base regression variables
  df_minmax <- df_filtered %>%
    mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))
  
  # 2. Load Lookup Table
  lookup_path <- file.path(data_out_dir, paste0("guild_predator_lookup_table_", cur_years, ".csv"))
  if (!file.exists(lookup_path)) {
    warning("Lookup table missing at: ", lookup_path, ". Skipping run.")
    next
  }
  guild_pred_lookup <- read_csv(lookup_path, show_col_types = FALSE)
  
  # 3. Filter Candidates according to Specific Domain Criteria
  if (cur_region == "NCC") {
    pred_lookup_filtered <- guild_pred_lookup %>% 
      filter(region == "NCC", stage_suffix == "smoltyr", !grepl("DFO", pred_data_col))
  } else {
    pred_lookup_filtered <- guild_pred_lookup %>% 
      filter(
        region == "AK" | 
          (region == "NCC" & stage_suffix == "adultyr") | 
          (region == "NCC" & stage_suffix == "smoltyr" & grepl("DFO", pred_data_col))
      )
  }
  
  # Apply blacklist filter
  target_preds <- sort(unique(pred_lookup_filtered$pred_data_col))
  target_preds <- target_preds[!sapply(target_preds, function(x) any(sapply(blacklist, grepl, x)))]
  target_preds <- intersect(target_preds, names(df_minmax))
  
  cat(sprintf("Identified %d candidate predators to evaluate.\n", length(target_preds)))
  if (length(target_preds) == 0) next
  
  out_csv <- file.path(data_out_dir, paste0("Goal3_Competition_Sweep_", cur_region, "_", cur_years, ".csv"))
  sweep_results <- list()
  start_time    <- Sys.time()
  
  # 4. Iterate Over Predators with 10-Item Checkpointing
  for (i in seq_along(target_preds)) {
    pred_name <- target_preds[i]
    pred_row  <- pred_lookup_filtered %>% filter(pred_data_col == pred_name) %>% slice(1)
    
    iter_start <- Sys.time()
    res <- run_single_predator_competition(
      pred_col        = pred_name,
      region_type     = cur_region,
      df_raw          = df_filtered,
      df_base         = df_minmax,
      lookup_row      = pred_row,
      max_combo_size  = 2
    )
    
    if (!is.null(res)) {
      sweep_results[[length(sweep_results) + 1]] <- res
      
      # Live Terminal Feedback
      elapsed_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
      rate_per_pred <- elapsed_total / i
      rem_preds     <- length(target_preds) - i
      eta_sec       <- round(rate_per_pred * rem_preds, 0)
      
      cat(sprintf("[%d/%d] %s | Best: %s (dAIC = %.2f) | Beta: %.3f | ETA: %ds\n",
                  i, length(target_preds), pred_name, res$abs_variant_type, res$abs_delta_aic, res$abs_beta_peff, eta_sec))
    }
    
    # Checkpoint to CSV every 10 predators
    if (i %% 10 == 0 && length(sweep_results) > 0) {
      checkpoint_df <- bind_rows(sweep_results) %>% arrange(abs_final_aic)
      write_csv(checkpoint_df, out_csv)
      cat(sprintf("--> [CHECKPOINT] Saved %d records to %s\n", nrow(checkpoint_df), basename(out_csv)))
    }
  }
  
  # Final Save for Completed Section
  if (length(sweep_results) > 0) {
    final_section_df <- bind_rows(sweep_results) %>% arrange(abs_final_aic)
    write_csv(final_section_df, out_csv)
    cat(sprintf("--> [COMPLETE] Successfully saved %d models to: %s\n\n", nrow(final_section_df), out_csv))
  }
}

cat("Master Competition Sweep Pipeline Finished Across All Datasets & Regions!\n")