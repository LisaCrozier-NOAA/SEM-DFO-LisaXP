


# ==============================================================================
# Script: 10e_run_ak_adult_full_competition_sweep_10212016.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Full 4-way competition sweep on SAR (Raw vs. SST Only vs. AltPrey Only vs. Joint):
#          CPUE ~ Hake_smoltyr
#          SAR  ~ CPUE + Peff
#          Fixed lower bound a = 0.2. Tests 1 & 2 prey pools with and without SST.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION & DATA LOADING
# ------------------------------------------------------------------------------
years_included <- "19yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]


if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))
guild_pred_lookup <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)

target_cpue <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar  <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

all_cols <- names(filtered_df)

blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

candidate_preds <- sort(all_cols[grepl("^(X10_|X15_)", all_cols) | all_cols == "X11_ssl_seak_pup_pred_smoltyr"])
candidate_preds <- candidate_preds[!sapply(candidate_preds, function(x) any(sapply(blacklist, grepl, x)))]

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))


# ------------------------------------------------------------------------------
# STEP 1: FULL 4-WAY COMPETITION SWEEP FUNCTION
# ------------------------------------------------------------------------------
test_ak_adult_full_competition <- function(pred_col, df_raw, df_base, max_combo_size = 2) {
  
  P_norm <- df_base[[pred_col]]
  
  row_data <- guild_pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  stage_suffix <- ifelse(grepl("adultyr$", pred_col), "adultyr", "smoltyr")
  region <- ifelse(grepl("^X15_", pred_col) | pred_col == "X11_ssl_seak_pup_pred_smoltyr", "AK", "NCC_adult")
  
  spatial_col <- case_when(
    stage_suffix == "smoltyr" ~ "X21_sst_egoa_junjulaug_smoltyr",
    stage_suffix == "adultyr" ~ "X21_sst_egoa_junjulaug_adultyr"
    # region == "AK" & stage_suffix == "smoltyr" ~ "X21_sst_egoa_junjulaug_smoltyr",
    # region == "AK" & stage_suffix == "adultyr" ~ "X21_sst_egoa_junjulaug_adultyr",
    # region == "NCC_adult" ~ paste0("X01_habCompInd_", stage_suffix)
  )
  if (!spatial_col %in% names(df_raw)) return(NULL)
  
  raw_env  <- df_raw[[spatial_col]]
  env_gate <- scale_custom_bounds(raw_env, a = 0.2, b = 1.0)
  P_temp   <- P_norm * env_gate
  
  all_cands <- unlist(strsplit(row_data$all_candidates_str, ", "))
  prey_candidates <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
  prey_candidates <- intersect(prey_candidates, names(df_base))
  
  # --- Raw Baseline Model (CPUE ~ Hake; SAR ~ CPUE + Raw Predator) ---
  base_model_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + ", pred_col)
  base_fit <- tryCatch(sem(base_model_str, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- tryCatch(parameterEstimates(base_fit, standardized = TRUE), error = function(e) NULL)
  if (is.null(base_pe) || !"rhs" %in% names(base_pe)) return(NULL)
  
  base_pe_sub <- base_pe %>% filter(op == "~" & lhs == target_sar & rhs == pred_col)
  if (nrow(base_pe_sub) == 0) return(NULL)
  
  base_est <- as.numeric(base_pe_sub$est[1])
  base_std <- as.numeric(base_pe_sub$std.all[1])
  base_p   <- as.numeric(base_pe_sub$pvalue[1])
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_sar]])
  
  best_abs     <- list(variant_type = "Raw Baseline", prey_combo = "None", weight_w = NA_real_,
                       aic = base_aic, delta_aic = 0, beta_peff = base_est, beta_std = base_std, pvalue = base_p, r2 = base_r2)
  best_abs_aic <- base_aic
  
  best_neg     <- if (base_est < 0) best_abs else NULL
  best_neg_aic <- if (base_est < 0) base_aic else Inf
  
  w_grid <- seq(0.1, 1.0, by = 0.1)
  
  # --- VARIANT 1: SST Only ---
  df_base$P_eff_test <- P_temp
  fit_sst <- tryCatch(sem(paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  if (!is.null(fit_sst) && lavInspect(fit_sst, "converged")) {
    pe <- tryCatch(parameterEstimates(fit_sst, standardized = TRUE), error = function(e) NULL)
    if (!is.null(pe) && "rhs" %in% names(pe)) {
      pe_sub <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_eff_test")
      if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
        curr_aic <- as.numeric(fitMeasures(fit_sst, "aic"))
        curr_est <- as.numeric(pe_sub$est[1])
        res_obj  <- list(variant_type = "SST Only", prey_combo = "None", weight_w = NA_real_, aic = curr_aic, delta_aic = curr_aic - base_aic,
                         beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]), pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_sst, "r2")[[target_sar]]))
        if (curr_aic < best_abs_aic) { best_abs_aic <- curr_aic; best_abs <- res_obj }
        if (curr_est < 0 && curr_aic < best_neg_aic) { best_neg_aic <- curr_aic; best_neg <- res_obj }
      }
    }
  }
  
  # --- VARIANTS 2 & 3: AltPrey Only (No SST) AND Joint (SST + AltPrey) ---
  if (length(prey_candidates) > 0) {
    max_k <- min(length(prey_candidates), max_combo_size)
    
    for (k in 1:max_k) {
      combos <- combn(prey_candidates, k, simplify = FALSE)
      
      for (cb in combos) {
        if (length(cb) == 1) {
          forage_pool <- df_base[[cb[1]]]
        } else {
          raw_sum <- rowSums(df_base[, cb, drop = FALSE], na.rm = TRUE)
          forage_pool <- scale_custom_bounds(raw_sum, 0.0, 1.0)
        }
        
        for (w in w_grid) {
          # --- A) AltPrey Only (No SST) ---
          P_alt_only <- P_norm * (1 - w * forage_pool)
          if (!is.na(sd(P_alt_only)) && sd(P_alt_only) >= 1e-6) {
            df_base$P_eff_test <- P_alt_only
            fit_alt <- tryCatch(sem(paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
            if (!is.null(fit_alt) && lavInspect(fit_alt, "converged")) {
              pe <- tryCatch(parameterEstimates(fit_alt, standardized = TRUE), error = function(e) NULL)
              if (!is.null(pe) && "rhs" %in% names(pe)) {
                pe_sub <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_eff_test")
                if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
                  curr_aic <- as.numeric(fitMeasures(fit_alt, "aic"))
                  curr_est <- as.numeric(pe_sub$est[1])
                  res_obj  <- list(variant_type = paste0("AltPrey Only (", k, " Prey)"), prey_combo = paste(cb, collapse = " + "), weight_w = w,
                                   aic = curr_aic, delta_aic = curr_aic - base_aic, beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]), pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_alt, "r2")[[target_sar]]))
                  if (curr_aic < best_abs_aic) { best_abs_aic <- curr_aic; best_abs <- res_obj }
                  if (curr_est < 0 && curr_aic < best_neg_aic) { best_neg_aic <- curr_aic; best_neg <- res_obj }
                }
              }
            }
          }
          
          # --- B) Joint SST + AltPrey ---
          P_joint <- P_temp * (1 - w * forage_pool)
          if (!is.na(sd(P_joint)) && sd(P_joint) >= 1e-6) {
            df_base$P_eff_test <- P_joint
            fit_jnt <- tryCatch(sem(paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_test"), data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
            if (!is.null(fit_jnt) && lavInspect(fit_jnt, "converged")) {
              pe <- tryCatch(parameterEstimates(fit_jnt, standardized = TRUE), error = function(e) NULL)
              if (!is.null(pe) && "rhs" %in% names(pe)) {
                pe_sub <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_eff_test")
                if (nrow(pe_sub) > 0 && !is.na(pe_sub$pvalue[1])) {
                  curr_aic <- as.numeric(fitMeasures(fit_jnt, "aic"))
                  curr_est <- as.numeric(pe_sub$est[1])
                  res_obj  <- list(variant_type = paste0("Joint SST + ", k, " AltPrey"), prey_combo = paste(cb, collapse = " + "), weight_w = w,
                                   aic = curr_aic, delta_aic = curr_aic - base_aic, beta_peff = curr_est, beta_std = as.numeric(pe_sub$std.all[1]), pvalue = as.numeric(pe_sub$pvalue[1]), r2 = as.numeric(inspect(fit_jnt, "r2")[[target_sar]]))
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
  
  # Diagnostic Evaluators
  abs_daic <- best_abs$delta_aic
  abs_beta <- best_abs$beta_peff
  
  return(tibble(
    predator            = pred_col,
    spatial_col         = spatial_col,
    base_beta           = base_est,
    base_beta_std       = base_std,
    base_pvalue         = base_p,
    base_aic            = base_aic,
    base_r2             = base_r2,
    
    # Absolute Best AIC Model Results
    abs_variant_type    = best_abs$variant_type,
    abs_prey            = best_abs$prey_combo,
    abs_weight_w        = best_abs$weight_w,
    abs_beta_peff       = best_abs$beta_peff,
    abs_beta_std        = best_abs$beta_std,
    abs_pvalue          = best_abs$pvalue,
    abs_final_aic       = best_abs$aic,
    abs_delta_aic       = best_abs$delta_aic,
    abs_final_r2        = best_abs$r2,
    
    # Diagnostic Evaluators
    Q1_aic_improved     = (best_abs$delta_aic < 0),
    Q2_beta_is_negative = (best_abs$beta_peff < 0),
    Q3_sign_flipped_pos_to_neg = (base_est > 0 && best_abs$beta_peff < 0),
    Q4_impact_increased = (best_abs$beta_peff < base_est),
    
    # Best Negative Beta Model Results
    neg_variant_type    = if (!is.null(best_neg)) best_neg$variant_type else "None",
    neg_prey            = if (!is.null(best_neg)) best_neg$prey_combo else "None",
    neg_weight_w        = if (!is.null(best_neg)) best_neg$weight_w else NA_real_,
    neg_beta_peff       = if (!is.null(best_neg)) best_neg$beta_peff else NA_real_,
    neg_pvalue          = if (!is.null(best_neg)) best_neg$pvalue else NA_real_,
    neg_final_aic       = if (!is.null(best_neg)) best_neg$aic else NA_real_,
    neg_delta_aic       = if (!is.null(best_neg)) best_neg$delta_aic else NA_real_,
    neg_final_r2        = if (!is.null(best_neg)) best_neg$r2 else NA_real_
  ))
}

# ------------------------------------------------------------------------------
# STEP 2: EXECUTE FULL COMPETITION SWEEP
# ------------------------------------------------------------------------------
cat(sprintf("Executing Full 4-Way Competition Sweep on %s dataset...\n", years_included))

sweep_results <- list()

for (i in seq_along(candidate_preds)) {
  pred_col <- candidate_preds[i]
  cat(sprintf("[%d/%d] Testing predator: %s | Time: %s\n", i, length(candidate_preds), pred_col, format(Sys.time(), "%H:%M:%S")))
  
  if (pred_col %in% names(minmax_base)) {
    res_df <- test_ak_adult_full_competition(pred_col = pred_col, df_raw = filtered_df, df_base = minmax_base, max_combo_size = 2)
    if (!is.null(res_df)) {
      sweep_results[[length(sweep_results) + 1]] <- res_df
      print(res_df %>% select(predator, base_beta, abs_variant_type, abs_prey, abs_beta_peff, abs_delta_aic, Q1_aic_improved))
    }
  }
}

full_comp_table <- bind_rows(sweep_results) %>% arrange(abs_final_aic)

# Summary
cat("\n==============================================================================\n")
cat(sprintf("   AK & ADULT FULL 4-WAY COMPETITION RESULTS (SAR ~ Peff) - %s DATASET         \n", years_included))
cat("==============================================================================\n\n")

print(
  full_comp_table %>% 
    select(predator,base_aic,abs_final_aic,base_r2,abs_final_r2,  base_beta, abs_variant_type, abs_prey, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved),
  n = 50
)

out_csv <- file.path(data_out_dir, paste0("Goal3_AK_Adult_Full_4Way_Competition_Sweep_", years_included, ".csv"))
write_csv(full_comp_table, out_csv)
cat(sprintf("\nSweep results successfully saved to: %s\n", out_csv))

print(
  full_comp_table %>%
    filter(Q1_aic_improved=="TRUE", )%>% 
    select(predator,base_aic,abs_final_aic,base_r2,abs_final_r2,  base_beta, abs_variant_type, abs_prey, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved),
  n = 50
)

#just wcvi hs and pollock still positive, but you wouldn't choose those models anyway:
# 1 X10_harbour_seal_DFO_Coastwide_smoltyr           -27.3         -33.1   0.506        0.607    -0.314 Joint SST + 1 AltPrey X05_DFA_abundSardine_smoltyr                          0.5        -0.998        -5.80  TRUE           
# 2 X15_DFA_sleeperSharkBSAI_predAK_smoltyr          -26.4         -29.8   0.481        0.556    -0.246 Joint SST + 1 AltPrey capelin_WGoA_smoltyr                                  1          -0.474        -3.40  TRUE           
# 3 X15_sablefishAK_biomass4plus_smoltyr             -25.6         -28.8   0.466        0.540     0.420 Joint SST + 1 AltPrey sitkaHerring_EGoA_smoltyr                             1          -1.62         -3.19  TRUE           
# 4 X15_sablefishAK_biomass3plus_smoltyr             -27.6         -28.8   0.519        0.534     0.504 Joint SST + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + sitkaHerr…          1          -1.93         -1.22  TRUE           
# 5 X10_harbour_seal_DFO_WCVI_smoltyr                -24.5         -28.4   0.432        0.532     0.240 AltPrey Only (2 Prey) biomassEuphShelfSum_smoltyr + X05_DFA_abund…          1           0.401        -3.94  TRUE           
# 6 X15_pollockGOA_SSB_lag1_smoltyr                  -25.5         -27.8   0.466        0.544     0.263 AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.8         0.455        -2.27  TRUE           
# 7 X15_DFA_halibutBiomassAge8plus_predAK_adultyr    -25.1         -27.1   0.457        0.502    -0.287 Joint SST + 1 AltPrey capelin_WGoA_adultyr                                  0.6        -0.443        -2.03  TRUE           
# 8 X11_ssl_seak_pup_pred_smoltyr                    -23.1         -27.1   0.393        0.491     0.122 Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1          -0.467        -4.03  TRUE           
# 9 X15_pollockGOA_SSB_smoltyr                       -24.1         -24.7   0.426        0.457     0.196 AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.7         0.290        -0.641 TRUE           
# 10 X15_salmonSharkGoA_predAK_smoltyr                -23.3         -24.3   0.410        0.432    -0.220 Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1          -0.403        -1.01  TRUE           


#RESULTS-------
# AK & ADULT FULL 4-WAY COMPETITION RESULTS (SAR ~ Peff) - 24yr DATASET         
# predator                                      base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_variant_type      abs_prey                                     abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved
# 1 X10_harbour_seal_DFO_Coastwide_smoltyr           -27.3         -33.1   0.506        0.607   -0.314  Joint SST + 1 AltPrey X05_DFA_abundSardine_smoltyr                          0.5       -0.998         -5.80  TRUE           
# 2 X15_DFA_sleeperSharkBSAI_predAK_smoltyr          -26.4         -29.8   0.481        0.556   -0.246  Joint SST + 1 AltPrey capelin_WGoA_smoltyr                                  1         -0.474         -3.40  TRUE           
# 3 X10_Harbour_s_WS_adultyr                         -29.0         -29.0   0.516        0.516   -0.530  Raw Baseline          None                                                 NA         -0.530          0     FALSE          
# 4 X15_sablefishAK_biomass4plus_smoltyr             -25.6         -28.8   0.466        0.540    0.420  Joint SST + 1 AltPrey sitkaHerring_EGoA_smoltyr                             1         -1.62          -3.19  TRUE           
# 5 X15_sablefishAK_biomass3plus_smoltyr             -27.6         -28.8   0.519        0.534    0.504  Joint SST + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + sitkaHerr…          1         -1.93          -1.22  TRUE           
# 6 X10_harbour_seal_DFO_WCVI_smoltyr                -24.5         -28.4   0.432        0.532    0.240  AltPrey Only (2 Prey) biomassEuphShelfSum_smoltyr + X05_DFA_abund…          1          0.401         -3.94  TRUE           
# 7 X15_pollockGOA_SSB_lag1_smoltyr                  -25.5         -27.8   0.466        0.544    0.263  AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.8        0.455         -2.27  TRUE           
# 8 X15_DFA_halibutBiomassAge8plus_predAK_adultyr    -25.1         -27.1   0.457        0.502   -0.287  Joint SST + 1 AltPrey capelin_WGoA_adultyr                                  0.6       -0.443         -2.03  TRUE           
# 9 X11_ssl_seak_pup_pred_smoltyr                    -23.1         -27.1   0.393        0.491    0.122  Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1         -0.467         -4.03  TRUE           
# 10 X10_Killer_whales_NR_BC_adultyr                  -25.4         -25.4   0.452        0.452    0.278  Raw Baseline          None                                                 NA          0.278          0     FALSE          
# 11 X15_pollockGOA_SSB_smoltyr                       -24.1         -24.7   0.426        0.457    0.196  AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.7        0.290         -0.641 TRUE           
# 12 X10_Northern_f_s_WS_adultyr                      -24.4         -24.4   0.435        0.435    0.325  Raw Baseline          None                                                 NA          0.325          0     FALSE          
# 13 X10_DFA_Steller_s_l_WS_adultyr                   -24.3         -24.3   0.430        0.430    0.192  Raw Baseline          None                                                 NA          0.192          0     FALSE          
# 14 X15_salmonSharkGoA_predAK_smoltyr                -23.3         -24.3   0.410        0.432   -0.220  Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1         -0.403         -1.01  TRUE           
# 15 X10_DFA_AllSeaLionsEMB_adultyr                   -22.6         -22.6   0.381        0.381    0.0972 Raw Baseline          None                                                 NA          0.0972         0     FALSE          


#19yr models w/ improved aic-----------
# predator                                      base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_variant_type      abs_prey                                     abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved
# 1 X10_harbour_seal_DFO_WCVI_smoltyr                -28.0         -40.3   0.763        0.883  0.319    AltPrey Only (2 Prey) X04_marketsquid_GAM_smoltyr + X05_eulachon_…          1          1.56         -12.4   TRUE           
# 2 X15_DFA_sleeperSharkBSAI_predAK_adultyr          -32.4         -33.8   0.816        0.832 -0.439    AltPrey Only (2 Prey) sitkaHerring_EGoA_adultyr + capelin_WGoA_ad…          1         -0.530         -1.38  TRUE           
# 3 X15_spinyDogfishBSAI_predAK_smoltyr              -23.4         -33.3   0.659        0.819  0.162    Joint SST + 2 AltPrey capelin_WGoA_smoltyr + X12_egoa_krill_smolt…          1         -1.07          -9.87  TRUE           
# 4 X15_DFA_sleeperSharkBSAI_predAK_smoltyr          -31.6         -33.0   0.808        0.822 -0.314    Joint SST + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + capelin_W…          1         -0.395         -1.33  TRUE           
# 5 X10_Harbour_s_WS_adultyr                         -30.7         -32.5   0.792        0.813 -0.462    Joint SST + 1 AltPrey X05_eulachon_during_chinook_adultyr                   0.8       -0.589         -1.80  TRUE           
# 6 X11_ssl_seak_pup_pred_smoltyr                    -30.7         -32.0   0.798        0.819  0.358    AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.4        0.440         -1.24  TRUE           
# 7 X15_DFA_halibutBiomassAge8plus_predAK_adultyr    -25.0         -31.5   0.709        0.822 -0.222    AltPrey Only (1 Prey) capelin_WGoA_adultyr                                  0.9       -0.354         -6.52  TRUE           
# 8 X15_pollockGOA_SSB_smoltyr                       -30.2         -31.4   0.795        0.819  0.344    AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.5        0.442         -1.19  TRUE           
# 9 X15_sablefishAK_biomass4plus_smoltyr             -24.9         -30.5   0.707        0.786 -0.193    Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1         -0.562         -5.56  TRUE           
# 10 X10_harbour_seal_DFO_Coastwide_smoltyr           -27.4         -30.0   0.754        0.787 -0.289    SST Only              None                                                 NA         -0.386         -2.59  TRUE           
# 11 X15_sablefishAK_biomass3plus_smoltyr             -23.2         -29.1   0.664        0.765 -0.118    Joint SST + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + sitkaHerr…          1         -0.858         -5.95  TRUE           
# 12 X15_pollockGOA_SSB_lag1_smoltyr                  -27.1         -28.7   0.745        0.781  0.256    AltPrey Only (1 Prey) X13_pollockGOA_age1and2_smoltyr                       0.6        0.387         -1.60  TRUE           
# 13 X15_salmonSharkBSAI_predAK_smoltyr               -26.1         -28.6   0.721        0.756 -0.362    Joint SST + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_sm…          1         -0.895         -2.47  TRUE           
# 14 X15_salmonSharkBSAI_predAK_adultyr               -22.8         -25.8   0.654        0.710  0.0270   AltPrey Only (1 Prey) sitkaHerring_EGoA_adultyr                             1         -0.572         -2.99  TRUE           
# 15 X15_salmonSharkGoA_predAK_smoltyr                -23.4         -24.1   0.664        0.677 -0.146    Joint SST + 1 AltPrey capelin_WGoA_smoltyr                                  1         -0.270         -0.733 TRUE           
# 16 X15_salmonSharkGoA_predAK_adultyr                -22.7         -22.9   0.653        0.661 -0.000524 AltPrey Only (1 Prey) capelin_WGoA_adultyr                                  1         -0.0643        -0.121 TRUE           
