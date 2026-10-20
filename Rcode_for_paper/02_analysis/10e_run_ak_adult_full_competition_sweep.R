# ==============================================================================
# Script: 10e_run_ak_adult_full_competition_sweep.R
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

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

if ("X01_habCompInd_smoltyr" %in% names(raw_df)) {
  raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr
}
if ("X01_habCompInd_adultyr" %in% names(raw_df)) {
  raw_df$X01_habCompInd_adultyr <- -1 * raw_df$X01_habCompInd_adultyr
}

if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))
guild_pred_lookup <- read_csv(file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")), show_col_types = FALSE)

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

all_cols <- names(filtered_df)

blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

candidate_preds <- all_cols[grepl("^(X10_|X15_)", all_cols) | all_cols == "X11_ssl_seak_pup_pred"]
candidate_preds <- candidate_preds[!grepl("^X10_.*smoltyr$", candidate_preds)]
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
  region <- ifelse(grepl("^X15_", pred_col) | pred_col == "X11_ssl_seak_pup_pred", "AK", "NCC_adult")
  
  spatial_col <- case_when(
    region == "AK" & stage_suffix == "smoltyr" ~ "X21_sst_egoa_junjulaug_smoltyr",
    region == "AK" & stage_suffix == "adultyr" ~ "X21_sst_egoa_junjulaug_adultyr",
    region == "NCC_adult" ~ paste0("X01_habCompInd_", stage_suffix)
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


#RESULTS-------
# AK & ADULT FULL 4-WAY COMPETITION RESULTS (SAR ~ Peff) - 24yr DATASET         
# predator                                          base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_variant_type      abs_prey        abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved
# 1 X10_Harbour_s_2yrLead_WS_adultyr                   -3.87         -9.39     0.443        0.570   -0.546  Joint SST + 2 AltPrey X12_DFA_biomas…          0.9        -2.32         -5.52  TRUE           
# 2 X10_Californian_s_l_2yrLead_WS_adultyr              0.507        -8.95     0.363        0.573    0.260  AltPrey Only (1 Prey) X09_DFA_HakeAg…          1           1.24         -9.46  TRUE           
# 3 X15_ArrowtoothFlounderBiomass_predAK_adultyr        0.992        -7.76     0.346        0.608   -0.187  AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.479        -8.75  TRUE           
# 4 X15_ArrowtoothFlounderBiomass_predAK_smoltyr        0.299        -5.93     0.369        0.503   -0.202  Joint SST + 1 AltPrey X13_capelin_WG…          1          -0.496        -6.23  TRUE           
# 5 X15_sablefishBiomass_predAK_smoltyr                 1.99         -5.63     0.316        0.504    0.131  Joint SST + 2 AltPrey X13_pollock_ag…          1          -2.35         -7.62  TRUE           
# 6 X15_DFA_sleeperSharks_adultyr                      -2.38         -5.29     0.453        0.529   -0.285  Joint SST + 2 AltPrey X13_pollock_ag…          1          -0.541        -2.92  TRUE           
# 7 X15_DFA_sleeperSharks_smoltyr                      -1.05         -4.82     0.409        0.495   -0.253  Joint SST + 1 AltPrey X13_capelin_WG…          1          -0.586        -3.77  TRUE           
# 8 X11_ssl_seak_pup_pred                               1.00         -4.08     0.354        0.458    0.183  Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.544        -5.08  TRUE           
# 9 X15_salmonSharkGoA_predAK_smoltyr                   2.46         -3.90     0.298        0.456    0.0565 Joint SST + 1 AltPrey X13_capelin_WG…          1          -1.63         -6.36  TRUE           
# 10 X10_DFA_ssl.est.wholerange_2yrLead_adultyr          1.94         -3.83     0.319        0.483    0.132  AltPrey Only (1 Prey) X09_DFA_HakeAg…          1           1.06         -5.77  TRUE           
# 11 X10_Northern_f_s_2yrLead_WS_adultyr                 0.211        -3.16     0.367        0.455    0.279  AltPrey Only (2 Prey) X05_DFA_abundS…          1           0.545        -3.37  TRUE           
# 12 X15_spinyDogfishGoA_predAK_smoltyr                  2.50         -1.59     0.301        0.384    0.0226 Joint SST + 2 AltPrey X13_pollock_ag…          1          -1.26         -4.08  TRUE           
# 13 X15_PacificCodBiomass_predAK_smoltyr                0.0557       -1.19     0.377        0.372    0.270  Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.706        -1.25  TRUE           
# 14 X15_sablefishBiomass_predAK_adultyr                 1.39         -1.11     0.333        0.432    0.228  AltPrey Only (2 Prey) X13_pollock_ag…          1          -1.32         -2.51  TRUE           
# 15 X10_AllSeaLionsEMB_adultyr                          2.49         -0.917    0.301        0.398    0.0301 AltPrey Only (2 Prey) X04_marketsqui…          1           1.23         -3.41  TRUE           
# 16 X10_AllSeaLionsEMB_2025_raw_sc_adultyr              2.49         -0.917    0.301        0.398    0.0301 AltPrey Only (2 Prey) X04_marketsqui…          1           1.23         -3.41  TRUE           
# 17 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr   1.96         -0.911    0.321        0.416   -0.138  Joint SST + 1 AltPrey X13_capelin_WG…          1          -0.493        -2.87  TRUE           
# 18 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr   2.33         -0.591    0.308        0.382   -0.0802 Joint SST + 1 AltPrey X13_capelin_WG…          1          -0.531        -2.92  TRUE           
# 19 X15_salmonSharkGoA_predAK_adultyr                   1.60          0.0163   0.328        0.399   -0.222  Joint SST + 2 AltPrey X13_pollock_ag…          0.9        -1.36         -1.58  TRUE           
# 20 X15_PacificCodBiomass_predAK_adultyr                2.46          0.211    0.302        0.401    0.0389 AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.242        -2.25  TRUE           
# 21 X15_sablefishRecruitment_predAK_smoltyr             2.36          1.12     0.306        0.339   -0.0674 Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.278        -1.24  TRUE           
# 22 X15_sablefishRecruitment_predAK_adultyr             2.32          2.00     0.305        0.324    0.0774 Joint SST + 2 AltPrey X13_pollock_ag…          1          -0.646        -0.319 TRUE           


#19yr-----------
# predator                                          base_aic abs_final_aic base_r2 abs_final_r2 base_beta abs_variant_type      abs_prey        abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved
#  1 X15_salmonSharkBSAI_predAK_smoltyr                   1.13        -20.0     0.311        0.778  -0.229   Joint SST + 2 AltPrey X13_pollock_ag…          1          -2.37        -21.2   TRUE           
# 2 X10_Californian_s_l_2yrLead_WS_adultyr              -0.730       -14.7     0.379        0.698   0.362   AltPrey Only (2 Prey) X09_DFA_HakeAg…          0.8         1.66        -14.0   TRUE           
# 3 X10_DFA_ssl.est.wholerange_2yrLead_adultyr           1.23        -12.8     0.313        0.678   0.198   AltPrey Only (2 Prey) X09_DFA_HakeAg…          1           2.04        -14.0   TRUE           
# 4 X15_DFA_sleeperSharks_smoltyr                       -8.80        -10.8     0.620        0.657  -0.613   Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.969        -2.01  TRUE           
# 5 X15_spinyDogfishBSAI_predAK_smoltyr                  1.34        -10.7     0.315        0.629   0.188   Joint SST + 2 AltPrey X13_capelin_WG…          1          -1.44        -12.0   TRUE           
# 6 X15_ArrowtoothFlounderBiomass_predAK_adultyr         0.869       -10.5     0.331        0.712  -0.188   AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.627       -11.4   TRUE           
# 7 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr   -1.44        -10.4     0.420        0.693  -0.377   AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.662        -8.98  TRUE           
# 8 X15_DFA_sleeperSharks_adultyr                       -5.28         -9.94    0.533        0.657  -0.412   AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.587        -4.66  TRUE           
# 9 X10_Harbour_s_2yrLead_WS_adultyr                    -3.77         -8.81    0.478        0.617  -0.578   Joint SST + 2 AltPrey X09_DFA_HakeAg…          0.7        -1.66         -5.04  TRUE           
# 10 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr    0.471        -8.49    0.347        0.591  -0.250   Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.814        -8.96  TRUE           
# 11 X15_ArrowtoothFlounderBiomass_predAK_smoltyr        -0.166        -7.50    0.373        0.565  -0.247   Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.675        -7.34  TRUE           
# 12 X11_ssl_seak_pup_pred                               -3.90         -7.13    0.502        0.577   0.457   AltPrey Only (1 Prey) X13_pollock_ag…          0.5         0.642        -3.23  TRUE           
# 13 X15_salmonSharkGoA_predAK_adultyr                    0.672        -6.14    0.318        0.563  -0.294   AltPrey Only (2 Prey) X13_pollock_ag…          1          -2.84         -6.81  TRUE           
# 14 X15_salmonSharkGoA_predAK_smoltyr                    2.08         -5.84    0.273        0.521   0.0260  Joint SST + 2 AltPrey X13_sitkaHerri…          1          -1.98         -7.91  TRUE           
# 15 X15_sablefishBiomass_predAK_smoltyr                  1.49         -5.82    0.300        0.530   0.150   Joint SST + 2 AltPrey X13_pollock_ag…          1          -2.36         -7.31  TRUE           
# 16 X15_PacificCodBiomass_predAK_smoltyr                -0.342        -4.16    0.370        0.453   0.293   Joint SST + 2 AltPrey X13_sitkaHerri…          1          -1.08         -3.82  TRUE           
# 17 X15_salmonSharkBSAI_predAK_adultyr                   1.42         -3.54    0.309        0.502  -0.190   AltPrey Only (2 Prey) X13_pollock_ag…          1          -0.612        -4.97  TRUE           
# 18 X10_Northern_f_s_2yrLead_WS_adultyr                  0.459        -2.69    0.332        0.437   0.301   AltPrey Only (2 Prey) X09_DFA_HakeAg…          1           0.677        -3.15  TRUE           
# 19 X10_AllSeaLionsEMB_2025_raw_sc_adultyr               2.05         -2.35    0.277        0.444   0.0450  AltPrey Only (2 Prey) X09_DFA_HakeAg…          1           2.92         -4.41  TRUE           
# 20 X10_AllSeaLionsEMB_adultyr                           2.05         -2.35    0.277        0.444   0.0450  AltPrey Only (2 Prey) X09_DFA_HakeAg…          1           2.92         -4.41  TRUE           
# 21 X15_PacificCodBiomass_predAK_adultyr                 2.09         -1.67    0.276        0.483   0.00597 AltPrey Only (2 Prey) X13_pollock_ag…          1          -0.562        -3.76  TRUE           
# 22 X15_sablefishBiomass_predAK_adultyr                  0.806        -1.14    0.326        0.430   0.258   AltPrey Only (2 Prey) X13_pollock_ag…          1          -1.41         -1.95  TRUE           
# 23 X15_spinyDogfishGoA_predAK_smoltyr                   2.09         -0.943   0.276        0.375  -0.00153 Joint SST + 2 AltPrey X13_sitkaHerri…          1          -0.451        -3.03  TRUE           
# 24 X15_spinyDogfishBSAI_predAK_adultyr                  1.96         -0.772   0.277        0.444   0.0886  AltPrey Only (1 Prey) X13_capelin_WG…          1          -0.371        -2.73  TRUE           
# 25 X15_sablefishRecruitment_predAK_smoltyr              1.81          0.240   0.287        0.345  -0.0994  Joint SST + 2 AltPrey X13_pollock_ag…          1          -0.603        -1.57  TRUE           
# 26 X15_sablefishRecruitment_predAK_adultyr              1.81          1.16    0.288        0.320   0.0994  AltPrey Only (1 Prey) X13_pollock_ag…          1           1.48         -0.650 TRUE           
# 
