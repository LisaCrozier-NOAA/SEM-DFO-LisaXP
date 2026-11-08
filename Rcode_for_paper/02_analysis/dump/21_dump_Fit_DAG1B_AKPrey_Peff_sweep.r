library(tidyverse)
library(lavaan)


#21_Fit_DAG1B_AKPrey_Peff_sweep.r -- REPLACED
# ==============================================================================
# 0. CONFIGURATION & DIRECTORIES-------
# ==============================================================================
proj_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_out_dir <- file.path(proj_dir, "Routput_for_paper/data")
raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")

raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Fixed anchors for DAG 1B
var_growth   <- "X06_DFA_IGF_mu_smoltyr"
var_prey_ncc <- "X05_DFA_abundSardine_smoltyr"
var_cpue     <- "X07_cpue_IntSprJunHW_smoltyr"
var_hake     <- "X09_DFA_HakeAge5Plus_smoltyr"
var_sar      <- "X16_DFA_sar_UCsprCh_smoltyr"
var_prey_ak  <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"

# ==============================================================================
# 1. EVALUATION FUNCTION FOR A SINGLE PREDATOR (DAG 1B TOPOLOGY)-------
# ==============================================================================
run_single_dag1b_sweep <- function(pred_col, df_raw, df_base, lookup_row, max_combo_size = 2) {
  
  if (is.na(pred_col) || !pred_col %in% names(df_base)) return(NULL)
  
  P_norm <- df_base[[pred_col]]
  
  # Thermal driver
  spatial_col <- trimws(as.character(lookup_row$spatial_overlap_col[1]))
  has_temp    <- (!is.na(spatial_col) && spatial_col != "" && spatial_col %in% names(df_raw))
  P_temp      <- if (has_temp) P_norm * scale_custom_bounds(df_raw[[spatial_col]], a = 0.2, b = 1.0) else NULL
  
  # Forage candidates
  all_cands <- unlist(strsplit(lookup_row$all_candidates_str[1], ", "))
  prey_cands <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
  prey_cands <- intersect(prey_cands, names(df_base))
  
  # Structural DAG 1B syntax generator
  build_syntax <- function(p_var) {
    paste0(
      var_growth, " ~ ", var_prey_ncc, "\n",
      var_cpue,   " ~ ", var_hake, "\n",
      var_sar,    " ~ ", var_growth, " + ", var_cpue, " + ", var_prey_ak, " + ", p_var
    )
  }
  
  # --- Raw Baseline ---
  base_fit <- tryCatch(sem(build_syntax(pred_col), data = df_base, missing = "ML", warn = FALSE), error = function(e) NULL)
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_pe  <- parameterEstimates(base_fit, standardized = TRUE)
  pe_pred  <- base_pe %>% filter(op == "~" & lhs == var_sar & rhs == pred_col)
  pe_grw   <- base_pe %>% filter(op == "~" & lhs == var_sar & rhs == var_growth)
  pe_pak   <- base_pe %>% filter(op == "~" & lhs == var_sar & rhs == var_prey_ak)
  if (nrow(pe_pred) == 0) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[var_sar]])
  
  best_abs <- list(variant = "Raw Baseline", prey = "None", aic = base_aic, delta_aic = 0,
                   b_pred = pe_pred$std.all[1], p_pred = pe_pred$pvalue[1],
                   b_growth = pe_grw$std.all[1], p_growth = pe_grw$pvalue[1],
                   b_preyak = pe_pak$std.all[1], p_preyak = pe_pak$pvalue[1], r2 = base_r2)
  best_abs_aic <- base_aic
  
  best_neg <- if (pe_pred$std.all[1] < 0) best_abs else NULL
  best_neg_aic <- if (pe_pred$std.all[1] < 0) base_aic else Inf
  
  # Internal evaluator for variant permutations
  eval_variant <- function(p_vector, var_type, prey_str) {
    if (is.na(sd(p_vector)) || sd(p_vector) < 1e-6) return()
    df_base$P_eff_cand <- p_vector
    fit_v <- tryCatch(sem(build_syntax("P_eff_cand"), data = df_base, missing = "ML", warn = FALSE), error = function(e) NULL)
    if (!is.null(fit_v) && lavInspect(fit_v, "converged")) {
      pe_v <- parameterEstimates(fit_v, standardized = TRUE)
      p_sub <- pe_v %>% filter(op == "~" & lhs == var_sar & rhs == "P_eff_cand")
      g_sub <- pe_v %>% filter(op == "~" & lhs == var_sar & rhs == var_growth)
      a_sub <- pe_v %>% filter(op == "~" & lhs == var_sar & rhs == var_prey_ak)
      
      if (nrow(p_sub) > 0 && !is.na(p_sub$pvalue[1])) {
        c_aic <- as.numeric(fitMeasures(fit_v, "aic"))
        res_obj <- list(variant = var_type, prey = prey_str, aic = c_aic, delta_aic = c_aic - base_aic,
                        b_pred = p_sub$std.all[1], p_pred = p_sub$pvalue[1],
                        b_growth = g_sub$std.all[1], p_growth = g_sub$pvalue[1],
                        b_preyak = a_sub$std.all[1], p_preyak = a_sub$pvalue[1],
                        r2 = as.numeric(inspect(fit_v, "r2")[[var_sar]]))
        
        if (c_aic < best_abs_aic) { best_abs_aic <<- c_aic; best_abs <<- res_obj }
        if (p_sub$std.all[1] < 0 && c_aic < best_neg_aic) { best_neg_aic <<- c_aic; best_neg <<- res_obj }
      }
    }
  }
  
  # --- Variant 1: Temp Only ---
  if (has_temp && !is.null(P_temp)) {
    eval_variant(P_temp, "Temp Only", "None")
  }
  
  # --- Variants 2 & 3: AltPrey Only and Joint (Temp + AltPrey) ---
  if (length(prey_cands) > 0) {
    max_k <- min(length(prey_cands), max_combo_size)
    for (k in 1:max_k) {
      combos <- combn(prey_cands, k, simplify = FALSE)
      for (cb in combos) {
        raw_f <- if (length(cb) == 1) df_base[[cb[1]]] else rowSums(df_base[, cb, drop = FALSE], na.rm = TRUE)
        gate_alt <- 1 - 1.0 * scale_custom_bounds(raw_f, a = 0.2, b = 1.0)
        
        # AltPrey Only
        eval_variant(P_norm * gate_alt, paste0("AltPrey Only (", k, ")"), paste(cb, collapse = " + "))
        
        # Joint Temp + AltPrey
        if (has_temp && !is.null(P_temp)) {
          eval_variant(P_temp * gate_alt, paste0("Joint Temp + AltPrey (", k, ")"), paste(cb, collapse = " + "))
        }
      }
    }
  }
  
  tibble(
    predator            = pred_col,
    spatial_col         = ifelse(has_temp, spatial_col, "None"),
    base_aic            = round(base_aic, 2),
    base_beta_pred      = round(pe_pred$std.all[1], 3),
    base_p_pred         = round(pe_pred$pvalue[1], 4),
    base_beta_growth    = round(pe_grw$std.all[1], 3),
    base_p_growth       = round(pe_grw$pvalue[1], 4),
    
    # Best Absolute Model
    abs_variant         = best_abs$variant,
    abs_prey            = best_abs$prey,
    abs_aic             = round(best_abs$aic, 2),
    abs_delta_aic       = round(best_abs$delta_aic, 2),
    abs_r2              = round(best_abs$r2, 3),
    abs_beta_pred       = round(best_abs$b_pred, 3),
    abs_p_pred          = round(best_abs$p_pred, 4),
    abs_beta_growth     = round(best_abs$b_growth, 3),
    abs_p_growth        = round(best_abs$p_growth, 4),
    abs_beta_preyak     = round(best_abs$b_preyak, 3),
    abs_p_preyak        = round(best_abs$p_preyak, 4),
    
    # Best Negative Beta Model
    neg_variant         = if (!is.null(best_neg)) best_neg$variant else "None",
    neg_prey            = if (!is.null(best_neg)) best_neg$prey else "None",
    neg_aic             = if (!is.null(best_neg)) round(best_neg$aic, 2) else NA_real_,
    neg_delta_aic       = if (!is.null(best_neg)) round(best_neg$delta_aic, 2) else NA_real_,
    neg_beta_pred       = if (!is.null(best_neg)) round(best_neg$b_pred, 3) else NA_real_,
    neg_p_pred          = if (!is.null(best_neg)) round(best_neg$p_pred, 4) else NA_real_,
    neg_beta_growth     = if (!is.null(best_neg)) round(best_neg$b_growth, 3) else NA_real_,
    neg_p_growth        = if (!is.null(best_neg)) round(best_neg$p_growth, 4) else NA_real_
  )
}

#Step 3: Master Execution Loop for 24-year and 19-year Runs-----------

# Execute sweeps for 24yr and 19yr AK predators under DAG 1B
for (cur_years in c("24yr", "19yr")) {
  
  cat(sprintf("\n>>> STARTING DAG 1B MODERATION SWEEP: %s <<<\n", cur_years))
  
  req_count  <- if (cur_years == "24yr") 24 else 19
  valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) >= req_count]
  
  df_sub <- raw_df_master %>% 
    select(all_of(valid_cols)) %>% 
    filter(complete.cases(.))
  
  df_scaled <- df_sub %>%
    mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))
  
  # Load lookup table
  lookup_file <- file.path(data_out_dir, paste0("guild_predator_lookup_table_", cur_years, ".csv"))
  if (!file.exists(lookup_file)) {
    warning("Lookup table not found: ", lookup_file)
    next
  }
  lookup_df <- read_csv(lookup_file, show_col_types = FALSE)
  
  # Filter candidate AK/Adult predators
  pred_candidates <- lookup_df %>%
    filter(region == "AK" | (region == "NCC" & stage_suffix == "adultyr") | grepl("DFO", pred_data_col)) %>%
    pull(pred_data_col) %>%
    unique()
  
  pred_candidates <- pred_candidates[!sapply(pred_candidates, function(x) any(sapply(blacklist, grepl, x)))]
  pred_candidates <- intersect(pred_candidates, names(df_scaled))
  
  sweep_list <- list()
  for (i in seq_along(pred_candidates)) {
    p_name <- pred_candidates[i]
    p_row  <- lookup_df %>% filter(pred_data_col == p_name) %>% slice(1)
    
    res <- run_single_dag1b_sweep(p_name, df_raw = df_sub, df_base = df_scaled, lookup_row = p_row)
    if (!is.null(res)) {
      sweep_list[[length(sweep_list) + 1]] <- res
      cat(sprintf("[%d/%d] %s | Best: %s (dAIC = %.2f) | Beta: %.3f\n", 
                  i, length(pred_candidates), p_name, res$abs_variant, res$abs_delta_aic, res$abs_beta_pred))
    }
  }
  
  if (length(sweep_list) > 0) {
    out_df <- bind_rows(sweep_list) %>% arrange(abs_aic)
    out_file <- file.path(data_out_dir, paste0("Goal3_DAG1B_Sweep_AK_", cur_years, ".csv"))
    write_csv(out_df, out_file)
    cat(sprintf("Saved %d DAG 1B sweep results to: %s\n", nrow(out_df), basename(out_file)))
  }
}


#Examine output--------
library(tidyverse)

# Path to the newly generated 24-year sweep results
sweep_out_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/data/Goal3_DAG1B_Sweep_AK_24yr.csv"

sweep_res <- read.csv(sweep_out_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. Top 15 Absolute Models (Ranked by Lowest AIC)
# ------------------------------------------------------------------------------
cat("=== TOP 15 ABSOLUTE MODELS (DAG 1B + PREY AK) ===\n")
top_abs <- sweep_res %>%
  arrange(abs_aic) %>%
  slice(1:15) %>%
  select(
    predator,
    abs_variant,
    abs_prey,
    abs_aic,
    abs_delta_aic,
    abs_r2,
    abs_beta_pred,
    abs_p_pred,
    abs_beta_growth,
    abs_p_growth,
    abs_beta_preyak,
    abs_p_preyak
  )

print(as.data.frame(top_abs))

reversed_signs <- sweep_res %>%
  filter(abs_beta_pred > 0 & abs_p_pred < 0.10) %>% # & abs_beta_growth > 0 & abs_p_growth < 0.05) %>%
  arrange(abs_aic) %>%
  select(
    predator,
    abs_variant,
    abs_prey,
    abs_aic,
    abs_delta_aic,
    abs_beta_pred,
    abs_p_pred,
    abs_beta_growth,
    abs_p_growth,
    abs_beta_preyak
  )

print(as.data.frame(reversed_signs))
#                                        predator              abs_variant                                               abs_prey abs_aic abs_delta_aic abs_beta_pred abs_p_pred abs_beta_growth abs_p_growth abs_beta_preyak
# 1 X15_DFA_halibutBiomassAge8plus_predAK_adultyr         AltPrey Only (1)                              sitkaHerring_EGoA_adultyr  -44.78         -3.34         0.667     0.0106           0.855        2e-04           0.722
# 2        X10_harbour_seal_DFO_Coastwide_smoltyr Joint Temp + AltPrey (2) X05_anchovy_GAM_smoltyr + X09_DFA_HakeAge5Plus_smoltyr  -43.73         -4.69         0.442     0.0224           0.701        4e-04           0.444
> 

# ------------------------------------------------------------------------------
# 2. Top Models that Fulfill Both Signs:
#    - Predator is negative (beta < 0, p < 0.10)
#    - Growth is positive (beta > 0, p < 0.05)
# ------------------------------------------------------------------------------
cat("\n=== MODELS WITH BOTH CORRECT SIGNS (Negative Pred & Positive Growth) ===\n")
valid_signs <- sweep_res %>%
  filter(abs_beta_pred < 0 & abs_p_pred < 0.10) %>% # & abs_beta_growth > 0 & abs_p_growth < 0.05) %>%
  arrange(abs_aic) %>%
  select(
    predator,
    abs_variant,
    abs_prey,
    abs_aic,
    abs_delta_aic,
    abs_beta_pred,
    abs_p_pred,
    abs_beta_growth,
    abs_p_growth,
    abs_beta_preyak
  )

print(as.data.frame(valid_signs))


#RESULTS-------
print(as.data.frame(valid_signs))
#                            predator              abs_variant                                           abs_prey abs_aic abs_delta_aic abs_beta_pred abs_p_pred abs_beta_growth abs_p_growth abs_beta_preyak
# 1     X11_ssl_seak_pup_pred_smoltyr         AltPrey Only (2) sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr  -51.36         -7.59        -0.514     0.0001           0.613       0.0000           0.088
# 2 X10_harbour_seal_DFO_WCVI_smoltyr         AltPrey Only (1)                            X05_anchovy_GAM_smoltyr  -45.41         -1.41        -0.444     0.0068           0.611       0.0001           0.146
# 3   X10_Killer_whales_NR_BC_adultyr Joint Temp + AltPrey (1)                X05_eulachon_during_chinook_adultyr  -45.27         -4.42        -0.536     0.0075           0.766       0.0001           0.220
# 4    X10_DFA_Steller_s_l_WS_adultyr Joint Temp + AltPrey (1)                X05_eulachon_during_chinook_adultyr  -42.82         -2.42        -0.478     0.0419           0.708       0.0013           0.138
# 5        X15_pollockGOA_SSB_smoltyr Joint Temp + AltPrey (2) sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr  -42.47         -3.36        -0.378     0.0536           0.562       0.0008           0.127

#What to Check in This Table
#abs_delta_aic: Negative values indicate by how many AIC units the moderation gate (alternate prey / temperature) improved fit over the unmoderated raw baseline.
#abs_beta_pred & abs_p_pred: Shows whether the predator retained or flipped to a negative top-down effect ($\beta < 0$).
#abs_beta_growth & abs_p_growth: Confirms whether the early-ocean juvenile growth pathway stayed robust and statistically significant ($p < 0.05$).
#abs_beta_preyak: Checks that subarctic ocean productivity continued to act as a significant positive bottom-up anchor.

#                                  predator              abs_variant                                                    abs_prey abs_aic abs_delta_aic abs_beta_pred abs_p_pred abs_beta_growth abs_p_growth abs_beta_preyak
# 1           X11_ssl_seak_pup_pred_smoltyr         AltPrey Only (2)          sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr  -51.36         -7.59        -0.514     0.0001           0.613       0.0000           0.088
# 2       X10_harbour_seal_DFO_WCVI_smoltyr         AltPrey Only (1)                                     X05_anchovy_GAM_smoltyr  -45.41         -1.41        -0.444     0.0068           0.611       0.0001           0.146
# 3         X10_Killer_whales_NR_BC_adultyr Joint Temp + AltPrey (1)                         X05_eulachon_during_chinook_adultyr  -45.27         -4.42        -0.536     0.0075           0.766       0.0001           0.220
# 4    X15_sablefishAK_biomass3plus_smoltyr Joint Temp + AltPrey (2) X13_pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr  -44.78         -2.84        -0.457     0.0106           0.105       0.5391           0.383
# 5    X15_sablefishAK_biomass4plus_smoltyr Joint Temp + AltPrey (1)                                   sitkaHerring_EGoA_smoltyr  -43.77         -1.88        -0.471     0.0218           0.009       0.9645           0.351
# 6          X10_DFA_Steller_s_l_WS_adultyr Joint Temp + AltPrey (1)                         X05_eulachon_during_chinook_adultyr  -42.82         -2.42        -0.478     0.0419           0.708       0.0013           0.138
# 7 X15_DFA_sleeperSharkBSAI_predAK_smoltyr Joint Temp + AltPrey (1)                                        capelin_WGoA_smoltyr  -42.61         -3.55        -0.548     0.0485          -0.134       0.6515           0.223
# 8              X15_pollockGOA_SSB_smoltyr Joint Temp + AltPrey (2)          sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr  -42.47         -3.36        -0.378     0.0536           0.562       0.0008           0.127
# 9                X10_Harbour_s_WS_adultyr             Raw Baseline                                                        None  -41.58          0.00        -0.285     0.0996           0.201       0.2785           0.265
# >
