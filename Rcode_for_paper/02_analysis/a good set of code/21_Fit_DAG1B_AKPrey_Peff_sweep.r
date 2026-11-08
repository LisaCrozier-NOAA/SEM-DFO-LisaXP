
#21_Fit_DAG1B_AKPrey_Peff_sweep.r -- NEW VERSION


library(tidyverse)
library(lavaan)

proj_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_out_dir <- file.path(proj_dir, "Routput_for_paper/data")
raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")

raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]
df_24         <- raw_df_master %>% filter(complete.cases(.))

# Blacklist general exclusions (leaving Killer Whales in)
blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")

scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

df_base <- df_24 %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# Structural anchors for DAG 1B
var_growth   <- "X06_DFA_IGF_mu_smoltyr"
var_prey_ncc <- "X05_DFA_abundSardine_smoltyr"
var_cpue     <- "X07_cpue_IntSprJunHW_smoltyr"
var_hake     <- "X09_DFA_HakeAge5Plus_smoltyr"
var_sar      <- "X16_DFA_sar_UCsprCh_smoltyr"
var_prey_ak  <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"

lookup_file <- file.path(data_out_dir, "guild_predator_lookup_table_24yr.csv")
lookup_df   <- read_csv(lookup_file, show_col_types = FALSE)

pred_candidates <- lookup_df %>%
  filter(region == "AK" | (region == "NCC" & stage_suffix == "adultyr") | grepl("DFO", pred_data_col)) %>%
  pull(pred_data_col) %>%
  unique()

pred_candidates <- pred_candidates[!sapply(pred_candidates, function(x) any(sapply(blacklist, grepl, x)))]
pred_candidates <- intersect(pred_candidates, names(df_base))

# ------------------------------------------------------------------------------
# Isolated Sweep Function with Killer Whale Specialization Exception
# ------------------------------------------------------------------------------
run_clean_predator_sweep <- function(pred_col, df_raw, df_scaled_in, lookup_row, max_combo_size = 2) {
  
  P_norm <- df_scaled_in[[pred_col]]
  
  # Check if salmon specialist (Killer Whales: Alternate prey = NULL)
  is_salmon_specialist <- grepl("Killer_whale", pred_col, ignore.case = TRUE)
  
  # Temperature overlap
  spatial_col <- trimws(as.character(lookup_row$spatial_overlap_col[1]))
  has_temp    <- (!is.na(spatial_col) && spatial_col != "" && spatial_col %in% names(df_raw))
  P_temp      <- if (has_temp) P_norm * scale_custom_bounds(df_raw[[spatial_col]], a = 0.2, b = 1.0) else NULL
  
  # Alternate prey candidates (suppressed for killer whales)
  if (is_salmon_specialist) {
    prey_cands <- character(0)
  } else {
    all_cands  <- unlist(strsplit(lookup_row$all_candidates_str[1], ", "))
    prey_cands <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
    prey_cands <- intersect(prey_cands, names(df_scaled_in))
  }
  
  syntax_template <- paste0(
    var_growth, " ~ ", var_prey_ncc, "\n",
    var_cpue,   " ~ ", var_hake, "\n",
    var_sar,    " ~ ", var_growth, " + ", var_cpue, " + ", var_prey_ak, " + TARGET_PRED"
  )
  
  # Evaluate inside an isolated data frame clone
  fit_variant <- function(p_vec, var_label, prey_label) {
    if (is.na(sd(p_vec)) || sd(p_vec) < 1e-6) return(NULL)
    
    d_fit <- df_scaled_in
    d_fit$TARGET_PRED <- p_vec
    
    fit <- tryCatch(
      sem(syntax_template, data = d_fit, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
    
    pe  <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
    rsq <- inspect(fit, "r2")
    
    p_est <- pe %>% filter(lhs == var_sar, rhs == "TARGET_PRED", op == "~")
    g_est <- pe %>% filter(lhs == var_sar, rhs == var_growth, op == "~")
    a_est <- pe %>% filter(lhs == var_sar, rhs == var_prey_ak, op == "~")
    
    if (nrow(p_est) == 0) return(NULL)
    
    list(
      variant    = var_label,
      prey       = prey_label,
      aic        = as.numeric(fitMeasures(fit, "aic")),
      r2         = as.numeric(rsq[[var_sar]]),
      b_pred     = p_est$std.all[1],
      p_pred     = p_est$pvalue[1],
      b_growth   = g_est$std.all[1],
      p_growth   = g_est$pvalue[1],
      b_preyak   = a_est$std.all[1],
      p_preyak   = a_est$pvalue[1]
    )
  }
  
  # 1. Raw Baseline
  base_res <- fit_variant(P_norm, "Raw Baseline", "None")
  if (is.null(base_res)) return(NULL)
  
  all_runs <- list(base_res)
  
  # 2. Temp Only
  if (has_temp && !is.null(P_temp)) {
    t_res <- fit_variant(P_temp, "Temp Only", "None")
    if (!is.null(t_res)) all_runs[[length(all_runs) + 1]] <- t_res
  }
  
  # 3. AltPrey Only & Joint Temp + AltPrey (skipped for killer whales)
  if (length(prey_cands) > 0) {
    max_k <- min(length(prey_cands), max_combo_size)
    for (k in 1:max_k) {
      combos <- combn(prey_cands, k, simplify = FALSE)
      for (cb in combos) {
        raw_f <- if (length(cb) == 1) df_scaled_in[[cb[1]]] else rowSums(df_scaled_in[, cb, drop = FALSE], na.rm = TRUE)
        gate_alt <- 1.0 - 1.0 * scale_custom_bounds(raw_f, a = 0.2, b = 1.0)
        
        ap_res <- fit_variant(P_norm * gate_alt, paste0("AltPrey Only (", k, ")"), paste(cb, collapse = " + "))
        if (!is.null(ap_res)) all_runs[[length(all_runs) + 1]] <- ap_res
        
        if (has_temp && !is.null(P_temp)) {
          jt_res <- fit_variant(P_temp * gate_alt, paste0("Joint Temp + AltPrey (", k, ")"), paste(cb, collapse = " + "))
          if (!is.null(jt_res)) all_runs[[length(all_runs) + 1]] <- jt_res
        }
      }
    }
  }
  
  df_runs <- bind_rows(all_runs) %>%
    mutate(delta_aic = round(aic - base_res$aic, 2))
  
  best_abs <- df_runs %>% arrange(aic) %>% slice(1)
  
  best_neg_df <- df_runs %>% filter(b_pred < 0) %>% arrange(aic)
  best_neg <- if (nrow(best_neg_df) > 0) best_neg_df %>% slice(1) else NULL
  
  tibble(
    predator         = pred_col,
    spatial_col      = ifelse(has_temp, spatial_col, "None"),
    base_aic         = round(base_res$aic, 2),
    base_beta_pred   = round(base_res$b_pred, 3),
    base_p_pred      = round(base_res$p_pred, 4),
    base_beta_growth = round(base_res$b_growth, 3),
    base_p_growth    = round(base_res$p_growth, 4),
    
    # Best Absolute Model
    abs_variant      = best_abs$variant,
    abs_prey         = best_abs$prey,
    abs_aic          = round(best_abs$aic, 2),
    abs_delta_aic    = best_abs$delta_aic,
    abs_r2           = round(best_abs$r2, 3),
    abs_beta_pred    = round(best_abs$b_pred, 3),
    abs_p_pred       = round(best_abs$p_pred, 4),
    abs_beta_growth  = round(best_abs$b_growth, 3),
    abs_p_growth     = round(best_abs$p_growth, 4),
    abs_beta_preyak  = round(best_abs$b_preyak, 3),
    abs_p_preyak     = round(best_abs$p_preyak, 4),
    
    # Best Negative Beta Model
    neg_variant      = if (!is.null(best_neg)) best_neg$variant else "None",
    neg_prey         = if (!is.null(best_neg)) best_neg$prey else "None",
    neg_aic          = if (!is.null(best_neg)) round(best_neg$aic, 2) else NA_real_,
    neg_delta_aic    = if (!is.null(best_neg)) best_neg$delta_aic else NA_real_,
    neg_beta_pred    = if (!is.null(best_neg)) round(best_neg$b_pred, 3) else NA_real_,
    neg_p_pred       = if (!is.null(best_neg)) round(best_neg$p_pred, 4) else NA_real_
  )
}

# ------------------------------------------------------------------------------
# Run Clean Sweep
# ------------------------------------------------------------------------------
clean_sweep_list <- list()

for (i in seq_along(pred_candidates)) {
  p_name <- pred_candidates[i]
  p_row  <- lookup_df %>% filter(pred_data_col == p_name) %>% slice(1)
  
  res <- run_clean_predator_sweep(p_name, df_raw = df_24, df_scaled_in = df_base, lookup_row = p_row)
  if (!is.null(res)) {
    clean_sweep_list[[length(clean_sweep_list) + 1]] <- res
    cat(sprintf("[%d/%d] %s | Variant: %s | dAIC = %.2f | Beta: %.3f (p = %.3f)\n", 
                i, length(pred_candidates), p_name, res$abs_variant, res$abs_delta_aic, res$abs_beta_pred, res$abs_p_pred))
  }
}

clean_sweep_df <- bind_rows(clean_sweep_list) %>% arrange(abs_aic)

# Print candidates with negative predator slope
valid_clean <- clean_sweep_df %>%
  filter(abs_beta_pred < 0) %>%
  arrange(abs_aic) %>%
  select(predator, abs_variant, abs_prey, abs_aic, abs_delta_aic, abs_beta_pred, abs_p_pred, abs_beta_growth, abs_p_growth, abs_beta_preyak)

cat("\n=== CLEAN SWEEP: NEGATIVE PREDATOR CANDIDATES (DAG 1B) ===\n")
print(as.data.frame(valid_clean))


reversed_signs <- clean_sweep_df %>%
  filter(abs_beta_pred > 0) %>%
  arrange(abs_aic) %>%
  select(predator, abs_variant, abs_prey, abs_aic, abs_delta_aic, abs_beta_pred, abs_p_pred, abs_beta_growth, abs_p_growth, abs_beta_preyak)

print(as.data.frame(reversed_signs))
#                                        predator      abs_variant                                                          abs_prey abs_aic abs_delta_aic abs_beta_pred abs_p_pred abs_beta_growth abs_p_growth abs_beta_preyak
# 1             X10_harbour_seal_DFO_WCVI_smoltyr AltPrey Only (2) X04_marketsquid_GAM_smoltyr + X05_eulachon_during_chinook_smoltyr  -49.57         -7.49         0.476     0.0010           0.072       0.6414           0.019
# 2        X10_harbour_seal_DFO_Coastwide_smoltyr AltPrey Only (2)                 X05_anchovy_GAM_smoltyr + X05_herring_NCC_smoltyr  -48.66         -6.10         0.426     0.0021           0.617       0.0000           0.348
# 3 X15_DFA_halibutBiomassAge8plus_predAK_adultyr AltPrey Only (1)                                         sitkaHerring_EGoA_adultyr  -45.06         -0.58         0.346     0.0358           0.736       0.0000           0.612

clean_sweep_df %>%
  filter(predator %in% c("X10_harbour_seal_DFO_WCVI_smoltyr","X10_harbour_seal_DFO_Coastwide_smoltyr","X15_DFA_halibutBiomassAge8plus_predAK_adultyr")) %>%
  select(predator, abs_variant, abs_prey, abs_aic,base_aic, neg_aic,neg_beta_pred,neg_delta_aic,neg_p_pred,neg_prey,neg_variant)
#   predator                                      abs_variant      abs_prey                                                          abs_aic base_aic neg_aic neg_beta_pred neg_delta_aic neg_p_pred neg_prey                     neg_variant     
# 1 X10_harbour_seal_DFO_WCVI_smoltyr             AltPrey Only (2) X04_marketsquid_GAM_smoltyr + X05_eulachon_during_chinook_smoltyr   -49.6    -42.1   -43.6        -0.44          -1.53     0.103  X05_DFA_abundSardine_smoltyr AltPrey Only (1)
# 2 X10_harbour_seal_DFO_Coastwide_smoltyr        AltPrey Only (2) X05_anchovy_GAM_smoltyr + X05_herring_NCC_smoltyr                   -48.7    -42.6    NA          NA             NA       NA      None                         None            
# 3 X15_DFA_halibutBiomassAge8plus_predAK_adultyr AltPrey Only (1) sitkaHerring_EGoA_adultyr                                           -45.1    -44.5   -44.4        -0.283          0.03     0.0562 capelin_WGoA_adultyr         AltPrey Only (1)

print(as.data.frame(valid_clean))
#                                   predator              abs_variant                                         abs_prey abs_aic abs_delta_aic abs_beta_pred abs_p_pred abs_beta_growth abs_p_growth abs_beta_preyak
# 1          X15_pollockGOA_SSB_lag1_smoltyr         AltPrey Only (2)    capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr  -47.90         -6.47        -0.489     0.0040           0.813       0.0000           0.398
# 2               X15_pollockGOA_SSB_smoltyr         AltPrey Only (2)    capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr  -47.34         -3.44        -0.463     0.0062           0.776       0.0000           0.409
# 3            X11_ssl_seak_pup_pred_smoltyr         AltPrey Only (2)    capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr  -46.04         -2.52        -0.366     0.0172           0.675       0.0000           0.350
# 4     X15_sablefishAK_biomass3plus_smoltyr Joint Temp + AltPrey (2) sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr  -45.95         -4.28        -0.401     0.0183           0.277       0.0456           0.475
# 5  X15_DFA_sleeperSharkBSAI_predAK_smoltyr             Raw Baseline                                             None  -45.89          0.00        -0.676     0.0192          -0.235       0.4574           0.253
# 6     X15_sablefishAK_biomass4plus_smoltyr Joint Temp + AltPrey (2) sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr  -45.24         -3.91        -0.456     0.0311           0.091       0.6662           0.364
# 7          X10_Killer_whales_NR_BC_adultyr                Temp Only                                             None  -43.75         -1.59        -0.411     0.0929           0.812       0.0014           0.344

#NS preds
# 8                 X10_Harbour_s_WS_adultyr Joint Temp + AltPrey (1)              X05_eulachon_during_chinook_adultyr  -42.58         -1.33        -0.200     0.2190           0.286       0.1228           0.164
# 9              X10_Northern_f_s_WS_adultyr Joint Temp + AltPrey (1)              X05_eulachon_during_chinook_adultyr  -41.98         -0.75        -0.101     0.3527           0.485       0.0002           0.339
# 10       X15_salmonSharkGoA_predAK_smoltyr Joint Temp + AltPrey (1)                             capelin_WGoA_smoltyr  -41.95         -0.55        -0.112     0.3599           0.433       0.0007           0.263
# 11          X10_DFA_Steller_s_l_WS_adultyr                Temp Only                                             None  -41.54         -0.36        -0.155     0.5216           0.558       0.0112           0.277
# 12          X10_DFA_AllSeaLionsEMB_adultyr         AltPrey Only (1)              X05_eulachon_during_chinook_adultyr  -41.30         -0.07        -0.057     0.6897           0.461       0.0004           0.291