# ==============================================================================
# Script: 06_postprocess_ak_logistic_models.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Post-process existing sweep outputs (Goal3_AK_Predators_Logistic_Switching_All.csv),
#          select top single model per predator, compute realized net impact metrics,
#          and evaluate sequential multi-prey logistic switching models.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

all_csv_file  <- file.path(data_out_dir, "Goal3_AK_Predators_Logistic_Switching_All.csv")
sem_data_file <- file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv")

if (!file.exists(all_csv_file))  stop("Missing previous output file: ", all_csv_file)
if (!file.exists(sem_data_file)) stop("Missing dataset file: ", sem_data_file)

all_sweep_df      <- read_csv(all_csv_file, show_col_types = FALSE)
sem_data_noNA     <- read.csv(sem_data_file, stringsAsFactors = FALSE)
sem_complete_data <- sem_data_noNA %>% filter(complete.cases(.))

# Apply uniform +5.0 domain shift
sem_complete_data <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"

min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

logistic_gate <- function(a_norm, gamma) {
  1 / (1 + exp(-gamma * (a_norm - 0.5)))
}

# ------------------------------------------------------------------------------
# STEP 1: CALCULATE REALIZED IMPACT METRICS FOR SINGLE BEST MODELS
# ------------------------------------------------------------------------------
cat("Extracting best models and calculating realized net impact metrics...\n")

best_single_candidates <- all_sweep_df %>%
  group_by(predator) %>%
  arrange(aic) %>%
  slice(1) %>%
  ungroup()

best_single_results <- list()

for (i in seq_len(nrow(best_single_candidates))) {
  row_i    <- best_single_candidates[i, ]
  p_col    <- row_i$predator
  mod_col  <- row_i$moderator
  w_val    <- as.numeric(row_i$best_w)
  g_val    <- as.numeric(row_i$best_gamma)
  b_est    <- as.numeric(row_i$pred_est)
  
  P_pos <- sem_complete_data[[p_col]]
  SAR   <- sem_complete_data[[target_sar]]
  
  if (mod_col == "None (Unattenuated)" || w_val == 0) {
    P_eff <- P_pos
    damping_pct <- 0.0
  } else {
    A_raw  <- sem_complete_data[[mod_col]]
    A_norm <- min_max_scale(A_raw)
    weight_mod <- 1 - (w_val * logistic_gate(A_norm, g_val))
    P_eff  <- P_pos * weight_mod
    damping_pct <- (1 - (sum(P_eff, na.rm = TRUE) / sum(P_pos, na.rm = TRUE))) * 100
  }
  
  net_impact_ts <- b_est * P_eff
  mean_impact   <- mean(net_impact_ts, na.rm = TRUE)
  std_beta      <- b_est * (sd(P_eff, na.rm = TRUE) / sd(SAR, na.rm = TRUE))
  
  best_single_results[[i]] <- row_i %>%
    mutate(
      pred_est             = as.numeric(pred_est),
      pvalue               = as.numeric(pvalue),
      best_w               = as.numeric(best_w),
      best_gamma           = as.numeric(best_gamma),
      aic                  = as.numeric(aic),
      delta_aic            = as.numeric(delta_aic),
      std_eff_beta         = as.numeric(std_beta),
      mean_net_impact      = as.numeric(mean_impact),
      realized_damping_pct = as.numeric(damping_pct),
      selection_type       = if_else(mod_col == "None (Unattenuated)", "Unattenuated Baseline", "Best Single Moderator")
    )
}

best_single_df <- map_dfr(best_single_results, ~ .x)

# ------------------------------------------------------------------------------
# STEP 2: EVALUATE SEQUENTIAL MULTI-PREY MODELS (Top 2 Candidates)
# ------------------------------------------------------------------------------
cat("Evaluating Sequential Multi-Prey Models (Combining Top 2 Moderators)...\n")

multi_prey_results <- list()

predators_with_multi <- all_sweep_df %>%
  filter(moderator != "None (Unattenuated)") %>%
  group_by(predator) %>%
  filter(n() >= 2) %>%
  pull(predator) %>%
  unique()

for (pred_col in predators_with_multi) {
  
  top2_mods <- all_sweep_df %>%
    filter(predator == pred_col & moderator != "None (Unattenuated)") %>%
    arrange(aic) %>%
    slice(1:2)
  
  m1_col   <- top2_mods$moderator[1]
  m2_col   <- top2_mods$moderator[2]
  m1_w     <- as.numeric(top2_mods$best_w[1])
  m1_g     <- as.numeric(top2_mods$best_gamma[1])
  m2_w     <- as.numeric(top2_mods$best_w[2])
  m2_g     <- as.numeric(top2_mods$best_gamma[2])
  
  P_pos <- sem_complete_data[[pred_col]]
  A1    <- sem_complete_data[[m1_col]]
  A2    <- sem_complete_data[[m2_col]]
  SAR   <- sem_complete_data[[target_sar]]
  
  A1_norm <- min_max_scale(A1)
  A2_norm <- min_max_scale(A2)
  
  # Sequential Multiplicative Gate: Gate 1 x Gate 2
  gate1 <- 1 - (m1_w * logistic_gate(A1_norm, m1_g))
  gate2 <- 1 - (m2_w * logistic_gate(A2_norm, m2_g))
  
  P_seq_eff <- P_pos * (gate1 * gate2)
  
  df_multi <- sem_complete_data
  df_multi$P_seq_eff <- P_seq_eff
  
  sem_syntax <- paste0(target_sar, " ~ ", target_cpue, " + P_seq_eff")
  fit <- tryCatch(sem(sem_syntax, data = df_multi, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (!is.null(fit) && lavInspect(fit, "converged")) {
    pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_seq_eff")
    
    if (nrow(pe) > 0 && !is.na(pe$pvalue[1]) && pe$est[1] < 0) {
      current_aic <- as.numeric(fitMeasures(fit, "aic"))
      b_est       <- as.numeric(pe$est[1])
      p_val       <- as.numeric(pe$pvalue[1])
      
      net_impact_ts <- b_est * P_seq_eff
      mean_impact   <- mean(net_impact_ts, na.rm = TRUE)
      damping_pct   <- (1 - (sum(P_seq_eff, na.rm = TRUE) / sum(P_pos, na.rm = TRUE))) * 100
      std_beta      <- b_est * (sd(P_seq_eff, na.rm = TRUE) / sd(SAR, na.rm = TRUE))
      
      base_aic <- as.numeric(min(all_sweep_df$aic[all_sweep_df$predator == pred_col & all_sweep_df$moderator == "None (Unattenuated)"], na.rm = TRUE))
      
      multi_prey_results[[length(multi_prey_results) + 1]] <- tibble(
        predator             = as.character(pred_col),
        moderator            = paste0("Sequential: ", m1_col, " x ", m2_col),
        best_w               = m1_w,
        best_gamma           = m1_g,
        pred_est             = b_est,
        pvalue               = p_val,
        aic                  = current_aic,
        delta_aic            = current_aic - base_aic,
        std_eff_beta         = as.numeric(std_beta),
        mean_net_impact      = as.numeric(mean_impact),
        realized_damping_pct = as.numeric(damping_pct),
        selection_type       = "Sequential Multi-Prey"
      )
    }
  }
}

multi_prey_df <- map_dfr(multi_prey_results, ~ .x)

# Ensure numeric types across both data frames before binding
best_single_df <- best_single_df %>%
  mutate(across(c(best_w, best_gamma, pred_est, pvalue, aic, delta_aic, std_eff_beta, mean_net_impact, realized_damping_pct), as.numeric))

if (nrow(multi_prey_df) > 0) {
  multi_prey_df <- multi_prey_df %>%
    mutate(across(c(best_w, best_gamma, pred_est, pvalue, aic, delta_aic, std_eff_beta, mean_net_impact, realized_damping_pct), as.numeric))
}

# ------------------------------------------------------------------------------
# STEP 3: CONSOLIDATE & EXPORT
# ------------------------------------------------------------------------------
final_selected_df <- bind_rows(best_single_df, multi_prey_df) %>%
  group_by(predator) %>%
  arrange(aic) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(aic)

cat("\n======================================================\n")
cat("      FINAL SELECTED ALASKA PREDATOR MODELS SUMMARY    \n")
cat("======================================================\n")
cat(sprintf("Total Alaska Predators Summarized: %d\n", nrow(final_selected_df)))
cat(sprintf("Predators with Substantial Attenuation Improvement (delta_aic < -2.0): %d\n", 
            sum(final_selected_df$delta_aic < -2.0)))
cat("======================================================\n\n")

cat("--- FINAL BEST MODEL PER ALASKA PREDATOR ---\n")
print(final_selected_df %>% 
        select(predator, moderator, selection_type, std_eff_beta, mean_net_impact, realized_damping_pct, delta_aic,aic,base_aic), 
      n = 30)

# Export Summary Tables
write_csv(best_single_df,     file.path(data_out_dir, "Goal3_AK_Predators_Best_Single_Moderator.csv"))
write_csv(final_selected_df,  file.path(data_out_dir, "Goal3_AK_Predators_Final_Selected_Models.csv"))

cat("\nOutputs saved to:\n")
cat(" - ", file.path(data_out_dir, "Goal3_AK_Predators_Best_Single_Moderator.csv"), "\n")
cat(" - ", file.path(data_out_dir, "Goal3_AK_Predators_Final_Selected_Models.csv"), "\n")



#RESULTS----------
#ordered by initial best loglik (aic/predator w/ no discounting for extra parameters)
# predator                                          moderator                                                                     selection_type        std_eff_beta mean_net_impact realized_damping_pct delta_aic   aic base_aic
# 1 X10_Californian_s_l_2yrLead_WS_smoltyr            None (Unattenuated)                                                           Unattenuated Baseline       0.442            2.21                   0    0         58.8     58.8
# 2 X15_spinyDogfishGoA_predAK_smoltyr                Sequential: X13_sitkaHerring_EGoA_smoltyr x X13_mid_il_capelin_smoltyr        Sequential Multi-Prey      -0.474           -1.02                  41.9 -4.49e+ 0  59.7     NA  

# 3 X11_ssl_seak_pup_pred                             Sequential: X13_mid_il_capelin_smoltyr x X12_egoa_krill_smoltyr               Sequential Multi-Prey      -0.460           -0.672                 40.6 -5.08e+ 0  59.9     NA  
# 4 X10_Harbour_s_2yrLead_WS_adultyr                  Sequential: X09_DFA_HakeAge5Plus x X05_DFA_abundSardine_adultyr               Sequential Multi-Prey      -0.410           -2.05                   0   -1.42e-14  60.1     NA  
# 5 X15_ArrowtoothFlounderBiomass_predAK_smoltyr      X13_mid_il_capelin_smoltyr                                                    Best Single Moderator      -0.417           -0.918                 29.6 -3.97e+ 0  60.3     64.2
# 6 X10_Harbor_seal_CR_2yrLead_adultyr                None (Unattenuated)                                                           Unattenuated Baseline       0.385            1.93                   0    0         60.9     60.9
# 7 X15_sablefishBiomass_predAK_smoltyr               Sequential: X12_egoa_krill_smoltyr x X13_pollock_age1plus_smoltyr             Sequential Multi-Prey      -0.385           -0.372                 64.9 -3.49e+ 0  61.1     NA  
# 8 X15_DFA_sleeperSharks_adultyr                     None (Unattenuated)                                                           Unattenuated Baseline      -0.362           -1.81                   0    0         61.6     61.6
# 9 X15_sablefishRecruitment_predAK_smoltyr           None (Unattenuated)                                                           Unattenuated Baseline      -0.585           -2.92                   0    0         62.7     62.7

# 10 X15_DFA_sleeperSharks_smoltyr                     Sequential: X13_pollock_age1plus_smoltyr x X13_sitkaHerring_EGoA_smoltyr      Sequential Multi-Prey      -0.311           -1.56                   0   -4.97e-14  63.0     NA  
# 11 X11_DFA_Harbour_p_WS                              None (Unattenuated)                                                           Unattenuated Baseline      -0.291           -1.46                   0    0         63.4     63.4
# 12 X15_ArrowtoothFlounderBiomass_predAK_adultyr      Sequential: X13_pollock_age1plus_adultyr x X13_sitkaHerring_EGoA_adultyr      Sequential Multi-Prey      -0.294           -0.288                 53.3 -8.39e- 1  63.5     NA  
# 13 X15_PacificCodBiomass_predAK_smoltyr              Sequential: X13_mid_il_capelin_smoltyr x X13_pollock_age1plus_smoltyr         Sequential Multi-Prey      -0.302           -0.662                 29.6 -4.67e- 1  63.6     NA  
# 14 X15_sablefishBiomass_predAK_adultyr               Sequential: X13_pollock_age1plus_adultyr x X13_sitkaHerring_EGoA_adultyr      Sequential Multi-Prey      -0.285           -0.276                 61.4 -6.96e- 1  63.7     NA  
# 15 X15_salmonSharkGoA_predAK_adultyr                 None (Unattenuated)                                                           Unattenuated Baseline      -0.281           -1.41                   0    0         64.0     64.0
# 16 X10_Northern_f_s_2yrLead_WS_adultyr               None (Unattenuated)                                                           Unattenuated Baseline       0.258            1.29                   0    0         64.2     64.2
# 17 X15_sablefishRecruitment_predAK_adultyr           Sequential: X13_pollock_age1plus_adultyr x X13_sitkaHerring_EGoA_adultyr      Sequential Multi-Prey      -0.247           -0.341                 53.9 -3.09e- 1  64.4     NA  
# 18 X10_Californian_s_l_2yrLead_WS_adultyr            None (Unattenuated)                                                           Unattenuated Baseline       0.237            1.19                   0    0         64.5     64.5
# 19 X15_salmonSharkGoA_predAK_smoltyr                 Sequential: X13_pollock_age1plus_smoltyr x X14_pinkSalmonNorthAmerica_smoltyr Sequential Multi-Prey      -0.238           -0.308                 50.7 -4.26e-14  64.7     NA  
# 20 X15_PacificCodBiomass_predAK_adultyr              Sequential: X13_pollock_age1plus_adultyr x X13_sitkaHerring_EGoA_adultyr      Sequential Multi-Prey      -0.225           -0.218                 57.5 -7.00e- 1  64.8     NA  
# 21 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr Sequential: X13_pollock_age1plus_adultyr x X13_sitkaHerring_EGoA_adultyr      Sequential Multi-Prey      -0.223           -0.202                 52.3 -3.92e- 1  64.8     NA  
# 22 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr Sequential: X13_pollock_age1plus_smoltyr x X13_sitkaHerring_EGoA_smoltyr      Sequential Multi-Prey      -0.219           -0.201                 52.0 -4.30e- 1  64.9     NA  
# 23 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr        None (Unattenuated)                                                           Unattenuated Baseline       0.216            1.08                   0    0         64.9     64.9
# 24 X15_spinyDogfishGoA_predAK_adultyr                Sequential: X13_sitkaHerring_EGoA_adultyr x X21_sst_egoa_junjulaug_smoltyr    Sequential Multi-Prey      -0.172           -0.227                 57.1 -1.42e-14  65.6     NA  
# 25 X10_Harbour_s_2yrLead_WS_smoltyr                  Sequential: X09_DFA_HakeAge5Plus x X05_DFA_abundSardine_smoltyr               Sequential Multi-Prey      -0.163           -0.815                  0   -4.26e-14  65.7     NA  
# 26 X10_DFA_ssl.est.wholerange_2yrLead_adultyr        None (Unattenuated)                                                           Unattenuated Baseline       0.131            0.656                  0    0         65.9     65.9
# 27 X10_Harbor_seal_CR_2yrLead_smoltyr                None (Unattenuated)                                                           Unattenuated Baseline       0.0370           0.185                  0    0         66.5     66.5
# 28 X10_Northern_f_s_2yrLead_WS_smoltyr               None (Unattenuated)                                                           Unattenuated Baseline       0.0232           0.116                  0    0         66.5     66.5
# > 