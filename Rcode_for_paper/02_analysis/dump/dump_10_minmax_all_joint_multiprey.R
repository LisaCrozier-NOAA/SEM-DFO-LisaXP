# ==============================================================================
# Script: 10_minmax_all_joint_multiprey.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Pure [0, 1] Min-Max scaling applied across ALL variables (predators & prey).
#          Joint multi-prey combinatorial sweep + sorted AIC table & 4-panel plots.
# Output: Rcode_for_paper/Routput_for_paper/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & RAW DATA (NO DOMAIN SHIFT)
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
plot_dir     <- file.path(master_out, "plots")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_file <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
lookup_file   <- file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv")
sem_orig     <- read.csv("data_Lisa/sem_master_data.csv", stringsAsFactors = FALSE)


# LOAD RAW UN-SHIFTED DATA
raw_sem_data <- read.csv(sem_data_file, stringsAsFactors = FALSE) %>% 
  filter(complete.cases(.))
sort(names(raw_sem_data))

raw_sem_data$X13_sitkaHerring_EGoA_smoltyr <- sem_orig$sitkaHerring_EGoA
raw_sem_data$X13_sitkaHerring_EGoA_adultyr <- sem_orig$sitkaHerring_EGoA[c(3:24,24,24)]


pred_lookup <- read.csv(lookup_file, stringsAsFactors = FALSE)
pred_lookup %>% select(pred_data_col) %>% unique()
pred_lookup <- pred_lookup %>% filter(!grepl("canaryRockfish", pred_data_col))

ncc_prey <- pred_lookup %>%
  filter(region=="NCC") %>%
  select(altprey1_data_col, altprey2_data_col, altprey3_data_col) %>%
  pivot_longer(cols = everything(), values_to = "target_var") %>%
  filter(!is.na(target_var) & target_var != "<NA>") %>%
  pull(target_var) %>%
  unique()


X08_prey<-ncc_prey[1:3]
X09_prey<-ncc_prey[1:5]
X10_prey<-ncc_prey[1:5]

NCCprey<-

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
temp_ak     <- "X21_sst_egoa_junjulaug_smoltyr"

# Universal Min-Max Helper Function [0, 1]
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Apply [0, 1] Min-Max Scaling across ALL numeric indicators EXCEPT targets & Year
minmax_complete_data <- raw_sem_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# ------------------------------------------------------------------------------
# STEP 1: JOINT MULTI-PREY SWEEP WITH PURE [0, 1] SCALED VARIABLES
# ------------------------------------------------------------------------------
find_best_joint_prey_minmax <- function(pred_col, df_model, max_combo_size = 2) {
  
  P_norm <- df_model[[pred_col]] # Now strictly in [0, 1]
  
  row_data <- pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  preys <- c(row_data$altprey1_data_col, row_data$altprey2_data_col, row_data$altprey3_data_col)
  preys <- unique(preys[!is.na(preys) & preys != ""])
  candidates <- intersect(unique(c(preys, temp_ak)), names(df_model))
  
  # Base Unattenuated Model
  base_fit <- tryCatch(
    sem(paste0(target_sar, " ~ ", target_cpue, " + ", pred_col), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), 
    error = function(e) NULL
  )
  
  if (is.null(base_fit) || !lavInspect(base_fit, "converged")) return(NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_pe  <- parameterEstimates(base_fit) %>% filter(op == "~" & rhs == pred_col)
  base_est <- as.numeric(base_pe$est[1])
  base_p   <- as.numeric(base_pe$pvalue[1])
  
  best_res <- list(
    predator = pred_col, selected_prey = "None (Unattenuated)",
    weights = c(), aic = base_aic, base_aic = base_aic, delta_aic = 0.0,
    beta = base_est, pval = base_p, P_eff = P_norm, P_eff_stage1 = P_norm
  )
  
  if (length(candidates) == 0) return(best_res)
  
  # Candidate matrices (already in [0, 1])
  Anorm_mat <- as.matrix(df_model[, candidates, drop = FALSE])
  w_grid    <- seq(0.1, 1.0, by = 0.1)
  
  # Joint search across 1-prey and 2-prey combinations
  for (k in 1:min(length(candidates), max_combo_size)) {
    combos <- combn(candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      w_combo_grid <- expand.grid(rep(list(w_grid), k))
      
      for (g_idx in 1:nrow(w_combo_grid)) {
        w_vals <- as.numeric(w_combo_grid[g_idx, ])
        
        w1_gate  <- 1 - w_vals[1] * Anorm_mat[, cb[1]]
        P_stage1 <- P_norm * w1_gate
        
        if (k > 1) {
          w2_gate <- 1 - w_vals[2] * Anorm_mat[, cb[2]]
          P_joint <- P_norm * (w1_gate * w2_gate)
        } else {
          P_joint <- P_stage1
        }
        
        df_model$P_joint_eff <- P_joint
        
        fit <- tryCatch(
          sem(paste0(target_sar, " ~ ", target_cpue, " + P_joint_eff"), data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), 
          error = function(e) NULL
        )
        
        if (!is.null(fit) && lavInspect(fit, "converged")) {
          pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_joint_eff")
          
          if (nrow(pe) > 0 && !is.na(pe$pvalue[1]) && pe$est[1] < 0) {
            curr_aic <- as.numeric(fitMeasures(fit, "aic"))
            
            if (curr_aic < best_res$aic) {
              best_res <- list(
                predator = pred_col, selected_prey = cb,
                weights = w_vals, aic = curr_aic, base_aic = base_aic,
                delta_aic = curr_aic - base_aic,
                beta = as.numeric(pe$est[1]), pval = as.numeric(pe$pvalue[1]),
                P_eff = P_joint, P_eff_stage1 = P_stage1
              )
            }
          }
        }
      }
    }
  }
  
  return(best_res)
}

# ------------------------------------------------------------------------------
# STEP 2: EXECUTE SWEEP & SORT BY AIC
# ------------------------------------------------------------------------------
cat("Executing [0, 1] Min-Max Joint Sweep across Alaska Predators...\n")

ak_candidates <- pred_lookup %>% 
  filter(grepl("^(X10_|X11_|X15_)", pred_data_col)) %>%
  filter(pred_data_col %in% names(minmax_complete_data)) %>%
  distinct(pred_data_col, .keep_all = TRUE) %>%
  pull(pred_data_col)

joint_sweep_results <- list()

for (pred_col in ak_candidates) {
  res <- find_best_joint_prey_minmax(pred_col, minmax_complete_data)
  if (!is.null(res)) {
    
    P_norm <- minmax_complete_data[[pred_col]]
    SAR    <- minmax_complete_data[[target_sar]]
    
    net_impact <- res$beta * res$P_eff
    mean_net   <- mean(net_impact, na.rm = TRUE)
    damp_pct   <- (1 - (sum(res$P_eff) / sum(P_norm))) * 100
    std_beta   <- res$beta * (sd(res$P_eff) / sd(SAR))
    
    formatted_weights <- if (length(res$weights) > 0 && is.numeric(res$weights)) {
      paste(round(res$weights, 2), collapse = "; ")
    } else {
      "None"
    }
    
    joint_sweep_results[[length(joint_sweep_results) + 1]] <- tibble(
      predator             = pred_col,
      selected_prey        = paste(res$selected_prey, collapse = " + "),
      fitted_weights       = formatted_weights,
      beta_peff            = res$beta,
      pvalue               = res$pval,
      std_eff_beta         = std_beta,
      mean_net_impact      = mean_net,
      realized_damping_pct = damp_pct,
      aic                  = res$aic,
      base_aic             = res$base_aic,
      delta_aic            = res$delta_aic
    )
  }
}

joint_sweep_table <- map_dfr(joint_sweep_results, ~ .x) %>%
  arrange(aic) # Sorted by AIC (Best fitting models top)

cat("\n==============================================================================\n")
cat("    UNIFORM [0, 1] MIN-MAX JOINT MODELS SORTED BY AIC (BEST FITS FIRST)       \n")
cat("==============================================================================\n\n")
print(joint_sweep_table %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, mean_net_impact, aic, delta_aic), n = 30)

write_csv(joint_sweep_table, file.path(data_out_dir, "Goal3_AK_Predators_PureMinMax_Sorted_AIC.csv"))

# ------------------------------------------------------------------------------
# STEP 3: REFINED 4-PANEL DIAGNOSTIC PLOT ON PURE [0, 1] SCALE
# ------------------------------------------------------------------------------
plot_pure_minmax_diagnostic <- function(pred_col, df_model) {
  
  res <- find_best_joint_prey_minmax(pred_col, df_model)
  if (is.null(res)) return(NULL)
  
  P_norm <- df_model[[pred_col]]
  Years  <- df_model$Year
  
  cpue_lm    <- lm(as.formula(paste(target_sar, "~", target_cpue)), data = df_model)
  sar_resids <- residuals(cpue_lm)
  
  selected_prey <- res$selected_prey
  weights       <- res$weights
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  # --- PANEL 1: RAW VS EFFECTIVE PREDATOR ([0, 1] BOUNDED) ---
  plot(Years, P_norm, type = "b", pch = 16, col = "black", ylim = c(0, 1.1),
       main = "A. Raw vs. Joint Effective Predator Index", xlab = "Year", ylab = "Predator Scale [0, 1]")
  lines(Years, res$P_eff, type = "b", pch = 17, col = "firebrick", lty = 2, lwd = 2)
  legend("topleft", legend = c("Raw Predator (P_norm)", "Joint Effective (P_eff)"),
         col = c("black", "firebrick"), pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.8)
  
  # --- PANEL 2: PREY ABUNDANCE VS WEIGHT GATE ---
  if (length(selected_prey) > 0 && selected_prey[1] != "None (Unattenuated)") {
    A1_norm <- df_model[[selected_prey[1]]]
    wt1     <- 1 - weights[1] * A1_norm
    
    plot(Years, A1_norm, type = "l", col = "steelblue", lty = 2, ylim = c(0, 1.1),
         main = "B. Prey Abundance & Weight Gate", xlab = "Year", ylab = "Proportion / Weight")
    lines(Years, wt1, col = "navy", lwd = 2)
    
    leg_txt <- c(paste("Raw", selected_prey[1]), paste("Weight W1 (w =", weights[1], ")"))
    leg_col <- c("steelblue", "navy")
    leg_lty <- c(2, 1)
    
    if (length(selected_prey) > 1) {
      A2_norm <- df_model[[selected_prey[2]]]
      wt2     <- 1 - weights[2] * A2_norm
      lines(Years, A2_norm, col = "darkorange", lty = 2)
      lines(Years, wt2, col = "sienna", lwd = 2)
      
      leg_txt <- c(leg_txt, paste("Raw", selected_prey[2]), paste("Weight W2 (w =", weights[2], ")"))
      leg_col <- c(leg_col, "darkorange", "sienna")
      leg_lty <- c(leg_lty, 2, 1)
    }
    legend("topright", legend = leg_txt, col = leg_col, lty = leg_lty, lwd = c(1,2,1,2)[1:length(leg_txt)], bty = "n", cex = 0.7)
  } else {
    plot(Years, rep(1, length(Years)), type = "l", col = "gray", ylim = c(0, 1.1),
         main = "B. No Alternate Prey Selected", xlab = "Year", ylab = "Weight")
  }
  
  # --- PANEL 3: STEPWISE ATTENUATION PROGRESSION ---
  plot(Years, P_norm, type = "l", col = "gray60", lwd = 1.5, ylim = c(0, 1.1),
       main = "C. Stepwise Attenuation Impact", xlab = "Year", ylab = "Predator Intensity [0, 1]")
  
  if (length(selected_prey) > 0 && selected_prey[1] != "None (Unattenuated)") {
    lines(Years, res$P_eff_stage1, col = "hotpink3", lwd = 1.5, lty = 2)
    lines(Years, res$P_eff, col = "darkred", lwd = 2)
    legend("topleft", legend = c("Unattenuated (P_norm)", "+ Attenuated by Prey 1", "+ Joint Attenuated (Prey 1 & 2)"),
           col = c("gray60", "hotpink3", "darkred"), lty = c(1, 2, 1), lwd = c(1.5, 1.5, 2), bty = "n", cex = 0.8)
  } else {
    legend("topleft", legend = "Unattenuated Baseline", col = "gray60", lty = 1, lwd = 1.5, bty = "n", cex = 0.8)
  }
  
  # --- PANEL 4: STANDARDIZED SAR RESIDUALS VS STANDARDIZED P_EFF ---
  sar_std  <- scale(sar_resids)
  peff_std <- scale(res$P_eff)
  
  plot(Years, -sar_std, type = "h", lwd = 3, col = ifelse(sar_std >= 0, "dodgerblue3", "darkorange2"),
       ylim = c(min(sar_std, peff_std, na.rm=T) - 0.5, max(sar_std, peff_std, na.rm=T) + 0.5),
       main = "D. SAR Residuals vs P_eff (Std Z-Scores)", xlab = "Year", ylab = "Z-Score")
  abline(h = 0, lty = 2, col = "gray40")
  lines(Years, peff_std, col = "darkred", lwd = 2, type = "b", pch = 18)
  legend("topleft", legend = c("SAR Residual (CPUE-Adjusted)", "P_eff (Std Effective Pred)"),
         col = c("dodgerblue3", "darkred"), pch = c(NA, 18), lty = 1, lwd = c(3, 2), bty = "n", cex = 0.8)
  
  title(main = paste("Predator:", pred_col, "\nSelected Gate:", paste(selected_prey, collapse = " + "), 
                     "| AIC =", round(res$aic, 2), "| dAIC =", round(res$delta_aic, 2)), outer = TRUE, cex.main = 1.1)
}

# Run diagnostic plot on Steller Sea Lion Pups
plot_pure_minmax_diagnostic("X11_ssl_seak_pup_pred", minmax_complete_data)
plot_pure_minmax_diagnostic("X15_spinyDogfishGoA_predAK_smoltyr", minmax_complete_data)
plot_pure_minmax_diagnostic("X10_Harbour_s_2yrLead_WS_adultyr", minmax_complete_data)
plot_pure_minmax_diagnostic("X15_ArrowtoothFlounderBiomass_predAK_smoltyr", minmax_complete_data)
plot_pure_minmax_diagnostic("X15_DFA_sleeperSharks_adultyr", minmax_complete_data)
plot_pure_minmax_diagnostic("X10_AllSeaLionsEMB_2025_raw_sc_adultyr", minmax_complete_data)


cat("\n==============================================================================\n")
cat("      JOINT MULTI-PREY MODELS SORTED BY AIC (BEST GLOBAL FITS FIRST)           \n")
cat("==============================================================================\n\n")
print(joint_sweep_table %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, mean_net_impact, aic, delta_aic), n = 30)

joint_sweep_table <- joint_sweep_table %>% filter(predator!="X10_Harbor_seal_CR_2yrLead_adultyr")


joint_sweep_table %>% filter(beta_peff>0, pvalue<0.05) # all of these are unattenuated, so they are the same ones that were positive before 
      #and only 2 are significant (X10_Californian_s_l_2yrLead_WS_smoltyr which shouldn't have been there, and X10_Harbor_seal_CR_2yrLead_adultyr -- which also shouldn't have been there -- so NO PROBLEM NOW!!!!
joint_sweep_table %>% filter(beta_peff<0, pvalue<0.05, !grepl("smoltyr",predator))

# predator                                     selected_prey                                              fitted_weights beta_peff  pvalue std_eff_beta mean_net_impact realized_damping_pct   aic base_aic delta_aic
# 1 X15_sablefishBiomass_predAK_adultyr          X13_pollock_age1plus_adultyr + X12_egoa_krill_smoltyr      1; 1              -2.52  0.00116       -0.463          -0.140                 74.6  1.21     8.83     -7.63
# 2 X10_Harbour_s_2yrLead_WS_adultyr             None (Unattenuated)                                        None              -0.546 0.00685       -0.410          -0.181                  0    3.57     3.57      0   
# 3 X15_ArrowtoothFlounderBiomass_predAK_adultyr X13_mid_il_capelin_smoltyr                                 1                 -0.320 0.0134        -0.378          -0.136                 27.7  4.51     8.43     -3.92
# 4 X15_DFA_sleeperSharks_adultyr                None (Unattenuated)                                        None              -0.285 0.0199        -0.362          -0.106                  0    5.07     5.07      0   
# 5 X11_ssl_seak_pup_pred                        X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr 1; 1              -0.366 0.0250        -0.378          -0.122                 50.1  5.39     8.44     -3.05
# 1 X10_Californian_s_l_2yrLead_WS_smoltyr None (Unatte…                                                    None               0.485 0.00250        0.442           0.284                    0  2.21     2.21         0
 

# Define output destination
mypred="X15_sablefishBiomass_predAK_adultyr"
        plot_dir <- file.path(getwd(), "Rcode_for_paper", "Routput_for_paper", "plots")
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        out_file <- file.path(plot_dir, paste0("altprey_plots_",mypred,".png"))
        png(filename = out_file, width = 10, height = 8, units = "in", res = 300)
        plot_pure_minmax_diagnostic(mypred, minmax_complete_data)
        dev.off()
