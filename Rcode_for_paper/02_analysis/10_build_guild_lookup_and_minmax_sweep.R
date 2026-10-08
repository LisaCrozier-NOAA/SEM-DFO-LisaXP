# ==============================================================================
# Script: 10_build_guild_lookup_and_minmax_sweep.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Build 5-Guild Predator-Prey Lookup Table directly from column names,
#          incorporating smolt vs. adult outmigration year matching & spatial 
#          overlap indicators, followed by pure [0, 1] Min-Max joint sweep.
# Output: Rcode_for_paper/Routput_for_paper/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & RAW DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
plot_dir     <- file.path(master_out, "plots")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_all_data <-read.csv( file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_24yr_noNA <-read.csv( file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_short_noNA <- sem_data_all_data %>% 
  filter(complete.cases(.))

summary(sem_data_24yr_noNA);ncol(sem_data_24yr_noNA) #121: 1998-2021
summary(sem_data_short_noNA);ncol(sem_data_short_noNA) #141: 2003-2021
raw_sem_data  <- sem_data_24yr_noNA

# sem_orig      <- read.csv("data_Lisa/sem_master_data.csv", stringsAsFactors = FALSE)
# cap<-grep("capelin",names(raw_sem_data),value=T)
# raw_sem_data[,c("Year",cap)]
# plot(raw_sem_data[,c("Year")],raw_sem_data[,c(cap[c(4)])],type='l')
# > lines(raw_sem_data[,c("Year")],raw_sem_data[,c(cap[c(3)])],col=2)
# >
# cap<-grep("capelin",names(sem_orig),value=T);cap
# raw_sem_data[,c("Year",cap)]
# plot(sem_orig[,c("year")],sem_orig[,c(cap[c(1)])],type='l')
# 
# # # Ensure Sitka Herring columns exist for both smoltyr and adultyr
# raw_sem_data$X13_capelin_WGoA_smoltyr <- sem_orig$capelin_WGoA
# raw_sem_data$X13_capelin_WGoA_adultyr <- sem_orig$capelin_WGoA[c(3:24, 23, 24)]
# # raw_sem_data$X13_sitkaHerring_EGoA_smoltyr <- sem_orig$sitkaHerring_EGoA
# # raw_sem_data$X13_sitkaHerring_EGoA_adultyr <- sem_orig$sitkaHerring_EGoA[c(3:24, 24, 24)]

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"

all_cols <- names(raw_sem_data)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE FUNCTIONAL PREY POOLS & SPATIAL INDICATORS BY GUILD
# ------------------------------------------------------------------------------
get_prey_pool <- function(guild_num, stage_suffix) {
  sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
  
  # Base Candidate Terms by Region & Guild
  m_squid     <- paste0("X04_marketsquid_GAM", sfx)
  anchovy     <- paste0("X05_anchovy_GAM", sfx)
  sardine     <- paste0("X05_DFA_abundSardine", sfx)
  herring_ncc <- paste0("X05_herring_GAM", sfx)
  eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
  hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
  
  pollock     <- paste0("X13_pollock_age1plus", sfx)
  herring_ak  <- paste0("X13_sitkaHerring_EGoA", sfx)
  capelin     <- paste0("X13_mid_il_capelin", sfx)
  krill_ak    <- paste0("X12_egoa_krill", sfx)
  pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
  
  prey_list <- switch(
    as.character(guild_num),
    "1" = c(m_squid, anchovy, sardine, herring_ncc, eulachon),         # NCC Small Avian
    "2" = c(hake, m_squid, anchovy, sardine, herring_ncc, eulachon),    # NCC Mammals & Large Fish
    "3" = c(pollock, herring_ak, capelin, krill_ak),                    # AK Krill-Eaters
    "4" = c(pollock, herring_ak, capelin),                              # AK Demersal Piscivores (No Krill)
    "5" = c(pink_ak, pollock, herring_ak, capelin),                     # AK Apex Sharks
    c()
  )
  
  # Keep only prey columns that exist in the raw dataset
  intersect(prey_list, all_cols)
}

# Assign designated spatial overlap indicator per predator
get_spatial_overlap_col <- function(pred_col, region, stage_suffix) {
  sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
  
  if (grepl("X09_DFA_HakeAge5Plus", pred_col)) {
    col <- paste0("X12_DFA_biomassEuphShelfSum", sfx)
  } else if (region == "NCC") {
    col <- paste0("X01_habCompInd", sfx)
  } else {
    col <- paste0("X21_sst_egoa_junjulaug", sfx)
  }
  
  if (col %in% all_cols) return(col) else return(NA_character_)
}

# ------------------------------------------------------------------------------
# STEP 2: DYNAMICALLY BUILD 5-GUILD LOOKUP TABLE
# ------------------------------------------------------------------------------
predator_cols <- all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)]

# Blacklist unwanted species / series
blacklist <- c("canaryRockfish", "Harbor_seal_CR")
predator_cols <- predator_cols[!sapply(predator_cols, function(x) any(sapply(blacklist, grepl, x)))]

lookup_rows <- list()

for (pred in predator_cols) {
  region <- if (grepl("^(X08_|X09_|X10_)", pred) || pred == "X11_DFA_Harbour_p_WS") "NCC" else "AK"
  
  # Explicitly preserve smoltyr vs adultyr matching (X10 naturally maps to either)
  stage_suffix <- if (grepl("smoltyr$", pred)) {
    "smoltyr"
  } else if (grepl("adultyr$", pred)) {
    "adultyr"
  } else {
    "smoltyr" # Default fallback for unsuffixed series
  }
  
  # Assign Functional Guild (1 through 5)
  guild_num <- case_when(
    region == "NCC" & grepl("^(X08_|X11_DFA_Harbour_p_WS)", pred) ~ 1, # NCC Small Avian / Porpoise
    region == "NCC" ~ 2,                                              # NCC Mammals & Large Fish (X09, X10)
    region == "AK" & grepl("(ssl|sablefish|spinyDogfish)", pred, ignore.case = TRUE) ~ 3, # AK Krill-Eaters
    region == "AK" & grepl("(halibut|PacificCod|Arrowtooth)", pred, ignore.case = TRUE) ~ 4, # AK Demersal Piscivores
    region == "AK" & grepl("shark", pred, ignore.case = TRUE) ~ 5,                          # AK Apex Sharks
    TRUE ~ 3
  )
  
  prey_candidates <- get_prey_pool(guild_num, stage_suffix)
  spatial_col     <- get_spatial_overlap_col(pred, region, stage_suffix)
  
  lookup_rows[[length(lookup_rows) + 1]] <- tibble(
    region              = region,
    guild_num           = guild_num,
    pred_data_col       = pred,
    stage_suffix        = stage_suffix,
    spatial_overlap_col = spatial_col,
    altprey1_data_col   = ifelse(length(prey_candidates) >= 1, prey_candidates[1], NA_character_),
    altprey2_data_col   = ifelse(length(prey_candidates) >= 2, prey_candidates[2], NA_character_),
    altprey3_data_col   = ifelse(length(prey_candidates) >= 3, prey_candidates[3], NA_character_),
    altprey4_data_col   = ifelse(length(prey_candidates) >= 4, prey_candidates[4], NA_character_),
    altprey5_data_col   = ifelse(length(prey_candidates) >= 5, prey_candidates[5], NA_character_),
    all_candidates_str  = paste(c(spatial_col, prey_candidates)[!is.na(c(spatial_col, prey_candidates))], collapse = ", ")
  )
}

guild_pred_lookup <- bind_rows(lookup_rows)

cat("\n==============================================================================\n")
cat("          5-GUILD PREDATOR-PREY LOOKUP TABLE GENERATED SUCCESSFULLY            \n")
cat("==============================================================================\n\n")
print(guild_pred_lookup %>% select(region, guild_num, pred_data_col, stage_suffix, spatial_overlap_col, all_candidates_str), n = 50)

write_csv(guild_pred_lookup, file.path(data_out_dir, "guild_predator_lookup_table.csv"))

# ------------------------------------------------------------------------------
# STEP 3: MIN-MAX SCALING [0, 1] & JOINT ATTENUATION SWEEP
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Scale all numeric indicators EXCEPT targets and Year into [0, 1]
minmax_complete_data <- raw_sem_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

find_best_joint_prey_minmax <- function(pred_col, df_model, max_combo_size = 2) {
  P_norm <- df_model[[pred_col]]
  row_data <- guild_pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  # Extract candidates (Spatial Overlap + Guild Prey Candidates)
  candidates <- unlist(strsplit(row_data$all_candidates_str, ", "))
  candidates <- intersect(candidates[!is.na(candidates) & candidates != ""], names(df_model))
  
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
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_sar]])
  
  
  best_res <- list(
    predator = pred_col, selected_prey = "None (Unattenuated)",
    weights = c(), aic = base_aic, base_aic = base_aic, delta_aic = 0.0,
    beta = base_est, pval = base_p, P_eff = P_norm, P_eff_stage1 = P_norm,
    sar_r2  <- base_r2,
    sar_r2_base  <- base_r2
    
  )
  
  if (length(candidates) == 0) return(best_res)
  
  Anorm_mat <- as.matrix(df_model[, candidates, drop = FALSE])
  w_grid    <- seq(0.1, 1.0, by = 0.1)
  
  # Evaluate 1-prey and 2-prey joint gate combinations
  for (k in 1:min(length(candidates), max_combo_size)) {
    combos <- combn(candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      w_combo_grid <- expand.grid(rep(list(w_grid), k))
      
      for (g_idx in 1:nrow(w_combo_grid)) {
        w_vals <- as.numeric(w_combo_grid[g_idx, ])

        # Gate 1
        w1_gate  <- 1 - w_vals[1] * Anorm_mat[, cb[1]]
        P_joint  <- P_norm * w1_gate
        
        # Gate 2 (if present)
        if (k >= 2) {
          w2_gate <- 1 - w_vals[2] * Anorm_mat[, cb[2]]
          P_joint <- P_joint * w2_gate
        }
        
        # Gate 3 (if present)
        if (k >= 3) {
          w3_gate <- 1 - w_vals[3] * Anorm_mat[, cb[3]]
          P_joint <- P_joint * w3_gate
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
                P_eff = P_joint, P_eff_stage1 =  P_norm * w1_gate,# P_eff_stage2 =   P_norm * w1_gate * w2_gate, P_eff_stage3 =   P_norm * w1_gate * w2_gate * w3_gate
                sar_r2  <- as.numeric(inspect(fit, "r2")[[target_sar]]),
                sar_r2_base  <- base_r2
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
# STEP 4: EXECUTE ATTENUATION SWEEP & SAVE OUTPUTS
# ------------------------------------------------------------------------------
cat("Executing Joint Attenuation Sweep across All 5-Guild Predators...\n")

# ------------------------------------------------------------------------------
# STEP 4 (UPDATED): EXECUTE SWEEP WITH PROGRESS LOGGING & CHECKPOINTS
# ------------------------------------------------------------------------------
cat("Executing Joint Attenuation Sweep across All 5-Guild Predators...\n")

sweep_predators     <- guild_pred_lookup$pred_data_col
total_preds         <- length(sweep_predators)
joint_sweep_results <- list()

checkpoint_file <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_CHECKPOINT.csv")

for (i in seq_along(sweep_predators[1])) {
  print(Sys.time())
  pred_col <- sweep_predators[i]
  
  # Print real-time progress counter
  cat(sprintf("[%d/%d] Testing predator: %s ...\n", i, total_preds, pred_col))
  
  if (pred_col %in% names(minmax_complete_data)) {
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
        delta_aic            = res$delta_aic,
        base_r2             = res$sar_r2_base ,
        final_r2            = res$sar_r2
      )
    }
  }
  
  # Save checkpoint every 5 predators (or on the final iteration)
  if (i %% 5 == 0 || i == total_preds) {
    if (length(joint_sweep_results) > 0) {
      chk_table <- map_dfr(joint_sweep_results, ~ .x) %>% arrange(aic)
      write_csv(chk_table, checkpoint_file)
      cat(sprintf("  --> Checkpoint saved to disc (%d models completed).\n", nrow(chk_table)))
    }
  }
}
print(Sys.time())
print(joint_sweep_results)

joint_sweep_table <- map_dfr(joint_sweep_results, ~ .x) %>% arrange(aic)
#write_csv(joint_sweep_table, file.path(data_out_dir, "Goal3_5predGuild_Joint_Sweep_Results_24yr_noNA.csv"))

cat("\n==============================================================================\n")
cat("       5-GUILD ATTENUATION SWEEP RESULTS SORTED BY AIC (BEST FITS TOP)        \n")
cat("==============================================================================\n\n")
print(joint_sweep_table %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, mean_net_impact, aic, delta_aic), n = 50)

write_csv(joint_sweep_table, file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA.csv"))

# ------------------------------------------------------------------------------
# STEP 5: 4-PANEL DIAGNOSTIC PLOTTING FUNCTION
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
  
  # PANEL 1: RAW VS EFFECTIVE PREDATOR
  plot(Years, P_norm, type = "b", pch = 16, col = "black", ylim = c(0, 1.1),
       main = "A. Raw vs. Joint Effective Predator Index", xlab = "Year", ylab = "Predator Scale [0, 1]")
  lines(Years, res$P_eff, type = "b", pch = 17, col = "firebrick", lty = 2, lwd = 2)
  legend("topleft", legend = c("Raw Predator (P_norm)", "Joint Effective (P_eff)"),
         col = c("black", "firebrick"), pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.8)
  
  # PANEL 2: PREY ABUNDANCE & WEIGHT GATE
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
    legend("topright", legend = leg_txt, col = leg_col, lty = leg_lty, lwd = c(1, 2, 1, 2)[1:length(leg_txt)], bty = "n", cex = 0.7)
  } else {
    plot(Years, rep(1, length(Years)), type = "l", col = "gray", ylim = c(0, 1.1),
         main = "B. No Alternate Prey Selected", xlab = "Year", ylab = "Weight")
  }
  
  # PANEL 3: STEPWISE ATTENUATION PROGRESSION
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
  
  # PANEL 4: SAR RESIDUALS VS EFFECTIVE PREDATOR Z-SCORES
  sar_std  <- scale(sar_resids)
  peff_std <- scale(res$P_eff)
  
  plot(Years, -sar_std, type = "h", lwd = 3, col = ifelse(sar_std >= 0, "dodgerblue3", "darkorange2"),
       ylim = c(min(sar_std, peff_std, na.rm = TRUE) - 0.5, max(sar_std, peff_std, na.rm = TRUE) + 0.5),
       main = "D. SAR Residuals vs P_eff (Std Z-Scores)", xlab = "Year", ylab = "Z-Score")
  abline(h = 0, lty = 2, col = "gray40")
  lines(Years, peff_std, col = "darkred", lwd = 2, type = "b", pch = 18)
  legend("topleft", legend = c("SAR Residual (CPUE-Adjusted)", "P_eff (Std Effective Pred)"),
         col = c("dodgerblue3", "darkred"), pch = c(NA, 18), lty = 1, lwd = c(3, 2), bty = "n", cex = 0.8)
  
  title(main = paste("Predator:", pred_col, "\nSelected Gate:", paste(selected_prey, collapse = " + "), 
                     "| AIC =", round(res$aic, 2), "| dAIC =", round(res$delta_aic, 2)), outer = TRUE, cex.main = 1.1)
}

# Example diagnostic plot generation for Sablefish
mypred   <- "X15_sablefishBiomass_predAK_adultyr"
plot_pure_minmax_diagnostic(mypred, minmax_complete_data)


out_file <- file.path(plot_dir, paste0("altprey_plots_", mypred, ".png"))
png(filename = out_file, width = 10, height = 8, units = "in", res = 300)
plot_pure_minmax_diagnostic(mypred, minmax_complete_data)
dev.off()

#RESULTS-----------
myresults<-(joint_sweep_table %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, mean_net_impact, aic, delta_aic), n = 50)

joint_sweep_table %>% filter(pvalue>0.05)
joint_sweep_table %>% filter(pvalue>0.05) %>% select(predator,selected_prey,fitted_weights,beta_peff,aic)
joint_sweep_table %>% filter(beta_peff>0) %>% select(predator,selected_prey,fitted_weights,beta_peff,aic)


#all pred w/ positive peff failed to select any altprey/sst
joint_sweep_table %>% filter(beta_peff>0) %>% select(predator,selected_prey,fitted_weights,beta_peff,aic)
# A tibble: 12 × 5
# predator                                   selected_prey       fitted_weights beta_peff   aic
# 1 X09_DFA_ChinAbundSnakeFall_adultyr         None (Unattenuated) None               0.678 -8.95
# 2 X10_Californian_s_l_2yrLead_WS_smoltyr     None (Unattenuated) None               0.588  1.93
# 3 X11_ssl_seak_pup_pred                      None (Unattenuated) None               0.443  5.37
# 4 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr None (Unattenuated) None               0.388  7.35
# 5 X15_PacificCodBiomass_predAK_smoltyr       None (Unattenuated) None               0.347  7.62
# 6 X10_Californian_s_l_2yrLead_WS_adultyr     None (Unattenuated) None               0.333  8.69
# 7 X09_DFA_ChinAbundSnakeFall_smoltyr         None (Unattenuated) None               0.278  8.81
# 8 X10_AllSeaLionsEMB_smoltyr                 None (Unattenuated) None               0.285  9.03
# 9 X10_AllSeaLionsEMB_adultyr                 None (Unattenuated) None               0.285  9.04
# 10 X10_Northern_f_s_2yrLead_WS_adultyr        None (Unattenuated) None               0.266  9.77
# 11 X10_DFA_ssl.est.wholerange_2yrLead_adultyr None (Unattenuated) None               0.174 10.3 
# 12 X15_spinyDogfishGoA_predAK_adultyr         None (Unattenuated) None               0.162 10.5 


print(joint_sweep_table %>% filter(beta_peff<0) %>% select(predator,selected_prey,fitted_weights,beta_peff,aic),n=Inf)
print(joint_sweep_table %>% filter(beta_peff<0) %>% select(predator,selected_prey,fitted_weights,beta_peff,aic,delta_aic),n=Inf)
# A tibble: 35 × 6
# predator                                          selected_prey                                                     fitted_weights beta_peff     aic delta_aic
# 1 X15_salmonSharkBSAI_predAK_smoltyr                X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         0.9; 1            -2.36  -9.87     -19.5  
# 2 X15_sablefishBiomass_predAK_smoltyr               X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -4.56  -7.10     -18.0  
# 3 X09_DFA_HakeAge5Plus_smoltyr                      X04_marketsquid_GAM_smoltyr + X05_anchovy_GAM_smoltyr             0.2; 0.4          -1.41  -5.37      -6.28 
# 4 X09_DFA_HakeAge5Plus_smoltyr                      X04_marketsquid_GAM_smoltyr + X05_anchovy_GAM_smoltyr             0.2; 0.4          -1.41  -5.37      -6.28 
# 5 X08_DFA_DC_corm_3_WS                              X01_habCompInd_smoltyr + X05_herring_GAM_smoltyr                  1; 0.8            -1.66  -4.46     -15.3  
# 6 X15_sharkCatchGoA_predAK_smoltyr                  X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -2.69  -4.02     -14.5  
# 7 X10_Harbour_s_2yrLead_WS_smoltyr                  X01_habCompInd_smoltyr + X05_herring_GAM_smoltyr                  1; 0.3            -1.12  -1.75     -10.1  
# 8 X11_DFA_Harbour_p_WS                              X01_habCompInd_smoltyr + X05_DFA_abundSardine_smoltyr             1; 0.5            -1.66  -1.00      -9.51 
# 9 X15_sablefishBiomass_predAK_adultyr               X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -2.27  -0.771    -11.4  
# 10 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr X21_sst_egoa_junjulaug_adultyr + X13_mid_il_capelin_adultyr       0.4; 1            -0.927 -0.712     -9.11 
# 11 X10_Harbour_s_2yrLead_WS_adultyr                  X01_habCompInd_adultyr + X05_eulachon_during_chinook_adultyr      0.6; 1            -0.874 -0.483     -6.36 
# 12 X09_chilipepper_adultyr                           X01_habCompInd_adultyr + X04_marketsquid_GAM_adultyr              0.9; 1            -4.43  -0.391     -9.23 
# 13 X15_DFA_sleeperSharks_adultyr                     X14_pinkSalmonNorthAmerica_adultyr + X13_mid_il_capelin_adultyr   0.4; 1            -0.715 -0.352     -4.79 
# 14 X15_ArrowtoothFlounderBiomass_predAK_adultyr      X21_sst_egoa_junjulaug_adultyr + X13_mid_il_capelin_adultyr       0.2; 1            -0.694  0.0275   -10.3  
# 15 X15_DFA_sleeperSharks_smoltyr                     X13_mid_il_capelin_smoltyr                                        1                 -0.623  0.918     -0.122
# 16 X09_DFA_HakeAge5Plus_adultyr                      X12_DFA_biomassEuphShelfSum_adultyr + X05_anchovy_GAM_adultyr     0.9; 0.8          -1.67   1.58      -7.66 
# 17 X15_salmonSharkGoA_predAK_smoltyr                 X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -2.28   2.59      -8.36 
# 18 X15_salmonSharkBSAI_predAK_adultyr                X14_pinkSalmonNorthAmerica_adultyr + X13_pollock_age1plus_adultyr 1; 1              -1.93   3.21      -6.04 
# 19 X15_sablefishRecruitment_predAK_adultyr           X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -4.00   3.30      -7.65 
# 20 X15_spinyDogfishBSAI_predAK_smoltyr               X13_mid_il_capelin_smoltyr + X12_egoa_krill_smoltyr               1; 1              -0.933  3.55      -5.71 
# 21 X15_PacificCodBiomass_predAK_adultyr              X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -0.794  4.37      -6.50 
# 22 X08_Loons_8_WS                                    X01_habCompInd_smoltyr + X05_herring_GAM_smoltyr                  1; 0.6            -1.33   4.66      -6.32 
# 23 X15_salmonSharkGoA_predAK_adultyr                 X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -1.76   5.11      -4.82 
# 24 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -0.465  5.14      -4.69 
# 25 X15_ArrowtoothFlounderBiomass_predAK_smoltyr      X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -0.402  5.46      -3.92 
# 26 X15_spinyDogfishBSAI_predAK_adultyr               X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -0.868  6.25      -3.96 
# 27 X10_Northern_f_s_2yrLead_WS_smoltyr               X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr              1; 0.5            -0.701  6.56      -4.24 
# 28 X15_spinyDogfishGoA_predAK_smoltyr                X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr         1; 1              -1.71   6.69      -4.30 
# 29 X08_Large_gulls_7_WS                              X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr              1; 0.5            -0.984  7.12      -3.40 
# 30 X10_AllSeaLionsEMB_2025_raw_sc_smoltyr            X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr              1; 1              -7.50   7.75      -2.43 
# 31 X15_sablefishRecruitment_predAK_smoltyr           X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr      1; 1              -5.01   8.89      -1.83 
# 32 X15_sharkCatchGoA_predAK_adultyr                  X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr         1; 1              -1.24   9.42      -1.53 
# 33 X08_commonMurre_JSOES                             X01_habCompInd_smoltyr + X05_DFA_abundSardine_smoltyr             1; 0.4            -1.89   9.54      -0.306
# 34 X09_chilipepper_smoltyr                           X01_habCompInd_smoltyr + X05_anchovy_GAM_smoltyr                  1; 1              -0.598 10.4       -0.584
# 35 X10_AllSeaLionsEMB_2025_raw_sc_adultyr            X01_habCompInd_adultyr + X05_eulachon_during_chinook_adultyr      1; 1              -0.116 10.8       -0.218
# > 


print(joint_sweep_table %>% filter(fitted_weights=="None") %>% select(predator,selected_prey,fitted_weights,beta_peff,aic,delta_aic),n=Inf)
# predator                                   selected_prey       fitted_weights beta_peff   aic delta_aic
# 1 X09_DFA_ChinAbundSnakeFall_adultyr         None (Unattenuated) None               0.678 -8.95         0
# 2 X10_Californian_s_l_2yrLead_WS_smoltyr     None (Unattenuated) None               0.588  1.93         0
# 3 X11_ssl_seak_pup_pred                      None (Unattenuated) None               0.443  5.37         0
# 4 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr None (Unattenuated) None               0.388  7.35         0
# 5 X15_PacificCodBiomass_predAK_smoltyr       None (Unattenuated) None               0.347  7.62         0
# 6 X10_Californian_s_l_2yrLead_WS_adultyr     None (Unattenuated) None               0.333  8.69         0
# 7 X09_DFA_ChinAbundSnakeFall_smoltyr         None (Unattenuated) None               0.278  8.81         0
# 8 X10_AllSeaLionsEMB_smoltyr                 None (Unattenuated) None               0.285  9.03         0
# 9 X10_AllSeaLionsEMB_adultyr                 None (Unattenuated) None               0.285  9.04         0
# 10 X10_Northern_f_s_2yrLead_WS_adultyr        None (Unattenuated) None               0.266  9.77         0
# 11 X10_DFA_ssl.est.wholerange_2yrLead_adultyr None (Unattenuated) None               0.174 10.3          0
# 12 X15_spinyDogfishGoA_predAK_adultyr         None (Unattenuated) None               0.162 10.5          0
# > 