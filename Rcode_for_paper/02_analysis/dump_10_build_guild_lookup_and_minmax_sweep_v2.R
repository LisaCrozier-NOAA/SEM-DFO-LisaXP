#output: Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA.csv

# ==============================================================================
# Script: 10_build_guild_lookup_and_minmax_sweep_v2.R
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

sem_data_all_data  <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_24yr_noNA <- read.csv(file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_short_noNA <- sem_data_all_data %>% 
  filter(complete.cases(.))

# Choose 19 vs 24-year complete dataset
#raw_sem_data <- sem_data_24yr_noNA
raw_sem_data <- sem_data_short_noNA;nrow(raw_sem_data);ncol(raw_sem_data)

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"

all_cols <- names(raw_sem_data)

# ------------------------------------------------------------------------------
# MIN-MAX SCALING [0, 1]
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_complete_data <- raw_sem_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))


# ------------------------------------------------------------------------------
# STEP 1: DEFINE FUNCTIONAL PREY POOLS & SPATIAL INDICATORS BY GUILD
# ------------------------------------------------------------------------------
get_prey_pool <- function(guild_num, stage_suffix) {
  sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
  
  m_squid     <- paste0("X04_marketsquid_GAM", sfx)
  anchovy     <- paste0("X05_anchovy_GAM", sfx)
  sardine     <- paste0("X05_DFA_abundSardine", sfx)
  herring_ncc <- paste0("X05_herring_GAM", sfx)
  eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
  hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
  
  pollock     <- paste0("X13_pollock_age1plus", sfx)
  herring_ak  <- paste0("X13_sitkaHerring_EGoA", sfx)
  capelin     <- paste0("X13_capelin_WGoA", sfx)
#  capelin     <- paste0("X13_mid_il_capelin", sfx)
  krill_ak    <- paste0("X12_egoa_krill", sfx)
  pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
  
  prey_list <- switch(
    as.character(guild_num),
    "1" = c(m_squid, anchovy, sardine, herring_ncc, eulachon),
    "2" = c(hake, m_squid, anchovy, sardine, herring_ncc, eulachon),
    "3" = c(pollock, herring_ak, capelin, krill_ak),
    "4" = c(pollock, herring_ak, capelin),
    "5" = c(pink_ak, pollock, herring_ak, capelin),
    c()
  )
  
  intersect(prey_list, all_cols)
}

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

blacklist <- c("canaryRockfish","chilipepper", "Harbor_seal_CR",
                "X09_DFA_ChinAbundSnakeFall_adultyr","X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA","X15_spinyDogfishGoA_predAK_adultyr")
predator_cols <- predator_cols[!sapply(predator_cols, function(x) any(sapply(blacklist, grepl, x)))]

lookup_rows <- list()

for (pred in predator_cols) {
  region <- if (grepl("^(X08_|X09_|X10_)", pred) || pred == "X11_DFA_Harbour_p_WS") "NCC" else "AK"
  
  stage_suffix <- if (grepl("smoltyr$", pred)) {
    "smoltyr"
  } else if (grepl("adultyr$", pred)) {
    "adultyr"
  } else {
    "smoltyr"
  }
  
  guild_num <- case_when(
    region == "NCC" & grepl("^(X08_|X11_DFA_Harbour_p_WS)", pred) ~ 1,
    region == "NCC" ~ 2,
    region == "AK" & grepl("(ssl|sablefish|spinyDogfish)", pred, ignore.case = TRUE) ~ 3,
    region == "AK" & grepl("(halibut|PacificCod|Arrowtooth)", pred, ignore.case = TRUE) ~ 4,
    region == "AK" & grepl("shark", pred, ignore.case = TRUE) ~ 5,
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
write_csv(guild_pred_lookup, file.path(data_out_dir, "guild_predator_lookup_table.csv"))

# ------------------------------------------------------------------------------
# STEP 3: JOINT ATTENUATION SWEEP-----
# ------------------------------------------------------------------------------

find_best_joint_prey_minmax <- function(pred_col, df_model, max_combo_size = 2) {
  P_norm <- df_model[[pred_col]]
  row_data <- guild_pred_lookup %>% filter(pred_data_col == pred_col)
  if (nrow(row_data) == 0) return(NULL)
  row_data <- row_data %>% slice(1)
  
  candidates <- unlist(strsplit(row_data$all_candidates_str, ", "))
  candidates <- intersect(candidates[!is.na(candidates) & candidates != ""], names(df_model))
  
  # Base Model
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
  
  # FIX: Correct key = value syntax in list definition
  best_res <- list(
    predator      = pred_col, 
    selected_prey = "None (Unattenuated)",
    weights       = c(), 
    aic           = base_aic, 
    base_aic      = base_aic, 
    delta_aic     = 0.0,
    beta          = base_est, 
    base_beta     = base_est,  # <--- PLACE 1: Store base_est here
    pval          = base_p, 
    P_eff         = P_norm, 
    P_eff_stage1  = P_norm,
    sar_r2        = base_r2,
    sar_r2_base   = base_r2
  )
  
  if (length(candidates) == 0) return(best_res)
  
  Anorm_mat <- as.matrix(df_model[, candidates, drop = FALSE])
  w_grid    <- seq(0.1, 1.0, by = 0.1)
  
  for (k in 1:min(length(candidates), max_combo_size)) {
    combos <- combn(candidates, k, simplify = FALSE)
    
    for (cb in combos) {
      w_combo_grid <- expand.grid(rep(list(w_grid), k))
      
      for (g_idx in 1:nrow(w_combo_grid)) {
        w_vals <- as.numeric(w_combo_grid[g_idx, ])
        
        w1_gate <- 1 - w_vals[1] * Anorm_mat[, cb[1]]
        P_joint <- P_norm * w1_gate
        
        if (k >= 2) {
          w2_gate <- 1 - w_vals[2] * Anorm_mat[, cb[2]]
          P_joint <- P_joint * w2_gate
        }
        
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
            # # --- CALCULATE PENALIZED AIC BEFORE MODEL COMPARISON ---
            # n_tuned_weights <- length(w_vals)
            # curr_aic_raw    <- as.numeric(fitMeasures(fit, "aic"))
            # curr_aic        <- curr_aic_raw + (2 * n_tuned_weights)
            
            if (curr_aic < best_res$aic) {
              
              # FIX: Correct key = value syntax in list definition
              best_res <- list(
                predator      = pred_col, 
                selected_prey = cb,
                weights       = w_vals, 
                aic           = curr_aic, 
                base_aic      = base_aic,
                delta_aic     = curr_aic - base_aic,
                beta          = as.numeric(pe$est[1]), 
                base_beta     = base_est,  # <--- PLACE 2: Keep base_est here
                pval          = as.numeric(pe$pvalue[1]),
                P_eff         = P_joint, 
                P_eff_stage1  = P_norm * w1_gate,
                sar_r2        = as.numeric(inspect(fit, "r2")[[target_sar]]),
                sar_r2_base   = base_r2
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
# STEP 4: EXECUTE ATTENUATION SWEEP ACROSS ALL PREDATORS
# ------------------------------------------------------------------------------
cat("Executing Joint Attenuation Sweep across All 5-Guild Predators...\n")

sweep_predators     <- guild_pred_lookup$pred_data_col
total_preds         <- length(sweep_predators)
joint_sweep_results <- list()

checkpoint_file <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_CHECKPOINT_19yr_2prey.csv")

# FIX: Loop over ALL predators in sweep_predators
for (i in seq_along(sweep_predators)) {
    print(Sys.time())
    pred_col <- sweep_predators[i]
  cat(sprintf("[%d/%d] Testing predator: %s ...\n", i, total_preds, pred_col))
  
  if (pred_col %in% names(minmax_complete_data)) {
    res <- find_best_joint_prey_minmax(pred_col=pred_col, df_model=minmax_complete_data, max_combo_size = 2)
    
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
        base_beta            = res$base_beta,  # <--- PLACE 3: Add column to tibble
        beta_peff            = res$beta,
        pvalue               = res$pval,
        std_eff_beta         = std_beta,
        mean_net_impact      = mean_net,
        realized_damping_pct = damp_pct,
        aic                  = res$aic,
        base_aic             = res$base_aic,
        delta_aic            = res$delta_aic,
        base_r2              = res$sar_r2_base,
        final_r2             = res$sar_r2
      )
    }
  }
  
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
print(Sys.time())

cat("\n==============================================================================\n")
cat("       5-GUILD ATTENUATION SWEEP RESULTS SORTED BY AIC (BEST FITS TOP)        \n")
cat("==============================================================================\n\n")
print(
  joint_sweep_table %>% 
    select(predator, selected_prey, fitted_weights,base_beta, beta_peff, pvalue, base_r2, final_r2, aic, delta_aic), 
  n = 50
)

#write_csv(joint_sweep_table, file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA_akferris.csv"))
write_csv(joint_sweep_table, file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_19yr_noNA_dougsmoothed.csv"))
