# ==============================================================================
# Script: 04_OceanPred_plus_EMBssl.r
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Test dual-arena predation on SAR:
#          Equation 1: CPUE ~ Hake_smoltyr
#          Equation 2: SAR  ~ CPUE + P_ocean_eff + P_inriver_adult
#          P_ocean_eff uses best candidate (a=0.2, w=1.0 for both SST & AltPrey)
#          P_inriver_adult uses raw baseline X10_DFA_AllSeaLionsEMB_adultyr
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper", "data")

raw_csv_path  <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Complete cases for 24-year series
valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) == 24]
df_raw <- raw_df_master %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# 0-1 Min-Max scaled base
df_base <- df_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue     <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar      <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col        <- "X09_DFA_HakeAge5Plus_smoltyr"
inriver_pred_col <- "X10_DFA_AllSeaLionsEMB_adultyr"

# Check in-river predictor presence
if (!inriver_pred_col %in% names(df_base)) {
  stop("In-river predator variable not found in 24-year complete dataset: ", inriver_pred_col)
}

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TOP OCEAN CANDIDATES TO TEST DUAL-ARENA HYPOTHESIS
# ------------------------------------------------------------------------------
ocean_candidates <- tribble(
  ~ocean_pred,                               ~spatial_col,                      ~prey_vars,
  "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", "X21_sst_egoa_junjulaug_smoltyr",  "capelin_WGoA_smoltyr",
  "X15_salmonSharkGoA_predAK_smoltyr",    "X21_sst_egoa_junjulaug_smoltyr",  "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr",
  "X15_sablefishAK_biomass3plus_smoltyr",    "X21_sst_egoa_junjulaug_smoltyr",  "X13_pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr",
  "X11_ssl_seak_pup_pred_smoltyr",           "X21_sst_egoa_junjulaug_smoltyr",  "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr",
  "X15_DFA_halibutBiomassAge8plus_predAK_adultyr", "X21_sst_egoa_junjulaug_adultyr", "capelin_WGoA_adultyr"
)

# ------------------------------------------------------------------------------
# STEP 2: FIT SINGLE vs. DUAL-ARENA MODELS
# ------------------------------------------------------------------------------
results_list <- list()

for (i in seq_len(nrow(ocean_candidates))) {
  cand <- ocean_candidates[i, ]
  
  # --- Construct P_ocean_eff with a = 0.2 and w = 1.0 ---
  P_norm <- df_base[[cand$ocean_pred]]
  
  # SST Gate (a = 0.2, w = 1.0)
  raw_temp  <- df_raw[[cand$spatial_col]]
  temp_gate <- scale_custom_bounds(raw_temp, a = 0.2, b = 1.0)
  P_temp    <- P_norm * temp_gate
  
  # Alternative Prey Gate (a = 0.2, w = 1.0)
  prey_cols <- trimws(unlist(strsplit(cand$prey_vars, "\\+")))
  if (length(prey_cols) == 1) {
    raw_forage <- df_base[[prey_cols[1]]]
  } else {
    raw_forage <- rowSums(df_base[, prey_cols, drop = FALSE], na.rm = TRUE)
  }
  forage_gate <- 1 - scale_custom_bounds(raw_forage, a = 0.2, b = 1.0)
  
  # Final Effective Ocean Predation
  df_base$P_ocean_eff <- P_temp * forage_gate
  
  # --- Model A: Single-Arena (Ocean Only) ---
  sem_single <- paste0(
    target_cpue, " ~ ", hake_col, "\n",
    target_sar,  " ~ ", target_cpue, " + P_ocean_eff"
  )
  fit_single <- sem(sem_single, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)
  
  aic_single <- as.numeric(fitMeasures(fit_single, "aic"))
  r2_single  <- as.numeric(inspect(fit_single, "r2")[[target_sar]])
  pe_single  <- parameterEstimates(fit_single, standardized = TRUE) %>%
    filter(op == "~" & lhs == target_sar & rhs == "P_ocean_eff")
  
  # --- Model B: Dual-Arena (Ocean Eff + In-River Adult Sea Lions) ---
  sem_dual <- paste0(
    target_cpue, " ~ ", hake_col, "\n",
    target_sar,  " ~ ", target_cpue, " + P_ocean_eff + ", inriver_pred_col
  )
  fit_dual <- sem(sem_dual, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)
  
  if (!is.null(fit_dual) && lavInspect(fit_dual, "converged")) {
    aic_dual <- as.numeric(fitMeasures(fit_dual, "aic"))
    r2_dual  <- as.numeric(inspect(fit_dual, "r2")[[target_sar]])
    pe_dual  <- parameterEstimates(fit_dual, standardized = TRUE) %>%
      filter(op == "~" & lhs == target_sar)
    
    pe_ocean   <- pe_dual %>% filter(rhs == "P_ocean_eff")
    pe_inriver <- pe_dual %>% filter(rhs == inriver_pred_col)
    
    results_list[[i]] <- tibble(
      ocean_predator    = cand$ocean_pred,
      alt_prey          = cand$prey_vars,
      
      # Single-Arena (Ocean Only)
      single_aic        = round(aic_single, 2),
      single_r2         = round(r2_single, 3),
      single_beta_ocean = round(pe_single$est[1], 3),
      single_std_ocean  = round(pe_single$std.all[1], 3),
      single_p_ocean    = round(pe_single$pvalue[1], 4),
      
      # Dual-Arena (Ocean + In-River)
      dual_aic          = round(aic_dual, 2),
      delta_aic         = round(aic_dual - aic_single, 2),
      dual_r2           = round(r2_dual, 3),
      r2_gain           = round(r2_dual - r2_single, 3),
      
      # Coefficients in Dual Model
      dual_std_ocean    = round(pe_ocean$std.all[1], 3),
      dual_p_ocean      = round(pe_ocean$pvalue[1], 4),
      dual_std_inriver  = round(pe_inriver$std.all[1], 3),
      dual_p_inriver    = round(pe_inriver$pvalue[1], 4)
    )
  }
}

dual_arena_table <- bind_rows(results_list) %>% arrange(dual_aic)

# ------------------------------------------------------------------------------
# STEP 3: DISPLAY RESULTS
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("          DUAL-ARENA PREDATION MODEL EVALUATION (24-YEAR DATASET)             \n")
cat("      SAR ~ CPUE + P_ocean_eff (Temp/Prey Gated) + InRiver_SeaLions (Raw)     \n")
cat("==============================================================================\n\n")

print(dual_arena_table, width = Inf)

out_csv <- file.path(data_out_dir, "Goal4_OceanPred_plus_EMBssl_24yr.csv")
write_csv(dual_arena_table, out_csv)
cat(sprintf("\nDual-arena evaluation saved to: %s\n", out_csv))


