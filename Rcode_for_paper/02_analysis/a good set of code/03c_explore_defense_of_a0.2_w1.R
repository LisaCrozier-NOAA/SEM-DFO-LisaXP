

#Defense of setting a=1:
#Why $R^2$, AIC, and $p$-values are Identical Despite Different a_preyLook closely at the numbers between the two runs:AIC, $\Delta\text{AIC}$, $R^2$, and $p$-value are identical to multiple decimal places ($-30.8$, $-4.41$, $0.578$, $0.00136$).beta_peff, however, changed:In the minimum AIC table (a_prey = 0.1): $\beta = -0.594$In the maximum $R^2$ table (a_prey = 0.5): $\beta = -1.07$This happens because altering $a_{\text{prey}}$ is performing an affine (linear) transformation on the predictor:$$P_{\text{eff}} = P_{\text{temp}} \times (1 - \text{prey}_{\text{scaled}})$$When you change the lower scaling bound $a_{\text{prey}}$ from $0.1$ to $0.5$, $\text{prey}_{\text{scaled}}$ undergoes a pure linear compression:$$\text{prey}_{\text{scaled}} = a + (1 - a) \times \frac{\text{raw} - \min}{\max - \min}$$In ordinary least squares and path regression (SEM without nonlinear latent interactions), multiplying a predictor by a scalar or shifting it by a constant does not change the model fit:The correlation between $P_{\text{eff}}$ and $\text{SAR}$ remains constant.$R^2$, log-likelihood, AIC, and test statistics ($t$/$z$ and $p$-values) remain unchanged.The regression coefficient $\beta$ rescales inversely ($\beta_{\text{new}} = \beta_{\text{old}} \times \text{scale factor}$) to absorb the change in variance.Because the optimizer found the exact same fit quality across multiple $a_{\text{prey}}$ values, slice_min and slice_max simply broke ties by selecting whichever row appeared first or matched tie-breaking criteria in the grid.

# ==============================================================================
# Script: 12_fine_grid_optimization_top_candidates_24yr.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: High-resolution parameter grid search (bounds 'a' and weights 'w')
#          on both SST and alternative prey for top 24yr candidates.
# Target: Maximize R2 and minimize AIC for SAR ~ CPUE + Peff
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA PREPARATION (24-YEAR DATASET)
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper", "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Filter complete 24-year series
valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) == 24]
df_raw <- raw_df_master %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

df_base <- df_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar  <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TOP CANDIDATE SPECIFICATIONS
# ------------------------------------------------------------------------------
# Top candidate configurations to optimize
top_candidates <- tribble(
  ~predator,                                 ~prey_vars,
  "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", "capelin_WGoA_smoltyr",
  "X15_sablefishAK_biomass4plus_smoltyr",    "sitkaHerring_EGoA_smoltyr",
  "X11_ssl_seak_pup_pred_smoltyr",           "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr"
  # Add candidates 4 & 5 here if tested:
  # "pred_col_4",                            "prey_col_4",
  # "pred_col_5",                            "prey_col_5"
)

# ------------------------------------------------------------------------------
# STEP 2: HIGH-RESOLUTION GRID EXPLORATION FUNCTION
# ------------------------------------------------------------------------------
run_fine_grid_optimization <- function(pred_col, prey_str, df_raw_data, df_base_data) {
  
  P_norm <- df_base_data[[pred_col]]
  
  # Determine spatial driver (standard EGoA summer SST for AK smolts)
  spatial_col <- "X21_sst_egoa_junjulaug_smoltyr"
  if (!spatial_col %in% names(df_raw_data)) {
    warning("Missing spatial col: ", spatial_col)
    return(NULL)
  }
  raw_temp <- df_raw_data[[spatial_col]]
  
  # Prepare unscaled forage pool
  prey_cols <- trimws(unlist(strsplit(prey_str, "\\+")))
  if (length(prey_cols) == 1) {
    raw_forage <- df_base_data[[prey_cols[1]]]
  } else {
    raw_forage <- rowSums(df_base_data[, prey_cols, drop = FALSE], na.rm = TRUE)
  }
  
  # --- Baseline Model (Unmodified Raw Predator) ---
  base_sem <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + ", pred_col)
  base_fit <- tryCatch(sem(base_sem, data = df_base_data, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  base_aic <- as.numeric(fitMeasures(base_fit, "aic"))
  base_r2  <- as.numeric(inspect(base_fit, "r2")[[target_sar]])
  base_pe  <- parameterEstimates(base_fit, standardized = TRUE) %>% 
    filter(op == "~" & lhs == target_sar & rhs == pred_col)
  base_beta <- as.numeric(base_pe$est[1])
  
  # Define Parameter Search Grids
  # Bounds 'a': scaling lower bound (0.0 = full range, 0.5 = mild compression)
  # Weights 'w': attenuation / gating strength (0.1 to 1.0)
  a_temp_grid <- seq(0.0, 0.5, by = 0.1)
  w_temp_grid <- seq(0.2, 1.0, by = 0.2)
  a_prey_grid <- seq(0.0, 0.5, by = 0.1)
  w_prey_grid <- seq(0.2, 1.0, by = 0.2)
  
  cat(sprintf("Optimizing: %s | Prey: %s\n", pred_col, prey_str))
  cat(sprintf("Total parameter space: %d combinations\n", 
              length(a_temp_grid) * length(w_temp_grid) * length(a_prey_grid) * length(w_prey_grid)))
  
  grid_results <- list()
  
  for (at in a_temp_grid) {
    temp_scaled <- scale_custom_bounds(raw_temp, a = at, b = 1.0)
    
    for (wt in w_temp_grid) {
      gate_temp <- 1 - wt * (1 - temp_scaled)
      P_temp    <- P_norm * gate_temp
      
      for (ap in a_prey_grid) {
        prey_scaled <- scale_custom_bounds(raw_forage, a = ap, b = 1.0)
        
        for (wp in w_prey_grid) {
          gate_prey <- 1 - wp * prey_scaled
          P_joint   <- P_temp * gate_prey
          
          if (is.na(sd(P_joint)) || sd(P_joint) < 1e-6) next
          
          df_base_data$P_fine_eff <- P_joint
          sem_model <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_fine_eff")
          
          fit <- tryCatch(
            sem(sem_model, data = df_base_data, std.lv = TRUE, missing = "ML", warn = FALSE), 
            error = function(e) NULL
          )
          
          if (!is.null(fit) && lavInspect(fit, "converged")) {
            pe <- parameterEstimates(fit, standardized = TRUE) %>% 
              filter(op == "~" & lhs == target_sar & rhs == "P_fine_eff")
            
            if (nrow(pe) > 0 && !is.na(pe$pvalue[1])) {
              curr_aic  <- as.numeric(fitMeasures(fit, "aic"))
              curr_r2   <- as.numeric(inspect(fit, "r2")[[target_sar]])
              curr_beta <- as.numeric(pe$est[1])
              curr_p    <- as.numeric(pe$pvalue[1])
              
              grid_results[[length(grid_results) + 1]] <- tibble(
                predator     = pred_col,
                prey_combo   = prey_str,
                a_temp       = at,
                w_temp       = wt,
                a_prey       = ap,
                w_prey       = wp,
                base_aic     = base_aic,
                base_r2      = base_r2,
                base_beta    = base_beta,
                aic          = curr_aic,
                delta_aic    = curr_aic - base_aic,
                r2           = curr_r2,
                r2_gain      = curr_r2 - base_r2,
                beta_peff    = curr_beta,
                beta_std     = as.numeric(pe$std.all[1]),
                pvalue       = curr_p
              )
            }
          }
        }
      }
    }
  }
  
  bind_rows(grid_results)
}

# ------------------------------------------------------------------------------
# STEP 3: EXECUTE OPTIMIZATION & SAVE DETAILED/SUMMARY RESULTS
# ------------------------------------------------------------------------------
all_grid_runs <- list()

for (i in seq_len(nrow(top_candidates))) {
  cand <- top_candidates[i, ]
  t0 <- Sys.time()
  
  res <- run_fine_grid_optimization(
    pred_col     = cand$predator,
    prey_str     = cand$prey_vars,
    df_raw_data  = df_raw,
    df_base_data = df_base
  )
  
  if (!is.null(res) && nrow(res) > 0) {
    all_grid_runs[[i]] <- res
    cat(sprintf("Candidate %d complete in %.1f seconds.\n", i, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

full_optimization_df <- bind_rows(all_grid_runs)

# Save the full parameter space
out_full_csv <- file.path(data_out_dir, "Goal3_TopCandidates_24yr_FineGrid_All.csv")
write_csv(full_optimization_df, out_full_csv)

# ------------------------------------------------------------------------------
# STEP 4: EXTRACT BEST MODELS (BY MAX R2 AND BY MIN AIC)
# ------------------------------------------------------------------------------
# Best by R2
best_by_r2 <- full_optimization_df %>%
  group_by(predator) %>%
  slice_max(order_by = r2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(r2))

# Best by AIC
best_by_aic <- full_optimization_df %>%
  group_by(predator) %>%
  slice_min(order_by = aic, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(aic)

cat("\n==============================================================================\n")
cat("      TOP CONFIGURATIONS OPTIMIZED BY MAXIMUM R2 (SAR EQUATION)               \n")
cat("==============================================================================\n")
print(best_by_r2 %>% select(predator, base_r2, r2, r2_gain, aic, delta_aic, beta_peff, pvalue, a_temp, w_temp, a_prey, w_prey), n = 10)

cat("\n==============================================================================\n")
cat("      TOP CONFIGURATIONS OPTIMIZED BY MINIMUM AIC                             \n")
cat("==============================================================================\n")
print(best_by_aic %>% select(predator, base_aic, aic, delta_aic, base_r2, r2, beta_peff, pvalue, a_temp, w_temp, a_prey, w_prey), n = 10)

write_csv(best_by_r2,  file.path(data_out_dir, "Goal3_TopCandidates_24yr_Best_R2.csv"))
write_csv(best_by_aic, file.path(data_out_dir, "Goal3_TopCandidates_24yr_Best_AIC.csv"))

# Inspect top models with an enforced ecological floor (a_temp >= 0.2)
full_optimization_df %>%
  filter(a_temp >= 0.2) %>%
  group_by(predator) %>%
  slice_max(order_by = r2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(r2)) %>%
  select(predator, base_r2, r2, r2_gain, aic, delta_aic, beta_peff, pvalue, a_temp, w_temp, a_prey, w_prey)
