# ==============================================================================
# Script: 07_create_ak_predator_comparison_table.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Build comparative summary table across all Alaska predators,
#          ordered by baseline magnitude, grouped by coefficient sign (+/-),
#          and tracking mean net impact shifts under alternate prey moderation.
# Output: Rcode_for_paper/Routput_for_paper/data/AK_Predator_AltPrey_Comparison_Table.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD PREVIOUS OUTPUTS
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

# Apply uniform +5.0 domain shift across numeric columns (excluding Year)
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
# STEP 1: EXTRACT BASELINE & BEST ALTPREY MODEL PER ALASKA PREDATOR
# ------------------------------------------------------------------------------
cat("Extracting baseline vs best alternate prey model comparisons...\n")

unique_predators <- unique(all_sweep_df$predator)
comparison_list  <- list()

for (pred_col in unique_predators) {
  
  P_pos <- sem_complete_data[[pred_col]]
  SAR   <- sem_complete_data[[target_sar]]
  
  # --- 1A. BASELINE UNMODERATED MODEL ---
  base_row <- all_sweep_df %>%
    filter(predator == pred_col & moderator == "None (Unattenuated)") %>%
    slice(1)
  
  base_coef <- as.numeric(base_row$pred_est)
  base_p    <- as.numeric(base_row$pvalue)
  base_aic  <- as.numeric(base_row$aic)
  
  base_net_ts     <- base_coef * P_pos
  base_mean_net   <- mean(base_net_ts, na.rm = TRUE)
  base_std_beta   <- base_coef * (sd(P_pos, na.rm = TRUE) / sd(SAR, na.rm = TRUE))
  
  # --- 1B. BEST ALTERNATE PREY / MODERATED MODEL (Lowest AIC regardless of dAIC threshold) ---
  alt_row <- all_sweep_df %>%
    filter(predator == pred_col & moderator != "None (Unattenuated)") %>%
    arrange(aic) %>%
    slice(1)
  
  if (nrow(alt_row) > 0) {
    alt_mod   <- alt_row$moderator
    alt_w     <- as.numeric(alt_row$best_w)
    alt_g     <- as.numeric(alt_row$best_gamma)
    alt_coef  <- as.numeric(alt_row$pred_est)
    alt_p     <- as.numeric(alt_row$pvalue)
    alt_aic   <- as.numeric(alt_row$aic)
    daic      <- as.numeric(alt_row$delta_aic)
    
    A_raw  <- sem_complete_data[[alt_mod]]
    A_norm <- min_max_scale(A_raw)
    weight_mod <- 1 - (alt_w * logistic_gate(A_norm, alt_g))
    P_eff  <- P_pos * weight_mod
    
    alt_net_ts     <- alt_coef * P_eff
    alt_mean_net   <- mean(alt_net_ts, na.rm = TRUE)
    alt_std_beta   <- alt_coef * (sd(P_eff, na.rm = TRUE) / sd(SAR, na.rm = TRUE))
    realized_damp  <- (1 - (sum(P_eff, na.rm = TRUE) / sum(P_pos, na.rm = TRUE))) * 100
  } else {
    alt_mod       <- "None Available"
    alt_w         <- 0; alt_g <- 0
    alt_coef      <- base_coef; alt_p <- base_p; alt_aic <- base_aic; daic <- 0
    alt_mean_net  <- base_mean_net; alt_std_beta <- base_std_beta; realized_damp <- 0
  }
  
  # --- 1C. NET IMPACT SHIFT METRICS ---
  net_shift      <- alt_mean_net - base_mean_net
  shift_direction <- case_when(
    abs(net_shift) < 0.001 ~ "No Change",
    net_shift > 0          ~ "Pushed Positive (Buffered Mortality)",
    net_shift < 0          ~ "Pushed Negative (Amplify Mortality)"
  )
  
  comparison_list[[length(comparison_list) + 1]] <- tibble(
    predator            = pred_col,
    base_group          = if_else(base_coef < 0, "Negative Base Effect (Mortality)", "Positive Base Effect (Co-varying)"),
    base_coef           = base_coef,
    base_abs_mag        = abs(base_coef),
    base_pvalue         = base_p,
    base_mean_net       = base_mean_net,
    best_altprey_mod    = alt_mod,
    altprey_w           = alt_w,
    altprey_gamma       = alt_g,
    altprey_coef        = alt_coef,
    altprey_pvalue      = alt_p,
    altprey_mean_net    = alt_mean_net,
    realized_damping_pct= realized_damp,
    delta_net_impact    = net_shift,
    impact_shift_dir    = shift_direction,
    base_aic            = base_aic,
    altprey_aic         = alt_aic,
    delta_aic           = daic
  )
}

comparison_df <- map_dfr(comparison_list, ~ .x)

# ------------------------------------------------------------------------------
# STEP 2: GROUP & ORDER BY ASCENDING MAGNITUDE OF BASELINE COEFFICIENT
# ------------------------------------------------------------------------------
sorted_table <- comparison_df %>%
  arrange(
    # 1. Group by Negative vs Positive Base Effect
    desc(base_group), 
    # 2. Order by ASCENDING magnitude of baseline coefficient |beta|
    base_abs_mag
  )

# ------------------------------------------------------------------------------
# STEP 3: DISPLAY & EXPORT COMPARISON TABLE
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("          ALASKA PREDATORS: BASELINE VS. ALTPREY COMPARISON TABLE             \n")
cat("==============================================================================\n\n")

cat("--- GROUP 1: NEGATIVE BASE EFFECT (ORDERED BY ASCENDING |COEF|) ---\n")
print(sorted_table %>% filter(base_group == "Negative Base Effect (Mortality)") %>%
        select(predator, base_coef, base_mean_net, best_altprey_mod, altprey_mean_net, delta_net_impact, impact_shift_dir, delta_aic),
      n = 30)

cat("\n--- GROUP 2: POSITIVE BASE EFFECT (ORDERED BY ASCENDING |COEF|) ---\n")
print(sorted_table %>% filter(base_group == "Positive Base Effect (Co-varying)") %>%
        select(predator, base_coef, base_mean_net, best_altprey_mod, altprey_mean_net, delta_net_impact, impact_shift_dir, delta_aic),
      n = 30)

out_table_path <- file.path(data_out_dir, "AK_Predator_AltPrey_Comparison_Table.csv")
write_csv(sorted_table, out_table_path)

cat("\nComparison table successfully exported to:\n", out_table_path, "\n")