# ==============================================================================
# Script: test_hake_hci_gating_standalone.R
# Purpose: Direct single-model comparison: Raw Hake vs. HCI-Gated Hake
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# --- 1. Load & Prep Data ---
years_included <- "24yr"
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper/data")

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Invert HCI (High = Warm / Compressed)
raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr

if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"
hci_col     <- "X01_habCompInd_smoltyr"

# --- 2. Raw Model (No Temperature) ---
raw_fit <- sem(paste0(target_cpue, " ~ ", hake_col), data = minmax_base, std.lv = TRUE, missing = "ML", warn = FALSE)
raw_aic <- as.numeric(fitMeasures(raw_fit, "aic"))
raw_pe  <- parameterEstimates(raw_fit, standardized = TRUE) %>% filter(op == "~" & rhs == hake_col)
raw_est <- as.numeric(raw_pe$est[1])
raw_p   <- as.numeric(raw_pe$pvalue[1])
raw_r2  <- as.numeric(inspect(raw_fit, "r2")[[target_cpue]])

# --- 3. Grid Search: HCI-Gated Hake ---
P_norm  <- minmax_base[[hake_col]]
raw_hci <- filtered_df[[hci_col]]

a_bounds <- seq(0.0, 0.5, by = 0.1)
w_grid   <- seq(0.1, 1.0, by = 0.1)

best_hci_fit <- NULL
best_hci_aic <- Inf
best_a       <- NA
best_w       <- NA

for (a in a_bounds) {
  hci_scaled <- scale_custom_bounds(raw_hci, a = a, b = 1.0)
  for (w in w_grid) {
    gate_factor <- 1 - w * (1 - hci_scaled)
    P_eff       <- P_norm * gate_factor
    
    if (is.na(sd(P_eff)) || sd(P_eff) < 1e-6) next
    
    df_temp <- minmax_base
    df_temp$P_hake_hci <- P_eff
    
    fit <- tryCatch(
      sem(paste0(target_cpue, " ~ P_hake_hci"), data = df_temp, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      curr_aic <- as.numeric(fitMeasures(fit, "aic"))
      if (curr_aic < best_hci_aic) {
        best_hci_aic <- curr_aic
        best_hci_fit <- fit
        best_a       <- a
        best_w       <- w
      }
    }
  }
}

best_pe  <- parameterEstimates(best_hci_fit, standardized = TRUE) %>% filter(op == "~" & rhs == "P_hake_hci")
best_est <- as.numeric(best_pe$est[1])
best_p   <- as.numeric(best_pe$pvalue[1])
best_r2  <- as.numeric(inspect(best_hci_fit, "r2")[[target_cpue]])

# --- 4. Print Side-by-Side Comparison ---
cat("\n==================================================================\n")
cat("      HAKE AGE 5+: RAW vs. HCI-GATED MODEL COMPARISON             \n")
cat("==================================================================\n\n")

comp_table <- tibble(
  Model           = c("Raw Hake (No Temperature)", "HCI-Gated Hake (Best Grid)"),
  MinMax_a        = c(NA, best_a),
  Weight_w        = c(NA, best_w),
  Beta_Hake       = c(raw_est, best_est),
  p_value         = c(raw_p, best_p),
  R2              = c(raw_r2, best_r2),
  AIC             = c(raw_aic, best_hci_aic),
  Delta_AIC       = c(0, best_hci_aic - raw_aic)
)

print(comp_table)
# Model                      MinMax_a Weight_w Beta_Hake  p_value    R2   AIC Delta_AIC
# 1 Raw Hake (No Temperature)      NA       NA      -0.555 0.000156 0.373 -7.44     0    
# 2 HCI-Gated Hake (Best Grid)      0.5      0.1    -0.550 0.000184 0.368 -7.25     0.196