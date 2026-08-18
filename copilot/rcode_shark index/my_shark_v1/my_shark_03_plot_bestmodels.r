#Best models for sharks


suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# -----------------------------------------------------------------------------
# 1. Load & Base Data Preparation
# -----------------------------------------------------------------------------
sem_master_data <- read.csv("data_Lisa/sem_master_data.csv", row.names = NULL) %>% clean_names()
ak_dat          <- read.csv("data_Lisa/data_all_tested_columns_annual.csv")%>% clean_names()
names(ak_dat)

names(sem_master_data)
data_base <- sem_master_data %>%
  select(year, x09_dfa_hake_age5plus) %>%
  left_join(ak_dat, by = "year") %>%
  dplyr::rename(
    Hake      = x09_dfa_hake_age5plus,
    CPUE_Juv  = x07_dfa_cpue_int_spr_jun_hw,
    SAR_Adult = x16_sar
  )



safe_scale <- function(x) as.vector(scale(x))

# -----------------------------------------------------------------------------
# 2. Fixed Parameters from Top Sweep Models
# -----------------------------------------------------------------------------
shark_col     <- "pacific_sleeper_shark_goa" # Top shark abundance column
Q10           <- 2.0
overlap_form  <- "logistic"
overlap_slope <- 0.8  # Using 0.8 from top model

# Parameter Combinations Requested (2 Temp x 3 Transformations = 6 Combinations)
combos <- tidyr::crossing(
  shark_temp_col  = c("enso_dj", "sst_wgoa_coastwatch_junjulaug"),
#  shark_transform = c("z", "log1p_z", "z_roll2")
  shark_transform = c("z_roll2")
)

# -----------------------------------------------------------------------------
# Subset Base Data to ONLY Tested Columns
# -----------------------------------------------------------------------------
tested_cols <- unique(c("year", "Hake", "CPUE_Juv", "SAR_Adult", shark_col, combos$shark_temp_col))

data_base <- data_base %>%
  select(all_of(tested_cols))
# -----------------------------------------------------------------------------
# 3. Calculate Weighted Shark Index & Fit Models Across Combinations
# -----------------------------------------------------------------------------
results_list <- list()

for (i in seq_len(nrow(combos))) {
  temp_col  <- combos$shark_temp_col[i]
  trans_type <- combos$shark_transform[i]
  scenario_id <- paste0(temp_col, " | ", trans_type)
  
  df <- data_base %>%
    filter(year >= 1998, year <= 2021) %>%
    mutate(
      Shark_raw = .data[[shark_col]],
      Temp_raw  = .data[[temp_col]]
    )
  
  # A. Calculate Abundance Transform (Shark_use)
  if (trans_type == "z") {
    df$Shark_use <- safe_scale(df$Shark_raw)
  } else if (trans_type == "log1p_z") {
    minv <- min(df$Shark_raw, na.rm = TRUE)
    df$Shark_use <- safe_scale(log1p(df$Shark_raw - minv))
  } else if (trans_type == "z_roll2") {
    z0 <- safe_scale(df$Shark_raw)
    df$Shark_use <- rowMeans(cbind(z0, dplyr::lag(z0, 1)), na.rm = TRUE)
  }
  
#   # Lead by 1 year: averages year t (smolt year) and year t+1 (ocean residency year)
#   zlead <- dplyr::lead(z0, 1) 
#   df$Shark_use <- rowMeans(cbind(z0, zlead), na.rm = TRUE)
# }
  
  # B. Calculate Metabolic Cost Factor (M_t) & Overlap Factor (O_t)
  T_scale <- safe_scale(df$Temp_raw)
  T_ref   <- mean(df$Temp_raw, na.rm = TRUE)
  
  df$M_t     <- Q10^((df$Temp_raw - T_ref) / 10)
  df$O_t     <- plogis(overlap_slope * T_scale) * 2  # Logistic form
  df$I_Shark <- df$Shark_use * df$M_t * df$O_t       # Weighted Index
  
  # C. Fit SEM with lavaan
  sem_spec <- '
    CPUE_Juv  ~ Hake
    SAR_Adult ~ CPUE_Juv + I_Shark
  '
  
  fit <- sem(sem_spec, data = df, missing = "ML", warn = FALSE)
  
  # D. Extract Parameter Estimates & Calculate Model Predictions
  pe <- parameterEstimates(fit)
  b0_sar   <- pe$est[pe$lhs == "SAR_Adult" & pe$op == "~1"]
  b_cpue   <- pe$est[pe$lhs == "SAR_Adult" & pe$op == "~" & pe$rhs == "CPUE_Juv"]
  b_shark  <- pe$est[pe$lhs == "SAR_Adult" & pe$op == "~" & pe$rhs == "I_Shark"]
  
  b0_cpue  <- pe$est[pe$lhs == "CPUE_Juv" & pe$op == "~1"]
  b_hake   <- pe$est[pe$lhs == "CPUE_Juv" & pe$op == "~" & pe$rhs == "Hake"]
  
  df <- df %>%
    mutate(
      Scenario    = scenario_id,
      Pred_CPUE   = b0_cpue + b_hake * Hake,
      Pred_SAR    = b0_sar + b_cpue * Pred_CPUE + b_shark * I_Shark
    )
  
  results_list[[i]] <- df
}

all_sims <- bind_rows(results_list)

write.csv(all_sims,file.path(out_dir, "shark_predictions.csv"))

# -----------------------------------------------------------------------------
# 4. Plot Initial Values, Weighted Index, and Predicted SAR
# -----------------------------------------------------------------------------

# Plot A: Initial Temperature & Shark Abundance
p1 <- all_sims %>%
  select(year, Scenario, Temp_raw, Shark_raw) %>%
  pivot_longer(cols = c(Temp_raw, Shark_raw), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = year, y = Value, color = Metric)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Scenario, scales = "free_y", ncol = 3) +
  labs(title = "Initial Data Inputs by Scenario", x = "Year", y = "Raw Value") +
  theme_minimal() +
  theme(legend.position = "top")

# Plot B: Calculated Weighted Shark Index (I_Shark)
p2 <- ggplot(all_sims, aes(x = year, y = I_Shark, color = Scenario)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario, ncol = 3) +
  labs(title = "Calculated Weighted Shark Index (I_Shark)", x = "Year", y = "I_Shark") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot C: Predicted vs Observed Salmon Response (SAR)
p3 <- all_sims %>%
  select(year, Scenario, Observed = SAR_Adult, Predicted = Pred_SAR) %>%
  pivot_longer(cols = c(Observed, Predicted), names_to = "Type", values_to = "SAR") %>%
  ggplot(aes(x = year, y = SAR, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario, ncol = 3) +
  labs(title = "Predicted vs Observed Salmon Response (SAR)", x = "Year", y = "SAR Value") +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "firebrick")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display Plots
print(p1)
print(p2)
print(p3)
plot(all_sims$I_Shark,all_sims$Pred_SAR)