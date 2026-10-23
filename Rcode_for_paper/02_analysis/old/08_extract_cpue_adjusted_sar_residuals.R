# ==============================================================================
# Script: 08_extract_cpue_adjusted_sar_residuals.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Fit baseline CPUE -> SAR model, extract residual time series,
#          and generate diagnostic plots to visualize variance to be explained.
# Output: Rcode_for_paper/Routput_for_paper/data/SAR_CPUE_Residuals_1998_2021.csv
# ==============================================================================

# When you fit sar~cpue, the residual represents the CPUE-adjusted marine survival anomaly for each year (1998–2021). 
# Positive residuals mean survival was higher than expected based on juvenile CPUE alone; 
# negative residuals mean survival was lower than expected (where predator mortality and environmental buffering play their role).

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
figure_out_dir <- file.path(master_out, "figures")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

sem_data_file <- file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv")
if (!file.exists(sem_data_file)) stop("Missing input file: ", sem_data_file)

sem_complete_data <- read.csv(sem_data_file, stringsAsFactors = FALSE) %>% 
  filter(complete.cases(.))

# Apply uniform +5.0 domain shift across numeric indicators (matching pipeline)
sem_plus5 <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"

# ------------------------------------------------------------------------------
# STEP 1: FIT BASELINE CPUE -> SAR MODEL & EXTRACT RESIDUALS
# ------------------------------------------------------------------------------
# Simple regression formulation matching the SEM stage link
baseline_lm <- lm(as.formula(paste(target_sar, "~", target_cpue)), data = sem_complete_data)

# Extract fitted values, raw residuals, and standard errors
res_df <- sem_complete_data %>%
  select(Year, SAR_observed = all_of(target_sar), CPUE_observed = all_of(target_cpue)) %>%
  mutate(
    SAR_fitted   = fitted(baseline_lm),
    SAR_residual = residuals(baseline_lm) # Observed SAR - Fitted SAR
  )

cat("Baseline Model Summary (SAR ~ CPUE):\n")
print(summary(baseline_lm))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.3824 -0.4573 -0.1088  0.5114  1.4446 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)              7.388e-17  1.746e-01   0.000  1.00000   
# X07_DFA_cpue_IntSprJunHW 5.479e-01  1.783e-01   3.072  0.00557 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8553 on 22 degrees of freedom
# Multiple R-squared:  0.3002,	Adjusted R-squared:  0.2684 
# F-statistic:  9.44 on 1 and 22 DF,  p-value: 0.005573

# ------------------------------------------------------------------------------
# STEP 2: PLOT RESIDUAL TIME SERIES
# ------------------------------------------------------------------------------
p_res <- ggplot(res_df, aes(x = Year, y = SAR_residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "navy", size = 2.5) +
  labs(
    title = "CPUE-Adjusted SAR Residuals (1998–2021)",
    subtitle = "Unexplained variance in salmon marine survival that predator & alternate prey models target",
    x = "Year",
    y = "SAR Residual (Observed - Predicted from CPUE)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_res)
ggsave(file.path(figure_out_dir, "SAR_CPUE_Residuals_1998_2021.png"),p_res)

# ------------------------------------------------------------------------------
# STEP 3: EXPORT RESIDUALS
# ------------------------------------------------------------------------------
out_res_path <- file.path(data_out_dir, "SAR_CPUE_Residuals_1998_2021.csv")
write_csv(res_df, out_res_path)

cat("\nExtracted residuals exported to:\n", out_res_path, "\n")
