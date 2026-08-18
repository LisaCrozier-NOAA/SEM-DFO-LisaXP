
#RESULTS -- EITHER IS FINE

#                     Quantile
# Numeric             HIGH LOW MODERATE MODERATE-LOW VERY HIGH
# HIGH RISK            3   0        1            1         1
# LOW-MODERATE RISK    0   0        0            3         0
# LOW RISK             0   9        1            1         0
# MODERATE RISK        0   0        2            0         0
# VERY HIGH RISK       0   0        0            0         5


# ==============================================================================
# Script: Compare Quantile vs. Numeric Risk Predictions (1998–2024)
# Evaluates side-by-side risk classification agreement across methods
# ==============================================================================

library(tidyverse)
library(janitor)
library(knitr)

output_dir <- "copilot/outputs_csl_cr"

# Safe Directory Check
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Load Both Prediction Datasets
# -----------------------------------------------------------------------------

estuary_risk_numeric  <- read.csv(file.path(output_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv"))
estuary_risk_quantile <- read.csv(file.path(output_dir, "chinook_quantile_csl_exposure_risk_1998_2024.csv"))

# -----------------------------------------------------------------------------
# 2. Join Predictions & Evaluate Agreement
# -----------------------------------------------------------------------------

# Detect risk column names automatically across both datasets
find_risk_col <- function(df) {
  cols <- colnames(df)
  risk_col <- cols[grepl("risk|tier", cols, ignore.case = TRUE) & !grepl("csl_|eulachon_|shad_", cols, ignore.case = TRUE)]
  if (length(risk_col) > 0) return(risk_col[1]) else return(cols[ncol(df)])
}

num_col_name <- find_risk_col(estuary_risk_numeric)
q_col_name   <- find_risk_col(estuary_risk_quantile)

comparison_df <- estuary_risk_numeric %>%
  select(year, numeric_risk = all_of(num_col_name)) %>%
  inner_join(
    estuary_risk_quantile %>% select(year, quantile_risk = all_of(q_col_name)),
    by = "year"
  ) %>%
  mutate(
    # Check exact match
    is_match = numeric_risk == quantile_risk,
    match_status = if_else(is_match, "Match", paste0("Mismatch (Numeric: ", numeric_risk, " vs Quantile: ", quantile_risk, ")"))
  )

# -----------------------------------------------------------------------------
# 3. Summary Statistics & Cross-Tabulation
# -----------------------------------------------------------------------------

total_years <- nrow(comparison_df)
matched_years <- sum(comparison_df$is_match)
pct_agreement <- round((matched_years / total_years) * 100, 1)

cat("\n========================================================\n")
cat("   RISK METHOD COMPARISON: QUANTILE VS. NUMERIC METHOD   \n")
cat("========================================================\n")
cat("Total Years Evaluated: ", total_years, "\n")
cat("Exact Predictions Matched: ", matched_years, " out of ", total_years, " (", pct_agreement, "% Agreement)\n", sep = "")
cat("========================================================\n\n")

# Confusion matrix / cross-tabulation table
cross_tab <- table(Numeric = comparison_df$numeric_risk, Quantile = comparison_df$quantile_risk)

cat("Cross-Tabulation Matrix (Numeric Rows vs Quantile Columns):\n")
print(cross_tab)
cat("\n")

#                     Quantile
# Numeric             HIGH LOW MODERATE MODERATE-LOW VERY HIGH
# HIGH RISK            3   0        1            1         1
# LOW-MODERATE RISK    0   0        0            3         0
# LOW RISK             0   9        1            1         0
# MODERATE RISK        0   0        2            0         0
# VERY HIGH RISK       0   0        0            0         5
# -----------------------------------------------------------------------------
# 4. Formatted Display Table (knitr::kable)
# -----------------------------------------------------------------------------

display_table <- comparison_df %>%
  select(
    Year = year,
    `Numeric Method Prediction` = numeric_risk,
    `Quantile Method Prediction` = quantile_risk,
    `Match Status` = match_status
  )

print(
  kable(
    display_table,
    caption = paste0("Year-by-Year Comparison of Numeric vs. Quantile Risk Predictions (Overall Agreement: ", pct_agreement, "%)"),
    align = c("c", "l", "l", "l")
  )
)

# Table: Year-by-Year Comparison of Numeric vs. Quantile Risk Predictions (Overall Agreement: 0%)
# 
# | Year |Numeric Method Prediction |Quantile Method Prediction |Match Status                                                    |
#   |:----:|:-------------------------|:--------------------------|:---------------------------------------------------------------|
#   | 1998 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 1999 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2000 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2001 |LOW RISK                  |MODERATE-LOW               |Mismatch (Numeric: LOW RISK vs Quantile: MODERATE-LOW)          |
#   | 2002 |LOW-MODERATE RISK         |MODERATE-LOW               |Mismatch (Numeric: LOW-MODERATE RISK vs Quantile: MODERATE-LOW) |
#   | 2003 |LOW-MODERATE RISK         |MODERATE-LOW               |Mismatch (Numeric: LOW-MODERATE RISK vs Quantile: MODERATE-LOW) |
#   | 2004 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2005 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2006 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2007 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2008 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2009 |LOW RISK                  |LOW                        |Mismatch (Numeric: LOW RISK vs Quantile: LOW)                   |
#   | 2010 |MODERATE RISK             |MODERATE                   |Mismatch (Numeric: MODERATE RISK vs Quantile: MODERATE)         |
#   | 2011 |LOW RISK                  |MODERATE                   |Mismatch (Numeric: LOW RISK vs Quantile: MODERATE)              |
#   | 2012 |MODERATE RISK             |MODERATE                   |Mismatch (Numeric: MODERATE RISK vs Quantile: MODERATE)         |
#   | 2013 |HIGH RISK                 |MODERATE                   |Mismatch (Numeric: HIGH RISK vs Quantile: MODERATE)             |
#   | 2014 |HIGH RISK                 |HIGH                       |Mismatch (Numeric: HIGH RISK vs Quantile: HIGH)                 |
#   | 2015 |HIGH RISK                 |VERY HIGH                  |Mismatch (Numeric: HIGH RISK vs Quantile: VERY HIGH)            |
#   | 2016 |VERY HIGH RISK            |VERY HIGH                  |Mismatch (Numeric: VERY HIGH RISK vs Quantile: VERY HIGH)       |
#   | 2017 |VERY HIGH RISK            |VERY HIGH                  |Mismatch (Numeric: VERY HIGH RISK vs Quantile: VERY HIGH)       |
#   | 2018 |VERY HIGH RISK            |VERY HIGH                  |Mismatch (Numeric: VERY HIGH RISK vs Quantile: VERY HIGH)       |
#   | 2019 |VERY HIGH RISK            |VERY HIGH                  |Mismatch (Numeric: VERY HIGH RISK vs Quantile: VERY HIGH)       |
#   | 2020 |VERY HIGH RISK            |VERY HIGH                  |Mismatch (Numeric: VERY HIGH RISK vs Quantile: VERY HIGH)       |
#   | 2021 |HIGH RISK                 |HIGH                       |Mismatch (Numeric: HIGH RISK vs Quantile: HIGH)                 |
#   | 2022 |LOW-MODERATE RISK         |MODERATE-LOW               |Mismatch (Numeric: LOW-MODERATE RISK vs Quantile: MODERATE-LOW) |
#   | 2023 |HIGH RISK                 |HIGH                       |Mismatch (Numeric: HIGH RISK vs Quantile: HIGH)                 |
#   | 2024 |HIGH RISK                 |MODERATE-LOW               |Mismatch (Numeric: HIGH RISK vs Quantile: MODERATE-LOW)         |
# -----------------------------------------------------------------------------
# 5. Export Results
# -----------------------------------------------------------------------------

write.csv(comparison_df, file.path(output_dir, "risk_prediction_comparison_quantile_vs_numeric.csv"), row.names = FALSE)

cat("\nPipeline Complete: Saved comparison output to 'risk_prediction_comparison_quantile_vs_numeric.csv' in", output_dir, "\n")
