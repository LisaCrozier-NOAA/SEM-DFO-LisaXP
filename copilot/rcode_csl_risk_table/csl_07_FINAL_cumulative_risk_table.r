#CUMULATIVE RISK TABLE
#write.csv(csl_risk_table_classified, file.path(output_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv"), row.names = FALSE)

# Table: Integrated CSL Predation Risk Table (1998–2024)
# 
# | Year |  CSL Count|Estuary Risk      | Ocean Buffer |Integrated Risk |        Hake Tier         |       Sardine Tier       |       Anchovy Tier       |        Sitka Tier        |       Herring Tier       |
#   |:----:|----------:|:-----------------|:------------:|:---------------|:------------------------:|:------------------------:|:------------------------:|:------------------------:|:------------------------:|
#   | 1998 |   581.8404|LOW RISK          |      3       |LOW RISK        |     Moderate Biomass     |  High Biomass (Buffer)   | Low Biomass (Upper Risk) | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |
#   | 1999 |   238.1267|LOW RISK          |      2       |LOW RISK        | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |
#   | 2000 |   280.3047|LOW RISK          |      2       |LOW RISK        | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) | Low Biomass (Upper Risk) |     Moderate Biomass     |
#   | 2001 |   865.6871|LOW RISK          |      4       |LOW RISK        |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     | Low Biomass (Upper Risk) |     Moderate Biomass     |
#   | 2002 |  1183.3938|LOW-MODERATE RISK |      4       |LOW RISK        |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     | Low Biomass (Upper Risk) |     Moderate Biomass     |
#   | 2003 |  1897.1685|LOW-MODERATE RISK |      5       |LOW RISK        |     Moderate Biomass     |  High Biomass (Buffer)   |     Moderate Biomass     |     Moderate Biomass     |  High Biomass (Buffer)   |
#   | 2004 |   652.0000|LOW RISK          |      5       |LOW RISK        |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |
#   | 2005 |   841.0000|LOW RISK          |      5       |LOW RISK        |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |     Moderate Biomass     |     Moderate Biomass     |
#   | 2006 |   535.0000|LOW RISK          |      4       |LOW RISK        |     Moderate Biomass     |  High Biomass (Buffer)   |  High Biomass (Buffer)   |     Moderate Biomass     | Low Biomass (Upper Risk) |
#   | 2007 |   804.0000|LOW RISK          |      4       |LOW RISK        | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     |
#   | 2008 |   822.0000|LOW RISK          |      3       |LOW RISK        | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |     Moderate Biomass     |
#   | 2009 |   761.0000|LOW RISK          |      2       |LOW RISK        | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |
#   | 2010 |  1450.0000|MODERATE RISK     |      2       |MODERATE RISK   | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |
#   | 2011 |   899.0000|LOW RISK          |      2       |LOW RISK        | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |
#   | 2012 |  1017.0000|MODERATE RISK     |      2       |MODERATE RISK   | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) |
#   | 2013 |  4179.0000|HIGH RISK         |      5       |HIGH RISK       |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     |
#   | 2014 |  4603.0000|HIGH RISK         |      4       |HIGH RISK       |     Moderate Biomass     | Low Biomass (Upper Risk) |     Moderate Biomass     |     Moderate Biomass     |     Moderate Biomass     |
#   | 2015 | 10640.0000|HIGH RISK         |      2       |HIGH RISK       |     Moderate Biomass     | Low Biomass (Upper Risk) |     Moderate Biomass     | Low Biomass (Upper Risk) | Low Biomass (Upper Risk) |
#   | 2016 |  7879.7396|VERY HIGH RISK    |      3       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |     Moderate Biomass     |     Moderate Biomass     | Low Biomass (Upper Risk) |
#   | 2017 |  6112.0000|VERY HIGH RISK    |      2       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) | Low Biomass (Upper Risk) |
#   | 2018 |  5226.8941|VERY HIGH RISK    |      3       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |
#   | 2019 |  5297.5871|VERY HIGH RISK    |      4       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |
#   | 2020 |  5197.3303|VERY HIGH RISK    |      4       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |
#   | 2021 |  4463.5882|HIGH RISK         |      4       |HIGH RISK       |  High Biomass (Buffer)   | Low Biomass (Upper Risk) |  High Biomass (Buffer)   |  High Biomass (Buffer)   |  High Biomass (Buffer)   |
#   > 




# ==============================================================================
# Script: CSL Risk Table Integration with Final 5 Ocean Prey Time Series
# Merges Estuary Chinook Exposure Risk with Quantile-Based Ocean Prey Tiers
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
# 1. Load Estuary Risk Data & Final 5 Ocean Prey Dataset
# -----------------------------------------------------------------------------

estuary_risk <- read.csv(file.path(output_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv"))
#estuary_risk_quantile <- read.csv(file.path(output_dir, "chinook_quantile_csl_exposure_risk_1998_2024.csv"))
ocean_prey   <- read.csv(file.path(output_dir, "csl_ocean_prey_5series_1998_2024.csv"))

# Merge estuary risk metrics with ocean prey metrics across 1998–2024
master_csl_risk_table <- estuary_risk %>%
  inner_join(ocean_prey, by = "year")

# -----------------------------------------------------------------------------
# 2. Apply Quantile Classification to Standardized Ocean Prey Metrics
# -----------------------------------------------------------------------------

# Function to categorize standardized continuous metrics into 3 Quantile Tiers
classify_quantile <- function(x) {
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  case_when(
    x <= q[2] ~ "Low Biomass (Upper Risk)",
    x <= q[3] ~ "Moderate Biomass",
    TRUE      ~ "High Biomass (Buffer)"
  )
}

csl_risk_table_classified <- master_csl_risk_table %>%
  mutate(
    # Quantile Tiers for Final 5 Ocean Prey Metrics
    hake_age5plus_tier = classify_quantile(x09_dfa_hake_age5plus_doug),
    sardine_tier       = classify_quantile(x05_dfa_abund_sardine_doug),
    sitka_herring_tier = classify_quantile(sitka_herring_e_go_a_doug),
    ncc_herring_tier   = classify_quantile(x05_herring_ncc_doug),
    ncc_anchovy_tier   = classify_quantile(x05_anchovy_gam_doug),
    
    # Composite Ocean Forage Buffer Score (Count of High/Moderate Ocean Prey Tiers)
    ocean_buffer_score = (hake_age5plus_tier != "Low Biomass (Upper Risk)") +
      (sardine_tier != "Low Biomass (Upper Risk)") +
      (sitka_herring_tier != "Low Biomass (Upper Risk)") +
      (ncc_herring_tier != "Low Biomass (Upper Risk)") +
      (ncc_anchovy_tier != "Low Biomass (Upper Risk)"),
    
    # Composite Multi-Prey Integrated Risk Tier
    integrated_csl_risk = case_when(
      overall_chinook_risk == "VERY HIGH RISK" & ocean_buffer_score <= 1 ~ "CRITICAL RISK",
      overall_chinook_risk == "VERY HIGH RISK" | overall_chinook_risk == "HIGH RISK" ~ "HIGH RISK",
      overall_chinook_risk == "MODERATE RISK" & ocean_buffer_score >= 4 ~ "LOW-MODERATE RISK",
      overall_chinook_risk == "MODERATE RISK" ~ "MODERATE RISK",
      TRUE ~ "LOW RISK"
    )
  )

# -----------------------------------------------------------------------------
# 3. Print Summary & Export Final Risk Table
# -----------------------------------------------------------------------------


cat("\n========================================================\n")
cat("   CSL PREDATION RISK TABLE WITH FINAL 5 OCEAN PREY     \n")
cat("========================================================\n")
# 1. Reorder risk_display_df columns to match your desired table headers
risk_display_df <- csl_risk_table_classified %>%
  select(
    year, 
    csl_during_chinook, 
    overall_chinook_risk,
    ocean_buffer_score,
    integrated_csl_risk,
    hake_age5plus_tier, 
    sardine_tier, 
    ncc_anchovy_tier,
    sitka_herring_tier, 
    ncc_herring_tier
  )

# 2. Print styled kable table with aligned column names
print(
  kable(
    risk_display_df,
    caption = "Integrated CSL Predation Risk Table (1998–2024)",
    col.names = c(
      "Year", "CSL Count", "Estuary Risk", 
      "Ocean Buffer", "Integrated Risk",
      "Hake Tier", "Sardine Tier", "Anchovy Tier", "Sitka Tier", "Herring Tier"
    ),
    align = c("c", "r", "l", "c", "l", "c", "c", "c", "c", "c")
  )
)

# Export complete table
write.csv(csl_risk_table_classified, file.path(output_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv"), row.names = FALSE)
out_dir <- "copilot/outputs_6"
write.csv(csl_risk_table_classified, file.path(out_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv"), row.names = FALSE)

cat("\nExport Complete: Saved 'master_csl_predation_risk_table_5ocean_prey_1998_2024.csv' to", output_dir, "\n")