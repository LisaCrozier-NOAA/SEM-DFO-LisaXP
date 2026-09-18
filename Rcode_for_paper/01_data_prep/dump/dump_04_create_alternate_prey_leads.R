


# ==============================================================================
# Script: 04_create_alternate_prey_leads.R
# Purpose: Derive 2-year lead variables for alternate prey & predator interaction
#          models while keeping baseline DFAs fixed on 1998–2021.
# ==============================================================================

library(tidyverse)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. Paths & Configuration
# ------------------------------------------------------------------------------
proj_dir   <- getwd()
input_dir  <- file.path(proj_dir, "output", "DFA")
output_dir <- file.path(proj_dir, "output", "interaction_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load master wide DFA matrix (contains data extended through 2023/2025)
# Note: Ensure this contains the extended years for the target columns
dfa_file <- file.path(input_dir, "clusDataDFA_wide_extended.rds") 
if (!file.exists(dfa_file)) {
  # Fallback to standard wide CSV if RDS is not present
  dfa_file <- file.path(input_dir, "completeness_check.csv")
}

clusData_lisa <- read.csv(dfa_file)

# Ensure Year column exists and is ordered chronologically
if ("date" %in% names(clusData_lisa)) {
  clusData_lisa$Year <- year(ymd(clusData_lisa$date))
}
clusData_lisa <- clusData_lisa %>% arrange(Year)

# ------------------------------------------------------------------------------
# 2. Target Columns for 2-Year Leads
# ------------------------------------------------------------------------------
prey_lead_cols <- c(
  "X04_marketsquid_GAM",
  "X05_DFA_abundSardine",         # Required for herring/sardine 2yrlead
  "X05_anchovy_GAM",
  "X09_DFA_ChinAbundSnakeFall",
  "X09_DFA_HakeAge5Plus",
  "X12_DFA_biomassEuphShelfSum",
  "X13_pollock_age1plus"
)

# Identify which specified columns exist in your dataset
target_cols <- intersect(prey_lead_cols, names(clusData_lisa))

cat("Found", length(target_cols), "of", length(prey_lead_cols), "target lead columns in dataset.\n")

# ------------------------------------------------------------------------------
# 3. Apply 2-Year Forward Lead Transformation
# ------------------------------------------------------------------------------
# Note: dplyr::lead(x, 2) shifts values back in time so that value at t+2 becomes value at t
clusData_led <- clusData_lisa %>%
  mutate(across(
    all_of(target_cols),
    ~ dplyr::lead(.x, 2),
    .names = "{.col}_2yrlead"
  ))

# ------------------------------------------------------------------------------
# 4. Construct Interaction Terms (Predator x Alternate Prey Lead)
# ------------------------------------------------------------------------------
# Example: Sea Lion * Capelin 2yrlead or Sea Lion * Herring 2yrlead
if ("PredAK_SeaLion" %in% names(clusData_led) && "X13_pollock_age1plus_2yrlead" %in% names(clusData_led)) {
  clusData_led <- clusData_led %>%
    mutate(
      SeaLion_x_Pollock_2yrlead = PredAK_SeaLion * X13_pollock_age1plus_2yrlead
    )
}

# ------------------------------------------------------------------------------
# 5. Trim Back to Baseline Evaluation Window (1998–2021) for SEM Analysis
# ------------------------------------------------------------------------------
# The leads were calculated using 2022–2023 data; now drop those outer years 
# so the evaluation matrix matches your main manuscript time frame.
final_sem_matrix <- clusData_led %>%
  filter(Year >= 1998, Year <= 2021)

# Save processed interaction matrix
saveRDS(final_sem_matrix, file.path(output_dir, "clusData_with_2yrleads_1998_2021.rds"))
write.csv(final_sem_matrix, file.path(output_dir, "clusData_with_2yrleads_1998_2021.csv"), row.names = FALSE)

cat("\nWorkflow Complete: 2-year leads constructed successfully.\n")
cat("Final matrix dimensions for SEM (1998-2021):", dim(final_sem_matrix)[1], "rows x", dim(final_sem_matrix)[2], "cols.\n")