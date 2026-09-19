
# --- SERIES WITH MISSING DATA IN TARGET WINDOW ---
#   # A tibble: 16 × 6
#   Column                               Total_Years Valid_Count NA_Count Percent_Complete Missing_Years
# 1 X10.PredMammalNCC_DFA1_smoltyr                27          26        1             96.3 1996         
# 2 X10.PredMammalNCC_0_smoothed_smoltyr          27          26        1             96.3 1996         
# 3 X10.PredMammalNCC_1_smoothed_smoltyr          27          26        1             96.3 1996         
# 4 X10.PredMammalNCC_2_smoothed_smoltyr          27          26        1             96.3 1996         
# 5 X10.PredMammalNCC_3_smoothed_smoltyr          27          26        1             96.3 1996         
# 6 X15.PredFishAK_0_smoothed_smoltyr             27          26        1             96.3 1996         
# 7 X15.PredFishAK_1_smoothed_smoltyr             27          26        1             96.3 1996         
# 8 X15.PredFishAK_2_smoothed_smoltyr             27          26        1             96.3 1996         
# 9 X15.PredFishAK_3_smoothed_smoltyr             27          26        1             96.3 1996         
# 10 X15.PredFishAK_4_smoothed_smoltyr             27          26        1             96.3 1996         
# 11 X15.PredFishAK_5_smoothed_smoltyr             27          26        1             96.3 1996         
# 12 X15.PredFishAK_b_DFA1_smoltyr                 27          26        1             96.3 1996         
# 13 X15.PredFishAK_b_0_smoothed_smoltyr           27          26        1             96.3 1996         
# 14 X15.PredFishAK_b_1_smoothed_smoltyr           27          26        1             96.3 1996         
# 15 X15.PredFishAK_b_2_smoothed_smoltyr           27          26        1             96.3 1996         
# 16 X15.PredFishAK_b_3_smoothed_smoltyr           27          26        1             96.3 1996   


# ==============================================================================
# Script: 03_audit_smoothed_extended_matrix.R
# Purpose: Check NA coverage across the final smoothed & shifted matrix
# Output: metadata/smoothed_extended_completeness_report.csv
# ==============================================================================

library(tidyverse)

rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
outputDir <- file.path(rootdir, "metadata")

matrix_file <- file.path(outputDir, "datWide_1996_2025_extended_shifted.csv")

if (!file.exists(matrix_file)) stop("Missing matrix file. Run script 02 first!")

dat <- read.csv(matrix_file, stringsAsFactors = FALSE)

# Define modeling target window (1996 through 2022/2023)
target_start <- 1996
target_end   <- 2022  # Set to 2022 for smolt year models, or 2023 for extended leads

dat_sub <- dat %>% filter(Year >= target_start & Year <= target_end)

# Audit each column
na_report <- lapply(setdiff(names(dat_sub), "Year"), function(col_name) {
  vec <- dat_sub[[col_name]]
  missing_years <- dat_sub$Year[is.na(vec)]
  
  tibble(
    Column = col_name,
    Total_Years = length(vec),
    Valid_Count = sum(!is.na(vec)),
    NA_Count = sum(is.na(vec)),
    Percent_Complete = round((sum(!is.na(vec)) / length(vec)) * 100, 1),
    Missing_Years = if (length(missing_years) > 0) paste(missing_years, collapse = ", ") else "None"
  )
}) %>% bind_rows()

cat("\n==============================================================================\n")
cat(sprintf("COMPLETENESS AUDIT FOR FINAL SMOOTHED MATRIX (%d - %d)\n", target_start, target_end))
cat("==============================================================================\n")
cat(sprintf("Total Series Evaluated: %d\n", nrow(na_report)))
cat(sprintf("100%% Complete Series: %d\n", sum(na_report$NA_Count == 0)))
cat(sprintf("Series with Missing Data: %d\n\n", sum(na_report$NA_Count > 0)))

incomplete_summary <- na_report %>% 
  filter(NA_Count > 0) %>% 
  arrange(desc(NA_Count))

if (nrow(incomplete_summary) > 0) {
  cat("--- SERIES WITH MISSING DATA IN TARGET WINDOW ---\n")
  print(incomplete_summary, n = 100)
} else {
  cat("SUCCESS: Every single smoothed series is 100% complete across the target window!\n")
}

# Export report
out_csv <- file.path(outputDir, "smoothed_extended_completeness_report.csv")
write.csv(na_report, out_csv, row.names = FALSE)
cat(sprintf("\nAudit report saved to: %s\n", out_csv))

#----------
library(tidyverse)

dat <- read.csv("metadata/datWide_1996_2025_extended_shifted.csv")
raw <- read.csv("metadata/datWide_1995_2025_reprocessed_raw.csv")

# Print the shifted predator column alongside raw unshifted values for 1998-2002
dat %>% 
#  select(Year, contains("10.PredMammalNCC_DFA1")) %>% 
  select(Year, contains("08.PredBirdNCC_b_0")) %>% 
  filter(Year >= 1998 & Year <= 2002)


library(tidyverse)

raw <- read.csv(file.path(rootdir, "metadata","datWide_1995_2025_reprocessed_raw.csv"))

# Check raw unshifted observations for 08.PredBirdNCC_b_0 across 1995-2002
raw %>% 
  select(Year, contains("murre")) %>% 
  filter(Year >= 1995 & Year <= 2002)
