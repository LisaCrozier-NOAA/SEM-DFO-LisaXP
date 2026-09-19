

# ==============================================================================
# Script: 03b_append_standalone_altprey_to_master.R
# Purpose: Merge verified standalone AltPrey constituents (X13_pollock_age1plus 
#          and X13_sitkaHerring_EGoA) into the master extended shifted LisaName matrix,
#          maintaining strict alphanumeric column ordering.
# Output: metadata/clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv
# ==============================================================================

library(tidyverse)

rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
outputDir <- file.path(rootdir, "metadata")

master_file   <- file.path(outputDir, "clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv")
altprey_file  <- file.path(outputDir, "standalone_altprey_constituents_1996_2025.csv")

if (!file.exists(master_file))  stop("Missing master file: ", master_file, ". Run Script 02 first!")
if (!file.exists(altprey_file)) stop("Missing standalone AltPrey file: ", altprey_file, ". Run Script 03a first!")

master_df  <- read.csv(master_file, stringsAsFactors = FALSE)
altprey_df <- read.csv(altprey_file, stringsAsFactors = FALSE)

# --- 1. REMOVE ANY PRE-EXISTING ALTPREY OVERLAPS & MERGE ---
altprey_cols <- setdiff(names(altprey_df), "Year")

updated_master <- master_df %>%
  select(-any_of(altprey_cols)) %>%
  left_join(altprey_df, by = "Year")

# --- 2. ENFORCE ALPHANUMERIC COLUMN ORDERING ---
# Year first, followed by X01_... through X16_SAR
lisa_sorted_cols <- c("Year", sort(setdiff(names(updated_master), "Year")))
updated_master   <- updated_master %>% select(all_of(lisa_sorted_cols))

# --- 3. EXPORT UPDATED MASTER MATRIX ---
write.csv(updated_master, master_file, row.names = FALSE)

cat("\n==============================================================================\n")
cat("SUCCESSFULLY MERGED STANDALONE ALTPREY CONSTITUENTS!\n")
cat("==============================================================================\n")
cat("Updated file:", master_file, "\n")
cat("Total Columns in Final Matrix:", ncol(updated_master), "\n")
cat("Added Columns:\n")
print(altprey_cols)
cat("\nFirst 10 Columns in Alphanumeric Order:\n")
print(names(updated_master)[1:10])
cat("==============================================================================\n")