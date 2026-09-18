library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Filter 1998-2021 window & evaluate NAs per column
# ---------------------------------------------------------------------------
df_window <- sem_altprey_data %>%
  filter(Year >= 1998 & Year <= 2021)

# Count NAs per column
na_counts <- colSums(is.na(df_window))

# Separate fully complete vs incomplete columns
complete_cols   <- names(na_counts[na_counts == 0])
incomplete_cols <- na_counts[na_counts > 0]

# ---------------------------------------------------------------------------
# Print Results
# ---------------------------------------------------------------------------
cat("=====================================================\n")
cat(sprintf("Total Columns Evaluated: %d (1998-2021, 24 Years)\n", ncol(df_window)))
cat(sprintf("Fully Complete (0 NAs):  %d columns\n", length(complete_cols)))
cat(sprintf("Incomplete (>0 NAs):     %d columns\n", length(incomplete_cols)))
cat("=====================================================\n\n")

cat("--- FULLY COMPLETE COLUMNS (1998-2021) ---\n")
print(complete_cols)

if (length(incomplete_cols) > 0) {
  cat("\n--- INCOMPLETE COLUMNS & NA COUNTS (1998-2021) ---\n")
  print(incomplete_cols)
}


# --- INCOMPLETE COLUMNS & NA COUNTS (1998-2021) ---
#   X06_Lmu_IntSprMayW                  X06_StomFull_May       X15_spinyDogfishBSAI_predAK 
# 1                                 4                                 5 
# X15_salmonSharkBSAI_predAK             X08_commonMurre_JSOES                X09_canaryRockfish 
# 5                                 5                                 5 
# X09_chilipepper X01_DFA_sumPreyOfPrey_planktonJun              X01_NHLlogSum_win_05 
# 5                                 1                                 3 
# X01_NHLlogSum_win_15              X01_NHLlogSum_win_25 
# 3                                 3 