# ==============================================================================
# Script: compare_clusData_to_doug_completeness.R
# Purpose: Compare local wide matrix against Doug's completeness.csv
# ==============================================================================

library(tidyverse)

# 1. Define File Paths
lisa_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_altprey/clusDataDFA_wide_1998_2021.csv"
doug_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/DAG1A_long/completeness.csv"

# 2. Read Files
df_lisa_raw <- read.csv(lisa_path, stringsAsFactors = FALSE)
df_doug_raw <- read.csv(doug_path, stringsAsFactors = FALSE)
df_doug_raw$year<-1998:2021

# 3. Standardize Year Keys & Clean Up Metadata Columns
# Lisa's format: Year is first
df_lisa <- df_lisa_raw %>%
  rename_with(~ "Year", matches("^year$", ignore.case = TRUE)) %>%
  arrange(Year)

# Doug's format: year is last or metadata present ("complete", date column)
df_doug <- df_doug_raw %>%
  rename_with(~ "Year", matches("^year$", ignore.case = TRUE)) %>%
  select(-any_of(c("complete", "date", "Date"))) %>%
  arrange(Year)

# 4. Define Explicit Exclusions
ignored_cols <- c(
  "X10.PredMammalNCC_1_smoothed",
  "X11.PredMammalSmolt_0_smoothed",
  "X06.Cond1NCC_2_smoothed"
)

# Filter out ignored columns
df_lisa_clean <- df_lisa %>% select(-any_of(ignored_cols))
df_doug_clean <- df_doug %>% select(-any_of(ignored_cols))

# 5. Align Shared Set of Shared Target Columns
shared_cols <- intersect(names(df_lisa_clean), names(df_doug_clean))
shared_cols <- setdiff(shared_cols, "Year")

cat("Comparing", length(shared_cols), "shared target columns across years 1998–2021...\n\n")

# Align years and columns strictly
lisa_sub <- df_lisa_clean %>% filter(Year >= 1998 & Year <= 2021) %>% select(Year, all_of(shared_cols))
doug_sub <- df_doug_clean %>% filter(Year >= 1998 & Year <= 2021) %>% select(Year, all_of(shared_cols))

# ------------------------------------------------------------------------------
# TEST A: NA POSITIONAL MATCH
# ------------------------------------------------------------------------------
na_mismatches <- list()

for (col in shared_cols) {
  lisa_na <- is.na(lisa_sub[[col]])
  doug_na <- is.na(doug_sub[[col]])
  
  if (!identical(lisa_na, doug_na)) {
    mismatch_years <- lisa_sub$Year[lisa_na != doug_na]
    na_mismatches[[col]] <- mismatch_years
  }
}

if (length(na_mismatches) == 0) {
  cat("[PASS] NA positions are 100% IDENTICAL across all remaining columns.\n")
} else {
  cat("[FAIL] NA positional mismatches detected in the following columns:\n")
  print(na_mismatches)
}

# ------------------------------------------------------------------------------
# TEST B: NUMERICAL VALUE EQUIVALENCE (Set tolerance here)
# ------------------------------------------------------------------------------
val_mismatches <- tibble()
tol <- 0.9  # Change tolerance in one place

for (col in shared_cols) {
  v_lisa <- lisa_sub[[col]]
  v_doug <- doug_sub[[col]]
  
  # Check only where both are non-NA
  valid_idx <- which(!is.na(v_lisa) & !is.na(v_doug))
  diffs     <- abs(v_lisa[valid_idx] - v_doug[valid_idx])
  
  # Identify indices exceeding tolerance
  exceed_mask <- diffs > tol
  bad_idx     <- valid_idx[exceed_mask]
  
  if (length(bad_idx) > 0) {
    tmp <- tibble(
      Year     = lisa_sub$Year[bad_idx],
      Column   = col,
      Lisa_Val = v_lisa[bad_idx],
      Doug_Val = v_doug[bad_idx],
      Diff     = diffs[exceed_mask]
    )
    val_mismatches <- bind_rows(val_mismatches, tmp)
  }
}

if (nrow(val_mismatches) == 0) {
  cat("[PASS] Value comparison passed: No entries exceed tolerance of", tol, "\n")
} else {
  cat("[FAIL] Found", nrow(val_mismatches), "mismatches exceeding tolerance of", tol, ":\n")
  print(val_mismatches, n = 50)
}


# ==============================================================================
# Script: check_pairwise_correlations.R
# Purpose: Evaluate correlation strength between Lisa and Doug time series
# ==============================================================================

library(tidyverse)

# Calculate Pearson correlation for each column
cor_results <- tibble()

for (col in shared_cols) {
  v_lisa <- lisa_sub[[col]]
  v_doug <- doug_sub[[col]]
  
  # Evaluate on complete pairs only
  valid_idx <- which(!is.na(v_lisa) & !is.na(v_doug))
  
  if (length(valid_idx) > 3) {
    r_val <- cor(v_lisa[valid_idx], v_doug[valid_idx], method = "pearson")
    max_diff <- max(abs(v_lisa[valid_idx] - v_doug[valid_idx]))
    
    cor_results <- bind_rows(cor_results, tibble(
      Column   = col,
      r_cor    = round(r_val, 5),
      Max_Diff = round(max_diff, 4),
      N_Years  = length(valid_idx)
    ))
  }
}

# ------------------------------------------------------------------------------
# REPORT SUMMARY
# ------------------------------------------------------------------------------
cor_results <- cor_results %>% arrange(r_cor)

cat("==============================================================================\n")
cat("PAIRWISE CORRELATION SUMMARY (Across", nrow(cor_results), "Variables)\n")
cat("==============================================================================\n\n")

# Highlight series with r < 0.99
low_cor <- cor_results %>% filter(r_cor < 0.99)

if (nrow(low_cor) == 0) {
  cat("[PASS] ALL variables have Pearson r >= 0.99!\n")
  cat("Conclusion: Discrepancies are purely scalar/additive shifts; relative trajectories are preserved.\n\n")
} else {
  cat("[WARNING]", nrow(low_cor), "variables have Pearson r < 0.99:\n")
  print(low_cor, n = 50)
  cat("\n")
}

# Print top 15 highest absolute difference series with their correlation
cat("Top 15 Series with Largest Numerical Differences:\n")
cor_results %>% arrange(desc(Max_Diff)) %>% head(15) %>% print()


# 3 variables have Pearson r < 0.99:
#   # A tibble: 3 × 4
#   Column                       r_cor Max_Diff N_Years
# 1 X12.ZooPreyAK_1_smoothed     0.974    0.324      24
# 2 X08.PredBirdNCC_b_2_smoothed 0.988    0.236      24
# 3 X10.PredMammalNCC_DFA1       0.990    1.19       24


library(tidyverse)

# Pull year-by-year values for X10.PredMammalNCC_DFA1
ssl_check <- tibble(
  Year = lisa_sub$Year,
  Lisa = lisa_sub[["X10.PredMammalNCC_DFA1"]],
  Doug = doug_sub[["X10.PredMammalNCC_DFA1"]]
) %>%
  mutate(
    Diff = Lisa - Doug,
    Abs_Diff = abs(Diff)
  )

cat("Year-by-Year Comparison for X10.PredMammalNCC_DFA1:\n")
print(ssl_check, n = 24)

# Plot trajectories together
ggplot(ssl_check, aes(x = Year)) +
  geom_line(aes(y = Lisa, color = "Lisa"), size = 1) +
  geom_line(aes(y = Doug, color = "Doug"), size = 1, linetype = "dashed") +
  labs(
    title = "Trajectory Comparison: X10.PredMammalNCC_DFA1",
    subtitle = "Checking for structural divergence vs scalar offset",
    y = "DFA Index Value"
  ) +
  theme_minimal()


#Re-standardize------
# ==============================================================================
# Script: check_standardized_correlations.R
# Purpose: Z-scale columns in Lisa and Doug datasets and re-evaluate differences
# ==============================================================================

library(tidyverse)

# 1. Z-score (standardize) columns independently for both datasets
lisa_z <- lisa_sub %>%
  mutate(across(all_of(shared_cols), ~ as.vector(scale(.x))))

doug_z <- doug_sub %>%
  mutate(across(all_of(shared_cols), ~ as.vector(scale(.x))))

# 2. Evaluate correlations and standardized max differences
z_cor_results <- tibble()

for (col in shared_cols) {
  v_lisa_z <- lisa_z[[col]]
  v_doug_z <- doug_z[[col]]
  
  valid_idx <- which(!is.na(v_lisa_z) & !is.na(v_doug_z))
  
  if (length(valid_idx) > 3) {
    r_val    <- cor(v_lisa_z[valid_idx], v_doug_z[valid_idx], method = "pearson")
    max_diff <- max(abs(v_lisa_z[valid_idx] - v_doug_z[valid_idx]))
    
    z_cor_results <- bind_rows(z_cor_results, tibble(
      Column       = col,
      r_cor        = round(r_val, 5),
      Max_Diff_SD  = round(max_diff, 4),
      N_Years      = length(valid_idx)
    ))
  }
}

# ------------------------------------------------------------------------------
# REPORT RESULTS
# ------------------------------------------------------------------------------
z_cor_results <- z_cor_results %>% arrange(desc(Max_Diff_SD))

cat("==============================================================================\n")
cat("STANDARDIZED (Z-SCORED) DIFFERENCE SUMMARY\n")
cat("==============================================================================\n\n")

# Top 15 series with the largest standardized differences
cat("Top 15 Series with Largest Z-Score Discrepancies:\n")
print(head(z_cor_results, 15))

cat("\n------------------------------------------------------------------------------\n")
cat("X10.PredMammalNCC_DFA1 Standardized Metrics:\n")
print(z_cor_results %>% filter(Column == "X10.PredMammalNCC_DFA1"))
