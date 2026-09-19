


# ==============================================================================
# Script: 03_append_all_altprey_and_supp_indicators.R
# Purpose: Extract standalone AltPrey, habCompInd temperature vector, and 
#          supplementary ak_yr/eulachon indicators (with _smoltyr and _adultyr 
#          2-year cohort leads), merge into master LisaName matrix, re-sort 
#          alphanumerically, and export 1998-2021 SEM modeling subsets.
# Outputs: 
#   - metadata/clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv
#   - metadata/sem_data_1998_2021.csv
#   - metadata/sem_altprey_data_complete_1998_2021.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(MARSS)

# --- Configuration & Paths ---
rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
outputDir <- file.path(rootdir, "metadata")

master_file  <- file.path(outputDir, "clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv")
raw_file     <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw_shortnames.csv")

if (!file.exists(raw_file)) {
  raw_file <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw.csv")
}

eulachon_file <- "copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv"
if (!file.exists(eulachon_file)) {
  eulachon_file <- file.path(dirname(rootdir), "copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv")
}

ak_yr_file <- "data_Lisa/ak_yr.csv"
if (!file.exists(ak_yr_file)) {
  ak_yr_file <- file.path(dirname(rootdir), "data_Lisa/ak_yr.csv")
}

# Pre-flight checks
if (!file.exists(master_file))   stop("Missing master LisaName file! Run Script 02 first.")
if (!file.exists(raw_file))      stop("Missing raw input file at: ", raw_file)
if (!file.exists(eulachon_file)) stop("Missing eulachon input file at: ", eulachon_file)
if (!file.exists(ak_yr_file))   stop("Missing ak_yr input file at: ", ak_yr_file)

master_df <- read.csv(master_file, stringsAsFactors = FALSE)
raw_df    <- read.csv(raw_file, stringsAsFactors = FALSE)

baseline_years <- 1998:2021
full_years     <- raw_df$Year

scale_baseline_lock <- function(vec, yrs) {
  base_idx  <- which(yrs %in% baseline_years)
  base_mean <- mean(vec[base_idx], na.rm = TRUE)
  base_sd   <- sd(vec[base_idx], na.rm = TRUE)
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

# ------------------------------------------------------------------------------
# 1. PROCESS RAW STRAGGLERS & ALTPREY (habCompInd, Pollock Age-1+, Sitka Herring)
# ------------------------------------------------------------------------------
cat("--- STEP 1: Processing Raw Stragglers & AltPrey Constituents ---\n")

# Single-indicator MARSS smoother helper
smooth_raw_straggler <- function(col_name) {
  matched_col <- names(raw_df)[names(raw_df) == col_name | 
                                 grep(gsub("_2025|_2026", "", col_name), names(raw_df), ignore.case = TRUE)][1]
  if (is.na(matched_col)) stop("Missing column in raw_df: ", col_name)
  
  vec_raw <- raw_df[[matched_col]]
  nas_raw <- raw_df$Year[is.na(vec_raw)]
  
  scaled_ext  <- scale_baseline_lock(vec_raw, full_years)
  scaled_base <- scaled_ext[full_years %in% baseline_years]
  
  fit_single <- MARSS(scaled_base, fit = FALSE, silent = TRUE)
  fit_single$par <- fit_single$start
  
  ext_proj <- MARSS(scaled_ext, fit = FALSE, silent = TRUE)
  ext_proj$par <- fit_single$par
  
  smoothed_trend <- as.numeric(t(MARSSkf(ext_proj)$xtT))
  if (length(nas_raw) > 0) smoothed_trend[full_years %in% nas_raw] <- NA_real_
  return(smoothed_trend)
}

# Smooth and generate _smoltyr and _adultyr 2-year leads
raw_stragglers_df <- tibble(
  Year = full_years,
  X01_habCompInd_smoltyr        = smooth_raw_straggler("habCompInd"),
  X13_pollock_age1plus_smoltyr  = smooth_raw_straggler("pollockBiomassAIage1plus_predAK_2026"),
  X13_sitkaHerring_EGoA_smoltyr = smooth_raw_straggler("sitkaHerring_EGoA")
) %>%
  mutate(
    X01_habCompInd_adultyr        = dplyr::lead(X01_habCompInd_smoltyr, 2),
    X13_pollock_age1plus_adultyr  = dplyr::lead(X13_pollock_age1plus_smoltyr, 2),
    X13_sitkaHerring_EGoA_adultyr = dplyr::lead(X13_sitkaHerring_EGoA_smoltyr, 2)
  )

# ------------------------------------------------------------------------------
# 2. PROCESS SUPPLEMENTARY INDICATORS (Eulachon & ak_yr.csv with Leads)
# ------------------------------------------------------------------------------
cat("--- STEP 2: Processing Eulachon & ak_yr Indicators (with Leads) ---\n")

eulachon_scaled <- read.csv(eulachon_file, stringsAsFactors = FALSE) %>%
  select(Year, eulachon_during_chinook, eulachon_annual_scaled) %>%
  rename(
    X05_eulachon_during_chinook_smoltyr = eulachon_during_chinook,
    X05_eulachon_annual_scaled_smoltyr   = eulachon_annual_scaled
  ) %>%
  mutate(
    X05_eulachon_during_chinook_adultyr = dplyr::lead(X05_eulachon_during_chinook_smoltyr, 2),
    X05_eulachon_annual_scaled_adultyr   = dplyr::lead(X05_eulachon_annual_scaled_smoltyr, 2)
  )

ak_yr_selected <- read.csv(ak_yr_file, row.names = NULL, stringsAsFactors = FALSE) %>%
  rename(Year = year) %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  select(
    Year,
    X11_ssl_seak_pup_pred        = ssl_seak_pup_pred,
    X12_egoa_krill_smoltyr       = secm_euph_dens,
    X13_mid_il_capelin_smoltyr   = mid_il_capelin,
    X21_sst_egoa_junjulaug_smoltyr = sst_egoa_coastwatch_junjulaug
  ) %>%
  mutate(
    X12_egoa_krill_adultyr         = dplyr::lead(X12_egoa_krill_smoltyr, 2),
    X13_mid_il_capelin_adultyr     = dplyr::lead(X13_mid_il_capelin_smoltyr, 2),
    X21_sst_egoa_junjulaug_adultyr = dplyr::lead(X21_sst_egoa_junjulaug_smoltyr, 2)
  )

# ------------------------------------------------------------------------------
# 3. MERGE ALL SUPPLEMENTARY INDICATORS & ENFORCE ALPHANUMERIC ORDER
# ------------------------------------------------------------------------------
cat("--- STEP 3: Merging & Sorting Master Matrix ---\n")

all_new_cols <- c(
  setdiff(names(raw_stragglers_df), "Year"),
  setdiff(names(eulachon_scaled), "Year"),
  setdiff(names(ak_yr_selected), "Year")
)

updated_master <- master_df %>%
  select(-any_of(all_new_cols)) %>%
  left_join(raw_stragglers_df, by = "Year") %>%
  left_join(eulachon_scaled,    by = "Year") %>%
  left_join(ak_yr_selected,     by = "Year") %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  arrange(Year)

# Enforce strict alphanumeric column ordering (Year, X01_..., X02_..., ..., X21_...)
sorted_cols    <- c("Year", sort(setdiff(names(updated_master), "Year")))
updated_master <- updated_master %>% select(all_of(sorted_cols))

# Export updated master matrix (1996-2025)
write.csv(updated_master, master_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. GENERATE 1998–2021 SEM MODELING SUBSETS
# ------------------------------------------------------------------------------
cat("--- STEP 4: Exporting 1998–2021 SEM Subsets ---\n")

df_1998_2021 <- updated_master %>% filter(Year >= 1998 & Year <= 2021)

valid_data_cols <- names(df_1998_2021)[colSums(!is.na(df_1998_2021)) > 0]

sem_complete_data <- df_1998_2021 %>%
  select(all_of(valid_data_cols)) %>%
  mutate(across(-Year, ~ as.vector(scale(.x))))

out_complete <- file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv")
out_raw_sem  <- file.path(outputDir, "sem_data_1998_2021.csv")

write.csv(sem_complete_data, out_complete, row.names = FALSE)
write.csv(df_1998_2021,       out_raw_sem,  row.names = FALSE)

# ------------------------------------------------------------------------------
# VERIFICATION SUMMARY
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("PROCESS COMPLETE! CONSOLIDATED SCRIPT 03 EXECUTED SUCCESSFULLY.\n")
cat("==============================================================================\n")
cat("1. Extended Master Matrix (1996-2025):\n   ", master_file, "\n")
cat("2. Unscaled SEM Dataset (1998-2021):\n   ", out_raw_sem, "\n")
cat("3. Scaled SEM Dataset (1998-2021):\n   ", out_complete, "\n")
cat(sprintf("\nTotal Columns in Final Master Matrix: %d\n", ncol(updated_master)))
cat("\nNewly Appended Indicator Series:\n")
print(all_new_cols)
cat("\nFirst 10 Columns in Alphanumeric Order:\n")
print(names(updated_master)[1:10])
cat("==============================================================================\n")

