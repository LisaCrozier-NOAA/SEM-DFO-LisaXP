

# ==============================================================================
# Script: 03a_extract_and_verify_altprey_constituents_FIXED.R
# Purpose: Extract exact pollock age1+ and Sitka herring columns from raw_df,
#          verify distinct values, run MARSS Kalman smoothing, and generate 2-year leads.
# ==============================================================================

library(tidyverse)
library(lubridate)
library(MARSS)

rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
outputDir <- file.path(rootdir, "metadata")

raw_file <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw_shortnames.csv")
if (!file.exists(raw_file)) raw_file <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw.csv")

raw_df <- read.csv(raw_file, stringsAsFactors = FALSE)

baseline_years <- 1998:2021
full_years     <- raw_df$Year

scale_baseline_lock <- function(vec, yrs) {
  base_idx  <- which(yrs %in% baseline_years)
  base_mean <- mean(vec[base_idx], na.rm = TRUE)
  base_sd   <- sd(vec[base_idx], na.rm = TRUE)
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

# --- EXACT MATCHES FROM RAW DATA ---
raw_pollock_col <- "pollockBiomassAIage1plus_predAK_2026"
raw_herring_col <- "sitkaHerring_EGoA"

cat(sprintf("Extracting Pollock Column: '%s'\n", raw_pollock_col))
cat(sprintf("Extracting Herring Column: '%s'\n\n", raw_herring_col))

# --- 1. PROCESS POLLOCK AGE 1+ ---
vec_pol_raw <- raw_df[[raw_pollock_col]]
nas_pol     <- raw_df$Year[is.na(vec_pol_raw)]

scaled_pol_ext  <- scale_baseline_lock(vec_pol_raw, full_years)
scaled_pol_base <- scaled_pol_ext[full_years %in% baseline_years]

fit_pol <- MARSS(scaled_pol_base, fit = FALSE, silent = TRUE)
fit_pol$par <- fit_pol$start
proj_pol <- MARSS(scaled_pol_ext, fit = FALSE, silent = TRUE)
proj_pol$par <- fit_pol$par
smooth_pol <- as.numeric(t(MARSSkf(proj_pol)$xtT))

# Preserve true boundary NAs from raw data
if (length(nas_pol) > 0) smooth_pol[full_years %in% nas_pol] <- NA_real_

# --- 2. PROCESS SITKA HERRING EGOA ---
vec_her_raw <- raw_df[[raw_herring_col]]
nas_her     <- raw_df$Year[is.na(vec_her_raw)]

scaled_her_ext  <- scale_baseline_lock(vec_her_raw, full_years)
scaled_her_base <- scaled_her_ext[full_years %in% baseline_years]

fit_her <- MARSS(scaled_her_base, fit = FALSE, silent = TRUE)
fit_her$par <- fit_her$start
proj_her <- MARSS(scaled_her_ext, fit = FALSE, silent = TRUE)
proj_her$par <- fit_her$par
smooth_her <- as.numeric(t(MARSSkf(proj_her)$xtT))

# Preserve true boundary NAs from raw data
if (length(nas_her) > 0) smooth_her[full_years %in% nas_her] <- NA_real_

# --- 3. COMBINE, LEAD(2), AND BUILD VERIFIED MATRIX ---
altprey_df <- tibble(
  Year = full_years,
  X13_pollock_age1plus_smoltyr  = smooth_pol,
  X13_sitkaHerring_EGoA_smoltyr = smooth_her
) %>%
  mutate(
    X13_pollock_age1plus_adultyr  = dplyr::lead(X13_pollock_age1plus_smoltyr, 2),
    X13_sitkaHerring_EGoA_adultyr = dplyr::lead(X13_sitkaHerring_EGoA_smoltyr, 2)
  ) %>%
  filter(Year >= 1996 & Year <= 2025)

# --- 4. VERIFICATION DIAGNOSTIC ---
are_identical <- identical(altprey_df$X13_pollock_age1plus_smoltyr, altprey_df$X13_sitkaHerring_EGoA_smoltyr)

cat("==============================================================================\n")
cat("VERIFICATION DIAGNOSTIC\n")
cat("==============================================================================\n")
if (are_identical) {
  stop("CRITICAL ERROR: Series are still identical! Check script execution.")
} else {
  cat("SUCCESS: Both standalone series successfully extracted and smoothed!\n\n")
  print(altprey_df, n = 15)
}

# --- 5. EXPORT STANDALONE MATRIX FOR STEP 03 MERGE ---
out_altprey_file <- file.path(outputDir, "standalone_altprey_constituents_1996_2025.csv")
write.csv(altprey_df, out_altprey_file, row.names = FALSE)
cat(sprintf("\nSaved standalone AltPrey matrix to: %s\n", out_altprey_file))