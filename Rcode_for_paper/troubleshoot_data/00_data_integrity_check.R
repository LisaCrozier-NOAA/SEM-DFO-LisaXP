# ==============================================================================
# Script: 00_data_integrity_check.R
# Directory: Rcode_for_paper/01_data_prep/
# Purpose: Detect corrupted, flatlined (zero-variance / straight line) time series
#          across all modeling indicators.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 1: LOAD RAW SEM DATA
# ------------------------------------------------------------------------------
proj_dir      <- file.path(getwd(), "Rcode_for_paper")
sem_data_file <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")

if (!file.exists(sem_data_file)) {
  stop("Cannot find data file at: ", sem_data_file)
}

raw_sem_data <- read.csv(sem_data_file, stringsAsFactors = FALSE)

cat(sprintf("Loaded dataset with %d rows (years) and %d columns.\n", nrow(raw_sem_data), ncol(raw_sem_data)))

# ------------------------------------------------------------------------------
# STEP 2: ISOLATE MODELING COLUMNS (FOCUS ON SMOLTYR / NON-DUPLICATED)
# ------------------------------------------------------------------------------
all_cols <- names(raw_sem_data)

# Extract modeling columns (X01 through X21)
model_cols <- all_cols[grepl("^X[0-9]{2}", all_cols)]

# Exclude 'adultyr' versions to keep the audit clean (focus on smoltyr and unsuffixed)
audit_cols <- model_cols[!grepl("adultyr$", model_cols)]


par(mfrow=c(5,5))
for (i in 1:length(audit_cols)){
  plot(raw_sem_data$Year,raw_sem_data[,audit_cols[i]],main=audit_cols[i])
}

cat(sprintf("Auditing %d primary modeling indicators...\n\n", length(audit_cols)))

# ------------------------------------------------------------------------------
# STEP 3: SCAN FOR ZERO VARIANCE / STRAIGHT LINES / CONSTANTS
# ------------------------------------------------------------------------------
corrupted_series <- list()
valid_series     <- list()

for (col in audit_cols) {
  vals <- raw_sem_data[[col]]
  
  # Remove NAs for variance check
  clean_vals <- vals[!is.na(vals)]
  
  if (length(clean_vals) == 0) {
    status <- "ALL NA"
    sd_val <- NA
    val_range <- NA
  } else {
    sd_val    <- sd(clean_vals)
    val_range <- max(clean_vals) - min(clean_vals)
    
    # Check if straight line (SD == 0 or Range == 0 or unique values <= 1)
    is_straight_line <- (sd_val == 0) || (val_range == 0) || (length(unique(clean_vals)) <= 1)
    
    status <- if (is_straight_line) "CORRUPTED (FLATLINE)" else "OK"
  }
  
  res_row <- tibble(
    column_name    = col,
    status         = status,
    n_unique       = length(unique(clean_vals)),
    std_dev        = round(sd_val, 4),
    val_range      = round(val_range, 4),
    mean_val       = round(mean(clean_vals, na.rm = TRUE), 4),
    matching_adult = if (paste0(col, "_adultyr") %in% all_cols) paste0(col, "_adultyr") else 
      if (gsub("smoltyr$", "adultyr", col) %in% all_cols) gsub("smoltyr$", "adultyr", col) else "None"
  )
  
  if (status != "OK") {
    corrupted_series[[length(corrupted_series) + 1]] <- res_row
  } else {
    valid_series[[length(valid_series) + 1]] <- res_row
  }
}

corrupted_df <- bind_rows(corrupted_series)
valid_df     <- bind_rows(valid_series)

# ------------------------------------------------------------------------------
# STEP 4: REPORT RESULTS
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("                       DATA INTEGRITY AUDIT RESULTS                           \n")
cat("==============================================================================\n\n")

if (nrow(corrupted_df) > 0) {
  cat(sprintf("⚠️ CRITICAL ALERT: Found %d CORRUPTED / FLATLINE series!\n\n", nrow(corrupted_df)))
  print(corrupted_df, n = 100)
  
  # Save corrupted column list for reference or automated dropping
  write_csv(corrupted_df, file.path(proj_dir, "metadata/CORRUPTED_FLATLINE_SERIES.csv"))
  cat("\nList of corrupted series saved to: metadata/CORRUPTED_FLATLINE_SERIES.csv\n")
} else {
  cat("✅ ALL CLEAR: No constant or zero-variance time series detected across audit columns!\n")
}

cat("\nSummary: ", nrow(valid_df), "valid series | ", nrow(corrupted_df), "corrupted series.\n")
