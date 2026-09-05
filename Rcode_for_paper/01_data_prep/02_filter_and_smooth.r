#Script 02: 01_data_prep/02_filter_and_smooth.R
#Consolidates analyzeIndicators.R, assessIndicators.R, functions.R, and indicator_exclude_draft.r. 
#Performs completeness filtering, NA linear interpolation, log transformations, and MARSS state-space Kalman smoothing for all time series.
#Exports allData.rds and datWide.csv

# ==============================================================================
# Script: 02_filter_and_smooth.R
# Purpose: QA/QC screening, log transformation, interpolation, and MARSS smoothing
# ==============================================================================

library(tidyverse)
library(lubridate)
library(imputeTS)
library(MARSS)

proj_dir     <- getwd()
data_dir     <- file.path(proj_dir, "data_processed")
output_dir   <- file.path(proj_dir, "output", "data_prep")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

indicators_meta <- read.csv(file.path(proj_dir, "metadata", "indicators.csv")) %>%
  filter(category != "", include == "Y")

min_year <- 1998
max_year <- 2021
min_frac_complete <- 0.50

# --- Filtering and Imputation Functions ---
impute_series <- function(vec, rule = "linear") {
  if (all(is.na(vec))) return(vec)
  if (rule == "linear") vec[vec <= 0] <- NA
  finite_idx <- which(is.finite(vec))
  if (length(finite_idx) < 2) return(vec)
  
  imp <- na_interpolation(vec)
  # Mask out extrapolated ends
  imp[1:max(1, finite_idx[1] - 1)] <- NA
  imp[min(max(finite_idx) + 1, length(imp)):length(imp)] <- NA
  return(imp)
}

smooth_marss <- function(vec) {
  scaled <- scale(vec) %>% as.vector()
  not_na <- which(!is.na(scaled))
  if (length(not_na) < 3) return(scaled)
  
  first_idx <- not_na[1]
  last_idx  <- not_na[length(not_na)]
  
  fit <- MARSS(scaled[first_idx:last_idx], fit = FALSE, silent = TRUE)
  fit$par <- fit$start
  kf <- MARSSkf(fit)
  
  res <- scaled
  res[first_idx:last_idx] <- as.numeric(t(kf$xtT))
  return(res)
}

# --- Main Iteration Loop ---
processed_list <- list()
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  short_name <- gsub(".csv", "", basename(f))
  meta <- indicators_meta %>% filter(shortName == short_name)
  if (nrow(meta) == 0) next
  
  raw_dat <- read.csv(f, skip = 4, header = TRUE, col.names = c("date", "value"), na.strings = c("NA", "null"))
  raw_dat$date <- ymd(raw_dat$date, truncated = 2L)
  raw_dat <- raw_dat %>% filter(year(date) >= min_year, year(date) <= max_year)
  
  # Completeness check
  frac_complete <- sum(!is.na(raw_dat$value)) / nrow(raw_dat)
  if (frac_complete < min_frac_complete && !startsWith(toupper(short_name), "SAR")) next
  
  # Impute & Transform
  raw_dat$imputed <- impute_series(raw_dat$value, meta$impute[1])
  if (!is.na(meta$logTransform[1]) && meta$logTransform[1] == "simple") {
    raw_dat$finalVal <- log(raw_dat$imputed)
  } else {
    raw_dat$finalVal <- raw_dat$imputed
  }
  
  # MARSS Smoothing
  raw_dat$smoothed  <- smooth_marss(raw_dat$finalVal)
  raw_dat$shortName <- short_name
  raw_dat$SEMlatent <- meta$SEMlatent[1]
  raw_dat$guild     <- meta$guild[1]
  
  processed_list[[short_name]] <- raw_dat
}

all_data <- bind_rows(processed_list)
saveRDS(all_data, file.path(output_dir, "allData.rds"))

# Output wide matrix for downstream validation
wide_data <- all_data %>%
  pivot_wider(id_cols = date, names_from = shortName, values_from = finalVal) %>%
  arrange(date)

write.csv(wide_data, file.path(output_dir, "datWide.csv"), row.names = FALSE)
message("Step 02 Complete: Indicators filtered, imputed, smoothed, and saved to allData.rds.")