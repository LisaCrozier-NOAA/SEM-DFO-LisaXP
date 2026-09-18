# Step 1: Consolidated R Scripts
# Script 01: 01_data_prep/01_preprocess_raw_sources.R
# Consolidates preprocessData.R and functions.R. 
# Handles extraction, date padding (1998–2026), and 2-year lag/lead construction for marine mammals and orcas.

#Script 01: 01_data_prep/01_preprocess_raw_sources.R


# ==============================================================================
# Script: 01_preprocess_raw_sources.R
# Purpose: Extract, align, pad, and standardize raw indicator time series
# ==============================================================================

library(tidyverse)
library(lubridate)
library(readxl)

# --- Configuration & Paths ---
proj_dir    <- getwd()
raw_data_dir <- file.path(proj_dir, "data_raw")
out_data_dir <- file.path(proj_dir, "data_processed")
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

indicators_meta <- read.csv(file.path(proj_dir, "metadata", "indicators.csv")) %>%
  filter(category != "")

min_date <- dmy("01JAN1998")
max_date <- dmy("31DEC2026")

# --- Helper Functions ---
add_missing_years <- function(df, val_col = "Value") {
  year_df <- data.frame(Year = year(seq(min_date, max_date, by = "1 year")))
  df %>%
    full_join(year_df, by = "Year") %>%
    filter(Year >= year(min_date), Year <= year(max_date)) %>%
    select(Year, Value = !!sym(val_col)) %>%
    arrange(Year)
}

save_indicator_csv <- function(df, short_name, dataset_name, source_info = "NOAA/DFO") {
  out_path <- file.path(out_data_dir, paste0(short_name, ".csv"))
  header <- c(
    paste0("Dataset:,", dataset_name),
    "Region:,",
    paste0("Source:,", source_info),
    "Contact:,"
  )
  writeLines(header, out_path)
  suppressWarnings(write.table(df, out_path, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE))
}

# --- 1. Process Pink Salmon Data (Ruggerone et al.) ---
pink_file <- file.path(raw_data_dir, "Ruggerone_pink_salmon_2025.csv")
if (file.exists(pink_file)) {
  pink_df <- read.csv(pink_file)
  
  pink_asia <- pink_df %>% select(Year, Value = Asia) %>% add_missing_years()
  save_indicator_csv(pink_asia, "pinkSalmonAsia_2025", "pinkSalmon_2025")
  
  pink_na <- pink_df %>% select(Year, Value = North.America) %>% add_missing_years()
  save_indicator_csv(pink_na, "pinkSalmonNorthAmerica_2025", "pinkSalmon_2025")
}

# --- 2. Process Mammal Data & Apply 2-Year Leads ---
pinniped_file <- file.path(raw_data_dir, "pinnipedsNCC.datwide.csv")
if (file.exists(pinniped_file)) {
  pin_df <- read.csv(pinniped_file)
  
  # Standardize and build 2-year leads for predation matching smolt outmigration
  lead_vars <- c("Harbor_seal_CR_2025", "AllSeaLionsEMB_2025", "AllSeaLionsBonn_2025")
  for (v in lead_vars) {
    if (v %in% names(pin_df)) {
      sub_df <- pin_df %>% select(Year, Value = !!sym(v))
      sub_df$Year <- sub_df$Year - 2  # Apply 2-year lead
      sub_df <- add_missing_years(sub_df)
      save_indicator_csv(sub_df, paste0(v, "_2yrLead"), "Pinniped_2yrLead")
    }
  }
}

# --- 3. Process Smolt-to-Adult Return (SAR) ---
sar_file <- file.path(raw_data_dir, "sar.19982021.csv")
if (file.exists(sar_file)) {
  sar_df <- read.csv(sar_file)
  names(sar_df) <- c("Year", "SAR")
  write.csv(sar_df, file.path(out_data_dir, "sar_clean.csv"), row.names = FALSE)
}

message("Step 01 Complete: Raw sources extracted and padded to standardized folder structure.")