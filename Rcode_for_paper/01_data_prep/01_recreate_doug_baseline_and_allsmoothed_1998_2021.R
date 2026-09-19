


# ==============================================================================
# Script: 01_reprocess_raw_data_scratch_1995_2025.R
# Purpose: Build un-truncated 1995-2025 raw data matrix directly from source files
#          without using Doug's legacy functions.R
# Output: metadata/datWide_1995_2025_reprocessed_raw.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(imputeTS)

# ------------------------------------------------------------------------------
# 1. SETUP UPDATED PATHS AND DIRECTORIES
# ------------------------------------------------------------------------------
rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
path      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"

dataDir   <- file.path(path, "data")
outputDir <- file.path(rootdir, "metadata")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# Target full timeline (1995 to 2025)
target_years <- data.frame(Year = 1995:2025)

# Subdirectories containing the source data files
subDirs <- c("Western_Aleutian_Islands", "Central_Aleutian_Islands", "Eastern_Aleutian_Islands",
             "Western_Gulf_of_Alaska", "Eastern_Gulf_of_Alaska", "CCIEA", "WCVI", "SEM_data_2024",
             "pinkSalmon", "BonnPinn", "ColumbiaRiverPinn", "SEM_data_2025", "pinkSalmon_2025",
             "data_2025", "NHL_2025", "planktonJuneNCC_2025", "plankton_2025b",
             "plankton_2026", "NHL_2026", "JSOES_2026", "Westport_2026", "WGOA_DFA_2026", "lingcodSA_2026",
             "predAK_2026", "orca_2026")

# ------------------------------------------------------------------------------
# 2. LOAD INDICATOR MANIFEST (indicators.csv)
# ------------------------------------------------------------------------------
indicators <- read.csv(file.path(path, "indicators.csv"), stringsAsFactors = FALSE) %>%
  filter(category != "", incl2026 == "Y")

# ------------------------------------------------------------------------------
# 3. CUSTOM IN-LINE INTERPOLATION & LOG-TRANSFORM FUNCTIONS
# ------------------------------------------------------------------------------
# Interpolate interior NAs ONLY; preserve leading/trailing NAs outside observed range
clean_interpolate <- function(vec, impute_type) {
  if (impute_type == "none" || all(is.na(vec))) return(vec)
  
  # Set non-positive values to NA if linear non-zero constraint is required
  if (impute_type == "linear") {
    vec[vec <= 0] <- NA
  } else if (impute_type == "linearKeepZeros") {
    vec[vec < 0]  <- NA
  }
  
  valid_idx <- which(!is.na(vec))
  if (length(valid_idx) < 2) return(vec)
  
  first_v <- min(valid_idx)
  last_v  <- max(valid_idx)
  
  # Linearly interpolate ONLY between the first and last valid data points
  interior <- vec[first_v:last_v]
  if (any(is.na(interior))) {
    vec[first_v:last_v] <- imputeTS::na_interpolation(interior, option = "linear")
  }
  return(vec)
}

# Apply simple natural log transformation if specified
clean_log_transform <- function(vec, log_setting) {
  if (log_setting == "simple") {
    return(log(vec))
  }
  return(vec)
}

# ------------------------------------------------------------------------------
# 4. INGEST AND PROCESS RAW DATA
# ------------------------------------------------------------------------------
dataList <- list()

for (subDir in subDirs) {
  subDirPath <- file.path(dataDir, subDir)
  if (!dir.exists(subDirPath)) next
  
  for (dataFile in list.files(subDirPath, pattern = "\\.csv$")) {
    thisFile <- file.path(subDirPath, dataFile)
    
    # Extract indicator name from header line 1
    thisHeader    <- readLines(thisFile, n = 1)
    thisIndicator <- trimws(unlist(strsplit(thisHeader, ","))[2])
    thisIndicator <- gsub('"', "", thisIndicator)
    
    # Match with indicators manifest
    ind_meta <- indicators %>% filter(indicator == thisIndicator, dataset == subDir)
    
    if (nrow(ind_meta) > 0) {
      raw_df <- read.csv(thisFile, skip = 4, colClasses = c("character", "numeric"), na.strings = c("null", "NA"))
      
      if (nrow(raw_df) == 0) next
      
      # Standardize Year column
      if ("Year" %in% names(raw_df)) {
        raw_df$Year <- suppressWarnings(as.numeric(raw_df$Year))
      } else if ("date" %in% names(raw_df)) {
        raw_df$Year <- year(ymd(raw_df$date))
      } else {
        raw_df$Year <- suppressWarnings(as.numeric(raw_df[[1]]))
      }
      
      val_col <- names(raw_df)[names(raw_df) %in% c("Value", "Index", "value")][1]
      if (is.na(val_col)) val_col <- names(raw_df)[2]
      
      df_clean <- raw_df %>%
        select(Year, Value = !!sym(val_col)) %>%
        filter(!is.na(Year), Year >= 1995, Year <= 2025) %>%
        group_by(Year) %>%
        summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop")
      
      # Merge onto 1995-2025 grid
      full_grid <- target_years %>% left_join(df_clean, by = "Year")
      
      # Apply interpolation and transformation
      impute_setting <- ind_meta$impute[1]
      log_setting    <- ind_meta$logTransform[1]
      
      vec_interp <- clean_interpolate(full_grid$Value, impute_setting)
      vec_final  <- clean_log_transform(vec_interp, log_setting)
      
      dataList[[length(dataList) + 1]] <- tibble(
        Year      = full_grid$Year,
        shortName = ind_meta$shortName[1],
        finalVal  = vec_final
      )
    }
  }
}

# ------------------------------------------------------------------------------
# 5. PIVOT WIDE AND EXPORT REPROCESSED MATRIX
# ------------------------------------------------------------------------------
allData_long <- bind_rows(dataList)

datWide_raw_full <- target_years %>%
  left_join(allData_long, by = "Year") %>%
  group_by(Year, shortName) %>%
  summarize(finalVal = mean(finalVal, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(shortName)) %>%
  pivot_wider(id_cols = Year, names_from = shortName, values_from = finalVal) %>%
  arrange(Year)

# Export complete un-truncated raw matrix
out_file <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw.csv")
write.csv(datWide_raw_full, out_file, row.names = FALSE)

cat("\n==============================================================================\n")
cat("RAW REPROCESSING COMPLETE (1995-2025)\n")
cat("Retained", ncol(datWide_raw_full) - 1, "indicators across 1995-2025.\n")
cat("Saved file to:", out_file, "\n")
cat("==============================================================================\n")

