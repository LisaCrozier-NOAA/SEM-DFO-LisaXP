

# ==============================================================================
# Script: 01_recreate_doug_baseline_and_allsmoothed_1998_2021.R
# Purpose: Recreate Doug's 1998-2021 DFA products + generate All-Smoothed matrix
# Output: copilot/outputs_altprey/
# ==============================================================================

library(tidyverse)
library(lubridate)
library(MARSS)
library(imputeTS)

rootdir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"
path    <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"

dataDir   <- file.path(path, "data")
outputDir <- file.path(rootdir, "copilot/outputs_altprey")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

source(file.path(path, "functions.R"))

screenStartDatetime <- dmy("01JAN1998")
screenEndDatetime   <- dmy("31DEC2021")
target_years        <- data.frame(date = ymd(paste0(1998:2021, "-01-01")))

# 1. Load Manifest & Process Raw Files
indicators <- read.csv(file.path(path, "indicators.csv")) %>%
  filter(category != "") %>%
  filter(incl2026 == "Y")

dataList <- list()
for (subDir in subDirs) {
  subDirPath <- file.path(dataDir, subDir)
  if (!dir.exists(subDirPath)) next
  
  for (dataFile in list.files(subDirPath, pattern = "\\.csv")) {
    thisFile <- file.path(subDirPath, dataFile)
    thisHeader <- readLines(thisFile, n = 4)
    thisIndicator <- trimws(unlist(strsplit(thisHeader[[1]], ","))[2])
    thisIndicator <- gsub('"', "", thisIndicator)
    
    if (thisIndicator %in% indicators$indicator) {
      thisData <- read.csv(thisFile, skip = 4, colClasses = c("character", "numeric"), na.strings = c("null", "NA"))
      
      if (names(thisData)[2] == "Index") {
        thisData$date <- ymd(thisData$Year, truncated = 2L)
        thisData <- thisData[, c("date", "Index")]
      } else if (names(thisData)[1] == "Year") {
        thisData$date <- ymd(thisData$Year, truncated = 2L)
        thisData <- thisData[, c("date", "Value")]
      } else {
        thisData$date <- ymd(thisData$date)
      }
      names(thisData) <- c("date", "value")
      
      thisInd  <- getInd(indicators, thisIndicator, subDir)
      thisData <- impute(thisInd, thisData)
      out      <- logTransform(thisInd, thisData)
      thisData <- out$thisData
      
      thisData$finalVal  <- if (out$transformed) thisData$logTransformed else thisData$imputed
      thisData$shortName <- thisInd$shortName
      
      dataList[[length(dataList) + 1]] <- thisData[, c("shortName", "date", "finalVal")]
    }
  }
}

allData <- bind_rows(dataList)

# 2. Filter 50% Completeness Criterion
yearDF <- data.frame(year = year(seq(screenStartDatetime, screenEndDatetime, by = "1 year")))

qualified_series <- allData %>%
  group_by(shortName) %>%
  group_modify(~ {
    thisScreenDat <- .x %>% filter(date >= screenStartDatetime & date <= screenEndDatetime)
    thisScreenDat$year <- year(thisScreenDat$date)
    thisYearDat <- thisScreenDat %>%
      group_by(year) %>%
      summarize(finalVal = mean(finalVal, na.rm = TRUE), .groups = "drop") %>%
      full_join(yearDF, by = "year")
    
    fracComplete <- sum(is.finite(thisYearDat$finalVal)) / nrow(thisYearDat)
    pass <- (fracComplete >= 0.5) | str_starts(toupper(.y$shortName), "SAR_")
    tibble(pass = pass)
  }) %>%
  filter(pass)

qualified_data <- allData %>% inner_join(qualified_series, by = "shortName")

# 3. Build Qualified Wide Matrix (1998-2021)
datWide_1998_2021 <- target_years %>%
  left_join(qualified_data, by = "date") %>%
  pivot_wider(id_cols = date, names_from = shortName, values_from = finalVal, values_fn = mean) %>%
  mutate(Year = year(date)) %>%
  select(Year, everything(), -date) %>%
  arrange(Year)

write.csv(datWide_1998_2021, file.path(outputDir, "datWide_1998_2021_qualified.csv"), row.names = FALSE)

# 4. Generate All-Smoothed Matrix (Single-Series MARSS Filter for ALL Variables)
all_marss_objects <- list()
allsmoothed_list  <- list()

for (col in setdiff(names(datWide_1998_2021), "Year")) {
  vec_raw <- datWide_1998_2021[[col]]
  valid_idx <- which(!is.na(vec_raw))
  
  if (length(valid_idx) > 3) {
    scaled_vec <- scale(vec_raw) %>% as.vector()
    
    # Fit MARSS univariate state-space model
    fit_single <- MARSS(scaled_vec, fit = FALSE, silent = TRUE)
    fit_single$par <- fit_single$start
    kfOut <- MARSSkf(fit_single)
    
    smoothed_vals <- as.numeric(t(kfOut$xtT))
    
    allsmoothed_list[[col]] <- smoothed_vals
    all_marss_objects[[paste0(col, "_smoothed")]] <- fit_single
  } else {
    allsmoothed_list[[col]] <- vec_raw
  }
}

datWide_1998_2021_allsmoothed <- as_tibble(allsmoothed_list) %>%
  mutate(Year = datWide_1998_2021$Year) %>%
  select(Year, everything())

# Export Containers
write.csv(datWide_1998_2021_allsmoothed, file.path(outputDir, "datWide_1998_2021_allsmoothed.csv"), row.names = FALSE)
saveRDS(all_marss_objects, file.path(outputDir, "all_marss_objects_1998_2021.rds"))

cat("\nPhase 1 Complete!\nExported: datWide_1998_2021_qualified.csv & datWide_1998_2021_allsmoothed.csv\n")