
#"C:\Users\Lisa.Crozier\Documents\Marine survival\SEM-DFO-LisaXP\copilot\rcode_alt_prey\doug_01_analyzeIndicators_extend_to_1996_2025_datwide.R"

#This script reads in data files in Doug's data folder that were previously processed in "preprocessData.R"
#it uses his screening criteria, except that I have gone back to the original requirement for minFracComplete <- 0.79, which excludes harbor seals from the CR
#I also extended the time frame 1996:2025 for the resulting data file to be used with 2yrlead requirements for alternate prey
#it creates a wide format data frame "copilot/outputs_8/datWide_1996_2025_qualified.csv"


library(tidyverse)
library(lubridate)
library(imputeTS)

# ---------------------------------------------------------------------------
# Setup & Directories
# ---------------------------------------------------------------------------
rootdir<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"

# Doug directory
# Change index to 1 for Lisa's path
path <- c("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results", "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/")[1]

workingDir <- file.path(path, "analyzeAKindices")


indicatorsFile <- file.path(workingDir, "indicators.csv")
dataDir <- file.path(workingDir, "data")
outputDir <- file.path(rootdir, "copilot/outputs_8")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# Source constants and functions
source(file.path(workingDir, "functions.R"))

scen <- "incl2026"
minFracComplete <- 0.79
screenStartDatetime <- dmy(screenStartDate) # 01JAN1998
screenEndDatetime <- dmy(screenEndDate)     # 31DEC2021

# Define full target output window (1996 to 2025)
target_years <- data.frame(date = ymd(paste0(1996:2025, "-01-01")))



# ---------------------------------------------------------------------------
# Load Indicators Manifest
# ---------------------------------------------------------------------------
indicators <- read.csv(indicatorsFile) %>% 
  filter(category != "")

indicators$include <- indicators[, scen]
indicators <- indicators %>% filter(include == "Y")

# ---------------------------------------------------------------------------
# Load and Process All Indicator Data
# ---------------------------------------------------------------------------
dataList <- list()

for (subDir in subDirs) {
  subDirPath <- file.path(dataDir, subDir)
  if (!dir.exists(subDirPath)) next
  
  thisDataFiles <- list.files(subDirPath, pattern = "\\.csv")
  
  for (dataFile in thisDataFiles) {
    thisFile <- file.path(subDirPath, dataFile)
    thisHeader <- readLines(thisFile, n = 4)
    thisIndicator <- trimws(unlist(strsplit(thisHeader[[1]], ","))[2])
    thisIndicator <- gsub('"', "", thisIndicator)
    
    if (thisIndicator %in% indicators$indicator) {
      cat("Processing:", dataFile, "-", thisIndicator, "\n")
      
      thisData <- read.csv(thisFile, skip = 4, colClasses = c("character", "numeric"), na.strings = c("null", "NA"))
      
      # Standardize date format
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
      
      # Imputation and Log Transformation
      thisInd <- getInd(indicators, thisIndicator, subDir)
      thisData <- impute(thisInd, thisData)
      out <- logTransform(thisInd, thisData)
      thisData <- out$thisData
      
      # Store metadata and finalVal
      thisData$finalVal <- if (out$transformed) thisData$logTransformed else thisData$imputed
      thisData$shortName <- thisInd$shortName
      thisData$indicator <- thisIndicator
      thisData$subDir <- subDir
      
      dataList[[length(dataList) + 1]] <- thisData[, c("indicator", "subDir", "shortName", "date", "finalVal")]
    }
  }
}

allData <- bind_rows(dataList)

# ---------------------------------------------------------------------------
# Screen Indicators (Grouped by indicator, subDir, and shortName)
# ---------------------------------------------------------------------------
yearDF <- data.frame(year = year(seq(screenStartDatetime, screenEndDatetime, by = "1 year")), ID = 1)

qualified_series <- allData %>% 
  group_by(indicator, subDir, shortName) %>% 
  group_modify(~ {
    thisDat <- .x
    thisScreenDat <- thisDat %>% filter(date >= screenStartDatetime & date <= screenEndDatetime)
    thisScreenDat$year <- year(thisScreenDat$date)
    
    thisYearDat <- thisScreenDat %>% 
      group_by(year) %>% 
      summarize(finalVal = mean(finalVal, na.rm = TRUE), .groups = "drop") %>% 
      full_join(yearDF, by = "year")
    
    fracComplete <- sum(is.finite(thisYearDat$finalVal)) / nrow(thisYearDat)
    passFracComplete <- fracComplete >= minFracComplete
    
    # Bypass for SAR variables
    if (str_starts(toupper(.y$shortName), "SAR_")) {
      passFracComplete <- TRUE
    }
    
    tibble(fracComplete = fracComplete, pass = passFracComplete)
  }) %>% 
  ungroup()

# Print excluded series for transparency
excluded <- qualified_series %>% filter(!pass)
if (nrow(excluded) > 0) {
  cat("\nExcluded Series (failed 79% completeness criterion):\n")
  print(excluded %>% select(indicator, subDir, shortName, fracComplete))
}

# Keep only the qualified series data
qualified_data <- allData %>% 
  inner_join(
    qualified_series %>% filter(pass) %>% select(indicator, subDir, shortName), 
    by = c("indicator", "subDir", "shortName")
  )

# ---------------------------------------------------------------------------
# Deduplicate shortName before Pivoting
# (If the same shortName passed screening from multiple subDirs, take the mean value)
# ---------------------------------------------------------------------------
final_long_data <- qualified_data %>% 
  group_by(shortName, date) %>% 
  summarize(finalVal = mean(finalVal, na.rm = TRUE), .groups = "drop") %>% 
  filter(!is.nan(finalVal))

# ---------------------------------------------------------------------------
# Format Wide Matrix (1996–2025)
# ---------------------------------------------------------------------------
datWide <- target_years %>% 
  left_join(final_long_data, by = "date") %>% 
  pivot_wider(
    id_cols = date, 
    names_from = shortName, 
    values_from = finalVal
  ) %>% 
  mutate(Year = year(date)) %>% 
  select(Year, everything(), -date) %>% 
  arrange(Year)

datWide$is_fully_complete <- complete.cases(datWide)

# Export wide data frame
write.csv(datWide, file.path(outputDir, "datWide_1996_2025_qualified.csv"), row.names = FALSE)

cat("\nPipeline Complete! Retained", ncol(datWide) - 1, "unique indicators across 1996-2025.\n")
head(datWide)