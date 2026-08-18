#Step 1. reconstruct eulachon annual abundance 1976-2024
rm(list=ls())


#Loads and cleans eulachon status review

# week_all<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),row.names=NULL)
# head(week_all)
# names(week_all)

# -----------------------------------------------------------------------------
# 1. Read, Clean, and Format Historical Annual Eulachon Data
# -----------------------------------------------------------------------------
library(tidyverse)
library(readxl)
path <- "data_Lisa/Eulachon_data_Gustavson_2010_status_review.xlsx"


# Helper function to convert commas, dashes, and "Unknown*" text into clean numeric / NA
clean_numeric <- function(x) {
  x_clean <- gsub(",", "", as.character(x))
  x_clean <- ifelse(grepl("Unknown|—|-|N/A", x_clean, ignore.case = TRUE), NA_character_, x_clean)
  suppressWarnings(as.numeric(x_clean))
}

eulachon_count_50yr <- read_xlsx(path = path, sheet = 2, skip = 1) %>%
  rename(
    year                             = Year,
    eulachon_cr_pounds               = `Total landings\r\n(pounds)`,
    eulachon_cr_nfish_10.8_per_pound = `Number of fish at\r\n10.8 per pound`,
    eulachon_cr_nfish_12.3_per_pound = `Number of fish at\r\n12.3 per pound`
  ) %>%
  mutate(
    year = as.integer(clean_numeric(year)),
    across(starts_with("eulachon_cr_"), clean_numeric)
  ) %>%
  # Filter to the last 50 years of data
  filter(!is.na(year), year >= (1976))

cat("--- Eulachon Annual Data (Last 50 Years) ---\n")
print(eulachon_count_50yr, n = 50)
