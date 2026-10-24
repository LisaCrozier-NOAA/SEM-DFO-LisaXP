

#Script: 00_create_master_crosswalk.R
#Step 1: Master Name Crosswalk Generator Script
#Below is the updated, clean version of your crosswalk creation script (00_create_master_crosswalk.R). 
#It standardizes column names, cleans up year tags (_2025/_2026), maps the 17 DFAs and 48 stragglers/smoothed series, and 
#exports master_name_crosswalk.csv.
#to be used with utility function source(file.path(proj_dir, "00b_crosswalk_utility_fxn.r"))


# ==============================================================================
# Script: 00_create_master_crosswalk.R
# Purpose: Generate the Master Variable Name Crosswalk (LisaName vs Doug dfa_cols)
# ==============================================================================

library(tidyverse)
library(stringr)

# --- Configuration & Paths ---
proj_dir     <- paste0(getwd(),"/Rcode_for_paper")
meta_dir     <- file.path(proj_dir, "metadata")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"
#doug_dir     <- "2026_10_21_SEM/shiftLisa_1_Lisa_21oct26"


rankedIndicators0 <- read.csv(file.path(doug_dir,"DFA", "rankedIndicators.csv"), row.names = NULL) 

rankedIndicators<-rankedIndicators0 %>%
  mutate(dfa_cols = paste0("X", DFAname))%>% 
  mutate(indicator1 = trimws(sub("-.*$", "", rankedIndicators)))%>%
  mutate(indicator1 = trimws(sub("_+$", "", indicator1)))%>%
  mutate(indicator1 = gsub("_2025|_2026", "", indicator1))%>%
  mutate(node=substr(dfa_cols, 1, 3)) %>%
  mutate(LisaName = case_when(
    type=="DFA" ~ paste(node,type,indicator1,sep="_"),
    TRUE ~ paste(node,indicator1,sep="_")
  )) %>%
  relocate(rankedIndicators,.after = last_col()) %>%
  as_tibble()


rankedIndicators %>% head()
print(rankedIndicators,n=Inf)

write.csv(rankedIndicators, file.path(meta_dir, "master_name_crosswalk.csv"), row.names = FALSE)
write.csv(rankedIndicators, file.path(meta_dir, "crosswalk_names.csv"), row.names = FALSE)

crosswalk_names<-read.csv(file.path(meta_dir, "crosswalk_names.csv"), row.names = NULL)

#00b_crosswalk_utility_fxn.r-------

# Utility function to automatically apply your standardized LisaName convention
apply_lisa_names <- function(df, crosswalk_path = "metadata/master_name_crosswalk.csv") {
  if (!file.exists(crosswalk_path)) {
    warning("Crosswalk file not found at: ", crosswalk_path, ". Returning original data.")
    return(df)
  }
  
  crosswalk <- read.csv(crosswalk_path)
  
  # Build lookup dictionary (maps dfa_cols or raw names -> LisaName)
  lookup_map <- setNames(crosswalk$LisaName, crosswalk$dfa_cols)
  
  # Also match raw DFAname strings without 'X' prefix if present
  lookup_map_raw <- setNames(crosswalk$LisaName, crosswalk$DFAname)
  lookup_map <- c(lookup_map, lookup_map_raw)
  
  # Rename matching columns
  current_cols <- names(df)
  new_cols <- current_cols
  
  for (i in seq_along(current_cols)) {
    col <- current_cols[i]
    if (col %in% names(lookup_map)) {
      new_cols[i] <- lookup_map[[col]]
    }
  }
  
  names(df) <- new_cols
  return(df)
}


