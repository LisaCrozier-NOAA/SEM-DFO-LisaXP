#Step 2: Universal Renaming Utility Function
#To make renaming automatic across all future scripts, add this short utility function into your pipeline helper (or at the top of your scripts). 
#It reads master_name_crosswalk.csv and renames dataset columns from dfa_cols (or raw indicator names) to your standardized LisaName.


#Step 3: Integrating the Naming Standardization Across Pipeline Modules
#Whenever you load wide DFA matrices or process results in downstream scripts (such as Lavaan SEM runs or variable importance tables), 
#call apply_lisa_names() right after loading the dataset.
# This ensures that:
# Every downstream script, figure, Lavaan model output, and table uses identical, readable column names.
# You never have to manually write inline case_when() statements or crosswalk logic again.
# If you decide to tweak a variable's label in the future, you only update master_name_crosswalk.csv once, and all scripts reflect the update automatically.


#00b_crosswalk_utility_fxn.r

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


