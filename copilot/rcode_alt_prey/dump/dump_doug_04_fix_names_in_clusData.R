


# ==============================================================================
# Script: doug_04_fix_names_in_clusData.R
# Purpose: Rename all headers in clusDataDFA_wide_1998_2021_altprey_extended.csv
#          to standardized LisaNames.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS
# ------------------------------------------------------------------------------
rootdir  <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"
meta_dir <- file.path(rootdir, "Rcode_for_paper/metadata")

target_csv_path <- file.path("copilot/outputs_altprey", "clusDataDFA_wide_1998_2021_altprey_extended.csv")
crosswalk_path  <- file.path(meta_dir, "master_name_crosswalk.csv")

if (!file.exists(target_csv_path)) stop("Missing target extended CSV at: ", target_csv_path)
if (!file.exists(crosswalk_path))  stop("Missing master crosswalk at: ", crosswalk_path)

# ------------------------------------------------------------------------------
# STEP 1: LOAD MATRIX AND CROSSWALK
# ------------------------------------------------------------------------------
ext_df <- read.csv(target_csv_path, stringsAsFactors = FALSE)
master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

current_cols <- names(ext_df)

# Build a lookup vector: DFAname -> LisaName
name_map <- setNames(master_crosswalk$LisaName, master_crosswalk$DFAname)

# ------------------------------------------------------------------------------
# STEP 2: MAP COLUMN HEADERS TO LISANAMES
# ------------------------------------------------------------------------------
new_cols <- sapply(current_cols, function(col) {
  if (col == "Year" || col == "year") return("Year")
  
  # 1. Check if column is already a valid LisaName (e.g. the 7 new 2yrLead cols or X16_SAR)
  if (col %in% master_crosswalk$LisaName || str_ends(col, "_2yrLead")) return(col)
  
  # 2. Strip R's auto-prepended 'X' from numeric DFA names (e.g., X07.Cond2NCC_DFA1 -> 07.Cond2NCC_DFA1)
  clean_col <- str_remove(col, "^X(?=\\d)")
  
  # 3. Match against crosswalk DFAname
  if (clean_col %in% names(name_map)) return(name_map[[clean_col]])
  
  # Return original if unmapped
  return(col)
})

# Verify mapping coverage
unmapped <- current_cols[new_cols == current_cols & current_cols != "Year"]
if (length(unmapped) > 0) {
  cat("Warning: Unmapped columns remaining:\n")
  print(unmapped)
} else {
  cat("SUCCESS: All columns successfully mapped to standardized LisaNames!\n")
}

# Assign standardized names
names(ext_df) <- unname(new_cols)

# ------------------------------------------------------------------------------
# STEP 3: OVERWRITE CSV IN METADATA DIRECTORY
# ------------------------------------------------------------------------------
write.csv(ext_df, file.path(target_csv_path), row.names = FALSE)

cat("\nSuccessfully updated:", target_csv_path, "\n")
cat("Total columns:", ncol(ext_df), "\n\n")
print(sort(names(ext_df)))
