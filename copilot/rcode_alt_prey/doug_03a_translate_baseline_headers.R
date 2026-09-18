



# ==============================================================================
# Script: doug_03a_translate_baseline_headers.R
# Purpose: Cleanly translate raw DFA header names in the baseline matrix
#          to official LisaNames using master_name_crosswalk.csv.
# Output: copilot/outputs_altprey/clusDataDFA_wide_1998_2021_LisaNames.csv
# ==============================================================================

library(dplyr)
library(stringr)

outputDir <- "copilot/outputs_altprey"
metaDir   <- "Rcode_for_paper/metadata"

baseline_path  <- file.path(outputDir, "clusDataDFA_wide_1998_2021.csv")
crosswalk_path <- file.path(metaDir,   "master_name_crosswalk.csv")

if (!file.exists(baseline_path))  stop("Missing baseline file: ", baseline_path)
if (!file.exists(crosswalk_path)) stop("Missing crosswalk file: ", crosswalk_path)

baseline  <- read.csv(baseline_path, stringsAsFactors = FALSE)
crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

clean_names <- sapply(names(baseline), function(col) {
  if (tolower(col) == "year") return("Year")
  if (col %in% crosswalk$LisaName) return(col)
  
  col_no_x <- str_remove(col, "^[Xx](?=\\d)")
  
  match_idx <- match(col_no_x, crosswalk$DFAname)
  if (!is.na(match_idx)) {
    return(crosswalk$LisaName[match_idx])
  }
  
  return(col)
})

names(baseline) <- unname(clean_names)

out_path  <- file.path(outputDir, "clusDataDFA_wide_1998_2021_LisaNames.csv")
meta_path <- file.path(metaDir,   "clusDataDFA_wide_1998_2021_LisaNames.csv")

write.csv(baseline, out_path,  row.names = FALSE)
write.csv(baseline, meta_path, row.names = FALSE)

cat("\nStandardized Baseline Created!\nSaved to:\n - ", out_path, "\n - ", meta_path, "\n")



