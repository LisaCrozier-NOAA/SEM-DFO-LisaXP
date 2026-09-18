library(dplyr)

# Path to metadata file
meta_file <- "Rcode_for_paper/metadata/sem_data_all_1998_2021.csv"
sem_df    <- read.csv(meta_file, stringsAsFactors = FALSE)

# Standardize X15 adult-timing names
sem_df_updated <- sem_df %>%
  rename_with(
    ~ paste0(gsub("_predAK$", "", .x), "_2yrLead_predAK"),
    .cols = setdiff(x15_adult_timing_cols, grep("sleeperSharks", x15_adult_timing_cols, value = TRUE))
  ) %>%
  rename_with(
    ~ paste0(.x, "_2yrLead"),
    .cols = intersect(names(.), "X15_DFA_sleeperSharks")
  )

# Overwrite metadata file with clean standardized headers
write.csv(sem_df_updated, meta_file, row.names = FALSE)
cat("\nUpdated X15 column headers saved to metadata!\n")