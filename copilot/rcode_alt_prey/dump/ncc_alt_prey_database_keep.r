library(dplyr)
library(stringr)

# ---------------------------------------------------------------------------
# 1. Build Master Alternate Prey Time Series preserving raw data column names
# ---------------------------------------------------------------------------

outputDir <- file.path(rootdir, "copilot/outputs_8")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))

# Load cleaned eulachon dataset
eulachon_clean <- read.csv(file.path(outputDir, "altprey_wide_1996_2025.csv")) %>%
  select(Year, eulachon_during_chinook = eulachon, eulachon_annual)

# Construct wide prey dataset with exact raw source column names and 2-year leads
altprey_master <- clusDataDFA_wide %>%
  select(
    Year,
    X09.PredFishNCC_DFA1,        # chin
    X09.PredFishNCC_b_DFA1,      # hake
    X12.ZooPreyAK_DFA1,          # krill
    X04.CompFishNCC_smoothed,    # market_squid
    X05.ForageFishNCC_0_smoothed,# anchovy
    X13.FishPreyAK_b_DFA1        # herring
  ) %>%
  left_join(eulachon_clean, by = "Year") %>%
  mutate(
    # Generate matching _2yrlead columns using exact source names
    X09.PredFishNCC_DFA1_2yrlead         = lead(X09.PredFishNCC_DFA1, 2),
    X09.PredFishNCC_b_DFA1_2yrlead       = lead(X09.PredFishNCC_b_DFA1, 2),
    X12.ZooPreyAK_DFA1_2yrlead           = lead(X12.ZooPreyAK_DFA1, 2),
    X04.CompFishNCC_smoothed_2yrlead     = lead(X04.CompFishNCC_smoothed, 2),
    X05.ForageFishNCC_0_smoothed_2yrlead = lead(X05.ForageFishNCC_0_smoothed, 2),
    X13.FishPreyAK_b_DFA1_2yrlead        = lead(X13.FishPreyAK_b_DFA1, 2),
    eulachon_during_chinook_2yrlead      = lead(eulachon_during_chinook, 2),
    eulachon_annual_2yrlead              = lead(eulachon_annual, 2)
  )

# Save master prey dataset with preserved source column names
write.csv(altprey_master, file.path(outputDir, "altprey_master_source_names_1996_2025.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Map exact data source column names into nccpred
# ---------------------------------------------------------------------------
nccpred <- read.csv(file.path(outputDir, "nccpred_altprey_assigned.csv"), row.names = NULL)

nccpred_mapped <- nccpred %>%
  mutate(
    # Check if indicator operates on 2-year lead
    is_2yrlead = str_detect(indicator, regex("_2yrlead|_2yrLead", ignore_case = TRUE)),
    
    # Map exact source column name for altprey1
    altprey1_data_col = case_when(
      altprey1 == "chin"         & !is_2yrlead ~ "X09.PredFishNCC_DFA1",
      altprey1 == "chin"         &  is_2yrlead ~ "X09.PredFishNCC_DFA1_2yrlead",
      
      altprey1 == "hake"         & !is_2yrlead ~ "X09.PredFishNCC_b_DFA1",
      altprey1 == "hake"         &  is_2yrlead ~ "X09.PredFishNCC_b_DFA1_2yrlead",
      
      altprey1 == "krill"        & !is_2yrlead ~ "X12.ZooPreyAK_DFA1",
      altprey1 == "krill"        &  is_2yrlead ~ "X12.ZooPreyAK_DFA1_2yrlead",
      
      altprey1 == "market_squid" & !is_2yrlead ~ "X04.CompFishNCC_smoothed",
      altprey1 == "market_squid" &  is_2yrlead ~ "X04.CompFishNCC_smoothed_2yrlead",
      
      altprey1 == "anchovy"      & !is_2yrlead ~ "X05.ForageFishNCC_0_smoothed",
      altprey1 == "anchovy"      &  is_2yrlead ~ "X05.ForageFishNCC_0_smoothed_2yrlead",
      
      altprey1 == "herring"      & !is_2yrlead ~ "X13.FishPreyAK_b_DFA1",
      altprey1 == "herring"      &  is_2yrlead ~ "X13.FishPreyAK_b_DFA1_2yrlead",
      
      altprey1 == "eulachon"     & !is_2yrlead ~ "eulachon_during_chinook",
      altprey1 == "eulachon"     &  is_2yrlead ~ "eulachon_during_chinook_2yrlead",
      
      TRUE ~ NA_character_
    ),
    
    # Map exact source column name for altprey2
    altprey2_data_col = case_when(
      altprey2 == "salmonids"                  ~ NA_character_, # Focal salmonids set to NA
      
      altprey2 == "anchovy"      & !is_2yrlead ~ "X05.ForageFishNCC_0_smoothed",
      altprey2 == "anchovy"      &  is_2yrlead ~ "X05.ForageFishNCC_0_smoothed_2yrlead",
      
      altprey2 == "herring"      & !is_2yrlead ~ "X13.FishPreyAK_b_DFA1",
      altprey2 == "herring"      &  is_2yrlead ~ "X13.FishPreyAK_b_DFA1_2yrlead",
      
      altprey2 == "eulachon"     & !is_2yrlead ~ "eulachon_during_chinook",
      altprey2 == "eulachon"     &  is_2yrlead ~ "eulachon_during_chinook_2yrlead",
      
      TRUE ~ NA_character_
    )
  )

# Export updated nccpred dataset
write.csv(nccpred_mapped, file.path(outputDir, "nccpred_mapped_source_names.csv"), row.names = FALSE)

# Display mapping verification
print(as_tibble(nccpred_mapped %>% select(indicator, altprey1, altprey1_data_col, altprey2, altprey2_data_col)), n = Inf)