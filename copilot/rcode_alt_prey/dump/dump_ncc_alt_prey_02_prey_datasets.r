

# This script processes eulachon annual and weekly datafiles from outputs_csl_cr
# then it joins that with the extended prey datafiles from doug_01 and doug_02 (clusDataDFA_wide_1996_2025.csv) and 

library(dplyr)
library(stringr)
library(imputeTS)

outputDir <- file.path( "copilot/outputs_8")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load and Process Eulachon Datasets
# ---------------------------------------------------------------------------

# Annual Run Index (Reconstructed Lbs)
eulachon_annual <- read.csv(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/eulachon_master_index_lbs_1993_2024.csv",
  row.names = NULL
) %>%
  select(year, eulachon_lbs_reconstructed) %>%
  arrange(year)

eulachon_annual$eulachon_lbs_interp <- na_interpolation(eulachon_annual$eulachon_lbs_reconstructed)

# Weekly Overlap Index (Active Chinook Window)
eulachon_weekly <- read.csv(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/master_estuary_2wklagweekly_all_species_1998_2024.csv",
  row.names = NULL
) %>%
  group_by(year) %>%
  summarise(eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE), .groups = "drop") %>%
  select(year, eulachon_during_chinook)

# Combine and apply Log + Z-score scaling
eulachon_clean <- full_join(eulachon_weekly, eulachon_annual %>% select(year, eulachon_lbs_interp), by = "year") %>%
  rename(Year = year) %>%
  mutate(
    # Primary index: Weekly overlap
    eulachon_weekly_scaled = as.vector(scale(log1p(eulachon_during_chinook))),
    # Secondary index: Annual reconstructed run size
    eulachon_annual_scaled = as.vector(scale(log1p(eulachon_lbs_interp)))
  )

write.csv(eulachon_clean, file.path(outputDir, "eulachon_scaled_annual_weekly_during_chinook.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Load DFA / Smoothed Prey Dataset and Merge
# ---------------------------------------------------------------------------

outputDir <- file.path("copilot/outputs_8")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))

# Load cleaned eulachon dataset
eulachon_clean <- read.csv(file.path(outputDir, "eulachon_scaled_annual_weekly_during_chinook.csv")) %>%
  select(Year, eulachon_during_chinook=eulachon_weekly_scaled, eulachon_annual=eulachon_annual_scaled)


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

