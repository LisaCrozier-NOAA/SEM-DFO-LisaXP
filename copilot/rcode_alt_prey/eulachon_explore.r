eulachon_weekly <- read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/master_estuary_2wklagweekly_all_species_1998_2024.csv",row.names = NULL) %>%
  group_by(year) %>%
  summarise(eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE)) %>%
  select(year, eulachon_during_chinook)

eulachon_weekly


eulachon_annual <- read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/eulachon_master_index_lbs_1993_2024.csv",row.names = NULL) %>%
 select(year, eulachon_lbs_reconstructed)

eulachon_annual


library(dplyr)
library(imputeTS)

# Load Annual Run Index (Reconstructed Lbs)
eulachon_annual <- read.csv(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/eulachon_master_index_lbs_1993_2024.csv",
  row.names = NULL
) %>%
  select(year, eulachon_lbs_reconstructed) %>%
  arrange(year)

# Interpolate the single missing year in the annual series
eulachon_annual$eulachon_lbs_interp <- na_interpolation(eulachon_annual$eulachon_lbs_reconstructed)

# Load Weekly Overlap Index (Eulachon Present During Active Chinook Migration)
eulachon_weekly <- read.csv(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/master_estuary_2wklagweekly_all_species_1998_2024.csv",
  row.names = NULL
) %>%
  group_by(year) %>%
  summarise(eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE), .groups = "drop") %>%
  select(year, eulachon_during_chinook)

# Merge both indices
eulachon_comparison <- inner_join(eulachon_annual, eulachon_weekly, by = "year")

# Calculate Correlations across overlapping years (1998-2024)
pearson_corr  <- cor(eulachon_comparison$eulachon_lbs_interp, eulachon_comparison$eulachon_during_chinook, method = "pearson", use = "pairwise.complete.obs")
spearman_corr <- cor(eulachon_comparison$eulachon_lbs_interp, eulachon_comparison$eulachon_during_chinook, method = "spearman", use = "pairwise.complete.obs")

cat("Pearson Correlation (Linear):", round(pearson_corr, 3), "\n")
cat("Spearman Rank Correlation:", round(spearman_corr, 3), "\n\n")

# Display comparison table
print(as_tibble(eulachon_comparison), n = Inf)

#final prey database------------
library(dplyr)
library(tidyr)
library(lubridate)
library(imputeTS)

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
clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))

# Build final target prey dataset
altprey_wide <- clusDataDFA_wide %>%
  select(
    Year,
    chin         = X09.PredFishNCC_DFA1,
    hake         = X09.PredFishNCC_b_DFA1,
    krill        = X12.ZooPreyAK_DFA1,
    market_squid = X04.CompFishNCC_smoothed,
    anchovy      = X05.ForageFishNCC_0_smoothed,
    herring      = X13.FishPreyAK_b_DFA1
  ) %>%
  left_join(
    eulachon_clean %>% select(Year, eulachon = eulachon_weekly_scaled, eulachon_annual = eulachon_annual_scaled),
    by = "Year"
  ) %>%
  arrange(Year)

# Export complete prey dataset
write.csv(altprey_wide, file.path(outputDir, "altprey_wide_1996_2025.csv"), row.names = FALSE)

head(altprey_wide)
