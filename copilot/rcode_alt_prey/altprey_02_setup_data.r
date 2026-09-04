

#This script adds 2yrlead columns for ncc altprey, eulachon data and single prey columns from ak.yr that were not there
# # Export Products
# write.csv(sem_altprey_data, 
#           file.path(new_outputDir, "sem_altprey_data.csv"), row.names = FALSE)
# 
# write.csv(sem_altprey_data_complete_1998_2021, 
#           file.path(new_outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)
# 
# write.csv(clusData_lisa, 
#           file.path(new_outputDir, "clusData_wide_Lisanames.csv"), row.names = FALSE)
# 
# cat("Scaled complete dataset (1998-2021) successfully exported to outputs_9!\n")
# 

library(dplyr)
library(tidyr)
library(stringr)

# Define Output Directory (Outputs 9)
new_outputDir <- file.path("copilot/outputs_9")
dir.create(new_outputDir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# SECTION 1: Load Primary Inputs & Crosswalk
# ===========================================================================
clusDataDFA_wide      <- read.csv("copilot/outputs_8/clusDataDFA_wide_1996_2025.csv")
datWide_qualified     <- read.csv("copilot/outputs_8/datWide_1996_2025_qualified.csv")
master_name_crosswalk <- read.csv(file.path(new_outputDir, "master_name_crosswalk.csv"))

# Build named vector: values = dfa_cols (old names), names = Lisaname (new names)
name_map <- setNames(master_name_crosswalk$dfa_cols, master_name_crosswalk$Lisaname)

# Explicitly map SAR_DFA1 to X16_SAR
name_map["X16_SAR"] <- "SAR_DFA1"

# ===========================================================================
# SECTION 2: Extract Age 1+ Pollock directly from datWide_qualified
# ===========================================================================
pollock_age1_series <- datWide_qualified %>%
  select(Year, X13_pollock_age1plus = pollockBiomassAIage1plus_predAK_2026) %>%
  filter(!is.na(Year))

# ===========================================================================
# SECTION 3: Rename clusDataDFA_wide to Lisanames & Join Age 1+ Pollock
# ===========================================================================
clusData_lisa <- clusDataDFA_wide %>%
  rename(any_of(name_map)) %>%
  left_join(pollock_age1_series, by = "Year")

# ===========================================================================
# SECTION 4: Create 2-Year Lead Versions for Specific Prey ONLY
# (No leads on Mammals or Guild 15)
# ===========================================================================
prey_lead_cols <- c(
  "X04_marketsquid_GAM",
  "X05_DFA_abundSardine",        # Added so herring 2yrlead exists
  "X05_anchovy_GAM",
  "X09_DFA_ChinAbundSnakeFall",
  "X09_DFA_HakeAge5Plus",
  "X12_DFA_biomassEuphShelfSum",
  "X13_pollock_age1plus"
)

clusData_led <- clusData_lisa %>%
  mutate(across(
    all_of(intersect(prey_lead_cols, names(clusData_lisa))),
    ~ dplyr::lead(.x, 2),
    .names = "{.col}_2yrlead"
  ))

# ===========================================================================
# SECTION 5: Integrate Supplemental Datasets (AK Prey & Eulachon)
# ===========================================================================
eulachon_scaled <- read.csv("copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv") %>%
  select(Year, eulachon_during_chinook, eulachon_annual_scaled) %>%
  mutate(
    eulachon_during_chinook_2yrlead = dplyr::lead(eulachon_during_chinook, 2),
    eulachon_annual_scaled_2yrlead  = dplyr::lead(eulachon_annual_scaled, 2)
  )

ak_yr_selected <- read.csv("data_Lisa/ak_yr.csv", row.names = NULL) %>%
  rename(Year = year) %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  select(
    Year,
    X10_ssl_seak_pup_pred = ssl_seak_pup_pred,
    X12_egoa.krill        = secm_euph_dens,
    X13_stka_herr_matbiom  = stka_herr_matbiom,
    X13_mid_il_capelin    = mid_il_capelin
  )

# ===========================================================================
# SECTION 6: Assemble Master SEM Datasets, Scale & Export
# ===========================================================================
# 1. Full Master Dataset
sem_altprey_data <- clusData_led %>%
  left_join(eulachon_scaled, by = "Year") %>%
  left_join(ak_yr_selected, by = "Year") %>%
  arrange(Year)

# 2. Filter 1998-2021 window and retain complete columns (0 NAs)
df_1998_2021 <- sem_altprey_data %>%
  filter(Year >= 1998 & Year <= 2021)

na_counts <- colSums(is.na(df_1998_2021))
complete_cols <- names(na_counts[na_counts == 0])

# 3. Select complete columns and z-score scale all indicators
sem_altprey_data_complete_1998_2021 <- df_1998_2021 %>%
  select(all_of(complete_cols)) %>%
  mutate(across(-Year, ~ as.numeric(scale(.x))))

# Export Products
write.csv(sem_altprey_data, 
          file.path(new_outputDir, "sem_altprey_data.csv"), row.names = FALSE)

write.csv(sem_altprey_data_complete_1998_2021, 
          file.path(new_outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)

write.csv(clusData_lisa, 
          file.path(new_outputDir, "clusData_wide_Lisanames.csv"), row.names = FALSE)

cat("Scaled complete dataset (1998-2021) successfully exported to outputs_9!\n")
