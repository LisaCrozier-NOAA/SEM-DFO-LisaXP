library(dplyr)
library(tidyr)
library(stringr)

# Define Output Directory (Outputs 9)
new_outputDir <- file.path("copilot/outputs_9")
dir.create(new_outputDir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# SECTION 1: Load Primary Inputs & Crosswalk (NO sem_data dependencies)
# ===========================================================================
clusDataDFA_wide      <- read.csv("copilot/outputs_8/clusDataDFA_wide_1996_2025.csv")
datWide_qualified     <- read.csv("copilot/outputs_8/datWide_1996_2025_qualified.csv")
master_name_crosswalk <- read.csv(file.path(new_outputDir, "master_name_crosswalk.csv"))

# Build inverted name map (c("dfa_cols" = "Lisaname") -> c("Lisaname" = "dfa_cols"))
# Ensures rename(any_of(name_map)) properly maps new_name = old_name
name_map <- setNames(master_name_crosswalk$dfa_cols, master_name_crosswalk$Lisaname)

# Explicitly map SAR_DFA1 to X16_SAR in case dfa_cols retained the extra X prefix
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
# SECTION 4: Create 2-Year Lead Versions for Mammals, Prey, & Guild 15
# ===========================================================================
# # Marine mammals shifted in-place
# mammal_lead_cols <- c(
#   "X10_DFA_ssl.est.wholerange_2yrLead",
#   "X10_Californian_s_l_2yrLead_WS",
#   "X10_Harbour_s_2yrLead_WS",
#   "X10_Northern_f_s_2yrLead_WS"
# )

# Base prey lead columns
prey_lead_cols <- c(
  "X09_DFA_ChinAbundSnakeFall",
  "X09_DFA_HakeAge5Plus",
  "X12_DFA_biomassEuphShelfSum",
  "X04_marketsquid_GAM",
  "X05_anchovy_GAM",
  "X13_pollock_age1plus"
)

# Dynamically extract all Guild 15 Alaska predator columns from Lisanames
#guild15_cols <- names(clusData_lisa)[grepl("^X15_", names(clusData_lisa))]
#cols_to_generate_lead <- unique(c(prey_lead_cols, guild15_cols))

clusData_led <- clusData_lisa %>%
  # 1. Shift marine mammals in-place
#  mutate(across(all_of(intersect(mammal_lead_cols, names(clusData_lisa))), ~ dplyr::lead(.x, 2))) %>%
  # 2. Generate explicit _2yrlead duplicate columns for prey and Guild 15 series
  mutate(across(
#    all_of(intersect(cols_to_generate_lead, names(clusData_lisa))),
    all_of(intersect(prey_lead_cols, names(clusData_lisa))),
    ~ dplyr::lead(.x, 2),
    .names = "{.col}_2yrlead"
  ))

# ===========================================================================
# SECTION 5: Integrate Supplemental Datasets (AK Prey & Eulachon)
# ===========================================================================
# Load Eulachon Index
eulachon_scaled <- read.csv("copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv") %>%
  select(Year, eulachon_during_chinook, eulachon_annual_scaled) %>%
  mutate(
    eulachon_during_chinook_2yrlead = dplyr::lead(eulachon_during_chinook, 2),
    eulachon_annual_scaled_2yrlead  = dplyr::lead(eulachon_annual_scaled, 2)
  )

# Load Alaska Series directly from ak_yr.csv
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
# SECTION 6: Assemble Master SEM Dataset & Export Final Products
# ===========================================================================
sem_altprey_data <- clusData_led %>%
  left_join(eulachon_scaled, by = "Year") %>%
  left_join(ak_yr_selected, by = "Year") %>%
  arrange(Year)

# Export Master Dataset & Individual Backup Files
write.csv(sem_altprey_data, file.path(new_outputDir, "sem_altprey_data.csv"), row.names = FALSE)
write.csv(clusData_lisa, file.path(new_outputDir, "clusData_wide_Lisanames.csv"), row.names = FALSE)
write.csv(datWide_qualified, file.path(new_outputDir, "all_individual_indicators_1996_2025.csv"), row.names = FALSE)

cat("Pipeline complete! All products generated with updated Lisanames and exported to copilot/outputs_9/sem_altprey_data.csv\n")

names(sem_altprey_data)
