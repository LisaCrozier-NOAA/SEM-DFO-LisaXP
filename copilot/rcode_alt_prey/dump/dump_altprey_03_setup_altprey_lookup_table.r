# ==============================================================================
# Script: altprey_03_setup_altprey_lookup_table.r
# Purpose: Build unified SEM dataset and construct AltPrey lookup table 
#          mapping adult predators (Guild 10 & 15) to _smoltYear columns.
# Output: copilot/outputs_altprey/sem_altprey_data_complete_1998_2021.csv
#         copilot/outputs_altprey/all_pred_dfa_altprey.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(janitor)

# ---------------------------------------------------------------------------
# Setup Directories
# ---------------------------------------------------------------------------
outputDir <- file.path("copilot/outputs_altprey")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
Rcode_for_paper_Dir <- file.path("Rcode_for_paper")

# ---------------------------------------------------------------------------
# 1. Load Extended Matrix & Additional Datasets
# ---------------------------------------------------------------------------
clusData_led <- read.csv(file.path(outputDir, "clusDataDFA_wide_1998_2021_altprey_extended_smoltyear_shifted_LisaNames.csv"))

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
    X11_ssl_seak_pup_pred = ssl_seak_pup_pred,
    X12_egoa.krill        = secm_euph_dens,
    X13_mid_il_capelin    = mid_il_capelin,
    X21_sst_egoa_junjulaug = sst_egoa_coastwatch_junjulaug
  )

# Assemble and export complete SEM dataset
sem_altprey_data <- clusData_led %>%
  left_join(eulachon_scaled, by = "Year") %>%
  left_join(ak_yr_selected, by = "Year") %>%
  arrange(Year) 

df_1998_2021 <- sem_altprey_data %>% filter(Year >= 1998 & Year <= 2021)
na_counts <- colSums(is.na(df_1998_2021))
complete_cols <- names(na_counts[na_counts == 0])

sem_complete_data <- df_1998_2021 %>%
  select(all_of(complete_cols)) %>%
  mutate(across(-Year, ~ as.numeric(scale(.x))))

# Export SEM Data Files
write.csv(sem_complete_data, file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)
write.csv(sem_complete_data, file.path(Rcode_for_paper_Dir, "metadata", "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)
write.csv(df_1998_2021, file.path(outputDir, "sem_data_1998_2021.csv"), row.names = FALSE)
write.csv(df_1998_2021, file.path(Rcode_for_paper_Dir, "metadata", "sem_data_1998_2021.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Build Unified Predator Assignment Matrix (NCC + AK)
# ---------------------------------------------------------------------------
excluded_preds <- c(
  "canaryrockfish", "chilipepper", "spinydogfishbsai", 
  "salmonsharkbsai", "pacificspinydogfish", "iphc20halibut", 
  "loon", "sharkcatchgoa", "harbor_seal_cr_2025"
)

combined_guildfiles <- read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/guildsWithExclude.csv", row.names = NULL) %>% 
  clean_names() %>% 
  rename(sem_name = se_mnode, short_name_lower = short_name) %>% 
  mutate(across(where(is.character), tolower)) %>% 
  filter(!is.na(latest_guild)) %>% 
  filter(!grepl(paste(excluded_preds, collapse = "|"), indicator))

# --- NCC PREDATOR ASSIGNMENTS ---
nccpred <- combined_guildfiles %>% 
  filter(grepl("08|09|10|11", latest_guild)) %>% 
  select(indicator, short_name_lower, sem_name, latest_guild) %>% 
  mutate(
    region = "NCC",
    altprey1 = case_when(
      grepl("killer.whale", indicator) ~ "chin",
      grepl("allsealionsbonn", indicator) ~ "eulachon",
      grepl("08\\.predbirdncc|peli|tern|murre|corm|alcid|shrw|shearwater|grebe", indicator) ~ "market_squid",
      grepl("hake|rockfish|mackerel|chinook abundance", indicator) ~ "krill",
      TRUE ~ "hake"
    ),
    altprey2 = case_when(
      grepl("californian_s_l|08\\.predbirdncc|peli|tern|murre|corm|alcid|shrw|shearwater|grebe|rockfish|mackerel|chinook abundance|hake", indicator) ~ "anchovy",
      TRUE ~ "herring"
    ),
    altprey3 = case_when(
      grepl("ssl.est.wholerange|allsealionsemb", indicator) ~ "eulachon",
      grepl("08\\.predbirdncc|peli|tern|corm|murre", indicator) ~ "herring",
      TRUE ~ NA_character_
    )
  )

# --- ALASKA PREDATOR ASSIGNMENTS ---
akpred <- combined_guildfiles %>% 
  filter(grepl("15", latest_guild)) %>% 
  select(indicator, short_name_lower, sem_name, latest_guild) %>% 
  mutate(
    region = "AK",
    altprey1 = case_when(
      grepl("spinydogfish", indicator) ~ "krill",
      TRUE ~ "pollock"
    ),
    altprey2 = case_when(
      grepl("salmonshark", indicator) ~ "pink_salmon",
      TRUE ~ "herring"
    ),
    altprey3 = case_when(
      grepl("sablefish", indicator) ~ "krill",
      grepl("pacificcod|arrowtooth|spinydogfish", indicator) ~ "capelin",
      TRUE ~ NA_character_
    )
  )

# Append Alaska SSL SEAK Pups as Guild 11
ssl_ak_row <- data.frame(
  indicator = "ssl_seak_pup_pred",
  short_name_lower = "ssl_seak_pup_pred",
  sem_name = "predak",
  latest_guild = "11",
  region = "AK",
  altprey1 = "herring",
  altprey2 = "capelin",
  altprey3 = "krill",
  stringsAsFactors = FALSE
)

all_preds_assigned <- bind_rows(nccpred, akpred, ssl_ak_row)

# ---------------------------------------------------------------------------
# 3. Map Exact Column Names (Both Original Adult & _smoltYear Variants)
# ---------------------------------------------------------------------------

# 1. Standard Base Mapping (Preserves original predator data column names)
all_preds_base <- all_preds_assigned %>%
  mutate(
    is_adult_predator = grepl("^10|^15", latest_guild),
    
    # --- Map Original Pred Data Column As-Is ---
    pred_data_col = case_when(
      grepl("ssl.est.wholerange", indicator)   ~ "X10_DFA_ssl.est.wholerange_2yrLead",
      grepl("californian_s_l", indicator)      ~ "X10_Californian_s_l_2yrLead_WS",
      grepl("harbour_s", indicator)            ~ "X10_Harbour_s_2yrLead_WS",
      grepl("northern_f_s", indicator)         ~ "X10_Northern_f_s_2yrLead_WS",
      grepl("ssl_seak_pup_pred", indicator)    ~ "X11_ssl_seak_pup_pred",
      grepl("harbour_p_ws", indicator)         ~ "X11_DFA_Harbour_p_WS",
      grepl("chinabundsnakefall", indicator)   ~ "X09_DFA_ChinAbundSnakeFall",
      grepl("hakeage5plus", indicator)         ~ "X09_DFA_HakeAge5Plus",
      grepl("08\\.predbirdncc|corm", indicator) ~ "X08_DFA_DC_corm_3_WS",
      grepl("gull", indicator)                 ~ "X08_Large_gulls_7_WS",
      grepl("arrowtooth", indicator)           ~ "X15_ArrowtoothFlounderBiomass_predAK",
      grepl("sablefishbiomass", indicator)     ~ "X15_sablefishBiomass_predAK",
      grepl("sablefishrecruitment", indicator) ~ "X15_sablefishRecruitment_predAK_2yrLead",
      grepl("pacificcod", indicator)           ~ "X15_PacificCodBiomass_predAK",
      grepl("spinydogfishgoa", indicator)      ~ "X15_spinyDogfishGoA_predAK_2yrLead",
      grepl("halibutbiomassage8plus", indicator)~ "X15_halibutBiomassAge8plus_2yrLead_predAK",
      grepl("sleepershark", indicator)         ~ "X15_DFA_sleeperSharks",
      grepl("salmonsharkgoa", indicator)       ~ "X15_salmonSharkGoA_predAK",
      TRUE ~ short_name_lower
    ),
    
    # --- Map Prey Columns for Original Predators (Using 2yrLead where applicable) ---
    altprey1_data_col = case_when(
      altprey1 == "chin"         ~ ifelse(is_adult_predator, "X09_DFA_ChinAbundSnakeFall_2yrLead", "X09_DFA_ChinAbundSnakeFall"),
      altprey1 == "hake"         ~ ifelse(is_adult_predator, "X09_DFA_HakeAge5Plus_2yrLead", "X09_DFA_HakeAge5Plus"),
      altprey1 == "eulachon"     ~ ifelse(is_adult_predator, "eulachon_during_chinook_2yrlead", "eulachon_during_chinook"),
      altprey1 == "krill"        ~ ifelse(region == "AK", "X12_egoa.krill", ifelse(is_adult_predator, "X12_DFA_biomassEuphShelfSum_2yrLead", "X12_DFA_biomassEuphShelfSum")),
      altprey1 == "market_squid" ~ ifelse(is_adult_predator, "X04_marketsquid_GAM_2yrLead", "X04_marketsquid_GAM"),
      altprey1 == "pollock"      ~ ifelse(is_adult_predator, "X13_pollock_age1plus_2yrLead", "X13_pollock_age1plus"),
      altprey1 == "herring"      ~ ifelse(region == "AK", 
                                          ifelse(is_adult_predator, "X13_sitkaHerring_EGoA_2yrLead", "X13_sitkaHerring_EGoA"), 
                                          ifelse(is_adult_predator, "X05_DFA_abundSardine_2yrLead", "X05_DFA_abundSardine")),
      TRUE ~ NA_character_
    ),
    
    altprey2_data_col = case_when(
      altprey2 == "anchovy"     ~ ifelse(is_adult_predator, "X05_anchovy_GAM_2yrLead", "X05_anchovy_GAM"),
      altprey2 == "herring"     ~ ifelse(region == "AK", 
                                         ifelse(is_adult_predator, "X13_sitkaHerring_EGoA_2yrLead", "X13_sitkaHerring_EGoA"), 
                                         ifelse(is_adult_predator, "X05_DFA_abundSardine_2yrLead", "X05_DFA_abundSardine")),
      altprey2 == "pink_salmon" ~ "X14_pinkSalmonNorthAmerica",
      altprey2 == "capelin"     ~ "X13_mid_il_capelin",
      TRUE ~ NA_character_
    ),
    
    altprey3_data_col = case_when(
      altprey3 == "eulachon"    ~ ifelse(is_adult_predator, "eulachon_during_chinook_2yrlead", "eulachon_during_chinook"),
      altprey3 == "anchovy"     ~ ifelse(is_adult_predator, "X05_anchovy_GAM_2yrLead", "X05_anchovy_GAM"),
      altprey3 == "herring"     ~ ifelse(region == "AK", 
                                         ifelse(is_adult_predator, "X13_sitkaHerring_EGoA_2yrLead", "X13_sitkaHerring_EGoA"), 
                                         ifelse(is_adult_predator, "X05_DFA_abundSardine_2yrLead", "X05_DFA_abundSardine")),
      altprey3 == "capelin"     ~ "X13_mid_il_capelin",
      altprey3 == "krill"       ~ ifelse(region == "AK", "X12_egoa.krill", ifelse(is_adult_predator, "X12_DFA_biomassEuphShelfSum_2yrLead", "X12_DFA_biomassEuphShelfSum")),
      TRUE ~ NA_character_
    )
  )

# 2. Duplicate adult predator rows specifically for the _smoltYear variants
all_preds_smoltYear <- all_preds_base %>%
  filter(is_adult_predator) %>%
  mutate(
    # Append _smoltYear to predator names
    pred_data_col = case_when(
      grepl("ssl.est.wholerange", indicator)   ~ "X10_DFA_ssl.est.wholerange_smoltYear",
      grepl("californian_s_l", indicator)      ~ "X10_Californian_s_l_WS_smoltYear",
      grepl("harbour_s", indicator)            ~ "X10_Harbour_s_WS_smoltYear",
      grepl("northern_f_s", indicator)         ~ "X10_Northern_f_s_WS_smoltYear",
      grepl("arrowtooth", indicator)           ~ "X15_ArrowtoothFlounderBiomass_predAK_smoltYear",
      grepl("sablefishbiomass", indicator)     ~ "X15_sablefishBiomass_predAK_smoltYear",
      grepl("sablefishrecruitment", indicator) ~ "X15_sablefishRecruitment_predAK_smoltYear",
      grepl("pacificcod", indicator)           ~ "X15_PacificCodBiomass_predAK_smoltYear",
      grepl("spinydogfishgoa", indicator)      ~ "X15_spinyDogfishGoA_predAK_smoltYear",
      grepl("halibutbiomassage8plus", indicator)~ "X15_halibutBiomassAge8plus_predAK_smoltYear",
      grepl("sleepershark", indicator)         ~ "X15_DFA_sleeperSharks_smoltYear",
      grepl("salmonsharkgoa", indicator)       ~ "X15_salmonSharkGoA_predAK_smoltYear",
      TRUE ~ paste0(pred_data_col, "_smoltYear")
    ),
    
    # Map prey to in-year smolt-entry versions
    altprey1_data_col = case_when(
      altprey1 == "chin"         ~ "X09_DFA_ChinAbundSnakeFall",
      altprey1 == "hake"         ~ "X09_DFA_HakeAge5Plus",
      altprey1 == "eulachon"     ~ "eulachon_during_chinook",
      altprey1 == "krill"        ~ ifelse(region == "AK", "X12_egoa.krill", "X12_DFA_biomassEuphShelfSum"),
      altprey1 == "market_squid" ~ "X04_marketsquid_GAM",
      altprey1 == "pollock"      ~ "X13_pollock_age1plus",
      altprey1 == "herring"      ~ ifelse(region == "AK", "X13_sitkaHerring_EGoA", "X05_DFA_abundSardine"),
      TRUE ~ NA_character_
    ),
    
    altprey2_data_col = case_when(
      altprey2 == "anchovy"     ~ "X05_anchovy_GAM",
      altprey2 == "herring"     ~ ifelse(region == "AK", "X13_sitkaHerring_EGoA", "X05_DFA_abundSardine"),
      altprey2 == "pink_salmon" ~ "X14_pinkSalmonNorthAmerica",
      altprey2 == "capelin"     ~ "X13_mid_il_capelin",
      TRUE ~ NA_character_
    ),
    
    altprey3_data_col = case_when(
      altprey3 == "eulachon"    ~ "eulachon_during_chinook",
      altprey3 == "anchovy"     ~ "X05_anchovy_GAM",
      altprey3 == "herring"     ~ ifelse(region == "AK", "X13_sitkaHerring_EGoA", "X05_DFA_abundSardine"),
      altprey3 == "capelin"     ~ "X13_mid_il_capelin",
      altprey3 == "krill"       ~ ifelse(region == "AK", "X12_egoa.krill", "X12_DFA_biomassEuphShelfSum"),
      TRUE ~ NA_character_
    )
  )

# 3. Combine both original and _smoltYear rows
all_pred_dfa_altprey <- bind_rows(all_preds_base, all_preds_smoltYear) %>%
  select(-is_adult_predator) %>%
  arrange(latest_guild, indicator)



names(all_pred_dfa_altprey)
all_pred_dfa_altprey %>% select(pred_data_col,altprey1_data_col,altprey2_data_col,altprey3_data_col)

# Export Final Lookup Table
out_lookup_path <- file.path(outputDir, "all_pred_dfa_altprey.csv")
write.csv(all_pred_dfa_altprey, out_lookup_path, row.names = FALSE)

cat("\nLookup Table Script Complete!")
cat("\nSaved updated lookup table to:", out_lookup_path, "\n")
