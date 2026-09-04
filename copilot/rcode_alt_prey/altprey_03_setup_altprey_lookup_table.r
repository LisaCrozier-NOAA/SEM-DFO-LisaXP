#This script produces 
# Export Final Lookup Table
#write.csv(all_pred_dfa_altprey_clean, file.path( "copilot/outputs_9/all_pred_dfa_altprey.csv"), row.names = FALSE)
#names shoild match 
#write.csv(sem_altprey_data_complete_1998_2021, 
#file.path(new_outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)

#sem_complete_data<- read.csv("copilot/outputs_9/sem_altprey_data_complete_1998_2021.csv", row.names = FALSE)


library(dplyr)
library(tidyr)
library(stringr)
library(janitor)

# ---------------------------------------------------------------------------
# Setup Directories
# ---------------------------------------------------------------------------
outputDir <- file.path("copilot/outputs_9")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Update Master Assembly: Rename ssl_seak_pup_pred to X11_
# ---------------------------------------------------------------------------
clusData_led <- read.csv(file.path(outputDir, "clusData_wide_Lisanames.csv"))

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
    X11_ssl_seak_pup_pred = ssl_seak_pup_pred,  # Renamed X10 -> X11
    X12_egoa.krill        = secm_euph_dens,
    X13_stka_herr_matbiom  = stka_herr_matbiom,
    X13_mid_il_capelin    = mid_il_capelin
  )

# Assemble and export updated complete dataset
sem_altprey_data <- clusData_led %>%
  left_join(eulachon_scaled, by = "Year") %>%
  left_join(ak_yr_selected, by = "Year") %>%
  arrange(Year)

df_1998_2021 <- sem_altprey_data %>% filter(Year >= 1998 & Year <= 2021)
na_counts <- colSums(is.na(df_1998_2021))
complete_cols <- names(na_counts[na_counts == 0])

which(na_counts>0)

sem_complete_data <- df_1998_2021 %>%
  select(all_of(complete_cols)) %>%
  mutate(across(-Year, ~ as.numeric(scale(.x))))

write.csv(sem_complete_data, file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Build Unified Predator Assignment Matrix (NCC + AK)
# ---------------------------------------------------------------------------
excluded_preds <- c(
  "canaryrockfish", "chilipepper", "spinydogfishbsai", 
  "salmonsharkbsai", "pacificspinydogfish", "iphc20halibut", 
  "commonmurre_jsoes", "loon", "sharkcatchgoa"
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
# 3. Map Exact Column Names from sem_altprey_data_complete_1998_2021.csv
# ---------------------------------------------------------------------------
all_pred_dfa_altprey <- all_preds_assigned %>%
  mutate(
    # --- Map Pred Name in Data ---
    pred_data_col = case_when(
      # Guild 10 Marine Mammals (NCC Adult Salmon Predators)
      grepl("ssl.est.wholerange", indicator)  ~ "X10_DFA_ssl.est.wholerange_2yrLead",
      grepl("californian_s_l", indicator)     ~ "X10_Californian_s_l_2yrLead_WS",
      grepl("harbour_s", indicator)           ~ "X10_Harbour_s_2yrLead_WS",
      grepl("northern_f_s", indicator)        ~ "X10_Northern_f_s_2yrLead_WS",
      
      # Guild 11 Smolt/Juvenile Mammal Predators
      grepl("ssl_seak_pup_pred", indicator)   ~ "X11_ssl_seak_pup_pred",
      grepl("harbour_p_ws", indicator)        ~ "X11_DFA_Harbour_p_WS",
      
      # Guild 09 Fish Predators
      grepl("chinabundsnakefall", indicator)  ~ "X09_DFA_ChinAbundSnakeFall",
      grepl("hakeage5plus", indicator)        ~ "X09_DFA_HakeAge5Plus",
      
      # Guild 08 Birds
      grepl("08\\.predbirdncc|corm", indicator)~ "X08_DFA_DC_corm_3_WS",
      grepl("gull", indicator)                ~ "X08_Large_gulls_7_WS",
      
      # Guild 15 AK Fish Predators
      grepl("arrowtooth", indicator)          ~ "X15_ArrowtoothFlounderBiomass_predAK",
      grepl("sablefishbiomass", indicator)    ~ "X15_sablefishBiomass_predAK",
      grepl("sablefishrecruitment", indicator)~ "X15_sablefishRecruitment_predAK",
      grepl("pacificcod", indicator)          ~ "X15_PacificCodBiomass_predAK",
      grepl("spinydogfishgoa", indicator)     ~ "X15_spinyDogfishGoA_predAK",
      grepl("halibutbiomassage8plus", indicator)~ "X15_halibutBiomassAge8plus_2yrLead_predAK",
      grepl("sleepershark", indicator)        ~ "X15_DFA_sleeperSharks",
      grepl("salmonsharkgoa", indicator)      ~ "X15_salmonSharkGoA_predAK",
      
      TRUE ~ short_name_lower
    ),
    
    # --- Robust Helper logic for Prey Lags ---
    # Catch any string starting with "10" for NCC region
    is_adult_predator = (region == "NCC" & grepl("^10", latest_guild)),
    
    # --- Map Altprey 1 Data Column ---
    altprey1_data_col = case_when(
      altprey1 == "chin"         ~ ifelse(is_adult_predator, "X09_DFA_ChinAbundSnakeFall_2yrlead", "X09_DFA_ChinAbundSnakeFall"),
      altprey1 == "hake"         ~ ifelse(is_adult_predator, "X09_DFA_HakeAge5Plus_2yrlead", "X09_DFA_HakeAge5Plus"),
      altprey1 == "eulachon"     ~ ifelse(is_adult_predator, "eulachon_during_chinook_2yrlead", "eulachon_during_chinook"),
      altprey1 == "krill"        ~ ifelse(region == "AK", "X12_egoa.krill", ifelse(is_adult_predator, "X12_DFA_biomassEuphShelfSum_2yrlead", "X12_DFA_biomassEuphShelfSum")),
      altprey1 == "market_squid" ~ ifelse(is_adult_predator, "X04_marketsquid_GAM_2yrlead", "X04_marketsquid_GAM"),
      altprey1 == "pollock"      ~ ifelse(is_adult_predator, "X13_pollock_age1plus_2yrlead", "X13_pollock_age1plus"),
      altprey1 == "herring"      ~ ifelse(region == "AK", "X13_stka_herr_matbiom", "X05_DFA_abundSardine"),
      TRUE ~ NA_character_
    ),
    
    # --- Map Altprey 2 Data Column ---
    altprey2_data_col = case_when(
      altprey2 == "anchovy"     ~ ifelse(is_adult_predator, "X05_anchovy_GAM_2yrlead", "X05_anchovy_GAM"),
      altprey2 == "herring"     ~ ifelse(region == "AK", "X13_stka_herr_matbiom", ifelse(is_adult_predator, "X05_DFA_abundSardine_2yrlead", "X05_DFA_abundSardine")),
      altprey2 == "pink_salmon" ~ "X14_pinkSalmonNorthAmerica",
      altprey2 == "capelin"     ~ "X13_mid_il_capelin",
      TRUE ~ NA_character_
    ),
    
    # --- Map Altprey 3 Data Column ---
    altprey3_data_col = case_when(
      altprey3 == "eulachon"    ~ ifelse(is_adult_predator, "eulachon_during_chinook_2yrlead", "eulachon_during_chinook"),
      altprey3 == "anchovy"     ~ ifelse(is_adult_predator, "X05_anchovy_GAM_2yrlead", "X05_anchovy_GAM"),
      altprey3 == "herring"     ~ ifelse(region == "AK", "X13_stka_herr_matbiom", "X05_DFA_abundSardine"),
      altprey3 == "capelin"     ~ "X13_mid_il_capelin",
      altprey3 == "krill"       ~ ifelse(region == "AK", "X12_egoa.krill", ifelse(is_adult_predator, "X12_DFA_biomassEuphShelfSum_2yrlead", "X12_DFA_biomassEuphShelfSum")),
      TRUE ~ NA_character_
    )
  ) %>% 
  select(-is_adult_predator)
# ---------------------------------------------------------------------------
# 4. Filter strictly to columns present in sem_altprey_data_complete_1998_2021
# ---------------------------------------------------------------------------
valid_cols <- names(sem_complete_data)

all_pred_dfa_altprey_clean <- all_pred_dfa_altprey %>% 
  filter(pred_data_col %in% valid_cols) %>% 
  distinct(pred_data_col, .keep_all = TRUE)

all_pred_dfa_altprey_clean %>% filter(grepl("2yrLead",pred_data_col),region=="NCC")

# Export Final Lookup Table
write.csv(all_pred_dfa_altprey_clean, file.path(outputDir, "all_pred_dfa_altprey.csv"), row.names = FALSE)

cat(sprintf("Lookup table complete! %d distinct predator nodes mapped to complete 1998-2021 data.\n", nrow(all_pred_dfa_altprey_clean)))