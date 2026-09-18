



# ==============================================================================
# Script: altprey_01_setup_altprey_lookup_table.r
# Purpose: Build unified SEM dataset and construct AltPrey lookup table 
#          mapping raw indicators to exact LisaNames via dual guild/indicator matching.
# Output: copilot/outputs_altprey/sem_altprey_data_complete_1998_2021.csv
#         copilot/outputs_altprey/all_pred_dfa_altprey.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(janitor)

outputDir <- file.path("copilot/outputs_altprey")
metaDir   <- file.path("Rcode_for_paper/metadata")

dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
dir.create(metaDir,   showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load Extended Dataset & Append Supplementary Indicators
# ---------------------------------------------------------------------------
clusData_led <- read.csv(
  file.path(outputDir, "clusDataDFA_wide_1998_2021_altprey_extended_smoltyear_shifted_LisaNames.csv"), 
  stringsAsFactors = FALSE
)

eulachon_scaled <- read.csv("copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv", stringsAsFactors = FALSE) %>%
  select(Year, eulachon_during_chinook, eulachon_annual_scaled) %>%
  mutate(
    eulachon_during_chinook_2yrlead = dplyr::lead(eulachon_during_chinook, 2),
    eulachon_annual_scaled_2yrlead  = dplyr::lead(eulachon_annual_scaled, 2)
  )

ak_yr_selected <- read.csv("data_Lisa/ak_yr.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  rename(Year = year) %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  select(
    Year,
    X11_ssl_seak_pup_pred = ssl_seak_pup_pred,
    X12_egoa.krill        = secm_euph_dens,
    X13_mid_il_capelin    = mid_il_capelin,
    X21_sst_egoa_junjulaug = sst_egoa_coastwatch_junjulaug
  )

sem_altprey_data <- clusData_led %>%
  left_join(eulachon_scaled, by = "Year") %>%
  left_join(ak_yr_selected, by = "Year") %>%
  arrange(Year)

df_1998_2021 <- sem_altprey_data %>% filter(Year >= 1998 & Year <= 2021)

# RETAIN ALL COLUMNS (Do not drop columns containing shifted edge NAs)
valid_data_cols <- names(df_1998_2021)[colSums(!is.na(df_1998_2021)) > 0]

sem_complete_data <- df_1998_2021 %>%
  select(all_of(valid_data_cols)) %>%
  mutate(across(-Year, ~ as.numeric(scale(.x))))

write.csv(sem_complete_data, file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)
write.csv(sem_complete_data, file.path(metaDir,   "sem_altprey_data_complete_1998_2021.csv"), row.names = FALSE)
write.csv(df_1998_2021,      file.path(outputDir, "sem_data_1998_2021.csv"), row.names = FALSE)
write.csv(df_1998_2021,      file.path(metaDir,   "sem_data_1998_2021.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Load Guild Assignments (NCC + AK)
# ---------------------------------------------------------------------------
excluded_preds <- c(
  "canaryrockfish", "chilipepper", "spinydogfishbsai", 
  "salmonsharkbsai", "pacificspinydogfish", "iphc20halibut", 
  "loon", "sharkcatchgoa", "harbor_seal_cr_2025"
)

combined_guildfiles <- read.csv(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/guildsWithExclude.csv", 
  row.names = NULL, stringsAsFactors = FALSE
) %>% 
  clean_names() %>% 
  rename(sem_name = se_mnode, short_name_lower = short_name) %>% 
  mutate(across(where(is.character), tolower)) %>% 
  filter(!is.na(latest_guild)) %>% 
  filter(!grepl(paste(excluded_preds, collapse = "|"), indicator))

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
# 3. Dynamic Crosswalk Mapping & Post-Processing Deduplication
# ---------------------------------------------------------------------------
crosswalk_path <- file.path(metaDir, "master_name_crosswalk.csv")
master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

lookup_pred_data_col <- function(ind_raw, short_raw, guild_raw) {
  g_target <- tolower(trimws(guild_raw))
  i_target <- tolower(trimws(ind_raw))
  s_target <- tolower(trimws(short_raw))
  
  match_row <- master_crosswalk %>%
    filter(
      tolower(trimws(guild)) == g_target & 
        (str_detect(tolower(rankedIndicators), fixed(i_target)) |
           str_detect(tolower(rankedIndicators), fixed(s_target)))
    )
  
  if (nrow(match_row) > 0) {
    return(match_row$LisaName[1])
  }
  
  if (grepl("ssl_seak_pup_pred", ind_raw, ignore.case = TRUE)) return("X11_ssl_seak_pup_pred")
  return(ind_raw)
}

all_preds_base <- all_preds_assigned %>%
  rowwise() %>%
  mutate(
    is_adult_predator = grepl("^10|^15", latest_guild),
    pred_data_col     = lookup_pred_data_col(indicator, short_name_lower, latest_guild),
    
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
  ) %>%
  ungroup()

all_preds_smoltYear <- all_preds_base %>%
  filter(is_adult_predator) %>%
  mutate(
    pred_data_col = ifelse(
      grepl("_2yrLead|_2yrlead", pred_data_col),
      str_replace(pred_data_col, "_2yrLead|_2yrlead", "_smoltYear"),
      paste0(pred_data_col, "_smoltYear")
    ),
    
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

all_pred_dfa_altprey <- bind_rows(all_preds_base, all_preds_smoltYear) %>%
  select(latest_guild, region, pred_data_col, altprey1_data_col, altprey2_data_col, altprey3_data_col) %>%
  distinct() %>%
  group_by(latest_guild, pred_data_col, altprey1_data_col, altprey2_data_col) %>%
  filter(n() == 1 | !is.na(altprey3_data_col)) %>%
  ungroup() %>%
  arrange(latest_guild, pred_data_col)

out_lookup_path  <- file.path(outputDir, "all_pred_dfa_altprey.csv")
meta_lookup_path <- file.path(metaDir,   "all_pred_dfa_altprey.csv")

write.csv(all_pred_dfa_altprey, out_lookup_path,  row.names = FALSE)
write.csv(all_pred_dfa_altprey, meta_lookup_path, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Pre-Flight Verification Report
# ---------------------------------------------------------------------------
sem_data <- read.csv(file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv"), stringsAsFactors = FALSE)
lookup   <- read.csv(out_lookup_path, stringsAsFactors = FALSE)

required_cols <- unique(c(
  lookup$pred_data_col,
  lookup$altprey1_data_col,
  lookup$altprey2_data_col,
  lookup$altprey3_data_col
))
required_cols <- required_cols[!is.na(required_cols) & required_cols != ""]

missing_cols <- setdiff(required_cols, names(sem_data))
present_cols <- intersect(required_cols, names(sem_data))

na_summary <- colSums(is.na(sem_data[, present_cols, drop = FALSE]))
cols_with_nas <- names(na_summary[na_summary > 0])

cat("\n======================================================\n")
cat("         PRE-FLIGHT SANITY CHECK RESULTS             \n")
cat("======================================================\n")

if (length(missing_cols) == 0) {
  cat("SUCCESS: All", length(required_cols), "referenced lookup columns exist in sem_altprey_data_complete_1998_2021.csv!\n")
} else {
  cat("WARNING: The following", length(missing_cols), "columns are missing from the dataset:\n")
  print(missing_cols)
}

if (length(cols_with_nas) == 0) {
  cat("SUCCESS: Zero NAs detected across all referenced lookup columns (1998–2021)!\n")
} else {
  cat("WARNING: NAs found in the following dataset columns:\n")
  print(na_summary[cols_with_nas])
}

cat("======================================================\n\n")

#SEM DATASET ORGANIZATION -------------

# ---------------------------------------------------------------------------
# Load Input Datasets
# ---------------------------------------------------------------------------
library(dplyr)
proj_dir     <- paste0(getwd(),"/Rcode_for_paper")
sem_data <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_complete_1998_2021.csv"), stringsAsFactors = FALSE)
altprey_lookup <- read.csv(file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv"), stringsAsFactors = FALSE)


# Filter explicitly to modeling window 1998-2021
sem_data_1998_2021 <- sem_data %>%
  filter(Year >= 1998 & Year <= 2021) %>%
  arrange(Year)

# ---------------------------------------------------------------------------
# Sort Column Names (Year first, followed by X01...X21 ordered)
# ---------------------------------------------------------------------------
all_other_cols <- setdiff(names(sem_data_1998_2021), "Year")

# Custom sorting helper so X01, X02, ..., X10, X21 sort numerically rather than strictly alphabetically
custom_sort_x <- function(cols) {
  cols[order(
    # Extract leading prefix number (e.g. X01 -> 1, X21 -> 21)
    as.numeric(gsub("^X(\\d+).*", "\\1", cols)),
    # Secondary alphabetical sort for suffixes/names within same prefix
    cols
  )]
}

sorted_indicator_cols <- custom_sort_x(all_other_cols)
target_col_order <- c("Year", sorted_indicator_cols)

# ---------------------------------------------------------------------------
# Version 1: All Data (Ordered Column Schema)
# ---------------------------------------------------------------------------
sem_data_all <- sem_data_1998_2021 %>%
  select(all_of(target_col_order))

# ---------------------------------------------------------------------------
# Version 2: Complete Cases Only (Zero NAs across 1998-2021)
# ---------------------------------------------------------------------------
complete_cols <- names(sem_data_all)[colSums(is.na(sem_data_all)) == 0]

sem_data_noNA <- sem_data_all %>%
  select(all_of(complete_cols))

# ---------------------------------------------------------------------------
# Export Both Versions to Data Output Directory
# ---------------------------------------------------------------------------
out_all_path  <- file.path(proj_dir, "metadata", "sem_data_all_1998_2021.csv")
out_noNA_path <- file.path(proj_dir, "metadata", "sem_data_noNA_1998_2021.csv")

write.csv(sem_data_all,  out_all_path,  row.names = FALSE)
write.csv(sem_data_noNA, out_noNA_path, row.names = FALSE)

# Print Summary
cat("\n======================================================\n")
cat("            SEM DATASET ORGANIZATION COMPLETE         \n")
cat("======================================================\n")
cat("Full dataset exported:", out_all_path, "\n")
cat("  - Rows:", nrow(sem_data_all), "| Columns:", ncol(sem_data_all), "\n\n")

cat("No-NA dataset exported:", out_noNA_path, "\n")
cat("  - Rows:", nrow(sem_data_noNA), "| Columns:", ncol(sem_data_noNA), "\n")
cat("======================================================\n\n")
