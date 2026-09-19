# ==============================================================================
# Script: 04_altprey_lookup_table.R
# Purpose: Build AltPrey lookup table mapping predator indicators to exact 
#          LisaNames matching master dataset headers (excluding killer whales),
#          and export SEM modeling subsets.
# Outputs: 
#   - metadata/all_pred_dfa_altprey.csv
#   - metadata/sem_data_all_1998_2021.csv
#   - metadata/sem_data_noNA_1998_2021.csv
# ==============================================================================

library(tidyverse)
library(janitor)

rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
outputDir <- file.path(rootdir, "metadata")
path      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"

# File paths
master_data_file <- file.path(outputDir, "sem_altprey_data_complete_1998_2021.csv")
crosswalk_path   <- file.path(outputDir, "master_name_crosswalk.csv")
guild_file       <- file.path(path, "guildsWithExclude.csv")

if (!file.exists(master_data_file)) stop("Missing master SEM data file: ", master_data_file)
if (!file.exists(crosswalk_path))   stop("Missing crosswalk file: ", crosswalk_path)
if (!file.exists(guild_file))       stop("Missing guild file at: ", guild_file)

sem_data     <- read.csv(master_data_file, stringsAsFactors = FALSE)
dataset_cols <- names(sem_data)

# ------------------------------------------------------------------------------
# 1. LOAD GUILD ASSIGNMENTS (NCC + AK, EXCLUDING KILLER WHALES)
# ------------------------------------------------------------------------------
excluded_preds <- c(
  "canaryrockfish", "chilipepper", "spinydogfishbsai", 
  "salmonsharkbsai", "pacificspinydogfish", "iphc20halibut", 
  "loon", "sharkcatchgoa", "harbor_seal_cr_2025", "killer.whale", "killerwhale"
)

combined_guildfiles <- read.csv(guild_file, row.names = NULL, stringsAsFactors = FALSE) %>% 
  clean_names() %>% 
  rename(sem_name = se_mnode, short_name_lower = short_name) %>% 
  mutate(across(where(is.character), tolower)) %>% 
  filter(!is.na(latest_guild)) %>% 
  filter(!grepl(paste(excluded_preds, collapse = "|"), indicator))

# NCC Predators (Guilds 08, 09, 10, 11)
nccpred <- combined_guildfiles %>% 
  filter(grepl("08|09|10|11", latest_guild)) %>% 
  select(indicator, short_name_lower, sem_name, latest_guild) %>% 
  mutate(
    region = "NCC",
    altprey1 = case_when(
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

# AK Predators (Guild 15)
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

# ------------------------------------------------------------------------------
# 2. DYNAMIC CROSSWALK LOOKUP & DATASET HEADER MATCHING
# ------------------------------------------------------------------------------
master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

lookup_pred_base <- function(ind_raw, short_raw, guild_raw) {
  g_target <- tolower(trimws(guild_raw))
  i_target <- tolower(trimws(ind_raw))
  s_target <- tolower(trimws(short_raw))
  
  if (grepl("ssl_seak_pup_pred", i_target)) return("X11_ssl_seak_pup_pred")
  
  match_row <- master_crosswalk %>%
    filter(
      tolower(trimws(guild)) == g_target & 
        (str_detect(tolower(rankedIndicators), fixed(i_target)) |
           str_detect(tolower(rankedIndicators), fixed(s_target)))
    )
  
  if (nrow(match_row) > 0) return(match_row$LisaName[1])
  return(ind_raw)
}

# Resolves dataset column name (checks suffix first, falls back to base)
resolve_dataset_col <- function(base_name, suffix_str) {
  if (is.na(base_name) || base_name == "") return(NA_character_)
  
  candidate_suffixed <- paste0(base_name, "_", suffix_str)
  if (candidate_suffixed %in% dataset_cols) return(candidate_suffixed)
  if (base_name %in% dataset_cols)          return(base_name)
  
  # Match prefix in dataset
  matched <- dataset_cols[grep(paste0("^", base_name), dataset_cols)][1]
  if (!is.na(matched)) return(matched)
  
  return(candidate_suffixed)
}

get_prey_base <- function(prey_key, region_val) {
  case_when(
    prey_key == "chin"         ~ "X09_DFA_ChinAbundSnakeFall",
    prey_key == "hake"         ~ "X09_DFA_HakeAge5Plus",
    prey_key == "eulachon"     ~ "X05_eulachon_during_chinook",
    prey_key == "krill"        ~ ifelse(region_val == "AK", "X12_egoa_krill", "X12_DFA_biomassEuphShelfSum"),
    prey_key == "market_squid" ~ "X04_marketsquid_GAM",
    prey_key == "pollock"      ~ "X13_pollock_age1plus",
    prey_key == "herring"      ~ ifelse(region_val == "AK", "X13_sitkaHerring_EGoA", "X05_DFA_abundSardine"),
    prey_key == "anchovy"      ~ "X05_anchovy_GAM",
    prey_key == "pink_salmon"  ~ "X14_pinkSalmonNorthAmerica",
    prey_key == "capelin"      ~ "X13_mid_il_capelin",
    TRUE                       ~ NA_character_
  )
}

# Build base lookup mapping
all_preds_base <- all_preds_assigned %>%
  rowwise() %>%
  mutate(
    base_pred_col = lookup_pred_base(indicator, short_name_lower, latest_guild),
    is_x10_x15    = grepl("^10|^15", latest_guild)
  ) %>%
  ungroup()

# Generate adultyr rows for X10 & X15 predators
preds_adultyr <- all_preds_base %>%
  filter(is_x10_x15) %>%
  rowwise() %>%
  mutate(
    pred_data_col     = resolve_dataset_col(base_pred_col, "adultyr"),
    altprey1_data_col = resolve_dataset_col(get_prey_base(altprey1, region), "adultyr"),
    altprey2_data_col = resolve_dataset_col(get_prey_base(altprey2, region), "adultyr"),
    altprey3_data_col = resolve_dataset_col(get_prey_base(altprey3, region), "adultyr")
  ) %>%
  ungroup()

# Generate smoltyr rows for all predators
preds_smoltyr <- all_preds_base %>%
  rowwise() %>%
  mutate(
    pred_data_col     = resolve_dataset_col(base_pred_col, "smoltyr"),
    altprey1_data_col = resolve_dataset_col(get_prey_base(altprey1, region), "smoltyr"),
    altprey2_data_col = resolve_dataset_col(get_prey_base(altprey2, region), "smoltyr"),
    altprey3_data_col = resolve_dataset_col(get_prey_base(altprey3, region), "smoltyr")
  ) %>%
  ungroup()

all_pred_dfa_altprey <- bind_rows(preds_adultyr, preds_smoltyr) %>%
  select(latest_guild, region, pred_data_col, altprey1_data_col, altprey2_data_col, altprey3_data_col) %>%
  distinct() %>%
  group_by(latest_guild, pred_data_col, altprey1_data_col, altprey2_data_col) %>%
  filter(n() == 1 | !is.na(altprey3_data_col)) %>%
  ungroup() %>%
  arrange(latest_guild, pred_data_col)

out_lookup_path <- file.path(outputDir, "all_pred_dfa_altprey.csv")
write.csv(all_pred_dfa_altprey, out_lookup_path, row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. DATASET ORGANIZATION & SANITY CHECK
# ------------------------------------------------------------------------------
custom_sort_x <- function(cols) {
  cols[order(
    as.numeric(gsub("^X(\\d+).*", "\\1", cols)),
    cols
  )]
}

all_other_cols        <- setdiff(names(sem_data), "Year")
sorted_indicator_cols <- custom_sort_x(all_other_cols)
target_col_order      <- c("Year", sorted_indicator_cols)

sem_data_all  <- sem_data %>% select(all_of(target_col_order))
complete_cols <- names(sem_data_all)[colSums(is.na(sem_data_all)) == 0]
sem_data_noNA <- sem_data_all %>% select(all_of(complete_cols))

out_all_path  <- file.path(outputDir, "sem_data_all_1998_2021.csv")
out_noNA_path <- file.path(outputDir, "sem_data_noNA_1998_2021.csv")

write.csv(sem_data_all,  out_all_path,  row.names = FALSE)
write.csv(sem_data_noNA, out_noNA_path, row.names = FALSE)

# Pre-flight SANITY CHECK
required_cols <- unique(c(
  all_pred_dfa_altprey$pred_data_col,
  all_pred_dfa_altprey$altprey1_data_col,
  all_pred_dfa_altprey$altprey2_data_col,
  all_pred_dfa_altprey$altprey3_data_col
))
required_cols <- required_cols[!is.na(required_cols) & required_cols != ""]

missing_cols <- setdiff(required_cols, names(sem_data_all))

cat("\n======================================================\n")
cat("          LOOKUP TABLE & SANITY CHECK RESULTS         \n")
cat("======================================================\n")
cat("Exported Lookup Table:", out_lookup_path, "\n")
if (length(missing_cols) == 0) {
  cat("SUCCESS: All", length(required_cols), "referenced lookup columns exist in the dataset!\n")
} else {
  cat("WARNING: The following", length(missing_cols), "columns are missing from the dataset:\n")
  print(missing_cols)
}
cat("======================================================\n")