library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)

# ==============================================================================
# Script: 02_analysis/07e_collect_DFA_constituents.r
# Step 1: Find raw data for the factors that were used to make DFAs 
# ==============================================================================


# ==============================================================================
# 0. DIRECTORY & DATA SETUP--------
# ==============================================================================
proj_dir      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
all_data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/allData.rds"
cw_path       <- file.path(proj_dir, "metadata/crosswalk_names.csv")
alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
lookup_path   <- file.path(proj_dir, "Routput_for_paper/data/guild_predator_lookup_table_24yr.csv")
fn_path       <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/functions/sem_output_summary_fxn.r"
save_rdata    <- file.path(proj_dir, "Routput_for_paper/models_DAG1B_latent_stages.RData")

if (file.exists(fn_path)) source(fn_path)

all_data   <- readRDS(all_data_path)
cw_df      <- read.csv(cw_path, stringsAsFactors = FALSE)
df_24_raw  <- read.csv(alt_data_path, stringsAsFactors = FALSE) %>% filter(complete.cases(.))
lookup_df  <- read.csv(lookup_path, stringsAsFactors = FALSE)


# ==============================================================================
# 1. PARSE & EXTRACT 100% COMPLETE CONSTITUENTS (1998-2021)-----
# ==============================================================================
parse_cw_tokens <- function(node_key, cw) {
  row_match <- cw %>%
    filter(LisaName == node_key | DFAname == node_key | dfa_cols == node_key)
  raw_tokens <- trimws(unlist(strsplit(row_match$rankedIndicators[1], "\\s*-\\s*|\\s*,\\s*")))
  raw_tokens <- raw_tokens[raw_tokens != ""]
  unique(gsub("_[0-9]{4}$", "", raw_tokens))
}

raw_x05_tokens <- parse_cw_tokens("X05_DFA_abundSardine", cw_df)
raw_x09_tokens <- parse_cw_tokens("X09_DFA_HakeAge5Plus", cw_df)
raw_x12_tokens <- parse_cw_tokens("X12_DFA_EGOA_DFA_lowerTrophic", cw_df)
raw_x13_tokens <- parse_cw_tokens("X13_DFA_WGOA_DFA_midTrophic", cw_df)
raw_x15_tokens <- parse_cw_tokens("X15_DFA_halibutBiomassAge8plus_predAK_adultyr", cw_df)

extract_complete_constituents <- function(all_df, token_list, node_prefix) {
  df_wide <- all_df %>%
    mutate(
      year            = year(ymd(date)),
      shortName_clean = gsub("_[0-9]{4}$", "", shortName)
    ) %>%
    filter(
      year >= 1998, year <= 2021,
      yearType == "smolt",
      shortName_clean %in% token_list
    ) %>%
    select(year, shortName = shortName_clean, finalVal) %>%
    distinct(year, shortName, .keep_all = TRUE) %>%
    pivot_wider(names_from = shortName, values_from = finalVal) %>%
    arrange(year)
  
  cand_cols     <- setdiff(names(df_wide), "year")
  col_nas       <- sapply(cand_cols, function(col) sum(is.na(df_wide[[col]])))
  complete_cols <- names(col_nas[col_nas == 0])
  
  df_wide %>%
    select(year, all_of(complete_cols)) %>%
    rename_with(~ paste0(node_prefix, .x), -year)
}

df_x05_comp <- extract_complete_constituents(all_data, raw_x05_tokens, "X05_")
df_x09_comp <- extract_complete_constituents(all_data, raw_x09_tokens, "X09_")
df_x12_comp <- extract_complete_constituents(all_data, raw_x12_tokens, "X12_")
df_x13_comp <- extract_complete_constituents(all_data, raw_x13_tokens, "X13_")
df_x15_comp <- extract_complete_constituents(all_data, raw_x15_tokens, "X15_")

dfa_constituents <- list(df_x05_comp, df_x09_comp, df_x12_comp, df_x13_comp, df_x15_comp) %>%
  reduce(full_join, by = "year") %>%
  arrange(year)

head(dfa_constituents)

write.csv(dfa_constituents, file.path(proj_dir, "metadata/sem_dfa_constituents.csv"),row.names = FALSE)
