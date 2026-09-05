
#Step 1: Master Name Crosswalk Generator Script
#Below is the updated, clean version of your crosswalk creation script (00_create_master_crosswalk.R). 
#It standardizes column names, cleans up year tags (_2025/_2026), maps the 17 DFAs and 48 stragglers/smoothed series, and 
#exports master_name_crosswalk.csv.
#to be used with utility function source(file.path(proj_dir, "00b_crosswalk_utility_fxn.r"))


# ==============================================================================
# Script: 00_create_master_crosswalk.R
# Purpose: Generate the Master Variable Name Crosswalk (LisaName vs Doug dfa_cols)
# ==============================================================================

library(tidyverse)
library(stringr)

# --- Configuration & Paths ---
proj_dir     <- getwd()
input_dir    <- file.path(proj_dir, "output", "DFA")
meta_dir     <- file.path(proj_dir, "metadata")
dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)

ranked_file  <- file.path(input_dir, "rankedIndicators.csv")

# Load Doug's ranked indicators file
if (!file.exists(ranked_file)) {
  stop("Missing required input: rankedIndicators.csv in ", input_dir)
}

rankedIndicators <- read.csv(ranked_file, row.names = NULL) %>%
  mutate(dfa_cols = paste0("X", DFAname))

# ------------------------------------------------------------------------------
# Construct Master Crosswalk Mapping
# ------------------------------------------------------------------------------
master_crosswalk <- rankedIndicators %>%
  mutate(
    # Clean off trailing _2025 / _2026 year tags
    clean_indicator_name = gsub("_2025|_2026", "", rankedIndicators),
    
    # Map Preferred Standardized Names (LisaName)
    LisaName = case_when(
      # --- Multi-Indicator DFA Overrides ---
      dfa_cols == "X01.ZooPreyNCC_JSOES_DFA1" ~ "X01_DFA_sumPreyOfPrey_planktonJun",
      dfa_cols == "X01.ZooPreyNCC_NHL_DFA1"   ~ "X01_DFA_NHL_krill",
      dfa_cols == "X02.CompJellyNCC_DFA1"     ~ "X02_DFA_SeaNettle",
      dfa_cols == "X03.FishPreyNCC_DFA1"      ~ "X03_DFA_comMurreDietHerrSard_Yaquina_NCC",
      dfa_cols == "X05.ForageFishNCC_DFA1"    ~ "X05_DFA_abundSardine",
      dfa_cols == "X06.Cond1NCC_DFA1"         ~ "X06_DFA_IGF_mu",
      dfa_cols == "X07.Cond2NCC_DFA1"         ~ "X07_DFA_cpue_IntSprJunHW",
      dfa_cols == "X08.PredBirdNCC_DFA1"      ~ "X08_DFA_DC_corm_3_WS",
      dfa_cols == "X09.PredFishNCC_DFA1"      ~ "X09_DFA_ChinAbundSnakeFall",
      dfa_cols == "X09.PredFishNCC_b_DFA1"    ~ "X09_DFA_HakeAge5Plus",
      dfa_cols == "X10.PredMammalNCC_DFA1"    ~ "X10_DFA_ssl.est.wholerange_2yrLead",
      dfa_cols == "X11.PredMammalSmolt_DFA1"  ~ "X11_DFA_Harbour_p_WS",
      dfa_cols == "X12.ZooPreyAK_DFA1"        ~ "X12_DFA_biomassEuphShelfSum",
      dfa_cols == "X13.FishPreyAK_DFA1"       ~ "X13_DFA_WGOA_DFA_midTrophic",
      dfa_cols == "X13.FishPreyAK_b_DFA1"     ~ "X13_DFA_WGOA_DFA_seabirds",
      dfa_cols == "X15.PredFishAK_b_DFA1"     ~ "X15_DFA_sleeperSharks",
      dfa_cols == "XSAR_DFA1"                 ~ "X16_SAR",
      
      # --- Single-Indicator Smoothed Nodes (Stragglers) ---
      TRUE ~ case_when(
        str_starts(clean_indicator_name, "X[0-9]{2}_") ~ clean_indicator_name,
        TRUE ~ paste0(substr(dfa_cols, 1, 3), "_", clean_indicator_name)
      )
    )
  ) %>%
  select(LisaName, dfa_cols, DFAname, guild, rankedIndicators) %>%
  as_tibble()

# ------------------------------------------------------------------------------
# Verification Output & Save
# ------------------------------------------------------------------------------
cat("Generated Crosswalk with", nrow(master_crosswalk), "total indicators.\n")
cat("DFAs mapped:", sum(grepl("DFA1", master_crosswalk$dfa_cols)), "\n")
cat("Stragglers mapped:", sum(!grepl("DFA1", master_crosswalk$dfa_cols)), "\n\n")

# Export Master Crosswalk
write.csv(master_crosswalk, file.path(meta_dir, "master_name_crosswalk.csv"), row.names = FALSE)
cat("Master Crosswalk saved to:", file.path(meta_dir, "master_name_crosswalk.csv"), "\n")