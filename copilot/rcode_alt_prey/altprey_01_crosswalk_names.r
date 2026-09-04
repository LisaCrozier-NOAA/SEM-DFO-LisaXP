library(dplyr)
library(stringr)

# Define Output Directory
new_outputDir <- file.path("copilot/outputs_9")
dir.create(new_outputDir, showWarnings = FALSE, recursive = TRUE)

# Load rankedIndicators from output directory
rankedIndicators <- read.csv("copilot/outputs_8/rankedIndicators.csv", row.names = NULL) %>%
  mutate(dfa_cols = paste0("X", DFAname))
smoothed_indiv_time_series_doug<-read.csv(file.path("copilot/outputs_8/datWide_1996_2025_qualified.csv"), row.names = NULL)
names(smoothed_indiv_time_series_doug)

# ---------------------------------------------------------------------------
# Construct Master Name Crosswalk (65 rows matching rankedIndicators)
# ---------------------------------------------------------------------------
master_name_crosswalk <- rankedIndicators %>%
  mutate(
    # Clean off trailing _2025 / _2026 year tags from rankedIndicators string
    clean_indicator_name = gsub("_2025|_2026", "", rankedIndicators),
    
    # Assign Preferred Target Names (Lisaname)
    Lisaname = case_when(
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
  select(Lisaname,dfa_cols,guild,rankedIndicators) %>%
  # Convert explicitly to a tibble to prevent print.default parser crashes
  as_tibble()

# Print DFAs cleanly with n = Inf
print(
  master_name_crosswalk %>% 
    filter(grepl("DFA1", dfa_cols)),
  n = Inf
)
#17 DFAs, including salmon

print(
  master_name_crosswalk %>% 
    filter(!grepl("DFA1", dfa_cols)),
  n = Inf
)

# Export to outputs_9
write.csv(master_name_crosswalk, file.path(new_outputDir, "master_name_crosswalk.csv"), row.names = FALSE)
