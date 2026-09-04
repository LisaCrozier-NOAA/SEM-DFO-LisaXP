library(dplyr)
library(stringr)
library(tidyr)

# ---------------------------------------------------------------------------
# 1. Load Files
# ---------------------------------------------------------------------------

#update Doug's DFAs and list of which indicators are in which DFAs

      source("copilot/rcode_alt_prey/doug_01_analyzeIndicators_extend_to_1996_2025_datwide.R")
      source("copilot/rcode_alt_prey/doug_02_createDFA.R")
      
      clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))
      all_dfa_cols     <- setdiff(names(clusDataDFA_wide), "Year")
      
      rankedIndicators <- read.csv(file.path(outputDir, "rankedIndicators.csv"), row.names = NULL) %>%
        mutate(dfa_cols=paste0("X",DFAname))%>% 
        relocate(dfa_cols,.after=DFAname) %>%
        filter(grepl("08|09|10|11",DFAname)) %>%
        select(dfa_cols,rankedIndicators)
      head(rankedIndicators)
      
      rankedIndicators_clean[1:2,"rankedIndicators_lower"]


nccpred_mapped <- read.csv(file.path(outputDir, "nccpred_mapped_source_names.csv"), row.names = NULL) %>%
  rename(short_name_lower=short_name) %>%
  select(short_name_lower,altprey1,altprey2,is_2yrlead,altprey1_data_col,altprey2_data_col) 

nccpred_mapped[1:2,]


# 3. Match via grepl (case-insensitive substring search inside the collapsed string)
nccpred_dfa_matched <- nccpred_mapped %>%
  # Join all combinations to test each short_name against each collapsed DFA string
  left_join(rankedIndicators, by = character(), relationship = "many-to-many") %>%
  # Keep rows where short_name_lower is found inside the raw collapsed string
  # (or keep as NA if an indicator isn't part of any DFA in rankedIndicators)
  filter(
    is.na(rankedIndicators) | 
      mapply(function(pattern, text) grepl(pattern, text, ignore.case = TRUE), 
             short_name_lower, rankedIndicators)
  ) %>%
  relocate(dfa_cols,.before=short_name_lower)

# Inspect top rows
head(nccpred_dfa_matched, 5)
names(nccpred_dfa_matched)

nccpred_dfa_matched[1,]
print(as_tibble(nccpred_dfa_matched %>% select(dfa_cols,short_name_lower,rankedIndicators)),n=Inf)


# ---------------------------------------------------------------------------
# 4. Filter DFAs with Conflicting Alternate Prey Assignments
# ---------------------------------------------------------------------------
prey_conflicts <- nccpred_dfa_matched %>%
  group_by(dfa_cols) %>%
  summarise(
    num_predators       = n(),
    predators_in_dfa    = paste(unique(short_name_lower), collapse = ", "),
    altprey1_choices    = paste(unique(altprey1[!is.na(altprey1)]), collapse = " | "),
    altprey2_choices    = paste(unique(altprey2[!is.na(altprey2)]), collapse = " | "),
    altprey1_conflict   = n_distinct(altprey1[!is.na(altprey1)]) > 1,
    altprey2_conflict   = n_distinct(altprey2[!is.na(altprey2)]) > 1,
    .groups = "drop"
  ) %>%
  filter(altprey1_conflict | altprey2_conflict)


prey_conflicts %>% select(dfa_cols,num_predators,altprey1_choices,altprey2_choices,predators_in_dfa)

nccpred_dfa_matched %>% 
  filter(dfa_cols=="X09.PredFishNCC_DFA1") %>% 
  select(dfa_cols,short_name_lower,altprey1,altprey2)

nccpred_dfa_matched %>% 
  filter(dfa_cols=="X09.PredFishNCC_b_DFA1") %>% 
  select(dfa_cols,short_name_lower,altprey1,altprey2)

nccpred_dfa_matched %>% 
  filter(dfa_cols=="X10.PredMammalNCC_DFA1") %>% 
  select(dfa_cols,short_name_lower,altprey1,altprey2)


library(dplyr)
library(stringr)

# ---------------------------------------------------------------------------
# 1. Update Sub-Groups for the 3 Conflicting DFAs------
# ---------------------------------------------------------------------------
nccpred_resolved <- nccpred_dfa_matched %>%
  mutate(
    # Recode dfa_cols and altprey for the sub-groups based on short_name_lower
    dfa_cols = case_when(
      # --- DFA 1: X09.PredFishNCC_DFA1 ---
      dfa_cols == "X09.PredFishNCC_DFA1" & grepl("halibut|lingcod|longnoseskate", short_name_lower, ignore.case = TRUE) ~ 
        "X09.PredFishNCC_DFA1_halibut_lingcod",
      dfa_cols == "X09.PredFishNCC_DFA1" & grepl("chin|blackrockfish|jackmack_ncc|stripetailrockfish_2025", short_name_lower, ignore.case = TRUE) ~ 
        "X09.PredFishNCC_DFA1_chin_blackrockfish",
      
      # --- DFA 2: X09.PredFishNCC_b_DFA1 ---
      dfa_cols == "X09.PredFishNCC_b_DFA1" & grepl("dogfish|bigskate", short_name_lower, ignore.case = TRUE) ~ 
        "X09.PredFishNCC_b_DFA1_dogfish",
      dfa_cols == "X09.PredFishNCC_b_DFA1" & grepl("hake|jackmackerel_gam_2025|pacificmackerel_gam_2025|darkblotchedrockfish|greestripedrockfish", short_name_lower, ignore.case = TRUE) ~ 
        "X09.PredFishNCC_b_DFA1_hake_mack",
      
      # --- DFA 3: X10.PredMammalNCC_DFA1 ---
      dfa_cols == "X10.PredMammalNCC_DFA1" & grepl("allsealionsbonn", short_name_lower, ignore.case = TRUE) ~ 
        "X10.PredMammalNCC_DFA1_bonn",
      dfa_cols == "X10.PredMammalNCC_DFA1" & grepl("steller|ssl.est.wholerange|allsealionsemb", short_name_lower, ignore.case = TRUE) ~ 
        "X10.PredMammalNCC_DFA1_ssl",
      
      TRUE ~ dfa_cols
    ),
    
    # Assign specific altprey1 & altprey2 choices to the new sub-DFAs
    altprey1 = case_when(
      dfa_cols == "X10.PredMammalNCC_DFA1_cslbonn"     ~ "eulachon",
      dfa_cols == "X10.PredMammalNCC_DFA1_ssl"         ~ "hake",
      dfa_cols == "X09.PredFishNCC_DFA1_halibut_lingcod" ~ "hake",
      dfa_cols == "X09.PredFishNCC_b_DFA1_dogfish"      ~ "hake",
      dfa_cols == "X09.PredFishNCC_DFA1_chin_blackrockfish" ~ "krill",
      dfa_cols == "X09.PredFishNCC_b_DFA1_hake_mack"     ~ "krill",
      TRUE ~ altprey1
    ),
    
    altprey2 = case_when(
      dfa_cols == "X09.PredFishNCC_DFA1_halibut_lingcod" ~ "herring",
      dfa_cols == "X09.PredFishNCC_DFA1_chin_blackrockfish" ~ "anchovy",
      dfa_cols == "X09.PredFishNCC_b_DFA1_dogfish"      ~ "herring",
      dfa_cols == "X09.PredFishNCC_b_DFA1_hake_mack"     ~ "anchovy",
      dfa_cols == "X10.PredMammalNCC_DFA1_cslbonn"     ~ "anchovy",
      dfa_cols == "X10.PredMammalNCC_DFA1_ssl"         ~ "herring",
      TRUE ~ altprey2
    ),
    altprey2_data_col = case_when(
      dfa_cols == "X10.PredMammalNCC_DFA1_cslbonn"     ~ "X05.ForageFishNCC_0_smoothed_2yrlead",
      TRUE ~ altprey2_data_col
      
  ))

nccpred_resolved %>% select(-rankedIndicators)
# ---------------------------------------------------------------------------
# 2. Collapse to Unique DFA Lines
# ---------------------------------------------------------------------------
nccpred_dfa_unique <- nccpred_resolved %>%
  filter(!is.na(dfa_cols)) %>%
  group_by(dfa_cols) %>%
  summarise(
    # Aggregate all constituent short names for reference
    constituent_predators = paste(unique(short_name_lower), collapse = " | "),
    altprey1              = unique(altprey1[!is.na(altprey1)])[1],
    altprey2              = unique(altprey2[!is.na(altprey2)])[1],
    is_2yrlead            = any(is_2yrlead, na.rm = TRUE),
    .groups = "drop"
  )

# View the final unconflicted DFA table
print(as_tibble(nccpred_dfa_unique), n = Inf)
