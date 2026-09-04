#Make a new lookup table with new DFA names since I re-ran Doug's script with a different threshold

outputDir <- file.path( "copilot/outputs_8")

rankedIndicators <- read.csv(file.path(outputDir, "rankedIndicators.csv"), row.names = NULL) %>%
  mutate(dfa_cols=paste0("X",DFAname)) %>% 
  relocate(dfa_cols,.after=DFAname) %>%
#  filter(grepl("12|13|14|15",DFAname)) %>%
  select(dfa_cols,rankedIndicators)

rankedIndicators[1,]
head(all_pred_dfa_altprey)

all_pred_dfa_altprey %>% filter(is_2yrlead==TRUE) %>% select(1:2)
all_pred_dfa_altprey %>% filter(is_2yrlead==FALSE) %>% select(1:2)

c(unique(all_pred_dfa_altprey %>% filter(is_2yrlead==TRUE) %>% pull(altprey1_data_col)),
unique(all_pred_dfa_altprey %>% filter(is_2yrlead==TRUE) %>% pull(altprey2_data_col)))

# all_pred_dfa_altprey %>%
#   mutate(Lisanames=)

head(var_lookup)


library(dplyr)
library(stringr)

# Read, create dfa_cols, and extract the 1st listed indicator into Lisaname
rankedIndicators <- read.csv(file.path(outputDir, "rankedIndicators.csv"), row.names = NULL) %>%
  mutate(Lisa_rankedIndicators = case_when(
    DFAname=="13.FishPreyAK_DFA1" ~  "pollockBiomassAIage1plus_predAK_2026 - WGOA_DFA_midTrophic_2026 - WGOA_DFA_lowerTrophic_2026 - rhinoAuk_EGoA - sitkaHerring_EGoA",
    DFAname=="13.FishPreyAK_b_DFA1" ~  "neg_capelin_WGoA - WGOA_DFA_seabirds_2026 - ammod_EAI - estAbundHerringRecruits - gadid_EAI",
.default = rankedIndicators
  )
    


  
  
  mutate(
    dfa_cols = paste0("X", DFAname),
    # Extract the first item before the " - " delimiter
    first_indicator = trimws(sapply(strsplit(as.character(rankedIndicators), " - "), `[`, 1)),
    # Combine dfa_cols with the first indicator name
    Lisaname = paste0(dfa_cols, "_", first_indicator)
  ) %>%
  select(-first_indicator) # Remove temporary column

# Inspect the result
head(rankedIndicators %>% select(DFAname, dfa_cols, Lisaname, rankedIndicators), 5)


write.csv(rankedIndicators,)
write.csv(rankedIndicators, file.path(outputDir, "rankedIndicators_Lisaname.csv"), row.names = FALSE)



library(dplyr)
library(tidyr)
library(stringr)

# ---------------------------------------------------------------------------
# 1. Identify all base prey data columns that require 2-year leads
# ---------------------------------------------------------------------------
prey_lead_cols <- all_pred_dfa_altprey %>%
  filter(is_2yrlead == TRUE) %>%
  select(altprey1_data_col, altprey2_data_col, altprey3_data_col) %>%
  pivot_longer(cols = everything(), values_to = "data_col") %>%
  filter(!is.na(data_col) & data_col != "") %>%
  # Remove existing '_2yrlead' suffix to ensure we match the base columns
  mutate(base_col = gsub("_2yrlead$", "", data_col)) %>%
  pull(base_col) %>%
  unique()

print(prey_lead_cols)
# ---------------------------------------------------------------------------
# 2. Extract matching rows from rankedIndicators and build _2yrlead rows
# ---------------------------------------------------------------------------
# Retrieve base rows matching the lead prey columns (checking dfa_cols or Lisaname)
base_prey_rows <- rankedIndicators %>%
  filter(dfa_cols %in% prey_lead_cols | Lisaname %in% prey_lead_cols)

# Create the new 2-year lead rows with updated dfa_cols, Lisaname, and metadata
new_lead_rows <- base_prey_rows %>%
  mutate(
    dfa_cols     = paste0(dfa_cols, "_2yrlead"),
    Lisaname     = paste0(Lisaname, "_2yrlead"),
    DFAname      = paste0(DFAname, "_2yrlead"),
    # Append lead suffix to the collapsed indicator manifest text if present
    rankedIndicators = if_else(
      !is.na(rankedIndicators) & rankedIndicators != "",
      paste0(rankedIndicators, "_2yrlead"),
      rankedIndicators
    )
  )

# ---------------------------------------------------------------------------
# 3. Append new lead rows back into rankedIndicators to make the master lookup
# ---------------------------------------------------------------------------
rankedIndicators_master <- rankedIndicators %>%
  bind_rows(new_lead_rows) %>%
  distinct(Lisaname, .keep_all = TRUE)

# Inspect the newly appended 2-year lead rows
print(
  rankedIndicators_master %>% 
    filter(grepl("_2yrlead$", Lisaname)) %>% 
    select(DFAname, dfa_cols, Lisaname)
)