# ---------------------------------------------------------------------------
# 2. Map exact data source column names into nccpred
# ---------------------------------------------------------------------------
nccpred_dfa <- read.csv(file.path(outputDir, "nccpred_dfa_altprey_assigned.csv"), row.names = NULL)

nccpred_dfa_prey_mapped <- nccpred_dfa %>%
  mutate(

    # Map exact source column name for altprey1
    altprey1_data_col = case_when(
      altprey1 == "chin"         & !is_2yrlead ~ "X09.PredFishNCC_DFA1",
      altprey1 == "chin"         &  is_2yrlead ~ "X09.PredFishNCC_DFA1_2yrlead",
      
      altprey1 == "hake"         & !is_2yrlead ~ "X09.PredFishNCC_b_DFA1",
      altprey1 == "hake"         &  is_2yrlead ~ "X09.PredFishNCC_b_DFA1_2yrlead",
      
      altprey1 == "krill"        & !is_2yrlead ~ "X12.ZooPreyAK_DFA1",
      altprey1 == "krill"        &  is_2yrlead ~ "X12.ZooPreyAK_DFA1_2yrlead",
      
      altprey1 == "market_squid" & !is_2yrlead ~ "X04.CompFishNCC_smoothed",
      altprey1 == "market_squid" &  is_2yrlead ~ "X04.CompFishNCC_smoothed_2yrlead",
      
      altprey1 == "anchovy"      & !is_2yrlead ~ "X05.ForageFishNCC_0_smoothed",
      altprey1 == "anchovy"      &  is_2yrlead ~ "X05.ForageFishNCC_0_smoothed_2yrlead",
      
      altprey1 == "herring"      & !is_2yrlead ~ "X13.FishPreyAK_b_DFA1",
      altprey1 == "herring"      &  is_2yrlead ~ "X13.FishPreyAK_b_DFA1_2yrlead",
      
      altprey1 == "eulachon"     & !is_2yrlead ~ "eulachon_during_chinook",
      altprey1 == "eulachon"     &  is_2yrlead ~ "eulachon_during_chinook_2yrlead",
      
      TRUE ~ NA_character_
    ),
    
    # Map exact source column name for altprey2
    altprey2_data_col = case_when(
      altprey2 == "salmonids"                  ~ NA_character_, # Focal salmonids set to NA
      
      altprey2 == "anchovy"      & !is_2yrlead ~ "X05.ForageFishNCC_0_smoothed",
      altprey2 == "anchovy"      &  is_2yrlead ~ "X05.ForageFishNCC_0_smoothed_2yrlead",
      
      altprey2 == "herring"      & !is_2yrlead ~ "X13.FishPreyAK_b_DFA1",
      altprey2 == "herring"      &  is_2yrlead ~ "X13.FishPreyAK_b_DFA1_2yrlead",
      
      altprey2 == "eulachon"     & !is_2yrlead ~ "eulachon_during_chinook",
      altprey2 == "eulachon"     &  is_2yrlead ~ "eulachon_during_chinook_2yrlead",
      
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(altprey2,altprey1)

unique(nccpred_dfa_prey_mapped$dfa_cols)
names(nccpred_dfa_prey_mapped)
head(nccpred_dfa_prey_mapped %>% filter(is_2yrlead=="TRUE"))
# Export updated nccpred dataset
write.csv(nccpred_dfa_prey_mapped, file.path(outputDir, "nccpred_dfa_prey_mapped.csv"), row.names = FALSE)

# Display mapping verification
print(as_tibble(nccpred_dfa_prey_mapped %>% 
                  select(constituent_predators, altprey1, altprey2,altprey1_data_col,  altprey2_data_col)), n = Inf)
