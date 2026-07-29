#To turn your detailed ecological diagnostics into a clean, publication-ready summary, we can write an R script that automates the filters, applies your trophic hypotheses, determines whether a variable is a DFA or a straggler, programmatically pulls the key nearest neighbors, and applies your custom explanations.
#The Summary Construction Pipeline
#This script runs the entire analysis in one clean pipe. 
#It dynamically uses your find_mds_neighbors function behind the scenes to extract the qualifying neighbors ($r > 0.5$ or $r < -0.2$) for every target variable and collapses them into a single scannable cell.


library(tidyverse)

# --- 1. Helper Function to Extract Neighbors by Correlation Direction ---
get_neighbor_summary <- function(target_var, mds_data, direction = "pos") {
  neighbors <- tryCatch({
    find_mds_neighbors(target_var, mds_data, n_neighbors = 10)
  }, error = function(e) NULL)
  
  if (is.null(neighbors) || nrow(neighbors) == 0) return("")
  
  # Filter based on the requested direction and thresholds
  if (direction == "pos") {
    filtered <- neighbors %>% filter(Pairwise_Cor_With_Target > 0.5)
  } else {
    filtered <- neighbors %>% filter(Pairwise_Cor_With_Target < -0.2)
  }
  
  if (nrow(filtered) == 0) return("")
  
  # Format as "Name (Value)"
  filtered %>%
    mutate(Formatted = sprintf("%s (%s)", Lisaname, round(Pairwise_Cor_With_Target, 2))) %>%
    pull(Formatted) %>%
    paste(collapse = ", ")
}

# --- 2. Build the Ecological Synthesis Table ---
ecological_synthesis <- top100_indicator_behavior %>%
  # Filter 1: Keep only indicators selected more than 10 times in the top 100 models
  filter(Times_Selected_In_Top100 > 10) %>%
  rename(Times_In_Top100 = Times_Selected_In_Top100) %>%
  
  # Filter 2: Exclude DAG1C for mediators & NCC nodes
  filter(!(rhs %in% c("Abundance", "Growth", "PreyNCC", "PredNCC") & str_detect(model_id, "DAG1C"))) %>%
  
  # Add columns for Hypothesis, Type (DFA vs. Straggler), and Confirmation
  mutate(
      Node2 = case_when(
        grepl("X04|X05|X14", Lisaname) ~ "Competitor",
        TRUE ~ as.character(rhs) # Force character type consistency
      ),
    # Define expected ecological relationship
    Hyp = case_when(
      grepl("X04|X05|X14", Lisaname) ~ "Negative",
      rhs %in% c("Abundance", "Growth", "PreyNCC", "PreyAK") ~ "Positive",
      rhs %in% c("PredNCC", "PredAK") ~ "Negative",
      TRUE                            ~ "Neutral"
    ),
    
    Type = if_else(str_detect(Lisaname, "DFA"), "DFA", "Straggler"),
    
    OK = case_when(
      Hyp == "Positive" & Pct_Positive > 50 ~ 1,
      Hyp == "Negative" & Pct_Negative > 50 ~ 1,
      TRUE                                  ~ 0
    )
  ) %>%
  
  # Fetch neighbors programmatically for both directions
  rowwise() %>%
  mutate(
    Pos_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "pos"),
    Neg_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "neg")
  ) %>%
  ungroup() %>%
  
  # Add Your Custom Explanations
  mutate(
    Explanation = case_when(
        # Core Salmon Diagnostics (Successes)
      rhs == "Abundance" & OK == 1 ~ "Better NCC survival carries through AK stage",
      rhs == "Growth"    & OK == 1 ~ "Bigger is better",
        
        # Core Salmon Diagnostics (Discrepancies)
      rhs == "Abundance" & OK == 0 ~ "Discrepancy: Survival signal decoupled",
      rhs == "Growth"    & OK == 0 ~ "Discrepancy: Growth signal decoupled",

      # PreyNCC explanations
            #Confirmations
            Lisaname == "X05_DFA_abundSardine" ~ "Competitor (tracking herring / fewer sardines = more food for salmon)",
            Lisaname == "X01_DFA_sumPreyOfPrey_planktonJun" ~ "Direct trophic support (more food = more growth)",
      
      # PredNCC explanations
            #Discrepancies
            Lisaname == "X08_commonMurre_JSOES" ~ "Ecosystem proxy (Murre ~ Habitat Compression Index)",
            Lisaname == "X08_Loons_8_WS" ~ "Coincidence (Loon population co-trending with CPUE)",
      
            #Confirmations
            Lisaname == "X09_DFA_HakeAge5Plus" ~ "Direct predator",
      
      # PredAK explanations
        #Discrepancies
          Lisaname == "X10_Harbor_seal_CR_2yrLead" ~ "Spurious (Noisy / high missing data; trend-rider on short series)",
          Lisaname == "X15_PacificCodBiomass_predAK" ~ "Ecosystem proxy (Pcod had negative response to heat wave, favored other predators)",
          Lisaname == "X15_halibutBiomassAge8plus_2yrLead_predAK" ~ "Ecosystem proxy (Halibut declines in warm water favor competitors)",
    
        #Confirmations
          Lisaname == "X10_Harbour_s_2yrLead_WS" ~ "Direct predator",
          Lisaname == "X10_DFA_ssl.est.wholerange_2yrLead" ~ "Direct predator (deepwater top-down consumer)",
          Lisaname == "X15_DFA_sleeperSharkBSAI_predAK" ~ "Direct predator (deepwater top-down consumer)",
          Lisaname == "X15_ArrowtoothFlounderBiomass_predAK" ~ "Direct predator (deepwater top-down consumer)",
      
      # PreyAK explanations
          #Discrepancies
          # Lisaname == "X12_copepodBiomass_WGoA" ~ "Ecosystem proxy (correlated with predators like Harbor Seals & Northern Fur Seals)",
          # Lisaname == "X12_copepodCom_WGoA" ~ "Ecosystem proxy (correlated with predators like Harbor Seals & Northern Fur Seals)",
          # Lisaname == "X13_hexagram_EAI" ~ "Ecosystem proxy (correlated with predators like Harbor Seals & Northern Fur Seals)",
          # Lisaname == "X12_DFA_biomassEuphShelfSum" ~ "Ecosystem proxy (associated with Salmon Shark abundance and Sea Nettle decline)",
          # Lisaname == "X13_DFA_WGOA_DFA_midTrophic" ~ "Ecosystem proxy (highly correlated with Sablefish recruitment and Shark catch)",
          # Lisaname == "X13_ammod_WAI" ~ "Ecosystem proxy (highly correlated with early Growth, fills in for growth path)",
          
          #Confirmations
          Lisaname == "X12_copepodCom_EGoA" ~ "Direct trophic support (cold-water lipid-rich copepod signal)",
          Lisaname == "X14_pinkSalmon" ~ "Competitor (more pink salmon = less food for Chinook salmon)",
      
      # Fallback default
      TRUE ~ "Investigate neighborhood relationships further"
    )
  ) %>%
  
  # Streamline, make Node an ordered factor, and add physical sorting column
  select(
    Model   = model_id, 
    Node    = rhs,
    Node2,
    Indicator = Lisaname, 
    Type, 
    Count   = Times_In_Top100, 
    Pct_Pos = Pct_Positive, 
    Pct_Neg = Pct_Negative, 
    Hyp, 
    OK, 
    Pos_Neigh, 
    Neg_Neigh, 
    Explanation
  ) %>%
  # 1. Set up the numeric sorting column and set factor level order
  mutate(
    order = case_when(
      Node == "Abundance" ~ 0,
      Node == "Growth"    ~ 0.5,
      Node == "PreyNCC"   ~ 1,
      Node == "PredNCC"   ~ 2,
      Node == "PreyAK"    ~ 3,
      Node == "PredAK"    ~ 4,
      TRUE                ~ 99
    ),
    Node = factor(Node, levels = c("Abundance", "Growth", "PreyNCC", "PredNCC", "PreyAK", "PredAK"))
  ) %>%
  # 2. Arrange exactly as requested: Trophic Order -> Success (Confirmation) -> Alphabetical Indicator
  arrange(order, desc(OK), Indicator)

# --- 3. Save & View the Summary Table ---
print(as_tibble(ecological_synthesis), n = Inf)
write_csv(ecological_synthesis, "LisaXP/outputs_4/Ecological_Synthesis_Summary_Table.csv")
