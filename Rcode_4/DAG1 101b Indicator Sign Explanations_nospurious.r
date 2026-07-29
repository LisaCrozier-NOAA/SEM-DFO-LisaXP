library(tidyverse)

cat("Running Clean Process 2: Individual Indicator Diagnostics & Synthesis...\n")

# ==============================================================================
# SECTION 1: LOAD METADATA & NEIGHBORHOOD DATA
# ==============================================================================
# Load MDS Plot Coordinates
plot_data_mds <- read.csv("LisaXP/outputs_4/mds_coordinates__guild.dfas1_data.csv", row.names = NULL)

# Correlation Matrix across all available years for neighbor lookups
full_data <- guild.dfas1 %>% select(-ncol(.), -year)
cor_matrix_full_years <- cor(full_data, use = "pairwise.complete.obs")

# Helper function to extract correlated neighborhood variables
get_neighbor_summary <- function(target_var, mds_data, direction = "pos") {
  neighbors <- tryCatch({
    find_mds_neighbors(target_var, mds_data, n_neighbors = 10)
  }, error = function(e) NULL)
  
  if (is.null(neighbors) || nrow(neighbors) == 0) return("")
  
  if (direction == "pos") {
    filtered <- neighbors %>% filter(Pairwise_Cor_With_Target > 0.5)
  } else {
    filtered <- neighbors %>% filter(Pairwise_Cor_With_Target < -0.2)
  }
  
  if (nrow(filtered) == 0) return("")
  
  filtered %>%
    mutate(Formatted = sprintf("%s (%s)", Lisaname, round(Pairwise_Cor_With_Target, 2))) %>%
    pull(Formatted) %>%
    paste(collapse = ", ")
}


# ==============================================================================
# SECTION 2: MAP ESTIMATES TO CLEAN TOP 100 MODELS (EXCLUDING SPURIOUS)
# ==============================================================================
# Filter estimates to only include the clean Top 100 model keys (top_100_keys_clean)
mapped_estimates_clean <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(rankings_mapped_clean, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0)
  ) %>%
  # Join strictly with top_100_keys_clean generated in Process 1
  inner_join(top_100_keys_clean, by = c("model_id", "modNum"))


# ==============================================================================
# SECTION 3: SUMMARIZE INDIVIDUAL INDICATOR PERFORMANCE
# ==============================================================================
top100_indicator_behavior_clean <- mapped_estimates_clean %>%
  group_by(model_id, rhs, Lisaname) %>%
  summarise(
    Times_Selected_In_Top100 = n(),
    Times_Positive           = sum(is_positive, na.rm = TRUE),
    Pct_Positive             = round((Times_Positive / Times_Selected_In_Top100) * 100, 1),
    Pct_Negative             = 100 - Pct_Positive,
    Times_Significant        = sum(is_significant, na.rm = TRUE),
    Pct_Significant          = round((Times_Significant / Times_Selected_In_Top100) * 100, 1),
    Mean_Estimate            = round(mean(est, na.rm = TRUE), 3),
    .groups                  = "drop"
  ) %>%
  arrange(model_id, rhs, desc(Times_Selected_In_Top100))

# Save individual indicator behavior table
write_csv(top100_indicator_behavior_clean, "LisaXP/outputs_4/Top_100_Indicator_Sign_and_Significance_nospurious.csv")

print(top100_indicator_behavior_clean,n=Inf)
# ==============================================================================
# SECTION 4: BUILD ECOLOGICAL SYNTHESIS TABLE
# ==============================================================================
ecological_synthesis_clean <- top100_indicator_behavior_clean %>%
  filter(Times_Selected_In_Top100 > 10) %>%
  rename(Times_In_Top100 = Times_Selected_In_Top100) %>%
  # Exclude DAG1C for mediators & NCC nodes
  filter(!(rhs %in% c("Abundance", "Growth", "PreyNCC", "PredNCC") & str_detect(model_id, "DAG1C"))) %>%
  mutate(
    # Node2 incorporates competitors cleanly
    Node2 = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Competitor",
      TRUE ~ as.character(rhs)
    ),
    # Hypotheses already account for Competitors vs. Direct Trophic links
    Hyp = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Negative",
      rhs %in% c("Abundance", "Growth", "PreyNCC", "PreyAK")     ~ "Positive",
      rhs %in% c("PredNCC", "PredAK")                           ~ "Negative",
      TRUE                                                      ~ "Neutral"
    ),
    Type = if_else(str_detect(Lisaname, "DFA"), "DFA", "Straggler"),
    OK   = case_when(
      Hyp == "Positive" & Pct_Positive > 50 ~ 1,
      Hyp == "Negative" & Pct_Negative > 50 ~ 1,
      TRUE                                  ~ 0
    )
  ) %>%
  rowwise() %>%
  mutate(
    Pos_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "pos"),
    Neg_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "neg")
  ) %>%
  ungroup() %>%
  mutate(
    Explanation = case_when(
      Lisaname == "X06_DFA_IGF_mu" ~ "Positive & significant (1B, 1C) or negative but NS (1A)",
      rhs == "Abundance" & OK == 1  ~ "Better NCC survival carries through AK stage",
      rhs == "Growth"    & OK == 1  ~ "Bigger is better",
      rhs == "Abundance" & OK == 0  ~ "Discrepancy: Survival signal decoupled",
      rhs == "Growth"    & OK == 0  ~ "Discrepancy: Growth signal decoupled",
      
      # Competitor confirmations
      Lisaname == "X05_DFA_abundSardine"                       ~ "Competitor (tracking herring / fewer sardines = more food for salmon)",
      Lisaname == "X14_pinkSalmon"                             ~ "Competitor (more pink salmon = less food for Chinook salmon)",
      Lisaname == "X13_DFA_WGOA_DFA_midTrophic"                ~ "MidTrophic = shrimp/jellyfish, potential competitors; neg loading of capelin & eulachon",
      
      # Direct Trophic confirmations
      Lisaname == "X01_DFA_sumPreyOfPrey_planktonJun"         ~ "Direct trophic support (more food = more growth)",
      Lisaname == "X09_DFA_HakeAge5Plus"                       ~ "Direct predator",
      Lisaname == "X10_Harbour_s_2yrLead_WS"                   ~ "Direct predator",
      Lisaname == "X10_DFA_ssl.est.wholerange_2yrLead"         ~ "Direct predator",
      Lisaname == "X10_Californian_s_l_2yrLead_WS"             ~ "Direct predator",
      Lisaname == "X12_copepodCom_EGoA"                        ~ "Direct trophic support (cold-water lipid-rich copepod signal)",
      Lisaname == "X15_DFA_sleeperSharkBSAI_predAK"            ~ "Direct predator (deepwater top-down consumer)",
      Lisaname == "X15_ArrowtoothFlounderBiomass_predAK"       ~ "Direct predator (deepwater top-down consumer)",
      
      # Ecosystem proxies
      Lisaname == "X04_marketsquid_GAM"                        ~ "Ecosystem proxy",
      Lisaname == "X12_copepodBiomass_WGoA"                    ~ "many small copepods, sum of small and large copepods S of Seward AK",
      Lisaname == "X12_copepodCom_WGoA"                        ~ "loads neg on lower trophic DFA, so pos on mid-trophic DFA, ratio of large copepods to (small+large)",
      Lisaname == "X12_DFA_biomassEuphShelfSum"                ~ "Ecosystem proxy (krill ~ pos PDO, Evans et al 2023)",
      Lisaname == "X13_ammod_WAI"                              ~ "Sandlance in puffin chick diet, negligible since 2010",
      Lisaname == "X13_gadid_WAI"                              ~ "Gadids in puffin chick diet",
      Lisaname == "X13_hexagram_EAI"                           ~ "Tufted puffin chick diet, similar to seabird DFA",
      Lisaname == "X13_pollockBiomassGoAage3plus_predAK"       ~ "No real reason limited to 3+, 1+ in WGOA_DFA_midTrophic",
      Lisaname == "X13_DFA_WGOA_DFA_seabirds"                  ~ "Birds have improved since ~2008",
      Lisaname == "X15_PacificCodBiomass_predAK"               ~ "Ecosystem proxy (Pcod had negative response to heat wave, favored other predators)",
      Lisaname == "X15_halibutBiomassAge8plus_2yrLead_predAK"  ~ "Ecosystem proxy (Halibut declines in warm water favor competitors)",
      Lisaname == "X15_sablefishBiomass_predAK"                ~ "Delay in responding to MHW, low in 2015",
      
      TRUE                                                     ~ "Investigate neighborhood relationships further"
    )
  ) %>%
  select(
    Model = model_id, Node = rhs, Node2, Indicator = Lisaname, Type, 
    Count = Times_In_Top100, Pct_Pos = Pct_Positive, Pct_Neg = Pct_Negative, 
    Hyp, OK, Pos_Neigh, Neg_Neigh, Explanation
  ) %>%
  mutate(
    order = case_when(
      Node2 == "Abundance"  ~ 0,
      Node2 == "Growth"     ~ 0.5,
      Node2 == "Competitor" ~ 0.75,
      Node2 == "PreyNCC"    ~ 1,
      Node2 == "PredNCC"    ~ 2,
      Node2 == "PreyAK"     ~ 3,
      Node2 == "PredAK"     ~ 4,
      TRUE                  ~ 99
    ),
    Node  = factor(Node,  levels = c("Abundance", "Growth", "PreyNCC", "PredNCC", "PreyAK", "PredAK")),
    Node2 = factor(Node2, levels = c("Abundance", "Growth", "Competitor", "PreyNCC", "PredNCC", "PreyAK", "PredAK"))
  ) %>%
  arrange(order, desc(OK), Indicator)


# ==============================================================================
# SECTION 5: SUMMARIZE ACROSS MODELS (FINAL OUTPUT)
# ==============================================================================
ecological_summary_clean <- ecological_synthesis_clean %>% 
  group_by(Indicator) %>%
  summarise(
    Node2         = first(Node2),
    Hyp           = first(Hyp),
    Type          = first(Type),
    Explanation   = first(Explanation),
    order         = first(order),
    Support       = mean(OK, na.rm = TRUE),
    mean_Positive = mean(Pct_Pos, na.rm = TRUE),
    Count         = sum(Count, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  relocate(Indicator, Node2, Hyp, Type, Support, mean_Positive, Count, Explanation, order) %>%
  arrange(order, Node2, desc(Support), Type)


print(ecological_summary_clean,n=Inf)

#explore
ecological_summary_clean %>% filter(Support > 0)
ecological_summary_clean %>% filter(Support == 0)

# --- Save Final Clean Tables ---
write_csv(ecological_synthesis_clean, "LisaXP/outputs_4/Indicator_Support_Table__nospurious.csv")
write_csv(ecological_summary_clean,   "LisaXP/outputs_4/Indicator_Summary_Table__nospurious.csv")

cat("Clean Process 2 complete! Summary saved to 'Indicator_Summary_Table_Clean.csv'.\n")