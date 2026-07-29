library(tidyverse)

cat("Running Spurious-Indicator Filtered Process 1...\n")

# --- 1. Map generic SEM nodes to actual indicators across all model rankings ---
rankings_mapped_clean <- rankings_all %>%
  select(
    model_id, modNum, AIC,
    PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
    AbundanceIndNames, PreyAKindNames, PredAKindNames
  ) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "actual_indicator"
  ) %>%
  mutate(generic_node = str_remove(generic_node, "(?i)indnames")) %>%
  left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = c("actual_indicator" = "var"))

# --- 2. Identify and Flag Models Containing Spurious Indicators ---
spurious_model_keys <- rankings_mapped_clean %>%
  filter(
    (generic_node == "PredNCC" & Lisaname %in% c("X08_Loons_8_WS", "X08_commonMurre_JSOES")) |
      (generic_node == "PredAK"  & Lisaname == "X10_Harbor_seal_CR_2yrLead")|
      (generic_node == "PredNCC"  & Lisaname == "X10_Harbor_seal_CR")
  ) %>%
  select(model_id, modNum) %>%
  distinct()

cat(sprintf("Identified and removing %d total model runs containing spurious indicators.\n", nrow(spurious_model_keys)))

# --- 3. Recalculate the CLEAN Top 100 Keys (Excluding HCI + Spurious Models) ---
top_100_keys_clean <- rankings_all %>%
  # Exclude Habitat Compression Index
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  # Exclude the 3 spurious models
  anti_join(spurious_model_keys, by = c("model_id", "modNum")) %>%
  # Slice the NEXT best 100 clean models per model run
  group_by(model_id) %>%
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
  select(model_id, modNum) %>%
  ungroup()

# --- 4. Join Estimates with Cleaned Model Keys & Evaluate ---
targeted_estimates_clean <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_mapped_clean, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0)
  ) %>%
  # Flag clean top 100 models
  left_join(
    top_100_keys_clean %>% mutate(is_top_100 = TRUE), 
    by = c("model_id", "modNum")
  ) %>%
  mutate(is_top_100 = replace_na(is_top_100, FALSE))

# --- 5. Per-Model Evaluation & Summary Table Generation ---
summary_table_clean_spurious <- targeted_estimates_clean %>%
  mutate(
    is_competitor = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname),
    
    Model_Prediction = case_when(
      is_competitor                                          ~ "Negative",
      rhs %in% c("Abundance", "Growth", "PreyNCC", "PreyAK") ~ "Positive",
      rhs %in% c("PredNCC", "PredAK")                        ~ "Negative",
      TRUE                                                  ~ "Positive"
    ),
    
    Model_Supports = case_when(
      Model_Prediction == "Positive" & is_positive == 1 ~ TRUE,
      Model_Prediction == "Negative" & is_negative == 1 ~ TRUE,
      TRUE                                              ~ FALSE
    )
  ) %>%
  group_by(model_id, path_label) %>%
  summarise(
    N_All         = n(),
    Sig_All       = round((sum(is_significant, na.rm = TRUE) / N_All) * 100, 1),
    Pos_All       = round((sum(is_positive, na.rm = TRUE) / N_All) * 100, 1),
    Neg_All       = 100 - Pos_All,
    
    # Clean Top 100 Model Stats
    N_100         = sum(is_top_100),
    Sig_100       = round((sum(is_significant[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    Pos_100       = round((sum(is_positive[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    Neg_100       = 100 - Pos_100,
    
    # % Supported among clean top 100 models
    Pct_Supported_100 = round((sum(Model_Supports[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    N_Competitors_100 = sum(is_competitor[is_top_100], na.rm = TRUE),
    
    .groups  = "drop"
  ) %>%
  mutate(
    H = case_when(
      path_label == "Abundance -> SAR"     ~ "H1",
      path_label == "Growth -> Abundance"  ~ "H2",
      path_label == "Growth -> SAR"        ~ "H3",
      path_label == "PreyNCC -> Growth"    ~ "H4",
      path_label == "PredNCC -> Abundance" ~ "H5",
      path_label == "PreyAK -> SAR"        ~ "H3_ak",
      path_label == "PredAK -> SAR"        ~ "H5",
      path_label == "PreyNCC -> Abundance" ~ "H6",
      path_label == "PreyNCC -> SAR"       ~ "H7",
      path_label == "PredNCC -> SAR"       ~ "H8",
      TRUE                                 ~ " "
    ),
    Hypothesis = case_when(
      path_label == "Abundance -> SAR"     ~ "Early marine survival predicts SAR",
      path_label == "Growth -> Abundance"  ~ "Bigger is better",
      path_label == "Growth -> SAR"        ~ "Bigger is better",
      path_label == "PreyNCC -> Growth"    ~ "More prey improves salmon growth",
      path_label == "PredNCC -> Abundance" ~ "NCC predators reduce early marine survival",
      path_label == "PreyAK -> SAR"        ~ "Bigger is better",
      path_label == "PredAK -> SAR"        ~ "AK predators reduce survival",
      path_label == "PreyNCC -> Abundance" ~ "Alternative prey improves salmon survival",
      path_label == "PreyNCC -> SAR"       ~ "PreyNCC have indirect effect on SAR",
      path_label == "PredNCC -> SAR"       ~ "NCC predators reduce survival beyond JSOES cpue",
      TRUE                                 ~ " "
    ),
    Base_Prediction = case_when(
      path_label %in% c("Abundance -> SAR", "Growth -> Abundance", "Growth -> SAR", 
                        "PreyNCC -> Growth", "PreyAK -> SAR", "PreyNCC -> Abundance", "PreyNCC -> SAR") ~ "Positive",
      path_label %in% c("PredNCC -> Abundance", "PredAK -> SAR", "PredNCC -> SAR")                       ~ "Negative",
      TRUE                                                                                               ~ " "
    ),
    Support = case_when(
      Sig_100 >= 50 & Pct_Supported_100 >= 50 ~ "Yes",
      Sig_100 >= 50 & Pct_Supported_100 < 50  ~ "Flipped Sign",
      Sig_100 < 50                            ~ "NS",
      TRUE                                    ~ " "
    )
  )

summary_table_no_spurious <- summary_table_clean_spurious %>%
  rename(Model = model_id, Path = path_label, Prediction = Base_Prediction) %>%
  select(H, Hypothesis, Path, Prediction, Support, Model, N_All, Sig_All, Pos_All, Neg_All, Sig_100, Pos_100, Neg_100, Pct_Supported_100, N_Competitors_100) %>%
  arrange(H, Path, Support, Model)


print(summary_table_no_spurious,n=Inf)

# --- 6. Save Clean Output ---
write_csv(summary_table_no_spurious, "LisaXP/outputs_4/MainRegression_Significance_Summary_NoSpuriousIndicators.csv")

cat("Clean summary table created and saved to 'MainRegression_Significance_Summary_NoSpuriousIndicators.csv'!\n")
