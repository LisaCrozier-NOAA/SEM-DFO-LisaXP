# ==============================================================================
# UPDATED PROCESS 1: MODEL-SPECIFIC COMPETITOR EVALUATION & ACCURATE AGGREGATION
# ==============================================================================
cat("Running Process 1 with True Per-Model Competitor Evaluation...\n")

# 1. Map generic SEM nodes to the actual indicators used in each model run
rankings_mapped_p1 <- rankings_all %>%
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
  mutate(generic_node = str_remove(generic_node, "(?i)indnames"))

# 2. Join estimates with mapped indicators to identify exact variable per model run
targeted_estimates <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_mapped_p1, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = c("actual_indicator" = "var")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0)
  ) %>%
  left_join(
    top_100_keys %>% mutate(is_top_100 = TRUE), 
    by = c("model_id", "modNum")
  ) %>%
  mutate(is_top_100 = replace_na(is_top_100, FALSE))

# 3. Evaluate EACH model run individually against its specific variable hypothesis
targeted_estimates_evaluated <- targeted_estimates %>%
  mutate(
    is_competitor = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname),
    
    # Specific prediction FOR THIS MODEL RUN ONLY
    Model_Prediction = case_when(
      is_competitor                                      ~ "Negative",
      rhs %in% c("Abundance", "Growth", "PreyNCC", "PreyAK") ~ "Positive",
      rhs %in% c("PredNCC", "PredAK")                        ~ "Negative",
      TRUE                                                  ~ "Positive"
    ),
    
    # Did THIS specific model support its specific indicator prediction?
    Model_Supports = case_when(
      Model_Prediction == "Positive" & is_positive == 1 ~ TRUE,
      Model_Prediction == "Negative" & is_negative == 1 ~ TRUE,
      TRUE                                              ~ FALSE
    )
  )

# 4. Summarize aggregate support accurately without altering baseline hypothesis labels
summary_table <- targeted_estimates_evaluated %>%
  group_by(model_id, path_label) %>%
  summarise(
    N_All         = n(),
    Sig_All       = round((sum(is_significant, na.rm = TRUE) / N_All) * 100, 1),
    Pos_All       = round((sum(is_positive, na.rm = TRUE) / N_All) * 100, 1),
    Neg_All       = 100 - Pos_All,
    
    # Top 100 Model Stats
    N_100         = sum(is_top_100),
    Sig_100       = round((sum(is_significant[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    Pos_100       = round((sum(is_positive[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    Neg_100       = 100 - Pos_100,
    
    # % of Top 100 models that supported their individual variable predictions
    Pct_Supported_100 = round((sum(Model_Supports[is_top_100], na.rm = TRUE) / N_100) * 100, 1),
    
    # Count of how many competitor variables were selected in Top 100
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
    # Baseline expected hypothesis remains pure and uncorrupted
    Base_Prediction = case_when(
      path_label %in% c("Abundance -> SAR", "Growth -> Abundance", "Growth -> SAR", 
                        "PreyNCC -> Growth", "PreyAK -> SAR", "PreyNCC -> Abundance", "PreyNCC -> SAR") ~ "Positive",
      path_label %in% c("PredNCC -> Abundance", "PredAK -> SAR", "PredNCC -> SAR")                       ~ "Negative",
      TRUE                                                                                               ~ " "
    ),
    # Final Support decision based on whether >50% of Top 100 models were significant AND aligned
    Support = case_when(
      Sig_100 >= 50 & Pct_Supported_100 >= 50 ~ "Yes",
      Sig_100 >= 50 & Pct_Supported_100 < 50  ~ "Flipped Sign",
      Sig_100 < 50                            ~ "NS",
      TRUE                                    ~ " "
    )
  )

summary_table_competitors <- summary_table %>%
  rename(Model = model_id, Path = path_label, Prediction = Base_Prediction) %>%
  select(H, Hypothesis, Path, Prediction, Support, Model, N_All, Sig_All, Pos_All, Neg_All, Sig_100, Pos_100, Neg_100, Pct_Supported_100, N_Competitors_100) %>%
  arrange(H, Path, Support, Model)


print(summary_table_competitors,n=Inf)


# Save corrected Process 1 output
write_csv(summary_table_competitors, "LisaXP/outputs_4/MainRegression_Significance_Summary_competitors.csv")
