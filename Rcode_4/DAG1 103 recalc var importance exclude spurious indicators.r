library(tidyverse)
library(ggh4x) # Required for facet_nested

cat("Running Full Ensemble Variable Importance Calculation & Barchart Plotting (nospurious)...\n")

# ==============================================================================
# STEP 0: DEFINE PATHS & TARGETS------------------------------------------------
# ==============================================================================
path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")

target_paths <- tribble(
  ~lhs_var,      ~rhs_var,   ~path_label,
  "logSAR",      "Abundance", "Abundance -> SAR",
  "Abundance",   "Growth",    "Growth -> Abundance",
  "logSAR",      "Growth",    "Growth -> SAR",
  "logSAR",      "PreyAK",    "PreyAK -> SAR",
  "Growth",      "PreyNCC",   "PreyNCC -> Growth",
  "Abundance",   "PredNCC",   "PredNCC -> Abundance",
  "logSAR",      "PredAK",    "PredAK -> SAR",
  "Abundance",   "PreyNCC",   "PreyNCC -> Abundance",
  "logSAR",      "PreyNCC",   "PreyNCC -> SAR",
  "logSAR",      "PredNCC",   "PredNCC -> SAR"
) %>% 
  distinct(lhs_var, rhs_var, .keep_all = TRUE)


# ==============================================================================
# STEP 0B: READ COMPLETENESS FILE FOR DATA AVAILABILITY STATUS-----------------
# ==============================================================================
cat("Loading indicator completeness lookup...\n")
completeness_file <- paste0(path, model_names[1], "/completeness.csv")

if (file.exists(completeness_file)) {
  completeness_raw <- read_csv(completeness_file, show_col_types = FALSE)
  
  # Wide format: Columns correspond to indicator 'var' names
  # Variable is complete ONLY if it contains 0 NA values
  completeness_df <- completeness_raw %>%
    select(-1) %>% # Drop date/year column
    summarise(across(everything(), ~ !any(is.na(.)))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "var",
      values_to = "complete"
    ) %>%
    left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = "var") %>%
    filter(!is.na(Lisaname))
} else {
  warning("completeness.csv not found! Check file path.")
}


# ==============================================================================
# STEP 1: LOAD, RESHAPE & MAP VAR -> LISANAME ACROSS ALL MODELS----------------
# ==============================================================================
cat("Loading all model rankings...\n")
rankings_all <- map_dfr(model_names, function(name) {
  file_path <- paste0(path, name, "/SEMresultsByClus.csv")
  if (!file.exists(file_path)) {
    warning(paste("File missing:", file_path))
    return(tibble())
  }
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(model_id = name)
})

# Unnest all indicator columns from SEMresultsByClus.csv into a single long frame
rankings_long <- rankings_all %>%
  # Exclude Habitat Compression Index
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  select(
    model_id, modNum, AIC,
    PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
    AbundanceIndNames, PreyAKindNames, PredAKindNames
  ) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(generic_node = str_remove(generic_node, "(?i)indnames")) %>%
  left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = c("var_name" = "var"))


# ==============================================================================
# STEP 2: IDENTIFY & EXCLUDE SPURIOUS MODELS GLOBALLY---------------------------
# ==============================================================================
# Identify EVERY model_id + modNum combination that ever contained any spurious variable
spurious_model_keys <- rankings_long %>%
  filter(Lisaname %in% c("X08_Loons_8_WS", 
                         "X08_commonMurre_JSOES", 
                         "X10_Harbor_seal_CR_2yrLead", 
                         "X10_Harbor_seal_CR")) %>%
  select(model_id, modNum) %>%
  distinct()

cat(sprintf("Globally purged %d model runs containing spurious indicators.\n", nrow(spurious_model_keys)))

# Keep ONLY models that contain zero spurious indicators anywhere
rankings_long_nospurious <- rankings_long %>%
  anti_join(spurious_model_keys, by = c("model_id", "modNum"))


# ==============================================================================
# STEP 3: CALCULATE AKAIKE WEIGHTS ACROSS THE ENTIRE ENSEMBLE------------------
# ==============================================================================
all_models_with_weights <- rankings_long_nospurious %>%
  select(model_id, modNum, AIC) %>%
  distinct() %>%
  group_by(model_id) %>%
  mutate(
    delta_AIC     = AIC - min(AIC, na.rm = TRUE),
    rel_lik       = exp(-0.5 * delta_AIC),
    akaike_weight = rel_lik / sum(rel_lik, na.rm = TRUE)
  ) %>%
  ungroup()


# ==============================================================================
# STEP 4: SUM VARIABLE IMPORTANCE SCORES BY LISANAME---------------------------
# ==============================================================================
variable_importance_scores <- all_models_with_weights %>%
  inner_join(rankings_long_nospurious, by = c("model_id", "modNum", "AIC")) %>%
  group_by(model_id, generic_node, Lisaname) %>%
  summarise(
    Total_Selections_In_Ensemble = n(),
    Akaike_Importance_Weight     = round(sum(akaike_weight, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(
    Node2 = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Competitor",
      TRUE ~ as.character(generic_node)
    ),
    order = case_when(
      Node2 == "Abundance"  ~ 0,
      Node2 == "Growth"     ~ 0.5,
      Node2 == "Competitor" ~ 0.75,
      Node2 == "PreyNCC"    ~ 1,
      Node2 == "PredNCC"    ~ 2,
      Node2 == "PreyAK"     ~ 3,
      Node2 == "PredAK"     ~ 4,
      TRUE                  ~ 99
    )
  ) %>%
  arrange(model_id, order, desc(Akaike_Importance_Weight))


# ==============================================================================
# STEP 5: OVERALL SUMMARY ACROSS ALL MODEL RUNS---------------------------------
# ==============================================================================
variable_importance_overall <- variable_importance_scores %>%
  group_by(Lisaname, Node2) %>%
  summarise(
    Total_Selections    = sum(Total_Selections_In_Ensemble),
    Mean_Akaike_Weight  = round(mean(Akaike_Importance_Weight), 4),
    Total_Akaike_Weight = round(sum(Akaike_Importance_Weight), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Akaike_Weight))

# Save Variable Importance Products-----
print(variable_importance_scores, n = Inf)
print(variable_importance_overall, n = Inf)

write_csv(variable_importance_scores, "LisaXP/outputs_4/Variable_Importance_By_ModelRun_nospurious.csv")
write_csv(variable_importance_overall, "LisaXP/outputs_4/Variable_Importance_Overall_nospurious.csv")

cat("Full ensemble variable importance calculation complete!\n")


# ==============================================================================
# ==============================================================================
# MAJOR BREAK: FORMING AND PLOTTING TOPVAR IMPORTANCE BARCHART------------------
# ==============================================================================
# ==============================================================================

cat("Building importance_topvar (dAIC <= 3, nospurious) and rendering barchart...\n")

# ==============================================================================
# STEP 6: CALCULATE IMPORTANCE FOR TOP VARS (dAIC <= 3, EXCLUDING SPURIOUS)-----
# ==============================================================================
topvar_models_with_weights <- rankings_long_nospurious %>%
  select(model_id, modNum, AIC) %>%
  distinct() %>%
  group_by(model_id) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) %>%
  filter(delta_AIC <= 3) %>%  # <--- RESTRICTED TO dAIC <= 3
  mutate(
    rel_lik       = exp(-0.5 * delta_AIC),
    akaike_weight = rel_lik / sum(rel_lik, na.rm = TRUE)
  ) %>%
  ungroup()

importance_scores_topvar <- topvar_models_with_weights %>%
  inner_join(rankings_long_nospurious, by = c("model_id", "modNum", "AIC")) %>%
  group_by(model_id, generic_node, Lisaname) %>%
  summarise(
    importance = round(sum(akaike_weight, na.rm = TRUE), 4),
    .groups = "drop"
  )


# ==============================================================================
# STEP 7: ATTACH ECOLOGICAL METADATA & BUILD WIDE FORMAT (`importance_topvar`)--
# ==============================================================================
importance_topvar <- importance_scores_topvar %>%
  mutate(
    Region = if_else(grepl("^X0[1-9]", Lisaname), "NCC", "AK"),
    
    Node2 = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Competitor",
      TRUE ~ as.character(generic_node)
    ),
    Trophic = case_when(
      Node2 == "Abundance"  ~ "Abundance",
      Node2 == "Growth"     ~ "Growth",
      Node2 == "Competitor" ~ "Competitor",
      Node2 == "PreyNCC"    ~ "Prey",
      Node2 == "PredNCC"    ~ "Predator",
      Node2 == "PreyAK"     ~ "Prey",
      Node2 == "PredAK"     ~ "Predator",
      TRUE                  ~ "Other"
    ),
    order = case_when(
      Node2 == "Abundance"  ~ 0,
      Node2 == "Growth"     ~ 0.5,
      Node2 == "Competitor" ~ 0.75,
      Node2 == "PreyNCC"    ~ 1,
      Node2 == "PredNCC"    ~ 2,
      Node2 == "PreyAK"     ~ 3,
      Node2 == "PredAK"     ~ 4,
      TRUE                  ~ 99
    )
  ) %>%
  pivot_wider(
    names_from = model_id,
    values_from = importance
  )


# ==============================================================================
# STEP 8: RESHAPE & ACCURATE MISSING DATA LABELING------------------------------
# ==============================================================================
df_long <- importance_topvar %>%
  mutate(
    Trophic = fct_reorder(Trophic, order),
    Region  = fct_reorder(Region, order),
    plot_names = fct_rev(fct_reorder(Lisaname, order))
  ) %>%
  pivot_longer(
    cols = any_of(c("DAG1A_long", "DAG1B_long", "DAG1A_short", "DAG1B_short", "DAG1C_long", "DAG1C_short")),
    names_to = c("variable", "type"),
    names_sep = "_",
    values_to = "importance"
  )

# Use completeness_df to flag TRUE missing data in long models
missing_labels <- df_long %>%
  filter(is.na(importance)) %>%
  left_join(completeness_df, by = c("plot_names" = "Lisaname")) %>%
  # An indicator genuinely has missing data ONLY if it's a 'long' run AND complete == FALSE
  filter(type == "long" & complete == FALSE) %>%
  group_by(plot_names, type, Trophic, Region) %>%
  summarize(label = "Missing data", .groups = 'drop')

# Render Plot
importance_topvar_barchart <- ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400", "DAG1C" = "#A3C1AD")) +
  geom_text(
    data = missing_labels, 
    aes(x = plot_names, y = 0, label = label), 
    inherit.aes = FALSE,
    hjust = -0.1, 
    color = "grey40", 
    fontface = "italic",
    size = 3
  ) +
  coord_flip() + 
  facet_nested(
    Region + Trophic ~ type, 
    scales = "free_y", 
    space = "free_y",
    nest_line = element_line(color = "black"),
    resect = unit(2, "mm")
  ) + 
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  labs(
    x = "Indicator", 
    y = "Importance (dAIC <= 3)",
    fill = NULL 
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    strip.text.y = element_text(angle = 0),
    panel.spacing.x = unit(2, "lines"), 
    panel.spacing.y = unit(1.5, "lines"),
    strip.background = element_blank()
  )

# Save Barchart Product-----
print(importance_topvar_barchart)

ggsave(
  "LisaXP/outputs_4/DAG1A_DAG1B_DAG1C_importance_topvar_barchart_nospurious.png",
  importance_topvar_barchart,
  width = 11, 
  height = 9, 
  dpi = 300
)

cat("Barchart plotting complete and saved to LisaXP/outputs_4/DAG1A_DAG1B_DAG1C_importance_topvar_barchart_nospurious.png!\n")