# ==============================================================================
# Script: 01_goal1_variable_importance_and_plot.R
# Directory: Rcode_for_paper/02_goal1_importance_analysis/
# Goal 1: Variable Importance Calculation & Faceted Barchart (DAG 1A & 1B)
# Output: Rcode_for_paper/Routput_for_paper/ (data, tables, and figures)
# ==============================================================================

library(tidyverse)
library(ggh4x) # Required for facet_nested
library(forcats)

# ------------------------------------------------------------------------------
# STEP 0: DEFINE PATHS & LOAD UTILITIES-------
# ------------------------------------------------------------------------------
proj_dir     <- getwd()
doug_dir     <- "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26"

input_dir    <- file.path(doug_dir,  "DFA")
#raw_path     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"

# Output Directories
master_out   <- file.path(proj_dir, "Rcode_for_paper", "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")
fig_out_dir  <- file.path(master_out, "figures")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# Source Universal Crosswalk Utility & Master Metadata
source(file.path("Rcode_for_paper", "01_data_prep", "00b_crosswalk_utility_fxn.r"))

crosswalk_path <- file.path("metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}
master_crosswalk <- read.csv(crosswalk_path)

# Models to Evaluate (Baseline linear models: DAG 1A & 1B)
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

cat("Running Goal 1 Ensemble Variable Importance Analysis & Plotting...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD COMPLETENESS LOOKUP --------
# ------------------------------------------------------------------------------
completeness_file <- file.path(doug_dir, model_names[1], "completeness.csv")

if (file.exists(completeness_file)) {
  completeness_raw <- read_csv(completeness_file, show_col_types = FALSE)
  
  completeness_df <- completeness_raw %>%
    select(-1) %>% # Drop date/year column
    summarise(across(everything(), ~ !any(is.na(.)))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "var",
      values_to = "complete"
    ) %>%
    left_join(master_crosswalk %>% select(DFAname, LisaName), by = c("var" = "DFAname")) %>%
    filter(!is.na(LisaName))
} else {
  warning("completeness.csv not found at: ", completeness_file)
}

# ------------------------------------------------------------------------------
# STEP 2: LOAD & MAP MODEL RANKINGS (SEMresultsByClus.csv) ------
# ------------------------------------------------------------------------------
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) {
    warning("File missing: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>%
    mutate(model_id = name)
})

# Unnest indicator columns into long format
rankings_long <- rankings_all %>%
  # Exclude Habitat Compression Index directly
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
  left_join(master_crosswalk %>% select(DFAname, LisaName), by = c("var_name" = "DFAname")) %>%
  # Fallback for any non-DFA stragglers
  mutate(LisaName = if_else(is.na(LisaName), var_name, LisaName))

# ------------------------------------------------------------------------------
# STEP 3: IDENTIFY & EXCLUDE SPECIFIC SPURIOUS MODELS --------
# ------------------------------------------------------------------------------
# Purge model combinations containing Harbor Seal CR variants
spurious_model_keys <- rankings_long %>%
  filter(LisaName %in% c("X10_Harbor_seal_CR_2yrLead", "X11_Harbor_seal_CR")) %>%
  select(model_id, modNum) %>%
  distinct()

cat(sprintf("Globally purged %d model runs containing Harbor Seal CR indicators.\n", nrow(spurious_model_keys)))

rankings_long_nospurious <- rankings_long %>%
  anti_join(spurious_model_keys, by = c("model_id", "modNum"))

# ------------------------------------------------------------------------------
# STEP 4: CALCULATE AKAIKE WEIGHTS & VARIABLE IMPORTANCE (non-spurious) -------
# ------------------------------------------------------------------------------
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

variable_importance_scores <- all_models_with_weights %>%
  inner_join(rankings_long_nospurious, by = c("model_id", "modNum", "AIC")) %>%
  group_by(model_id, generic_node, LisaName) %>%
  summarise(
    Total_Selections_In_Ensemble = n(),
    Akaike_Importance_Weight     = round(sum(akaike_weight, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(
    Node2 = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName) ~ "Competitor",
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

variable_importance_overall <- variable_importance_scores %>%
  group_by(LisaName, Node2) %>%
  summarise(
    Total_Selections    = sum(Total_Selections_In_Ensemble),
    Mean_Akaike_Weight  = round(mean(Akaike_Importance_Weight), 4),
    Total_Akaike_Weight = round(sum(Akaike_Importance_Weight), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Akaike_Weight))

# Save overall importance outputs----------
write_csv(variable_importance_scores, file.path(data_out_dir, "Variable_Importance_By_ModelRun.csv"))
write_csv(variable_importance_overall, file.path(data_out_dir, "Variable_Importance_Overall.csv"))

# ------------------------------------------------------------------------------
# STEP 5: CALCULATE IMPORTANCE FOR TOP VARS (dAIC <= 3) -----
# ------------------------------------------------------------------------------
topvar_models_with_weights <- rankings_long_nospurious %>%
  select(model_id, modNum, AIC) %>%
  distinct() %>%
  group_by(model_id) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) %>%
  filter(delta_AIC <= 3) %>%
  mutate(
    rel_lik       = exp(-0.5 * delta_AIC),
    akaike_weight = rel_lik / sum(rel_lik, na.rm = TRUE)
  ) %>%
  ungroup()

importance_scores_topvar <- topvar_models_with_weights %>%
  inner_join(rankings_long_nospurious, by = c("model_id", "modNum", "AIC")) %>%
  group_by(model_id, generic_node, LisaName) %>%
  summarise(
    importance = round(sum(akaike_weight, na.rm = TRUE), 4),
    .groups = "drop"
  )

importance_topvar <- importance_scores_topvar %>%
  mutate(
    Region = if_else(grepl("^X0[1-9]", LisaName), "NCC", "AK"),
    Node2 = case_when(
      grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName) ~ "Competitor",
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

#save topvar-------
write_csv(importance_topvar, file.path(data_out_dir, "importance_topvar_daic3.csv"))

# ------------------------------------------------------------------------------
# STEP 6: PLOT BARCHART---------
# ------------------------------------------------------------------------------
df_long <- importance_topvar %>%
  filter(order > 0) %>% # Omit primary salmon response anchors
  mutate(
    Trophic    = fct_reorder(Trophic, order),
    Region     = fct_reorder(Region, order),
    plot_names = fct_rev(fct_reorder(LisaName, order))
  ) %>%
  pivot_longer(
    cols      = any_of(c("DAG1A_long", "DAG1B_long", "DAG1A_short", "DAG1B_short")),
    names_to  = c("variable", "type"),
    names_sep = "_",
    values_to = "importance"
  )

# Flag TRUE missing data for long time series runs
missing_labels <- df_long %>%
  filter(is.na(importance)) %>%
  left_join(completeness_df, by = c("plot_names" = "LisaName")) %>%
  filter(type == "long" & complete == FALSE) %>%
  group_by(plot_names, type, Trophic, Region) %>%
  summarize(label = "Missing data", .groups = "drop")

importance_topvar_barchart <- ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(
    values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400"),
    labels = c("DAG 1A (Sequential)", "DAG 1B (Growth Direct)")
  ) +
  geom_text(
    data        = missing_labels, 
    aes(x = plot_names, y = 0, label = label), 
    inherit.aes = FALSE,
    hjust       = -0.1, 
    color       = "grey40", 
    fontface    = "italic",
    size        = 3
  ) +
  coord_flip() + 
  facet_nested(
    Region + Trophic ~ type, 
    scales     = "free_y", 
    space      = "free_y",
    nest_line  = element_line(color = "black"),
    resect     = unit(2, "mm")
  ) + 
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1.10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x    = "Ecosystem Indicator", 
    y    = "AIC Importance Weight (dAIC <= 3)",
    fill = "DAG Model:" 
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    strip.text.y       = element_text(angle = 0, face = "bold"),
    strip.text.x       = element_text(face = "bold"),
    panel.spacing.x    = unit(1.5, "lines"), 
    panel.spacing.y    = unit(1.0, "lines"),
    strip.background   = element_blank(),
    legend.position    = "bottom"
  )

# Export Final Plot Products--------
print(importance_topvar_barchart)
ggsave(
  file.path(fig_out_dir, "DAG1A_DAG1B_importance_topvar_barchart.png"),
  importance_topvar_barchart,
  width = 10, 
  height = 8, 
  dpi = 300
)

ggsave(
  file.path(fig_out_dir, "DAG1A_DAG1B_importance_topvar_barchart.pdf"),
  importance_topvar_barchart,
  width = 10, 
  height = 8
)

cat("Goal 1 Analysis Complete! Figures and tables written to: ", master_out, "\n")
