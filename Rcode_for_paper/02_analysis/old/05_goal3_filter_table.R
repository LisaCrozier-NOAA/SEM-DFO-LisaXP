library(tidyverse)

#Extract SST models------

# 1. Load the Goal 3 Comparison Summaries
summaries <- read.csv("Rcode_for_paper/Routput_for_paper/data/Goal3_Model_Summaries_Comparison.csv")

# 2. Filter exclusively for Thermal Moderation models
thermal_models <- summaries %>%
  filter(model_type == "Thermal Moderation (Pred*SST)") %>%
  arrange(aic)

# Display Top 10 Thermal Models
print(head(thermal_models, 10))

# 3. Load Parameter Estimates for the Top Thermal Model
params <- read.csv("Rcode_for_paper/Routput_for_paper/data/Goal3_All_Parameter_Estimates.csv")

top_thermal_params <- params %>%
  filter(
    model_type == "Thermal Moderation (Pred*SST)",
    pair_id == thermal_models$pair_id[1]
  )

print(top_thermal_params)

# ==============================================================================
# Script: 02_goal3_generate_table5_from_saved_data.R
# Directory: Rcode_for_paper/04_goal3_interaction_analysis/
# Purpose: Generate Table 5 including Thermal Moderation Models
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(stringr)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD goal 3 model results--------
# ------------------------------------------------------------------------------
proj_dir     <- getwd()
master_out   <- file.path(proj_dir, "Rcode_for_paper", "Routput_for_paper")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

summaries_path <- file.path(data_out_dir, "Goal3_Model_Summaries_Comparison.csv")
params_path    <- file.path(data_out_dir, "Goal3_All_Parameter_Estimates.csv")

if (!file.exists(summaries_path) || !file.exists(params_path)) {
  stop("Missing required input CSVs in: ", data_out_dir)
}

master_crosswalk <- read.csv(file.path("metadata", "master_name_crosswalk.csv"), stringsAsFactors = FALSE)

lookup_dict     <- setNames(master_crosswalk$LisaName, tolower(master_crosswalk$dfa_cols))
lookup_dict_raw <- setNames(master_crosswalk$LisaName, tolower(master_crosswalk$DFAname))
lookup_dict     <- c(lookup_dict, lookup_dict_raw)

get_lisa_name <- function(raw_col) {
  raw_lower <- tolower(raw_col)
  if (raw_lower %in% names(lookup_dict)) return(lookup_dict[[raw_lower]])
  return(raw_col)
}

cat("Loading pre-calculated Goal 3 datasets...\n")

summary_results <- read_csv(summaries_path, show_col_types = FALSE)
all_params_df   <- read_csv(params_path, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: FILTERING LOGIC----
# For Main Effects models: requires beta <= 0 --------
# For Interaction models: int <= 0 -------
#retains models where the interaction term provides top-down moderation
# ------------------------------------------------------------------------------
# Extract interaction parameter signs
interaction_slopes <- all_params_df %>%
  filter(term_type %in% c("Thermal Interaction (Pred*SST)", "AltPrey Interaction (Pred*Prey)")) %>%
  group_by(pair_id, model_type) %>%
  summarise(
    has_negative_int = any(est < 0),
    .groups = "drop"
  )

# Combine and Filter
summary_filtered <- summary_results %>%
  filter(global_p>0.05) %>%
  left_join(interaction_slopes, by = c("pair_id", "model_type")) %>%
  filter(
    # Keep baseline models only if main effect is top-down (beta <= 0)
    (model_type == "Baseline DAG 1A (Main Effects)" & beta_ncc_pred <= 0 & beta_ak_pred <= 0) |
      # Keep interaction models if thermal or altprey interaction term is active
      (model_type != "Baseline DAG 1A (Main Effects)")
  )

# ------------------------------------------------------------------------------
# STEP 2: TABLE PREPARED -------
#CONSTRUCT EXPLICIT TREATMENT COLUMNS & MAP LISANAME
# ------------------------------------------------------------------------------
table5_prepared <- summary_filtered %>%
  mutate(
    ncc_pred_clean = sapply(ncc_predator, get_lisa_name),
    ak_pred_clean  = sapply(ak_predator,  get_lisa_name),
    
    ncc_treatment = case_when(
      model_type == "Thermal Moderation (Pred*SST)"  ~ paste0(ncc_pred_clean, " × SST"),
      model_type == "AltPrey Moderation (Pred*Prey)" ~ paste0(ncc_pred_clean, " × AltPrey"),
      TRUE                                           ~ "None (Main Effect Only)"
    ),
    
    ak_treatment = case_when(
      model_type == "Thermal Moderation (Pred*SST)"  ~ paste0(ak_pred_clean, " × SST"),
      model_type == "AltPrey Moderation (Pred*Prey)" ~ paste0(ak_pred_clean, " × AltPrey"),
      TRUE                                           ~ "None (Main Effect Only)"
    )
  ) %>%
  group_by(model_type) %>%
  slice_min(order_by = aic, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(aic) %>%
  mutate(delta_aic_overall = round(aic - min(aic), 2))

# ------------------------------------------------------------------------------
# STEP 3: RENDER GT TABLE-----------
# ------------------------------------------------------------------------------
table5_gt <- table5_prepared %>%
  select(
    model_type, 
    ncc_pred_clean, 
    ncc_treatment, 
    ak_pred_clean, 
    ak_treatment, 
    aic, 
    delta_aic_overall, 
    cfi, 
    cpue_r2, 
    sar_r2
  ) %>%
  gt(groupname_col = "model_type") %>%
  tab_header(
    title = md("**Table 5. Model Evaluation: Baseline DAG 1A vs. Non-Linear Interaction Structures**"),
    subtitle = "Standardized evaluation on complete 1998–2021 observations"
  ) %>%
  cols_label(
    ncc_pred_clean    = md("**NCC Predator**"),
    ncc_treatment     = md("**NCC Moderation Term**"),
    ak_pred_clean     = md("**AK Predator**"),
    ak_treatment      = md("**AK Moderation Term**"),
    aic               = md("**AIC**"),
    delta_aic_overall = md("**ΔAIC**"),
    cfi               = md("**CFI**"),
    cpue_r2           = md("**R² (CPUE)**"),
    sar_r2            = md("**R² (SAR)**")
  ) %>%
  cols_align(align = "center", columns = c(aic, delta_aic_overall, cfi, cpue_r2, sar_r2)) %>%
  cols_align(align = "left", columns = c(ncc_pred_clean, ncc_treatment, ak_pred_clean, ak_treatment)) %>%
  fmt_number(columns = c(aic, delta_aic_overall, cfi, cpue_r2, sar_r2), decimals = 2) %>%
  tab_options(
    table.font.size = px(9),
    heading.title.font.size = px(11.5),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold"
  )

# Save Outputs-----------
gtsave(table5_gt, file.path(tbl_out_dir, "Table5_Model_Comparison_Baseline_vs_Interactions.html"))
gtsave(
  data     = table5_gt,
  filename = file.path(tbl_out_dir, "Table5_Model_Comparison_Baseline_vs_Interactions.png"),
  vwidth   = 1250,
  vheight  = 650
)

message("Table 5 successfully restored with Thermal Moderation models! Saved to: ", tbl_out_dir)
