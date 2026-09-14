# ==============================================================================
# Script: 02_goal2_table4_indicator_diagnostics.R
# Directory: Rcode_for_paper/03_goal2_hypothesis_testing/
# Purpose: Goal 2 - Compact Table 4 with Standardized LisaNames via Crosswalk
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

library(tidyverse)
library(gt)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD MASTER CROSSWALK UTILITY
# ------------------------------------------------------------------------------
proj_dir <- getwd()
doug_dir <- "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26"

master_out   <- file.path(proj_dir, "Rcode_for_paper", "Routput_for_paper")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

# Source Universal Crosswalk Utility
source(file.path("Rcode_for_paper", "01_data_prep", "00b_crosswalk_utility_fxn.r"))

# Load Master Crosswalk directly
crosswalk_path <- file.path("metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}

master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# Active Baseline Models
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

# Target Structural Paths
target_paths <- tribble(
  ~lhs_var,      ~rhs_var,     ~path_label,             ~H_code,
  "Growth",      "PreyNCC",    "PreyNCC -> Growth",     "H1a",
  "Abundance",   "Growth",     "Growth -> CPUE",        "H1b",
  "logSAR",      "Growth",     "Growth -> SAR",         "H1c",
  "logSAR",      "PreyAK",     "PreyAK -> SAR",         "H1e",
  "Abundance",   "PredNCC",    "PredNCC -> CPUE",       "H2a",
  "logSAR",      "PredAK",     "PredAK -> SAR",         "H2c",
  "logSAR",      "Abundance",  "CPUE -> SAR",           "H4a"
)

cat("Building Compact Table 4 with Crosswalk LisaName Mapping...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD RANKINGS, APPLY CROSSWALK & PURGE HCI + CR SEALS
# ------------------------------------------------------------------------------
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

# Build comprehensive lookup vector matching DFAname, dfa_cols, and raw indicator names
lookup_dict <- setNames(master_crosswalk$LisaName, master_crosswalk$dfa_cols)
lookup_dict_raw <- setNames(master_crosswalk$LisaName, master_crosswalk$DFAname)
lookup_dict <- c(lookup_dict, lookup_dict_raw)

# Unnest generic node indicator columns and map to LisaName
rankings_long <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>% # Exclude HCI
  select(model_id, modNum, AIC, ends_with("indNames")) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(
    generic_node = str_remove(generic_node, "(?i)indnames"),
    # Convert 'var_name' or 'X' + 'var_name' to LisaName using crosswalk lookup
    LisaName = case_when(
      var_name %in% names(lookup_dict) ~ lookup_dict[var_name],
      paste0("X", var_name) %in% names(lookup_dict) ~ lookup_dict[paste0("X", var_name)],
      TRUE ~ var_name
    )
  )

# Purge CR Harbor Seal Runs
cr_seal_keys <- rankings_long %>%
  filter(LisaName %in% c("X10_Harbor_seal_CR_2yrLead", "X10_Harbor_seal_CR", "11.PredMammalSmolt_0_smoothed")) %>%
  select(model_id, modNum) %>%
  distinct()

rankings_long_clean <- rankings_long %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum"))

top_100_keys_clean <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
  select(model_id, modNum) %>%
  ungroup()

# ------------------------------------------------------------------------------
# STEP 2: LOAD PARAMETER ESTIMATES & EVALUATE MATCH STATUS
# ------------------------------------------------------------------------------
estimates_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "parameterEstimates.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

indicator_analysis <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_long_clean, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0),
    
    # Classify Series Structure
    Series_Type    = if_else(grepl("DFA", LisaName) | grepl("DFA", var_name), "Multi-Taxa DFA", "Single-Series Straggler"),
    
    # Classify Competitor Species
    is_competitor  = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    
    # Expected Ecological Direction
    Expected_Sign  = case_when(
      is_competitor                   ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK") ~ "Negative (-)",
      TRUE                            ~ "Positive (+)"
    ),
    
    # Match Status Evaluation
    Match_Status = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 & is_significant == 1 ~ "Supported",
      Expected_Sign == "Positive (+)" & is_negative == 1 & is_significant == 1 ~ "Flipped",
      Expected_Sign == "Negative (-)" & is_negative == 1 & is_significant == 1 ~ "Supported",
      Expected_Sign == "Negative (-)" & is_positive == 1 & is_significant == 1 ~ "Flipped",
      TRUE                                                                    ~ "NS"
    )
  ) %>%
  inner_join(top_100_keys_clean, by = c("model_id", "modNum"))

# ------------------------------------------------------------------------------
# STEP 3: AGGREGATE BY NODE & SERIES TYPE USING LISANAME
# ------------------------------------------------------------------------------
table4_compact <- indicator_analysis %>%
  group_by(H_code, path_label, rhs, Series_Type) %>%
  summarise(
    N_Indicators   = n_distinct(LisaName),
    N_Evaluations  = n(),
    Pct_Sig        = round((sum(is_significant) / N_Evaluations) * 100, 1),
    Mean_Beta      = round(mean(est, na.rm = TRUE), 2),
    
    # Concatenate unique LisaNames for Supported and Flipped categories
    Supported_Vars = paste(sort(unique(LisaName[Match_Status == "Supported"])), collapse = ", "),
    Flipped_Vars   = paste(sort(unique(LisaName[Match_Status == "Flipped"])), collapse = ", "),
    
    .groups = "drop"
  ) %>%
  mutate(
    Supported_Vars = if_else(Supported_Vars == "", "None", Supported_Vars),
    Flipped_Vars   = if_else(Flipped_Vars == "", "None", Flipped_Vars)
  ) %>%
  rename(SEM_Node = rhs) %>%
  arrange(H_code, SEM_Node, Series_Type)

# Save Raw CSV Summary
write_csv(table4_compact, file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.csv"))

# # ------------------------------------------------------------------------------
# # STEP 4: RENDER FORMATTED GT PUBLICATION TABLE
# # ------------------------------------------------------------------------------
# table4_gt <- table4_compact %>%
#   select(H_code, SEM_Node, Series_Type, N_Indicators, Pct_Sig, Mean_Beta, Supported_Vars, Flipped_Vars) %>%
#   gt(groupname_col = "SEM_Node") %>%
#   tab_header(
#     title = md("**Table 4. Diagnostic Summary of Indicator Drivers Across SEM Nodes**"),
#     subtitle = "Performance comparison of Multi-Taxa DFAs vs. Single-Series Stragglers using standardized crosswalk labels"
#   ) %>%
#   cols_label(
#     H_code         = md("**Code**"),
#     Series_Type    = md("**Series Structure**"),
#     N_Indicators   = md("**N Indicators**"),
#     Pct_Sig        = md("**% Sig**"),
#     Mean_Beta      = md("**Mean β**"),
#     Supported_Vars = md("**Supported Indicators (LisaName)**"),
#     Flipped_Vars   = md("**Flipped-Sign Indicators (LisaName)**")
#   ) %>%
#   cols_align(align = "center", columns = c(H_code, Series_Type, N_Indicators, Pct_Sig, Mean_Beta)) %>%
#   cols_align(align = "left", columns = c(Supported_Vars, Flipped_Vars)) %>%
#   tab_options(
#     table.font.size = px(10),
#     heading.title.font.size = px(12),
#     column_labels.font.weight = "bold",
#     row_group.font.weight = "bold"
#   )
# 
# # Save HTML and High-Res PNG
# gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.html"))
# gtsave(
#   data     = table4_gt,
#   filename = file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.png"),
#   vwidth   = 1200,
#   vheight  = 700
# )
# 
# message("Table 4 with standardized LisaName labels generated successfully! Saved to: ", tbl_out_dir)

# ------------------------------------------------------------------------------
# STEP 4: RENDER COMPACT GT PUBLICATION TABLE (TIGHT CELL PADDING FOR PPT)
# ------------------------------------------------------------------------------
table4_gt <- table4_compact %>%
  select(H_code, SEM_Node, Series_Type, N_Indicators, Pct_Sig, Mean_Beta, Supported_Vars, Flipped_Vars) %>%
  gt(groupname_col = "SEM_Node") %>%
  tab_header(
    title = md("**Table 4. Diagnostic Summary of Indicator Drivers Across SEM Nodes**"),
    subtitle = "Performance comparison of Multi-Taxa DFAs vs. Single-Series Stragglers using standardized crosswalk labels"
  ) %>%
  cols_label(
    H_code         = md("**Code**"),
    Series_Type    = md("**Series Structure**"),
    N_Indicators   = md("**N**"),
    Pct_Sig        = md("**% Sig**"),
    Mean_Beta      = md("**Mean β**"),
    Supported_Vars = md("**Supported Indicators (LisaName)**"),
    Flipped_Vars   = md("**Flipped-Sign Indicators (LisaName)**")
  ) %>%
  cols_align(align = "center", columns = c(H_code, Series_Type, N_Indicators, Pct_Sig, Mean_Beta)) %>%
  cols_align(align = "left", columns = c(Supported_Vars, Flipped_Vars)) %>%
  # Set specific column widths to force text wrapping on long variable lists
  cols_width(
    H_code ~ px(55),
    Series_Type ~ px(100),
    N_Indicators ~ px(40),
    Pct_Sig ~ px(55),
    Mean_Beta ~ px(65),
    Supported_Vars ~ px(400),
    Flipped_Vars ~ px(450)
  ) %>%
  # Compact Styling Options
  tab_options(
    table.font.size = px(9),                        # Slightly smaller font for tight layouts
    heading.title.font.size = px(12),
    heading.subtitle.font.size = px(10),
    column_labels.font.size = px(9.5),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.font.size = px(10),
    data_row.padding = px(2),                       # Tightens row padding (default is ~8px)
    row_group.padding = px(3),                      # Tightens header group padding
    column_labels.padding = px(4),                  # Tightens column label padding
    table.width = pct(100)
  ) %>%
  # Add CSS styling for tight line heights in multiline text cells
  opt_css(
    css = "
    .gt_row {
      line-height: 1.15 !important;
      vertical-align: middle !important;
    }
    .gt_group_heading {
      background-color: #f4f4f4 !important;
      padding-top: 3px !important;
      padding-bottom: 3px !important;
    }
    "
  )

# Save HTML and High-Res Compressed PNG for PowerPoint
gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.html"))

gtsave(
  data     = table4_gt,
  filename = file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.png"),
  vwidth   = 1200,   # High resolution width
  vheight  = 550     # Reduced vertical height to eliminate bottom empty space
)

message("Tightened Table 4 generated successfully! Saved to: ", tbl_out_dir)
