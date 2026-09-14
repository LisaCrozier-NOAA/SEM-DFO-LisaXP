# ==============================================================================
# Script: 02_goal2_table4_indicator_diagnostics.R
# Directory: Rcode_for_paper/03_goal2_hypothesis_testing/
# Purpose: Goal 2 - Compact Table 4 (All Indicators Retained + Tie-Breaker)
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

library(tidyverse)
library(gt)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DEFINE HYPOTHESES-------
# ------------------------------------------------------------------------------
proj_dir <- getwd()
doug_dir <- "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26"

master_out   <- file.path(proj_dir, "Rcode_for_paper", "Routput_for_paper")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path("Rcode_for_paper", "01_data_prep", "00b_crosswalk_utility_fxn.r"))

crosswalk_path <- file.path("metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}

master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

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

cat("Building Compact Table 4 (All Indicators Retained with Tie-Breaker Logic)...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD RANKINGS, MAP LISANAME & PURGE HCI + CR SEALS----------
# ------------------------------------------------------------------------------
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

lookup_dict     <- setNames(master_crosswalk$LisaName, master_crosswalk$dfa_cols)
lookup_dict_raw <- setNames(master_crosswalk$LisaName, master_crosswalk$DFAname)
lookup_dict     <- c(lookup_dict, lookup_dict_raw)

rankings_long <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  select(model_id, modNum, AIC, ends_with("indNames")) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(
    generic_node = str_remove(generic_node, "(?i)indnames"),
    LisaName = case_when(
      var_name %in% names(lookup_dict) ~ lookup_dict[var_name],
      paste0("X", var_name) %in% names(lookup_dict) ~ lookup_dict[paste0("X", var_name)],
      TRUE ~ var_name
    )
  )

cr_seal_keys <- rankings_long %>%
  filter(LisaName %in% c("X10_Harbor_seal_CR_2yrLead", "X11_Harbor_seal_CR", "11.PredMammalSmolt_0_smoothed")) %>%
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
# STEP 2: LOAD PARAMETER ESTIMATES & EVALUATE MATCH STATUS----------
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
    Series_Type    = if_else(grepl("DFA", LisaName) | grepl("DFA", var_name), "Multi-Taxa DFA", "Single-Series Straggler"),
    is_competitor  = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    
    Expected_Sign = case_when(
      LisaName %in% c("05.ForageFishNCC_DFA1", "X05_DFA_abundSardine") ~ "Negative (-)",
      is_competitor                                                   ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK")                                 ~ "Negative (-)",
      TRUE                                                            ~ "Positive (+)"
    ),
    
    Match_Status = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 & is_significant == 1 ~ "Supported",
      Expected_Sign == "Positive (+)" & is_negative == 1 & is_significant == 1 ~ "Flipped",
      Expected_Sign == "Negative (-)" & is_negative == 1 & is_significant == 1 ~ "Supported",
      Expected_Sign == "Negative (-)" & is_positive == 1 & is_significant == 1 ~ "Flipped",
      TRUE                                                                    ~ "NS"
    )
  ) %>%
  inner_join(top_100_keys_clean, by = c("model_id", "modNum"))

head(indicator_analysis)
# Save indicator_analysis-------
write_csv(indicator_analysis, file.path(data_out_dir, "indicator_analysis.csv"))

# ------------------------------------------------------------------------------
# STEP 3: AGGREGATE ALL INDICATORS (APPLY TIE-BREAKER ONLY TO MIXED RESULTS)---------
# ------------------------------------------------------------------------------
table4_compact <- indicator_analysis %>%
  mutate(SEM_Node = rhs) %>%
  
  # Step 3A: Evaluate each individual indicator across all model runs
  group_by(H_code, path_label, SEM_Node, Series_Type, LisaName) %>%
  summarise(
    N_Runs        = n(),
    N_Supported   = sum(Match_Status == "Supported", na.rm = TRUE),
    N_Flipped     = sum(Match_Status == "Flipped", na.rm = TRUE),
    N_NS          = sum(Match_Status == "NS", na.rm = TRUE),
    Ind_Mean_Beta = mean(est, na.rm = TRUE),
    
    # Categorization Rules:
    # 1. Direct assignment if strictly Supported or strictly Flipped
    # 2. Tie-breaker rule if present in both: assign to majority directional count
    Final_Category = case_when(
      N_Supported > 0 & N_Flipped == 0 ~ "Supported",
      N_Flipped > 0 & N_Supported == 0 ~ "Flipped",
      N_Supported > 0 & N_Flipped > 0 & N_Supported >= N_Flipped ~ "Supported",
      N_Supported > 0 & N_Flipped > 0 & N_Flipped > N_Supported   ~ "Flipped",
      TRUE                             ~ "NS"
    ),
    .groups = "drop"
  ) %>%
  
  # Step 3B: Aggregate into Table 4 categories
  group_by(H_code, path_label, SEM_Node, Series_Type) %>%
  summarise(
    N_Indicators   = n_distinct(LisaName),
    N_Evaluations  = sum(N_Runs),
    Pct_Sig        = round((sum(Final_Category != "NS") / n()) * 100, 1),
    Mean_Beta      = round(mean(Ind_Mean_Beta, na.rm = TRUE), 2),
    
    Supported_Vars = paste(sort(unique(LisaName[Final_Category == "Supported"])), collapse = ", "),
    Flipped_Vars   = paste(sort(unique(LisaName[Final_Category == "Flipped"])), collapse = ", "),
    
    .groups        = "drop"
  ) %>%
  mutate(
    Supported_Vars = if_else(Supported_Vars == "", "None", Supported_Vars),
    Flipped_Vars   = if_else(Flipped_Vars == "", "None", Flipped_Vars)
  ) %>%
  arrange(H_code, SEM_Node, Series_Type)

# Save Summary CSV-------
write_csv(table4_compact, file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.csv"))

# ------------------------------------------------------------------------------
# STEP 4: RENDER GT PUBLICATION TABLE ---------
# ------------------------------------------------------------------------------
table4_gt <- table4_compact %>%
  select(H_code, SEM_Node, Series_Type, N_Indicators, Pct_Sig, Mean_Beta, Supported_Vars, Flipped_Vars) %>%
  gt(groupname_col = "SEM_Node") %>%
  tab_header(
    title = md("**Table 4. Diagnostic Summary of Indicator Drivers Across SEM Nodes**"),
    subtitle = "Complete evaluation of Multi-Taxa DFAs vs. Single-Series Stragglers with tie-breaker classification"
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
  cols_width(
    H_code ~ px(55),
    Series_Type ~ px(100),
    N_Indicators ~ px(40),
    Pct_Sig ~ px(55),
    Mean_Beta ~ px(65),
    Supported_Vars ~ px(400),
    Flipped_Vars ~ px(450)
  ) %>%
  tab_options(
    table.font.size = px(9),
    heading.title.font.size = px(12),
    heading.subtitle.font.size = px(10),
    column_labels.font.size = px(9.5),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.font.size = px(10),
    data_row.padding = px(2),
    row_group.padding = px(3),
    column_labels.padding = px(4),
    table.width = pct(100)
  ) %>%
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

# Save Outputs---------
gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.html"))
gtsave(
  data     = table4_gt,
  filename = file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.png"),
  vwidth   = 1200,
  vheight  = 550
)

message("Table 4 complete with all indicators retained! Saved to: ", tbl_out_dir)