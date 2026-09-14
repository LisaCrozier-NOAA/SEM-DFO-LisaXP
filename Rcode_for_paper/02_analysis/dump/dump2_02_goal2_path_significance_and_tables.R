# ==============================================================================
# Script: 01_goal2_path_significance_and_tables.R
# Directory: Rcode_for_paper/03_goal2_hypothesis_testing/
# Purpose: Goal 2 - Individual Diagnostics & Manuscript-Ready Table 3
# Output: Rcode_for_paper/Routput_for_paper/ (data and tables)
# ==============================================================================

library(tidyverse)
library(gt)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD MASTER CROSSWALK ONLY
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

# Single Source of Truth
crosswalk_path <- file.path("metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}

master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# Active Baseline Models (DAG 1A & 1B across long and short windows)
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

# Target Hypothesized Structural Paths
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

cat("Running Goal 2 Pipeline (Diagnostics & Table 3)...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD MODEL RANKINGS & PURGE HCI + CR SEALS
# ------------------------------------------------------------------------------
cat("Loading model rankings...\n")
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) {
    warning("File missing: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

# Unnest generic node indicator columns and map exclusively via master_crosswalk
rankings_long <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>% # Exclude HCI
  select(model_id, modNum, AIC, ends_with("indNames")) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(generic_node = str_remove(generic_node, "(?i)indnames")) %>%
  left_join(master_crosswalk %>% select(dfa_cols, LisaName), by = c("var_name" = "dfa_cols")) %>%
  mutate(LisaName = if_else(is.na(LisaName), var_name, LisaName))

# Identify & globally purge model runs containing CR Harbor Seals
cr_seal_keys <- rankings_long %>%
  filter(LisaName %in% c("X10_Harbor_seal_CR_2yrLead", "X10_Harbor_seal_CR", "11.PredMammalSmolt_0_smoothed")) %>%
  select(model_id, modNum) %>%
  distinct()

cat(sprintf("Purged %d model runs containing CR Harbor Seal indicators.\n", nrow(cr_seal_keys)))

rankings_long_clean <- rankings_long %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum"))

# Extract clean Top 100 models per model run
top_100_keys_clean <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
  select(model_id, modNum) %>%
  ungroup()

# ------------------------------------------------------------------------------
# STEP 2: LOAD PARAMETER ESTIMATES & EVALUATE COMPETITOR SIGNS
# ------------------------------------------------------------------------------
cat("Loading parameter estimates...\n")
estimates_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "parameterEstimates.csv")
  if (!file.exists(f_path)) {
    warning("File missing: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

targeted_estimates <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_long_clean, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0),
    is_competitor  = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    
    # Model-specific prediction based on species ecological role
    Expected_Sign = case_when(
      is_competitor                   ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK") ~ "Negative (-)",
      TRUE                            ~ "Positive (+)"
    ),
    
    # Check if estimated sign matches expected sign
    Model_Supports = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 ~ TRUE,
      Expected_Sign == "Negative (-)" & is_negative == 1 ~ TRUE,
      TRUE                                               ~ FALSE
    )
  ) %>%
  inner_join(top_100_keys_clean, by = c("model_id", "modNum"))

# ------------------------------------------------------------------------------
# OUTPUT 1: INDIVIDUAL INDICATOR DIAGNOSTICS CSV
# ------------------------------------------------------------------------------
indicator_diagnostic_summary <- targeted_estimates %>%
  group_by(H_code, path_label, model_id, rhs, LisaName, Expected_Sign) %>%
  summarise(
    N_Top_Runs      = n(),
    Pct_Significant = round((sum(is_significant) / N_Top_Runs) * 100, 1),
    Pct_Positive    = round((sum(is_positive) / N_Top_Runs) * 100, 1),
    Pct_Negative    = round((sum(is_negative) / N_Top_Runs) * 100, 1),
    Pct_Supported   = round((sum(Model_Supports) / N_Top_Runs) * 100, 1),
    Mean_Beta       = round(mean(est, na.rm = TRUE), 3),
    .groups         = "drop"
  ) %>%
  arrange(H_code, path_label, model_id, desc(N_Top_Runs))

write_csv(indicator_diagnostic_summary, file.path(data_out_dir, "Goal2_Individual_Indicator_Diagnostics.csv"))
cat("Diagnostic CSV saved to: Goal2_Individual_Indicator_Diagnostics.csv\n")

# ------------------------------------------------------------------------------
# OUTPUT 2: MANUSCRIPT-READY TABLE 3 (WIDE FORMAT)
# ------------------------------------------------------------------------------
path_summary_by_model <- targeted_estimates %>%
  group_by(model_id, H_code, path_label) %>%
  summarise(
    N_Top100      = n(),
    Pct_Sig       = round((sum(is_significant) / N_Top100) * 100, 1),
    Pct_Supported = round((sum(Model_Supports) / N_Top100) * 100, 1),
    .groups       = "drop"
  ) %>%
  mutate(
    Baseline_Expectation = case_when(
      H_code %in% c("H2a", "H2c") ~ "Negative (-)",
      TRUE                        ~ "Positive (+)"
    ),
    Verdict = case_when(
      Pct_Sig >= 50 & Pct_Supported >= 50 ~ "Supported",
      Pct_Sig >= 50 & Pct_Supported < 50  ~ "Flipped Sign",
      TRUE                                ~ "NS"
    ),
    Cell_Display = paste0(Verdict, "\n(", Pct_Sig, "% Sig)")
  )

# Reshape into Wide Format (Hypothesis Rows x Model Columns)
table3_wide <- path_summary_by_model %>%
  select(H_code, path_label, Baseline_Expectation, model_id, Cell_Display) %>%
  pivot_wider(
    names_from  = model_id,
    values_from = Cell_Display,
    values_fill = "Path Excluded"
  ) %>%
  arrange(H_code)

write_csv(table3_wide, file.path(tbl_out_dir, "Table3_Hypothesis_Evaluation_Summary.csv"))

# Render Formatted gt Table
table3_gt <- table3_wide %>%
  gt() %>%
  tab_header(
    title = md("**Table 3. Empirical Evaluation of Structural Path Hypotheses**"),
    subtitle = "Path significance (% Sig at p < 0.05) and directional support across top-ranked candidate models"
  ) %>%
  cols_label(
    H_code               = md("**Code**"),
    path_label           = md("**Target Structural Path**"),
    Baseline_Expectation = md("**Baseline Expectation**"),
    DAG1A_long           = md("**DAG 1A (Long)**"),
    DAG1A_short          = md("**DAG 1A (Short)**"),
    DAG1B_long           = md("**DAG 1B (Long)**"),
    DAG1B_short          = md("**DAG 1B (Short)**")
  ) %>%
  cols_align(align = "center", columns = c(H_code, Baseline_Expectation, DAG1A_long, DAG1A_short, DAG1B_long, DAG1B_short)) %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(13),
    column_labels.font.weight = "bold"
  )

# Save Outputs
gtsave(table3_gt, file.path(tbl_out_dir, "Table3_Hypothesis_Evaluation_Summary.html"))
gtsave(
  data     = table3_gt,
  filename = file.path(tbl_out_dir, "Table3_Hypothesis_Evaluation_Summary.png"),
  vwidth   = 1200,
  vheight  = 700
)

message("Goal 2 Pipeline Complete! Table 3 and diagnostics exported to: ", master_out)