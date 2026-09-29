# ==============================================================================
# Script: 02_goal2_table4_indicator_diagnostics.R
# Directory: Rcode_for_paper/03_goal2_hypothesis_testing/
# Purpose: Goal 2 - Table 4 Indicator Diagnostics (dAIC <= 3, Single Indicator per Row)
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

library(tidyverse)
library(gt)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DEFINE HYPOTHESES
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26"

master_out   <- file.path(proj_dir, "Routput_for_paper")
tbl_out_dir  <- file.path(proj_dir, "tables")
data_out_dir <- file.path(proj_dir, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

# Source Universal Crosswalk Utility
source(file.path("Rcode_for_paper", "01_data_prep", "00b_crosswalk_utility_fxn.r"))

# Single Source of Truth
crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")
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

cat("Building Table 4 (dAIC <= 3 Criterion, Individual Indicators per Row)...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD RANKINGS, APPLY dAIC <= 3 & PURGE CR SEALS / HCI
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

# Filter by dAIC <= 3 criterion per model run
daic3_keys_clean <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  filter(dAIC <= 3.0) %>%
  select(model_id, modNum, dAIC) %>%
  ungroup()

cat(sprintf("Retained %d model evaluations across all runs under dAIC <= 3.\n", nrow(daic3_keys_clean)))

# ------------------------------------------------------------------------------
# STEP 2: LOAD ESTIMATES & CALCULATE DIRECTIONAL SUPPORT METRICS
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
  inner_join(daic3_keys_clean, by = c("model_id", "modNum")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0),
    is_competitor  = grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    
    Expected_Sign = case_when(
      LisaName %in% c("05.ForageFishNCC_DFA1", "X05_DFA_abundSardine") ~ "Negative (-)",
      is_competitor                                                   ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK")                                 ~ "Negative (-)",
      TRUE                                                            ~ "Positive (+)"
    ),
    
    # Classify whether run directionally supports or flips expectation
    Direction_Status = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 ~ "Supported",
      Expected_Sign == "Negative (-)" & is_negative == 1 ~ "Supported",
      TRUE                                               ~ "Flipped"
    )
  )

# Save indicator-level evaluation CSV
write_csv(indicator_analysis, file.path(data_out_dir, "Goal2_Indicator_Evaluations_dAIC3.csv"))

# ------------------------------------------------------------------------------
# STEP 3: AGGREGATE PER INDIVIDUAL INDICATOR (SEPARATE SUPPORTED VS FLIPPED)
# ------------------------------------------------------------------------------
table4_data <- indicator_analysis %>%
  group_by(H_code, path_label, rhs, LisaName) %>%
  summarise(
    N_Models = n(),
    
    # Metrics for Supported Direction
    N_Supp       = sum(Direction_Status == "Supported"),
    Pct_Sig_Supp = if_else(N_Supp > 0, round((sum(Direction_Status == "Supported" & is_significant == 1) / N_Supp) * 100, 1), NA_real_),
    Beta_Supp    = if_else(N_Supp > 0, round(mean(est[Direction_Status == "Supported"], na.rm = TRUE), 3), NA_real_),
    
    # Metrics for Flipped Direction
    N_Flip       = sum(Direction_Status == "Flipped"),
    Pct_Sig_Flip = if_else(N_Flip > 0, round((sum(Direction_Status == "Flipped" & is_significant == 1) / N_Flip) * 100, 1), NA_real_),
    Beta_Flip    = if_else(N_Flip > 0, round(mean(est[Direction_Status == "Flipped"], na.rm = TRUE), 3), NA_real_),
    
    .groups = "drop"
  ) %>%
  rename(SEM_Node = rhs, Indicator = LisaName) %>%
  arrange(H_code, SEM_Node, Indicator)

write_csv(table4_data, file.path(tbl_out_dir, "Table4_Indicator_Diagnostics_dAIC3.csv"))

# ------------------------------------------------------------------------------
# STEP 4: RENDER PUBLICATION GT TABLE
# ------------------------------------------------------------------------------
table4_gt <- table4_data %>%
  select(H_code, Indicator, N_Models, Pct_Sig_Supp, Beta_Supp, Pct_Sig_Flip, Beta_Flip) %>%
  gt(groupname_col = "H_code") %>%
  tab_header(
    title = md("**Table 4. Individual Indicator Diagnostics Across Structural Paths**"),
    subtitle = "Models evaluated under ΔAIC ≤ 3. % Sig and Mean β split separately by direction."
  ) %>%
  tab_spanner(
    label = md("**Supported Direction**"),
    columns = c(Pct_Sig_Supp, Beta_Supp)
  ) %>%
  tab_spanner(
    label = md("**Flipped Sign Direction**"),
    columns = c(Pct_Sig_Flip, Beta_Flip)
  ) %>%
  cols_label(
    Indicator    = md("**Indicator (LisaName)**"),
    N_Models     = md("**N**"),
    Pct_Sig_Supp = md("**% Sig**"),
    Beta_Supp    = md("**Mean β**"),
    Pct_Sig_Flip = md("**% Sig**"),
    Beta_Flip    = md("**Mean β**")
  ) %>%
  cols_align(align = "left", columns = c(Indicator)) %>%
  cols_align(align = "center", columns = c(N_Models, Pct_Sig_Supp, Beta_Supp, Pct_Sig_Flip, Beta_Flip)) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  
  # Highlight high significance in supported direction (soft green)
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(
      columns = c(Pct_Sig_Supp, Beta_Supp),
      rows = !is.na(Pct_Sig_Supp) & Pct_Sig_Supp >= 50
    )
  ) %>%
  
  # Highlight high significance in flipped direction (soft red)
  tab_style(
    style = list(cell_fill(color = "#FFEBEE"), cell_text(color = "#C62828", weight = "bold")),
    locations = cells_body(
      columns = c(Pct_Sig_Flip, Beta_Flip),
      rows = !is.na(Pct_Sig_Flip) & Pct_Sig_Flip >= 50
    )
  ) %>%
  
  cols_width(
    Indicator ~ px(280),
    N_Models ~ px(50),
    everything() ~ px(85)
  ) %>%
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = px(12),
    heading.subtitle.font.size = px(10),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    data_row.padding = px(3),
    table.width = px(700)
  )

print(table4_gt)

# Save Outputs
gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Indicator_Diagnostics_dAIC3.html"))
gtsave(
  data     = table4_gt,
  filename = file.path(tbl_out_dir, "Table4_Indicator_Diagnostics_dAIC3.png"),
  vwidth   = 800,
  vheight  = 750
)

message("Table 4 complete! Saved to: ", tbl_out_dir)
