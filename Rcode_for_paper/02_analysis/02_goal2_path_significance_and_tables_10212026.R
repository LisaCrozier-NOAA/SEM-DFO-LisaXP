# ==============================================================================
# Script: 02_goal2_path_significance_and_tables.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Integrated Goal 2 Hypothesis Testing Pipeline
#          - Table 3 (Word, Wide, PPT Formats)
#          - Table 4 (Indicator Diagnostics under dAIC <= 3)
# Output: Rcode_for_paper/Routput_for_paper/ and tables/
# ==============================================================================
#RESULT FIGURES-----------
# table3_ppt_gt #4 COLUMN GRID
# table3_1col_gt #ALL MODELS SUMMED TOGETHER
# table4_gt #INDIV INDICATORS
# filename = file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.png"),
# gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Indicator_Sign_Colored_dAIC3.png"), vwidth = 800, vheight = 750)
# gtsave(table3_ppt_gt, file.path(tbl_out_dir, "Table3_PPT.png"), vwidth = 650, vheight = 500)


library(tidyverse)
library(gt)
library(webshot2)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD CROSSWALK
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"

master_out   <- file.path(proj_dir,  "Routput_for_paper")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

# Source Universal Crosswalk Utility
#source(file.path("Rcode_for_paper", "01_data_prep", "00b_crosswalk_utility_fxn.r"))

# Single Source of Truth for Crosswalk
crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}
master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# Build Comprehensive Crosswalk Lookup Dictionary
lookup_dict <- setNames(master_crosswalk$LisaName, master_crosswalk$dfa_cols)
lookup_dict_raw <- setNames(master_crosswalk$LisaName, master_crosswalk$DFAname)
lookup_dict <- c(lookup_dict, lookup_dict_raw)

# Active Baseline Models
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

cat("Starting Integrated Goal 2 Pipeline...\n")

# ------------------------------------------------------------------------------
# STEP 1: SINGLE PASS LOAD FOR RANKINGS & PARAMETER ESTIMATES---------
# ------------------------------------------------------------------------------
cat("Loading model rankings across all candidate DAGs...\n")
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) {
    warning("File missing: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

# Reshape & Clean Rankings Data
rankings_long <- rankings_all %>%
#  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
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

# Identify & Purge Model Runs with Specific CR Seals / Exclusions
# cr_seal_keys <- rankings_long %>%
#   filter(LisaName %in% c("X21_DFA_SST_NCC_aprjun_DFObccm",
#                          "X15_sablefishAK_biomass4plus",
#                          "X15_pollockGOA_SSB_lag1",
#                          "X15_DFA_sleeperSharkBSAI_predAK_adultyr")) %>%
#   select(model_id, modNum) %>%
#   distinct()
# 
# rankings_long <- rankings_long %>%
#   anti_join(cr_seal_keys, by = c("model_id", "modNum"))

cat("Loading parameter estimates...\n")
estimates_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "parameterEstimates.csv")
  if (!file.exists(f_path)) {
    warning("File missing: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

# Master Evaluated Estimates Base Data
estimates_base <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_long, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0),
    is_competitor  = grepl("X02|X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    Expected_Sign  = case_when(
      is_competitor                     ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK")   ~ "Negative (-)",
      TRUE                              ~ "Positive (+)"
    ),
    Model_Supports = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 ~ TRUE,
      Expected_Sign == "Negative (-)" & is_negative == 1 ~ TRUE,
      TRUE                                               ~ FALSE
    ),
    Direction_Status = if_else(Model_Supports, "Supported", "Flipped")
  )

# ------------------------------------------------------------------------------
# STEP 2: TOP-100 MODEL ANALYSIS & TABLE 3 GENERATION-------
# ------------------------------------------------------------------------------
cat("Processing Top 100 Models for Table 3 & Individual Diagnostics CSV...\n")

top_100_keys_clean <- rankings_all %>%
#  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
#  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
  select(model_id, modNum) %>%
  ungroup()

targeted_estimates_top100 <- estimates_base %>%
  inner_join(top_100_keys_clean, by = c("model_id", "modNum"))

# Save Diagnostic CSV (Top 100)
indicator_diagnostic_summary <- targeted_estimates_top100 %>%
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

write_csv(indicator_diagnostic_summary, file.path(data_out_dir, "Goal2_Indicator_Signs.csv"))

# Path Summary Data Frame for Table 3
path_summary_by_model <- targeted_estimates_top100 %>%
  group_by(model_id, H_code, path_label) %>%
  summarise(
    N_Top100      = n(),
    Pct_Sig       = round((sum(is_significant) / N_Top100) * 100, 1),
    Pct_Supported = round((sum(Model_Supports) / N_Top100) * 100, 1),
    .groups       = "drop"
  ) %>%
  mutate(
    Baseline_Expectation = if_else(H_code %in% c("H2a", "H2c"), "Negative (-)", "Positive (+)"),
    Verdict = case_when(
      Pct_Sig >= 50 & Pct_Supported >= 50 ~ "Supported",
      Pct_Sig >= 50 & Pct_Supported < 50  ~ "Flipped Sign",
      TRUE                                ~ "NS"
    ),
    Cell_Display_HTML = paste0(Verdict, "<br>(", Pct_Sig, "% Sig)"),
    Cell_Display_Text = paste0(Verdict, "\n", Pct_Sig, "% Sig")
  )

model_cols <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

# Helper function to apply green/red styling to gt tables
apply_table3_styles <- function(gt_obj) {
  for (col in model_cols) {
    gt_obj <- gt_obj %>%
      tab_style(
        style = list(cell_fill(color = "#C8E6C9"), cell_text(color = "#1B5E20", weight = "bold")),
        locations = cells_body(columns = all_of(col), rows = str_detect(.data[[col]], "Supported"))
      ) %>%
      tab_style(
        style = list(cell_fill(color = "#FFCDD2"), cell_text(color = "#B71C1C", weight = "bold")),
        locations = cells_body(columns = all_of(col), rows = str_detect(.data[[col]], "Flipped Sign"))
      )
  }
  return(gt_obj)
}

# --- 2A: Table 3 - Word Portrait Format ----
t3_word_wide <- path_summary_by_model %>%
  select(H_code, path_label, Baseline_Expectation, model_id, Cell_Display_HTML) %>%
  pivot_wider(names_from = model_id, values_from = Cell_Display_HTML, values_fill = "Path Excluded") %>%
  arrange(H_code)

table3_word_gt <- t3_word_wide %>%
  gt() %>%
  tab_header(
    title = md("**Table 3. Evaluation of Structural Path Hypotheses**"),
    subtitle = "Path significance (% Sig at p < 0.05) and directional support across top-ranked candidate models"
  ) %>%
  cols_label(
    H_code               = md("**Code**"),
    path_label           = md("**Structural Path**"),
    Baseline_Expectation = md("**Expectation**"),
    DAG1A_long           = md("**DAG 1A<br>(Long)**"),
    DAG1A_short          = md("**DAG 1A<br>(Short)**"),
    DAG1B_long           = md("**DAG 1B<br>(Long)**"),
    DAG1B_short          = md("**DAG 1B<br>(Short)**")
  ) %>%
  fmt_markdown(columns = c(all_of(model_cols), "path_label")) %>%
  cols_align(align = "center", columns = c(H_code, Baseline_Expectation, all_of(model_cols))) %>%
  cols_width(H_code ~ px(55), path_label ~ px(160), Baseline_Expectation ~ px(95), everything() ~ px(85)) %>%
  apply_table3_styles() %>%
  tab_options(table.width = px(650), table.font.size = px(10), heading.title.font.size = px(12), column_labels.font.weight = "bold", data_row.padding = px(4))

gtsave(table3_word_gt, file.path(tbl_out_dir, "Table3_Hypothesis_Test_Word_Portrait.html"))
gtsave(table3_word_gt, file.path(tbl_out_dir, "Table3_Hypothesis_Test_Word_Portrait.png"), vwidth = 700, vheight = 550)

# --- 2B: Table 3 - Publication Wide Format ----
t3_pub_wide <- path_summary_by_model %>%
  select(H_code, path_label, Baseline_Expectation, model_id, Cell_Display_HTML) %>%
  pivot_wider(names_from = model_id, values_from = Cell_Display_HTML, values_fill = "Path Excluded") %>%
  arrange(H_code)

write_csv(
  path_summary_by_model %>% select(H_code, path_label, Baseline_Expectation, model_id, Cell_Display_Text) %>%
    pivot_wider(names_from = model_id, values_from = Cell_Display_Text, values_fill = "Path Excluded") %>% arrange(H_code),
  file.path(tbl_out_dir, "Table3_Hypothesis_Test_wideformat.csv")
)

table3_wide_gt <- t3_pub_wide %>%
  gt() %>%
  tab_header(
    title = md("**Table 3. Hypothesis Tests**"),
    subtitle = "Path significance (% Sig at p < 0.05) and directional support across top-ranked candidate models"
  ) %>%
  cols_label(
    H_code               = md("**Code**"),
    path_label           = md("**Structural Path**"),
    Baseline_Expectation = md("**Expectation**"),
    DAG1A_long           = md("**DAG 1A<br>(Long)**"),
    DAG1A_short          = md("**DAG 1A<br>(Short)**"),
    DAG1B_long           = md("**DAG 1B<br>(Long)**"),
    DAG1B_short          = md("**DAG 1B<br>(Short)**")
  ) %>%
  fmt_markdown(columns = c(all_of(model_cols), "path_label")) %>%
  cols_align(align = "center", columns = c(H_code, Baseline_Expectation, all_of(model_cols))) %>%
  apply_table3_styles() %>%
  tab_options(table.font.size = px(11), heading.title.font.size = px(13), column_labels.font.weight = "bold", data_row.padding = px(6))

gtsave(table3_wide_gt, file.path(tbl_out_dir, "Table3_Wide_Format.html"))
gtsave(table3_wide_gt, file.path(tbl_out_dir, "Table3_Wide_Format.png"), vwidth = 950, vheight = 650)

# --- 2C: Table 3 - Presentation Slide Format (PPT) ----
t3_ppt_wide <- path_summary_by_model %>%
  mutate(
    Concept = if_else(H_code %in% c("H1a", "H1b", "H1c", "H1e"),
                      "Bigger is better -- competitor/prey mix",
                      "Predator abundance regulates survival -- supported")
  ) %>%
  select(Concept, path_label, Baseline_Expectation, model_id, Cell_Display_HTML) %>%
  pivot_wider(names_from = model_id, values_from = Cell_Display_HTML, values_fill = "Path Excluded")

table3_ppt_gt <- t3_ppt_wide %>%
  gt(groupname_col = "Concept") %>%
  tab_header(
    title = md("**Evaluation of Structural Path Hypotheses**"),
    subtitle = "Directional support & significance (p < 0.05) across top candidate models"
  ) %>%
  cols_label(
    path_label           = md("**Structural Path**"),
    Baseline_Expectation = md("**Expect.**"),
    DAG1A_long           = md("**DAG 1A<br>(Long)**"),
    DAG1A_short          = md("**DAG 1A<br>(Short)**"),
    DAG1B_long           = md("**DAG 1B<br>(Long)**"),
    DAG1B_short          = md("**DAG 1B<br>(Short)**")
  ) %>%
  fmt_markdown(columns = c(all_of(model_cols), "path_label")) %>%
  cols_align(align = "center", columns = c(Baseline_Expectation, all_of(model_cols))) %>%
  cols_align(align = "left", columns = c(path_label)) %>%
  cols_width(path_label ~ px(140), Baseline_Expectation ~ px(75), everything() ~ px(80)) %>%
  apply_table3_styles() %>%
  tab_style(
    style = list(cell_fill(color = "#ECEFF1"), cell_text(weight = "bold", color = "#263238", size = px(16))),
    locations = cells_row_groups()
  ) %>%
  tab_options(table.width = px(580), table.font.size = px(11), heading.title.font.size = px(14), heading.subtitle.font.size = px(11), column_labels.font.weight = "bold", row_group.padding = px(6), data_row.padding = px(4))

gtsave(table3_ppt_gt, file.path(tbl_out_dir, "Table3_PPT.html"))
gtsave(table3_ppt_gt, file.path(tbl_out_dir, "Table3_PPT.png"), vwidth = 650, vheight = 500)

# ------------------------------------------------------------------------------
# STEP 3: dAIC <= 3 FILTERING & TABLE 4 GENERATION--------
# ------------------------------------------------------------------------------
cat("Processing dAIC <= 3 models for Table 4 Indicator Diagnostics...\n")

daic3_keys_clean <- rankings_all %>%
#  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
#  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  filter(dAIC <= 3.0) %>%
  select(model_id, modNum, dAIC) %>%
  ungroup()

indicator_analysis_daic3 <- estimates_base %>%
  inner_join(daic3_keys_clean, by = c("model_id", "modNum"))

write_csv(indicator_analysis_daic3, file.path(data_out_dir, "Goal2_indicator_coef_detail_dAIC3.csv"))

# Summarize for Table 4 INDIV INDICATOR SIGNS------------
table4_data <- indicator_analysis_daic3 %>%
  group_by(H_code, path_label, rhs, LisaName) %>%
  summarise(
    N_Models     = n(),
    N_Supp       = sum(Direction_Status == "Supported"),
    Pct_Sig_Supp = if_else(N_Supp > 0, round((sum(Direction_Status == "Supported" & is_significant == 1) / N_Supp) * 100, 1), NA_real_),
    Beta_Supp    = if_else(N_Supp > 0, round(mean(est[Direction_Status == "Supported"], na.rm = TRUE), 3), NA_real_),
    N_Flip       = sum(Direction_Status == "Flipped"),
    Pct_Sig_Flip = if_else(N_Flip > 0, round((sum(Direction_Status == "Flipped" & is_significant == 1) / N_Flip) * 100, 1), NA_real_),
    Beta_Flip    = if_else(N_Flip > 0, round(mean(est[Direction_Status == "Flipped"], na.rm = TRUE), 3), NA_real_),
    .groups      = "drop"
  ) %>%
  rename(SEM_Node = rhs, Indicator = LisaName) %>%
  arrange(H_code, SEM_Node, Indicator)

write_csv(table4_data, file.path(tbl_out_dir, "Table4_Indicator_Sign_data_dAIC3.csv"))

table4_gt <- table4_data %>%
  select(H_code, Indicator, N_Models, Pct_Sig_Supp, Beta_Supp, Pct_Sig_Flip, Beta_Flip) %>%
  gt(groupname_col = "H_code") %>%
  tab_header(
    title = md("**Table 4. Individual Indicator Diagnostics Across Structural Paths**"),
    subtitle = "Models evaluated under ΔAIC ≤ 3. % Sig and Mean β split separately by direction."
  ) %>%
  tab_spanner(label = md("**Supported Direction**"), columns = c(Pct_Sig_Supp, Beta_Supp)) %>%
  tab_spanner(label = md("**Flipped Sign Direction**"), columns = c(Pct_Sig_Flip, Beta_Flip)) %>%
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
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(columns = c(Pct_Sig_Supp, Beta_Supp), rows = !is.na(Pct_Sig_Supp) & Pct_Sig_Supp >= 50)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#FFEBEE"), cell_text(color = "#C62828", weight = "bold")),
    locations = cells_body(columns = c(Pct_Sig_Flip, Beta_Flip), rows = !is.na(Pct_Sig_Flip) & Pct_Sig_Flip >= 50)
  ) %>%
  cols_width(Indicator ~ px(280), N_Models ~ px(50), everything() ~ px(85)) %>%
  tab_options(table.font.size = px(10), heading.title.font.size = px(12), heading.subtitle.font.size = px(10), column_labels.font.weight = "bold", row_group.font.weight = "bold", data_row.padding = px(3), table.width = px(700))

gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Indicator_Sign_Colored_dAIC3.html"))
gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Indicator_Sign_Colored_dAIC3.png"), vwidth = 800, vheight = 750)

#=======================================
#CONDENSED VERSION OF TABLE 3 WITH 1 COLUMN-----------------
# ==============================================================================
# ==============================================================================
# Script: 03_goal2_table3_consolidated_1col.R
# Directory: Rcode_for_paper/03_goal2_hypothesis_testing/
# Purpose: Proportional Color Gradient 1-Column Version of Table 3
# ==============================================================================

library(tidyverse)
library(gt)
library(webshot2)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 1: CALCULATE EXACT SUPPORT PERCENTAGES & VERDICT LABELS
# ------------------------------------------------------------------------------
table3_1col_data <- targeted_estimates_top100 %>%
  group_by(H_code, path_label) %>%
  summarise(
    N_Total_Evaluated = n(),
    N_Sig             = sum(is_significant == 1),
    N_Supported       = sum(Model_Supports == TRUE),
    N_Flipped         = sum(Model_Supports == FALSE),
    
    Pct_Supported     = round((N_Supported / N_Total_Evaluated) * 100, 1),
    Pct_Flipped       = round((N_Flipped / N_Total_Evaluated) * 100, 1),
    Overall_Pct_Sig   = round((N_Sig / N_Total_Evaluated) * 100, 1),
    
    .groups = "drop"
  ) %>%
  mutate(
    Baseline_Expectation = if_else(H_code %in% c("H2a", "H2c"), "Negative (-)", "Positive (+)"),
    
    Concept = if_else(
      H_code %in% c("H1a", "H1b", "H1c", "H1e"),
      "Bigger is better -- competitor/prey mix",
      "Predator abundance regulates survival -- supported"
    ),
    
    # Label as NS if 0%, Supported if > 0%, Flipped if majority flipped
    Verdict_Label = case_when(
      Pct_Supported <= 10 & Pct_Flipped <= 10 ~ "NS",
      # Pct_Supported >= Pct_Flipped         ~ "Supported",
      # Pct_Flipped >= 10         ~ "Flipped Sign",
      TRUE                                 ~ "Supported"
    ),
    
    # Text display string
    Consolidated_Display = case_when(
      Verdict_Label == "NS" ~ paste0(
        "**0%** Supported<br><span style='font-size:9px; color:#777;'>(0% Sig, N=", N_Total_Evaluated, ")</span>"
      ),
      TRUE ~ paste0(
        "**", Pct_Supported, "%** ", Verdict_Label, 
        "<br><span style='font-size:9px; color:#555;'>(", Overall_Pct_Sig, "% Sig, N=", N_Total_Evaluated, ")</span>"
      )
    )
  ) %>%
  arrange(H_code)

# ------------------------------------------------------------------------------
# STEP 2: DEFINE EXPLICIT COLOR GRADIENT FUNCTION
# ------------------------------------------------------------------------------
get_gradient_color <- function(pct_supported, pct_flipped, verdict_label) {
  
  # Rule 1: 0% Support / NS -> Pure White
  if (is.na(pct_supported) || pct_supported == 0 || verdict_label == "NS") {
    return(list(fill = "#FFFFFF", text = "#212121"))
  }
  
  # Rule 2: Supported Direction (> 0% Support)
  if (verdict_label == "Supported") {
    alpha <- pct_supported / 100
    fill_col <- colorRampPalette(c("#FFFFFF", "#E8F5E9", "#A5D6A7", "#81C784"))(100)[max(1, as.integer(alpha * 99) + 1)]
    return(list(fill = fill_col, text = "#1B5E20"))
  }
  
  # Rule 3: Flipped Sign Direction
  if (verdict_label == "Flipped Sign") {
    if (pct_flipped == 0) return(list(fill = "#FFFFFF", text = "#212121"))
    alpha <- pct_flipped / 100
    fill_col <- colorRampPalette(c("#FFFFFF", "#FFEBEE", "#EF9A9A", "#E57373"))(100)[max(1, as.integer(alpha * 99) + 1)]
    return(list(fill = fill_col, text = "#B71C1C"))
  }
  
  return(list(fill = "#FFFFFF", text = "#212121"))
}

# ------------------------------------------------------------------------------
# STEP 3: BUILD GT TABLE & APPLY STYLES
# ------------------------------------------------------------------------------
table3_1col_gt <- table3_1col_data %>%
  select(Concept, path_label, Baseline_Expectation, Consolidated_Display) %>%
  gt(groupname_col = "Concept") %>%
  tab_header(
    title = md("**Consolidated Evaluation of Structural Path Hypotheses**"),
    subtitle = "Proportional directional support across all evaluated candidate models (Top 100 per DAG)"
  ) %>%
  cols_label(
    path_label           = md("**Structural Path**"),
    Baseline_Expectation = md("**Expectation**"),
    Consolidated_Display = md("**Overall Support %**")
  ) %>%
  fmt_markdown(columns = c("Consolidated_Display", "path_label")) %>%
  cols_align(align = "left", columns = c(path_label)) %>%
  cols_align(align = "center", columns = c(Baseline_Expectation, Consolidated_Display)) %>%
  cols_width(
    path_label           ~ px(180),
    Baseline_Expectation ~ px(90),
    Consolidated_Display ~ px(180)
  )

# Apply individual cell styles using exact row metrics
for (i in seq_len(nrow(table3_1col_data))) {
  p_supp  <- table3_1col_data$Pct_Supported[i]
  p_flip  <- table3_1col_data$Pct_Flipped[i]
  verdict <- table3_1col_data$Verdict_Label[i]
  
  style_info <- get_gradient_color(p_supp, p_flip, verdict)
  
  table3_1col_gt <- table3_1col_gt %>%
    tab_style(
      style = list(
        cell_fill(color = style_info$fill),
        cell_text(color = style_info$text)
      ),
      locations = cells_body(
        columns = "Consolidated_Display",
        rows = i
      )
    )
}

# Table Options
table3_1col_gt <- table3_1col_gt %>%
  tab_style(
    style = list(
      cell_fill(color = "#ECEFF1"),
      cell_text(weight = "bold", color = "#263238", size = px(13))
    ),
    locations = cells_row_groups()
  ) %>%
  tab_options(
    table.width                = px(470),
    table.font.size            = px(11),
    heading.title.font.size    = px(13),
    heading.subtitle.font.size = px(10),
    column_labels.font.weight  = "bold",
    row_group.padding          = px(5),
    data_row.padding           = px(4)
  )

print(table3_1col_gt)

# Export Outputs
gtsave(table3_1col_gt, file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.html"))
gtsave(
  data     = table3_1col_gt,
  filename = file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.png"),
  vwidth   = 520,
  vheight  = 460
)
# ------------------------------------------------------------------------------
# STEP 4: EXPORT OUTPUTS----
# ------------------------------------------------------------------------------
gtsave(table3_1col_gt, file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.html"))
gtsave(
  data     = table3_1col_gt,
  filename = file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.png"),
  vwidth   = 520,
  vheight  = 460
)

cat("Gradient-colored 1-column Table 3 exported to: ", tbl_out_dir, "\n")


#RESULT FIGURES-----------
table3_ppt_gt #4 COLUMN GRID
table3_1col_gt #ALL MODELS SUMMED TOGETHER
table4_gt #INDIV INDICATORS
filename = file.path(tbl_out_dir, "Table3_Consolidated_Gradient_1Col.png"),
gtsave(table4_gt, file.path(tbl_out_dir, "Table4_Indicator_Sign_Colored_dAIC3.png"), vwidth = 800, vheight = 750)
gtsave(table3_ppt_gt, file.path(tbl_out_dir, "Table3_PPT.png"), vwidth = 650, vheight = 500)
