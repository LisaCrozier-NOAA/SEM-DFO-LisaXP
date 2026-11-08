# ==============================================================================
# Script: 03_figs_tables/04_Table_var_importance.R
# Purpose: Filter duplicate variables, renormalize importance, and export Word table
# Output: Rcode_for_paper/Routput_for_paper/tables/Table_Variable_Importance_daic3_filtered.docx
# ==============================================================================
# ==============================================================================
# Script: 03_figs_tables/03_format_importance_table_word.R
# Purpose: Filter duplicates, sort within node by max importance, and export compact Word table
# Output: Rcode_for_paper/Routput_for_paper/tables/Table_Variable_Importance_daic3_filtered.docx
# ==============================================================================

library(tidyverse)
library(gt)

# ------------------------------------------------------------------------------
# 1. Paths & Setup
# ------------------------------------------------------------------------------
csv_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/data/importance_topvar_daic3.csv"
tbl_dir  <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/tables"
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

df_raw <- read.csv(csv_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. Filter Duplicates, Recalculate Weights & Sort by Max Importance
# ------------------------------------------------------------------------------
spurious_vars <- c("X15_pollockGOA_SSB_lag1", "X15_sablefishAK_biomass4plus")
model_cols    <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")
node_order    <- c("Abundance", "Growth", "PreyNCC", "PredNCC", "PreyAK", "PredAK")

df_clean <- df_raw %>%
  select(-order, -Node2) %>%
  filter(!LisaName %in% spurious_vars) %>%
  # Renormalize weights within generic_node so each active set sums to 1.000
  group_by(generic_node) %>%
  mutate(across(all_of(model_cols), ~ {
    col_sum <- sum(.x, na.rm = TRUE)
    if (col_sum > 0) .x / col_sum else .x
  })) %>%
  ungroup() %>%
  # Compute max importance across the 4 models for sorting
  rowwise() %>%
  mutate(max_imp = max(c_across(all_of(model_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  # Preserve natural ecological node ordering, then sort within node by max_imp descending
  mutate(node_fac = factor(generic_node, levels = node_order)) %>%
  arrange(node_fac, desc(max_imp)) %>%
  # Clean variable names for publication
  mutate(
    Variable = str_replace(LisaName, "^X[0-9]{2}_", ""),
    Variable = str_replace_all(Variable, "_", " ")
  ) %>%
  select(
    Node = generic_node,
    Variable,
    Region,
    Role = Trophic,
    `1A (Long)`  = DAG1A_long,
    `1A (Short)` = DAG1A_short,
    `1B (Long)`  = DAG1B_long,
    `1B (Short)` = DAG1B_short
  )

# ------------------------------------------------------------------------------
# 3. Build gt Table (Narrow Footprint for Portrait Word Document)
# ------------------------------------------------------------------------------
importance_gt <- df_clean %>%
  gt(groupname_col = "Node") %>%
  tab_header(
    title = md("**Table X. Variable Relative Importance Across Top SEM Models (&Delta;AIC &le; 3)**"),
    subtitle = "Ranked within node by maximum relative importance weight across model sweeps"
  ) %>%
  tab_spanner(
    label = "DAG 1A",
    columns = c(`1A (Long)`, `1A (Short)`)
  ) %>%
  tab_spanner(
    label = "DAG 1B",
    columns = c(`1B (Long)`, `1B (Short)`)
  ) %>%
  cols_label(
    Variable     = md("**Variable**"),
    Region       = md("**Region**"),
    Role         = md("**Trophic Role**"),
    `1A (Long)`  = md("**Long**"),
    `1A (Short)` = md("**Short**"),
    `1B (Long)`  = md("**Long**"),
    `1B (Short)` = md("**Short**")
  ) %>%
  fmt_number(
    columns = c(`1A (Long)`, `1A (Short)`, `1B (Long)`, `1B (Short)`),
    decimals = 3
  ) %>%
  sub_missing(
    columns = c(`1A (Long)`, `1A (Short)`, `1B (Long)`, `1B (Short)`),
    missing_text = "—"
  ) %>%
  cols_align(
    align = "center",
    columns = c(Region, Role, `1A (Long)`, `1A (Short)`, `1B (Long)`, `1B (Short)`)
  ) %>%
  cols_align(
    align = "left",
    columns = Variable
  ) %>%
  # Compact column widths (Total width ~470px = ~4.9 inches; comfortably inside 6.5" page margins)
  cols_width(
    Variable     ~ px(165),
    Region       ~ px(45),
    Role         ~ px(70),
    `1A (Long)`  ~ px(45),
    `1A (Short)` ~ px(45),
    `1B (Long)`  ~ px(45),
    `1B (Short)` ~ px(45)
  ) %>%
  tab_options(
    table.font.size = px(9.5),
    heading.title.font.size = px(11.5),
    heading.subtitle.font.size = px(9.5),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    data_row.padding = px(2.5),
    table.width = px(470)
  )

# ------------------------------------------------------------------------------
# 4. Save DOCX, PNG, and CSV
# ------------------------------------------------------------------------------
docx_out <- file.path(tbl_dir, "Table_Variable_Importance_daic3_filtered.docx")
png_out  <- file.path(tbl_dir, "Table_Variable_Importance_daic3_filtered.png")
csv_out  <- file.path(tbl_dir, "importance_topvar_daic3_filtered.csv")

gtsave(importance_gt, docx_out)
gtsave(importance_gt, png_out, vwidth = 520, vheight = 850)
write.csv(df_clean, csv_out, row.names = FALSE)

message("Success! Saved compact sorted table to:\n  Word: ", docx_out, "\n  PNG:  ", png_out)