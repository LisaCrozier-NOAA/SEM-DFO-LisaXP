# ==============================================================================
# Script: 06_generate_table2_hypotheses.R
# Purpose: Generate Table 2 (A Priori Hypotheses) with Category 4 (H4 Regional)
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

library(gt)
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Master Output Directory Setup
# ------------------------------------------------------------------------------
out_dir     <- "Rcode_for_paper/Routput_for_paper"
tbl_out_dir <- file.path(out_dir, "tables")
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Build Hypotheses Data Frame (Includes H4 Regional Bottleneck)
# ------------------------------------------------------------------------------
hypotheses_df <- tribble(
  ~Code, ~Category, ~Ecological_Mechanism, ~Target_Path, ~Expected_Sign, ~Evaluated_DAGs, ~Primary_Reference,
  
  # Category 1: Bottom-Up
  "H1a", "H1: Bottom-up", "Early marine prey availability promotes juvenile growth", "PreyNCC -> Growth", "+", "DAG1A, DAG1B", "Crozier et al. (2021)",
  "H1b", "H1: Bottom-up", "Early marine growth enhances early ocean abundance (CPUE)", "Growth -> CPUE", "+", "DAG1A", "Wells et al. (2020)",
  "H1c", "H1: Bottom-up", "Early marine growth directly increases smolt-to-adult survival", "Growth -> SAR", "+", "DAG1B", "Crozier et al. (2021)",
  "H1e", "H1: Bottom-up", "Subarctic marine productivity enhances late marine survival", "PreyAK -> SAR", "+", "DAG1A, DAG1B", "Wells et al. (2023)",
  
  # Category 2: Top-Down
  "H2a", "H2: Top-down", "Early ocean predator pressure reduces juvenile abundance", "PredNCC -> CPUE", "-", "DAG1A, DAG1B", "Wells et al. (2025)",
  "H2c", "H2: Top-down", "Late marine predation acts as a bottleneck on adult survival", "PredAK -> SAR", "-", "DAG1A, DAG1B", "Crozier et al. (2025)",
  
  # Category 3: Indirect & Moderated
  "H3a", "H3: Top-down Moderated", "Alternate prey buffer early juvenile salmon from NCC predators", "PredNCC * AltPreyNCC -> Abundance", "+", "DAG1D (Two-Stage)", "Wells et al. (2023)",
  "H3b", "H3: Top-down Moderated", "Buffer prey induce functional switching by sea lions/sharks", "PredAK * AltPreyAK -> SAR", "+", "DAG1D (Two-Stage)", "Wells et al. (2023)",
  "H3c", "H3: Top-down Moderated", "Thermal conditions alter spatial overlap and predation intensity", "Pred * SST -> Survival", "+/-", "DAG1D (Two-Stage)", "Wells et al. (2025)",
  
  # Category 4: Regional Life-Stage Dominance
  "H4a", "H4: Regional Bottleneck", "Early ocean residence (NCC) sets survival trajectory prior to AK stage", "CPUE -> SAR", "+ (Rel. Magnitude)", "DAG1A, DAG1B, DAG1D", "Wells et al. (2020)"
)

# ------------------------------------------------------------------------------
# 3. Generate Formatted gt Table
# ------------------------------------------------------------------------------
table2_gt <- hypotheses_df %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 2. Proposed Ecological & Regional Hypotheses**"),
    subtitle = "Categorized by Trophic Control (H1-H3) and Regional Life-Stage Bottlenecks (H4)"
  ) %>%
  cols_label(
    Code = md("**Code**"),
    Ecological_Mechanism = md("**Ecological Mechanism**"),
    Target_Path = md("**Target Structural Path**"),
    Expected_Sign = md("**Expected Sign**"),
    Evaluated_DAGs = md("**Target DAGs**"),
    Primary_Reference = md("**Supporting Literature**")
  ) %>%
  cols_align(align = "center", columns = c(Code, Expected_Sign, Evaluated_DAGs)) %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(13),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold"
  )

# ------------------------------------------------------------------------------
# 4. Save HTML and Image (PNG) Outputs to Rcode_for_paper/Routput_for_paper/tables/
# ------------------------------------------------------------------------------
# Save as HTML
gtsave(table2_gt, file.path(tbl_out_dir, "Table2_Hypotheses_Matrix.html"))

# Save as PNG (for PowerPoint / Slides)
gtsave(
  data = table2_gt, 
  filename = file.path(tbl_out_dir, "Table2_Hypotheses_Matrix.png"),
  vwidth = 1200,
  vheight = 800
)

message("Table 2 exported to: ", tbl_out_dir)