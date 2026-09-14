library(gt)
library(tidyverse)
install.packages("webshot2", type = "binary")
library(webshot2)
# ------------------------------------------------------------------------------
# 1. Directories & Setup
# ------------------------------------------------------------------------------
out_dir <- "Routput_for_paper"
table_out_dir <- file.path(out_dir, "tables")
dir.create(table_out_dir, showWarnings = FALSE, recursive = TRUE)

#To convert your gt table directly into a high-resolution PNG image that you can copy and paste into PowerPoint, you can use the gtsave() function.

#Behind the scenes, saving a gt table as an image requires a web rendering engine (like chromote or webshot2).

# ------------------------------------------------------------------------------
# 2. Build Hypotheses Data Frame
# ------------------------------------------------------------------------------
hypotheses_df <- tribble(
  ~Code, ~Category, ~Ecological_Mechanism, ~Target_Path, ~Expected_Sign, ~Evaluated_DAGs, ~Primary_Reference,
  
  # Category 1: Bottom-Up
  "H1a", "H1: Bottom-up", "Early marine prey availability promotes juvenile growth", "PreyNCC -> Growth", "+", "DAG1A, DAG1B", "Crozier et al. (2021)",
  "H1b", "H1: Bottom-up", "Early marine growth enhances early ocean abundance (CPUE)", "Growth -> CPUE", "+", "DAG1A", "Wells et al. (2020)",
  "H1c", "H1: Bottom-up", "Early marine growth directly increases smolt-to-adult survival", "Growth -> SAR", "+", "DAG1B", "Crozier et al. (2021)",
  "H1d", "H1: Bottom-up", "Early marine abundance dictates adult return strength", "CPUE -> SAR", "+", "DAG1A, DAG1B", "Wells et al. (2020)",
  "H1e", "H1: Bottom-up", "Subarctic marine productivity enhances late marine survival", "PreyAK -> SAR", "+", "DAG1A, DAG1B", "Wells et al. (2023)",
  
  # Category 2: Top-Down
  "H2a", "H2: Top-down", "Early ocean predator pressure reduces juvenile abundance", "PredNCC -> CPUE", "-", "DAG1A, DAG1B", "Wells et al. (2025)",
  "H2c", "H2: Top-down", "Late marine predation acts as a bottleneck on adult survival", "PredAK -> SAR", "-", "DAG1A, DAG1B", "Crozier et al. (2025)",
  
  # Category 3: Indirect & Moderated
  "H3a", "H3: Top-down Moderated", "Alternate prey buffer early juvenile salmon from NCC predators", "PredNCC * AltPreyNCC -> Abundance", "+", "DAG1D (Two-Stage)", "Wells et al. (2023)",
  "H3b", "H3: Top-down Moderated", "Buffer prey induce functional switching by sea lions/sharks", "PredAK * AltPreyAK -> SAR", "+", "DAG1D (Two-Stage)", "Wells et al. (2023)",
  "H3c", "H3: Top-down Moderated", "Thermal conditions alter spatial overlap and predation intensity", "Pred * SST -> Survival", "+/-", "DAG1D (Two-Stage)", "Wells et al. (2025)"
)

# ------------------------------------------------------------------------------
# 3. Generate Formatted gt Table
# ------------------------------------------------------------------------------
table2_gt <- hypotheses_df %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 2. Proposed Ecological Hypotheses and Path Expectations**"),
    subtitle = "Tested across baseline SEM configurations (DAG 1A & 1B) and the two-stage interaction model (DAG 1D)"
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
# 4. Save Outputs to Routput_for_paper/tables/
# ------------------------------------------------------------------------------
gtsave(table2_gt, file.path(table_out_dir, "Table2_Hypotheses_Matrix.html"))
 gtsave(table2_gt, file.path(table_out_dir, "Table2_Hypotheses_Matrix.docx")) # Direct Word export

# Save as a high-resolution PNG for PowerPoint
gtsave(
  data = table2_gt, 
  filename = file.path(table_out_dir, "Table2_Hypotheses_Matrix.png"),
  vwidth = 1200,   # Set custom width for crisp PowerPoint slides
  vheight = 800     # Set custom height
)

message("Table PNG exported to: ", file.path(table_out_dir, "Table2_Hypotheses_Matrix.png"))
message("Table 2 saved to: ", file.path(table_out_dir, "Table2_Hypotheses_Matrix.html"))
