# ==============================================================================
# Script: 03_figs_tables/02_Fig_DAG_Tab_Hypo_SinglePage.R
# Purpose: Clean Spacing for DAGs (1A, 1B, 1D Full Width) + Table 2 Portrait Page
# Output: Rcode_for_paper/figures/Figure_DAG_and_Table2_Portrait_Page.png
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(gt)
library(patchwork)
library(grid)
library(png)
library(stringr)

# ------------------------------------------------------------------------------
# 1. Directory & Color Setup
# ------------------------------------------------------------------------------
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
tbl_out_dir <- file.path(proj_dir, "tables")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

color_salmon <- "#F8AFA6"
color_prey   <- "#A8E6CF"
color_pred   <- "#D7C4B7"

# ------------------------------------------------------------------------------
# 2. Table 2: Compact Hypothesis Table for Portrait Display
# ------------------------------------------------------------------------------
hypotheses_df <- tribble(
  ~Code, ~Category, ~Ecological_Mechanism, ~Target_Path, ~Expected_Sign, ~Evaluated_DAGs, ~Primary_Reference,
  
  # Category 1: Bottom-Up
  "H1a", "H1: Bottom-up", "Early-marine prey increases juvenile growth", "PreyNCC -> Growth", "+", "DAG1A, DAG1B", "Crozier et al. (2021)",
  "H1b", "H1: Bottom-up", "Early-marine growth increases early-ocean survival", "Growth -> CPUE", "+", "DAG1A", "Wells et al. (2020)",
  "H1c", "H1: Bottom-up", "Early-marine growth increases smolt-to-adult survival directly", "Growth -> SAR", "+", "DAG1B", "Crozier et al. (2021)",
  "H1e", "H1: Bottom-up", "Subarctic marine productivity increases later-marine survival", "PreyAK -> SAR", "+", "DAG1A, DAG1B", "Wells et al. (2023)",
  
  # Category 2: Top-Down
  "H2a", "H2: Top-down", "Early-ocean predators reduce early-ocean survival", "PredNCC -> CPUE", "-", "DAG1A, DAG1B", "Wells et al. (2025)",
  "H2c", "H2: Top-down", "Later-marine predators reduce later-marine survival", "PredAK -> SAR", "-", "DAG1A, DAG1B", "Crozier et al. (2025)",
  
  # Category 3: Indirect & Moderated
  "H3a", "H3: Top-down Moderated", "Alternate prey buffer early-ocean survival", "PredNCC ~ AltPreyNCC -> CPUE", "+", "DAG1D", "Wells et al. (2023)",
  "H3b", "H3: Top-down Moderated", "Alternate prey buffer later-ocean survival", "PredAK ~ AltPreyAK -> SAR", "+", "DAG1D", "Wells et al. (2023)",
  "H3c", "H3: Top-down Moderated", "Ocean temperatures alter predator impacts on survival", "Pred ~ SST -> Survival", "+/-", "DAG1D", "Wells et al. (2025)",
  
  # Category 4: Regional Life-Stage
  "H4a", "H4: Regional Bottleneck", "Early-ocean survival regulates smolt-to-adult survival", "CPUE -> SAR", "+", "DAG1A, DAG1B, DAG1D", "Wells et al. (2020)",
  "H4b", "H4: Regional Bottleneck", "Later-ocean survival regulates smolt-to-adult survival", "CPUE + PredAK -> SAR", "- (Pred)", "DAG1A, DAG1B, DAG1D", "Wells et al. (2020)"
)

table2_gt <- hypotheses_df %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 2. Ecological & Regional Hypotheses**"),
    subtitle = "Categorized by Trophic Control (H1-H3) and Regional Bottlenecks (H4)"
  ) %>%
  cols_label(
    Code = md("**Code**"),
    Ecological_Mechanism = md("**Ecological Mechanism**"),
    Target_Path = md("**Path**"),
    Expected_Sign = md("**Sign**"),
    Evaluated_DAGs = md("**Target DAGs**"),
    Primary_Reference = md("**Literature**")
  ) %>%
  cols_align(align = "center", columns = c(Code, Expected_Sign, Evaluated_DAGs)) %>%
  cols_width(
    Code ~ px(48),
    Ecological_Mechanism ~ px(240),
    Target_Path ~ px(125),
    Expected_Sign ~ px(55),
    Evaluated_DAGs ~ px(95),
    Primary_Reference ~ px(110)
  ) %>%
  tab_options(
    table.font.size = px(8.5),
    heading.title.font.size = px(10.5),
    heading.subtitle.font.size = px(9.0),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    data_row.padding = px(1.8),
    table.width = px(675)
  )

tbl_png_path <- file.path(tbl_out_dir, "Table2_Compact_Hypotheses.png")
gtsave(table2_gt, tbl_png_path, vwidth = 700, vheight = 510)

tbl_grob <- rasterGrob(readPNG(tbl_png_path), interpolate = TRUE)

# ------------------------------------------------------------------------------
# 3. Clean Standard Rectangular Node DAG Engine (1A & 1B)
# ------------------------------------------------------------------------------
render_dag_panel <- function(nodes_df, edges_df, panel_title) {
  
  box_w <- 1.10
  box_h <- 0.48
  
  if (!"x_mid_override" %in% names(edges_df)) edges_df$x_mid_override <- NA_real_
  if (!"y_mid_override" %in% names(edges_df)) edges_df$y_mid_override <- NA_real_
  
  nodes_df <- nodes_df %>%
    mutate(label_wrapped = str_wrap(label, width = 11))
  
  edge_coords <- edges_df %>%
    left_join(nodes_df %>% select(node, x_from = x, y_from = y), by = c("from" = "node")) %>%
    left_join(nodes_df %>% select(node, x_to = x, y_to = y), by = c("to" = "node")) %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      dist = sqrt(dx^2 + dy^2),
      x_start = x_from + (dx / dist) * (box_w / 1.7),
      y_start = y_from + (dy / dist) * (box_h / 1.7),
      x_end   = x_to - (dx / dist) * (box_w / 1.7),
      y_end   = y_to - (dy / dist) * (box_h / 1.7),