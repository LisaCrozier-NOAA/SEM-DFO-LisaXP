# ==============================================================================
# Script: 03_figs_tables/01_conceptual_dags.R
# Purpose: Render revised DAG 1A, 1B, & 1D + Table 2 (Optimized for MS Word Portrait)
# Output: Rcode_for_paper/
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(gt)
library(patchwork)
library(grid)
library(png)
library(stringr)

# ------------------------------------------------------------------------------
# 1. Directory & Color Palette Setup
# ------------------------------------------------------------------------------
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
tbl_out_dir <- file.path(proj_dir, "tables")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# Color Palette
color_salmon <- "#F8AFA6"  # Soft Pink
color_prey   <- "#A8E6CF"  # Soft Light Green
color_pred   <- "#D7C4B7"  # Soft Light Brown

# ------------------------------------------------------------------------------
# 2. Table 2: A Priori Hypotheses (PNG Export & Screen Print)
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
  "H4b", "H4: Regional Bottleneck", "Later-ocean survival regulates smolt-to-adult survival", "CPUE + PredAK -> SAR", "- (Predators)", "DAG1A, DAG1B, DAG1D", "Wells et al. (2020)"
)

table2_gt <- hypotheses_df %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 2. Ecological & Regional Hypotheses**"),
    subtitle = "Categorized by Trophic Control (H1-H3) and Regional/Life-Stage Influences (H4)"
  ) %>%
  cols_label(
    Code = md("**Code**"),
    Ecological_Mechanism = md("**Ecological Mechanism**"),
    Target_Path = md("**Structural Path**"),
    Expected_Sign = md("**Expected Sign**"),
    Evaluated_DAGs = md("**Target DAGs**"),
    Primary_Reference = md("**Supporting Literature**")
  ) %>%
  cols_align(align = "center", columns = c(Code, Expected_Sign, Evaluated_DAGs)) %>%
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(13),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    table.width = px(850)
  )

# Save Table 2 as PNG
tbl_png_path <- file.path(tbl_out_dir, "Table2_Hypotheses_Matrix.png")
gtsave(table2_gt, tbl_png_path, vwidth = 1000, vheight = 850)

# Render Table to screen/viewer device
tbl_img <- readPNG(tbl_png_path)
grid.newpage()
grid.raster(tbl_img)

# ------------------------------------------------------------------------------
# 3. Standard Rectangular Node DAG Engine (DAG 1A & 1B)
# ------------------------------------------------------------------------------
render_dag_panel <- function(nodes_df, edges_df, panel_title) {
  
  box_w <- 0.95
  box_h <- 0.55
  
  # Ensure override columns exist even if omitted in edges_df
  if (!"x_mid_override" %in% names(edges_df)) edges_df$x_mid_override <- NA_real_
  if (!"y_mid_override" %in% names(edges_df)) edges_df$y_mid_override <- NA_real_
  
  # Wrap text inside node labels
  nodes_df <- nodes_df %>%
    mutate(label_wrapped = str_wrap(label, width = 10))
  
  edge_coords <- edges_df %>%
    left_join(nodes_df %>% select(node, x_from = x, y_from = y), by = c("from" = "node")) %>%
    left_join(nodes_df %>% select(node, x_to = x, y_to = y), by = c("to" = "node")) %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      dist = sqrt(dx^2 + dy^2),
      # Truncate arrows at rect borders cleanly
      x_start = x_from + (dx / dist) * (box_w / 1.8),
      y_start = y_from + (dy / dist) * (box_h / 1.8),
      x_end   = x_to - (dx / dist) * (box_w / 1.8),
      y_end   = y_to - (dy / dist) * (box_h / 1.8),
      # Explicit x and y overrides with fallback to midpoints
      x_mid   = ifelse(!is.na(x_mid_override), x_mid_override, (x_start + x_end) / 2),
      y_mid   = ifelse(!is.na(y_mid_override), y_mid_override, (y_start + y_end) / 2),
      # Offset labels vertically from lines
      y_label = y_mid + ifelse(dy == 0, 0.22, 0.05 * sign(dy))
    )
  
  p <- ggplot() +
    # Edge Path Arrows
    geom_segment(
      data = edge_coords,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
      color = "gray25", linewidth = 1.1
    ) +
    # Path Hypothesis Labels
    geom_label(
      data = edge_coords,
      aes(x = x_mid, y = y_label, label = hyp_code),
      fill = alpha("white", 0.85), color = "navy", fontface = "bold", size = 4.2,
      label.padding = unit(0.12, "lines"), label.size = NA
    ) +
    # Rectangular Variable Nodes
    geom_tile(
      data = nodes_df,
      aes(x = x, y = y, fill = fill_color),
      width = box_w, height = box_h, color = "black", linewidth = 1.0
    ) +
    # Wrapped Node Text Labels
    geom_text(
      data = nodes_df,
      aes(x = x, y = y, label = label_wrapped),
      color = "black", fontface = "bold", size = 4.0, lineheight = 0.88
    ) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.1, 6.3)) +
    scale_y_continuous(limits = c(0.2, 3.8)) +
    labs(title = panel_title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(8, 8, 8, 8),
      panel.background = element_rect(fill = "white", color = "gray80", linewidth = 0.8)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# 4. DAG 1A & 1B Node and Edge Definitions
# ------------------------------------------------------------------------------
# DAG 1A: Nodes pulled wider apart horizontally
nodes_1a <- tribble(
  ~node,      ~x,   ~y,   ~label,               ~fill_color,
  "PreyNCC",  0.7,  3.2,  "Prey (NCC)",         color_prey,
  "PredNCC",  0.7,  0.8,  "Predators (NCC)",    color_pred,
  "Growth",   2.1,  2.0,  "Salmon Growth",      color_salmon,
  "CPUE",     3.5,  2.0,  "Salmon Abundance",   color_salmon,
  "SAR",      4.9,  2.0,  "Salmon SAR",         color_salmon,
  "PreyAK",   5.7,  3.2,  "Prey (AK)",          color_prey,
  "PredAK",   5.7,  0.8,  "Predators (AK)",     color_pred
)

edges_1a <- tribble(
  ~from,      ~to,      ~hyp_code,  ~x_mid_override, ~y_mid_override,
  "PreyNCC",  "Growth", "H1a (+)",  1.10,            2.85,
  "PredNCC",  "CPUE",   "H2a (-)",  1.30,            1.15,
  "Growth",   "CPUE",   "H1b (+)",  NA,              NA,
  "CPUE",     "SAR",    "H4a (+)",  NA,              NA,
  "PreyAK",   "SAR",    "H1e (+)",  5.30,            2.85,
  "PredAK",   "SAR",    "H2c (-)",  5.30,            1.15
)

# DAG 1B: Refined edge positioning for steep vertical paths (H1a and H2a offset cleanly)
nodes_1b <- tribble(
  ~node,      ~x,   ~y,   ~label,               ~fill_color,
  "PreyNCC",  0.7,  3.2,  "Prey (NCC)",         color_prey,
  "Growth",   2.2,  2.7,  "Salmon Growth",      color_salmon,
  "PredNCC",  0.7,  0.8,  "Predators (NCC)",    color_pred,
  "CPUE",     2.2,  1.3,  "Salmon Abundance",   color_salmon,
  "SAR",      4.3,  2.0,  "Salmon SAR",         color_salmon,
  "PreyAK",   5.7,  3.2,  "Prey (AK)",          color_prey,
  "PredAK",   5.7,  0.8,  "Predators (AK)",     color_pred
)

edges_1b <- tribble(
  ~from,      ~to,      ~hyp_code,  ~x_mid_override, ~y_mid_override,
  "PreyNCC",  "Growth", "H1a (+)",  1.45,            2.90,
  "PredNCC",  "CPUE",   "H2a (-)",  1.45,            1.10,
  "Growth",   "SAR",    "H1c (+)",  NA,              NA,
  "CPUE",     "SAR",    "H4a (+)",  NA,              NA,
  "PreyAK",   "SAR",    "H1e (+)",  5.20,            2.85,
  "PredAK",   "SAR",    "H2c (-)",  5.20,            1.15
)

p_1a <- render_dag_panel(nodes_1a, edges_1a, "DAG 1A: Horizontal Main Line Model")
p_1b <- render_dag_panel(nodes_1b, edges_1b, "DAG 1B: Converging Model")

# ------------------------------------------------------------------------------
# 5. DAG 1D Panel: Moderation Model with "Moderated" Hypothesis Text
# ------------------------------------------------------------------------------
render_dag1d_panel <- function() {
  
  box_w <- 1.05
  box_h <- 0.55
  
  nodes_1d <- tribble(
    ~node,         ~x,    ~y,   ~label,                  ~fill_color,
    "PredNCC",     0.85,  0.8,  "Predators\n(NCC)",      color_pred,
    "AltPreyNCC",  2.15,  0.8,  "Alt Prey /\nSST (NCC)", color_prey,
    "CPUE",        1.50,  2.7,  "Salmon\nAbundance",     color_salmon,
    "SAR",         4.50,  2.7,  "Salmon\nSAR",           color_salmon,
    "PredAK",      3.85,  0.8,  "Predators\n(AK)",       color_pred,
    "AltPreyAK",   5.15,  0.8,  "Alt Prey /\nSST (AK)",  color_prey
  )
  
  tilde_nodes <- tribble(
    ~id,     ~x,    ~y,   ~to_node, ~hyp_code,
    "NCC_m", 1.50,  0.8,  "CPUE",   "H3a (Moderated)",
    "AK_m",  4.50,  0.8,  "SAR",    "H3b (Moderated)"
  )
  
  main_edge <- tribble(
    ~from,  ~to,  ~hyp_code,  ~x_from, ~y_from, ~x_to, ~y_to,
    "CPUE", "SAR", "H4a (+)", 1.50 + box_w/2, 2.7, 4.50 - box_w/2, 2.7
  )
  
  p <- ggplot() +
    # Main CPUE -> SAR Path Edge
    geom_segment(
      data = main_edge,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
      color = "gray25", linewidth = 1.1
    ) +
    geom_label(
      data = main_edge,
      aes(x = (x_from + x_to)/2, y = y_from + 0.22, label = hyp_code),
      fill = alpha("white", 0.85), color = "navy", fontface = "bold", size = 4.2,
      label.padding = unit(0.12, "lines"), label.size = NA
    ) +
    # Dashed lines leading from boxes into ~ junction
    geom_segment(aes(x = 0.85 + box_w/2, y = 0.8, xend = 1.50 - 0.22, yend = 0.8), color = "gray40", linetype = "dashed", linewidth = 0.8) +
    geom_segment(aes(x = 2.15 - box_w/2, y = 0.8, xend = 1.50 + 0.22, yend = 0.8), color = "gray40", linetype = "dashed", linewidth = 0.8) +
    geom_segment(aes(x = 3.85 + box_w/2, y = 0.8, xend = 4.50 - 0.22, yend = 0.8), color = "gray40", linetype = "dashed", linewidth = 0.8) +
    geom_segment(aes(x = 5.15 - box_w/2, y = 0.8, xend = 4.50 + 0.22, yend = 0.8), color = "gray40", linetype = "dashed", linewidth = 0.8) +
    
    # Moderated Upward Path Arrows (from ~ junction to response node)
    geom_segment(
      data = tilde_nodes,
      aes(x = x, y = y + 0.22, xend = x, yend = 2.7 - box_h/2),
      arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
      color = "navy", linewidth = 1.1
    ) +
    geom_label(
      data = tilde_nodes,
      aes(x = x + 0.70, y = (y + 0.22 + 2.7 - box_h/2)/2, label = hyp_code),
      fill = alpha("white", 0.85), color = "navy", fontface = "bold", size = 4.2,
      label.padding = unit(0.12, "lines"), label.size = NA
    ) +
    # Rectangular Structural Variable Nodes
    geom_tile(
      data = nodes_1d,
      aes(x = x, y = y, fill = fill_color),
      width = box_w, height = box_h, color = "black", linewidth = 1.0
    ) +
    geom_text(
      data = nodes_1d,
      aes(x = x, y = y, label = label),
      color = "black", fontface = "bold", size = 4.0, lineheight = 0.88
    ) +
    # Tilde Node Junction Circle
    geom_point(
      data = tilde_nodes,
      aes(x = x, y = y),
      shape = 21, fill = "white", color = "navy", size = 9, stroke = 1.3
    ) +
    geom_text(
      data = tilde_nodes,
      aes(x = x, y = y, label = "~"),
      color = "navy", fontface = "bold", size = 7, vjust = 0.35
    ) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.1, 6.3)) +
    scale_y_continuous(limits = c(0.2, 3.8)) +
    labs(title = "DAG 1D: Weighted Moderation Model [Predators ~ (Alt Prey / SST)]") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(8, 8, 8, 8),
      panel.background = element_rect(fill = "white", color = "gray80", linewidth = 0.8)
    )
}

p_1d <- render_dag1d_panel()

# ------------------------------------------------------------------------------
# 6. Combine All 3 DAGs into 1 Single Publication Figure & Print
# ------------------------------------------------------------------------------
combined_dag <- p_1a / p_1b / p_1d +
  plot_annotation(
    title = "Conceptual Directed Acyclic Graphs (DAGs) for Salmon Survival Hypotheses",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  )

# Output sized specifically for portrait mode in MS Word (7.5" x 10.5" @ 300 DPI)
fig_png_path <- file.path(fig_out_dir, "Figure_Conceptual_DAGs_Combined.png")
ggsave(
  filename = fig_png_path,
  plot = combined_dag,
  width = 7.5,
  height = 10.5,
  units = "in",
  dpi = 300
)

# Render Combined DAG Figure directly to screen/viewer
print(combined_dag)