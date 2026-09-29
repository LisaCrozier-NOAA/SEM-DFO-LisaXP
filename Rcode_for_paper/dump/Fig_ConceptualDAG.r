# ==============================================================================
# Script: 05_render_conceptual_dags_custom.R
# Purpose: Render revised DAG 1A, 1B, & 1D + Table 2 with Product Interaction Nodes
# Output: Rcode_for_paper/Routput_for_paper/
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(gt)

# ------------------------------------------------------------------------------
# 1. Directory & Color Palette Setup
# ------------------------------------------------------------------------------
out_dir     <- "Rcode_for_paper/Routput_for_paper"
fig_out_dir <- file.path(out_dir, "figures")
tbl_out_dir <- file.path(out_dir, "tables")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# Color Palette
color_salmon <- "#F8AFA6"  # Soft Pink
color_prey   <- "#A8E6CF"  # Soft Light Green
color_pred   <- "#D7C4B7"  # Soft Light Brown

# ------------------------------------------------------------------------------
# 2. Table 2: A Priori Hypotheses (Including Regional Category H4)
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
  
  # Category 4: Regional Life-Stage
  "H4a", "H4: Regional Bottleneck", "Early ocean residence (NCC) sets survival trajectory prior to AK stage", "CPUE -> SAR", "+ (Rel. Magnitude)", "DAG1A, DAG1B, DAG1D", "Wells et al. (2020)"
)

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

gtsave(table2_gt, file.path(tbl_out_dir, "Table2_Hypotheses_Matrix.html"))

# ------------------------------------------------------------------------------
# 3. Standard DAG Plotting Engine (DAG 1A & 1B)
# ------------------------------------------------------------------------------
render_dag_plot <- function(nodes_df, edges_df, dag_title, filename) {
  
  edge_coords <- edges_df %>%
    left_join(nodes_df %>% select(node, x_from = x, y_from = y), by = c("from" = "node")) %>%
    left_join(nodes_df %>% select(node, x_to = x, y_to = y), by = c("to" = "node")) %>%
    mutate(
      x_mid = (x_from + x_to) / 2,
      y_mid = (y_from + y_to) / 2
    )
  
  p <- ggplot() +
    # Draw Straight Path Arrows
    geom_segment(
      data = edge_coords,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      color = "gray30", linewidth = 0.85
    ) +
    # Draw Path Hypothesis Labels
    geom_label(
      data = edge_coords,
      aes(x = x_mid, y = y_mid, label = hyp_code),
      fill = "lightyellow", color = "darkblue", fontface = "bold", size = 3.2,
      label.padding = unit(0.18, "lines")
    ) +
    # Draw Nodes
    geom_point(
      data = nodes_df, 
      aes(x = x, y = y, fill = fill_color), 
      shape = 21, color = "gray20", size = 23, stroke = 1.1
    ) +
    # Draw Node Labels
    geom_text(
      data = nodes_df, 
      aes(x = x, y = y, label = label), 
      color = "black", fontface = "bold", size = 2.5
    ) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.2, 5.8)) +
    scale_y_continuous(limits = c(0.2, 3.8)) +
    labs(title = dag_title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 1),
      panel.background = element_rect(fill = "white", color = "gray85")
    )
  
  ggsave(file.path(fig_out_dir, paste0(filename, ".png")), plot = p, width = 9, height = 6, dpi = 300)
  ggsave(file.path(fig_out_dir, paste0(filename, ".pdf")), plot = p, width = 9, height = 6)
  
  return(p)
}

# ------------------------------------------------------------------------------
# 4. DAG 1A: Horizontal Main Line Layout
# ------------------------------------------------------------------------------
nodes_1a <- tribble(
  ~node,         ~x,   ~y,   ~label,                      ~fill_color,
  "PreyNCC",     1.0,  3.1,  "Prey\n(NCC)",               color_prey,
  "PredNCC",     1.0,  0.9,  "Predators\n(NCC)",          color_pred,
  "Growth",      2.2,  2.0,  "Salmon\nGrowth",            color_salmon,
  "CPUE",        3.5,  2.0,  "Salmon\nAbundance",         color_salmon,
  "SAR",         4.8,  2.0,  "Salmon\nSAR",               color_salmon,
  "PreyAK",      5.5,  3.1,  "Prey\n(AK)",                color_prey,
  "PredAK",      5.5,  0.9,  "Predators\n(AK)",           color_pred
)

edges_1a <- tribble(
  ~from,        ~to,      ~hyp_code,
  "PreyNCC",    "Growth", "H1a \n(+)",
  "PredNCC",    "CPUE",   "H2a \n(-)",
  "Growth",     "CPUE",   "H1b \n(+)",
  "CPUE",       "SAR",    "H4a \n(Regional)",
  "PreyAK",     "SAR",    "H1e \n(+)",
  "PredAK",     "SAR",    "H2c \n(-)"
)

p_1a <- render_dag_plot(nodes_1a, edges_1a, "DAG 1A: Horizontal Main Line Model", "DAG1A_Conceptual")

# ------------------------------------------------------------------------------
# 5. DAG 1B: Converging Layout
# ------------------------------------------------------------------------------
nodes_1b <- tribble(
  ~node,         ~x,   ~y,   ~label,                      ~fill_color,
  "PreyNCC",     0.9,  3.3,  "Prey\n(NCC)",               color_prey,
  "Growth",      2.2,  2.8,  "Salmon\nGrowth",            color_salmon,
  "PredNCC",     0.9,  0.7,  "Predators\n(NCC)",          color_pred,
  "CPUE",        2.2,  1.2,  "Salmon\nAbundance",         color_salmon,
  "SAR",         4.0,  2.0,  "Salmon\nSAR",               color_salmon,
  "PreyAK",      5.2,  3.1,  "Prey\n(AK)",                color_prey,
  "PredAK",      5.2,  0.9,  "Predators\n(AK)",           color_pred
)

edges_1b <- tribble(
  ~from,        ~to,      ~hyp_code,
  "PreyNCC",    "Growth", "H1a\n (+)",
  "PredNCC",    "CPUE",   "H2a \n(-)",
  "Growth",     "SAR",    "H1c\n (+)",
  "CPUE",       "SAR",    "H4a\n (Regional)",
  "PreyAK",     "SAR",    "H1e\n (+)",
  "PredAK",     "SAR",    "H2c\n (-)"
)

p_1b <- render_dag_plot(nodes_1b, edges_1b, "DAG 1B: Converging ( > ) Model Layout", "DAG1B_Conceptual")

# ------------------------------------------------------------------------------
# 6. DAG 1D Engine: Interaction Model with Touching Circles & Product Node
# ------------------------------------------------------------------------------
render_dag1d_interaction_plot <- function() {
  
  # Node Definitions for DAG 1D
  # Circles overlap slightly: delta_x = 0.52 (size 23 ~ 0.55 width units)
  nodes_1d <- tribble(
    ~node,         ~x,     ~y,   ~label,                      ~fill_color,
    "AltPreyNCC",  1.24,   0.9,  "Alt Prey\n(NCC)",           color_prey,
    "PredNCC",     1.76,   0.9,  "Predators\n(NCC)",          color_pred,
    "CPUE",        1.50,   2.7,  "Salmon\nAbundance",         color_salmon,
    "SAR",         3.80,   2.7,  "Salmon\nSAR",               color_salmon,
    "AltPreyAK",   3.54,   0.9,  "Alt Prey\n(AK)",            color_prey,
    "PredAK",      4.06,   0.9,  "Predators\n(AK)",           color_pred
  )
  
  # Product Junction Coordinates (Point where circles kiss/intersect)
  mult_junctions <- tribble(
    ~prod_id, ~x,     ~y,   ~to_node, ~hyp_code,
    "NCC_x",  1.50,   0.9,  "CPUE",   "H3a \n(xAltPrey)",
    "AK_x",   3.80,   0.9,  "SAR",    "H3b \n(xAltPrey)"
  )
  
  # Standard Path Edge (CPUE -> SAR)
  main_edge <- tribble(
    ~from,   ~to,   ~hyp_code, ~x_from, ~y_from, ~x_to, ~y_to,
    "CPUE",  "SAR", "H4a \n(Regional)", 1.50, 2.7, 3.80, 2.7
  )
  
  p <- ggplot() +
    # 1. Draw Main Line Edge (CPUE -> SAR)
    geom_segment(
      data = main_edge,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      color = "gray30", linewidth = 0.85
    ) +
    geom_label(
      data = main_edge,
      aes(x = (x_from + x_to)/2, y = (y_from + y_to)/2, label = hyp_code),
      fill = "lightyellow", color = "darkblue", fontface = "bold", size = 3.2,
      label.padding = unit(0.18, "lines")
    ) +
    # 2. Draw Interaction Product Arrows (from 'x' junction to Response)
    geom_segment(
      data = mult_junctions,
      aes(x = x, y = y + 0.15, xend = x, yend = 2.45),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      color = "darkblue", linewidth = 0.95
    ) +
    geom_label(
      data = mult_junctions,
      aes(x = x, y = (y + 0.15 + 2.45)/2, label = hyp_code),
      fill = "lightyellow", color = "darkblue", fontface = "bold", size = 3.2,
      label.padding = unit(0.18, "lines")
    ) +
    # 3. Draw Kissing Circles for Trophic Pairs
    geom_point(
      data = nodes_1d, 
      aes(x = x, y = y, fill = fill_color), 
      shape = 21, color = "gray20", size = 23, stroke = 1.1
    ) +
    geom_text(
      data = nodes_1d, 
      aes(x = x, y = y, label = label), 
      color = "black", fontface = "bold", size = 2.5
    ) +
    # 4. Draw Small Product Multiplication Node ('x') at Boundary Intersection
    geom_point(
      data = mult_junctions,
      aes(x = x, y = y),
      shape = 21, fill = "white", color = "darkblue", size = 6, stroke = 1.2
    ) +
    geom_text(
      data = mult_junctions,
      aes(x = x, y = y, label = "×"),
      color = "darkblue", fontface = "bold", size = 4
    ) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.2, 5.2)) +
    scale_y_continuous(limits = c(0.2, 3.5)) +
    labs(title = "DAG 1D: Two-Stage Interaction Model (Top-Down Moderated)") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 1),
      panel.background = element_rect(fill = "white", color = "gray85")
    )
  
  filename <- "DAG1D_Conceptual"
  ggsave(file.path(fig_out_dir, paste0(filename, ".png")), plot = p, width = 9, height = 6, dpi = 300)
  ggsave(file.path(fig_out_dir, paste0(filename, ".pdf")), plot = p, width = 9, height = 6)
  
  return(p)
}

p_1d <- render_dag1d_interaction_plot()

message("Table 2 and revised DAG figures (1A, 1B, 1D) saved to: ", out_dir)