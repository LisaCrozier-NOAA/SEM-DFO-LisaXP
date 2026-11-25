# ==============================================================================
# Script: 02_analysis/00_Fig_DAG_Tab_Hypo.R
# Purpose: Render revised DAG 1A, 1B, & 1C + Table 2 (Optimized for MS Word Portrait)
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
master_out  <- file.path(proj_dir, "Outputs")
fig_out_dir <- file.path(master_out, "figures")
tbl_out_dir <- file.path(master_out, "tables")

dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# Color Palette
color_salmon <- "#F8AFA6"  # Soft Pink
color_prey   <- "#A8E6CF"  # Soft Light Green
color_pred   <- "#D7C4B7"  # Soft Light Brown

# ------------------------------------------------------------------------------
# 2. Table 2: A Priori Hypotheses (Without Literature Column)
# ------------------------------------------------------------------------------
hypotheses_df <- tribble(
  ~Code, ~Category, ~Ecological_Mechanism, ~Target_Path, ~Expected_Sign, ~Evaluated_DAGs,
  
  # Category 1: Bottom-Up
  "H1a", "H1: Bottom-up", "Early-marine prey increases juvenile growth", "PreyNCC -> Growth", "+", "DAG1A, DAG1B, DAG1C",
  "H1b", "H1: Bottom-up", "Early-marine growth increases early-ocean survival", "Growth -> CPUE", "+", "DAG1A",
  "H1c", "H1: Bottom-up", "Early-marine growth increases smolt-to-adult survival directly", "Growth -> SAR", "+", "DAG1B, DAG1C",
  "H1e", "H1: Bottom-up", "Subarctic marine productivity increases later-marine survival", "PreyAK -> SAR", "+", "DAG1A, DAG1B, DAG1C",
  
  # Category 2: Top-Down
  "H2a", "H2: Top-down", "Early-ocean predators reduce early-ocean survival", "PredNCC -> CPUE", "-", "DAG1A, DAG1B, DAG1C",
  "H2c", "H2: Top-down", "Later-marine predators reduce later-marine survival", "PredAK -> SAR", "-", "DAG1A, DAG1B, DAG1C",
  
  # Category 3: Indirect & Moderated
  "H3a", "H3: Top-down Mediated", "Alternate prey buffer early-ocean predation", "AltPreyNCC -> PredNCC", "-", "DAG1C",
  "H3b", "H3: Top-down Mediated", "Alternate prey buffer later-ocean predation", "AltPreyAK -> PredAK", "-", "DAG1C",
  "H3c", "H3: Ocean Conditions", "Ocean temperatures influence predator presence/metabolism", "SST -> Pred", "+/-", "DAG1C",
  
  # Category 4: Regional Life-Stage
  "H4a", "H4: Regional Bottleneck", "Early-ocean survival regulates smolt-to-adult survival", "CPUE -> SAR", "+", "DAG1A, DAG1B, DAG1C",
  "H4b", "H4: Regional Bottleneck", "Later-ocean survival regulates smolt-to-adult survival", "CPUE + PredAK -> SAR", "- (Predators)", "DAG1A, DAG1B, DAG1C"
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
    Evaluated_DAGs = md("**Target DAGs**")
  ) %>%
  cols_align(align = "center", columns = c(Code, Expected_Sign, Evaluated_DAGs)) %>%
  cols_align(align = "left", columns = c(Ecological_Mechanism, Target_Path)) %>%
  cols_width(
    Code ~ px(60),
    Ecological_Mechanism ~ px(340),
    Target_Path ~ px(180),
    Expected_Sign ~ px(100),
    Evaluated_DAGs ~ px(160)
  ) %>%
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(13),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    table.width = px(840),
    data_row.padding = px(4)
  )

# Save Table 2 as PNG and PDF (tighter canvas: 860px width)
tbl_png_path <- file.path(tbl_out_dir, "obj0_table2_hypo_matrix.png")
tbl_pdf_path <- file.path(tbl_out_dir, "obj0_table2_hypo_matrix.pdf")
gtsave(table2_gt, tbl_png_path, vwidth = 860, vheight = 760)
gtsave(table2_gt, tbl_pdf_path, vwidth = 860, vheight = 760)

# Render Table to screen/viewer device
tbl_img <- readPNG(tbl_png_path)
grid.newpage()
grid.raster(tbl_img)
# ------------------------------------------------------------------------------
# 3. Standard Rectangular Node DAG Engine
# ------------------------------------------------------------------------------
render_dag_panel <- function(nodes_df, edges_df, panel_title) {
  
  box_w <- 0.95
  box_h <- 0.55
  
  # Ensure all override columns exist
  for (col in c("x_mid_override", "y_mid_override", 
                "x_start_override", "y_start_override", 
                "x_end_override", "y_end_override")) {
    if (!col %in% names(edges_df)) edges_df[[col]] <- NA_real_
  }
  
  nodes_df <- nodes_df %>%
    mutate(label_wrapped = str_wrap(label, width = 10))
  
  edge_coords <- edges_df %>%
    left_join(nodes_df %>% select(node, x_from = x, y_from = y), by = c("from" = "node")) %>%
    left_join(nodes_df %>% select(node, x_to = x, y_to = y), by = c("to" = "node")) %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      dist = sqrt(dx^2 + dy^2),
      
      # Default radial truncation at rectangle boundaries
      x_start_calc = x_from + (dx / dist) * (box_w / 1.8),
      y_start_calc = y_from + (dy / dist) * (box_h / 1.8),
      x_end_calc   = x_to - (dx / dist) * (box_w / 1.8),
      y_end_calc   = y_to - (dy / dist) * (box_h / 1.8),
      
      # Use explicit coordinate overrides if provided
      x_start = ifelse(!is.na(x_start_override), x_start_override, x_start_calc),
      y_start = ifelse(!is.na(y_start_override), y_start_override, y_start_calc),
      x_end   = ifelse(!is.na(x_end_override),   x_end_override,   x_end_calc),
      y_end   = ifelse(!is.na(y_end_override),   y_end_override,   y_end_calc),
      
      # Midpoint label positioning
      x_mid   = ifelse(!is.na(x_mid_override), x_mid_override, (x_start + x_end) / 2),
      y_mid   = ifelse(!is.na(y_mid_override), y_mid_override, (y_start + y_end) / 2),
      y_label = y_mid + ifelse(dy == 0, 0.22, 0.05 * sign(dy))
    )
  
  p <- ggplot() +
    # Edge Path Arrows
    geom_segment(
      data = edge_coords,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      arrow = arrow(length = unit(0.30, "cm"), type = "closed"),
      color = "gray25", linewidth = 1.1
    ) +
    # Path Hypothesis Labels
    geom_label(
      data = edge_coords,
      aes(x = x_mid, y = y_label, label = hyp_code),
      fill = alpha("white", 0.90), color = "navy", fontface = "bold", size = 2.9,
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
      color = "black", fontface = "bold", size = 3.0, lineheight = 0.88
    ) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(-0.15, 6.35), ylim = c(0.12, 3.75), clip = "off") +
    labs(title = panel_title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(6, 12, 6, 12),
      panel.background = element_rect(fill = "white", color = "gray80", linewidth = 0.8)
    )
  
  return(p)
}
# ------------------------------------------------------------------------------
# 4. DAG 1A & 1B Node and Edge Definitions
# ------------------------------------------------------------------------------
# DAG 1A: Horizontal Main Line Model
nodes_1a <- tribble(
  ~node,      ~x,   ~y,   ~label,               ~fill_color,
  "PreyNCC",  0.7,  3.2,  "Prey/Competitor (NCC)",         color_prey,
  "PredNCC",  0.7,  0.8,  "Predators (NCC)",    color_pred,
  "Growth",   2.1,  2.0,  "Salmon Growth",      color_salmon,
  "CPUE",     3.5,  2.0,  "Salmon Abundance",   color_salmon,
  "SAR",      4.9,  2.0,  "Salmon SAR",         color_salmon,
  "PreyAK",   5.7,  3.2,  "Prey/Competitor (AK)",          color_prey,
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

# DAG 1B: Converging Model
nodes_1b <- tribble(
  ~node,      ~x,   ~y,   ~label,               ~fill_color,
  "PreyNCC",  0.7,  3.2,  "Prey/Competitor (NCC)",         color_prey,
  "Growth",   2.2,  2.7,  "Salmon Growth",      color_salmon,
  "PredNCC",  0.7,  0.8,  "Predators (NCC)",    color_pred,
  "CPUE",     2.2,  1.3,  "Salmon Abundance",   color_salmon,
  "SAR",      4.3,  2.0,  "Salmon SAR",         color_salmon,
  "PreyAK",   5.7,  3.2,  "Prey/Competitor (AK)",          color_prey,
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
# 5. DAG 1C Panel: Converging Model with Environmental/Alt Prey Conditioning
# ------------------------------------------------------------------------------
# Same structure as DAG 1B, but Predators are functions of Alt Prey / SST
# ------------------------------------------------------------------------------
# 5. DAG 1C Panel: Shifted Bottom Row + Top-of-Box PredAK Arrow
# ------------------------------------------------------------------------------
nodes_1c <- tribble(
  ~node,         ~x,    ~y,   ~label,               ~fill_color,
  "PreyNCC",     0.65,  3.20, "Prey/Competitor (NCC)",    color_prey,
  "Growth",      2.25,  2.70, "Salmon Growth",      color_salmon,
  "AltPreyNCC",  0.55,  0.55, "Alt Prey / SST NCC", color_prey,   # Moved down to y = 0.55
  "PredNCC",     1.75,  0.55, "Predators (NCC)",    color_pred,   # Moved down to y = 0.55
  "CPUE",        2.85,  1.30, "Salmon Abundance",   color_salmon, # Now has clear margin above y = 0.55
  "SAR",         4.20,  2.00, "Salmon SAR",         color_salmon,
  "PreyAK",      5.75,  3.20, "Prey/Competitor (AK)",     color_prey,
  "AltPreyAK",   4.55,  0.55, "Alt Prey / SST AK",  color_prey,   # Moved down to y = 0.55
  "PredAK",      5.75,  0.55, "Predators (AK)",     color_pred    # Moved down to y = 0.55
)

# PredAK box top is y = 0.55 + 0.55/2 = 0.825; box center x = 5.75
edges_1c <- tribble(
  ~from,         ~to,        ~hyp_code,   ~x_mid_override, ~y_mid_override, ~x_start_override, ~y_start_override,
  "PreyNCC",     "Growth",   "H1a (+)",   1.40,            2.95,            NA_real_,          NA_real_,
  "AltPreyNCC",  "PredNCC",  "H3a/c",     1.15,            0.80,            NA_real_,          NA_real_,
  "PredNCC",     "CPUE",     "H2a (-)",   2.25,            0.95,            NA_real_,          NA_real_, # Clear of node edges
  "Growth",      "SAR",      "H1c (+)",   3.20,            2.48,            NA_real_,          NA_real_,
  "CPUE",        "SAR",      "H4a (+)",   3.55,            1.72,            NA_real_,          NA_real_,
  "PreyAK",      "SAR",      "H1e (+)",   5.00,            2.70,            NA_real_,          NA_real_,
  "AltPreyAK",   "PredAK",   "H3b/c",     5.15,            0.80,            NA_real_,          NA_real_,
  # Explicitly start at the top center of PredAK box (x = 5.75, y = 0.83)
  "PredAK",      "SAR",      "H2c (-)",   5.15,            1.45,            5.75,              0.83
)

p_1c <- render_dag_panel(nodes_1c, edges_1c, "DAG 1C: Converging Model with Alt Prey & SST Driven Predators")

# ------------------------------------------------------------------------------
# 6. Combine All 3 DAGs into 1 Single Publication Figure & Print
# ------------------------------------------------------------------------------
combined_dag <- p_1a / p_1b / p_1c +
  plot_annotation(
#    title = "Conceptual Directed Acyclic Graphs (DAGs) for Salmon Survival Hypotheses",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  )
combined_dag

# Output sized specifically for portrait mode in MS Word (7.5" x 10.5" @ 300 DPI)
fig_png_path <- file.path(fig_out_dir, "obj0_fig0_Conceptual_DAGs.png")
fig_pdf_path <- file.path(fig_out_dir, "obj0_fig0_Conceptual_DAGs.pdf")
ggsave(
  filename = fig_png_path,
  plot = combined_dag,
  width = 6.5,
  height = 8.5,
  units = "in",
  dpi = 300
)
ggsave(
  filename = fig_pdf_path,
  plot = combined_dag,
  width = 6.5,
  height = 8.5,
  units = "in",
  dpi = 300
)


# Render Combined DAG Figure directly to screen/viewer
print(combined_dag)
