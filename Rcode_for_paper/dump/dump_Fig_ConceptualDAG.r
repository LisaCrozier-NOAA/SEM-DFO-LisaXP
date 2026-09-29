# ==============================================================================
# Script: 05_render_conceptual_dags.R
# Purpose: Generate publication-ready conceptual DAG diagrams with H-label paths
# ==============================================================================

library(tidyverse)
library(lavaan)
library(tidySEM)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Output Setup
# ------------------------------------------------------------------------------
proj_dir   <- getwd()
output_dir <- file.path(proj_dir, "Routput_for_paper", "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Node display labels
node_pretty_labels <- c(
  "PreyNCC"   = "Prey\n(NCC)",
  "PredNCC"   = "Predators\n(NCC)",
  "Growth"    = "Salmon\nGrowth",
  "CPUE"      = "Salmon Abundance\n(CPUE)",
  "PreyAK"    = "Prey\n(AK)",
  "PredAK"    = "Predators\n(AK)",
  "SAR"       = "Salmon\nSAR"
)

# Shared 5x5 Layout Matrix
dag_layout_mat <- matrix(c(
  "PreyNCC", ""      , ""   , ""    , "PredNCC",
  ""       , "Growth", ""   , "CPUE", ""       ,
  ""       , ""      , "SAR", ""    , ""       ,
  "PreyAK" , ""      , ""   , ""    , "PredAK" ,
  ""       , ""      , ""   , ""    , ""
), nrow = 5, ncol = 5, byrow = TRUE)

# Helper to render a conceptual DAG graph using tidySEM
render_conceptual_dag <- function(model_syntax, layout_mat, edge_labels, dag_title, filename) {
  
  # 1. Dummy fit object to generate graph structure
  # Generate synthetic dataset matching variable names
  set.seed(123)
  vars <- c("PreyNCC", "PredNCC", "Growth", "CPUE", "PreyAK", "PredAK", "SAR")
  dummy_data <- as.data.frame(matrix(rnorm(100 * length(vars)), ncol = length(vars)))
  names(dummy_data) <- vars
  
  fit_dummy <- sem(model_syntax, data = dummy_data, fixed.x = FALSE)
  
  # 2. Prepare graph via tidySEM
  g <- prepare_graph(fit_dummy, layout = layout_mat, intercepts = FALSE, residuals = FALSE)
  
  # 3. Swap node labels for clean multi-line display names
  g$nodes$label <- node_pretty_labels[g$nodes$name]
  
  # 4. Map Hypothesis Codes (e.g. H1a, H2a) onto Structural Edges
  g$edges$label <- ""
  for (i in seq_len(nrow(g$edges))) {
    from_node <- g$edges$from[i]
    to_node   <- g$edges$to[i]
    key <- paste0(from_node, "->", to_node)
    
    if (key %in% names(edge_labels)) {
      g$edges$label[i] <- edge_labels[[key]]
    }
  }
  
  # Styling Options
  g$nodes$fill <- "aliceblue"
  g$nodes$color <- "skyblue4"
  g$nodes$label_color <- "black"
  g$edges$color <- "gray30"
  g$edges$label_color <- "darkblue"
  
  # Build ggplot object
  p <- plot(g) +
    labs(title = dag_title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      text = element_text(family = "sans")
    )
  
  # Save PNG/PDF
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot = p, width = 8, height = 6, bg = "white")
  
  return(p)
}

# ------------------------------------------------------------------------------
# 2. DAG 1A: Sequential Marine Stage Model
# ------------------------------------------------------------------------------
syntax_1a <- "
  Growth ~ PreyNCC
  CPUE   ~ Growth + PredNCC
  SAR    ~ CPUE + PreyAK + PredAK
"

labels_1a <- c(
  "PreyNCC->Growth" = "H1a (+)",
  "Growth->CPUE"    = "H1b (+)",
  "PredNCC->CPUE"   = "H2a (-)",
  "CPUE->SAR"       = "H1d (+)",
  "PreyAK->SAR"     = "H1e (+)",
  "PredAK->SAR"     = "H2c (-)"
)

p_1a <- render_conceptual_dag(
  syntax_1a, dag_layout_mat, labels_1a, 
  "DAG 1A: Sequential Marine Stage Model", "DAG1A_Conceptual"
)

# ------------------------------------------------------------------------------
# 3. DAG 1B: Direct Growth-to-SAR Link Model
# ------------------------------------------------------------------------------
syntax_1b <- "
  Growth ~ PreyNCC
  CPUE   ~ PredNCC
  SAR    ~ Growth + CPUE + PreyAK + PredAK
"

labels_1b <- c(
  "PreyNCC->Growth" = "H1a (+)",
  "PredNCC->CPUE"   = "H2a (-)",
  "Growth->SAR"     = "H1c (+)",
  "CPUE->SAR"       = "H1d (+)",
  "PreyAK->SAR"     = "H1e (+)",
  "PredAK->SAR"     = "H2c (-)"
)

p_1b <- render_conceptual_dag(
  syntax_1b, dag_layout_mat, labels_1b, 
  "DAG 1B: Direct Growth-to-SAR Link Model", "DAG1B_Conceptual"
)

# ------------------------------------------------------------------------------
# 4. DAG 1C: Direct Early Marine Predation & Saturated Path Model
# ------------------------------------------------------------------------------
syntax_1c <- "
  CPUE ~ PreyNCC + PredNCC
  SAR  ~ PreyNCC + Growth + CPUE + PredNCC + PreyAK + PredAK
"

labels_1c <- c(
  "PreyNCC->CPUE"   = "H1a* (+)",
  "PredNCC->CPUE"   = "H2a (-)",
  "PreyNCC->SAR"    = "H1a (+)",
  "Growth->SAR"     = "H1c (+)",
  "CPUE->SAR"       = "H1d (+)",
  "PredNCC->SAR"    = "H2b (-)",
  "PreyAK->SAR"     = "H1e (+)",
  "PredAK->SAR"     = "H2c (-)"
)

p_1c <- render_conceptual_dag(
  syntax_1c, dag_layout_mat, labels_1c, 
  "DAG 1C: Direct Early Marine Predation Model", "DAG1C_Conceptual"
)

cat("Conceptual DAG rendering complete. Diagrams saved to 'output/figures/'.\n")
