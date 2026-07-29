library(Cairo)
library(cowplot)
library(ggplot2)

CairoPDF("LisaXP/outputs_4/DAG1_TopModels_noHCI_Plots.pdf", width = 20, height = 14)

# Define a clean line style to reuse
divider_style <- list(color = "black", size = 1.5, linetype = "solid")

# ==========================================
# PAGE 1: LONG MODELS
# ==========================================

row1_long <- plot_grid(
  DAG1A_long_Dougtopmodel_plot, DAG1A_long_Dougtopmodel_reduced_plot, 
  labels = c("DAG1A Long: Full", "DAG1A Long: Reduced"), label_size = 14, hjust = -0.1
)

row2_long <- plot_grid(
  DAG1B_long_Dougtopmodel_plot, DAG1B_long_Dougtopmodel_reduced_plot, 
  labels = c("DAG1B Long: Full", "DAG1B Long: Reduced"), label_size = 14, hjust = -0.1
)

page1_grid <- plot_grid(row1_long, row2_long, ncol = 1)

# Turn the grid into a canvas and draw crosshairs directly in the middle (0.5)
ggdraw(page1_grid) +
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey70", size = 1) + # Horizontal Divider
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "grey70", size = 1)   # Vertical Divider


# ==========================================
# PAGE 2: SHORT MODELS
# ==========================================

row1_short <- plot_grid(
  DAG1A_short_Dougtopmodel_plot, DAG1A_short_Dougtopmodel_reduced_plot, 
  labels = c("DAG1A Short: Full", "DAG1A Short: Reduced"), label_size = 14, hjust = -0.1
)

row2_short <- plot_grid(
  DAG1B_short_Dougtopmodel_plot, DAG1B_short_Dougtopmodel_reduced_plot, 
  labels = c("DAG1B Short: Full", "DAG1B Short: Reduced"), label_size = 14, hjust = -0.1
)

page2_grid <- plot_grid(row1_short, row2_short, ncol = 1)

# Apply the same crosshair dividers to page 2
ggdraw(page2_grid) +
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey70", size = 1) + # Horizontal Divider
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "grey70", size = 1)   # Vertical Divider

# ==============================================================================
# PAGE 3: DAG1C MODELS (LONG & SHORT)
# ==============================================================================

row1_c <- plot_grid(
  DAG1C_long_Dougtopmodel_plot, DAG1C_long_Dougtopmodel_reduced_plot, 
  labels = c("DAG1C Long: Full", "DAG1C Long: Reduced"), label_size = 14, hjust = -0.1
)

row2_c <- plot_grid(
  DAG1C_short_noHCI_Dougtopmodel_plot, DAG1C_short_noHCI_Dougtopmodel_reduced_plot, 
  labels = c("DAG1C Short: Full", "DAG1C Short: Reduced"), label_size = 14, hjust = -0.1
)

page3_grid <- plot_grid(row1_c, row2_c, ncol = 1)

# Render Page 3 Canvas + Crosshair Dividers
print(
  ggdraw(page3_grid) +
    draw_line(x = c(0, 1), y = c(0.5, 0.5), color = divider_style$color, size = divider_style$size) + 
    draw_line(x = c(0.5, 0.5), y = c(0, 1), color = divider_style$color, size = divider_style$size)
)

dev.off()