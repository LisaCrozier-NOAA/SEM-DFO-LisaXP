# ==============================================================================
# Script: plot_peff_formula_presentation.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Generate an updated high-resolution slide visual explaining the Effective
#          Predator Load (Peff) mathematical formulation for Goal 3.
# Output: Rcode_for_paper/Routput_for_paper/plots/Peff_Formula_Slide.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# --- Create Presentation Canvas ---
png("Rcode_for_paper/Routput_for_paper/plots/Peff_Formula_Slide.png", 
    width = 2400, height = 1350, res = 300) # 16:9 Widescreen Ratio

grid.newpage()

# Slide Title
grid.text(
  "Goal 3: Effective Predator Load (Peff) Formulation",
  x = 0.05, y = 0.90, just = "left",
  gp = gpar(fontsize = 18, fontface = "bold", col = "#1A2530")
)

# Subtitle
grid.text(
  "Accounting for Spatial/Thermal Gating (SST) and Alternative Prey Buffering (Forage Pool)",
  x = 0.05, y = 0.84, just = "left",
  gp = gpar(fontsize = 12, fontface = "italic", col = "#555555")
)

# Main Formula Box
grid.rect(
  x = 0.5, y = 0.68, width = 0.9, height = 0.18,
  gp = gpar(fill = "#F4F6F8", col = "#2C3E50", lwd = 2)
)

# Main Formula Text
grid.text(
  expression(P[eff] == P[norm] %*% SST[scaled] %*% (1 - w[prey] %*% F[pool])),
  x = 0.5, y = 0.68,
  gp = gpar(fontsize = 20, fontface = "bold", col = "#2C3E50")
)

# Variable Breakdown Header
grid.text(
  "Model Components & Parameter Specifications:",
  x = 0.08, y = 0.52, just = "left",
  gp = gpar(fontsize = 13, fontface = "bold", col = "#1A2530")
)

# Component 1: P_norm
grid.text(
  expression(bold(P[norm]) * " : Min-Max Normalized Predator Abundance " * group("[", list(0, 1), "]")),
  x = 0.10, y = 0.44, just = "left",
  gp = gpar(fontsize = 11, col = "#E69F00")
)

# Component 2: S_scaled
grid.text(
  expression(bold(SST[scaled]) * " : Thermal/Spatial Overlap Gate " * group("(", list(w[SST] == 1.0), ")") * " scaled to " * group("[", list(a, 1), "]") * " where " * a == 0.2),
  x = 0.10, y = 0.36, just = "left",
  gp = gpar(fontsize = 11, col = "#56B4E9")
)

# Component 3: F_pool
grid.text(
  expression(bold(F[pool]) * " : Alternative Prey Buffering Pool " * group("(", list(sum(Prey[i], i==1, N)), ")") * " scaled to " * group("[", list(0, 1), "]")),
  x = 0.10, y = 0.28, just = "left",
  gp = gpar(fontsize = 11, col = "#009E73")
)

# Component 4: w_prey
grid.text(
  expression(bold(w[prey]) * " : Forage Pool Attenuation Weight Grid Search Parameter " * group("[", list(0.1, 1.0), "]") * " (saved as " * bold(abs_weight_w) * ")"),
  x = 0.10, y = 0.20, just = "left",
  gp = gpar(fontsize = 11, col = "#D55E00")
)

# Bottom Callout / Annotation
grid.text(
  "Note: Life stage suffixes (_smoltyr vs _adultyr) are strictly matched across all terms for a given model.",
  x = 0.5, y = 0.08,
  gp = gpar(fontsize = 9, fontface = "italic", col = "#777777")
)

dev.off()

message("Updated presentation slide graphic exported to Routput_for_paper/plots/Peff_Formula_Slide.png")