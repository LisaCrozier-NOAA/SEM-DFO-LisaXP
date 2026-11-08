
"03_figs_tables/Weighting_slide_ppt.r"

library(ggplot2)
library(grid)

# Create a clean canvas
p <- ggplot() + 
  theme_void() + 
  xlim(0, 10) + ylim(0, 6) +
  
  # --- Equation Header Box ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 4.8, ymax = 5.8, fill = "#1E293B", color = NA) +
  annotate("text", x = 5, y = 5.3, label = "Min-Max Joint Weighting Formula for Predators in SEM", 
           color = "white", fontface = "bold", size = 5) +
  
  # --- Component Boxes ---
  # Raw Predator
  annotate("rect", xmin = 0.5, xmax = 2.8, ymin = 3.2, ymax = 4.2, fill = "#EFF6FF", color = "#2563EB", linewidth = 1) +
  annotate("text", x = 1.65, y = 3.8, label = "Raw Predator (P)", fontface = "bold", color = "#1E40AF", size = 4) +
  annotate("text", x = 1.65, y = 3.4, label = "Scaled [0, 1]  \n(x - min) / (max - min)", color = "#475569", size = 3) +
  
  # Multiply Symbol 1
  annotate("text", x = 3.25, y = 3.7, label = "×", fontface = "bold", size = 7, color = "#0F172A") +
  
  # Spatial Overlap Gate
  annotate("rect", xmin = 3.7, xmax = 6.3, ymin = 3.2, ymax = 4.2, fill = "#FEF2F2", color = "#DC2626", linewidth = 1) +
  annotate("text", x = 5.0, y = 3.8, label = "(1 - w1 × HCI/SST)", fontface = "bold", color = "#991B1B", size = 4) +
  annotate("text", x = 5.0, y = 3.4, label = "Habitat Compression", color = "#475569", size = 3) +
  
  # Multiply Symbol 2
  annotate("text", x = 6.75, y = 3.7, label = "×", fontface = "bold", size = 7, color = "#0F172A") +
  
  # Alternate Prey Gate
  annotate("rect", xmin = 7.2, xmax = 9.5, ymin = 3.2, ymax = 4.2, fill = "#ECFDF5", color = "#059669", linewidth = 1) +
  annotate("text", x = 8.35, y = 3.8, label = "(1 - w2 × Alt Prey)", fontface = "bold", color = "#065F46", size = 4) +
  annotate("text", x = 8.35, y = 3.4, label = "Forage fish / Squid / Hake-pollock", color = "#475569", size = 3) +
  
  # Downward Arrows
  annotate("segment", x = 1.65, xend = 5.0, y = 3.0, yend = 2.2, arrow = arrow(length = unit(0.2, "cm")), color = "#475569", linewidth = 1) +
  annotate("segment", x = 5.0, xend = 5.0, y = 3.0, yend = 2.2, arrow = arrow(length = unit(0.2, "cm")), color = "#475569", linewidth = 1) +
  annotate("segment", x = 8.35, xend = 5.0, y = 3.0, yend = 2.2, arrow = arrow(length = unit(0.2, "cm")), color = "#475569", linewidth = 1) +
  
  # Result Box
  annotate("rect", xmin = 2.5, xmax = 7.5, ymin = 1.2, ymax = 2.1, fill = "#F8FAFC", color = "#0F172A", linewidth = 1.5) +
  annotate("text", x = 5.0, y = 1.7, label = "P_eff (Joint Effective Predator Intensity)", fontface = "bold", color = "#0F172A", size = 4.5) +
  annotate("text", x = 5.0, y = 1.4, label = "SEM: SAR ~ CPUE + P_eff  ==>  Yields negative path (β < 0)", color = "#1E293B", size = 3.2) +
  
  # Summary Banner at Bottom
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 0.1, ymax = 0.8, fill = "#F1F5F9", color = NA) +
  annotate("text", x = 5.0, y = 0.5, label = "Key Mechanism: Abundant preferred prey or cool habitat drops the multiplier toward 0 (w->1), reducing predator pressure on salmon.\nWhen buffering is low (w->0), effective predation equals raw predator intensity (Peff = 1.0 * Praw).", color = "#334155", size = 3, fontface = "italic")

# Print plot in RStudio viewer
print(p)

# Save high-res PNG for PowerPoint
ggsave(file.path("Rcode_for_paper","figures","weighting_function_slide.png"), plot = p, width = 10, height = 6, dpi = 300)
       