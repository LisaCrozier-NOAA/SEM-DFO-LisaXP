library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS
# ------------------------------------------------------------------------------
figures_dir <- file.path("Rcode_for_paper", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PREY SPECIES & CARD COORDINATES
# ------------------------------------------------------------------------------
prey_df <- data.frame(
  group  = c("Market Squid",
             "Northern Anchovy",
             "Pacific Sardine",
             "Eulachon",
             "Pacific Herring (NCC & AK)",
             "Walleye Pollock / Hake",
             "Capelin",
             "EGoA & WCVI Krill (Euphausiids)",
             "Pink Salmon (North America)"),
  region = c("NCC (Guild 1 & 2)", "NCC (Guild 1 & 2)", "NCC (Guild 1 & 2)","NCC (Guild 1 & 2)", 
             "NCC & AK (Guilds 1-5)",  "AK / NCC (Guilds 2-5)", 
             "AK (Guilds 3-5)", "AK (Guild 3)", "AK (Guild 5)"),
  x = rep(c(1.8, 5.0, 8.2), length.out = 9),
  y = c(5.8, 5.8, 5.8, 3.8, 3.8, 3.8, 1.8, 1.8, 1.8),
  type = c("squid", "forage_fish", "forage_fish", "herring", "forage_fish", "gadid", "capelin", "krill", "salmon")
)

# ------------------------------------------------------------------------------
# STEP 2: BUILD REFINED ANATOMICAL VECTOR POLYGONS
# ------------------------------------------------------------------------------
get_prey_polygons <- function(df) {
  poly_list <- list()
  
  for (i in 1:nrow(df)) {
    cx <- df$x[i]
    cy <- df$y[i] + 0.05
    tp <- df$type[i]
    
    if (tp == "squid") {
      # Mantle / Head
      m_x <- cx + c(-0.1, -0.3, -0.1, 0.25, 0.35, 0.25) * 0.75
      m_y <- cy + c(0.18, 0.0, -0.18, -0.12, 0.0, 0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = m_x, y = m_y, id = paste0(i, "_mantle"), fill = "#C2410C")
      
      # Tail Fins (Rhomboid)
      f_x <- cx + c(-0.3, -0.45, -0.3, -0.2) * 0.75
      f_y <- cy + c(0.0, 0.18, 0.0, -0.18) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = f_x, y = f_y, id = paste0(i, "_fin"), fill = "#9A3412")
      
      # Tentacles (Trailing right)
      t_x <- cx + c(0.35, 0.55, 0.52, 0.35) * 0.75
      t_y <- cy + c(0.05, 0.12, 0.08, 0.02) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = t_x, y = t_y, id = paste0(i, "_t1"), fill = "#EA580C")
      
      t2_x <- cx + c(0.35, 0.58, 0.54, 0.35) * 0.75
      t2_y <- cy + c(-0.05, -0.14, -0.09, -0.02) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = t2_x, y = t2_y, id = paste0(i, "_t2"), fill = "#EA580C")
      
    } else if (tp == "krill") {
      # Curved Euphausiid Segmented Body
      b_x <- cx + c(-0.35, -0.2, 0.0, 0.2, 0.35, 0.25, 0.0, -0.2) * 0.75
      b_y <- cy + c(0.12, 0.18, 0.15, 0.08, -0.02, 0.02, 0.06, 0.05) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#DC2626")
      
      # Cephalothorax / Head
      h_x <- cx + c(-0.35, -0.48, -0.42, -0.3) * 0.75
      h_y <- cy + c(0.12, 0.18, 0.02, 0.02) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = h_x, y = h_y, id = paste0(i, "_head"), fill = "#B91C1C")
      
      # Feathery Legs (Thoracopods)
      l_x <- cx + c(-0.25, -0.28, -0.1, -0.13, 0.05, 0.02) * 0.75
      l_y <- cy + c(0.02, -0.2, 0.04, -0.18, 0.02, -0.15) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = l_x, y = l_y, id = paste0(i, "_legs"), fill = "#EF4444")
      
    } else if (tp %in% c("forage_fish", "herring", "capelin", "gadid", "salmon")) {
      # Sleek Streamlined Fish Body
      # Color palette by type
      col_main <- case_when(
        tp == "herring"     ~ "#0284C7",
        tp == "capelin"     ~ "#0D9488",
        tp == "gadid"       ~ "#475569",
        tp == "salmon"      ~ "#1E3A8A",
        TRUE                ~ "#2563EB"
      )
      col_dark <- case_when(
        tp == "herring"     ~ "#0369A1",
        tp == "capelin"     ~ "#0F766E",
        tp == "gadid"       ~ "#334155",
        tp == "salmon"      ~ "#1E293B",
        TRUE                ~ "#1D4ED8"
      )
      
      # Fusiform Body
      px <- cx + c(-0.45, -0.2, 0.1, 0.35, 0.45, 0.35, 0.1, -0.2) * 0.75
      py <- cy + c(0.0, 0.15, 0.12, 0.05, 0.0, -0.05, -0.12, -0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = px, y = py, id = paste0(i, "_body"), fill = col_main)
      
      # Dorsal Fin
      df_x <- cx + c(-0.05, 0.05, 0.12) * 0.75
      df_y <- cy + c(0.12, 0.28, 0.08) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = col_dark)
      
      # Forked Caudal Tail
      tf_x <- cx + c(-0.45, -0.6, -0.52, -0.6, -0.42) * 0.75
      tf_y <- cy + c(0.0, 0.22, 0.0, -0.22, 0.0) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = tf_x, y = tf_y, id = paste0(i, "_tail"), fill = col_dark)
      
      # Pectoral Fin
      pf_x <- cx + c(0.18, 0.28, 0.22) * 0.75
      pf_y <- cy + c(-0.05, -0.18, -0.02) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = pf_x, y = pf_y, id = paste0(i, "_pfin"), fill = col_dark)
    }
  }
  
  bind_rows(poly_list)
}

polys <- get_prey_polygons(prey_df)

# ------------------------------------------------------------------------------
# STEP 3: ASSEMBLE GGPLOT CANVAS & SAVE SLIDE
# ------------------------------------------------------------------------------
p <- ggplot() + 
  theme_void() + 
  xlim(0, 10) + ylim(0, 7.2) +
  
  # Slide Header
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 6.4, ymax = 7.0, fill = "#0F172A", color = NA) +
  annotate("text", x = 5, y = 6.7, label = "Candidate Alternate Prey Species in SEM Gating Models", 
           color = "white", fontface = "bold", size = 5)

# Render Card Boxes and Text Labels
for (i in 1:nrow(prey_df)) {
  cx <- prey_df$x[i]
  cy <- prey_df$y[i]
  
  p <- p + 
    annotate("rect", xmin = cx - 1.45, xmax = cx + 1.45, ymin = cy - 0.75, ymax = cy + 0.65, 
             fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.8) +
    annotate("rect", xmin = cx - 1.45, xmax = cx + 1.45, ymin = cy + 0.38, ymax = cy + 0.65, 
             fill = "#334155", color = NA) +
    annotate("text", x = cx, y = cy + 0.51, label = prey_df$region[i], 
             color = "white", fontface = "bold", size = 2.9) +
    annotate("text", x = cx, y = cy - 0.55, label = prey_df$group[i], 
             color = "#0F172A", fontface = "bold", size = 3.0)
}

# Add Vector Icon Polygons to Plot
p <- p + 
  geom_polygon(data = polys, aes(x = x, y = y, group = id, fill = fill), color = NA) +
  scale_fill_identity()

# Print in RStudio
print(p)

# Save High-Resolution Image to target figures directory
out_png <- file.path(figures_dir, "prey_icons_slide.png")
ggsave(out_png, plot = p, width = 11, height = 7.2, dpi = 300)

cat(sprintf("\nPrey icons slide generated successfully and saved to: %s\n", out_png))