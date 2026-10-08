library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS
# ------------------------------------------------------------------------------
figures_dir <- file.path("Rcode_for_paper", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PREDATOR GROUPS
# ------------------------------------------------------------------------------
predators <- data.frame(
  group = c("Common Murre / Cormorant / Gulls / Loons",
            "Chilipepper Rockfish",
            "Adult Chinook Salmon (Cannibalism/Competition)",
            "Pacific Hake",
            "Harbor Seals, Fur Seals & Sea Lions",
            "Harbor Porpoise",
            "Arrowtooth Flounder & Halibut",
            "Pacific Cod & Sablefish",
            "Spiny Dogfish",
            "Sleeper & Salmon Sharks"),
  region = c("NCC (Avian)", "NCC (Fish)", "NCC (Salmonid)", "NCC (Gadid)", 
             "NCC / AK (Pinniped)", "NCC (Cetacean)", "AK (Flatfish)", 
             "AK (Gadid/Groundfish)", "AK (Elasmobranch)", "AK (Apex Shark)"),
  x = rep(c(1.8, 5.0, 8.2), length.out = 10),
  y = c(6.2, 6.2, 6.2, 4.5, 4.5, 4.5, 2.8, 2.8, 2.8, 1.1),
  shape_type = c("bird", "rockfish", "salmon", "hake", "seal", "porpoise", "flatfish", "cod", "dogfish", "shark")
)

# ------------------------------------------------------------------------------
# STEP 2: BUILD DETAILED MULTI-PART ANATOMICAL POLYGONS
# ------------------------------------------------------------------------------
get_predator_polygons <- function(df) {
  poly_list <- list()
  
  for (i in 1:nrow(df)) {
    cx <- df$x[i]
    cy <- df$y[i] + 0.05
    st <- df$shape_type[i]
    
    if (st == "bird") {
      # Cormorant / Seabird Body (Sleek diving bird)
      b_x <- cx + c(-0.3, -0.1, 0.15, 0.35, 0.2, -0.05) * 0.75
      b_y <- cy + c(0.05, 0.18, 0.12, 0.0, -0.15, -0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#1E293B")
      
      # Wings (Swept back)
      w_x <- cx + c(-0.05, -0.25, -0.05) * 0.75
      w_y <- cy + c(0.12, 0.32, 0.1) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = w_x, y = w_y, id = paste0(i, "_wing"), fill = "#0F172A")
      
      # Beak / Head
      bk_x <- cx + c(0.35, 0.52, 0.38) * 0.75
      bk_y <- cy + c(0.02, 0.0, -0.04) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = bk_x, y = bk_y, id = paste0(i, "_beak"), fill = "#D97706")
      
    } else if (st == "seal") {
      # Seal / Sea Lion Torso
      b_x <- cx + c(-0.4, -0.15, 0.15, 0.35, 0.25, -0.05, -0.25) * 0.75
      b_y <- cy + c(-0.08, 0.18, 0.22, 0.05, -0.12, -0.2, -0.15) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#334155")
      
      # Front Flipper
      f_x <- cx + c(-0.05, 0.1, 0.02) * 0.75
      f_y <- cy + c(-0.15, -0.32, -0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = f_x, y = f_y, id = paste0(i, "_flipper"), fill = "#1E293B")
      
      # Hind Tail Flippers
      tf_x <- cx + c(0.35, 0.52, 0.48) * 0.75
      tf_y <- cy + c(0.05, 0.12, -0.08) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = tf_x, y = tf_y, id = paste0(i, "_tflipper"), fill = "#1E293B")
      
    } else if (st == "porpoise") {
      # Harbor Porpoise Body
      b_x <- cx + c(-0.45, -0.2, 0.1, 0.38, 0.25, -0.1) * 0.75
      b_y <- cy + c(-0.02, 0.16, 0.18, 0.0, -0.14, -0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#0F172A")
      
      # Triangular Dorsal Fin
      df_x <- cx + c(-0.05, 0.05, 0.08) * 0.75
      df_y <- cy + c(0.17, 0.35, 0.15) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = "#0F172A")
      
      # Tail Fluke
      fl_x <- cx + c(0.38, 0.52, 0.45, 0.52) * 0.75
      fl_y <- cy + c(0.0, 0.15, 0.0, -0.15) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = fl_x, y = fl_y, id = paste0(i, "_fluke"), fill = "#0F172A")
      
    } else if (st == "flatfish") {
      # Flatfish Asymmetric Body (Flounder / Halibut)
      b_x <- cx + c(-0.4, -0.15, 0.15, 0.38, 0.15, -0.15) * 0.75
      b_y <- cy + c(0.0, 0.28, 0.24, 0.0, -0.24, -0.28) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#475569")
      
      # Tail Fin
      t_x <- cx + c(-0.4, -0.52, -0.45, -0.52) * 0.75
      t_y <- cy + c(0.0, 0.18, 0.0, -0.18) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = t_x, y = t_y, id = paste0(i, "_tail"), fill = "#334155")
      
    } else if (st %in% c("shark", "dogfish")) {
      # Shark / Dogfish Sleek Streamlined Body
      col_main <- if (st == "dogfish") "#64748B" else "#881337"
      col_dark <- if (st == "dogfish") "#475569" else "#4C0519"
      
      b_x <- cx + c(-0.45, -0.15, 0.15, 0.42, 0.2, -0.15) * 0.75
      b_y <- cy + c(0.0, 0.14, 0.1, 0.0, -0.1, -0.08) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = col_main)
      
      # Sharp Dorsal Fin
      df_x <- cx + c(-0.1, 0.02, 0.1) * 0.75
      df_y <- cy + c(0.12, 0.35, 0.08) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = col_dark)
      
      # Heterocercal Shark Tail
      st_x <- cx + c(0.38, 0.55, 0.42, 0.48) * 0.75
      st_y <- cy + c(0.0, 0.28, -0.02, -0.18) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = st_x, y = st_y, id = paste0(i, "_stail"), fill = col_dark)
      
    } else if (st %in% c("rockfish", "salmon", "hake", "cod")) {
      # Spiny Rockfish / Deep Cod / Salmon / Hake
      col_main <- case_when(
        st == "rockfish" ~ "#C2410C",
        st == "salmon"   ~ "#1E3A8A",
        st == "cod"      ~ "#3B82F6",
        TRUE             ~ "#475569"
      )
      col_dark <- case_when(
        st == "rockfish" ~ "#9A3412",
        st == "salmon"   ~ "#1E293B",
        st == "cod"      ~ "#1D4ED8",
        TRUE             ~ "#334155"
      )
      
      # Body
      b_x <- cx + c(-0.4, -0.15, 0.15, 0.38, 0.2, -0.15) * 0.75
      b_y <- cy + c(0.0, 0.18, 0.14, 0.0, -0.14, -0.15) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = col_main)
      
      # Spiny / Pronounced Dorsal Fin
      df_x <- cx + c(-0.15, -0.05, 0.05, 0.18) * 0.75
      df_y <- cy + c(0.15, 0.32, 0.28, 0.12) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = col_dark)
      
      # Tail
      tf_x <- cx + c(-0.4, -0.55, -0.48, -0.55) * 0.75
      tf_y <- cy + c(0.0, 0.2, 0.0, -0.2) * 0.75
      poly_list[[length(poly_list) + 1]] <- data.frame(x = tf_x, y = tf_y, id = paste0(i, "_tail"), fill = col_dark)
    }
  }
  
  bind_rows(poly_list)
}

polys <- get_predator_polygons(predators)

# ------------------------------------------------------------------------------
# STEP 3: BUILD GGPLOT CANVAS & SAVE SLIDE
# ------------------------------------------------------------------------------
p <- ggplot() + 
  theme_void() + 
  xlim(0, 10) + ylim(0, 7.5) +
  
  # Slide Header
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 6.8, ymax = 7.4, fill = "#0F172A", color = NA) +
  annotate("text", x = 5, y = 7.1, label = "Predator Taxonomic Groups in SEM Analysis", 
           color = "white", fontface = "bold", size = 5)

# Add Cards & Text
for (i in 1:nrow(predators)) {
  cx <- predators$x[i]
  cy <- predators$y[i]
  
  p <- p + 
    annotate("rect", xmin = cx - 1.45, xmax = cx + 1.45, ymin = cy - 0.65, ymax = cy + 0.65, 
             fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.8) +
    annotate("rect", xmin = cx - 1.45, xmax = cx + 1.45, ymin = cy + 0.38, ymax = cy + 0.65, 
             fill = "#334155", color = NA) +
    annotate("text", x = cx, y = cy + 0.51, label = predators$region[i], 
             color = "white", fontface = "bold", size = 3) +
    annotate("text", x = cx, y = cy - 0.45, label = predators$group[i], 
             color = "#0F172A", fontface = "bold", size = 2.7)
}

# Add Polygon Icons Directly to ggplot Object
p <- p + 
  geom_polygon(data = polys, aes(x = x, y = y, group = id, fill = fill), color = NA) +
  scale_fill_identity()

# Print to RStudio Screen
print(p)

# Save high-resolution PNG for PowerPoint
out_png <- file.path(figures_dir, "predator_icons_slide.png")
ggsave(out_png, plot = p, width = 11, height = 7.5, dpi = 300)

cat(sprintf("\nPredator icons slide updated successfully and saved to: %s\n", out_png))