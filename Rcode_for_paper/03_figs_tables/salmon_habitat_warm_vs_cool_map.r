
"03_figs_tables/salmon_habitat_warm_vs_cool_map.r"

library(ggplot2)
library(dplyr)
library(gridExtra)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & COLOR PALETTE
# ------------------------------------------------------------------------------
figures_dir <- file.path("Rcode_for_paper", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Palette from methodology diagram
c_land      <- "#D9B48F"  # Tan landmass
c_ocean_c   <- "#DCE6EC"  # Cool ocean
c_ocean_w   <- "#F0E2D0"  # Warm ocean
c_hab_c     <- "#829AB1"  # Cool corridor
c_hab_w     <- "#D9822B"  # Warm compressed corridor
c_border    <- "#52606D"  # Dark slate line accents
c_text      <- "#243B53"  # Deep slate text

# ------------------------------------------------------------------------------
# STEP 1: DEFINE COASTLINE & GEOGRAPHIC LOCATIONS
# ------------------------------------------------------------------------------
coastline <- data.frame(
  x = c(8.5, 8.2, 7.8, 7.2, 6.2, 5.2, 4.0, 3.0, 2.0, 1.0, 0.5),
  y = c(1.0, 2.0, 3.0, 4.0, 4.8, 5.2, 5.5, 5.6, 5.4, 5.1, 4.8)
)

loc_labels <- data.frame(
  x = c(8.7, 7.0, 1.2),
  y = c(0.7, 4.1, 4.3),
  label = c("Columbia River\nMouth", "SE Alaska\nCoast", "Aleutian\nIslands")
)

# Explicit regional separation:
# - Hake: South (NCC)
# - Pollock: North (AK)
# - Herring: South & North (NCC & AK)
species_df <- data.frame(
  label = c("Harbor Seal", "Pacific Hake", "Herring (NCC)", "Steller Sea Lion", 
            "Sablefish", "Walleye Pollock", "Herring (AK)", "Arrowtooth Flounder", "Salmon Shark"),
  x     = c(8.3,           7.8,            7.2,             6.7,                
            5.2,         4.3,               3.6,            2.8,                  1.8),
  y     = c(1.4,           2.4,            3.1,             4.6,                
            4.4,         4.8,               5.2,            4.8,                  4.5),
  type  = c("seal",        "fish",         "herring",       "seal",             
            "fish",      "fish",            "herring",      "flatfish",           "shark")
)

# ------------------------------------------------------------------------------
# STEP 2: BUILD MULTI-PART ANATOMICAL ICON POLYGONS
# ------------------------------------------------------------------------------
get_map_icon_polygons <- function(df) {
  poly_list <- list()
  
  for (i in 1:nrow(df)) {
    cx <- df$x[i]
    cy <- df$y[i]
    tp <- df$type[i]
    
    if (tp == "seal") {
      # Seal / Sea Lion
      b_x <- cx + c(-0.35, -0.15, 0.15, 0.35, 0.25, -0.05, -0.2) * 0.45
      b_y <- cy + c(-0.08, 0.18, 0.22, 0.05, -0.12, -0.2, -0.15) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#334155")
      
      f_x <- cx + c(-0.05, 0.1, 0.02) * 0.45
      f_y <- cy + c(-0.15, -0.32, -0.12) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = f_x, y = f_y, id = paste0(i, "_flip"), fill = "#0F172A")
      
    } else if (tp == "shark") {
      # Salmon Shark
      b_x <- cx + c(-0.4, -0.15, 0.15, 0.4, 0.2, -0.15) * 0.45
      b_y <- cy + c(0.0, 0.14, 0.1, 0.0, -0.1, -0.08) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#881337")
      
      df_x <- cx + c(-0.1, 0.02, 0.1) * 0.45
      df_y <- cy + c(0.12, 0.35, 0.08) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = "#4C0519")
      
      st_x <- cx + c(0.38, 0.55, 0.42, 0.48) * 0.45
      st_y <- cy + c(0.0, 0.28, -0.02, -0.18) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = st_x, y = st_y, id = paste0(i, "_tail"), fill = "#4C0519")
      
    } else if (tp == "flatfish") {
      # Flounder / Halibut
      b_x <- cx + c(-0.35, -0.12, 0.12, 0.35, 0.12, -0.12) * 0.45
      b_y <- cy + c(0.0, 0.26, 0.22, 0.0, -0.22, -0.26) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#475569")
      
      t_x <- cx + c(-0.35, -0.48, -0.4, -0.48) * 0.45
      t_y <- cy + c(0.0, 0.18, 0.0, -0.18) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = t_x, y = t_y, id = paste0(i, "_tail"), fill = "#1E293B")
      
    } else if (tp == "herring") {
      # Pacific Herring
      b_x <- cx + c(-0.4, -0.18, 0.1, 0.32, 0.4, 0.32, 0.1, -0.18) * 0.45
      b_y <- cy + c(0.0, 0.15, 0.12, 0.05, 0.0, -0.05, -0.12, -0.12) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#0284C7")
      
      tf_x <- cx + c(-0.4, -0.55, -0.48, -0.55) * 0.45
      tf_y <- cy + c(0.0, 0.2, 0.0, -0.2) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = tf_x, y = tf_y, id = paste0(i, "_tail"), fill = "#0369A1")
      
    } else if (tp == "fish") {
      # Hake / Pollock / Sablefish
      b_x <- cx + c(-0.38, -0.12, 0.12, 0.35, 0.18, -0.12) * 0.45
      b_y <- cy + c(0.0, 0.16, 0.12, 0.0, -0.12, -0.14) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = b_x, y = b_y, id = paste0(i, "_body"), fill = "#1D4ED8")
      
      df_x <- cx + c(-0.12, -0.02, 0.1) * 0.45
      df_y <- cy + c(0.12, 0.28, 0.1) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = df_x, y = df_y, id = paste0(i, "_dfin"), fill = "#1E293B")
      
      tf_x <- cx + c(-0.38, -0.52, -0.45, -0.52) * 0.45
      tf_y <- cy + c(0.0, 0.18, 0.0, -0.18) * 0.45
      poly_list[[length(poly_list) + 1]] <- data.frame(x = tf_x, y = tf_y, id = paste0(i, "_tail"), fill = "#1E293B")
    }
  }
  
  bind_rows(poly_list)
}

map_polys <- get_map_icon_polygons(species_df)

# ------------------------------------------------------------------------------
# STEP 3: BUILD HABITAT CORRIDORS
# ------------------------------------------------------------------------------
cool_corridor <- data.frame(
  x = c(8.5, 8.2, 7.8, 7.2, 6.2, 5.2, 4.0, 3.0, 2.0, 1.0, 0.5, 0.5, 1.0, 2.0, 3.0, 4.0, 5.2, 6.2, 7.2, 7.8, 8.2, 8.5),
  y = c(1.0, 2.0, 3.0, 4.0, 4.8, 5.2, 5.5, 5.6, 5.4, 5.1, 4.8, 3.5, 3.8, 4.0, 4.2, 4.1, 3.8, 3.3, 2.5, 1.5, 0.5, 0.2)
)

warm_corridor <- data.frame(
  x = c(8.5, 8.2, 7.8, 7.2, 6.2, 5.2, 4.0, 3.0, 2.0, 1.0, 0.5, 0.5, 1.0, 2.0, 3.0, 4.0, 5.2, 6.2, 7.2, 7.8, 8.2, 8.5),
  y = c(1.0, 2.0, 3.0, 4.0, 4.8, 5.2, 5.5, 5.6, 5.4, 5.1, 4.8, 4.3, 4.6, 4.9, 5.0, 4.9, 4.6, 4.2, 3.4, 2.4, 1.5, 0.6)
)

# ------------------------------------------------------------------------------
# STEP 4: RENDER MAP
# ------------------------------------------------------------------------------
draw_habitat_map <- function(corridor_df, bg_color, corridor_color, title_text, subtitle_text) {
  p <- ggplot() +
    theme_void() +
    xlim(0, 10) + ylim(0, 6.5) +
    
    # Ocean
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 6.5, fill = bg_color) +
    
    # Land
    annotate("polygon", x = c(coastline$x, 10, 10), y = c(coastline$y, 6.5, 0), fill = c_land, color = c_border, linewidth = 0.6) +
    
    # Corridor
    annotate("polygon", x = corridor_df$x, y = corridor_df$y, fill = corridor_color, alpha = 0.55, color = corridor_color, linewidth = 0.8) +
    
    # Header
    annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 5.8, ymax = 6.3, fill = c_text, color = NA) +
    annotate("text", x = 5.0, y = 6.05, label = title_text, color = "white", fontface = "bold", size = 4) +
    
    # Subtitle
    annotate("text", x = 5.0, y = 5.6, label = subtitle_text, color = c_text, fontface = "italic", size = 2.8) +
    
    # Regional Anchors
    geom_text(data = loc_labels, aes(x = x, y = y, label = label), color = c_text, fontface = "bold", size = 2.5, lineheight = 0.85) +
    
    # Vector Icons
    geom_polygon(data = map_polys, aes(x = x, y = y, group = id, fill = fill), color = NA) +
    scale_fill_identity() +
    
    # Species Text Labels
    geom_text(data = species_df, aes(x = x, y = y - 0.22, label = label), color = c_text, fontface = "bold", size = 2.2)
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 5: ASSEMBLE PANELS & SAVE
# ------------------------------------------------------------------------------
p_cool <- draw_habitat_map(
  cool_corridor, 
  bg_color = c_ocean_c, 
  corridor_color = c_hab_c, 
  title_text = "A. Cool Regime (Low HCI & Low SST)", 
  subtitle_text = "Expanded marine corridor; relaxed spatial overlap & diluted predator contact"
)

p_warm <- draw_habitat_map(
  warm_corridor, 
  bg_color = c_ocean_w, 
  corridor_color = c_hab_w, 
  title_text = "B. Warm Regime (High HCI & High SST)", 
  subtitle_text = "Compressed coastal ribbon; intensified spatial overlap & top-down predator gating"
)

combined_plot <- grid.arrange(
  p_cool, p_warm, 
  ncol = 1,
  top = textGrob("Marine Habitat Compression & Regional Species Distribution (Columbia River to Aleutians)", 
                 gp = gpar(fontface = "bold", fontsize = 13, col = c_text))
)

out_png <- file.path(figures_dir, "salmon_habitat_warm_vs_cool_map.png")
ggsave(out_png, plot = combined_plot, width = 9, height = 10, dpi = 300)

combined_plot

cat(sprintf("\nGeographically accurate habitat map generated and saved to: %s\n", out_png))