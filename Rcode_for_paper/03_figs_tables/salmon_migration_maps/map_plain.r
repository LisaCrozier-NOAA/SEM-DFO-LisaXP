# ==============================================================================
# Script: 04_visuals/01_oceanic_salmon_base_map.R
# Purpose: Generate Exact High-Quality PNW & GOA Base Map (Slide Layer)
# Output: map_layer.png (Columbia River Labels Removed)
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)

# Set up visual directory
proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Fetch geographic layers
land <- ne_countries(scale = "medium", returnclass = "sf")

# Define coordinate boundaries precisely mapping PNW to Aleutian Islands
map_bounds <- coord_sf(
  xlim = c(-170, -119), 
  ylim = c(43, 62), 
  expand = FALSE
)

# Place labels inspired by original map (omitting CR and Columbia River)
labels_df <- tribble(
  ~label,              ~x,      ~y,      ~fontface,   ~size,
  "Aleutian \nIslands",  -167.0,  55.5,    "italic",    4.5,
#  "GOA",               -143.5,  58.5,    "bold",      5.5,
  "Alaska",            -143.0,  61.5,    "bold",      5.2,
#  "SEAK",              -132.5,  56.0,    "plain",     4.5,
  "British\nColumbia", -123.5,  52.5,    "plain",     4.5,
  "WCVI",              -126.0,  49.5,    "plain",     4.2,
  "Columbia \nRiver",        -121.25,  45.5,    "plain",     4.5,
"USA",        -120.75,  47.5,    "plain",     4.5,
"Pacific Ocean",     -158.0,  45.5,    "plain",     5.5
)

# 2. Build the exact base map without CR text
base_map <- ggplot() +
  # Add landmass with smooth earthy tones
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  
  # Custom ocean styling matching the deep teal/blue transitions
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#1d5873", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  
  # Add custom geographic labels
  geom_text(
    data = labels_df,
    aes(x = x, y = y, label = label),
    color = "white",
    fontface = labels_df$fontface,
    size = labels_df$size,
    lineheight = 0.9,
    hjust = 0.5
  ) +
  
  # Inject boundaries and projection
  map_bounds

# Render & Print Map Layer to view
print(base_map)

# Save high-resolution mapping layer for slides (Slide 1 Layer)
ggsave(
  filename = file.path(fig_out_dir, "map_layer.png"),
  plot = base_map,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Precise map_layer.png saved successfully to: ", fig_out_dir)