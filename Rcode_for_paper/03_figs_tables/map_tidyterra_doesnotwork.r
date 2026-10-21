# ==============================================================================
# Script: 04_visuals/01_oceanic_salmon_base_map.R
# Purpose: Generate Exact High-Quality PNW & GOA Base Map with Faint Relief
# Output: map_layer.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(maptiles)
library(tidyterra)
library(terra)

# Set up visual directory
proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Fetch geographic vector boundaries
land <- ne_countries(scale = "medium", returnclass = "sf")

# Define spatial bounding box in WGS84 (EPSG:4326)
bbox_sf <- st_bbox(c(xmin = -170, xmax = -119, ymin = 43, ymax = 62), crs = st_crs(4326)) %>%
  st_as_sfc()

# Download faint relief raster tiles (Esri World Ocean / World Hillshade)
# Note: "Esri.WorldOceanBase" or "Esri.WorldImagery" can also be used
bg_tiles <- get_tiles(bbox_sf, provider = "Esri.WorldPhysical", zoom = 5, crop = TRUE)

# Place labels (adjusted coordinate placement)
labels_df <- tribble(
  ~label,               ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   61.5,     "bold",      5.2,
  "British\nColumbia",  -123.5,   52.5,     "plain",     4.5,
  "WCVI",               -126.0,   49.5,     "plain",     4.2,
  "Columbia \nRiver",   -121.25,  45.5,     "plain",     4.5,
  "USA",                -120.75,  47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5
)

# 2. Build map with faint bathymetry & topographic relief
base_map <- ggplot() +
  # Layer 1: Render background terrain/bathymetry raster faintly (alpha = 0.45)
  geom_spatraster_rgb(data = bg_tiles, alpha = 0.45) +
  
  # Layer 2: Semi-transparent land overlay to blend colors smoothly
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", alpha = 0.5, linewidth = 0.4) +
  
  # Coordinate limits matching exact study bounds
  coord_sf(
    xlim = c(-170, -119), 
    ylim = c(43, 62), 
    expand = FALSE
  ) +
  
  # Styling & clean background theme
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#1d5873", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  
  # Layer 3: Place geographic labels on top of relief layer
  geom_text(
    data = labels_df,
    aes(x = x, y = y, label = label),
    color = "white",
    fontface = labels_df$fontface,
    size = labels_df$size,
    lineheight = 0.9,
    hjust = 0.5
  )

# Render & save base map layer
print(base_map)

ggsave(
  filename = file.path(fig_out_dir, "map_layer.png"),
  plot = base_map,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Relief map_layer.png saved successfully to: ", fig_out_dir)