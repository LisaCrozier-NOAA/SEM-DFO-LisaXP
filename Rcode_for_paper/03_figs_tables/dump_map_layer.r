# ==============================================================================
# Script: 04_visuals/02_oceanic_salmon_migration_layer.R
# Purpose: Add 8-Stage Salmon Migration Loop to Plain Base Map
# Output: map_layer_with_salmon.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic) # For Chinook salmon silhouette

# Set up visual directory
proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. FETCH BASE VECTOR LANDMASS
# ------------------------------------------------------------------------------
land <- ne_countries(scale = "medium", returnclass = "sf")

# Fetch Chinook Salmon PhyloPic UUID
# Uses standard Chinook salmon silhouette
chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha", img_index = 1)

# ------------------------------------------------------------------------------
# 2. DEFINE MIGRATION STAGES & TEARDROP GEOMETRY
# ------------------------------------------------------------------------------
# Stage details matching growth trajectory and bathymetric depth profile
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,     ~lat,    ~length_mm, ~depth_m, ~angle,
  1,      "CR Exit",           -124.2,   45.5,    100,        10,       75,
  2,      "WCVI Coastal",      -126.8,   49.2,    150,        30,       55,
  3,      "SEAK Shelf",        -136.0,   56.5,    200,        40,       30,
  4,      "GOA Top Shelf",     -146.0,   59.2,    300,        50,       -10,
  5,      "GOA Upper Offshore",-153.0,   56.0,    400,        150,      -110,
  6,      "GOA West Offshore", -147.0,   51.5,    450,        120,      -145,
  7,      "GOA Swing Back",    -138.0,   48.5,    475,        80,       -165,
  8,      "Return WCVI Offsh", -129.5,   46.5,    500,        30,       -175
) %>%
  mutate(
    # Scale icon sizes proportionally based on 100mm baseline
    icon_size = (length_mm / 100) * 0.045
  )

# Base Map Geographic Labels (from your updated script)
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

# Color Gradient Palette: Light Blue (Surface/10m) -> Deep Navy (150m Depth)
depth_palette <- colorRampPalette(c("#B2EBF2", "#4DD0E1", "#0288D1", "#01579B", "#0A192F"))

# Assign depth colors across the 8 stages
stage_colors <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# ------------------------------------------------------------------------------
# 3. BUILD ggplot MAP WITH SALMON OVERLAY
# ------------------------------------------------------------------------------
p <- ggplot() +
  # Base Landmass
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  
  # Oceanic Styling
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#1d5873", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  
  # Geographic Labels
  geom_text(
    data = labels_df,
    aes(x = x, y = y, label = label),
    color = "white",
    fontface = labels_df$fontface,
    size = labels_df$size,
    lineheight = 0.9,
    hjust = 0.5
  ) +
  
  # Teardrop Migration Path Line (Dashed)
  geom_path(
    data = salmon_stages,
    aes(x = lon, y = lat),
    color = "white",
    linetype = "dashed",
    linewidth = 0.6,
    alpha = 0.7
  )

# Add Individual Salmon Icons iteratively (Scaled, Colored, and Rotated along trajectory)
for (i in 1:nrow(salmon_stages)) {
  stg <- salmon_stages[i, ]
  p <- p + add_phylopic(
    uuid = chinook_uuid,
    x = stg$lon,
    y = stg$lat,
    ysize = stg$icon_size,
    color = stage_colors[i],
    angle = stg$angle
  )
}

# Apply Coordinate Bounding Box
p <- p + coord_sf(
  xlim = c(-170, -119),
  ylim = c(43, 62),
  expand = FALSE
)

# ------------------------------------------------------------------------------
# 4. RENDER AND SAVE
# ------------------------------------------------------------------------------
print(p)

ggsave(
  filename = file.path(fig_out_dir, "map_layer_with_salmon.png"),
  plot = p,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Salmon migration overlay saved to: ", file.path(fig_out_dir, "map_layer_with_salmon.png"))