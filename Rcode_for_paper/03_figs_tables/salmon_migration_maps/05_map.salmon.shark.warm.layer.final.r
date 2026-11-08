# ==============================================================================
# Script: 04_visuals/03_oceanic_salmon_and_sharks_layer_v2.R
# Purpose: Version 2 - Projection-Corrected Shark Sizes Along AK Coastal Arc
# Output: map_layer_with_salmon_and_sharks_warm.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass & PhyloPic setup
land <- ne_countries(scale = "medium", returnclass = "sf")

chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
shark_uuid   <- get_uuid(name = "sharks")

# 3. Salmon Migration Stages (Offshore & Deeper AK Juveniles)
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,     ~lat,    ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",           -125.5,   46.8,    100,        10,       40,     FALSE,
  2,      "WCVI Coastal",      -128.2,   49.5,    150,        30,       30,     FALSE,
  3,      "BC Coast",          -132.2,   52.2,    150,        50,       30,     FALSE,
  4,      "SEAK Shelf",        -137.8,   55.8,    200,        80,       15,     FALSE,
  5,      "GOA Top Shelf",     -147.5,   58.2,    120,       100,        0,     FALSE,
  6,      "GOA West Offshore", -153.5,   55.5,    310,       120,     -35,     FALSE,
  7,      "GOA Central",       -147.0,   51.5,    360,       150,      20,     TRUE,
  8,      "GOA Swing Back",    -137.5,   48.2,    410,       180,      10,     TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Version 2 Sharks: 3 Light Grey Shallow Sharks (Scaled for Latitude)
# Reduced height for northernmost shark (id 2) to compensate for latitudinal distortion
shark_stages_v2 <- tribble(
  ~shark_id, ~location,           ~lon,     ~lat,    ~height, ~color,    ~angle, ~flip,
  1,         "SEAK Juvenile Zone",-136.2,   54.8,    1.15,    "#B5B9CC",  20,    FALSE,
  2,         "GOA Northern Arc",  -145.2,   57.5,    0.85,    "#C2C6D8",  -5,    FALSE, # Scaled down for lat 57.5 N
  3,         "Western AK Arc",    -151.8,   54.2,    1.20,    "#A8ACBF", -30,    FALSE
)

# Base Map Geographic Labels
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

# Salmon depth color gradient
depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# 5. Build Map Layer
p <- ggplot() +
  geom_sf(data = land, fill = "#d2d8b0", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#849e92", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  geom_text(
    data = labels_df,
    aes(x = x, y = y, label = label),
    color = "white",
    fontface = labels_df$fontface,
    size = labels_df$size,
    lineheight = 0.9,
    hjust = 0.5
  )

# Add Salmon Silhouettes
for (i in 1:nrow(salmon_stages)) {
  stg <- salmon_stages[i, ]
  p <- p + add_phylopic(
    uuid = chinook_uuid,
    x = stg$lon,
    y = stg$lat,
    height = stg$icon_size,
    fill = stage_colors[i],
    color = stage_colors[i],
    angle = stg$angle,
    horizontal = stg$flip
  )
}

# Add Version 2 Shark Silhouettes
for (i in 1:nrow(shark_stages_v2)) {
  shk <- shark_stages_v2[i, ]
  p <- p + add_phylopic(
    uuid = shark_uuid,
    x = shk$lon,
    y = shk$lat,
    height = shk$height,
    fill = shk$color,
    color = shk$color,
    angle = shk$angle,
    horizontal = shk$flip
  )
}

# Crop map bounds
p <- p + coord_sf(
  xlim = c(-170, -119),
  ylim = c(43, 62),
  expand = FALSE
)

# 6. Save Output Image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "map_layer_with_salmon_and_sharks_warm.png"),
  plot = p,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Success! Saved map with resized top shark.")