# ==============================================================================
# Script: 03_figs_tables/03_oceanic_salmon_and_sharks_layer.R
# Purpose: Add 3 Deep Dark Grey Sharks to Central GOA Inside Migration Loop
# Output: map_layer_with_salmon_and_sharks_v1.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "03_figs_tables","salmon_migration_maps")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass & PhyloPic setup
land <- ne_countries(scale = "medium", returnclass = "sf")

# Retrieve silhouettes for Chinook salmon and Salmon Shark (Lamna ditropis)
chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
shark_uuid   <- get_uuid(name = "sharks") # Salmon shark silhouette

# 3. Salmon Migration Stages
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,     ~lat,    ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",           -125.5,   46.8,    100,        10,       40,     FALSE,
  2,      "WCVI Coastal",      -128.2,   49.5,    150,        30,       30,     FALSE,
  3,      "BC Coast",          -132.2,   52.2,    150,        50,       30,     FALSE,
  4,      "SEAK Shelf",        -136.5,   56.2,    200,        70,       15,     FALSE,
  5,      "GOA Top Shelf",     -146.0,   59.2,    120,        90,        0,     FALSE,
  6,      "GOA West Offshore", -153.5,   55.5,    310,       120,     -35,     FALSE,
  7,      "GOA Central",       -147.0,   51.5,    360,       150,      20,     TRUE,
  8,      "GOA Swing Back",    -137.5,   48.2,    410,       180,      10,     TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Version 1: 3 Deep-Water Sharks (Central GOA Inside Loop)
# Dark shades of grey (#2B2D42, #1A1A24, #3A3D52) representing deep depths
shark_stages_v1 <- tribble(
  ~shark_id, ~location,          ~lon,     ~lat,    ~height, ~color,    ~angle, ~flip,
  1,         "Central GOA North",-144.5,   55.2,    1.30,    "#2B2D42", -20,    FALSE,
  2,         "Central GOA West", -149.0,   53.8,    1.45,    "#1A1A24",  15,    TRUE,
  3,         "Central GOA East", -141.5,   52.5,    1.25,    "#3A3D52", -10,    FALSE
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

# Salmon depth color gradient (Light Cyan -> Deep Cobalt Blue)
depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# 5. Build Base Map Layer
p <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#6fa2bb", color = NA),
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

# Add Version 1 Shark Silhouettes (Dark Grey, Central GOA)
for (i in 1:nrow(shark_stages_v1)) {
  shk <- shark_stages_v1[i, ]
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

# 6. Save Version 1 Image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "map_layer_with_salmon_and_sharks_cool.png"),
  plot = p,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Success! Saved Version 1 (Central Sharks) to: ", 
        file.path(fig_out_dir, "map_layer_with_salmon_and_sharks_v1.png"))