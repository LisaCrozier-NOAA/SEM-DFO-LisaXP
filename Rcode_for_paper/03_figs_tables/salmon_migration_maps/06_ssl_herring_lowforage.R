# ==============================================================================
# Script: 03_figs_tables/salmon_migration_maps/05_ssl_herring_lowforage.R
# Purpose: Offline Map using local RDS icon cache.
#          Steller Sea Lions scaled up and shifted slightly inshore from salmon track.
# Output: map_layer_with_forage_and_ssl_lowforage.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "03_figs_tables", "salmon_migration_maps", "figs")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass & Local RDS Cache Retrieval (100% Offline)
land <- ne_countries(scale = "medium", returnclass = "sf")

cache_file <- file.path(proj_dir, "03_figs_tables", "phylopic_cache", "phylopic_icons.rds")
if (!file.exists(cache_file)) {
  stop(sprintf("Icon cache not found at: %s. Run the caching script first.", cache_file))
}
icons <- readRDS(cache_file)

# 3. Salmon Migration Stages
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",          -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",     -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",         -132.2,   52.2,   150,        50,        30,   FALSE,
  4,      "SEAK Shelf",       -137.8,   55.8,   200,        80,        15,   FALSE,
  5,      "GOA Top Shelf",    -147.5,   58.2,   240,       100,         0,   FALSE,
  6,      "GOA West Offshore",-153.5,   55.5,   310,       120,       -35,   FALSE,
  7,      "GOA Central",      -147.0,   51.5,   360,       150,        20,   TRUE,
  8,      "GOA Swing Back",   -137.5,   48.2,   410,       180,        10,   TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Reduced Forage Fish (Low forage condition: 2 fish per school)
forage_schools <- list(
  # WCVI / Haida Gwaii Herring (2 fish)
  list(img = icons$herring, lon = -129.8, lat = 51.2, height = 0.45, angle = 15),
  list(img = icons$herring, lon = -129.4, lat = 51.0, height = 0.40, angle = -10),
  
  # Sitka Herring (2 fish)
  list(img = icons$herring, lon = -135.9, lat = 57.1, height = 0.45, angle = 20),
  list(img = icons$herring, lon = -135.7, lat = 56.9, height = 0.40, angle = -5),
  
  # Middleton Capelin (2 fish)
  list(img = icons$capelin, lon = -146.3, lat = 59.2, height = 0.40, angle = -15),
  list(img = icons$capelin, lon = -146.0, lat = 59.0, height = 0.35, angle = 10)
)

# 5. Scaled-Up Sea Lions Positioned Inshore of the Salmon Track
# Shifted eastward toward the coastline so they do not overlap the smolts
ssl_inshore <- list(
  # Inshore of Stage 2 (WCVI Coastal: salmon at -128.2, 49.5) -> shifted to -127.1, 49.9
  list(img = icons$ssl, lon = -127.1, lat = 49.9, height = 1.45, angle = -15, flip = TRUE),
  
  # Inshore of Stage 4 (SEAK Shelf: salmon at -137.8, 55.8) -> shifted to -136.3, 56.4
  list(img = icons$ssl, lon = -136.3, lat = 56.4, height = 1.45, angle = -20, flip = TRUE),
  
  # Inshore/North of Stage 5 (GOA Top Shelf: salmon at -147.5, 58.2) -> shifted to -146.8, 59.5
  list(img = icons$ssl, lon = -146.8, lat = 59.5, height = 1.40, angle = 10, flip = FALSE)
)

# Base Map Geographic Labels
labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   61.5,     "bold",      5.2,
  "British\nColumbia",  -123.5,   52.5,     "plain",     4.5,
  "WCVI",               -125.8,   48.8,     "plain",     4.2,
  "Columbia \nRiver",   -121.25,  45.5,     "plain",     4.5,
  "USA",                -120.75,  47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5
)

# Salmon depth color gradient
depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors  <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# 6. Build Baseline Ocean Map
p <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    # Standard baseline light ocean blue (#6fa2bb)
    panel.background = element_rect(fill = "#6fa2bb", color = NA),
    plot.margin      = margin(0, 0, 0, 0)
  ) +
  geom_text(
    data       = labels_df,
    aes(x = x, y = y, label = label),
    color      = "white",
    fontface   = labels_df$fontface,
    size       = labels_df$size,
    lineheight = 0.9,
    hjust      = 0.5
  )

# Add Salmon Silhouettes from local cache
for (i in seq_len(nrow(salmon_stages))) {
  stg <- salmon_stages[i, ]
  p <- p + add_phylopic(
    img        = icons$chinook,
    x          = stg$lon,
    y          = stg$lat,
    height     = stg$icon_size,
    fill       = stage_colors[i],
    color      = stage_colors[i],
    angle      = stg$angle,
    horizontal = stg$flip
  )
}

# Add Forage Fish Silhouettes from local cache
for (fg in forage_schools) {
  p <- p + add_phylopic(
    img    = fg$img,
    x      = fg$lon,
    y      = fg$lat,
    height = fg$height,
    fill   = "#D1D5DB",
    color  = "#9CA3AF",
    angle  = fg$angle
  )
}

# Add Scaled-Up Inshore Sea Lions from local cache
for (sl in ssl_inshore) {
  p <- p + add_phylopic(
    img        = sl$img,
    x          = sl$lon,
    y          = sl$lat,
    height     = sl$height,
    fill       = "#372E29",
    color      = "#1F1A17",
    angle      = sl$angle,
    horizontal = sl$flip
  )
}

# Crop map bounds
p <- p + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# 7. Save Output Image
print(p)

file_out <- file.path(fig_out_dir, "map_layer_with_forage_and_ssl_lowforage.png")
ggsave(
  filename = file_out,
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

cat(sprintf("\n[Success] Offline map saved to: %s\n", file_out))