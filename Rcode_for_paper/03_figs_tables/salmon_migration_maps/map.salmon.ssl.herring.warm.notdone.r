# ==============================================================================
# Script: 04_visuals/05_oceanic_salmon_forage_and_ssl_warm.R
# Purpose: Warm Map - Subspecies Sea Lions Further Offshore Intersecting Salmon
# Output: map_layer_with_forage_and_ssl_warm.png
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

# 2. Base vector landmass & PhyloPic UUID retrieval
land <- ne_countries(scale = "medium", returnclass = "sf")

chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
# Specific subspecies for sea lion silhouette
ssl_uuid     <- get_uuid(name = "Zalophus californianus wollebaeki")
herring_uuid <- get_uuid(name = "Clupea pallasii")
capelin_uuid <- get_uuid(name = "Mallotus villosus")

# 3. Salmon Migration Stages (Shifted offshore & deeper in Alaskan waters)
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,     ~lat,    ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",           -125.5,   46.8,    100,        10,       40,     FALSE,
  2,      "WCVI Coastal",      -128.2,   49.5,    150,        30,       30,     FALSE,
  3,      "BC Coast",          -132.2,   52.2,    150,        50,       30,     FALSE,
  4,      "SEAK Shelf",        -137.8,   55.8,    200,        80,       15,     FALSE, # Shifted offshore & deeper
  5,      "GOA Top Shelf",     -147.5,   58.2,    120,       100,        0,     FALSE, # Shifted offshore & deeper
  6,      "GOA West Offshore", -153.5,   55.5,    310,       120,     -35,     FALSE,
  7,      "GOA Central",       -147.0,   51.5,    360,       150,      20,     TRUE,
  8,      "GOA Swing Back",    -137.5,   48.2,    410,       180,      10,     TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Reduced Forage Fish (2 fish per school)
forage_schools <- tribble(
  ~species,  ~uuid,        ~lon,     ~lat,    ~height, ~angle,
  # WCVI / Haida Gwaii Herring (2 fish)
  "Herring", herring_uuid, -129.8,   51.2,    0.45,     15,
  "Herring", herring_uuid, -129.4,   51.0,    0.40,    -10,
  
  # Sitka Herring (2 fish)
  "Herring", herring_uuid, -135.9,   57.1,    0.45,     20,
  "Herring", herring_uuid, -135.7,   56.9,    0.40,     -5,
  
  # Middleton Capelin (2 fish)
  "Capelin", capelin_uuid, -146.3,   59.2,    0.40,    -15,
  "Capelin", capelin_uuid, -146.0,   59.0,    0.35,     10
)

# 5. Offshore Sea Lions (Positioned further out than salmon, swimming inward/overlapping)
ssl_offshore <- tribble(
  ~location,                  ~uuid,    ~lon,     ~lat,    ~height, ~angle, ~flip,
  "Far Offshore WCVI",        ssl_uuid, -130.5,   48.8,    1.00,    -15,    TRUE,   # Outside Stage 2, facing inward
  "Far Offshore SEAK",        ssl_uuid, -139.8,   54.8,    0.95,    -25,    TRUE,   # Outside Stage 4, swimming toward salmon
  "Far Offshore GOA North",   ssl_uuid, -149.2,   57.0,    0.80,     15,    FALSE   # Outside Stage 5, overlapping curve
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

# 6. Build Warm Map Layer
p <- ggplot() +
  geom_sf(data = land, fill = "#d2d8b0", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    # Warm ocean tone (#849e92)
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

# Add Forage Fish Silhouettes (Silver-grey #D1D5DB)
for (i in 1:nrow(forage_schools)) {
  fg <- forage_schools[i, ]
  p <- p + add_phylopic(
    uuid = fg$uuid,
    x = fg$lon,
    y = fg$lat,
    height = fg$height,
    fill = "#D1D5DB",
    color = "#9CA3AF",
    angle = fg$angle
  )
}

# Add Far-Offshore Sea Lions (Espresso #372E29)
for (i in 1:nrow(ssl_offshore)) {
  sl <- ssl_offshore[i, ]
  p <- p + add_phylopic(
    uuid = sl$uuid,
    x = sl$lon,
    y = sl$lat,
    height = sl$height,
    fill = "#372E29",
    color = "#1F1A17",
    angle = sl$angle,
    horizontal = sl$flip
  )
}

# Crop map bounds
p <- p + coord_sf(
  xlim = c(-170, -119),
  ylim = c(43, 62),
  expand = FALSE
)

# 7. Save Output Image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "map_layer_with_forage_and_ssl_warm.png"),
  plot = p,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Success! Saved warm ocean map with subspecies sea lions positioned further offshore.")