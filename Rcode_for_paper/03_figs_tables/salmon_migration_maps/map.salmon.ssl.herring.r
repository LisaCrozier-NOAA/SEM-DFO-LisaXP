# ==============================================================================
# Script: 04_visuals/04_oceanic_salmon_forage_and_ssl_cool.R
# Purpose: Cool Map - Salmon Loop + Expanded Forage Schools + Scaled-Up SSLs
# Output: map_layer_with_forage_and_ssl_cool.png
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
#ssl_uuid     <- get_uuid(name = "Eumetopias jubatus")  # Steller sea lion
ssl_uuid <- get_uuid(name = "Zalophus californianus wollebaeki")
herring_uuid <- get_uuid(name = "Clupea pallasii")     # Pacific herring
capelin_uuid <- get_uuid(name = "Mallotus villosus")    # Capelin

# 3. Salmon Migration Stages
salmon_stages <- tribble(
  ~stage, ~location,           ~lon,     ~lat,    ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",           -125.5,   46.8,    100,        10,       40,     FALSE,
  2,      "WCVI Coastal",      -128.2,   49.5,    150,        30,       30,     FALSE,
  3,      "BC Coast",          -132.2,   52.2,    150,        50,       30,     FALSE,
  4,      "SEAK Shelf",        -136.5,   56.2,    200,        70,       15,     FALSE,
  5,      "GOA Top Shelf",     -146.0,   58.2,    120,        90,        0,     FALSE,
  6,      "GOA West Offshore", -153.5,   55.5,    310,       120,     -35,     FALSE,
  7,      "GOA Central",       -147.0,   51.5,    360,       150,      20,     TRUE,
  8,      "GOA Swing Back",    -137.5,   48.2,    410,       180,      10,     TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Forage Fish Schools (5 Fish per School)
forage_schools <- tribble(
  ~species,  ~uuid,        ~lon,     ~lat,    ~height, ~angle,
  # WCVI / Haida Gwaii Herring School (5 fish)
  "Herring", herring_uuid, -129.8,   51.2,    0.45,     15,
  "Herring", herring_uuid, -129.4,   51.0,    0.40,    -10,
  "Herring", herring_uuid, -130.1,   51.4,    0.42,      5,
  "Herring", herring_uuid, -129.2,   51.3,    0.38,     20,
  "Herring", herring_uuid, -129.7,   50.8,    0.40,     -5,
  
  # Sitka Herring School (5 fish)
  "Herring", herring_uuid, -135.9,   57.1,    0.45,     20,
  "Herring", herring_uuid, -135.7,   56.9,    0.40,     -5,
  "Herring", herring_uuid, -136.3,   57.3,    0.42,     10,
  "Herring", herring_uuid, -135.4,   57.1,    0.38,     15,
  "Herring", herring_uuid, -136.1,   56.7,    0.40,      0,
  
  # Middleton Island Capelin School (5 fish)
  "Capelin", capelin_uuid, -146.3,   59.2,    0.40,    -15,
  "Capelin", capelin_uuid, -146.0,   59.0,    0.35,     10,
  "Capelin", capelin_uuid, -146.6,   59.3,    0.38,     -5,
  "Capelin", capelin_uuid, -145.8,   59.2,    0.35,      5,
  "Capelin", capelin_uuid, -146.4,   58.8,    0.37,    -10
)

# 5. Steller Sea Lions (SSL) Sitting on Coastal Land (SCALED UP SIZES)
ssl_haulouts <- tribble(
  ~location,            ~uuid,    ~lon,     ~lat,    ~height, ~angle, ~flip,
  "WCVItop",            ssl_uuid, -128.2,   50.5,    1.50,    -10,    TRUE,   # Increased height to 1.50
  "Sitka Coastal Land", ssl_uuid, -136.1,   58.3,    1.35,    -10,    TRUE,   # Increased height to 1.35
  "Middleton Is. Land", ssl_uuid, -146.3,   60.0,    1.30,      0,    FALSE   # Increased height to 1.30
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

# 6. Build Cool Map Layer
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

# Add Forage Fish Silhouettes
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

# Add Scaled-Up Steller Sea Lion Silhouettes
for (i in 1:nrow(ssl_haulouts)) {
  sl <- ssl_haulouts[i, ]
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
  filename = file.path(fig_out_dir, "map_layer_with_forage_and_ssl_cool.png"),
  plot = p,
  width = 10,
  height = 7.5,
  dpi = 300
)

message("Success! Saved cool map with larger Steller Sea Lions.")