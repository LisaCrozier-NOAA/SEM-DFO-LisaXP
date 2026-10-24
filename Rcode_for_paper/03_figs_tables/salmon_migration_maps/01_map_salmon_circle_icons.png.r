# ==============================================================================
# Script: 03_figs_tables/salmon_migration_maps/01_map_salmon_circle_icons.png.r
# Purpose: High-Contrast Light Ocean Base with Expanded Multi-Year Rearing Gyre
# Output: map_layer_with_salmon.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "03_figs_tables", "salmon_migration_maps")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass & PhyloPic setup
land         <- ne_countries(scale = "medium", returnclass = "sf")
chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")

# 3. Original 8 stages preserved intact + 4 additional immature gyre stages
# Stages 1-8 are identical to your original coordinate configuration
salmon_stages <- tribble(
  ~stage, ~location,                 ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  # --- Original 8 stages ---
  1,      "CR Exit",                -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",           -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",               -132.2,   52.2,   150,        50,        30,   FALSE,
  4,      "SEAK Shelf",             -136.5,   56.2,   200,        70,        15,   FALSE,
  5,      "GOA Top Shelf",          -146.0,   59.2,   120,        90,         0,   FALSE,
  6,      "GOA West Offshore",      -153.5,   55.5,   310,       120,       -35,   FALSE,
  7,      "GOA Central",            -147.0,   51.5,   510,       150,        20,   TRUE,
  8,      "GOA Swing Back",         -137.5,   48.2,   410,       180,        10,   TRUE,
  
  # --- 4 Added stages completing the multi-year counterclockwise gyre circuit ---
  9,      "Gyre East Bound",        -141.5,   52.8,   510,       130,       -75,   TRUE,
  10,     "Gyre North Re-entry",    -148.5,   54.8,   330,       140,       -20,   FALSE,
  11,     "Gyre Southwest Sweep",   -156.0,   52.5,   400,       160,       -75,   FALSE,
  12,     "Gyre South Basin",       -150.0,   48.8,   500,       170,        30,   TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# Base Map Geographic Labels
labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   61.5,     "bold",      5.2,
  "British\nColumbia",  -123.5,   52.5,     "plain",     4.5,
  "WCVI",               -126.0,   49.5,     "plain",     4.2,
  "Columbia \nRiver",   -121.25,  45.5,     "plain",     4.5,
  "USA",                -120.75,  47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5
)

# COLOR GRADIENT: Light Cyan (Surface Smolts) -> Deep Navy Blue (Deep Immatures/Adults)
depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors  <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# 4. Build Map Layer (Pure icon layout, no line underlays)
p <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
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

# Render salmon silhouettes
for (i in seq_len(nrow(salmon_stages))) {
  stg <- salmon_stages[i, ]
  p <- p + add_phylopic(
    uuid       = chinook_uuid,
    x          = stg$lon,
    y          = stg$lat,
    height     = stg$icon_size,
    fill       = stage_colors[i],
    color      = stage_colors[i],
    angle      = stg$angle,
    horizontal = stg$flip
  )
}

# Crop map bounds
p <- p + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# 5. Save output image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "map_layer_with_salmon_circle_icons.png"),
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

message("Success! Saved clean icon map with original 8 points preserved plus 4 gyre rearing stages.")