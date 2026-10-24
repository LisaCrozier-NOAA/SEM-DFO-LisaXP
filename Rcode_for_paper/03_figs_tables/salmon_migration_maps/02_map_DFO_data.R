# ==============================================================================
# Script: 03_figs_tables/02_map_DFO_data.R
# Purpose: Map showing 1st-year smolt outmigration + returning adult, plus
#          forage schools (herring, capelin, krill), marine predators 
#          (sea lions, harbor seals, killer whales), and a BC maple leaf.
# Output: map_layer_with_forage_and_ssl_cool.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "03_figs_tables","salmon_migration_maps")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass & PhyloPic UUID retrieval
land <- ne_countries(scale = "medium", returnclass = "sf")

chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
ssl_uuid     <- get_uuid(name = "Zalophus californianus wollebaeki") # Sea lion
herring_uuid <- get_uuid(name = "Clupea pallasii")                   # Pacific herring
capelin_uuid <- get_uuid(name = "Mallotus villosus")                 # Capelin
orca_uuid    <- get_uuid(name = "Orcinus orca")                      # Killer whale
seal_uuid    <- get_uuid(name = "Phoca vitulina")                    # Harbour seal
krill_uuid   <- get_uuid(name = "Euphausia pacifica")                # Krill
maple_uuid   <- get_uuid(name = "Acer")                              # Maple leaf

# 3. Salmon Migration Stages: 1st Year (Stages 1-5) + Returning Adult (Stage 8)
salmon_stages <- tribble(
  ~stage, ~location,                 ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",                -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",           -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",               -132.2,   52.2,   150,        50,        30,   FALSE,
#  4,      "SEAK Shelf",             -136.5,   56.2,   200,        70,        15,   FALSE,
#  5,      "GOA Top Shelf",          -146.0,   58.2,   240,        90,         0,   FALSE,
  8,      "GOA Swing Back (Adult)", -130.5,   45.2,   450,       180,        10,   TRUE
) %>%
  mutate(
    icon_size = 0.60 + (length_mm - 100) * (0.80 / 350)
  )

# 4. Forage Fish & Krill Schools
forage_schools <- tribble(
  ~species,  ~uuid,        ~lon,    ~lat,    ~height, ~angle, ~fill,     ~color,
  # WCVI / Haida Gwaii Herring School
  "Herring", herring_uuid, -129.8,  51.2,    0.45,     15,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -129.4,  51.0,    0.40,    -10,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -130.1,  51.4,    0.42,      5,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -129.2,  51.3,    0.38,     20,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -129.7,  50.8,    0.40,     -5,    "#D1D5DB", "#9CA3AF",
  
  # Sitka Herring School
  "Herring", herring_uuid, -135.9,  57.1,    0.45,     20,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -135.7,  56.9,    0.40,     -5,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -136.3,  57.3,    0.42,     10,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -135.4,  57.1,    0.38,     15,    "#D1D5DB", "#9CA3AF",
  "Herring", herring_uuid, -136.1,  56.7,    0.40,      0,    "#D1D5DB", "#9CA3AF",
  
  # Middleton Island Capelin School
  "Capelin", capelin_uuid, -146.3,  59.2,    0.40,    -15,    "#E0F2FE", "#7DD3FC",
  "Capelin", capelin_uuid, -146.0,  59.0,    0.35,     10,    "#E0F2FE", "#7DD3FC",
  "Capelin", capelin_uuid, -146.6,  59.3,    0.38,     -5,    "#E0F2FE", "#7DD3FC",
  "Capelin", capelin_uuid, -145.8,  59.2,    0.35,      5,    "#E0F2FE", "#7DD3FC",
  "Capelin", capelin_uuid, -146.4,  58.8,    0.37,    -10,    "#E0F2FE", "#7DD3FC",
  
  # WCVI Krill Cluster (Offshore shelf break)
  "Krill",   krill_uuid,   -127.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -130.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -132.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -131.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444"
)

# 5. Marine Mammals (Steller Sea Lions, Harbour Seals, Killer Whales)
marine_mammals <- tribble(
  ~group,         ~uuid,      ~lon,    ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  # Steller Sea Lions (Haulouts)
  "SSL",          ssl_uuid,   -128.2,  50.5,    1.50,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -136.1,  58.3,    1.35,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -146.3,  60.0,    1.30,      0,    FALSE, "#372E29", "#1F1A17",
  
  # Harbour Seals (BC Coast / Coastal Vancouver Island)
  "Harbour Seal", seal_uuid,  -125.4,  49.8,    0.85,     15,    FALSE, "#4B5563", "#1F2937",
  "Harbour Seal", seal_uuid,  -130.6,  53.0,    0.80,    -10,    TRUE,  "#4B5563", "#1F2937",
  
  # Killer Whales (Northern Resident / BC Offshore Shelf)
  "Orca",         orca_uuid,  -133.8,  51.6,    2.25,     25,    FALSE, "#111827", "#000000",
  "Orca",         orca_uuid,  -136.5,  53.7,    2.15,     30,    FALSE, "#111827", "#000000"
)

# 6. Maple Leaf Badge on British Columbia / Vancouver Island
maple_leaf_df <- tibble(
  lon    = -124.5,
  lat    =  51.2,
  height =  1.20
)
maple_leaf_df2 <- tibble(
  lon    = -126.5,
  lat    =  53.2,
  height =  1.20
)
maple_leaf_df3 <- tibble(
  lon    = -132,
  lat    =  54,
  height =  1.20
)


# Base Map Geographic Labels
labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   61.5,     "bold",      5.2,
  "British\nColumbia",  -122.5,   53.5,     "plain",     4.5,
  "WCVI",               -125.8,   48.8,     "plain",     4.2,
  "Columbia \nRiver",   -121.25,  45.5,     "plain",     4.5,
  "USA",                -120.75,  47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5
)

# Salmon depth color gradient
depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors  <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# 7. Build Map Layer
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

# Add Maple Leaf Symbol on BC
p <- p + add_phylopic(
  uuid   = maple_uuid,
  x      = maple_leaf_df$lon,
  y      = maple_leaf_df$lat,
  height = maple_leaf_df$height,
  fill   = "#B91C1C", # Canadian flag crimson
  color  = "#7F1D1D",
  alpha  = 0.85
)
p <- p + add_phylopic(
  uuid   = maple_uuid,
  x      = maple_leaf_df2$lon,
  y      = maple_leaf_df2$lat,
  height = maple_leaf_df2$height,
  fill   = "#B91C1C", # Canadian flag crimson
  color  = "#7F1D1D",
  alpha  = 0.85
)
p <- p + add_phylopic(
  uuid   = maple_uuid,
  x      = maple_leaf_df3$lon,
  y      = maple_leaf_df3$lat,
  height = maple_leaf_df3$height,
  fill   = "#B91C1C", # Canadian flag crimson
  color  = "#7F1D1D",
  alpha  = 0.85
)

# Add Salmon Silhouettes (Smolts 1-5 + Returning Adult 8)
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

# Add Forage Species (Herring, Capelin, Krill)
for (i in seq_len(nrow(forage_schools))) {
  fg <- forage_schools[i, ]
  p <- p + add_phylopic(
    uuid   = fg$uuid,
    x      = fg$lon,
    y      = fg$lat,
    height = fg$height,
    fill   = fg$fill,
    color  = fg$color,
    angle  = fg$angle
  )
}

# Add Marine Mammals (SSL, Harbour Seals, Orcas)
for (i in seq_len(nrow(marine_mammals))) {
  mm <- marine_mammals[i, ]
  p <- p + add_phylopic(
    uuid       = mm$uuid,
    x          = mm$lon,
    y          = mm$lat,
    height     = mm$height,
    fill       = mm$fill,
    color      = mm$color,
    angle      = mm$angle,
    horizontal = mm$flip
  )
}

# Crop map bounds
p <- p + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# 8. Save Output Image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "02_map_DFO_data.png"),
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

message("Success! Saved updated migration map with mammals, forage, and BC maple leaf.")