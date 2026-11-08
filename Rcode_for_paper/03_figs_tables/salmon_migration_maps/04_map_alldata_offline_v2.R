# ==============================================================================
# Script: 03_figs_tables/04_map_alldata_offline.R
# Purpose: Comprehensive Ecosystem & Monitoring Map (100% Offline via local RDS cache).
#          - Krill cluster pushed far offshore into deep ocean (completely off land).
#          - Seward Line copepod and capelin pushed south into open GOA waters.
#          - SEAK copepod pushed offshore clear of Baranof/Chichagof islands.
# Output: 03_figs_tables/salmon_migration_maps/figs/04_map_alldata_offline.png
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

# 2. Base vector landmass
land <- ne_countries(scale = "medium", returnclass = "sf")

# 3. Load Offline PhyloPic RDS Cache
cache_file <- file.path(proj_dir, "03_figs_tables", "phylopic_cache", "phylopic_icons.rds")
if (!file.exists(cache_file)) {
  stop(sprintf("Icon cache not found at: %s. Run your caching script first.", cache_file))
}
icons <- readRDS(cache_file)

# ==============================================================================
# DATA DEFINITIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# A. Salmon Stages (Smolts 1-3 + Returning Adult 8)
# ------------------------------------------------------------------------------
salmon_stages <- tribble(
  ~stage, ~location,                 ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",                -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",           -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",               -132.2,   52.2,   150,        50,        30,   FALSE,
  8,      "GOA Swing Back (Adult)", -130.5,   45.2,   450,       180,        10,   TRUE
) %>%
  mutate(icon_size = 0.60 + (length_mm - 100) * (0.80 / 350))

depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors  <- depth_palette(100)[as.numeric(cut(salmon_stages$depth_m, breaks = 100))]

# ------------------------------------------------------------------------------
# B. Canadian Maple Leaves
# ------------------------------------------------------------------------------
maple_leaves <- list(
  list(lon = -124.5, lat = 51.2, height = 1.20),
  list(lon = -126.5, lat = 53.2, height = 1.20),
  list(lon = -132.0, lat = 54.0, height = 1.20)
)

# ------------------------------------------------------------------------------
# C. Northern Mammals
# ------------------------------------------------------------------------------
marine_mammals <- list(
  # Steller Sea Lions
  list(img = icons$ssl,  lon = -127.5, lat = 50.8, height = 1.45, angle = -10, flip = TRUE,  fill = "#372E29", color = "#1F1A17"),
  list(img = icons$ssl,  lon = -135.5, lat = 58.6, height = 1.35, angle = -10, flip = TRUE,  fill = "#372E29", color = "#1F1A17"),
  list(img = icons$ssl,  lon = -146.3, lat = 60.3, height = 1.30, angle =   0, flip = FALSE, fill = "#372E29", color = "#1F1A17"),
  
  # Harbour Seals (Outer Pacific coast, facing west)
  list(img = icons$seal, lon = -127.2, lat = 49.6, height = 0.85, angle =  15, flip = FALSE, fill = "#4B5563", color = "#1F2937"),
  list(img = icons$seal, lon = -132.5, lat = 52.6, height = 0.85, angle = -15, flip = FALSE, fill = "#4B5563", color = "#1F2937"),
  
  # Killer Whales
  list(img = icons$orca, lon = -133.8, lat = 51.6, height = 2.25, angle =  25, flip = FALSE, fill = "#111827", color = "#000000"),
  list(img = icons$orca, lon = -136.5, lat = 53.7, height = 2.15, angle =  30, flip = FALSE, fill = "#111827", color = "#000000")
)

# ------------------------------------------------------------------------------
# D. Coastal & Gulf Forage Schools
# ------------------------------------------------------------------------------
forage_schools <- list(
  # Queen Charlotte Sound Herring School
  list(img = icons$herring, lon = -131.2, lat = 51.4, height = 0.45, angle =  15, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -130.8, lat = 51.1, height = 0.40, angle = -10, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -131.5, lat = 51.5, height = 0.42, angle =   5, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -130.6, lat = 51.3, height = 0.38, angle =  20, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -131.0, lat = 50.8, height = 0.40, angle =  -5, fill = "#D1D5DB", color = "#9CA3AF"),
  
  # Sitka Herring School
  list(img = icons$herring, lon = -137.0, lat = 56.8, height = 0.45, angle =  20, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -136.6, lat = 56.6, height = 0.40, angle =  -5, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -137.3, lat = 57.0, height = 0.42, angle =  10, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -136.3, lat = 56.8, height = 0.38, angle =  15, fill = "#D1D5DB", color = "#9CA3AF"),
  list(img = icons$herring, lon = -136.8, lat = 56.4, height = 0.40, angle =   0, fill = "#D1D5DB", color = "#9CA3AF"),
  
  # Middleton Island Capelin School (Shifted south into open water 57.2N - 57.6N)
  list(img = icons$capelin, lon = -146.5, lat = 58.5, height = 0.40, angle = -15, fill = "#E0F2FE", color = "#7DD3FC"),
  list(img = icons$capelin, lon = -146.1, lat = 58.3, height = 0.35, angle =  10, fill = "#E0F2FE", color = "#7DD3FC"),
  list(img = icons$capelin, lon = -146.8, lat = 58.6, height = 0.38, angle =  -5, fill = "#E0F2FE", color = "#7DD3FC"),
  list(img = icons$capelin, lon = -145.9, lat = 58.4, height = 0.35, angle =   5, fill = "#E0F2FE", color = "#7DD3FC"),
  list(img = icons$capelin, lon = -146.4, lat = 58.1, height = 0.37, angle = -10, fill = "#E0F2FE", color = "#7DD3FC"),
  # list(img = icons$capelin, lon = -146.5, lat = 57.5, height = 0.40, angle = -15, fill = "#E0F2FE", color = "#7DD3FC"),
  # list(img = icons$capelin, lon = -146.1, lat = 57.3, height = 0.35, angle =  10, fill = "#E0F2FE", color = "#7DD3FC"),
  # list(img = icons$capelin, lon = -146.8, lat = 57.6, height = 0.38, angle =  -5, fill = "#E0F2FE", color = "#7DD3FC"),
  # list(img = icons$capelin, lon = -145.9, lat = 57.4, height = 0.35, angle =   5, fill = "#E0F2FE", color = "#7DD3FC"),
  # list(img = icons$capelin, lon = -146.4, lat = 57.1, height = 0.37, angle = -10, fill = "#E0F2FE", color = "#7DD3FC"),
  
  # WCVI Krill Cluster (Shifted WEST to -130.5W to -133.5W: COMPLETELY IN BLUE WATER)
  list(img = icons$krill,   lon = -130.8, lat = 48.6, height = 0.32, angle =  35, fill = "#FCA5A5", color = "#EF4444"),
  list(img = icons$krill,   lon = -131.4, lat = 48.3, height = 0.28, angle =  15, fill = "#FCA5A5", color = "#EF4444"),
  list(img = icons$krill,   lon = -131.0, lat = 48.8, height = 0.30, angle =  50, fill = "#FCA5A5", color = "#EF4444"),
  list(img = icons$krill,   lon = -132.5, lat = 48.5, height = 0.32, angle =  35, fill = "#FCA5A5", color = "#EF4444"),
  list(img = icons$krill,   lon = -133.5, lat = 48.2, height = 0.28, angle =  15, fill = "#FCA5A5", color = "#EF4444"),
  list(img = icons$krill,   lon = -133.0, lat = 48.6, height = 0.30, angle =  50, fill = "#FCA5A5", color = "#EF4444")
)

# ------------------------------------------------------------------------------
# E. Ferriss et al. 2025 Regional ESR Arenas (Shifted Safely into Open Waters)
# ------------------------------------------------------------------------------
ferriss_arenas <- list(
  # WGOA Low-Trophic: Seward Line (Pushed south to 58.2N, well off Kenai Peninsula)
  list(img = icons$copepod, lon = -149.2, lat = 58.2, height = 0.48, angle =  25, flip = FALSE, fill = "#FEF08A", color = "#CA8A04"),
  # WGOA Mid-Trophic: Shelikof / Kodiak shelf (Pushed south to 56.2N)
  list(img = icons$capelin, lon = -153.5, lat = 56.2, height = 0.48, angle = -15, flip = FALSE, fill = "#7DD3FC", color = "#0284C7"),
  # WGOA Seabirds: Middleton Island
  list(img = icons$murre,   lon = -145.4, lat = 58.8, height = 0.65, angle =   5, flip = FALSE, fill = "#0F172A", color = "#020617"),
  
  # EGOA Low-Trophic: SECM Shelf (Pushed offshore to -137.8W, 57.2N clear of land)
  list(img = icons$copepod, lon = -137.8, lat = 57.2, height = 0.48, angle =  15, flip = FALSE, fill = "#FEF08A", color = "#CA8A04"),
  # EGOA Seabirds: St. Lazaria Island
  list(img = icons$murre,   lon = -135.2, lat = 56.3, height = 0.65, angle = -10, flip = FALSE, fill = "#0F172A", color = "#020617")
)

# ------------------------------------------------------------------------------
# F. DFO Line P Transect
# ------------------------------------------------------------------------------
line_p_zooplankton <- list(
  list(img = icons$copepod,     lon = -130.0, lat = 48.9, height = 0.40, angle =  10, fill = "#FDE047", color = "#CA8A04"),
  list(img = icons$krill,       lon = -134.5, lat = 49.2, height = 0.38, angle =  25, fill = "#F87171", color = "#DC2626"),
  list(img = icons$larval_fish, lon = -139.0, lat = 49.5, height = 0.40, angle = -15, fill = "#FDE68A", color = "#D97706"),
  list(img = icons$crab,        lon = -143.5, lat = 49.8, height = 0.38, angle =   0, fill = "#FB923C", color = "#C2410C")
)

# ------------------------------------------------------------------------------
# G. Gulf of Alaska Demersal Predators & Sharks
# ------------------------------------------------------------------------------
goa_predators <- list(
  # Walleye Pollock
  list(img = icons$pollock,   lon = -149.5, lat = 57.3, height = 0.75, angle =  10, flip = TRUE,  fill = "#64748B", color = "#334155"),
  list(img = icons$pollock,   lon = -148.0, lat = 57.8, height = 0.65, angle = -15, flip = TRUE,  fill = "#64748B", color = "#334155"),
  list(img = icons$pollock,   lon = -146.5, lat = 56.2, height = 0.65, angle = -15, flip = TRUE,  fill = "#64748B", color = "#334155"),
  
  # Cod, Halibut, Flounder, Sablefish, Sharks
  list(img = icons$cod,       lon = -154.5, lat = 54.5, height = 0.85, angle =  20, flip = FALSE, fill = "#52525B", color = "#27272A"),
  list(img = icons$halibut,   lon = -147.2, lat = 55.0, height = 1.05, angle = -25, flip = FALSE, fill = "#44403C", color = "#1C1917"),
  list(img = icons$flounder,  lon = -143.5, lat = 56.5, height = 0.80, angle =  15, flip = TRUE,  fill = "#78716C", color = "#292524"),
  list(img = icons$sablefish, lon = -143.0, lat = 55.2, height = 0.85, angle = -30, flip = TRUE,  fill = "#2B547E", color = "#1E3A5F"),
  list(img = icons$sablefish, lon = -141.0, lat = 55.2, height = 0.85, angle = -30, flip = TRUE,  fill = "#2B547E", color = "#1E3A5F"),
  list(img = icons$shark,     lon = -145.0, lat = 52.8, height = 1.30, angle = -10, flip = FALSE, fill = "#1E293B", color = "#0F172A"),
  list(img = icons$shark,     lon = -151.0, lat = 52.0, height = 1.15, angle =  20, flip = TRUE,  fill = "#334155", color = "#0F172A")
)

# ------------------------------------------------------------------------------
# H. Northern California Current
# ------------------------------------------------------------------------------
nh_zooplankton <- list(
  list(img = icons$copepod,     lon = -125.8, lat = 44.6, height = 0.45, angle = 20, fill = "#FDE047", color = "#CA8A04"),
  list(img = icons$crab,        lon = -127.0, lat = 44.6, height = 0.45, angle =  0, fill = "#FB923C", color = "#C2410C"),
  list(img = icons$larval_fish, lon = -128.2, lat = 44.6, height = 0.45, angle = 15, fill = "#FDE68A", color = "#D97706"),
  list(img = icons$krill,       lon = -129.4, lat = 44.6, height = 0.42, angle = 30, fill = "#F87171", color = "#DC2626")
)

ncc_shelf_fauna <- list(
  list(img = icons$hake,        lon = -127.2, lat = 47.0, height = 0.85, angle =  20, flip = FALSE, fill = "#38BDF8", color = "#0284C7"),
  list(img = icons$rockfish,    lon = -126.2, lat = 46.2, height = 0.70, angle = -15, flip = FALSE, fill = "#FB7185", color = "#E11D48"),
  list(img = icons$squid,       lon = -128.2, lat = 46.1, height = 0.65, angle =  35, flip = FALSE, fill = "#C084FC", color = "#9333EA"),
  list(img = icons$sea_nettle,  lon = -127.2, lat = 45.3, height = 0.80, angle = -10, flip = FALSE, fill = "#F472B6", color = "#DB2777"),
  list(img = icons$forage_fish, lon = -126.0, lat = 47.6, height = 0.40, angle =  15, flip = FALSE, fill = "#94A3B8", color = "#475569"),
  list(img = icons$forage_fish, lon = -126.4, lat = 47.4, height = 0.35, angle = -10, flip = FALSE, fill = "#94A3B8", color = "#475569")
)

ncc_seabirds <- list(
  list(img = icons$murre,       lon = -128.2, lat = 47.8, height = 0.75, angle =  10, flip = FALSE, fill = "#1E293B", color = "#0F172A"),
  list(img = icons$shearwater,  lon = -129.5, lat = 46.8, height = 0.75, angle =  30, flip = FALSE, fill = "#44403C", color = "#1C1917"),
  list(img = icons$pelican,     lon = -128.8, lat = 45.5, height = 0.95, angle =  15, flip = FALSE, fill = "#78716C", color = "#292524"),
  list(img = icons$murre,       lon = -127.8, lat = 44.0, height = 0.75, angle =   5, flip = FALSE, fill = "#1E293B", color = "#0F172A")
)

columbia_pinnipeds <- list(
  list(img = icons$seal, lon = -124.0, lat = 46.3, height = 0.60, angle = 10, flip = FALSE, fill = "#475569", color = "#1E293B"),
  list(img = icons$ssl,  lon = -123.0, lat = 46.0, height = 0.80, angle = -5, flip = FALSE, fill = "#334155", color = "#0F172A"),
  list(img = icons$ssl,  lon = -122.0, lat = 45.7, height = 0.90, angle =  0, flip = FALSE, fill = "#1E293B", color = "#020617")
)

labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   62.0,     "bold",      5.2,
  "British\nColumbia",  -122.5,   53.5,     "plain",     4.5,
  "WCVI",               -125.8,   48.8,     "plain",     4.2,
  "Columbia \nRiver",   -120.2,   45.8,     "plain",     4.2,
  "USA",                -120.5,   47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5,
  "NH Line",            -130.6,   44.6,     "italic",    3.4,
  "Line P",             -138.5,   48.8,     "italic",    3.4
)

# ==============================================================================
# ASSEMBLE MAP
# ==============================================================================

p <- ggplot() +
  # Basemap
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#6fa2bb", color = NA),
    plot.margin      = margin(0, 0, 0, 0)
  ) +
  
  # Transects
  annotate("segment", x = -124.8, xend = -129.8, y = 44.6, yend = 44.6, 
           color = "#FDE047", linewidth = 0.8, linetype = "dashed", alpha = 0.7) +
  annotate("segment", x = -126.5, xend = -145.0, y = 48.6, yend = 50.0, 
           color = "#FFFFFF", linewidth = 0.8, linetype = "dashed", alpha = 0.7) +
  
  # Labels
  geom_text(
    data       = labels_df,
    aes(x = x, y = y, label = label),
    color      = "white",
    fontface   = labels_df$fontface,
    size       = labels_df$size,
    lineheight = 0.9,
    hjust      = 0.5
  )

# Canadian Maple Leaves
for (ml in maple_leaves) {
  p <- p + add_phylopic(
    img = icons$maple, x = ml$lon, y = ml$lat, height = ml$height,
    fill = "#B91C1C", color = "#7F1D1D", alpha = 0.85
  )
}

# Salmon Stages
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

# Forage Schools & Krill
for (fg in forage_schools) {
  p <- p + add_phylopic(
    img    = fg$img,
    x      = fg$lon,
    y      = fg$lat,
    height = fg$height,
    fill   = fg$fill,
    color  = fg$color,
    angle  = fg$angle
  )
}

# Northern Mammals
for (mm in marine_mammals) {
  p <- p + add_phylopic(
    img        = mm$img,
    x          = mm$lon,
    y          = mm$lat,
    height     = mm$height,
    fill       = mm$fill,
    color      = mm$color,
    angle      = mm$angle,
    horizontal = mm$flip
  )
}

# GOA Groundfish & Sharks
for (gp in goa_predators) {
  p <- p + add_phylopic(
    img        = gp$img,
    x          = gp$lon,
    y          = gp$lat,
    height     = gp$height,
    fill       = gp$fill,
    color      = gp$color,
    angle      = gp$angle,
    horizontal = gp$flip
  )
}

# DFO Line P Zooplankton
for (lp in line_p_zooplankton) {
  p <- p + add_phylopic(
    img    = lp$img,
    x      = lp$lon,
    y      = lp$lat,
    height = lp$height,
    fill   = lp$fill,
    color  = lp$color,
    angle  = lp$angle
  )
}

# Ferriss et al. 2025 Arenas
for (fa in ferriss_arenas) {
  p <- p + add_phylopic(
    img        = fa$img,
    x          = fa$lon,
    y          = fa$lat,
    height     = fa$height,
    fill       = fa$fill,
    color      = fa$color,
    angle      = fa$angle,
    horizontal = fa$flip
  )
}

# Newport Line Zooplankton
for (zp in nh_zooplankton) {
  p <- p + add_phylopic(
    img    = zp$img,
    x      = zp$lon,
    y      = zp$lat,
    height = zp$height,
    fill   = zp$fill,
    color  = zp$color,
    angle  = zp$angle
  )
}

# NCC Shelf Fauna
for (sf in ncc_shelf_fauna) {
  p <- p + add_phylopic(
    img        = sf$img,
    x          = sf$lon,
    y          = sf$lat,
    height     = sf$height,
    fill       = sf$fill,
    color      = sf$color,
    angle      = sf$angle,
    horizontal = sf$flip
  )
}

# Offshore Seabirds
for (sb in ncc_seabirds) {
  p <- p + add_phylopic(
    img        = sb$img,
    x          = sb$lon,
    y          = sb$lat,
    height     = sb$height,
    fill       = sb$fill,
    color      = sb$color,
    angle      = sb$angle,
    horizontal = sb$flip
  )
}

# Yaquina Head Murre with Prey in Bill
p <- p +
  add_phylopic(
    img = icons$murre, x = -124.3, y = 44.8, height = 0.80, horizontal = FALSE,
    fill = "#1E293B", color = "#0F172A"
  ) +
  add_phylopic(
    img = icons$forage_fish, x = -123.9, y = 45.0, height = 0.22, angle = -35,
    fill = "#38BDF8", color = "#0284C7"
  )

# Columbia River Pinnipeds
for (cr in columbia_pinnipeds) {
  p <- p + add_phylopic(
    img        = cr$img,
    x          = cr$lon,
    y          = cr$lat,
    height     = cr$height,
    fill       = cr$fill,
    color      = cr$color,
    angle      = cr$angle,
    horizontal = cr$flip
  )
}

# Map Extents
p <- p + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# ------------------------------------------------------------------------------
# RENDER & SAVE
# ------------------------------------------------------------------------------
print(p)

file_out <- file.path(fig_out_dir, "04_map_alldata_offline.v2.png")
ggsave(
  filename = file_out,
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

cat(sprintf("\n[Map Complete] Saved clean offline map to:\n%s\n", file_out))