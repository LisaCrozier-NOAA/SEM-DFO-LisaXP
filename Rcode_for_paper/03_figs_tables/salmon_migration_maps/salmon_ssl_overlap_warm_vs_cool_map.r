# ==============================================================================
# Script: salmon_ssl_overlap_warm_vs_cool_map.r
# Purpose: Replicate Santora et al. (2020) Fig 1c style map comparing 2013 (Cool)
#          vs 2014 (Warm) SST habitat compression off WA/OR coast with PhyloPic 
#          predator and salmon silhouettes.
# Output: salmon_habitat_warm_vs_cool_map.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(rphylopic) # Icon library for silhouettes
})
pick_phylopic(name="steller sea lion")
pick_phylopic(name="chinook salmon")

# --- 1. Fetch PhyloPic UUIDs for Sea Lion and Salmon ---
cat("Fetching animal silhouettes from PhyloPic database...\n")
salmon_uuid  <- get_uuid(name = "chinook salmon")
sealion_uuid <- get_uuid(name = "steller sea lion", n = 1)

# Fetch SVG objects
img_salmon  <- get_phylopic(uuid = salmon_uuid)
img_sealion <- get_phylopic(uuid = sealion_uuid)

# --- 2. Generate Synthetic Coastal SST Grid (WA/OR Coast: 42N to 48.5N) ---
grid_lat <- seq(42, 48.5, by = 0.1)
grid_lon <- seq(-126, -123.5, by = 0.1)
map_grid <- expand.grid(lon = grid_lon, lat = grid_lat)

map_grid <- map_grid %>%
  mutate(dist_from_coast = abs(lon - (-124.0)))

# 2013 (Cool Year): Cold upwelling band (< 12°C) extends far offshore
map_grid_2013 <- map_grid %>%
  mutate(
    Year = "2013 (Cool / Expanded Shelf)",
    SST = 8 + 3 * dist_from_coast + rnorm(n(), 0, 0.2)
  )

# 2014 (Warm Year): Cold upwelling compressed into tight nearshore strip
map_grid_2014 <- map_grid %>%
  mutate(
    Year = "2014 (Warm / Compressed Shelf)",
    SST = 11.5 + 6 * dist_from_coast + rnorm(n(), 0, 0.2)
  )

full_grid <- rbind(map_grid_2013, map_grid_2014)

# --- 3. Define Animal Coordinates (5 Sea Lions, 5 Salmon) ---
# 2013: Sea Lions spread wide across the expansive blue cool-water shelf
csl_2013 <- data.frame(
  lon = c(-125.4, -124.9, -124.6, -125.6, -124.4),
  lat = c(43.5, 45.0, 46.2, 47.1, 48.0),
  Type = "Sea Lion"
)

# 2014: Sea Lions compressed nearshore tracking the shrunken blue shelf
csl_2014 <- data.frame(
  lon = c(-124.3, -124.2, -124.35, -124.25, -124.15),
  lat = c(43.5, 45.0, 46.2, 47.1, 48.0),
  Type = "Sea Lion"
)

# Salmon: Always restricted nearshore (Columbia River plume / coastal migration)
salmon_base <- data.frame(
  lon = c(-124.18, -124.12, -124.22, -124.08, -124.15),
  lat = c(44.0, 45.2, 46.1, 47.0, 47.8),
  Type = "Spring Chinook"
)

# --- 4. Plot Builder Function using rphylopic ---
build_map_panel <- function(df_grid, df_csl, df_salmon, title_text) {
  
  p <- ggplot() +
    # SST Heatmap raster
    geom_raster(data = df_grid, aes(x = lon, y = lat, fill = SST)) +
    scale_fill_gradientn(
      colors = c("#1E3A8A", "#3B82F6", "#93C5FD", "#FDE047", "#EF4444"),
      limits = c(8, 18),
      name = "SST (°C)"
    ) +
    # 12°C Isotherm / Habitat Compression Boundary
    geom_contour(data = df_grid, aes(x = lon, y = lat, z = SST), breaks = 12, color = "white", size = 1.2) +
    coord_quickmap(xlim = c(-126, -123.5), ylim = c(42, 48.5)) +
    labs(title = title_text, x = "Longitude (°W)", y = "Latitude (°N)") +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "gray90"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # Add Sea Lion PhyloPic Silhouettes
  for (i in 1:nrow(df_csl)) {
    p <- p + add_phylopic(
      img = img_sealion,
      x = df_csl$lon[i],
      y = df_csl$lat[i],
      height = 0.35,
      fill = "#F97316" # Orange silhouette
    )
  }
  
  # Add Salmon PhyloPic Silhouettes
  for (i in 1:nrow(df_salmon)) {
    p <- p + add_phylopic(
      img = img_salmon,
      x = df_salmon$lon[i],
      y = df_salmon$lat[i],
      height = 0.25,
      fill = "#10B981" # Green silhouette
    )
  }
  
  return(p)
}

# --- 5. Generate Panels & Combine ---
p_2013 <- build_map_panel(map_grid_2013, csl_2013, salmon_base, "2013 (Cool / Expanded Shelf)")
p_2014 <- build_map_panel(map_grid_2014, csl_2014, salmon_base, "2014 (Warm / Compressed Shelf)")

final_map <- p_2013 + p_2014 + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Thermal Habitat Compression Driven Spatial Overlap",
    subtitle = "White line = 12°C Isotherm. Green = Spring Chinook smolts; Orange = Sea Lions (PhyloPic icons).\nIn 2014 (warm), predators compress into the narrow nearshore salmon zone.",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

final_map

ggsave("salmon_habitat_warm_vs_cool_map.png", final_map, width = 12, height = 8, dpi = 300)
cat("Habitat map successfully saved to salmon_habitat_warm_vs_cool_map.png\n")