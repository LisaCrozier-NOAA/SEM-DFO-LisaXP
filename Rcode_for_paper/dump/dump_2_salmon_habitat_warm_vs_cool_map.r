# ==============================================================================
# Script: salmon_habitat_warm_vs_cool_map.r
# Purpose: Replicate Santora et al. (2020) Fig 1c style map comparing 2013 (Cool)
#          vs 2014 (Warm) SST habitat compression off WA/OR coast with predator 
#          and salmon distribution icons.
# Output: salmon_habitat_warm_vs_cool_map.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(patchwork)
})

# --- 1. Generate Synthetic Coastal SST Grid (WA/OR Coast: 42N to 48.5N) ---
grid_lat <- seq(42, 48.5, by = 0.1)
grid_lon <- seq(-126, -123.5, by = 0.1)
map_grid <- expand.grid(lon = grid_lon, lat = grid_lat)

# Calculate Distance from Coast (Approximate Coastline at lon = -124.0)
map_grid <- map_grid %>%
  mutate(dist_from_coast = abs(lon - (-124.0)))

# 2013 (Cool Year): Cold upwelling band (< 12°C) extends far offshore (~1.5 degrees lon)
map_grid_2013 <- map_grid %>%
  mutate(
    Year = "2013 (Cool / Expanded Shelf)",
    SST = 8 + 3 * dist_from_coast + rnorm(n(), 0, 0.3),
    Cool_Habitat = ifelse(SST <= 12, "Cool Upwelling (<12°C)", "Warm Offshore (>12°C)")
  )

# 2014 (Warm Year / Blob): Cold upwelling compressed into tight nearshore strip (<0.4 degrees lon)
map_grid_2014 <- map_grid %>%
  mutate(
    Year = "2014 (Warm / Compressed Shelf)",
    SST = 11.5 + 6 * dist_from_coast + rnorm(n(), 0, 0.3),
    Cool_Habitat = ifelse(SST <= 12, "Cool Upwelling (<12°C)", "Warm Offshore (>12°C)")
  )

full_grid <- rbind(map_grid_2013, map_grid_2014)

# --- 2. Fixed Animals (5 Sea Lions, 5 Salmon) ---
# Sea Lions (CSL): Spread across the ENTIRE blue cool-water zone in each year
csl_2013 <- data.frame(
  lon = c(-125.2, -124.8, -124.5, -125.5, -124.2),
  lat = c(43.5, 45.0, 46.2, 47.1, 48.0),
  Type = "Sea Lion", Year = "2013 (Cool / Expanded Shelf)"
)

# In 2014, blue area shrinks nearshore -> Sea Lions COMPRESSED nearshore
csl_2014 <- data.frame(
  lon = c(-124.3, -124.2, -124.35, -124.25, -124.15),
  lat = c(43.5, 45.0, 46.2, 47.1, 48.0),
  Type = "Sea Lion", Year = "2014 (Warm / Compressed Shelf)"
)

# Salmon: Always restricted nearshore (Columbia River plume / coastal migration)
salmon_fixed <- data.frame(
  lon = c(-124.2, -124.15, -124.25, -124.1, -124.2),
  lat = c(44.0, 45.2, 46.1, 47.0, 47.8),
  Type = "Spring Chinook"
)

salmon_2013 <- salmon_fixed %>% mutate(Year = "2013 (Cool / Expanded Shelf)")
salmon_2014 <- salmon_fixed %>% mutate(Year = "2014 (Warm / Compressed Shelf)")

animals_df <- rbind(csl_2013, csl_2014, salmon_2013, salmon_2014)

# --- 3. Plot Map Panels ---
plot_habitat_map <- function(year_label) {
  sub_grid <- full_grid %>% filter(Year == year_label)
  sub_anim <- animals_df %>% filter(Year == year_label)
  
  ggplot() +
    # SST Heatmap raster
    geom_raster(data = sub_grid, aes(x = lon, y = lat, fill = SST)) +
    scale_fill_gradientn(
      colors = c("#1E3A8A", "#3B82F6", "#93C5FD", "#FDE047", "#EF4444"),
      limits = c(8, 18),
      name = "SST (°C)"
    ) +
    # 12°C Isotherm / Habitat Compression Boundary Line
    geom_contour(data = sub_grid, aes(x = lon, y = lat, z = SST), breaks = 12, color = "white", size = 1.2, linetype = "solid") +
    # Animal Overlays
    geom_point(data = sub_anim, aes(x = lon, y = lat, shape = Type, color = Type), size = 5, stroke = 1.5) +
    scale_shape_manual(values = c("Sea Lion" = 17, "Spring Chinook" = 19)) +
    scale_color_manual(values = c("Sea Lion" = "#F97316", "Spring Chinook" = "#10B981")) +
    coord_quickmap(xlim = c(-126, -123.5), ylim = c(42, 48.5)) +
    labs(
      title = year_label,
      x = "Longitude (°W)",
      y = "Latitude (°N)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "gray90"),
      plot.title = element_text(face = "bold", size = 14)
    )
}

p_2013 <- plot_habitat_map("2013 (Cool / Expanded Shelf)")
p_2014 <- plot_habitat_map("2014 (Warm / Compressed Shelf)")

# Combine Both Panels Side-by-Side
final_map <- p_2013 + p_2014 + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Thermal Habitat Compression Driven Spatial Overlap",
    subtitle = "White line = 12°C Isotherm. In 2014 (warm), predators (triangles) compress into the nearshore salmon zone (circles).",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("salmon_habitat_warm_vs_cool_map.png", final_map, width = 12, height = 8, dpi = 300)
cat("Habitat map successfully saved to salmon_habitat_warm_vs_cool_map.png\n")