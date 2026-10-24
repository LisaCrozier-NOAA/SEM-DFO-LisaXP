library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)

proj_dir <- file.path(getwd(), "Rcode_for_paper")
fig_out_dir <- file.path(proj_dir, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

land <- ne_countries(scale = "medium", returnclass = "sf")

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

base_map <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  
  coord_sf(xlim = c(-170, -119), ylim = c(43, 62), expand = FALSE) +
  
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#1d5873", color = NA),
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

print(base_map)

ggsave(
  filename = file.path(fig_out_dir, "map_layer2.png"),
  plot = base_map,
  width = 10,
  height = 7.5,
  dpi = 300
)