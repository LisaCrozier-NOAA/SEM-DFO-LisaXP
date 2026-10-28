# ==============================================================================
# Script: 03_figs_tables/06_map_summary_model2_modulated_triad.R
# Purpose: Summary Slide Map B: Mechanistic Triad & Dual-Arena Ecosystem.
#          Displays the full Chinook migration loop, early coastal bottlenecks
#          (Hake + Bonneville SSL), and the 3 modulated candidate predators
#          (Sleeper Shark, Sablefish, SEAK SSL) paired with their respective
#          forage buffers (Capelin, Pollock/Herring, Sitka Herring) and
#          subarctic thermal gating (EGoA SST).
# Output: Rcode_for_paper/Routput_for_paper/plots/Summary_Map_Model2_Triad.png
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

# 1. Output directory setup
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
plot_out_dir <- file.path(proj_dir, "Routput_for_paper", "plots")
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Base vector landmass
land <- ne_countries(scale = "medium", returnclass = "sf")

# 3. Verified PhyloPic UUIDs
chinook_uuid   <- get_uuid(name = "Oncorhynchus tshawytscha")
ssl_uuid       <- "21db0060-c29a-4870-96a4-629d7424647f"
hake_uuid      <- "9e87f9e3-a1a6-49e0-9037-924d81ab53a3"
shark_uuid     <- get_uuid(name = "sharks")
sablefish_uuid <- get_uuid(name = "Anoplopoma fimbria")
pollock_uuid   <- get_uuid(name = "Gadus chalcogrammus")
herring_uuid   <- get_uuid(name = "Clupea pallasii")
capelin_uuid   <- get_uuid(name = "Mallotus villosus")

# ------------------------------------------------------------------------------
# 4. SALMON MIGRATION LOOP (Juvenile Smolts -> Immature Rearing -> Returning Adult)
# ------------------------------------------------------------------------------
salmon_loop <- tribble(
  ~stage, ~location,                 ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",                -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",           -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",               -132.2,   52.2,   150,        50,        30,   FALSE,
  4,      "SEAK Shelf",             -136.5,   56.2,   200,        70,        15,   FALSE,
  5,      "GOA Top Shelf",          -146.0,   58.2,   100,        90,         0,   FALSE,
  6,      "Immature Rearing W",     -154.5,   55.8,   330,       130,       -35,   FALSE,
  7,      "Immature Rearing S",     -146.5,   50.5,   370,       150,        25,   TRUE,
  8,      "Returning Adult Swing",  -133.0,   46.5,   450,       180,        12,   TRUE
) %>%
  mutate(icon_size = 0.58 + (length_mm - 100) * (0.80 / 350))

depth_palette <- colorRampPalette(c("#B2EBF2", "#00ACC1", "#0288D1", "#1565C0", "#0D47A1"))
stage_colors  <- depth_palette(100)[as.numeric(cut(salmon_loop$depth_m, breaks = 100))]

# ------------------------------------------------------------------------------
# 5. EARLY COASTAL BOTTLENECK ARENAS (Shared Across All 3 Models)
# ------------------------------------------------------------------------------
early_arenas <- tribble(
  ~taxa,          ~uuid,       ~lon,     ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  "Pacific Hake", hake_uuid,   -126.8,   45.8,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
  "Pacific Hake", hake_uuid,   -126.8,   46.1,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
  "Pacific Hake", hake_uuid,   -126.2,   45.8,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
  "Bonneville SSL", ssl_uuid,  -122.0,   45.7,    0.95,      0,    FALSE, "#1E293B", "#020617"
)

# ------------------------------------------------------------------------------
# 6. THE 3 MODULATED CANDIDATE PREDATORS (Unified Midnight Navy #08306B)
# ------------------------------------------------------------------------------
triad_predators <- tribble(
  ~candidate,       ~uuid,          ~lon,     ~lat,    ~height, ~angle, ~flip,
  # Candidate 1: Pacific Sleeper Shark (WGoA / Pelagic Slope)
  "Sleeper Shark",  shark_uuid,     -149.0,   53.5,    1.50,    -10,    FALSE,
  # Candidate 2: Alaska Sablefish 3+ (Deep GOA Slope / Subarctic Basin)
  "Sablefish",      sablefish_uuid, -143.0,   55.5,    1.10,    -25,    TRUE,
  # Candidate 3: Steller Sea Lion SEAK Pups (Coastal SEAK Haulout)
  "SEAK SSL",       ssl_uuid,       -136.2,   58.5,    1.45,    -10,    TRUE
)

# ------------------------------------------------------------------------------
# 7. MATCHED ALTERNATIVE PREY BUFFERS (Deep Forest Green #1B7837)
# ------------------------------------------------------------------------------
triad_prey <- tribble(
  ~target_pred,    ~uuid,        ~lon,     ~lat,    ~height, ~angle,
  # Capelin Buffer for Sleeper Shark (WGoA)
  "Sleeper Shark", capelin_uuid, -152.0,   56.2,    0.45,    -15,
  "Sleeper Shark", capelin_uuid, -151.3,   55.8,    0.40,     10,
  "Sleeper Shark", capelin_uuid, -152.6,   56.5,    0.42,     -5,
  
  # Pollock + Sitka Herring Buffer for Sablefish (EGoA / Slope)
  "Sablefish",     pollock_uuid, -145.2,   57.5,    0.70,     10,
  "Sablefish",     pollock_uuid, -144.2,   57.1,    0.70,     10,
  "Sablefish",     herring_uuid, -140.0,   57.2,    0.45,     15,
  "Sablefish",     herring_uuid, -141.0,   56.8,    0.45,     15,
  
  # Sitka Herring + Capelin Buffer for SEAK Sea Lions
  "SEAK SSL",      herring_uuid, -136.5,   57.3,    0.45,     20,
  "SEAK SSL",      herring_uuid, -136.1,   57.0,    0.40,     -5,
  "SEAK SSL",      capelin_uuid, -137.2,   56.8,    0.38,     10
)

# ------------------------------------------------------------------------------
# 8. ANNOTATIONS & LABELS
# ------------------------------------------------------------------------------
labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   62.0,     "bold",      5.2,
  "British\nColumbia",  -122.5,   53.5,     "plain",     4.5,
  "Columbia \nRiver",   -120.2,   45.8,     "plain",     4.2,
  "USA",                -120.5,   47.5,     "plain",     4.5,
  "Pacific Ocean",      -160.0,   46.0,     "plain",     5.5
)

# ------------------------------------------------------------------------------
# 9. ASSEMBLE MAP
# ------------------------------------------------------------------------------
p_map2 <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#6fa2bb", color = NA),
    plot.margin      = margin(0, 0, 0, 0)
  ) +
  
  # # Thermal Gating Zone (Summer EGoA SST contour footprint)
  # annotate(
  #   "rect", xmin = -153.5, xmax = -135.0, ymin = 51.5, ymax = 59.5,
  #   fill = "#C00000", alpha = 0.08, color = "#C00000", linetype = "dashed", linewidth = 0.7
  # ) +
  # annotate(
  #   "text", x = -144.0, y = 52.2,
  #   label = "Thermal Gating Arena (Summer EGoA SST)\nControls Predator Spatial-Thermal Overlap",
  #   color = "#990000", size = 3.6, fontface = "bold.italic", lineheight = 0.95
  # ) +
  # 
  # Geographic text labels
  geom_text(
    data       = labels_df,
    aes(x = x, y = y, label = label),
    color      = "white",
    fontface   = labels_df$fontface,
    size       = labels_df$size,
    lineheight = 0.9,
    hjust      = 0.5
  )

# Add Early Coastal Bottleneck Arena icons
for (i in seq_len(nrow(early_arenas))) {
  ea <- early_arenas[i, ]
  p_map2 <- p_map2 + add_phylopic(
    uuid       = ea$uuid,
    x          = ea$lon,
    y          = ea$lat,
    height     = ea$height,
    fill       = ea$fill,
    color      = ea$color,
    angle      = ea$angle,
    horizontal = ea$flip
  )
}

# Early Arena text callouts
p_map2 <- p_map2 +
  annotate("text", x = -128.5, y = 45.4, label = "Pacific Hake",
           size = 3.3, fontface = "bold", color = "white", lineheight = 0.9)# +
#   annotate("text", x = -120.0, y = 45.0, label = "Bonneville SSL",
#            size = 3.1, fontface = "italic", color = "white", lineheight = 0.9)
# # 
# Add Salmon Lifecycle Sequence
for (i in seq_len(nrow(salmon_loop))) {
  stg <- salmon_loop[i, ]
  p_map2 <- p_map2 + add_phylopic(
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

# Add Alternative Forage Buffers (Deep Green)
for (i in seq_len(nrow(triad_prey))) {
  pr <- triad_prey[i, ]
  p_map2 <- p_map2 + add_phylopic(
    uuid   = pr$uuid,
    x      = pr$lon,
    y      = pr$lat,
    height = pr$height,
    fill   = "#1B7837",
    color  = "#0F4A21",
    angle  = pr$angle
  )
}

# Add 3 Modulated Predators (Midnight Navy)
for (i in seq_len(nrow(triad_predators))) {
  tp <- triad_predators[i, ]
  p_map2 <- p_map2 + add_phylopic(
    uuid       = tp$uuid,
    x          = tp$lon,
    y          = tp$lat,
    height     = tp$height,
    fill       = "#08306B",
    color      = "#021226",
    angle      = tp$angle,
    horizontal = tp$flip
  )
}

# Predator & Prey identification callouts
p_map2 <- p_map2 +
  annotate("text", x = -149.0, y = 52.3, label = "Sleeper Shark\n+ Capelin Buffer",
           size = 3.2, fontface = "bold", color = "#08306B", lineheight = 0.9) +
  annotate("text", x = -141.5, y = 54.9, label = "Sablefish\n+ Pollock/Herring Buffer",
           size = 3.2, fontface = "bold", color = "#08306B", lineheight = 0.9) +
  annotate("text", x = -131.0, y = 59.4, label = "SEAK SSL\n+ Herring Buffer",
           size = 3.2, fontface = "bold", color = "#08306B", lineheight = 0.9)

# Crop bounds
p_map2 <- p_map2 + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# ------------------------------------------------------------------------------
# 10. SAVE OUTPUT (16:9 Slide Dimensions)
# ------------------------------------------------------------------------------
print(p_map2)

file_out <- file.path(plot_out_dir, "Summary_Map_Model2_Triad.png")
ggsave(
  filename = file_out,
  plot     = p_map2,
  width    = 13.33,
  height   = 7.5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

cat(sprintf("\n[Saved] Summary Slide Map B saved to: %s\n", file_out))