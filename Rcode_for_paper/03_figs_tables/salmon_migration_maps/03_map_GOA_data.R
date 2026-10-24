# ==============================================================================
# Script: 03_figs_tables/03_map_GOA_data.R
# Purpose: Comprehensive GOA & Coastal Map showing 1st-year smolt outmigration,
#          returning adult, forage schools (herring, capelin, krill), coastal
#          marine mammals (SSL, seals, orcas), and the Gulf of Alaska demersal
#          and pelagic predator guild (pollock, arrowtooth, sablefish, halibut,
#          Pacific cod, sharks).
# Output: figures/02_map_DFO_data.png
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

# Existing taxa
chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
ssl_uuid     <- get_uuid(name = "Zalophus californianus wollebaeki") # Sea lion
herring_uuid <- get_uuid(name = "Clupea pallasii")                   # Pacific herring
capelin_uuid <- get_uuid(name = "Mallotus villosus")                 # Capelin
orca_uuid    <- get_uuid(name = "Orcinus orca")                      # Killer whale
seal_uuid    <- get_uuid(name = "Phoca vitulina")                    # Harbour seal
krill_uuid   <- get_uuid(name = "Euphausia pacifica")                # Krill
maple_uuid   <- get_uuid(name = "Acer")                              # Maple leaf

# Added GOA Groundfish and Shark taxa
pollock_uuid    <- get_uuid(name = "Gadus chalcogrammus")            # Walleye pollock
flounder_uuid   <- get_uuid(name = "righteye flounders")            # Arrowtooth flounder
#flounder_uuid   <- get_uuid(name = "Atheresthes stomias")            # Arrowtooth flounder
sablefish_uuid  <- get_uuid(name = "Anoplopoma fimbria")             # Sablefish / Blackcod
halibut_uuid    <- get_uuid(name = "halibuts")        # Pacific halibut
#halibut_uuid    <- get_uuid(name = "Hippoglossus stenolepis")        # Pacific halibut
cod_uuid        <- get_uuid(name = "cods")            # Pacific cod
#cod_uuid        <- get_uuid(name = "Gadus macrocephalus")            # Pacific cod
shark_uuid      <- get_uuid(name = "sharks")                         # Pelagic / sleeper shark

# 3. Salmon Migration Stages: 1st Year (Stages 1-3) + Returning Adult (Stage 8)
salmon_stages <- tribble(
  ~stage, ~location,                 ~lon,    ~lat,  ~length_mm, ~depth_m, ~angle, ~flip,
  1,      "CR Exit",                -125.5,   46.8,   100,        10,        40,   FALSE,
  2,      "WCVI Coastal",           -128.2,   49.5,   150,        30,        30,   FALSE,
  3,      "BC Coast",               -132.2,   52.2,   150,        50,        30,   FALSE,
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
  
  # WCVI Krill Cluster
  "Krill",   krill_uuid,   -127.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -130.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -132.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -131.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444"
)

# 5. Marine Mammals
marine_mammals <- tribble(
  ~group,         ~uuid,      ~lon,    ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  # Steller Sea Lions (Haulouts)
  "SSL",          ssl_uuid,   -128.2,  50.5,    1.50,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -136.1,  58.3,    1.35,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -146.3,  60.0,    1.30,      0,    FALSE, "#372E29", "#1F1A17",
  
  # Harbour Seals
  "Harbour Seal", seal_uuid,  -125.4,  49.8,    0.85,     15,    FALSE, "#4B5563", "#1F2937",
  "Harbour Seal", seal_uuid,  -130.6,  53.0,    0.80,    -10,    TRUE,  "#4B5563", "#1F2937",
  
  # Killer Whales
  "Orca",         orca_uuid,  -133.8,  51.6,    2.25,     25,    FALSE, "#111827", "#000000",
  "Orca",         orca_uuid,  -136.5,  53.7,    2.15,     30,    FALSE, "#111827", "#000000"
)

# 6. Gulf of Alaska Groundfish, Flatfish, and Pelagic Sharks
goa_predators <- tribble(
  ~species,              ~uuid,          ~lon,    ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  # Walleye Pollock (Mid-shelf GOA)
  "Walleye Pollock",     pollock_uuid,   -149.5,  58.0,    0.75,     10,    TRUE,  "#64748B", "#334155",
  "Walleye Pollock",     pollock_uuid,   -148.0,  58.4,    0.65,    -15,    TRUE,  "#64748B", "#334155",
  "Walleye Pollock",     pollock_uuid,   -147.0,  57.6,    0.65,    -15,    TRUE,  "#64748B", "#334155",
  
  # Pacific Cod (Kodiak / Western GOA Shelf)
  "Pacific Cod",         cod_uuid,       -153.5,  54.8,    0.85,     20,    FALSE, "red", "#27272A",
#  "Pacific Cod",         cod_uuid,       -153.5,  56.8,    0.85,     20,    FALSE, "#52525B", "#27272A",
  
  # Pacific Halibut (Outer continental shelf break)
  "Pacific Halibut",     halibut_uuid,   -147.2,  55.5,    1.05,    -25,    FALSE, "green", "#1C1917",
#  "Pacific Halibut",     halibut_uuid,   -147.2,  55.5,    1.05,    -25,    FALSE, "#44403C", "#1C1917",
  
  # Arrowtooth Flounder (GOA Shelf Slope)
  "Arrowtooth Flounder", flounder_uuid,  -143.5,  57.2,    0.80,     15,    TRUE,  "yellow", "#292524",
"Arrowtooth Flounder", flounder_uuid,  -143.5,  67.2,    0.80,     15,    TRUE,  "yellow", "#292524",
#  "Arrowtooth Flounder", flounder_uuid,  -143.5,  57.2,    0.80,     15,    TRUE,  "#78716C", "#292524",
  
  # Sablefish (Deep bathyal slope / Eastern GOA)
  "Sablefish",           sablefish_uuid, -143.0,  55.8,    0.85,    -30,    TRUE,  "#2B547E", "#1E3A5F",
"Sablefish",           sablefish_uuid, -141.0,  55.8,    0.85,    -30,    TRUE,  "#2B547E", "#1E3A5F",

  # Sharks (Pacific Sleeper & Salmon Shark - Offshore Pelagic/Slope)
  "Sleeper Shark",       shark_uuid,     -145.0,  52.8,    1.30,    -10,    FALSE, "#1E293B", "#0F172A",
  "Salmon Shark",        shark_uuid,     -151.0,  52.0,    1.15,     20,    TRUE,  "#334155", "#0F172A"
)

# 7. Maple Leaf Badges on BC
maple_leaves <- tribble(
  ~lon,   ~lat,  ~height,
  -124.5, 51.2,  1.20,
  -126.5, 53.2,  1.20,
  -132.0, 54.0,  1.20
)

# Geographic Labels
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

# 8. Build Map Layer
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

# Add Maple Leaves on BC
for (i in seq_len(nrow(maple_leaves))) {
  ml <- maple_leaves[i, ]
  p <- p + add_phylopic(
    uuid   = maple_uuid,
    x      = ml$lon,
    y      = ml$lat,
    height = ml$height,
    fill   = "#B91C1C",
    color  = "#7F1D1D",
    alpha  = 0.85
  )
}

# Add GOA Groundfish and Sharks
for (i in seq_len(nrow(goa_predators))) {
  gp <- goa_predators[i, ]
  p <- p + add_phylopic(
    uuid       = gp$uuid,
    x          = gp$lon,
    y          = gp$lat,
    height     = gp$height,
    fill       = gp$fill,
    color      = gp$color,
    angle      = gp$angle,
    horizontal = gp$flip
  )
}

# Add Salmon Silhouettes
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

# 9. Save Output Image
print(p)

ggsave(
  filename = file.path(fig_out_dir, "03_map_GOA_data.png"),
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

message("Success! Saved comprehensive map including GOA groundfish, halibut, cod, and sharks.")