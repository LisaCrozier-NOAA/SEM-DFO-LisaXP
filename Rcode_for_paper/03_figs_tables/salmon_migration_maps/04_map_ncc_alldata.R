# ==============================================================================
# Script: 03_figs_tables/04_map_alldata.R
# Purpose: Ocean-basin scale map restoring exact Canadian layout
#          (WCVI krill cluster, Canadian maple leaves, mammals, herring)
#          while keeping spaced-out NCC fauna and GOA layout.
# Output: 03_figs_tables/salmon_migration_maps/figs/04_map_alldata_spaced.png
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

# Taxa UUIDs (verified)
chinook_uuid   <- get_uuid(name = "Oncorhynchus tshawytscha")
ssl_uuid       <- get_uuid(name = "Zalophus californianus wollebaeki")
seal_uuid      <- get_uuid(name = "Phoca vitulina")
herring_uuid   <- get_uuid(name = "Clupea pallasii")
capelin_uuid   <- get_uuid(name = "Mallotus villosus")
orca_uuid      <- get_uuid(name = "Orcinus orca")
krill_uuid     <- get_uuid(name = "Euphausia pacifica")
maple_uuid     <- get_uuid(name = "Acer")

pollock_uuid   <- get_uuid(name = "Gadus chalcogrammus")
flounder_uuid  <- get_uuid(name = "righteye flounders")
sablefish_uuid <- get_uuid(name = "Anoplopoma fimbria")
halibut_uuid   <- get_uuid(name = "halibuts")
cod_uuid       <- get_uuid(name = "cods")
shark_uuid     <- get_uuid(name = "sharks")

copepod_uuid   <- "be9c4ede-4f33-4800-a202-4c73f493b5e6"
crab_uuid      <- "035b0041-ec5b-4a9c-ac8e-5a3cc3ee02a1"
larval_uuid    <- "7360374f-8860-44c5-b74b-af53ca1a8af4"
sea_nettle_uuid<- "ebed2a7d-0291-41bf-97be-c961fb28f3c9"
squid_uuid     <- "68922290-f613-4409-af5d-d7fefbb055b1"
murre_uuid     <- "822f0866-55b5-4b8b-aa0e-df6cf2678865"
shearw_uuid    <- "a8ecc1e9-f339-4d85-84c3-aa2fff47b2c1"
pelican_uuid   <- "2a168f51-4016-4a32-ba72-64a53bee31ca"
hake_uuid      <- "9e87f9e3-a1a6-49e0-9037-924d81ab53a3"
rockfish_uuid  <- "590dd989-103a-4b22-8788-81471386d860"
forage_uuid    <- "5d05c64d-9232-45a4-84f5-6f6519c4a6fb"

# ------------------------------------------------------------------------------
# 3. SALMON STAGES (Smolts 1-3 + Returning Adult 8 at -130.5, 45.2)
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
# 4. EXACT CANADIAN & NORTHERN FORAGE SCHOOLS (Including 6-krill WCVI cluster)
# ------------------------------------------------------------------------------
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
  
  # WCVI Krill Cluster (Full 6 krill cluster)
  "Krill",   krill_uuid,   -127.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -128.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -130.8,  48.7,    0.32,     35,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -132.3,  48.5,    0.28,     15,    "#FCA5A5", "#EF4444",
  "Krill",   krill_uuid,   -131.0,  48.9,    0.30,     50,    "#FCA5A5", "#EF4444"
)

# ------------------------------------------------------------------------------
# 5. EXACT CANADIAN & NORTHERN MARINE MAMMALS
# ------------------------------------------------------------------------------
marine_mammals <- tribble(
  ~group,         ~uuid,      ~lon,    ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  # Steller Sea Lions
  "SSL",          ssl_uuid,   -128.2,  50.5,    1.50,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -136.1,  58.3,    1.35,    -10,    TRUE,  "#372E29", "#1F1A17",
  "SSL",          ssl_uuid,   -146.3,  60.0,    1.30,      0,    FALSE, "#372E29", "#1F1A17",
  
  # Harbour Seals
  "Harbour Seal", seal_uuid,  -125.4,  49.8,    0.85,     15,    FALSE, "#4B5563", "#1F2937",
  "Harbour Seal", seal_uuid,  -130.6,  53.0,    0.80,    -10,    TRUE,  "#4B5563", "#1F2937",
  
  # Killer Whales (Larger sizes: 2.25 & 2.15)
  "Orca",         orca_uuid,  -133.8,  51.6,    2.25,     25,    FALSE, "#111827", "#000000",
  "Orca",         orca_uuid,  -136.5,  53.7,    2.15,     30,    FALSE, "#111827", "#000000"
)

# ------------------------------------------------------------------------------
# 6. EXACT CANADIAN MAPLE LEAVES (BC Coast & Haida Gwaii/Vancouver Island)
# ------------------------------------------------------------------------------
maple_leaves <- tribble(
  ~lon,   ~lat,  ~height,
  -124.5, 51.2,  1.20,
  -126.5, 53.2,  1.20,
  -132.0, 54.0,  1.20
)

# ------------------------------------------------------------------------------
# 7. GULF OF ALASKA PREDATORS
# ------------------------------------------------------------------------------
goa_predators <- tribble(
  ~species,              ~uuid,          ~lon,    ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  "Walleye Pollock",     pollock_uuid,   -149.5,  58.0,    0.75,     10,    TRUE,  "#64748B", "#334155",
  "Walleye Pollock",     pollock_uuid,   -148.0,  58.4,    0.65,    -15,    TRUE,  "#64748B", "#334155",
  "Walleye Pollock",     pollock_uuid,   -147.0,  57.6,    0.65,    -15,    TRUE,  "#64748B", "#334155",
  "Pacific Cod",         cod_uuid,       -153.5,  54.8,    0.85,     20,    FALSE, "#52525B", "#27272A",
  "Pacific Halibut",     halibut_uuid,   -147.2,  55.5,    1.05,    -25,    FALSE, "#44403C", "#1C1917",
  "Arrowtooth Flounder", flounder_uuid,  -143.5,  57.2,    0.80,     15,    TRUE,  "#78716C", "#292524",
  "Sablefish",           sablefish_uuid, -143.0,  55.8,    0.85,    -30,    TRUE,  "#2B547E", "#1E3A5F",
  "Sablefish",           sablefish_uuid, -141.0,  55.8,    0.85,    -30,    TRUE,  "#2B547E", "#1E3A5F",
  "Sleeper Shark",       shark_uuid,     -145.0,  52.8,    1.30,    -10,    FALSE, "#1E293B", "#0F172A",
  "Salmon Shark",        shark_uuid,     -151.0,  52.0,    1.15,     20,    TRUE,  "#334155", "#0F172A"
)

# ------------------------------------------------------------------------------
# 8. SPATIALLY EXPANDED NCC SUITE (Offshore & Fanned Out)
# ------------------------------------------------------------------------------
nh_zooplankton <- tribble(
  ~taxa,          ~uuid,           ~lon,     ~lat,    ~height, ~angle, ~fill,     ~color,
  "Copepod",      copepod_uuid,    -125.8,   44.6,    0.45,     20,    "#FDE047", "#CA8A04",
  "Crab Megalopa",crab_uuid,       -127.0,   44.6,    0.45,      0,    "#FB923C", "#C2410C",
  "Larval Fish",  larval_uuid,     -128.2,   44.6,    0.45,     15,    "#FDE68A", "#D97706",
  "Krill",        krill_uuid,      -129.4,   44.6,    0.42,     30,    "#F87171", "#DC2626"
)

ncc_shelf_fauna <- tribble(
  ~taxa,          ~uuid,           ~lon,     ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  "Pacific Hake", hake_uuid,       -127.2,   47.0,    0.85,      20,   FALSE, "#38BDF8", "#0284C7",
  "Rockfish",     rockfish_uuid,   -126.2,   46.2,    0.70,     -15,   FALSE, "#FB7185", "#E11D48",
  "Market Squid", squid_uuid,      -128.2,   46.1,    0.65,      35,   FALSE, "#C084FC", "#9333EA",
  "Sea Nettle",   sea_nettle_uuid, -127.2,   45.3,    0.80,     -10,   FALSE, "#F472B6", "#DB2777",
  "Forage Fish",  forage_uuid,     -126.0,   47.6,    0.40,      15,   FALSE, "#94A3B8", "#475569",
  "Forage Fish",  forage_uuid,     -126.4,   47.4,    0.35,     -10,   FALSE, "#94A3B8", "#475569"
)

ncc_seabirds <- tribble(
  ~taxa,              ~uuid,        ~lon,     ~lat,    ~height, ~angle, ~fill,     ~color,
  "Common Murre",     murre_uuid,   -128.2,   47.8,    0.75,     10,    "#1E293B", "#0F172A",
  "Sooty Shearwater", shearw_uuid,  -129.5,   46.8,    0.75,     30,    "#44403C", "#1C1917",
  "Brown Pelican",    pelican_uuid, -128.8,   45.5,    0.95,     15,    "#78716C", "#292524",
  "Common Murre",     murre_uuid,   -127.8,   44.0,    0.75,      5,    "#1E293B", "#0F172A"
)

yaquina_murre <- tibble(lon = -124.3, lat = 44.8, height = 0.80)
yaquina_fish  <- tibble(lon = -123.9, lat = 45.0, height = 0.22, angle = -35)

columbia_pinnipeds <- tribble(
  ~taxa,          ~uuid,       ~lon,     ~lat,    ~height, ~angle, ~flip, ~fill,     ~color,
  "Harbor Seal",  seal_uuid,   -124.0,   46.3,    0.60,     10,    FALSE, "#475569", "#1E293B",
  "Sea Lion",     ssl_uuid,    -123.0,   46.0,    0.80,     -5,    FALSE, "#334155", "#0F172A",
  "Sea Lion",     ssl_uuid,    -122.0,   45.7,    0.90,      0,    FALSE, "#1E293B", "#020617"
)

# ------------------------------------------------------------------------------
# 9. ASSEMBLE MAP
# ------------------------------------------------------------------------------
labels_df <- tribble(
  ~label,                ~x,       ~y,       ~fontface,   ~size,
  "Aleutian \nIslands", -167.0,   55.5,     "italic",    4.5,
  "Alaska",             -143.0,   61.5,     "bold",      5.2,
  "British\nColumbia",  -122.5,   53.5,     "plain",     4.5,
  "WCVI",               -125.8,   48.8,     "plain",     4.2,
  "Columbia \nRiver",   -120.2,   45.8,     "plain",     4.2,
  "USA",                -120.5,   47.5,     "plain",     4.5,
  "Pacific Ocean",      -158.0,   45.5,     "plain",     5.5,
  "NH Line",            -130.6,   44.6,     "italic",    3.6
)

p <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#6fa2bb", color = NA),
    plot.margin      = margin(0, 0, 0, 0)
  ) +
  
  # Dashed NH Line across the zooplankton transect
  annotate("segment", x = -124.8, xend = -129.8, y = 44.6, yend = 44.6, 
           color = "#FDE047", linewidth = 0.8, linetype = "dashed", alpha = 0.7) +
  
  # Geographic text
  geom_text(
    data       = labels_df,
    aes(x = x, y = y, label = label),
    color      = "white",
    fontface   = labels_df$fontface,
    size       = labels_df$size,
    lineheight = 0.9,
    hjust      = 0.5
  )

# Add All 3 Canadian Maple Leaves
for (i in seq_len(nrow(maple_leaves))) {
  p <- p + add_phylopic(
    uuid   = maple_uuid,
    x      = maple_leaves$lon[i],
    y      = maple_leaves$lat[i],
    height = maple_leaves$height[i],
    fill   = "#B91C1C",
    color  = "#7F1D1D",
    alpha  = 0.85
  )
}

# Add Forage Species (Herring, Capelin, full 6-krill WCVI cluster)
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

# Add Marine Mammals (SSL, Harbour Seals, Orcas at 2.25/2.15 size)
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

# Add Newport Hydrographic Line Zooplankton
for (i in seq_len(nrow(nh_zooplankton))) {
  zp <- nh_zooplankton[i, ]
  p <- p + add_phylopic(
    uuid   = zp$uuid,
    x      = zp$lon,
    y      = zp$lat,
    height = zp$height,
    fill   = zp$fill,
    color  = zp$color,
    angle  = zp$angle
  )
}

# Add Columbia River Pinnipeds
for (i in seq_len(nrow(columbia_pinnipeds))) {
  cr <- columbia_pinnipeds[i, ]
  p <- p + add_phylopic(
    uuid       = cr$uuid,
    x          = cr$lon,
    y          = cr$lat,
    height     = cr$height,
    fill       = cr$fill,
    color      = cr$color,
    angle      = cr$angle,
    horizontal = cr$flip
  )
}

# Add NCC Shelf Pelagics & Invertebrates
for (i in seq_len(nrow(ncc_shelf_fauna))) {
  sf <- ncc_shelf_fauna[i, ]
  p <- p + add_phylopic(
    uuid       = sf$uuid,
    x          = sf$lon,
    y          = sf$lat,
    height     = sf$height,
    fill       = sf$fill,
    color      = sf$color,
    angle      = sf$angle,
    horizontal = sf$flip
  )
}

# Add Offshore Seabirds
for (i in seq_len(nrow(ncc_seabirds))) {
  sb <- ncc_seabirds[i, ]
  p <- p + add_phylopic(
    uuid   = sb$uuid,
    x      = sb$lon,
    y      = sb$lat,
    height = sb$height,
    fill   = sb$fill,
    color  = sb$color,
    angle  = sb$angle
  )
}

# Add Murre carrying Forage Fish in Bill
p <- p + 
  add_phylopic(
    uuid = murre_uuid, x = yaquina_murre$lon, y = yaquina_murre$lat,
    height = yaquina_murre$height, fill = "#1E293B", color = "#0F172A"
  ) +
  add_phylopic(
    uuid = forage_uuid, x = yaquina_fish$lon, y = yaquina_fish$lat,
    height = yaquina_fish$height, angle = yaquina_fish$angle,
    fill = "#38BDF8", color = "#0284C7"
  )

# Set map limits
p <- p + coord_sf(
  xlim   = c(-170, -119),
  ylim   = c(43, 62),
  expand = FALSE
)

# ------------------------------------------------------------------------------
# 10. SAVE OUTPUT
# ------------------------------------------------------------------------------
print(p)

ggsave(
  filename = file.path(fig_out_dir, "04_map_alldata.png"),
  plot     = p,
  width    = 10,
  height   = 7.5,
  dpi      = 300
)

message("Success! Saved map with exact Canadian configuration and spaced NCC suite.")