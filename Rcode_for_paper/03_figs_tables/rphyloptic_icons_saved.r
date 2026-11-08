# ==============================================================================
# Script: rphylopic_icons_saved.R
# Purpose: Download and save all project PhyloPic icons locally for offline use.
# ==============================================================================

# Step 2: Use the Local Icons in Any Script Without Internet
# In your mapping and plotting scripts, replace live get_uuid() calls with readRDS(), and pass img = icons$name instead of uuid = ... into add_phylopic():

library(rphylopic)

cache_dir <- file.path(getwd(), "Rcode_for_paper", "03_figs_tables", "phylopic_cache")
# cache_dir <- file.path(getwd(), "Rcode_for_paper", "assets", "phylopic_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Master UUID list used across your maps and stepped plots
icon_uuids <- list(
  chinook       = get_uuid(name = "Oncorhynchus tshawytscha"),
  ssl           = "21db0060-c29a-4870-96a4-629d7424647f",
  seal          = "1925a9aa-bf33-4270-824c-b8f781f8ec82",
  herring       = get_uuid(name = "Clupea pallasii"),
  capelin       = get_uuid(name = "Mallotus villosus"),
  orca          = get_uuid(name = "Orcinus orca"),
  krill         = "44a3628d-aafd-45cc-97a6-1cb74bd43dec",
  maple         = get_uuid(name = "Acer"),
  pollock       = get_uuid(name = "Gadus chalcogrammus"),
  flounder      = get_uuid(name = "righteye flounders"),
  sablefish     = get_uuid(name = "Anoplopoma fimbria"),
  halibut       = get_uuid(name = "halibuts"),
  cod           = get_uuid(name = "cods"),
  shark         = get_uuid(name = "sharks"),
  copepod       = "be9c4ede-4f33-4800-a202-4c73f493b5e6",
  crab          = "035b0041-ec5b-4a9c-ac8e-5a3cc3ee02a1",
  larval_fish   = "7360374f-8860-44c5-b74b-af53ca1a8af4",
  sea_nettle    = "ebed2a7d-0291-41bf-97be-c961fb28f3c9",
  squid         = "68922290-f613-4409-af5d-d7fefbb055b1",
  murre         = "822f0866-55b5-4b8b-aa0e-df6cf2678865",
  shearwater    = "a8ecc1e9-f339-4d85-84c3-aa2fff47b2c1",
  pelican       = "2a168f51-4016-4a32-ba72-64a53bee31ca",
  hake          = "9e87f9e3-a1a6-49e0-9037-924d81ab53a3",
  rockfish      = "590dd989-103a-4b22-8788-81471386d860",
  forage_fish   = "5d05c64d-9232-45a4-84f5-6f6519c4a6fb"
)

cat("Downloading icon vector objects...\n")
local_icons <- lapply(icon_uuids, function(id) {
  # get_phylopic() downloads the parsed picture object directly
  get_phylopic(uuid = id)
})

# Save everything into a single lightweight binary file
rds_path <- file.path(cache_dir, "phylopic_icons.rds")
saveRDS(local_icons, file = rds_path)
cat(sprintf("Success! Saved %d icons locally to:\n%s\n", length(local_icons), rds_path))


# 1. Load your local icon bank (No internet needed)
icons <- readRDS("Rcode_for_paper/assets/phylopic_cache/phylopic_icons.rds")

# 2. Add icons directly from memory using the `img` parameter
p <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32") +
  theme_void() +
  
  # Chinook Salmon
  add_phylopic(
    img        = icons$chinook,
    x          = -125.5,
    y          = 46.8,
    height     = 0.65,
    fill       = "#00ACC1",
    angle      = 40
  ) +
  
  # Sleeper Shark
  add_phylopic(
    img        = icons$shark,
    x          = -149.0,
    y          = 53.5,
    height     = 1.40,
    fill       = "#08306B"
  ) +
  
  # Steller Sea Lion
  add_phylopic(
    img        = icons$ssl,
    x          = -128.8,
    y          = 49.8,
    height     = 1.50,
    fill       = "#372E29",
    horizontal = TRUE
  )