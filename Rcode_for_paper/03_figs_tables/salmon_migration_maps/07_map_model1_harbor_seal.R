# ==============================================================================
# Script: map_summary_model1_harbor_seal.R
# Purpose: Minimalist map showing ONLY the Hake + Harbor Seal + CR SSL paradigm
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(rphylopic)

fig_out_dir <- file.path(proj_dir, "03_figs_tables", "salmon_migration_maps", "figs")

land <- ne_countries(scale = "medium", returnclass = "sf")

chinook_uuid <- get_uuid(name = "Oncorhynchus tshawytscha")
hake_uuid    <- "9e87f9e3-a1a6-49e0-9037-924d81ab53a3"
#seal_uuid    <- "1925a9aa-bf33-4270-824c-b8f781f8ec82"
seal_uuid    <-  get_uuid(name = "Phoca vitulina")
ssl_uuid     <- "21db0060-c29a-4870-96a4-629d7424647f"

p_map1 <- ggplot() +
  geom_sf(data = land, fill = "#c5cc9d", color = "#404a32", linewidth = 0.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#6fa2bb", color = NA)) +
  
  # 1. Early smolt outmigration
  add_phylopic(uuid = chinook_uuid, x = -125.5, y = 46.8, height = 0.65, fill = "#B2EBF2", color = "#B2EBF2", angle = 40) +
  add_phylopic(uuid = chinook_uuid, x = -128.2, y = 49.5, height = 0.85, fill = "#00ACC1", color = "#00ACC1", angle = 30) +
  
  # 2. Hake arena (NCC Shelf)
  add_phylopic(uuid = hake_uuid, x = -126.8, y = 46.1, height = 0.85, fill = "#38BDF8", color = "#0284C7", angle = 20) +
#   "Pacific Hake", hake_uuid,   -126.8,   45.8,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
# "Pacific Hake", hake_uuid,   -126.8,   46.1,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
# "Pacific Hake", hake_uuid,   -126.2,   45.8,    0.85,     -50,    TRUE, "#38BDF8", "#0284C7",
annotate("text", x = -128.8, y = 46.5, label = "Pacific Hake", size = 3.5, fontface = "bold", color = "white") +
  
  # 3. Harbour Seals (BC Coast / Strait)
  add_phylopic(uuid = seal_uuid, x = -127.5, y = 50.8, height = 1.10, fill = "#1F2937", color = "#111827", angle = 15) +
  add_phylopic(uuid = seal_uuid, x = -125.4, y = 49.8, height = 0.85, fill = "#1F2937", color = "#111827", angle = 15, horizontal = FALSE) +
  add_phylopic(uuid = seal_uuid, x = -132.0, y = 53.0, height = 1.10, fill = "#1F2937", color = "#111827", angle = 10, horizontal = TRUE) +
  add_phylopic(uuid = seal_uuid, x = -124.0, y = 46.3, height = 0.85, fill = "#1F2937", color = "#111827", angle = 15, horizontal = FALSE) +
#  "Harbor Seal",  seal_uuid,   -124.0,   46.3,    0.60,     10,    FALSE, "#475569", "#1E293B",
#  "Harbour Seal", seal_uuid,  -125.4,  49.8,    0.85,     15,    FALSE, "#4B5563", "#1F2937",
#  add_phylopic(uuid = seal_uuid, x = -131.0, y = 53.0, height = 1.10, fill = "#1F2937", color = "#111827", angle = -10, horizontal = TRUE) +
  annotate("text", x = -133.0, y = 51.5, label = "Coastal Predation\n(Harbour Seals)", size = 3.5, fontface = "bold", color = "white") +
  
  # 4. In-River Sea Lions (Bonneville bottleneck)
  add_phylopic(uuid = ssl_uuid, x = -122.0, y = 45.9, height = 1.00, fill = "#1F1A17", color = "#1F1A17") +
  annotate("text", x = -120.0, y = 45.2, label = "Adult Bottleneck\n(Bonneville SSL)", size = 3.2, fontface = "italic", color = "white") +
  
  coord_sf(xlim = c(-137, -119), ylim = c(43, 56), expand = FALSE) +
  labs(
    title = "Paradigm 1: Direct Coastal Predation (Unmodulated)",
    subtitle = "Hake (Smolt CPUE) + DFO Coastwide Harbour Seal + Bonneville Dam Sea Lions | R² ≈ 0.51"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black", margin = margin(t = 10, b = 4)),
    plot.subtitle = element_text(size = 11, face = "italic", color = "gray25", margin = margin(b = 10)),
    plot.margin   = margin(10, 15, 10, 15)
  )

ggsave(  filename = file.path(fig_out_dir, "07_map_harbor_seal_model_CR_WCVI.png"),
      p_map1, width = 8, height = 6.5, dpi = 300, bg = "white")

#old: "Rcode_for_paper/Routput_for_paper/plots/Summary_Map_Model1_HarborSeal.png", 
