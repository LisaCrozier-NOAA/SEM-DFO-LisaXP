library(ggplot2)
library(grid)

# Create canvas
p <- ggplot() + 
  theme_void() + 
  xlim(0, 10) + ylim(0, 7.5) +
  
  # --- Slide Header ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 6.6, ymax = 7.3, fill = "#1E293B", color = NA) +
  annotate("text", x = 5, y = 6.95, label = "5-Guild Predator Classification & Functional Prey Pools", 
           color = "white", fontface = "bold", size = 5) +
  
  # --- Table Column Headers ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 6.0, ymax = 6.5, fill = "#334155", color = NA) +
  annotate("text", x = 0.9, y = 6.25, label = "Guild", color = "white", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = 1.8, y = 6.25, label = "Region", color = "white", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 6.25, label = "Predator Group / Key Species", color = "white", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = 6.3, y = 6.25, label = "Candidate Alternate Prey Pool", color = "white", fontface = "bold", size = 3.5, hjust = 0) +
  
  # --- Row 1: Guild 1 ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 4.9, ymax = 5.9, fill = "#F8FAFC", color = "#E2E8F0") +
  annotate("rect", xmin = 0.5, xmax = 1.6, ymin = 4.9, ymax = 5.9, fill = "#EFF6FF", color = NA) +
  annotate("text", x = 0.7, y = 5.4, label = "Guild 1", fontface = "bold", color = "#1E40AF", size = 3.8, hjust = 0) +
  annotate("text", x = 1.8, y = 5.4, label = "NCC", fontface = "bold", color = "#475569", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 5.6, label = "Avian Predators & Small Porpoise", fontface = "bold", color = "#0F172A", size = 3.4, hjust = 0) +
  annotate("text", x = 2.8, y = 5.2, label = "Cormorants, Harbor Porpoise", color = "#64748B", size = 3.0, hjust = 0) +
  annotate("text", x = 6.3, y = 5.4, label = "Market Squid, Anchovy, Sardine, Herring, Eulachon", color = "#1E293B", size = 3.2, hjust = 0) +
  
  # --- Row 2: Guild 2 ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 3.8, ymax = 4.8, fill = "#FFFFFF", color = "#E2E8F0") +
  annotate("rect", xmin = 0.5, xmax = 1.6, ymin = 3.8, ymax = 4.8, fill = "#F0FDF4", color = NA) +
  annotate("text", x = 0.7, y = 4.3, label = "Guild 2", fontface = "bold", color = "#166534", size = 3.8, hjust = 0) +
  annotate("text", x = 1.8, y = 4.3, label = "NCC", fontface = "bold", color = "#475569", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 4.5, label = "Mammals & Demersal Fish", fontface = "bold", color = "#0F172A", size = 3.4, hjust = 0) +
  annotate("text", x = 2.8, y = 4.1, label = "Pinnipeds (Harbor/Fur Seals, Sea Lions), Hake", color = "#64748B", size = 3.0, hjust = 0) +
  annotate("text", x = 6.3, y = 4.3, label = "Pacific Hake + Guild 1 Forage Pool (Squid, Anchovy, etc.)", color = "#1E293B", size = 3.2, hjust = 0) +
  
  # --- Row 3: Guild 3 ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 2.7, ymax = 3.7, fill = "#F8FAFC", color = "#E2E8F0") +
  annotate("rect", xmin = 0.5, xmax = 1.6, ymin = 2.7, ymax = 3.7, fill = "#FFF7ED", color = NA) +
  annotate("text", x = 0.7, y = 3.2, label = "Guild 3", fontface = "bold", color = "#C2410C", size = 3.8, hjust = 0) +
  annotate("text", x = 1.8, y = 3.2, label = "AK", fontface = "bold", color = "#475569", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 3.4, label = "Mammals & Fish (krill eaters)", fontface = "bold", color = "#0F172A", size = 3.4, hjust = 0) +
  annotate("text", x = 2.8, y = 3.0, label = "Steller Sea Lions, Sablefish, Spiny Dogfish", color = "#64748B", size = 3.0, hjust = 0) +
  annotate("text", x = 6.3, y = 3.2, label = "Walleye Pollock (Age 1+), Sitka Herring, Capelin, Krill", color = "#1E293B", size = 3.2, hjust = 0) +
  
  # --- Row 4: Guild 4 ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 1.6, ymax = 2.6, fill = "#FFFFFF", color = "#E2E8F0") +
  annotate("rect", xmin = 0.5, xmax = 1.6, ymin = 1.6, ymax = 2.6, fill = "#FAF5FF", color = NA) +
  annotate("text", x = 0.7, y = 2.1, label = "Guild 4", fontface = "bold", color = "#7E22CE", size = 3.8, hjust = 0) +
  annotate("text", x = 1.8, y = 2.1, label = "AK", fontface = "bold", color = "#475569", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 2.3, label = "Demersal Piscivores (Gadids & Flatfish)", fontface = "bold", color = "#0F172A", size = 3.4, hjust = 0) +
  annotate("text", x = 2.8, y = 1.9, label = "Pacific Halibut, Pacific Cod, Arrowtooth Flounder", color = "#64748B", size = 3.0, hjust = 0) +
  annotate("text", x = 6.3, y = 2.1, label = "Walleye Pollock (Age 1+), Sitka Herring, Capelin", color = "#1E293B", size = 3.2, hjust = 0) +
  
  # --- Row 5: Guild 5 ---
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 0.5, ymax = 1.5, fill = "#F8FAFC", color = "#E2E8F0") +
  annotate("rect", xmin = 0.5, xmax = 1.6, ymin = 0.5, ymax = 1.5, fill = "#FEF2F2", color = NA) +
  annotate("text", x = 0.7, y = 1.0, label = "Guild 5", fontface = "bold", color = "#B91C1C", size = 3.8, hjust = 0) +
  annotate("text", x = 1.8, y = 1.0, label = "AK", fontface = "bold", color = "#475569", size = 3.5, hjust = 0) +
  annotate("text", x = 2.8, y = 1.2, label = "Apex Pelagic Sharks", fontface = "bold", color = "#0F172A", size = 3.4, hjust = 0) +
  annotate("text", x = 2.8, y = 0.8, label = "Salmon Shark, Sleeper Shark", color = "#64748B", size = 3.0, hjust = 0) +
  annotate("text", x = 6.3, y = 1.0, label = "Pink Salmon, Walleye Pollock, Sitka Herring, Capelin", color = "#1E293B", size = 3.2, hjust = 0) +
  
  # --- Footer Note ---
  annotate("text", x = 5.0, y = 0.2, label = "Note: Spatial overlap indicators (Habitat Suitability / Euphotic Shelf Biomass in NCC; SST in AK) enter as additional gating candidates.", 
           color = "#475569", size = 2.8, fontface = "italic")

# Render in RStudio
print(p)

# Save high-res PNG for PowerPoint
ggsave(file.path("Rcode_for_paper","figures","predator_guild_assignments_slide.png"), plot = p, width = 11, height = 7, dpi = 300)
