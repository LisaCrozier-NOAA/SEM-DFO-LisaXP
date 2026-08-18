# ==============================================================================
# Script: Ocean Components Breakdown (Guaranteed Custom Stack Heights)
# Stack Order (Bottom -> Top): Hake -> Anchovy -> Sardine -> NCC Herring -> Sitka Herring
# ==============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Prepare Data & Compute Explicit Rectangle Coordinates per Species
# -----------------------------------------------------------------------------

# Desired visual order from bottom to top
species_order <- c("Hake (Age 5+)", "Anchovy", "Sardine", "NCC Herring", "Sitka Herring")

ocean_components_fixed <- csl_risk_table_classified %>%
  filter(year <= 2021) %>%
  select(
    year, 
    `Hake (Age 5+)` = hake_age5plus_tier,
    `Anchovy`       = ncc_anchovy_tier,
    `Sardine`       = sardine_tier,
    `NCC Herring`   = ncc_herring_tier,
    `Sitka Herring` = sitka_herring_tier,
    ocean_buffer_score
  ) %>%
  pivot_longer(
    cols = -c(year, ocean_buffer_score),
    names_to = "Prey_Component",
    values_to = "Tier_Status"
  ) %>%
  filter(Tier_Status != "Low Biomass (Upper Risk)") %>%
  mutate(
    Prey_Component = factor(Prey_Component, levels = species_order)
  ) %>%
  arrange(year, Prey_Component) %>%
  group_by(year) %>%
  mutate(
    # Compute physical bar positions so level 1 is ALWAYS on the ground (0 to 1)
    stack_pos = row_number(),
    ymin = stack_pos - 1,
    ymax = stack_pos,
    xmin = year - 0.35,
    xmax = year + 0.35
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2. Render Plot with geom_rect
# -----------------------------------------------------------------------------

p_ocean_components_guaranteed <- ggplot() +
  # Draw explicit stacked rectangles (Hake always stays at the base)
  geom_rect(
    data = ocean_components_fixed,
    aes(
      xmin = xmin, xmax = xmax, 
      ymin = ymin, ymax = ymax, 
      fill = Prey_Component
    ),
    color = "white",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  # Line overlay showing total Ocean Buffer Score
  geom_line(
    data = csl_risk_table_classified %>% filter(year <= 2021),
    aes(x = year, y = ocean_buffer_score),
    color = "#1e3a8a",
    linewidth = 1.2
  ) +
  geom_point(
    data = csl_risk_table_classified %>% filter(year <= 2021),
    aes(x = year, y = ocean_buffer_score),
    color = "#1e3a8a",
    size = 2.5
  ) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(breaks = 0:5, limits = c(0, 5.5)) +
  scale_fill_manual(
    values = c(
      "Hake (Age 5+)" = "#1e40af", # Dark Blue (Bottom)
      "Anchovy"       = "#059669", # Green
      "Sardine"       = "#d97706", # Amber/Orange
      "NCC Herring"   = "#7c3aed", # Purple
      "Sitka Herring" = "#ec4899"  # Pink (Top)
    )
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title = "Ocean Forage Buffer Breakdown (1998–2021)",
    subtitle = "Guaranteed Stack Order (Bottom to Top): Hake → Anchovy → Sardine → NCC Herring → Sitka Herring\nDark Blue Line = Composite Ocean Buffer Score (0–5)",
    x = "Year",
    y = "Ocean Buffer Score (Count of Species Available)",
    fill = "Prey Component"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12, color = "#1e3a8a"),
    panel.grid.minor = element_blank()
  )

print(p_ocean_components_guaranteed)

# Save output
ggsave(
  "copilot/outputs_csl_cr/ocean_buffer_components_guaranteed_hake_bottom.png", 
  p_ocean_components_guaranteed, 
  width = 12, height = 6.5, dpi = 300
)