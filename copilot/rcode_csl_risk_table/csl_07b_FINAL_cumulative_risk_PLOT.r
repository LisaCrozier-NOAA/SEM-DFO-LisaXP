
# ==============================================================================
# Script: Ranked Risk Trajectory Plot (1998–2021)
# Compares CSL Count, Estuary Risk, Ocean Risk (Inverted Buffer), & Integrated Risk
# ==============================================================================

library(tidyverse)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Recreate/Load Data & Standardize Ordinal Risk Levels
# -----------------------------------------------------------------------------

# Mapping Risk Tiers to a numeric scale: 0 (Lowest Risk) -> 4 (Critical Risk)
risk_tier_map <- c(
  "LOW RISK"          = 0,
  "LOW-MODERATE RISK" = 1,
  "MODERATE RISK"     = 2,
  "HIGH RISK"         = 3,
  "VERY HIGH RISK"    = 4,
  "CRITICAL RISK"     = 4
)

plot_data <- csl_risk_table_classified %>%
  filter(year <= 2021) %>%
  mutate(
    # Invert Ocean Buffer Score to represent OCEAN RISK (5 Buffer = 0 Risk; 0 Buffer = 5 Risk)
    ocean_risk_numeric = 5 - ocean_buffer_score,
    
    # Map text risk categories to numeric ranks (0–4)
    estuary_risk_rank    = risk_tier_map[overall_chinook_risk],
    integrated_risk_rank = risk_tier_map[integrated_csl_risk],
    
    # Normalize CSL Count onto a comparable 0–4 scale for overlay
    csl_count_scaled = (csl_during_chinook - min(csl_during_chinook, na.rm = TRUE)) / 
      (max(csl_during_chinook, na.rm = TRUE) - min(csl_during_chinook, na.rm = TRUE)) * 4
  )

# -----------------------------------------------------------------------------
# 2. Build Panel 1: CSL Count Trajectory
# -----------------------------------------------------------------------------

p1 <- ggplot(plot_data, aes(x = year, y = csl_during_chinook)) +
  geom_line(color = "#1e3a8a", linewidth = 1.2) +
  geom_point(color = "#1e3a8a", size = 2.5) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "1. Raw Sea Lion Exposure (CSL Count)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11, color = "#1e3a8a"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 3. Build Panel 2: Estuary Risk Tier vs Ocean Risk Tier (Inverted Buffer)
# -----------------------------------------------------------------------------

p2_data <- plot_data %>%
  select(year, Estuary = estuary_risk_rank, `Ocean (Inverted Buffer)` = ocean_risk_numeric) %>%
  pivot_longer(-year, names_to = "Risk_Type", values_to = "Risk_Rank")

p2 <- ggplot(p2_data, aes(x = year, y = Risk_Rank, color = Risk_Type, linetype = Risk_Type)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Estuary" = "#d97706", "Ocean (Inverted Buffer)" = "#2563eb")) +
  scale_linetype_manual(values = c("Estuary" = "solid", "Ocean (Inverted Buffer)" = "dashed")) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(
    breaks = 0:4,
    labels = c("Low (0)", "Low-Mod (1)", "Mod (2)", "High (3)", "Very High (4)")
  ) +
  labs(
    title = "2. Estuary Risk vs Ocean Risk (Inverted Buffer: 0 = Low Risk, 5 = High Risk)",
    y = "Risk Level"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11, color = "#2563eb"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 4. Build Panel 3: Integrated Risk Tier Trajectory
# -----------------------------------------------------------------------------

p3 <- ggplot(plot_data, aes(x = year, y = integrated_risk_rank)) +
  geom_line(color = "#dc2626", linewidth = 1.3) +
  geom_point(color = "#dc2626", size = 2.8) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(
    breaks = 0:4,
    labels = c("LOW", "LOW-MOD", "MODERATE", "HIGH", "VERY HIGH/CRIT")
  ) +
  labs(
    title = "3. Integrated CSL Predation Risk Tier",
    x = "Year",
    y = "Integrated Risk Tier"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = "#dc2626"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 5. Combine & Render Multi-Panel Figure
# -----------------------------------------------------------------------------

p_final_comparison <- (p1 / p2 / p3) +
  plot_annotation(
    title = "Rank-Ordered CSL Predation Risk Trajectory (1998–2021)",
    subtitle = "Comparing Raw CSL Counts, Estuary Risk, Inverted Ocean Buffer Risk, and Integrated Risk Tiers"
  )

print(p_final_comparison)

# Save high-resolution graphic
ggsave("copilot/outputs_csl_cr/csl_risk_trajectory_rank_ordered.png", p_final_comparison, width = 12, height = 10, dpi = 300)
