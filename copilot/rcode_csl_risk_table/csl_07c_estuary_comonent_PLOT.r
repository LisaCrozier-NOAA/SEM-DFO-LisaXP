# ==============================================================================
# Script: Estuary Components Breakdown (1998–2021)
# ==============================================================================

# Mapping Estuary Risk Tiers to a numeric scale
risk_tier_map <- c(
  "LOW RISK"          = 0,
  "LOW-MODERATE RISK" = 1,
  "MODERATE RISK"     = 2,
  "HIGH RISK"         = 3,
  "VERY HIGH RISK"    = 4
)

estuary_plot_data <- csl_risk_table_classified %>%
  filter(year <= 2021) %>%
  mutate(estuary_risk_rank = risk_tier_map[overall_chinook_risk])

# Panel A: CSL Predator Exposure
p_est_csl <- ggplot(estuary_plot_data, aes(x = year, y = csl_during_chinook)) +
  geom_area(fill = "#ef4444", alpha = 0.2) +
  geom_line(color = "#dc2626", linewidth = 1.1) +
  geom_point(color = "#dc2626", size = 2) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "A. CSL Predator Exposure During Chinook Migration Window",
    y = "CSL Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank())

# Panel B: In-River Alternate Prey Buffers (Eulachon & Shad)
p_est_buffers <- ggplot(estuary_plot_data, aes(x = year)) +
  geom_line(aes(y = eulachon_during_chinook, color = "Eulachon (Early Buffer)"), linewidth = 1) +
  geom_point(aes(y = eulachon_during_chinook, color = "Eulachon (Early Buffer)"), size = 2) +
  geom_line(aes(y = shad_during_chinook, color = "Shad (Late Buffer)"), linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = shad_during_chinook, color = "Shad (Late Buffer)"), size = 2) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = c("Eulachon (Early Buffer)" = "#059669", "Shad (Late Buffer)" = "#2563eb")) +
  labs(
    title = "B. Estuary Alternate Prey Biomass Buffers (Log Scale)",
    y = "Passage Count (Log Scale)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel C: Resulting Estuary Risk Tier
p_est_risk <- ggplot(estuary_plot_data, aes(x = year, y = estuary_risk_rank)) +
  geom_line(color = "#d97706", linewidth = 1.2) +
  geom_point(color = "#d97706", size = 2.5) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  scale_y_continuous(
    breaks = 0:4,
    labels = c("LOW", "LOW-MOD", "MODERATE", "HIGH", "VERY HIGH")
  ) +
  labs(
    title = "C. Baseline Estuary Predation Risk Tier",
    x = "Year",
    y = "Estuary Risk Tier"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

# Combine into 3-Panel Stack
p_estuary_components <- (p_est_csl / p_est_buffers / p_est_risk) +
  plot_annotation(
    title = "Estuary Risk Component Breakdown (1998–2021)",
    subtitle = "Comparing CSL Predator Numbers, Eulachon/Shad Estuary Buffers, and the Final Estuary Risk Tier"
  )

print(p_estuary_components)

ggsave("copilot/outputs_csl_cr/estuary_risk_components_breakdown.png", p_estuary_components, width = 11, height = 9, dpi = 300)
