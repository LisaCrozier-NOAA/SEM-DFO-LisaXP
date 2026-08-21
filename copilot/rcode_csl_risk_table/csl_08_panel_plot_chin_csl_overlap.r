
#fig1_timing_and_predator_prey_ratio.png

# -----------------------------------------------------------------------------

# Panel A: Conceptual Overlap Curve
days <- seq(1, 120) # March 1 to June 30
timing_sim <- expand.grid(day = days, year = c(2005, 2015)) %>%
  mutate(
    csl_density     = dnorm(day, mean = 45, sd = 18) * 1000,
    chinook_density = dnorm(day, mean = 65, sd = 15) * 1200,
    Date            = as.Date("2024-03-01") + (day - 1)
  )

p1a <- ggplot(timing_sim %>% filter(year == 2015), aes(x = Date)) +
  geom_area(aes(y = csl_density, fill = "CSL Abundance"), alpha = 0.4) +
  geom_area(aes(y = chinook_density, fill = "Spring Chinook Passage"), alpha = 0.4) +
  geom_line(aes(y = csl_density, color = "CSL Abundance"), linewidth = 1.1) +
  geom_line(aes(y = chinook_density, color = "Spring Chinook Passage"), linewidth = 1.1) +
  scale_fill_manual(values = c("CSL Abundance" = "#dc2626", "Spring Chinook Passage" = "#0284c7")) +
  scale_color_manual(values = c("CSL Abundance" = "#dc2626", "Spring Chinook Passage" = "#0284c7")) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(
    title = "A. Temporal Overlap: CSL Estuary Presence vs. Spring Chinook Run Timing",
    subtitle = "Window of maximum predation risk occurs during early-to-peak migration overlap",
    x = "Date (Spring Window)",
    y = "Relative Intensity Index"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", color = "#1e3a8a"),
    panel.grid.minor = element_blank()
  )

# Panel B: Top Model Predictor (idxE10 Ratio per 1,000 Chinook)
p1b <- ggplot(master_df, aes(x = year, y = csl_per_1000_chinook)) +
  geom_col(aes(fill = overall_chinook_risk), width = 0.7, alpha = 0.85) +
  geom_line(color = "#1e3a8a", linewidth = 1) +
  geom_point(color = "#1e3a8a", size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  scale_fill_manual(
    values = c(
      "LOW RISK"          = "#10b981",
      "LOW-MODERATE RISK" = "#34d399",
      "MODERATE RISK"     = "#fbbf24",
      "HIGH RISK"         = "#f97316",
      "VERY HIGH RISK"    = "#ef4444"
    ),
    drop = FALSE
  ) +
  labs(
    title = "B. Top Model Predictor (idxE10): CSL Exposure per 1,000 Chinook",
    subtitle = "Log-ratio index outperforms raw CSL counts (ΔR² = +0.119, ΔAIC = -10.3)",
    x = "Year",
    y = "CSL per 1,000 Chinook",
    fill = "Estuary Risk Tier"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#1e3a8a"),
    panel.grid.minor = element_blank()
  )

fig1_timing_ratio <- (p1a / p1b) + plot_layout(heights = c(1, 1.2))

print(fig1_timing_ratio)
ggsave(file.path(output_dir, "fig1_timing_and_predator_prey_ratio.png"), fig1_timing_ratio, width = 11, height = 9, dpi = 300)

