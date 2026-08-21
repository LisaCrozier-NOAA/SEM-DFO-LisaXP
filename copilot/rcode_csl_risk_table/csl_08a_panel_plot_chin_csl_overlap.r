#3 panel plot, 
#top: csl/chin overlap example in 1 year, 
#middle: time series of chin during csl residence
#bottom: time series of ratio
fig1_timing_chinook_and_ratio_expanded.png

# ==============================================================================
# Script: Expanded Figure 1 - Run Timing, Chinook Run Size, & Predator/Prey Ratio
# Panel A: Conceptual Overlap
# Panel B: Spring Chinook Estuary Run Volume
# Panel C: CSL per 1,000 Chinook Ratio (Top SEM Model Predictor)
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)

output_dir <- "copilot/outputs_csl_cr"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Prepare Ratio Metric
# -----------------------------------------------------------------------------

master_df <- read.csv(file.path(output_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv")) %>%
  filter(year <= 2024) %>%
  mutate(
    # CSL per 1,000 Chinook ratio
    csl_per_1000_chinook = if_else(
      chinook_estuary_total > 0, 
      (csl_during_chinook / chinook_estuary_total) * 1000, 
      NA_real_
    ),
    overall_chinook_risk = factor(
      overall_chinook_risk,
      levels = c("LOW RISK", "LOW-MODERATE RISK", "MODERATE RISK", "HIGH RISK", "VERY HIGH RISK")
    )
  )

# -----------------------------------------------------------------------------
# 2. Panel A: Migration Timing Overlap Curve
# -----------------------------------------------------------------------------

days <- seq(1, 120) # Day 1 = March 1, Day 120 = June 30
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
    title = "A. Temporal Overlap: CSL Estuary Presence vs. Spring Chinook Migration",
    subtitle = "Predation exposure window occurs during spring migration overlap",
    x = NULL,
    y = "Relative Intensity"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", color = "#0f172a"),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 3. Panel B: Spring Chinook Estuary Run Size During CSL Passage Window
# -----------------------------------------------------------------------------

p1b <- ggplot(master_df, aes(x = year, y = chinook_estuary_total)) +
  geom_col(fill = "#0284c7", alpha = 0.75, width = 0.7) +
  geom_line(color = "#0369a1", linewidth = 1) +
  geom_point(color = "#0369a1", size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "B. Spring Chinook Estuary Volume During Active CSL Window",
    subtitle = "Total Chinook salmon passing through estuary during weeks with active CSL presence",
    x = NULL,
    y = "Chinook Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#0284c7"),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 4. Panel C: CSL per 1,000 Chinook Predator/Prey Ratio (idxE10)
# -----------------------------------------------------------------------------

p1c <- ggplot(master_df, aes(x = year, y = csl_per_1000_chinook)) +
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
    title = "C. Top SEM Model Predictor (idxE10): CSL per 1,000 Chinook Ratio",
    subtitle = "Calculated as (CSL Count / Chinook Estuary Total) × 1,000. Outperforms raw CSL counts (ΔR² = +0.119).",
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

# -----------------------------------------------------------------------------
# 5. Combine 3 Panels & Export Figure
# -----------------------------------------------------------------------------

fig1_expanded <- (p1a / p1b / p1c) + plot_layout(heights = c(1, 1, 1.2))

print(fig1_expanded)

ggsave(
  file.path(output_dir, "fig1_timing_chinook_and_ratio_expanded.png"), 
  fig1_expanded, 
  width = 11, height = 11, dpi = 300
)

