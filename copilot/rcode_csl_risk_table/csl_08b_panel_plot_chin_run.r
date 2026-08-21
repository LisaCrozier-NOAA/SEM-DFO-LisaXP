#top: run exposed during csl window in time series
#in pct (almost always 100%)


# ==============================================================================
# Script: Chinook Run Exposure % During CSL Occupancy Window (1998–2024)
# Calculates total vs window-restricted Chinook and plots % run exposure
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)

output_dir <- "copilot/outputs_csl_cr"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Calculate Total Annual Chinook Run vs CSL Window Run
# -----------------------------------------------------------------------------

# Load raw weekly estuary dataset to compute total annual run size
master_estuary <- read.csv(file.path(output_dir, "master_estuary_2wklagweekly_all_species_1998_2024.csv"))

annual_chinook_overlap <- master_estuary %>%
  group_by(year) %>%
  mutate(
    chin_est_max   = max(chin_estuary_count, na.rm = TRUE),
    chin_est_prop  = if_else(is.na(chin_est_max) | chin_est_max <= 0, 0, chin_estuary_count / chin_est_max),
    is_csl_active  = chin_est_prop >= 0.10
  ) %>%
  summarise(
    total_chinook_run   = sum(chin_estuary_count, na.rm = TRUE),
    csl_window_chinook  = sum(chin_estuary_count[is_csl_active], na.rm = TRUE),
    csl_during_chinook  = sum(csl_final[is_csl_active], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_run_exposed = (csl_window_chinook / total_chinook_run) * 100,
    pct_run_escaped = 100 - pct_run_exposed
  ) %>%
  filter(year <= 2024)

# -----------------------------------------------------------------------------
# 2. Panel A: Total Chinook Run vs. Run Exposed During CSL Window
# -----------------------------------------------------------------------------

p_chinook_volume <- ggplot(annual_chinook_overlap, aes(x = year)) +
  geom_col(aes(y = total_chinook_run / 1000, fill = "Total Spring Chinook Run"), alpha = 0.4, width = 0.75) +
  geom_col(aes(y = csl_window_chinook / 1000, fill = "Run Exposed During CSL Window"), alpha = 0.85, width = 0.75) +
  geom_line(aes(y = total_chinook_run / 1000), color = "#0284c7", linewidth = 1) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c(
      "Total Spring Chinook Run"        = "#93c5fd",
      "Run Exposed During CSL Window"  = "#1d4ed8"
    )
  ) +
  labs(
    title = "A. Annual Spring Chinook Estuary Volume: Total vs. CSL Co-Occurrence",
    subtitle = "Dark blue bars represent the portion of the run passing while CSLs actively occupy the estuary.",
    x = NULL,
    y = "Chinook Count (Thousands)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 3. Panel B: Percentage of Total Chinook Run Encountering CSL Window
# -----------------------------------------------------------------------------

p_pct_exposed <- ggplot(annual_chinook_overlap, aes(x = year, y = pct_run_exposed)) +
  geom_col(fill = "#1d4ed8", alpha = 0.85, width = 0.7) +
  geom_line(color = "#1e3a8a", linewidth = 1) +
  geom_point(color = "#1e3a8a", size = 2.2) +
  geom_hline(
    yintercept = mean(annual_chinook_overlap$pct_run_exposed, na.rm = TRUE), 
    linetype = "dashed", color = "#dc2626", linewidth = 0.8
  ) +
  annotate(
    "text", x = 1998.5, 
    y = mean(annual_chinook_overlap$pct_run_exposed, na.rm = TRUE) + 3,
    label = paste0("27-Year Mean Exposure: ", round(mean(annual_chinook_overlap$pct_run_exposed, na.rm = TRUE), 1), "%"),
    color = "#dc2626", fontface = "bold", size = 3.2, hjust = 0
  ) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(
    title = "B. Percentage of Total Spring Chinook Run Exposed to CSLs",
    subtitle = "Calculated as: (Chinook Count During CSL Window / Total Annual Chinook Run) × 100",
    x = "Year",
    y = "% Chinook Run Exposed"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#1e3a8a"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 4. Combine & Save Dedicated Figure
# -----------------------------------------------------------------------------

fig_chinook_exposure_pct <- (p_chinook_volume / p_pct_exposed) + plot_layout(heights = c(1.1, 1))

print(fig_chinook_exposure_pct)

ggsave(
  file.path(output_dir, "fig_chinook_run_exposure_pct_1998_2024.png"), 
  fig_chinook_exposure_pct, 
  width = 11, height = 9, dpi = 300
)