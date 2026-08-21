# ==============================================================================
# Script: Figure 2 - Dual Log-Axis Plot (Eulachon Left, CSL Right)
# Stable log-space transformation eliminating sec_axis NaN calculation errors
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)

output_dir <- "copilot/outputs_csl_cr"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load & Transform Data into Log10 Space
# -----------------------------------------------------------------------------

master_df <- read.csv(file.path(output_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv")) %>%
  filter(year <= 2024) %>%
  mutate(
    # Log10 transformation (+1 for 0 values)
    log_eul = log10(if_else(is.na(eulachon_during_chinook) | eulachon_during_chinook < 1, 1, eulachon_during_chinook)),
    log_csl = log10(if_else(is.na(csl_during_chinook) | csl_during_chinook < 1, 1, csl_during_chinook)),
    
    overall_chinook_risk = factor(
      overall_chinook_risk,
      levels = c("LOW RISK", "LOW-MODERATE RISK", "MODERATE RISK", "HIGH RISK", "VERY HIGH RISK")
    )
  )

# Set fixed visual ranges on log scale
# Eulachon Left Axis: 0 (1 count) to 7 (10 Million)
eul_min <- 0
eul_max <- 7

# CSL Right Axis: 2 (100 count) to 4.5 (~30,000 count)
csl_min <- 2
csl_max <- 4.5

# Linear mapping parameters between Log(CSL) and Log(Eulachon) visual spaces
slope <- (eul_max - eul_min) / (csl_max - csl_min)
intercept <- eul_min - slope * csl_min

# Map Log(CSL) into the Log(Eulachon) coordinate space
master_df <- master_df %>%
  mutate(csl_mapped = log_csl * slope + intercept)

# Map thresholds to visual coordinates
eul_thresh_val <- log10(50000)                   # ~4.70
csl_thresh_val <- log10(3000) * slope + intercept # ~4.36

# -----------------------------------------------------------------------------
# 2. Panel A: Dual Log-Axis Time Series
# -----------------------------------------------------------------------------

p2a_dual_log <- ggplot(master_df, aes(x = year)) +
  # Eulachon Line (Left Axis, Green)
  geom_line(aes(y = log_eul, color = "Eulachon Index"), linewidth = 1.1) +
  geom_point(aes(y = log_eul, color = "Eulachon Index"), size = 2.2) +
  
  # CSL Line (Right Axis mapped through linear transformation, Red)
  geom_line(aes(y = csl_mapped, color = "CSL Count"), linewidth = 1.1) +
  geom_point(aes(y = csl_mapped, color = "CSL Count"), size = 2.2) +
  
  # Dotted Threshold Lines
  geom_hline(yintercept = eul_thresh_val, linetype = "dashed", color = "#059669", linewidth = 0.8) +
  annotate(
    "text", x = 1998.5, y = eul_thresh_val + 0.35, 
    label = "50K Eulachon Buffer Threshold", 
    color = "#059669", fontface = "bold", size = 3, hjust = 0
  ) +
  
  geom_hline(yintercept = csl_thresh_val, linetype = "dashed", color = "#dc2626", linewidth = 0.8) +
  annotate(
    "text", x = 1998.5, y = csl_thresh_val - 0.35, 
    label = "3K CSL High Exposure Threshold", 
    color = "#dc2626", fontface = "bold", size = 3, hjust = 0
  ) +
  
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  
  # Left Axis: Log(Eulachon) mapped back to true counts
  # Right Axis: Log(CSL) mapped back to true counts using the inverse linear formula
  scale_y_continuous(
    name = "Eulachon Index (Log Scale)",
    limits = c(eul_min, eul_max),
    breaks = log10(c(1, 100, 1000, 10000, 50000, 100000, 1000000, 8000000)),
    labels = comma(c(1, 100, 1000, 10000, 50000, 100000, 1000000, 8000000)),
    sec.axis = sec_axis(
      transform = ~ (. - intercept) / slope,
      name = "CSL Count (Log Scale)",
      breaks = log10(c(100, 300, 1000, 3000, 10000)),
      labels = comma(c(100, 300, 1000, 3000, 10000))
    )
  ) +
  scale_color_manual(
    values = c(
      "Eulachon Index" = "#059669", 
      "CSL Count"      = "#dc2626"
    )
  ) +
  labs(
    title = "A. Eulachon Prey Availability vs. CSL Predator Exposure (Dual Log Axes)",
    subtitle = "Left Axis = Eulachon Index; Right Axis = CSL Count. Dotted lines mark critical risk thresholds.",
    x = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    axis.title.y.left = element_text(color = "#059669", face = "bold"),
    axis.title.y.right = element_text(color = "#dc2626", face = "bold"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 3. Panel B: Mechanism Boxplot (Eulachon Deficit Across Risk Tiers)
# -----------------------------------------------------------------------------

p2b_dual_log <- ggplot(
  master_df %>% filter(!is.na(overall_chinook_risk)), 
  aes(x = overall_chinook_risk, y = 10^log_eul)
) +
  geom_boxplot(aes(fill = overall_chinook_risk), alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.6) +
  geom_hline(yintercept = 50000, linetype = "dashed", color = "#059669", linewidth = 0.8) +
  annotate(
    "text", x = 1.3, y = 110000, 
    label = "50K Buffer Threshold", 
    color = "#059669", fontface = "bold", size = 3
  ) +
  scale_y_log10(
    labels = comma,
    breaks = c(1, 100, 1000, 10000, 50000, 100000, 1000000, 8000000)
  ) +
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
    title = "B. Mechanism: Risk Escalation Corresponds to Eulachon Deficits (<50K)",
    subtitle = "HIGH and VERY HIGH Risk years systematically fall below the 50K Eulachon buffer line.",
    x = "Estuary Risk Tier",
    y = "Eulachon Index (Log Scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 4. Combine & Save Figure
# -----------------------------------------------------------------------------

fig2_eulachon_risk_dual_log <- (p2a_dual_log / p2b_dual_log) + plot_layout(heights = c(1.2, 1))

print(fig2_eulachon_risk_dual_log)

ggsave(
  file.path(output_dir, "fig2_eulachon_heavy_risk_dynamics_dual_log.png"), 
  fig2_eulachon_risk_dual_log, 
  width = 11, height = 9, dpi = 300
)
