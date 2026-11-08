# ==============================================================================
# Script: 01b_goal1_variable_importance_pointwhisker_plots.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Point & Whisker Rank-Ordered Frequency Distribution of Variable Importance
#          (Extends 01_goal1_variable_importance_and_plot.R)
# ==============================================================================

library(tidyverse)
library(ggrepel)

# ------------------------------------------------------------------------------
# STEP 1: LOAD & PREPARE ALL-MODEL IMPORTANCE SCORES
# ------------------------------------------------------------------------------
# Load full-ensemble score data generated in Step 4 of script 01
by_model_path <- file.path(data_out_dir, "Variable_Importance_By_ModelRun.csv")

if (!file.exists(by_model_path)) {
  stop("Missing input file: ", by_model_path, ". Please run Script 01 first.")
}

importance_by_model <- read_csv(by_model_path, show_col_types = FALSE)

# Filter out salmon variables, assign regional panels, and summarize min/mean/max
importance_summary <- importance_by_model %>%
  # Exclude salmon response nodes
  filter(!generic_node %in% c("Abundance", "Growth")) %>%
  mutate(
    # Classify into 4 subpanels
    Region = if_else(grepl("NCC", generic_node) | grepl("^X0[1-9]", LisaName), "NCC", "AK"),
    Trophic = case_when(
      grepl("Prey", generic_node) ~ "Prey",
      grepl("Pred", generic_node) ~ "Predator",
      TRUE ~ "Other"
    ),
    Subpanel = factor(
      paste0(Region, " ", Trophic),
      levels = c("NCC Prey", "NCC Predator", "AK Prey", "AK Predator")
    )
  ) %>%
  filter(Trophic %in% c("Prey", "Predator")) %>%
  
  # Group across all model types per indicator
  group_by(Subpanel, LisaName) %>%
  summarise(
    Mean_Imp = mean(Akaike_Importance_Weight, na.rm = TRUE),
    Min_Imp  = min(Akaike_Importance_Weight, na.rm = TRUE),
    Max_Imp  = max(Akaike_Importance_Weight, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  
  # Rank order variables within each subpanel (1 = Most Important)
  group_by(Subpanel) %>%
  arrange(desc(Mean_Imp)) %>%
  mutate(
    Rank = row_number(),
    # Top 3 flag for labeling variant
    Is_Top3 = if_else(Rank <= 3, LisaName, NA_character_)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# STEP 2: DEFINE CUSTOM COLOR PALETTE PER SUBPANEL
# ------------------------------------------------------------------------------
panel_colors <- c(
  "NCC Prey"     = "#E65100", # Warm Coral / Orange
  "NCC Predator" = "#2E7D32", # Green
  "AK Prey"      = "#0288D1", # Teal / Blue
  "AK Predator"  = "#7B1FA2"  # Purple
)

# Base ggplot theme matching reference panel layout
theme_publication_panels <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(1.2, "lines"),
      strip.background   = element_rect(fill = "grey92", color = "black", linewidth = 0.8),
      strip.text         = element_text(face = "bold", size = 11),
      axis.title         = element_text(face = "bold", size = 11),
      axis.text.x        = element_blank(), # Hidden X text for clean point distribution
      axis.ticks.x       = element_blank(),
      legend.position    = "none"
    )
}

# ------------------------------------------------------------------------------
# STEP 3: VARIANT A — NO VARIABLE NAMES (CLEAN FREQUENCY DISTRIBUTION)----------
# ------------------------------------------------------------------------------
plot_nolabels <- ggplot(
  importance_summary, 
  aes(x = Rank, y = Mean_Imp, color = Subpanel)
) +
  # Whiskers: Min to Max Range across model runs
  geom_errorbar(
    aes(ymin = Min_Imp, ymax = Max_Imp), 
    width = 0.2, 
    linewidth = 0.6,
    alpha = 0.7
  ) +
  # Points: Mean Importance Weight
  geom_point(size = 2.5) +
  # Connect rank distribution with thin line
  geom_line(linewidth = 0.4, alpha = 0.4) +
  facet_wrap(~ Subpanel, scales = "free_x", nrow = 1) +
  scale_color_manual(values = panel_colors) +
  scale_y_continuous(
    limits = c(0, 1.05), 
    breaks = seq(0, 1, by = 0.2),
    expand = c(0.02, 0)
  ) +
  labs(
    x = "Rank-Ordered Indicators (Most to Least Important)",
    y = "Mean Akaike Importance Weight (Range across models)",
    title = "Ecosystem Indicator Variable Importance Distribution"
  ) +
  theme_publication_panels()

# Save Variant A
ggsave(
  file.path(fig_out_dir, "Goal1_Importance_PointWhisker_NoLabels.png"),
  plot_nolabels, width = 10, height = 4.5, dpi = 300
)
ggsave(
  file.path(fig_out_dir, "Goal1_Importance_PointWhisker_NoLabels.pdf"),
  plot_nolabels, width = 10, height = 4.5
)

# ------------------------------------------------------------------------------
# STEP 4: VARIANT B — PRECISE MANUAL & DYNAMIC LABEL POSITIONING
# ------------------------------------------------------------------------------

# 1. Labels for NCC Prey (Dynamic repulsion)
labels_ncc_prey <- importance_summary %>%
  filter(Subpanel == "NCC Prey" & Rank <= 4)

# 2. Fixed Manual Coordinates for NCC Predator (Fixes point overlap & low axis collision)
labels_ncc_pred <- importance_summary %>%
  filter(Subpanel == "NCC Predator" & Rank <= 2) %>%
  mutate(
    x_pos = case_when(
      Rank == 1 ~ 1.35,  # Shifted slightly right of Rank 1 point
      Rank == 2 ~ 2.00   # Centered directly over Rank 2 point
    ),
    y_pos = case_when(
      Rank == 1 ~ 0.96,  # Positioned just below top point (y = 1.0)
      Rank == 2 ~ 0.05   # Elevated safely above y = 0.0 dots and axis line
    )
  )

# 3. Fixed Manual Coordinates for AK Panels (Staggered staircase above points)
labels_ak <- importance_summary %>%
  filter(str_detect(Subpanel, "^AK") & Rank <= 3) %>%
  mutate(
    y_pos = case_when(
      Rank == 1 ~ 0.45,
      Rank == 2 ~ 0.40,
      Rank == 3 ~ 0.35
    )
  )

plot_top3labels <- plot_nolabels +
  # --- Panel 1: NCC Prey (Repelled labels) ---
  geom_text_repel(
    data = labels_ncc_prey,
    aes(x = Rank, y = Mean_Imp, label = LisaName, color = Subpanel),
    size = 2.8,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  # --- Panel 2: NCC Predator (Manual positioning) ---
  geom_text(
    data = labels_ncc_pred,
    aes(x = x_pos, y = y_pos, label = LisaName, color = Subpanel),
    hjust = 0,         # Left-aligned for Rank 1, clean fit for Rank 2
    size = 2.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # --- Panels 3 & 4: AK Prey & AK Predator (High staggered staircase) ---
  geom_text(
    data = labels_ak,
    aes(x = Rank, y = y_pos, label = paste0(Rank, ". ", LisaName), color = Subpanel),
    hjust = 0,
    size = 2.7,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_segment(
    data = labels_ak,
    aes(x = Rank, xend = Rank, y = y_pos - 0.02, yend = Max_Imp + 0.02, color = Subpanel),
    linetype = "dotted",
    linewidth = 0.4,
    show.legend = FALSE
  )

# Save Updated Variant B
ggsave(
  file.path(fig_out_dir, "Goal1_Importance_PointWhisker_Top3Labeled_StaggeredAK.png"),
  plot_top3labels, width = 11, height = 5.2, dpi = 300
)
ggsave(
  file.path(fig_out_dir, "Goal1_Importance_PointWhisker_Top3Labeled_StaggeredAK.pdf"),
  plot_top3labels, width = 11, height = 5.2
)

print(plot_top3labels)