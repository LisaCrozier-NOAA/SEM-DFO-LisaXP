# ==============================================================================
# Script: Merge All Species (Eulachon, CSL, Chinook, Shad) & Plot Overlap
# Multi-Species Standardized Alignment (1998–2024, Weeks 10–26)
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# -----------------------------------------------------------------------------
# 1. Read Clean Individual Datasets
# -----------------------------------------------------------------------------

# Reconstructed Eulachon (1998–2024)
eulachon_df <- read.csv(file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"))

# Reconstructed CSL (1998–2024)
csl_df <- read.csv(file.path(output_dir, "csl_reconstructed_weekly_1998_2024.csv"))

# Interpolated Chinook Salmon (1998–2024)
chinook_df <- read.csv(file.path(output_dir, "bonn_chinook_weekly_interpolated_1998_2024.csv"))

# Interpolated American Shad (1998–2024, from standalone script)
shad_df<-read.csv(file.path(output_dir, "bonn_shad_weekly_interpolated_1998_2024.csv"), row.names = NULL)

# Helper function to prevent max() from throwing -Inf warnings on all-NA years
safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}
# -----------------------------------------------------------------------------
# 2. Build Master Combined Grid & Calculate Proportions (0 to 1 Scale)
# -----------------------------------------------------------------------------

master_all_species <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(eulachon_df %>% select(year, week, eulachon_final), by = c("year", "week")) %>%
  left_join(csl_df %>% select(year, week, csl_final), by = c("year", "week")) %>%
  left_join(chinook_df %>% select(year, week, Chin_weekly_count), by = c("year", "week")) %>%
  left_join(shad_df %>% select(year, week, Shad_weekly_count), by = c("year", "week")) %>%
  # Fill missing zeros for passage counts
  mutate(
    Chin_weekly_count = replace_na(Chin_weekly_count, 0),
    Shad_weekly_count = replace_na(Shad_weekly_count, 0)
  ) %>%
  group_by(year) %>%
  mutate(
    # Eulachon Proportions
    eul_max  = safe_max(eulachon_final),
    eul_prop = if_else(is.na(eul_max) | eul_max <= 0, 0, eulachon_final / eul_max),
    
    # CSL Proportions
    csl_max  = safe_max(csl_final),
    csl_prop = if_else(is.na(csl_max) | csl_max <= 0, 0, csl_final / csl_max),
    
    # Chinook Proportions
    chin_max  = safe_max(Chin_weekly_count),
    chin_prop = if_else(is.na(chin_max) | chin_max <= 0, 0, Chin_weekly_count / chin_max),
    
    # Shad Proportions
    shad_max  = safe_max(Shad_weekly_count),
    shad_prop = if_else(is.na(shad_max) | shad_max <= 0, 0, Shad_weekly_count / shad_max)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. Multi-Species Overlap Plots (Faceted 1998–2024)
# -----------------------------------------------------------------------------

# Plot 1: Full Faceted Multi-Species Overlap Series (1998–2024)
p_multispecies_full <- ggplot(master_all_species, aes(x = week)) +
  geom_line(aes(y = eul_prop, color = "Eulachon (Reconstructed)"), linewidth = 0.9) +
  geom_line(aes(y = csl_prop, color = "CSL (Reconstructed)"), linewidth = 0.9) +
  geom_line(aes(y = chin_prop, color = "Chinook Salmon (BON)"), linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(y = shad_prop, color = "American Shad (BON)"), linewidth = 0.8, linetype = "dotted") +
  facet_wrap(~ year, ncol = 4) +
  scale_color_manual(
    values = c(
      "Eulachon (Reconstructed)" = "seagreen4",
      "CSL (Reconstructed)"      = "steelblue",
      "Chinook Salmon (BON)"     = "firebrick",
      "American Shad (BON)"      = "goldenrod1"
    )
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
  labs(
    title = "Standardized Weekly Multi-Species Alignment (1998–2024, Weeks 10–26)",
    subtitle = "Comparing predator (CSL) and prey (Eulachon, Chinook, Shad) seasonal phenology",
    x = "Week (10–26)",
    y = "Proportion of Annual Peak",
    color = "Species Series"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_multispecies_full)

# Plot 2: Recent Training Period Focused View (2011–2024)
p_multispecies_recent <- ggplot(master_all_species %>% filter(year >= 2011), aes(x = week)) +
  geom_line(aes(y = eul_prop, color = "Eulachon (Reconstructed)"), linewidth = 1) +
  geom_line(aes(y = csl_prop, color = "CSL (Reconstructed)"), linewidth = 1) +
  geom_line(aes(y = chin_prop, color = "Chinook Salmon (BON)"), linewidth = 0.9, linetype = "dashed") +
  geom_line(aes(y = shad_prop, color = "American Shad (BON)"), linewidth = 0.9, linetype = "dotted") +
  facet_wrap(~ year, ncol = 4) +
  scale_color_manual(
    values = c(
      "Eulachon (Reconstructed)" = "seagreen4",
      "CSL (Reconstructed)"      = "steelblue",
      "Chinook Salmon (BON)"     = "firebrick",
      "American Shad (BON)"      = "goldenrod1"
    )
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
  labs(
    title = "Standardized Weekly Multi-Species Alignment (2011–2024, Weeks 10–26)",
    subtitle = "Proportion of annual peak biomass/count within each year",
    x = "Week (10–26)",
    y = "Proportion of Annual Peak",
    color = "Species Series"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_multispecies_recent)

# -----------------------------------------------------------------------------
# 4. Save Master Outputs
# -----------------------------------------------------------------------------

write.csv(master_all_species, file.path(output_dir, "master_weekly_all_species_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "multispecies_weekly_standardized_1998_2024.png"), p_multispecies_full, width = 13, height = 9)
ggsave(file.path(output_dir, "multispecies_weekly_standardized_2011_2024.png"), p_multispecies_recent, width = 12, height = 8)
