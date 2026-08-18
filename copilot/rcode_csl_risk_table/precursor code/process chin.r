# ==============================================================================
# Dedicated Script: Process & Interpolate Bonneville Daily Chinook Passage
# ==============================================================================

library(tidyverse)
library(readxl)
library(lubridate)
library(zoo) # For linear interpolation

output_dir <- "copilot/outputs_csl_cr"

# -----------------------------------------------------------------------------
# 1. Read Raw Excel File & Clean Chinook Daily Data
# -----------------------------------------------------------------------------

excel_path <- "data_Lisa/Adult_BONpassage_1976_2024.xlsx"
bonn_daily <- read_excel(excel_path, sheet = 1)

chin_daily_clean <- bonn_daily %>%
  filter(parameter == "Chin") %>%
  mutate(
    date_clean  = as.Date(`mm-dd`),
    week_num    = week(date_clean),
    
    # Strip non-numeric characters (commas, spaces)
    count_clean = gsub("[^0-9.]", "", trimws(count)),
    
    # Convert empty strings or non-numerics to true NAs for interpolation
    count_num   = if_else(count_clean == "", NA_real_, suppressWarnings(as.numeric(count_clean)))
  ) %>%
  filter(
    year >= 1998, year <= 2024,
    week_num >= 10, week_num <= 26
  )

# -----------------------------------------------------------------------------
# 2. Linearly Interpolate Missing Daily Counts Within Each Year
# -----------------------------------------------------------------------------

chin_daily_interpolated <- chin_daily_clean %>%
  arrange(year, date_clean) %>%
  group_by(year) %>%
  mutate(
    # Interpolate interior NAs; rule=2 carries boundary values if endpoints are NA
    count_interp = na.approx(count_num, x = date_clean, na.rm = FALSE, rule = 2),
    # Default to 0 only if an entire year has no data
    count_interp = replace_na(count_interp, 0)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. Aggregate Daily Interpolated Counts to Weekly Sums & Means
# -----------------------------------------------------------------------------

chin_weekly_sum <- chin_daily_interpolated %>%
  group_by(year, week = week_num) %>%
  summarise(
    Chin_weekly_count = sum(count_interp, na.rm = TRUE),
    Chin_weekly_mean  = mean(count_interp, na.rm = TRUE),
    .groups           = "drop"
  )

# -----------------------------------------------------------------------------
# 4. Construct Full Grid (Weeks 10–26, 1998–2024)
# -----------------------------------------------------------------------------

chin_weekly_grid <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(chin_weekly_sum, by = c("year", "week")) %>%
  mutate(
    Chin_weekly_count = replace_na(Chin_weekly_count, 0),
    Chin_weekly_mean  = replace_na(Chin_weekly_mean, 0)
  )

# -----------------------------------------------------------------------------
# 5. Diagnostic Check
# -----------------------------------------------------------------------------

cat("--- Interpolated Chinook Weekly Dataset Summary (1998–2024) ---\n")
print(head(chin_weekly_grid, 10))
cat("\nSummary of Weekly Interpolated Counts:\n")
print(summary(chin_weekly_grid$Chin_weekly_count))

# Save standalone clean interpolated Chinook file
write.csv(chin_weekly_grid, file.path(output_dir, "bonn_chinook_weekly_interpolated_1998_2024.csv"), row.names = FALSE)



# ==============================================================================
# Script: Faceted Diagnostic Plots for Interpolated Weekly Chinook Passage
# Matching exact project plotting style (1998–2024, Weeks 10–26)
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# -----------------------------------------------------------------------------
# 1. Read Clean Interpolated Chinook Data & Prep 0–1 Standardized Proportions
# -----------------------------------------------------------------------------

chin_weekly_grid <- read.csv(file.path(output_dir, "bonn_chinook_weekly_interpolated_1998_2024.csv"))

chin_plot_prep <- chin_weekly_grid %>%
  group_by(year) %>%
  mutate(
    chin_peak = max(Chin_weekly_count, na.rm = TRUE),
    chin_prop = if_else(chin_peak == 0 | is.na(chin_peak), 0, Chin_weekly_count / chin_peak)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2. Plot 1: Faceted Absolute Weekly Chinook Counts (1998–2024)
# -----------------------------------------------------------------------------

p_chin_absolute <- ggplot(chin_plot_prep, aes(x = week, y = Chin_weekly_count)) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  geom_point(color = "firebrick", size = 1.5) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Weekly Interpolated Chinook Salmon Passage (1998–2024, Weeks 10–26)",
    subtitle = "Bonneville Dam passage totals; faceted by year with free y-axes",
    x = "Week (10–26)",
    y = "Weekly Chinook Passage Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_chin_absolute)

# -----------------------------------------------------------------------------
# 3. Plot 2: Faceted Standardized Seasonal Shape (0 to 1 Scale)
# -----------------------------------------------------------------------------

p_chin_standardized <- ggplot(chin_plot_prep, aes(x = week, y = chin_prop)) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  geom_point(color = "firebrick", size = 1.5) +
  facet_wrap(~ year, ncol = 4) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
  labs(
    title = "Standardized Weekly Chinook Seasonal Profile (1998–2024, Weeks 10–26)",
    subtitle = "Proportion of annual peak weekly count within each year",
    x = "Week (10–26)",
    y = "Proportion of Annual Peak"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_chin_standardized)

# -----------------------------------------------------------------------------
# 4. Save Image Outputs
# -----------------------------------------------------------------------------

ggsave(file.path(output_dir, "chinook_weekly_absolute_faceted_1998_2024.png"), p_chin_absolute, width = 12, height = 8)
ggsave(file.path(output_dir, "chinook_weekly_standardized_faceted_1998_2024.png"), p_chin_standardized, width = 12, height = 8)