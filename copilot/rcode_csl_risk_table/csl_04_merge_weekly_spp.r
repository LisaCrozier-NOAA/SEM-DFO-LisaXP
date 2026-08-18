



#Start processing Chinook---------

#process chin.r
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


#End processing Chinook-----------------



#start processing shad------------
#process shar.r
# ==============================================================================
# Dedicated Script: Clean & Interpolate Daily Bonneville Shad Passage
# ==============================================================================

library(tidyverse)
library(lubridate)
library(zoo) # For linear interpolation (na.approx)

# -----------------------------------------------------------------------------
# 1. Read Raw Wide File & Drop Feb 29
# -----------------------------------------------------------------------------

shad_raw <- read.csv("data_Lisa/SHAD_bonn_19382026.csv", row.names = NULL, stringsAsFactors = FALSE)

# Filter out Feb 29 upfront
shad_clean_wide <- shad_raw %>%
  filter(!grepl("29-Feb|29.Feb|Feb.29", mm.dd, ignore.case = TRUE))

# -----------------------------------------------------------------------------
# 2. Pivot to Long Format & Standardize Column Names
# -----------------------------------------------------------------------------

shad_long <- shad_clean_wide %>%
  pivot_longer(
    cols      = -mm.dd,
    names_to  = "year_col",
    values_to = "raw_count"
  ) %>%
  mutate(
    # Parse 4-digit year from column names (e.g. "X2024.BON..." -> 2024)
    year = parse_number(year_col),
    
    # Strip spaces and commas from count values
    count_clean = gsub("[^0-9.]", "", trimws(raw_count)),
    
    # Convert empty strings or non-numerics to NA for interpolation
    count_num   = if_else(count_clean == "", NA_real_, as.numeric(count_clean))
  ) %>%
  filter(year >= 1998, year <= 2024)

# -----------------------------------------------------------------------------
# 3. Construct Safe Dates & Linearly Interpolate Missing Counts
# -----------------------------------------------------------------------------

shad_interpolated <- shad_long %>%
  # Force standard year-day string format
  mutate(
    date_str   = paste(mm.dd, year, sep = "-"),
    date_clean = dmy(date_str)
  ) %>%
  # Drop any residual unparseable dates
  filter(!is.na(date_clean)) %>%
  arrange(year, date_clean) %>%
  group_by(year) %>%
  mutate(
    # Linear interpolation for interior NAs; rule=2 carries boundary values
    count_interp = na.approx(count_num, x = date_clean, na.rm = FALSE, rule = 2),
    # If a full year is NA, default to 0
    count_interp = replace_na(count_interp, 0),
    
    week_num = week(date_clean)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 4. Aggregate Weekly Sums (Weeks 10–26)
# -----------------------------------------------------------------------------

shad_weekly_final <- shad_interpolated %>%
  filter(week_num >= 10, week_num <= 26) %>%
  group_by(year, week = week_num) %>%
  summarise(
    Shad_weekly_count = sum(count_interp, na.rm = TRUE),
    Shad_weekly_mean  = mean(count_interp, na.rm = TRUE),
    .groups           = "drop"
  )

# Ensure full grid for weeks 10–26 across all years
shad_weekly_grid <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(shad_weekly_final, by = c("year", "week")) %>%
  mutate(
    Shad_weekly_count = replace_na(Shad_weekly_count, 0),
    Shad_weekly_mean  = replace_na(Shad_weekly_mean, 0)
  )

cat("--- Processed Shad Dataset Summary (1998–2024) ---\n")
print(head(shad_weekly_grid))
print(summary(shad_weekly_grid$Shad_weekly_count))


write.csv(shad_weekly_grid, file.path(output_dir, "bonn_shad_weekly_interpolated_1998_2024.csv"), row.names = FALSE)

#end processing shad-------------

#combine all spp in master dataset --------
#combine_weekly_counts.r-----
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

#lag for estuary--------
# Helper function to avoid -Inf warnings on all-NA years
safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA_real_)
  max(x_clean)
}

master_estuary <- master_all %>%
  select(year, week, csl_final, eulachon_final, Chin_weekly_count, Shad_weekly_count) %>%
  group_by(year) %>%
  mutate(
    # 2-Week Shift for Estuary Arrival Timing
    chin_estuary_count = lead(Chin_weekly_count, 2, default = 0),
    shad_estuary_count = lead(Shad_weekly_count, 2, default = 0)
  ) %>%
  # Define within-year proportions to identify main run timing
  mutate(
    chin_max  = safe_max(chin_estuary_count),
    chin_prop = if_else(is.na(chin_max) | chin_max <= 0, 0, chin_estuary_count / chin_max),
    
    # Flag weeks where Chinook presence is >= 10% of seasonal peak (Active Run Window)
    is_chinook_active = chin_prop >= 0.10
  ) %>%
  ungroup()
# -----------------------------------------------------------------------------
# 4. Save Master Outputs----------
# -----------------------------------------------------------------------------

write.csv(master_all_species, file.path(output_dir, "master_weekly_all_species_1998_2024.csv"), row.names = FALSE)
write.csv(master_estuary, file.path(output_dir, "master_estuary_2wklagweekly_all_species_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "multispecies_weekly_standardized_1998_2024.png"), p_multispecies_full, width = 13, height = 9)
ggsave(file.path(output_dir, "multispecies_weekly_standardized_2011_2024.png"), p_multispecies_recent, width = 12, height = 8)

#redo master files--------
# ==============================================================================
# Script: Merge All Species & Create 2-Week Estuary Lagged Dataset
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# Helper function to prevent max() from throwing -Inf warnings on all-NA years
safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA_real_)
  max(x_clean, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Read Clean Individual Datasets
# -----------------------------------------------------------------------------

eulachon_df <- read.csv(file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"))
csl_df      <- read.csv(file.path(output_dir, "csl_reconstructed_weekly_1998_2024.csv"))
chinook_df  <- read.csv(file.path(output_dir, "bonn_chinook_weekly_interpolated_1998_2024.csv"))
shad_df     <- read.csv(file.path(output_dir, "bonn_shad_weekly_interpolated_1998_2024.csv"), row.names = NULL)

# -----------------------------------------------------------------------------
# 2. Build Master Combined Grid at Bonneville Dam Timing
# -----------------------------------------------------------------------------

master_all_species <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(eulachon_df %>% select(year, week, eulachon_final), by = c("year", "week")) %>%
  left_join(csl_df %>% select(year, week, csl_final), by = c("year", "week")) %>%
  left_join(chinook_df %>% select(year, week, Chin_weekly_count), by = c("year", "week")) %>%
  left_join(shad_df %>% select(year, week, Shad_weekly_count), by = c("year", "week")) %>%
  mutate(
    Chin_weekly_count = replace_na(Chin_weekly_count, 0),
    Shad_weekly_count = replace_na(Shad_weekly_count, 0)
  ) %>%
  group_by(year) %>%
  mutate(
    # Standardized 0-1 proportions at dam timing
    eul_max  = safe_max(eulachon_final),
    eul_prop = if_else(is.na(eul_max) | eul_max <= 0, 0, eulachon_final / eul_max),
    
    csl_max  = safe_max(csl_final),
    csl_prop = if_else(is.na(csl_max) | csl_max <= 0, 0, csl_final / csl_max),
    
    chin_max  = safe_max(Chin_weekly_count),
    chin_prop = if_else(is.na(chin_max) | chin_max <= 0, 0, Chin_weekly_count / chin_max),
    
    shad_max  = safe_max(Shad_weekly_count),
    shad_prop = if_else(is.na(shad_max) | shad_max <= 0, 0, Shad_weekly_count / shad_max)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. Create Estuary-Shifted Dataset (2-Week Advance for Chinook & Shad)
# -----------------------------------------------------------------------------

# Bonneville Week 12 passage = Estuary Week 10 presence
master_estuary <- master_all_species %>%
  select(year, week, csl_final, eulachon_final, Chin_weekly_count, Shad_weekly_count) %>%
  group_by(year) %>%
  mutate(
    # Pull passage data 2 weeks earlier to represent Estuary timing
    chin_estuary_count = lead(Chin_weekly_count, 2, default = 0),
    shad_estuary_count = lead(Shad_weekly_count, 2, default = 0),
    
    # Recalculate 0-1 proportions for Estuary-shifted run curves
    chin_est_max  = safe_max(chin_estuary_count),
    chin_est_prop = if_else(is.na(chin_est_max) | chin_est_max <= 0, 0, chin_estuary_count / chin_est_max),
    
    shad_est_max  = safe_max(shad_estuary_count),
    shad_est_prop = if_else(is.na(shad_est_max) | shad_est_max <= 0, 0, shad_estuary_count / shad_est_max),
    
    # Flag active Chinook window in the estuary (>= 10% of peak estuary arrival)
    is_chinook_active = chin_est_prop >= 0.10
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 4. Save Final CSV Files & Figures
# -----------------------------------------------------------------------------

write.csv(master_all_species, file.path(output_dir, "master_weekly_all_species_1998_2024.csv"), row.names = FALSE)
write.csv(master_estuary, file.path(output_dir, "master_estuary_2wklagweekly_all_species_1998_2024.csv"), row.names = FALSE)

cat("--- Master Datasets Successfully Created & Saved ---\n")
print(head(master_estuary))
