


library(tidyverse)
library(readxl)
library(scales)


# -----------------------------------------------------------------------------
# 1. Historical Landings (1976–2009) -> Convert Pounds to Grams-------
# -----------------------------------------------------------------------------

path <- "data_Lisa/Eulachon_data_Gustavson_2010_status_review.xlsx"

# Helper function to convert commas, dashes, and "Unknown*" text into clean numeric / NA
clean_numeric <- function(x) {
  x_clean <- gsub(",", "", as.character(x))
  x_clean <- ifelse(grepl("Unknown|—|-|N/A", x_clean, ignore.case = TRUE), NA_character_, x_clean)
  suppressWarnings(as.numeric(x_clean))
}

eulachon_count_50yr <- read_xlsx(path = path, sheet = 2, skip = 1) %>%
  rename(
    year                             = Year,
    eulachon_cr_pounds               = `Total landings\r\n(pounds)`,
    eulachon_cr_nfish_10.8_per_pound = `Number of fish at\r\n10.8 per pound`,
    eulachon_cr_nfish_12.3_per_pound = `Number of fish at\r\n12.3 per pound`
  ) %>%
  mutate(
    year = as.integer(clean_numeric(year)),
    across(starts_with("eulachon_cr_"), clean_numeric)
  ) %>%
  # Filter to the last 50 years of data
  filter(!is.na(year), year >= (1976))

cat("--- Eulachon Annual Data (Last 50 Years) ---\n")
print(eulachon_count_50yr, n = 50)


landings_grams_df <- eulachon_count_50yr %>%
  filter(year >= 1976, year <= 2009) %>%
  select(year, eulachon_cr_pounds) %>%
  mutate(
    # 1 lb = 453.59237 grams
    eulachon_grams = eulachon_cr_pounds * 453.59237,
    data_type      = "Commercial Landings"
  )

# -----------------------------------------------------------------------------
# 2. Modern Egg-Based SSB (2011–2024) -> Convert Pounds to Grams---------
# -----------------------------------------------------------------------------
# ssb=plankton outflow/(eggs/fish)/(fish/lb); eggs/fish = 16,293; fish/lb = 11.16 
#per email from Laura Heironimus WDFW Jun 10, 2025

library(tidyverse)

jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
lre_dat_yearly %>% select(year,csl_nonpup_total_emb)


ssb_grams_df <- lre_dat_yearly %>%
  filter(year >= 2011, year <= 2024) %>%
  select(year, eulachon_ssb_est) %>%
  mutate(
    # Formula confirmed output is in lbs -> convert to grams
    eulachon_grams = eulachon_ssb_est * 453.59237,
    data_type      = "Egg SSB"
  )

# -----------------------------------------------------------------------------
# 3. Combine into Single 50-Year Biomass Series (Grams)--------
# -----------------------------------------------------------------------------
eulachon_50yr_grams <- bind_rows(
  landings_grams_df %>% select(year, eulachon_grams, data_type),
  ssb_grams_df %>% select(year, eulachon_grams, data_type)
) %>%
  arrange(year)

# -----------------------------------------------------------------------------
# 4. Summary & Scale Comparison in Grams----------
# -----------------------------------------------------------------------------
cat("--- EULACHON ANNUAL BIOMASS (GRAMS) SUMMARY ---\n")
eulachon_50yr_grams %>%
  group_by(data_type) %>%
  summarise(
    years_count  = n(),
    mean_grams   = mean(eulachon_grams, na.rm = TRUE),
    median_grams = median(eulachon_grams, na.rm = TRUE),
    min_grams    = min(eulachon_grams, na.rm = TRUE),
    max_grams    = max(eulachon_grams, na.rm = TRUE)
  ) %>%
  print()

# Save corrected dataset to output folder
write.csv(
  eulachon_50yr_grams,
  file.path(output_dir, "eulachon_annual_biomass_grams_1976_2024.csv"),
  row.names = FALSE
)


library(tidyverse)
library(scales)

# -----------------------------------------------------------------------------
# 1. Combine Annual Eulachon Grams Time Series------
# -----------------------------------------------------------------------------
eulachon_annual_grams_plot <- eulachon_50yr_grams %>%
  mutate(
    # Convert to metric megagrams (10^6 g) or metric tons for cleaner axis display
    eulachon_metric_tons = eulachon_grams / 1e6
  )

# -----------------------------------------------------------------------------
# 2. Plot Both Eras (1976–2009 Landings vs. 2011–2024 SSB)---------
# -----------------------------------------------------------------------------
p_eulachon_history <- ggplot(eulachon_annual_grams_plot, aes(x = year, y = eulachon_grams, fill = data_type)) +
  geom_col(width = 0.75, alpha = 0.85) +
  scale_fill_manual(
    values = c("Commercial Landings" = "darkgoldenrod3", "Egg SSB" = "steelblue")
  ) +
  scale_x_continuous(breaks = seq(1976, 2024, by = 4)) +
  scale_y_continuous(labels = label_scientific()) +
  labs(
    title = "Eulachon Annual Biomass History (1976–2024)",
    subtitle = "Commercial Landings (1976–2009) and Modern Egg-Based SSB Estimates (2011–2024)",
    x = "Year",
    y = "Annual Biomass (Grams)",
    fill = "Data Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_eulachon_history)

# Save plot to output folder
ggsave(
  filename = file.path(output_dir, "eulachon_annual_grams_1976_2024.png"),
  plot     = p_eulachon_history,
  width    = 10, height = 5
)

#Plot in Pounds--------
# -----------------------------------------------------------------------------
# 1. Combine Annual Eulachon Series in Pounds (1976–2024)----
# -----------------------------------------------------------------------------
# A. Historical Landings (1976–2009)
landings_lbs_df <- eulachon_count_50yr %>%
  filter(year >= 1976, year <= 2009) %>%
  select(year, eulachon_lbs = eulachon_cr_pounds) %>%
  mutate(data_type = "Commercial Landings")

# B. Modern Egg-Based SSB (2011–2024)
ssb_lbs_df <- lre_dat_yearly %>%
  filter(year >= 2011, year <= 2024) %>%
  select(year, eulachon_lbs = eulachon_ssb_est) %>%
  mutate(data_type = "Egg SSB")

# -----------------------------------------------------------------------------
# 2. Merge, Interpolate 2010 Gap, and Finalize Annual Series-----
# -----------------------------------------------------------------------------
eulachon_annual_lbs_50yr <- bind_rows(landings_lbs_df, ssb_lbs_df) %>%
  # Force complete 1976–2024 grid (explicitly creates empty row for 2010)
  complete(year = 1976:2024) %>%
  arrange(year) %>%
  mutate(
    # Label 2010 explicitly
    data_type    = if_else(year == 2010, "Interpolated Gap (2010)", data_type),
    # Linearly interpolate 2010 between 2009 landings and 2011 SSB
    eulachon_lbs = zoo::na.approx(eulachon_lbs, na.rm = FALSE, rule = 2)
  )

# Save continuous annual pounds series to output directory
write.csv(
  eulachon_annual_lbs_50yr,
  file.path(output_dir, "eulachon_annual_lbs_1976_2024.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Summary & Visual Inspection of Annual Series in Pounds-----------
# -----------------------------------------------------------------------------
cat("--- EULACHON ANNUAL BIOMASS (POUNDS) SUMMARY (1976–2024) ---\n")
print(summary(eulachon_annual_lbs_50yr$eulachon_lbs))

p_lbs_history <- ggplot(eulachon_annual_lbs_50yr, aes(x = year, y = eulachon_lbs, fill = data_type)) +
  geom_col(width = 0.75, alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "Commercial Landings"     = "darkgoldenrod3",
      "Interpolated Gap (2010)" = "gray50",
      "Egg SSB"                 = "steelblue"
    )
  ) +
  scale_x_continuous(breaks = seq(1976, 2024, by = 4)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Eulachon Annual Biomass History in Pounds (1976–2024)",
    subtitle = "Commercial Landings (1976–2009), Interpolated Gap (2010), and Egg SSB (2011–2024)",
    x = "Year",
    y = "Annual Biomass (Pounds)",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_lbs_history)

# Save plot to outputs directory
ggsave(
  filename = file.path(output_dir, "eulachon_annual_lbs_1976_2024.png"),
  plot     = p_lbs_history,
  width    = 10, height = 5
)

#Scaling & Sine Wave Analysis of Commercial Landings (1976–2009)-------------
library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Filter & Scale Historical Commercial Landings Data (1976–2009)
# -----------------------------------------------------------------------------
landings_scaled <- eulachon_count_50yr %>%
  filter(year >= 1976, year <= 2009, !is.na(eulachon_cr_pounds)) %>%
  select(year, eulachon_cr_pounds) %>%
  mutate(
    # Z-score standardization (Mean = 0, SD = 1)
    landings_z = (eulachon_cr_pounds - mean(eulachon_cr_pounds)) / sd(eulachon_cr_pounds),
    
    # Min-Max Scaling (0 to 1 domain)
    landings_minmax = (eulachon_cr_pounds - min(eulachon_cr_pounds)) / 
      (max(eulachon_cr_pounds) - min(eulachon_cr_pounds)),
    
    # Continuous time index for multi-year harmonic modeling
    time_index = year - min(year) + 1
  )

# -----------------------------------------------------------------------------
# 2. Fit Multi-Year Harmonic Sine Model across Time Series
# -----------------------------------------------------------------------------
# Test for periodic cycle (e.g., 5-year, 8-year, or 10-year ocean productivity cycle)
# We test a 10-year harmonic cycle angular domain: 2*pi * time / 10

landings_scaled <- landings_scaled %>%
  mutate(
    # 10-year harmonic domain (radians)
    rad_10yr = 2 * pi * time_index / 10
  )

# Fit linear model with sine and cosine terms
sine_fit_10yr <- lm(landings_z ~ sin(rad_10yr) + cos(rad_10yr), data = landings_scaled)
summary(sine_fit_10yr)

# Add fitted sine values to dataframe
landings_scaled$sine_pred_z <- predict(sine_fit_10yr)

# -----------------------------------------------------------------------------
# 3. Plot Scaled Landings with Fitted Sine Wave Overlaid
# -----------------------------------------------------------------------------
p_landings_sine <- ggplot(landings_scaled, aes(x = year)) +
  # Raw Z-Score Points and Line
  geom_point(aes(y = landings_z), size = 2.5, color = "black") +
  geom_line(aes(y = landings_z), color = "gray40", linetype = "solid", linewidth = 0.8) +
  
  # Fitted Sine Wave Model (Red)
  geom_line(
    aes(y = sine_pred_z, color = "Fitted Sine Wave (10-Yr Cycle)"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  
  # Reference line at mean (0 Z-score)
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkred", alpha = 0.7) +
  
  scale_x_continuous(breaks = seq(1976, 2009, by = 4)) +
  scale_color_manual(values = c("Fitted Sine Wave (10-Yr Cycle)" = "firebrick")) +
  labs(
    title = "Historical Commercial Eulachon Landings (1976–2009)",
    subtitle = "Standardized Z-Scores (Mean = 0, SD = 1) with Overlaid Multi-Year Sine Wave Model",
    x = "Year",
    y = "Standardized Landings (Z-Score)",
    color = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_landings_sine)

# Save plot to output directory
ggsave(
  filename = file.path(output_dir, "eulachon_landings_scaled_sine_fit.png"),
  plot     = p_landings_sine,
  width    = 9, height = 5
)


#1993+ -----------------
library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Filter & Scale Landings Data Starting Post-Crash (1993–2009)
# -----------------------------------------------------------------------------
landings_post_crash <- eulachon_count_50yr %>%
  filter(year >= 1993, year <= 2009, !is.na(eulachon_cr_pounds)) %>%
  select(year, eulachon_cr_pounds) %>%
  mutate(
    # Z-score standardization calculated purely on the 1993-2009 window
    landings_z = (eulachon_cr_pounds - mean(eulachon_cr_pounds)) / sd(eulachon_cr_pounds),
    
    # Time index starting at 1 for 1993
    time_index = year - 1993 + 1
  )

# -----------------------------------------------------------------------------
# 2. Fit Multi-Year Sine Wave Model Across Post-Crash Window
# -----------------------------------------------------------------------------
# Test a harmonic cycle (~6-year ocean periodicity across the 17-year window)
landings_post_crash <- landings_post_crash %>%
  mutate(
    rad_cycle = 2 * pi * time_index / 6
  )

sine_post_crash_mod <- lm(landings_z ~ sin(rad_cycle) + cos(rad_cycle), data = landings_post_crash)
summary(sine_post_crash_mod)

# Add fitted sine predictions to dataframe
landings_post_crash$sine_pred_z <- predict(sine_post_crash_mod)

# -----------------------------------------------------------------------------
# 3. Plot Scaled Post-Crash Landings (1993–2009)
# -----------------------------------------------------------------------------
p_post_crash_sine <- ggplot(landings_post_crash, aes(x = year)) +
  # Observed Points and Trajectory
  geom_point(aes(y = landings_z), size = 3, color = "black") +
  geom_line(aes(y = landings_z), color = "gray30", linewidth = 0.9) +
  
  # Overlaid Sine Model (Red Dashed)
  geom_line(
    aes(y = sine_pred_z, color = "Fitted Sine Wave (~6-Yr Cycle)"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  
  # Mean reference line at 0 Z-score
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkred", alpha = 0.7) +
  
  scale_x_continuous(breaks = 1993:2009) +
  scale_color_manual(values = c("Fitted Sine Wave (~6-Yr Cycle)" = "firebrick")) +
  labs(
    title = "Post-Crash Commercial Eulachon Landings (1993–2009)",
    subtitle = "Standardized Z-Scores (Mean = 0, SD = 1) starting after the 1992–1993 step change",
    x = "Year",
    y = "Standardized Landings (Z-Score)",
    color = "Model Fit"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_post_crash_sine)

# Save plot to output directory
ggsave(
  filename = file.path(output_dir, "eulachon_landings_post_crash_1993_2009_sine.png"),
  plot     = p_post_crash_sine,
  width    = 9, height = 5
)

#Wrong sine wave period was used. refit------------
library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Scale Each Era Independently (Z-scores per era)
# -----------------------------------------------------------------------------
landings_era <- eulachon_count_50yr %>%
  filter(year >= 1993, year <= 2009, !is.na(eulachon_cr_pounds)) %>%
  select(year, pounds = eulachon_cr_pounds) %>%
  mutate(
    z_score   = (pounds - mean(pounds)) / sd(pounds),
    data_era  = "Post-Crash Landings (1993-2009)"
  )

ssb_era <- lre_dat_yearly %>%
  filter(year >= 2011, year <= 2024, !is.na(eulachon_ssb_est)) %>%
  select(year, pounds = eulachon_ssb_est) %>%
  mutate(
    z_score   = (pounds - mean(pounds)) / sd(pounds),
    data_era  = "Modern Egg SSB (2011-2024)"
  )

# Combine scaled eras
eul_two_eras_scaled <- bind_rows(landings_era, ssb_era) %>%
  arrange(year) %>%
  mutate(t = year - 1993 + 1) # Continuous year index starting at 1993

# -----------------------------------------------------------------------------
# 2. Fit Unconstrained Sine Model via Non-Linear Least Squares (nls)
# -----------------------------------------------------------------------------
# Model: z = A * sin(2*pi*t / T + phi)
# Initial guesses: Amplitude A ~ 1, Period T ~ 10-12 years, Phase phi ~ 0
sine_unconstrained <- nls(
  z_score ~ A * sin((2 * pi * t / T) + phi),
  data = eul_two_eras_scaled,
  start = list(A = 1, T = 11, phi = 0),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

summary(sine_unconstrained)

# Extract estimated period length
est_period <- coef(sine_unconstrained)["T"]
cat("--- ESTIMATED SINE PERIOD LENGTH (YEARS) ---\n")
print(round(est_period, 2))

# -----------------------------------------------------------------------------
# 3. Generate Continuous Sine Curve Predictions across 1993–2024
# -----------------------------------------------------------------------------
grid_pred <- tibble(
  year = seq(1993, 2024, by = 0.1)
) %>%
  mutate(
    t = year - 1993 + 1,
    sine_fit = predict(sine_unconstrained, newdata = tibble(t = t))
  )

# -----------------------------------------------------------------------------
# 4. Plot Both Scaled Eras with Unconstrained Sine Fit
# -----------------------------------------------------------------------------
p_two_eras_sine <- ggplot() +
  # Raw data points scaled to self
  geom_point(
    data = eul_two_eras_scaled,
    aes(x = year, y = z_score, color = data_era),
    size = 3
  ) +
  geom_line(
    data = eul_two_eras_scaled,
    aes(x = year, y = z_score, group = data_era, color = data_era),
    linewidth = 0.8, alpha = 0.7
  ) +
  # Fitted Sine Curve across full window
  geom_line(
    data = grid_pred,
    aes(x = year, y = sine_fit),
    color = "black", linetype = "dashed", linewidth = 1.1
  ) +
  # Reference zero line
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  
  scale_color_manual(
    values = c(
      "Post-Crash Landings (1993-2009)" = "darkgoldenrod3",
      "Modern Egg SSB (2011-2024)"     = "steelblue"
    )
  ) +
  scale_x_continuous(breaks = seq(1993, 2024, by = 2)) +
  labs(
    title = "Eulachon Biomass Oscillations (1993–2024)",
    subtitle = paste0(
      "Each era scaled to its own Z-score (Mean=0, SD=1). Black dashed line = fitted sine wave (Period = ", 
      round(est_period, 1), " years)"
    ),
    x = "Year",
    y = "Standardized Anomaly (Z-score)",
    color = "Dataset Era"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_two_eras_sine)

# Save plot
ggsave(
  filename = file.path(output_dir, "eulachon_independently_scaled_sine_alignment.png"),
  plot     = p_two_eras_sine,
  width    = 10, height = 5
)


#Detect Empirical Peaks in the Scaled Data------------
library(tidyverse)
library(slider)

# -----------------------------------------------------------------------------
# 1. Scale Each Era Independently to Z-Scores
# -----------------------------------------------------------------------------
landings_era <- eulachon_count_50yr %>%
  filter(year >= 1993, year <= 2009, !is.na(eulachon_cr_pounds)) %>%
  select(year, pounds = eulachon_cr_pounds) %>%
  mutate(
    z_score  = (pounds - mean(pounds)) / sd(pounds),
    data_era = "Post-Crash Landings (1993-2009)"
  )

ssb_era <- lre_dat_yearly %>%
  filter(year >= 2011, year <= 2024, !is.na(eulachon_ssb_est)) %>%
  select(year, pounds = eulachon_ssb_est) %>%
  mutate(
    z_score  = (pounds - mean(pounds)) / sd(pounds),
    data_era = "Modern Egg SSB (2011-2024)"
  )

eul_scaled_combined <- bind_rows(landings_era, ssb_era) %>%
  arrange(year)

# -----------------------------------------------------------------------------
# 2. Identify "2 out of 3 Positive Years" Clusters via slider::slide_dbl()
# -----------------------------------------------------------------------------
eul_peak_detection <- eul_scaled_combined %>%
  group_by(data_era) %>%
  mutate(
    is_pos = as.numeric(z_score > 0),
    # 3-year rolling sum centered on current year (1 before, 1 after)
    pos_count_3yr = slide_dbl(is_pos, sum, .before = 1, .after = 1, .complete = FALSE)
  ) %>%
  ungroup()

cat("--- YEARS WITH AT LEAST 2 OUT OF 3 POSITIVE YEARS ---\n")
eul_peak_detection %>%
  filter(pos_count_3yr >= 2) %>%
  select(year, z_score, pos_count_3yr, data_era) %>%
  print(n = 30)

# -----------------------------------------------------------------------------
# 3. Fit Unconstrained Sine/Cosine Harmonic Regression
# -----------------------------------------------------------------------------
grid_search_period <- tibble(T_candidate = seq(5, 18, by = 0.1)) %>%
  mutate(
    fit = map(T_candidate, function(T_val) {
      lm(
        z_score ~ sin(2 * pi * year / T_val) + cos(2 * pi * year / T_val),
        data = eul_scaled_combined
      )
    }),
    r_squared = map_dbl(fit, ~ summary(.x)$r.squared)
  )

# Extract best period T determined entirely by data fit
best_T     <- grid_search_period %>% arrange(desc(r_squared)) %>% slice(1) %>% pull(T_candidate)
best_model <- grid_search_period %>% arrange(desc(r_squared)) %>% slice(1) %>% pull(fit) %>% pluck(1)

cat("\n--- DATA-DRIVEN OPTIMAL SINE PERIOD ---\n")
cat("Best Period Length (T):", best_T, "years\n")
cat("R-Squared:", round(summary(best_model)$r.squared, 3), "\n")

# -----------------------------------------------------------------------------
# 4. Predict & Plot Unconstrained Data-Driven Fit
# -----------------------------------------------------------------------------
pred_grid <- tibble(year = seq(1993, 2024, by = 0.1))
pred_grid$sine_fit <- predict(best_model, newdata = pred_grid)

# Calculate empirical positive cluster highlights for plotting background
pos_clusters <- eul_peak_detection %>%
  filter(pos_count_3yr >= 2)

p_empirical_sine <- ggplot() +
  # Highlight 2-out-of-3 positive clusters with vertical shading
  geom_tile(
    data = pos_clusters,
    aes(x = year, y = 0, height = Inf),
    fill = "gold", alpha = 0.25
  ) +
  # Raw era points and lines
  geom_point(
    data = eul_scaled_combined,
    aes(x = year, y = z_score, color = data_era),
    size = 3
  ) +
  geom_line(
    data = eul_scaled_combined,
    aes(x = year, y = z_score, group = data_era, color = data_era),
    linewidth = 0.8, alpha = 0.7
  ) +
  # Data-driven Sine fit
  geom_line(
    data = pred_grid,
    aes(x = year, y = sine_fit),
    color = "black", linetype = "dashed", linewidth = 1.2
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  
  scale_color_manual(
    values = c(
      "Post-Crash Landings (1993-2009)" = "darkgoldenrod3",
      "Modern Egg SSB (2011-2024)"     = "steelblue"
    )
  ) +
  scale_x_continuous(breaks = seq(1993, 2024, by = 2)) +
  labs(
    title = "Eulachon Multi-Year Cycle (Empirical Data-Driven Fit)",
    subtitle = paste0(
      "Yellow bands = empirical 2-of-3 positive year clusters; Black dashed line = fitted harmonic sine (Period T = ", 
      best_T, " yrs)"
    ),
    x = "Year",
    y = "Standardized Anomaly (Z-score)",
    color = "Dataset Era"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_empirical_sine)

# Save plot
ggsave(
  filename = file.path(output_dir, "eulachon_empirical_data_driven_sine.png"),
  plot     = p_empirical_sine,
  width    = 10, height = 5
)
