

#RESULTS:
#eulachon_index_1993_2024<-read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"),row.names = FALSE)
      #  names(eulachon_index_1993_2024)
      # [1] "year"                       "index_source"               "z_score_mashup"            
      # [4] "eulachon_lbs_reconstructed" "sine_z_pred"                "sine_lbs_pred"             
      
      #plots of z-score and raw reconstructed pounds
      # ggsave(
      #   filename = file.path(output_dir, "eulachon_empirical_data_driven_sine.png"),
      #   plot     = p_empirical_sine,
      #   width    = 10, height = 5
      # )
      # ggsave(
      #   filename = file.path(output_dir, "eulachon_master_reconstructed_index_lbs.png"),
      #   plot     = p_master_index,
      #   width    = 10, height = 5.5)

# print(p_empirical_sine)
# print(p_master_index)


# -----------------------------------------------------------------------------
# Setup & Package Loading--------
# -----------------------------------------------------------------------------
output_dir <- "outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

library(tidyverse)
library(readxl)
library(slider)   # For slide_dbl rolling windows
library(scales)   # For comma number formatting
library(zoo)      # For na.approx gap interpolation

# -----------------------------------------------------------------------------
# Part 1: Read & Clean Raw Data--------
# -----------------------------------------------------------------------------
# Part 1A: Jake's yearly (1976-2024) and weekly (2011-2024) LRE Data------
# jake_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
# load(paste0(jake_path, "lre_dat_yearly.RData"), verbose = TRUE)
# 
# write.csv(
#   lre_dat_yearly,
#   file.path(output_dir, "jake_lre_dat_yearly.csv"),
#   row.names = FALSE
# )
# week_all <- read.csv(
#   paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
#   row.names = NULL
# )
# 
# write.csv(
#   week_all,
#   file.path(output_dir, "jake_week_all.csv"),
#   row.names = FALSE
# )
# Part 1B. Read & Clean Historical Eulachon Excel Data--------
# path <- "data_Lisa/Eulachon_data_Gustavson_2010_status_review.xlsx"
# 
# clean_numeric <- function(x) {
#   x_clean <- gsub(",", "", as.character(x))
#   x_clean <- ifelse(grepl("Unknown|—|-|N/A", x_clean, ignore.case = TRUE), NA_character_, x_clean)
#   suppressWarnings(as.numeric(x_clean))
# }
# 
# eulachon_count_50yr <- read_xlsx(path = path, sheet = 2, skip = 1) %>%
#   rename(
#     year                             = Year,
#     eulachon_cr_pounds               = `Total landings\r\n(pounds)`,
#     eulachon_cr_nfish_10.8_per_pound = `Number of fish at\r\n10.8 per pound`,
#     eulachon_cr_nfish_12.3_per_pound = `Number of fish at\r\n12.3 per pound`
#   ) %>%
#   mutate(
#     year = as.integer(clean_numeric(year)),
#     across(starts_with("eulachon_cr_"), clean_numeric)
#   ) %>%
#   filter(!is.na(year), year >= 1976)
# 
# cat("--- Eulachon Annual Data (1976–2024) ---\n")
# print(eulachon_count_50yr, n = 50)
# 
# write.csv(
#   eulachon_count_50yr,
#   file.path(output_dir, "status_review_eulachon_count_50yr.csv"),
#   row.names = FALSE
# )

# Part 1C. Read cleaned eulachon data-------
lre_dat_yearly<- read.csv(file.path(output_dir, "jake_lre_dat_yearly.csv"),row.names = NULL)
week_all<- read.csv(file.path(output_dir, "jake_week_all.csv"),row.names = NULL)
eulachon_count_50yr<- read.csv(file.path(output_dir, "status_review_eulachon_count_50yr.csv"),row.names = NULL)
#per Laura's email, "don't use 2020 data -- season was curtailed")
lre_dat_yearly <- lre_dat_yearly %>%
  mutate(eulachon_ssb_est_corrected=case_when(
    year==2020 ~ NA,
    year==2010 ~ NA,
    TRUE ~ eulachon_ssb_est))

lre_dat_yearly %>% select(year,eulachon_ssb_est,eulachon_ssb_est_corrected)
# -----------------------------------------------------------------------------
# Part 2: Scale 2 Datasets & Fit Sine Curve--------
# -----------------------------------------------------------------------------
landings_era <- eulachon_count_50yr %>%
  filter(year >= 1993, year <= 2009, !is.na(eulachon_cr_pounds)) %>%
  select(year, pounds = eulachon_cr_pounds) %>%
  mutate(
    z_score  = (pounds - mean(pounds)) / sd(pounds),
    data_era = "Post-Crash Landings (1993-2009)"
  )

ssb_era <- lre_dat_yearly %>%
  filter(year >= 2011, year <= 2024, !is.na(eulachon_ssb_est_corrected)) %>%
  select(year, pounds = eulachon_ssb_est_corrected) %>%
  mutate(
    z_score  = (pounds - mean(pounds)) / sd(pounds),
    data_era = "Modern Egg SSB (2011-2024)"
  )

eul_scaled_combined <- bind_rows(landings_era, ssb_era) %>%
  arrange(year)

# Detect 2-out-of-3 positive years across the timeline
eul_peak_detection <- eul_scaled_combined %>%
  mutate(
    is_pos = as.numeric(z_score > 0),
    pos_count_3yr = slide_dbl(is_pos, sum, .before = 1, .after = 1, .complete = FALSE)
  )

cat("--- YEARS WITH AT LEAST 2 OUT OF 3 POSITIVE YEARS ---\n")
eul_peak_detection %>%
  filter(pos_count_3yr >= 2) %>%
  select(year, z_score, pos_count_3yr, data_era)
#  print(n = 30)

# Unconstrained Grid Search for Harmonic Sine Wave Period Length T
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

best_T     <- grid_search_period %>% arrange(desc(r_squared)) %>% slice(1) %>% pull(T_candidate)
best_model <- grid_search_period %>% arrange(desc(r_squared)) %>% slice(1) %>% pull(fit) %>% pluck(1)

cat("\n--- DATA-DRIVEN OPTIMAL SINE PERIOD ---\n")
cat("Best Period Length (T):", best_T, "years\n")
cat("R-Squared:", round(summary(best_model)$r.squared, 3), "\n")

# Predict & Plot Empirical Sine Fit
pred_grid <- tibble(year = seq(1993, 2024, by = 0.1))
pred_grid$sine_fit <- predict(best_model, newdata = pred_grid)

pos_clusters <- eul_peak_detection %>% filter(pos_count_3yr >= 2)

p_empirical_sine <- ggplot() +
  geom_tile(
    data = pos_clusters,
    aes(x = year, y = 0, height = Inf),
    fill = "gold", alpha = 0.25
  ) +
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

ggsave(
  filename = file.path(output_dir, "eulachon_empirical_data_driven_sine.png"),
  plot     = p_empirical_sine,
  width    = 10, height = 5
)

# -----------------------------------------------------------------------------
# Part 3: Mash Up & Interpolate 2010 & 2020 Gap--------
# -----------------------------------------------------------------------------
ssb_mean <- mean(ssb_era$pounds, na.rm = TRUE)
ssb_sd   <- sd(ssb_era$pounds, na.rm = TRUE)

cat("--- MODERN EGG SSB SCALING PARAMETERS (POUNDS) ---\n")
cat("Mean SSB (2011-2024):", scales::comma(ssb_mean), "lbs\n")
cat("SD SSB   (2011-2024):", scales::comma(ssb_sd), "lbs\n\n")

# Build complete 1993-2024 timeline and interpolate missing 2010 and 2020 gap
eulachon_index_1993_2024 <- tibble(year = 1993:2024) %>%
  left_join(eul_scaled_combined, by = "year") %>%
  mutate(
    # Back-transform Z-scores into pounds using Modern SSB scale
    eulachon_lbs_raw = (z_score * ssb_sd) + ssb_mean,
    
    # Linearly interpolate 2010 gap across z-scores and pounds
    z_score_mashup             = zoo::na.approx(z_score, na.rm = FALSE),
    eulachon_lbs_reconstructed = zoo::na.approx(eulachon_lbs_raw, na.rm = FALSE),
    
    # Label source for clarity
    index_source = case_when(
      year >= 1993 & year <= 2009 ~ "Back-Transformed Landings (1993-2009)",
      year == 2010                ~ "Interpolated Gap (2010)",
      year == 2020                ~ "Interpolated Gap (2020)",
      year >= 2011 & year <= 2024 ~ "Observed Egg SSB (2011-2024)"
    )
  )

print(eulachon_index_1993_2024 ,n=Inf)

# Predict sine wave across entire timeline
sine_preds_z <- predict(best_model, newdata = tibble(year = 1993:2024))

eulachon_index_1993_2024 <- eulachon_index_1993_2024 %>%
  mutate(
    sine_z_pred   = sine_preds_z,
    sine_lbs_pred = pmax(0, (sine_z_pred * ssb_sd) + ssb_mean)
  ) %>%
  select(
    year, 
    index_source,
    z_score_mashup,
    eulachon_lbs_reconstructed,
    sine_z_pred,
    sine_lbs_pred
  )

# Verify zero missing values remain
cat("--- CHECK FOR NA VALUES IN MASTER INDEX ---\n")
print(colSums(is.na(eulachon_index_1993_2024)))

cat("\n--- MASTER EULACHON INDEX IN POUNDS (1993–2024) ---\n")
print(eulachon_index_1993_2024, n = 32)

write.csv(
  eulachon_index_1993_2024,
  file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 4: Plot Reconstructed Master Index--------
# -----------------------------------------------------------------------------
p_master_index <- ggplot(eulachon_index_1993_2024, aes(x = year)) +
  geom_point(
    aes(y = eulachon_lbs_reconstructed, color = index_source),
    size = 3.2
  ) +
  geom_line(
    aes(y = eulachon_lbs_reconstructed),
    color = "gray30", linewidth = 0.8, linetype = "dotted"
  ) +
  geom_line(
    aes(y = sine_lbs_pred, linetype = "Sine Model Estimate"),
    color = "black", linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Back-Transformed Landings (1993-2009)" = "darkgoldenrod3",
      "Interpolated Gap (2010)"                = "firebrick",
      "Observed Egg SSB (2011-2024)"          = "steelblue"
    )
  ) +
  scale_linetype_manual(values = c("Sine Model Estimate" = "dashed")) +
  scale_x_continuous(breaks = seq(1993, 2024, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Reconstructed Eulachon Annual Index in Pounds (1993–2024)",
    subtitle = "1993–2009 landings back-transformed; 2010 linearly interpolated; Sine wave provides smooth trough baseline",
    x = "Year",
    y = "Eulachon Biomass Index (Pounds)",
    color = "Data Source",
    linetype = "Model Signal"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_master_index)

# -----------------------------------------------------------------------------
# Part 5: Save reconstructed index and plot --------
# -----------------------------------------------------------------------------
ggsave(
  filename = file.path(output_dir, "eulachon_master_reconstructed_index_lbs.png"),
  plot     = p_master_index,
  width    = 10, height = 5.5
)