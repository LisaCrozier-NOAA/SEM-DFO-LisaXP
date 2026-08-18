output_dir <- "outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}


library(tidyverse)
library(readxl)

#Part 1: read in data========
jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
head(lre_dat_yearly)
names(lre_dat_yearly)
#re-save in working directory
write.csv(
  lre_dat_yearly,
  file.path(output_dir, "jake_lre_dat_yearly.csv"),
  row.names = FALSE
)



week_all<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),row.names=NULL)
head(week_all)
write.csv(
  week_all,
  file.path(output_dir, "jake_week_all.csv"),
  row.names = FALSE
)

#week_final_scaled<-read.csv("csl_week_filled_in_2011_2024.csv",row.names = NULL)

# -----------------------------------------------------------------------------
# 1. Read, Clean, and Format Historical Annual Eulachon Data
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


#Part 1: scale 2 datasets, fit sine curve================
          # -----------------------------------------------------------------------------
          # 1. Scale Each Era Independently to Z-Scores-------
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
          # 2. Identify "2 out of 3 Positive Years" Clusters via slider::slide_dbl()------
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
          # 3. Fit Unconstrained Sine/Cosine Harmonic Regression---------
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
          # 4. Predict & Plot Unconstrained Data-Driven Fit-------
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
          

#Part 2: Mash Up =================================
# -----------------------------------------------------------------------------
# 1. Extract Modern SSB Scaling Parameters (Mean & SD)------
# -----------------------------------------------------------------------------
ssb_mean <- mean(ssb_era$pounds, na.rm = TRUE)
ssb_sd   <- sd(ssb_era$pounds, na.rm = TRUE)

cat("--- MODERN EGG SSB SCALING PARAMETERS (POUNDS) ---\n")
cat("Mean SSB (2011-2024):", scales::comma(ssb_mean), "lbs\n")
cat("SD SSB   (2011-2024):", scales::comma(ssb_sd), "lbs\n\n")

# -----------------------------------------------------------------------------
# 2. Build Unified Master Table (1993–2024)--------
# -----------------------------------------------------------------------------
# We expand the sequence to cover 1993:2024 (including the 2010 gap)
eulachon_index_1993_2024 <- tibble(year = 1993:2024) %>%
  left_join(eul_scaled_combined, by = "year") %>%
  mutate(
    # A. Unified Z-score (Mash-up of separately scaled series)
    z_score_mashup = z_score,
    
    # B. Back-transform Z-scores into pounds using the Modern SSB scale
    # Formula: Pounds = (Z_score * Modern_SSB_SD) + Modern_SSB_Mean
    eulachon_lbs_reconstructed = (z_score_mashup * ssb_sd) + ssb_mean,
    
    # Label source for clarity
    index_source = case_when(
      year >= 1993 & year <= 2009 ~ "Back-Transformed Landings (1993-2009)",
      year == 2010                ~ "Gap Year (2010)",
      year >= 2011 & year <= 2024 ~ "Observed Egg SSB (2011-2024)"
    )
  )

# -----------------------------------------------------------------------------
# 3. Predict & Back-Transform Sine Estimate Across Full Timeline--------
# -----------------------------------------------------------------------------
# Predict Z-score sine wave across 1993-2024 (including 2010)
sine_preds_z <- predict(best_model, newdata = tibble(year = 1993:2024))

# Back-transform sine Z-scores into pounds using the modern SSB scale
eulachon_index_1993_2024 <- eulachon_index_1993_2024 %>%
  mutate(
    sine_z_pred   = sine_preds_z,
    sine_lbs_pred = pmax(0, (sine_z_pred * ssb_sd) + ssb_mean) # Truncate negative tail at 0
  ) %>%
  select(
    year, 
    index_source,
    z_score_mashup,
    eulachon_lbs_reconstructed,
    sine_z_pred,
    sine_lbs_pred
  )

# View head and tail of new master table
cat("--- MASTER EULACHON INDEX IN POUNDS (1993–2024) ---\n")
print(eulachon_index_1993_2024, n = 32)

# Save to output directory
write.csv(
  eulachon_index_1993_2024,
  file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Plot Visual Comparison: Reconstructed Raw Index vs. Sine Curve in Pounds--------
# -----------------------------------------------------------------------------
p_master_index <- ggplot(eulachon_index_1993_2024, aes(x = year)) +
  # Back-transformed raw data points
  geom_point(
    aes(y = eulachon_lbs_reconstructed, color = index_source),
    size = 3.2, na.rm = TRUE
  ) +
  geom_line(
    aes(y = eulachon_lbs_reconstructed),
    color = "gray30", linewidth = 0.8, linetype = "dotted", na.rm = TRUE
  ) +
  # Back-transformed Sine Wave Estimate (Black Dashed)
  geom_line(
    aes(y = sine_lbs_pred, linetype = "Sine Model Estimate"),
    color = "black", linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Back-Transformed Landings (1993-2009)" = "darkgoldenrod3",
      "Gap Year (2010)"                        = "transparent",
      "Observed Egg SSB (2011-2024)"          = "steelblue"
    )
  ) +
  scale_linetype_manual(values = c("Sine Model Estimate" = "dashed")) +
  scale_x_continuous(breaks = seq(1993, 2024, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Reconstructed Eulachon Annual Index in Pounds (1993–2024)",
    subtitle = "1993–2009 landings back-transformed to modern Egg SSB scale; Sine wave provides smooth trough/peak baseline",
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

# Save plot to outputs directory
ggsave(
  filename = file.path(output_dir, "eulachon_master_reconstructed_index_lbs.png"),
  plot     = p_master_index,
  width    = 10, height = 5.5
)
