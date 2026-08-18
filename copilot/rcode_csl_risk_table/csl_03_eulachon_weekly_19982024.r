


#RESULTS
# write.csv(eulachon_full_reconstructed, file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"), row.names = FALSE)
# ggsave(file.path(output_dir, "eulachon_weekly_reconstructed_faceted_1998_2024.png"), p_eul_full_faceted, width = 12, height = 8)
# ggsave(file.path(output_dir, "eulachon_weekly_reconstructed_continuous_1998_2024.png"), p_eul_full_continuous, width = 12, height = 6)

# -----------------------------------------------------------------------------
# 1. Read Data and Coerce Factor Column Explicitly
# -----------------------------------------------------------------------------
week_final_scaled <- read.csv(file.path(output_dir, "csl_week_filled_in_2011_2024.csv"), row.names = NULL) %>%
  mutate(
    year_factor = factor(year) # <-- FIX 1: Explicitly create factor from year
  )
week_raw<-week_final_scaled %>% select(year,week,month,year_factor,date,Spring_Achin_bonn_pass,csl_nonpup_total_emb,eulachon_ssb_4week_est,eulachon_input)
head(week_raw)


# -----------------------------------------------------------------------------
# 2. Standardize Weekly Data Within Each Year (0 to 1 Scale)
# -----------------------------------------------------------------------------
week_standardized_prep <- week_raw %>%
  group_by(year) %>%
  mutate(
    year_peak = max(eulachon_input, na.rm = TRUE),
    year_peak_clean = if_else(is.infinite(year_peak) | year_peak <= 0, NA_real_, year_peak),
    eul_prop = eulachon_input / year_peak_clean
  ) %>%
  ungroup()

# Extract observed annual peaks for training years (2011–2024)
observed_peaks_training <- week_standardized_prep %>%
  group_by(year) %>%
  summarise(observed_peak = unique(year_peak_clean), .groups = "drop")

# Join observed peaks with ALL years in eulachon_annual (1998–2024)
annual_peaks_full <- eulachon_annual %>%
  left_join(observed_peaks_training, by = "year")

# -----------------------------------------------------------------------------
# 3. Fit Peak Model & Predict Peak Scalars for ALL Years (1998–2024)
# -----------------------------------------------------------------------------
# Fit linear model on training years where observed_peak is available
peak_model <- lm(observed_peak ~ eulachon_lbs_reconstructed, 
                 data = annual_peaks_full %>% filter(!is.na(observed_peak)))

summary(peak_model)

# Predict expected peak height for EVERY year from 1998 to 2024
annual_peaks_calibrated <- annual_peaks_full %>%
  mutate(
    predicted_peak = predict(peak_model, newdata = .),
    predicted_peak = pmax(0, predicted_peak), # Ensure non-negative
    
    # Use actual observed peak for 2011–2024 if present; use predicted peak for 1998–2010
    final_peak_scalar = if_else(!is.na(observed_peak), observed_peak, predicted_peak)
  )

# Master Grid across ALL years (1998–2024, excl. 2020)
master_grid_eul <- expand_grid(
  year = setdiff(1998:2024, 2020),
  week = 10:26
) %>%
  left_join(week_standardized_prep, by = c("year", "week")) %>%
  left_join(annual_peaks_calibrated %>% select(year, final_peak_scalar), by = "year") %>%
  arrange(year, week) %>%
  mutate(time_index = year + (week - 10) / 17)

# -----------------------------------------------------------------------------
# 4. Fit Seasonal Shape GAM on Standardized (0–1) Data
# -----------------------------------------------------------------------------
gam_train_shape <- master_grid_eul %>%
  filter(
    year >= 2011, year <= 2024,
    !is.na(eul_prop)
  )

shape_gam <- gam(
  eul_prop ~ s(week, k = 7),
  data = gam_train_shape,
  family = quasibinomial(link = "logit"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# 5. Predict Shape & Rescale (Now Valid for 1998–2024)
# -----------------------------------------------------------------------------
shape_preds <- predict(shape_gam, newdata = master_grid_eul, type = "link", se.fit = TRUE)

master_eul_reconstructed <- master_grid_eul %>%
  mutate(
    prop_fit     = shape_gam$family$linkinv(shape_preds$fit),
    prop_lwr     = shape_gam$family$linkinv(shape_preds$fit - 1.96 * shape_preds$se.fit),
    prop_upr     = shape_gam$family$linkinv(shape_preds$fit + 1.96 * shape_preds$se.fit),
    
    eul_gam_pred = prop_fit * final_peak_scalar,
    eul_lwr_95   = prop_lwr * final_peak_scalar,
    eul_upr_95   = prop_upr * final_peak_scalar,
    
    data_source  = if_else(!is.na(eulachon_input), "Observed", "Predicted"),
    eulachon_final = if_else(!is.na(eulachon_input),
                             as.numeric(eulachon_input),
                             eul_gam_pred)
  )


# -----------------------------------------------------------------------------
# 6. Diagnostic Plots
# -----------------------------------------------------------------------------

# Plot 1: Standardized Shape GAM (Excl. 2020)
p_shape_fit <- ggplot(master_eul_reconstructed %>% filter(year >= 2011), aes(x = week)) +
  geom_ribbon(aes(ymin = prop_lwr, ymax = prop_upr), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = prop_fit, color = "GAM Standardized Curve"), linewidth = 1) +
  geom_point(aes(y = eul_prop, color = "Observed Proportions"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4) +
  scale_color_manual(values = c("Observed Proportions" = "black", "GAM Standardized Curve" = "blue")) +
  labs(
    title = "Standardized Weekly Eulachon Shape GAM (0–1 Scale, Excl. 2020)",
    subtitle = "Fitted across valid training years with equal weighting",
    x = "Week (10–26)",
    y = "Proportion of Annual Peak",
    color = "Legend"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_shape_fit)

# Plot 2: Calibrated Rescaled Absolute Biomass Predictions (2011–2024, Excl. 2020)
p_eul_recent_rescaled <- ggplot(master_eul_reconstructed %>% filter(year >= 2011), aes(x = week)) +
  geom_ribbon(aes(ymin = eul_lwr_95, ymax = eul_upr_95), fill = "seagreen", alpha = 0.2) +
  geom_line(aes(y = eul_gam_pred, color = "Rescaled GAM Prediction"), linewidth = 0.8) +
  geom_point(aes(y = eulachon_input, color = "Observed Data"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Data" = "black", "Rescaled GAM Prediction" = "seagreen4")) +
  labs(
    title = "Peak-Calibrated Reconstructed Weekly Eulachon Biomass (2011–2024)",
    subtitle = "Excludes 2020 COVID gap; Shape rescaled directly to observed/predicted peak height",
    x = "Week (10–26)",
    y = "Eulachon Biomass / Count",
    color = "Legend"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_eul_recent_rescaled)


#back to 1998-------
cat("\n========================================================\n")
cat("   GENERATING FULL EULACHON RECONSTRUCTION (1998–2024)   \n")
cat("========================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Format Final Dataset
# -----------------------------------------------------------------------------
eulachon_full_reconstructed <- master_eul_reconstructed %>%
  select(
    year, 
    week, 
    month, 
    date, 
    eulachon_input, 
    prop_fit, 
    prop_lwr, 
    prop_upr, 
    final_peak_scalar, 
    eul_gam_pred, 
    eul_lwr_95, 
    eul_upr_95, 
    data_source, 
    eulachon_final,
    time_index
  ) %>%
  arrange(year, week)

# -----------------------------------------------------------------------------
# 2. Diagnostic Plots (Full Eulachon Series 1998–2024)
# -----------------------------------------------------------------------------

# Plot A: Full Faceted Series (1998–2024, Weeks 10–26)
p_eul_full_faceted <- ggplot(eulachon_full_reconstructed, aes(x = week)) +
  geom_ribbon(aes(ymin = eul_lwr_95, ymax = eul_upr_95), fill = "seagreen", alpha = 0.2) +
  geom_line(aes(y = eul_gam_pred, color = "GAM Reconstructed Fit"), linewidth = 0.8) +
  geom_point(aes(y = eulachon_input, color = "Observed Data"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Data" = "black", "GAM Reconstructed Fit" = "seagreen4")) +
  labs(
    title = "Reconstructed Weekly Eulachon Biomass (1998–2024, Weeks 10–26)",
    subtitle = "Faceted by year (excl. 2020); Black points = Observed, Green line = Peak-calibrated GAM",
    x = "Week (10–26)",
    y = "Eulachon Biomass / Count",
    color = "Legend"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_eul_full_faceted)

# Plot B: Continuous Timeline across all years (1998–2024)
p_eul_full_continuous <- ggplot(eulachon_full_reconstructed, aes(x = time_index)) +
  geom_ribbon(aes(ymin = eul_lwr_95, ymax = eul_upr_95), fill = "seagreen", alpha = 0.15) +
  geom_line(aes(y = eulachon_final), color = "gray40", linewidth = 0.5) +
  geom_point(aes(y = eulachon_final, color = data_source), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "firebrick")) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Continuous Reconstructed Eulachon Time Series (1998–2024, Weeks 10–26)",
    subtitle = "Black = Observed weekly data | Red = Peak-calibrated GAM predictions",
    x = "Year",
    y = "Eulachon Biomass / Count",
    color = "Data Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_eul_full_continuous)

# -----------------------------------------------------------------------------
# 3. Save Output Files
# -----------------------------------------------------------------------------

write.csv(eulachon_full_reconstructed, file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "eulachon_weekly_reconstructed_faceted_1998_2024.png"), p_eul_full_faceted, width = 12, height = 8)
ggsave(file.path(output_dir, "eulachon_weekly_reconstructed_continuous_1998_2024.png"), p_eul_full_continuous, width = 12, height = 6)
