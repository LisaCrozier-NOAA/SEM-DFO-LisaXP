# ==============================================================================
# CSL & Eulachon Weekly Reconstruction Pipeline (1998–2024, Weeks 10–26)
# ==============================================================================
rm(list=ls())

library(tidyverse)
library(mgcv)
library(zoo)

output_dir <- "copilot/outputs_csl_cr"

# -----------------------------------------------------------------------------
# 1. Load & Clean Raw Data
# -----------------------------------------------------------------------------

# CSL Weekly Census Data
week_all_census <- read.csv(file.path(output_dir, "csl_week_all_census.csv"), row.names = NULL) %>%
  filter(year >= 1998, week >= 10, week <= 26)

# CSL Annual Baseline Data (1993–2024)
csl_annual <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"), row.names = NULL) %>%
  filter(year >= 1998) %>%
  select(year, csl_annual_mean)

# Eulachon Weekly Data (2011–2024)
eulachon_week_all <- read.csv(file.path(output_dir, "jake_week_all.csv"), row.names = NULL) %>%
  filter(year >= 1998, week >= 10, week <= 26) %>%
  select(year, week, eulachon_ssb_4week_est)

# Eulachon Annual Master Index
eulachon_annual <- read.csv(file.path(output_dir, "eulachon_reconstructed_index_lbs_1993_2024.csv")) %>%
  filter(year >= 1998) %>%
  select(year, eulachon_lbs_reconstructed)

# -----------------------------------------------------------------------------
# 2. Build Seasonal Mean Shapes & Scaling Factors (2011–2024 Baseline)
# -----------------------------------------------------------------------------

# A. Average 2011–2024 Weekly Eulachon Seasonal Profile
eulachon_seasonal_shape <- eulachon_week_all %>%
  filter(year >= 2011, year <= 2024) %>%
  group_by(week) %>%
  summarise(eulachon_mean_shape = mean(eulachon_ssb_4week_est, na.rm = TRUE), .groups = "drop")

# B. Baseline Mean Annual Eulachon LBS (2011–2024)
eulachon_lbs_baseline_mean <- eulachon_annual %>%
  filter(year >= 2011, year <= 2024) %>%
  summarise(mean_lbs = mean(eulachon_lbs_reconstructed, na.rm = TRUE)) %>%
  pull(mean_lbs)

# C. Calculate Annual CSL Mean & Log Annual Scale Factor
csl_annual_prep <- csl_annual %>%
  mutate(
    # Fallback to mean if csl_annual_mean is missing or zero
    csl_annual_mean_clean = if_else(is.na(csl_annual_mean) | csl_annual_mean <= 0, 
                                    mean(csl_annual_mean, na.rm = TRUE), 
                                    csl_annual_mean),
    log_annual_scale = log(csl_annual_mean_clean)
  )

# -----------------------------------------------------------------------------
# 3. Construct Complete Master Grid (1998–2024, Weeks 10–26)
# -----------------------------------------------------------------------------

master_grid <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(week_all_census, by = c("year", "week")) %>%
  left_join(csl_annual_prep, by = "year") %>%
  left_join(eulachon_week_all, by = c("year", "week")) %>%
  left_join(eulachon_annual, by = "year") %>%
  left_join(eulachon_seasonal_shape, by = "week") %>%
  arrange(year, week) %>%
  mutate(
    # Reconstruct Weekly Eulachon:
    # Use observed eulachon_ssb_4week_est if present; otherwise scale 2011-2024 mean shape
    eulachon_reconstructed = if_else(
      !is.na(eulachon_ssb_4week_est),
      eulachon_ssb_4week_est,
      eulachon_mean_shape * (eulachon_lbs_reconstructed / eulachon_lbs_baseline_mean)
    ),
    # Continuous time index for multi-year plotting
    time_index = year + (week - 10) / 17
  )

# -----------------------------------------------------------------------------
# 4. Fit GAM Model on High-Quality Observed Data (2011–2024)
# -----------------------------------------------------------------------------

# Filter training set: 2011-2024 with observed CSL and valid Eulachon
gam_train_data <- master_grid %>%
  filter(
    year >= 2011, year <= 2024,
    !is.na(csl_nonpup_total_emb_lisa),
    !is.na(eulachon_reconstructed)
  )

csl_gam <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, k = 5) + 
    s(eulachon_reconstructed, k = 5) + 
    offset(log_annual_scale),
  data = gam_train_data,
  family = quasipoisson(link = "log"),
  method = "REML"
)

summary(csl_gam)

# -----------------------------------------------------------------------------
# 5. Predict & Construct Final Time Series (1998–2024)
# -----------------------------------------------------------------------------

# Predict across full dataset
gam_preds <- predict(csl_gam, newdata = master_grid, type = "link", se.fit = TRUE)

master_reconstructed <- master_grid %>%
  mutate(
    fit_link     = gam_preds$fit,
    se_link      = gam_preds$se.fit,
    
    # Scale predictions back to count scale
    csl_gam_pred = exp(fit_link),
    csl_lwr_95   = exp(fit_link - 1.96 * se_link),
    csl_upr_95   = exp(fit_link + 1.96 * se_link),
    
    # Hierarchical selection: Observed > Predicted
    data_source  = if_else(!is.na(csl_nonpup_total_emb_lisa), "Observed", "Predicted"),
    csl_final    = if_else(!is.na(csl_nonpup_total_emb_lisa), 
                           as.numeric(csl_nonpup_total_emb_lisa), 
                           csl_gam_pred)
  )

# -----------------------------------------------------------------------------
# 6. Visualization & Diagnostics
# -----------------------------------------------------------------------------

# Plot A: Recent Years (2011–2024) - Observed Points vs GAM Curve
p_recent <- ggplot(master_reconstructed %>% filter(year >= 2011), aes(x = week)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = csl_gam_pred, color = "GAM Model Fit"), linewidth = 0.8) +
  geom_point(aes(y = csl_nonpup_total_emb_lisa, color = "Observed Data"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Data" = "black", "GAM Model Fit" = "steelblue")) +
  labs(
    title = "Weekly CSL Phenology Model Fit (2011–2024)",
    subtitle = "Evaluates GAM predictions across observed and missing weeks",
    x = "Week (10–26)",
    y = "CSL Non-Pup Count",
    color = "Data Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_recent)

# Plot B: Continuous Full Reconstructed Series (1998–2024)
p_full_series <- ggplot(master_reconstructed, aes(x = time_index)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.15) +
  geom_line(aes(y = csl_final), color = "gray40", linewidth = 0.5) +
  geom_point(aes(y = csl_final, color = data_source), size = 1.3) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "firebrick")) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Reconstructed CSL Weekly Time Series (1998–2024, Weeks 10–26)",
    subtitle = "Black points indicate observed surveys; red points indicate GAM reconstructed values",
    x = "Year",
    y = "CSL Count",
    color = "Series Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_full_series)

# ==============================================================================
# Dual-Axis Plot: Eulachon (Observed vs Predicted) & CSL Overlay (1998–2024)---------
# ==============================================================================

# 1. Calculate dynamic scaling factor to align Eulachon and CSL on dual axes
max_csl <- max(master_reconstructed$csl_final, na.rm = TRUE)
max_eulachon <- max(master_reconstructed$eulachon_reconstructed, na.rm = TRUE)
eulachon_scale_factor <- max_csl / max_eulachon

# Add indicator for observed vs predicted eulachon
plot_data <- master_reconstructed %>%
  mutate(
    eulachon_source = if_else(!is.na(eulachon_ssb_4week_est), "Eulachon (Observed)", "Eulachon (Reconstructed)")
  )

# 2. Build Multi-Panel Faceted Plot (by Year) for Detailed Inspection
p_overlay_faceted <- ggplot(plot_data %>% filter(year >= 2011), aes(x = week)) +
  # Eulachon Layer (Area / Line)
  geom_area(aes(y = eulachon_reconstructed * eulachon_scale_factor, fill = eulachon_source), 
            alpha = 0.35) +
  geom_line(aes(y = eulachon_reconstructed * eulachon_scale_factor, color = eulachon_source), 
            linewidth = 0.8, linetype = "dashed") +
  
  # CSL GAM Confidence Interval Ribbon
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "darkblue", alpha = 0.15) +
  
  # CSL Layer (Line + Points)
  geom_line(aes(y = csl_gam_pred, color = "CSL (GAM Predicted)"), linewidth = 0.9) +
  geom_point(aes(y = csl_final, color = data_source), size = 2) +
  
  # Axes & Faceting
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_y_continuous(
    name = "CSL Count",
    sec.axis = sec_axis(~ . / eulachon_scale_factor, name = "Eulachon Biomass Index")
  ) +
  scale_color_manual(
    name = "Series",
    values = c(
      "CSL (GAM Predicted)"    = "darkblue",
      "Observed"              = "black",
      "Predicted"             = "firebrick",
      "Eulachon (Observed)"     = "seagreen4",
      "Eulachon (Reconstructed)" = "darkgoldenrod2"
    )
  ) +
  scale_fill_manual(
    name = "Series",
    values = c(
      "Eulachon (Observed)"     = "seagreen3",
      "Eulachon (Reconstructed)" = "goldenrod1"
    )
  ) +
  labs(
    title = "Eulachon Run Dynamics vs. CSL Phenology (2011–2024)",
    subtitle = "Dashed area = Eulachon biomass (Primary scale) | Solid lines/points = CSL counts (Secondary scale)",
    x = "Week (10–26)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.title.y.right = element_text(color = "seagreen4")
  )

print(p_overlay_faceted)

# 3. Continuous Multi-Year Overlay (1998–2024)
p_overlay_continuous <- ggplot(plot_data, aes(x = time_index)) +
  # Eulachon Area Background
  geom_area(aes(y = eulachon_reconstructed * eulachon_scale_factor, fill = eulachon_source), 
            alpha = 0.3) +
  # CSL Reconstructed Line
  geom_line(aes(y = csl_final, color = "CSL Final (Obs + Pred)"), linewidth = 0.7) +
  geom_point(aes(y = csl_final, color = data_source), size = 1.1) +
  
  scale_y_continuous(
    name = "CSL Count",
    sec.axis = sec_axis(~ . / eulachon_scale_factor, name = "Eulachon Biomass Index")
  ) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  scale_color_manual(
    name = "CSL Data",
    values = c(
      "CSL Final (Obs + Pred)" = "navy",
      "Observed"               = "black",
      "Predicted"              = "firebrick"
    )
  ) +
  scale_fill_manual(
    name = "Eulachon Data",
    values = c(
      "Eulachon (Observed)"     = "seagreen3",
      "Eulachon (Reconstructed)" = "goldenrod1"
    )
  ) +
  labs(
    title = "Full Reconstructed Time Series: Eulachon & CSL Overlay (1998–2024)",
    subtitle = "Shaded background indicates weekly eulachon biomass; points/lines represent weekly CSL counts",
    x = "Year"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_overlay_continuous)

# Save Overlay Plots
ggsave(file.path(output_dir, "eulachon_csl_overlay_faceted.png"), p_overlay_faceted, width = 12, height = 9)
ggsave(file.path(output_dir, "eulachon_csl_overlay_continuous.png"), p_overlay_continuous, width = 13, height = 6)

# -----------------------------------------------------------------------------
# 7. Save Outputs
# -----------------------------------------------------------------------------

write.csv(master_reconstructed, file.path(output_dir, "csl_reconstructed_weekly_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "csl_weekly_fit_2011_2024.png"), p_recent, width = 11, height = 8)
ggsave(file.path(output_dir, "csl_weekly_reconstructed_1998_2024.png"), p_full_series, width = 12, height = 6)