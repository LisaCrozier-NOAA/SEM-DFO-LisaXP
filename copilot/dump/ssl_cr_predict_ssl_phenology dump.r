
#True Hazard Overlap: 
#The actual predation overlap index is the weekly product of Predicted Sea Lions * Estuary Chinook:
#Overlap(Salmon\ssl) = cumulative sum over all weeks in Apr–Jun} {ssl_predicted * Chinook_Estuary}


library(tidyverse)
library(mgcv)
library(zoo)

# -----------------------------------------------------------------------------
# 1. Read Data & Back-Track Chinook (2 Weeks)
# -----------------------------------------------------------------------------
week_all <- read.csv(
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
  row.names = NULL
)

# Back-track Chinook passage by 2 weeks (Estuary arrival timing)
chinook_estuary_df <- week_all %>%
  select(year, week, Spring_Achin_bonn_pass) %>%
  mutate(week_estuary = week - 2) %>%
  filter(week_estuary >= 1) %>%
  select(year, week = week_estuary, chinook_estuary_pass = Spring_Achin_bonn_pass)

# -----------------------------------------------------------------------------
# 2. Build Dataset, Expand 2019 and 2020 Grid via complete(), & Interpolate
# -----------------------------------------------------------------------------
# A. Compute annual ssl mean scalar (incorporating missing 2019 & 2020)
ssl_annual_mean <- week_all %>%
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  group_by(year) %>%
  summarise(
    ssl_annual_mean = if_else(all(is.na(ssl_nonpup_total_sj)), NA_real_, mean(ssl_nonpup_total_sj, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # Force 2019 and 2020 into the annual summary
  complete(year = 2011:2024) %>%
  mutate(
    ssl_annual_mean = zoo::na.approx(ssl_annual_mean, na.rm = FALSE, rule = 2)
  )

# B. Assemble full weekly dataset and expand grid for 2020
week_scaled_prep <- week_all %>%
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  # FORCE MISSING YEARS/WEEKS (2020) INTO THE DATASET GRID
  complete(year = 2011:2024, week = full_seq(week, 1)) %>%
  left_join(ssl_annual_mean, by = "year") %>%
  left_join(chinook_estuary_df, by = c("year", "week")) %>%
  arrange(year, week) %>%
  mutate(
    # UNGROUPED INTERPOLATION: Fills Eulachon continuously across newly created 2020 rows
    eulachon_input       = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
    year_factor          = factor(year),
    log_annual_scale     = log(pmax(1, ssl_annual_mean)),
    chinook_estuary_pass = if_else(is.na(chinook_estuary_pass), 0, chinook_estuary_pass)
  )

# VERIFY THAT 2020 ROWS NOW EXIST
cat("--- Checking 2020 Rows in Console ---\n")
week_scaled_prep %>%
  filter(year == 2020) %>%
  select(year, week, ssl_annual_mean, log_annual_scale, eulachon_input) %>%
  print(n = 13)

# -----------------------------------------------------------------------------
# 3. Fit GAM on Observed Data
# -----------------------------------------------------------------------------
ssl_scaled_gam <- gam(
  ssl_nonpup_total_sj ~ s(week, k = 5) + 
    s(eulachon_input, k = 5) + 
    log_annual_scale,
  data = week_scaled_prep %>% filter(!is.na(ssl_nonpup_total_sj) & !is.na(eulachon_input)),
  family = quasipoisson(link = "log"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# 4. Predict Reconstructed ssl Phenology (2019 & 2020 Fully Reconstructed)--------
# -----------------------------------------------------------------------------
scaled_preds_link <- predict(
  ssl_scaled_gam,
  newdata = week_scaled_prep,
  type = "link",
  se.fit = TRUE
)

week_final_scaled <- week_scaled_prep %>%
  mutate(
    fit_link     = scaled_preds_link$fit,
    se_link      = scaled_preds_link$se.fit,
    
    ssl_gam_pred = exp(fit_link),
    ssl_lwr_95   = exp(fit_link - 1.96 * se_link),
    ssl_upr_95   = exp(fit_link + 1.96 * se_link),
    
    ssl_final    = if_else(is.na(ssl_nonpup_total_sj), ssl_gam_pred, as.numeric(ssl_nonpup_total_sj)),
    ssl_final    = pmax(0, ssl_final),
    
    # Create a continuous time variable for multi-year time series plotting
    # Assuming standard calendar week spacing
    time_index   = year + (week - 1) / 52,
    
    weekly_ssl_chinook_overlap = ssl_final * chinook_estuary_pass
  )

# -----------------------------------------------------------------------------
# 5. Extract Annual Integrated Predation Overlap Index--------
# -----------------------------------------------------------------------------
safe_sum <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

salmon_ssl_annual_index <- week_final_scaled %>%
  group_by(year) %>%
  summarise(
    weeks_evaluated        = n(),
    ssl_obs_weeks          = sum(!is.na(ssl_nonpup_total_sj)),
    total_estuary_chinook  = safe_sum(chinook_estuary_pass),
    total_ssl_exposure     = safe_sum(ssl_final),
    
    i_ssl_chinook_overlap  = safe_sum(weekly_ssl_chinook_overlap)
  )

cat("\n--- ANNUAL SALMON x ssl PREDATION OVERLAP INDEX (2011-2024) ---\n")
print(salmon_ssl_annual_index, n = 20)

# -----------------------------------------------------------------------------
# 6. Plot Gap Inspection (2018–2021)---------
# -----------------------------------------------------------------------------
ssl_weekly_fill_2011_2024_plot <-
ggplot(week_final_scaled , aes(x = week)) +
#  ggplot(week_final_scaled %>% filter(year %in% c(2018, 2019, 2020, 2021)), aes(x = week)) +
  geom_ribbon(aes(ymin = ssl_lwr_95, ymax = ssl_upr_95), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = ssl_gam_pred, color = "Reconstructed Curve"), linewidth = 1) +
  geom_point(aes(y = ssl_nonpup_total_sj, color = "Observed Points"), size = 2.5, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Points" = "black", "Reconstructed Curve" = "steelblue")) +
  labs(
    title = "Reconstructed Sea Lion Phenology Across COVID Gap (2018–2021)",
    subtitle = "2020 explicitly expanded via complete() and predicted from GAM",
    x = "Week Number (Apr–Jun)",
    y = "ssl Non-Pup Count",
    color = "Legend"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# -----------------------------------------------------------------------------
# 3. Plot 2: Continuous Time Series Plot (Overall Model Performance Across Years)
# -----------------------------------------------------------------------------
p_continuous <- ggplot(week_final_scaled, aes(x = time_index)) +
  # 95% Confidence Interval Ribbon
  geom_ribbon(aes(ymin = ssl_lwr_95, ymax = ssl_upr_95), fill = "steelblue", alpha = 0.25) +
  # GAM Prediction Line
  geom_line(aes(y = ssl_gam_pred, color = "GAM Predicted"), size = 0.9) +
  # Observed Data Points
  geom_point(aes(y = ssl_nonpup_total_sj, color = "Observed Data"), size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    values = c("Observed Data" = "black", "GAM Predicted" = "steelblue")
  ) +
  scale_x_continuous(breaks = 2011:2024) +
  labs(
    title = "Continuous Time Series of ssl Phenology Model Fits (2011–2024)",
    subtitle = "Demonstrates multi-year fit and 2019 back-prediction based on Eulachon biomass",
    x = "Year",
    y = "ssl Non-Pup Count",
    color = "Series"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_continuous)


#--------------------
#Save weekly predictions and plot--------
#-------------------------

write.csv(week_final_scaled,"ssl_week_filled_in_2011_2024.csv")

ggsave(file="ssl_weekly_fill_2011_2024_plot.png",ssl_weekly_fill_2011_2024_plot)
ggsave(file="ssl_weekly_continuous_2011_2024_plot.png",ssl_weekly_fill_2011_2024_plot)


