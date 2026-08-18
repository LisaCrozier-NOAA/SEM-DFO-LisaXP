


#RESULTS
# baseline seasonal signal to new outputs directory
output_dir <- "outputs_csl_cr"
eulachon_seasonal_signal<-read.csv(file.path(output_dir, "eulachon_baseline_seasonal_signal.csv"),row.names = NULL)
head(eulachon_seasonal_signal)




library(tidyverse)
library(mgcv)
library(zoo)

#note, eulachon_input was created in ssl_cr_predict_csl_phenology v2.r to fill in missing data:
# UNGROUPED INTERPOLATION: Fills Eulachon continuously across newly created 2020 rows
#eulachon_input       = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
week_final_scaled<-read.csv("csl_week_filled_in_2011_2024.csv",row.names = NULL)

week_raw<-week_final_scaled %>% select(year,week,month,date,Spring_Achin_bonn_pass,csl_nonpup_total_emb,eulachon_ssb_4week_est)
head(week_raw)
# -----------------------------------------------------------------------------
output_dir <- "outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}


# -----------------------------------------------------------------------------
# 1. Fit Tweedie GAM to Weekly Eulachon Biomass
# -----------------------------------------------------------------------------
# family = tw() natively handles exact zero-biomass weeks without errors
eulachon_gam <- gam(
  eulachon_input ~ s(week, k = 6) + s(year_factor, bs = "re"),
  data = week_final_scaled %>% filter(!is.na(eulachon_input)),
  family = tw(link = "log"),
  method = "REML"
)

summary(eulachon_gam)

# -----------------------------------------------------------------------------
# 2. Build Consistent Prediction Grid
# -----------------------------------------------------------------------------
# We supply year_factor in the grid so predict.gam() is satisfied, 
# while using exclude = "s(year_factor)" to compute the population-wide average.
eul_pred_grid <- tibble(
  week        = seq(min(week_final_scaled$week), max(week_final_scaled$week), length.out = 100),
  year_factor = factor("2011", levels = levels(week_final_scaled$year_factor))
)

# Predict population-level average run timing
gam_eul_preds <- predict(
  eulachon_gam,
  newdata = eul_pred_grid,
  type    = "link",
  se.fit  = TRUE,
  exclude = "s(year_factor)"
)

# Format population-wide predictions for plotting
eul_pop_phenology <- eul_pred_grid %>%
  mutate(
    fit_link = gam_eul_preds$fit,
    se_link  = gam_eul_preds$se.fit,
    
    # Back-transform from log scale
    fit_eul  = exp(fit_link),
    lwr_95   = pmax(0, exp(fit_link - 1.96 * se_link)),
    upr_95   = exp(fit_link + 1.96 * se_link)
  )

# -----------------------------------------------------------------------------
# 3. Plot 1: Absolute Biomass Consistency (Annuals vs. Population Average)
# -----------------------------------------------------------------------------
p_consistency <- ggplot() +
  # Individual year weekly trajectories
  geom_line(
    data = week_final_scaled,
    aes(x = week, y = eulachon_input, group = year, color = factor(year)),
    alpha = 0.5, linewidth = 0.8
  ) +
  # Overall fitted average 95% CI ribbon
  geom_ribbon(
    data = eul_pop_phenology,
    aes(x = week, ymin = lwr_95, ymax = upr_95),
    fill = "black", alpha = 0.2
  ) +
  # Overall fitted average curve
  geom_line(
    data = eul_pop_phenology,
    aes(x = week, y = fit_eul),
    color = "black", linewidth = 1.5, linetype = "solid"
  ) +
  labs(
    title = "Eulachon Run Timing & Biomass Consistency (2011–2024)",
    subtitle = "Thin colored lines = observed/interpolated annual runs; Heavy black curve = GAM baseline average (± 95% CI)",
    x = "Week Number (Apr–Jun)",
    y = "Eulachon 4-Week SSB Estimate",
    color = "Year"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p_consistency)

# -----------------------------------------------------------------------------
# 4. Plot 2: Proportional Run Timing (Standardized Peak Alignment)
# -----------------------------------------------------------------------------
# Standardize each year's run from 0 to 1 to isolate timing from run magnitude
week_standardized <- week_final_scaled %>%
  group_by(year) %>%
  mutate(
    eul_max = max(eulachon_input, na.rm = TRUE),
    eul_prop = if_else(eul_max == 0, 0, eulachon_input / eul_max)
  ) %>%
  ungroup()

p_proportional <- ggplot(week_standardized, aes(x = week, y = eul_prop, group = year, color = factor(year))) +
  geom_line(linewidth = 0.9, alpha = 0.7) +
  stat_summary(
    aes(group = 1), 
    fun = mean, geom = "line", color = "black", linewidth = 1.8, linetype = "solid"
  ) +
  labs(
    title = "Standardized Eulachon Run Timing Alignment",
    subtitle = "Proportion of annual peak biomass (Heavy black line = mean proportional timing profile)",
    x = "Week Number (Apr–Jun)",
    y = "Proportion of Peak Biomass",
    color = "Year"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p_proportional)


library(tidyverse)
library(readxl)

# -----------------------------------------------------------------------------
# 0. Create Output Directory
# -----------------------------------------------------------------------------
output_dir <- "outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# -----------------------------------------------------------------------------
# 1. Extract Dominant Seasonal Eulachon Phenology Signal-----
# -----------------------------------------------------------------------------
# Standardize the GAM baseline curve so weekly values sum to 1 across weeks
eulachon_seasonal_signal <- eul_pop_phenology %>%
  select(week, fit_eul) %>%
  mutate(
    # Weekly proportion of annual total run
    eul_prop_signal = fit_eul / sum(fit_eul)
  )

# Save baseline seasonal signal to new outputs directory
write.csv(
  eulachon_seasonal_signal, 
  file.path(output_dir, "eulachon_baseline_seasonal_signal.csv"), 
  row.names = FALSE
)

cat("Extracted and saved baseline seasonal Eulachon signal across weeks.\n")

