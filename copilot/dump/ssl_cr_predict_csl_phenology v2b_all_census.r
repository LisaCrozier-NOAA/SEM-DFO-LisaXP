
#this version uses visual, uav and airplane surveys for csl instead of just visual
output_dir <- "copilot/outputs_csl_cr"




library(tidyverse)
library(mgcv)
library(zoo)

# -----------------------------------------------------------------------------
# 1. Read Data, Expand 2019 and 2020 Grid via complete(), & Interpolate--------
# -----------------------------------------------------------------------------
      week_all_census <- read.csv(file.path(output_dir,"csl_week_all_census.csv"),row.names=NULL)
      # > names(week_all_census)
      # [1] "year"                      "week"                     
      # [3] "csl_nonpup_total_emb_lisa" "csl_nobs_emb_lisa"        
      
      csl_annual<- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),row.names = NULL) %>%
              filter(year>=1993) %>%
              select(year,csl_annual_mean,log_annual_scale)
      # > names(csl_annual)
      # [1] "year"             "csl_annual_mean"  "eulachon_ssb_est" "log_annual_scale"
      
      eulachon_week_all <- read.csv(file.path(output_dir,"jake_week_all.csv"),row.names=NULL) %>%
        select(year,week,month,date,eulachon_ssb_4week_est ,Spring_Achin_bonn_pass,total_bonn  )
      
      eulachon_annual<- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv")) %>%
        select(year,eulachon_lbs_reconstructed)
    
     
      head(week_all_census)
      
# A. Compute annual CSL mean scalar (incorporating missing 2019 & 2020)--------
csl_annual_from_week <- week_all_census %>%
  filter(week %in% 10:26, year >= 2011, year <= 2024) %>%
#  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  group_by(year) %>%
  summarise(
    csl_annual_mean = if_else(all(is.na(csl_nonpup_total_emb_lisa )), NA_real_, mean(csl_nonpup_total_emb_lisa , na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # Force 2019 and 2020 into the annual summary
  complete(year = 2011:2024) %>%
  mutate(
    csl_annual_mean = zoo::na.approx(csl_annual_mean, na.rm = FALSE, rule = 2)
  )

# test<-left_join(csl_annual_from_week,csl_annual, join_by(year) )
# plot(csl_annual$year,csl_annual$csl_annual_mean,type='l')
# lines(csl_annual_from_week$year,csl_annual_from_week$csl_annual_mean,col=2)

# B. Assemble full weekly dataset and expand grid for 2020---------
week_scaled_prep <- week_all_census %>%
  #filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  # FORCE MISSING YEARS/WEEKS (2020) INTO THE DATASET GRID
  complete(year = 2011:2024, week = full_seq(week, 1)) %>%
  left_join(csl_annual_from_week, by = "year") %>%
  left_join(eulachon_week_all, by = c("year", "week")) %>%
  arrange(year, week) %>%
  mutate(
    # UNGROUPED INTERPOLATION: Fills Eulachon continuously across newly created 2020 rows
    eulachon_input       = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
    year_factor          = factor(year),
    log_annual_scale     = log(pmax(1, csl_annual_mean)),
#    chinook_estuary_pass = if_else(is.na(chinook_estuary_pass), 0, chinook_estuary_pass)
  )

# VERIFY THAT 2020 ROWS NOW EXIST
cat("--- Checking 2020 Rows in Console ---\n")
week_scaled_prep %>%
  filter(year == 2020,week %in% 10:26) %>%
  select(year, week, csl_annual_mean, log_annual_scale, eulachon_input) %>%
  print(n = 13)


# -----------------------------------------------------------------------------
# 3. Fit GAM on Observed Data-------
# -----------------------------------------------------------------------------
csl_scaled_gam <- gam(
  csl_nonpup_total_emb_lisa  ~ s(week, k = 5) + 
    s(eulachon_input, k = 5) + 
    log_annual_scale,
  data = week_scaled_prep %>% 
    filter(!is.na(csl_nonpup_total_emb_lisa ) & !is.na(eulachon_input)),
  family = quasipoisson(link = "log"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# 4. Predict Reconstructed CSL Phenology (2019 & 2020 Fully Reconstructed)--------
# -----------------------------------------------------------------------------
scaled_preds_link <- predict(
  csl_scaled_gam,
  newdata = week_scaled_prep,
  type = "link",
  se.fit = TRUE
)

week_final_scaled <- week_scaled_prep %>%
  mutate(
    fit_link     = scaled_preds_link$fit,
    se_link      = scaled_preds_link$se.fit,
    
    csl_gam_pred = exp(fit_link),
    csl_lwr_95   = exp(fit_link - 1.96 * se_link),
    csl_upr_95   = exp(fit_link + 1.96 * se_link),
    
    csl_final    = if_else(is.na(csl_nonpup_total_emb_lisa ), csl_gam_pred, as.numeric(csl_nonpup_total_emb_lisa )),
    csl_final    = pmax(0, csl_final),
    
    # Create a continuous time variable for multi-year time series plotting
    # Assuming standard calendar week spacing
    time_index   = year + (week - 1) / 52
  )


# -----------------------------------------------------------------------------
# 6. Plot Gap Inspection (2018–2021)---------
# -----------------------------------------------------------------------------
csl_weekly_fill_2011_2024_plot <-
ggplot(week_final_scaled %>% filter(year>=2011), aes(x = week)) +
#  ggplot(week_final_scaled %>% filter(year %in% c(2018, 2019, 2020, 2021)), aes(x = week)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = csl_gam_pred, color = "Reconstructed Curve"), linewidth = 1) +
  geom_point(aes(y = csl_nonpup_total_emb_lisa , color = "Observed Points"), size = 2.5, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Points" = "black", "Reconstructed Curve" = "steelblue")) +
  labs(
    title = "Reconstructed Sea Lion Phenology Across COVID Gap (2018–2021)",
    subtitle = "2020 explicitly expanded via complete() and predicted from GAM",
    x = "Week Number (Apr–Jun)",
    y = "CSL Non-Pup Count",
    color = "Legend"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(csl_weekly_fill_2011_2024_plot)
# -----------------------------------------------------------------------------
# 3. Plot 2: Continuous Time Series Plot (Overall Model Performance Across Years)
# -----------------------------------------------------------------------------
p_continuous <- ggplot(week_final_scaled %>% filter(year>=2011), aes(x = time_index)) +
  # 95% Confidence Interval Ribbon
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.25) +
  # GAM Prediction Line
  geom_line(aes(y = csl_gam_pred, color = "GAM Predicted"), size = 0.9) +
  # Observed Data Points
  geom_point(aes(y = csl_nonpup_total_emb_lisa , color = "Observed Data"), size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    values = c("Observed Data" = "black", "GAM Predicted" = "steelblue")
  ) +
  scale_x_continuous(breaks = 2011:2024) +
  labs(
    title = "Continuous Time Series of CSL Phenology Model Fits (2011–2024)",
    subtitle = "Demonstrates multi-year fit and 2019 back-prediction based on Eulachon biomass",
    x = "Year",
    y = "CSL Non-Pup Count",
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
write.csv(week_final_scaled,"csl_week_filled_in_2011_2024.csv",row.names=FALSE)

ggsave(file="csl_weekly_fill_2011_2024_plot.png",csl_weekly_fill_2011_2024_plot)
ggsave(file="csl_weekly_continuous_2011_2024_plot.png",p_continuous)


