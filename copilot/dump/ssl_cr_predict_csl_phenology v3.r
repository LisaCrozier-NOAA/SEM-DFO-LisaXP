#52-week version-- to be dumped


output_dir <- "outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

library(tidyverse)
library(mgcv)
library(readxl)
library(lubridate)

# -----------------------------------------------------------------------------
# 1. Load Long-Term Annual Sea Lion Dataset (1976–2024)
# -----------------------------------------------------------------------------

jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
lre_dat_yearly %>% select(year,csl_nonpup_total_emb)

library(lubridate)
path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/LisaDataProcessScripts 2025/mammals/ODFW atlas count Columbia River v 20250218.xlsx"

odfw <- read_xlsx(path = path, sheet = 1)  %>%
  filter(spp == "ZC", 
         #    type == "VISUAL", 
         location == "COLUMBIA RIVER-EAST MOORING BASIN") %>%
  mutate(
    date_parsed = ymd(datemil),
    year = year(date_parsed),
    week = week(date_parsed)   # Use isoweek(date_parsed) if aligned to ISO weeks
  )

csl_year <- odfw  %>% 
  filter(week %in% 10:26) %>% #jake used 14:25 and just visual surveys
  group_by(year) %>% 
  summarize(
    csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb_lisa = n(),  # Number of observations
    .groups = "drop"
  ) 

csl_week <- odfw  %>% 
  filter(week %in% 10:26) %>% #jake used 14:25 and just visual surveys
  group_by(year,week) %>% 
  summarize(
    csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb_lisa = n(),  # Number of observations
    .groups = "drop"
  ) 

lre_dat_yearly_lisa<-left_join(lre_dat_yearly,csl_year,join_by(year))

# -----------------------------------------------------------------------------
# 2. Extract & Format 50-Year Annual CSL Baseline (1976–2024)
# -----------------------------------------------------------------------------
csl_annual_baseline_50yr <- lre_dat_yearly_lisa %>%
  select(year, csl_nonpup_total_emb_lisa, eulachon_ssb_est) %>%
  filter(year >= 1976, year <= 2024) %>%
  rename(csl_annual_mean = csl_nonpup_total_emb_lisa) %>%
  # Ensure all years from 1976 to 2024 exist in the sequence
  complete(year = 1976:2024) %>%
  arrange(year) %>%
  mutate(
    # Interpolate missing years (e.g., 2019/2020 COVID survey gaps) across time series
    csl_annual_mean = zoo::na.approx(csl_annual_mean, na.rm = FALSE, rule = 2),
    # Log scalar for downstream GAM / regression models
    log_annual_scale = log(pmax(1, csl_annual_mean))
  )

# Inspect long-term sea lion baseline summary
cat("--- Long-Term CSL Annual Baseline (1976–2024) ---\n")
print(csl_annual_baseline_50yr, n = 50)

# Save to output directory
write.csv(
  csl_annual_baseline_50yr,
  file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 1. Expand Data Grid to Full Calendar Year (Weeks 1 to 52)
# -----------------------------------------------------------------------------

week_all<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),row.names=NULL)
week_all_lisa<-left_join(week_all,csl_week,join_by(year,week))

# Expand dataset to cover all 52 weeks per year
week_full_year <- week_all_lisa %>%
  complete(year = 2011:2024, week = 1:52) %>%
  left_join(csl_annual_baseline_50yr, by = "year") %>%
  arrange(year, week) %>%
  mutate(
    # Full annual circular domain in radians (2*pi across 52 weeks)
    week_rad         = 2 * pi * (week - 1) / 52,
    year_factor      = factor(year),
    eulachon_input   = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
    log_annual_scale = log(pmax(1, csl_annual_mean))
  )

# -----------------------------------------------------------------------------
# 2. Fit Full Annual Cyclic Models (Weeks 1–52)
# -----------------------------------------------------------------------------
week_final_scaled<-read.csv("csl_week_filled_in_2011_2024.csv",row.names = NULL)

# A. Annual Sine-Cosine Harmonic Model (Parametric)
eulachon_annual_sine <- glm(
  eulachon_input ~ sin(week_rad) + cos(week_rad) + 
    sin(2 * week_rad) + cos(2 * week_rad) + # 2nd harmonic to allow skewness
    log_annual_scale,
  data   = week_full_year %>% filter(!is.na(eulachon_input)),
  family = quasipoisson(link = "log")
)

# B. Annual Cyclic GAM (bs = "cc" forces smooth tail to wrap from week 52 -> week 1)
eulachon_annual_gam <- gam(
  eulachon_input ~ s(week, bs = "cc", k = 12) + s(year_factor, bs = "re"),
  data    = week_full_year %>% filter(!is.na(eulachon_input)),
  family  = tw(link = "log"),
  knots   = list(week = c(1, 52)),
  method  = "REML"
)

# -----------------------------------------------------------------------------
# 3. Predict & Evaluate Across All 52 Weeks
# -----------------------------------------------------------------------------
annual_eval_grid <- tibble(
  week        = 1:52,
  year_factor = factor("2011", levels = levels(week_full_year$year_factor)),
  log_annual_scale = mean(week_full_year$log_annual_scale, na.rm = TRUE)
) %>%
  mutate(
    week_rad = 2 * pi * (week - 1) / 52
  )

# Extract predictions across the full year
annual_eval_grid$sine_pred <- predict(eulachon_annual_sine, newdata = annual_eval_grid, type = "response")
annual_eval_grid$gam_pred  <- exp(predict(eulachon_annual_gam, newdata = annual_eval_grid, type = "link", exclude = "s(year_factor)"))

# -----------------------------------------------------------------------------
# 4. Plot Comparison: 52-Week Sine/Harmonic vs. Cyclic GAM
# -----------------------------------------------------------------------------
p_annual_comp <- ggplot() +
  # Raw data across all weeks
  geom_point(
    data = week_full_year,
    aes(x = week, y = eulachon_input),
    alpha = 0.2, color = "darkgray"
  ) +
  # Parametric Harmonic Sine Model (Red Dashed)
  geom_line(
    data = annual_eval_grid,
    aes(x = week, y = sine_pred, color = "Harmonic Sine Model (2-Harmonics)"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  # Cyclic GAM Curve (Blue Solid)
  geom_line(
    data = annual_eval_grid,
    aes(x = week, y = gam_pred, color = "Full Year Cyclic GAM"),
    linewidth = 1.2, linetype = "solid"
  ) +
  scale_x_continuous(breaks = seq(4, 52, by = 4), name = "Calendar Week (1 to 52)") +
  scale_color_manual(
    values = c("Harmonic Sine Model (2-Harmonics)" = "firebrick", "Full Year Cyclic GAM" = "steelblue")
  ) +
  labs(
    title = "Full Calendar Year Eulachon Phenology (Weeks 1–52)",
    subtitle = "Comparing annual cyclic curves for disaggregating annual landing totals into weekly biomass",
    y = "Eulachon Biomass / Run Index",
    color = "Model Fit"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p_annual_comp)

# Save full 52-week plot to outputs directory
ggsave(
  filename = file.path(output_dir, "eulachon_52week_sine_vs_gam.png"),
  plot     = p_annual_comp,
  width    = 10, height = 5
)