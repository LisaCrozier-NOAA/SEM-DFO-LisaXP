#weeks 10-26 fleshed out

# -----------------------------------------------------------------------------
# Setup & Package Loading--------
# -----------------------------------------------------------------------------
output_dir <- "copilot/outs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

library(tidyverse)
library(mgcv)
library(readxl)
library(lubridate)
library(zoo)
library(scales)

# -----------------------------------------------------------------------------
# Part 1: Load Raw Inputs & Calculate Lisa's Spring CSL Metrics (Weeks 10–26)--------
# -----------------------------------------------------------------------------
jake_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path, "lre_dat_yearly.RData"), verbose = TRUE)

odfw_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/LisaDataProcessScripts 2025/mammals/ODFW atlas count Columbia River v 20250218.xlsx"

odfw <- read_xlsx(path = odfw_path, sheet = 1) %>%
  filter(
    spp == "ZC", 
    location == "COLUMBIA RIVER-EAST MOORING BASIN"
  ) %>%
  mutate(
    date_parsed = ymd(datemil),
    year        = year(date_parsed),
    week        = week(date_parsed)
  )

# Calculate Lisa's Spring Annual CSL Mean (Weeks 10–26)
csl_year <- odfw %>% 
  filter(week %in% 10:26) %>%
  group_by(year) %>% 
  summarize(
    csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb_lisa         = n(),
    .groups                   = "drop"
  ) 

# Calculate Lisa's Spring Weekly CSL Mean (Weeks 10–26)
csl_week <- odfw %>% 
  filter(week %in% 10:26) %>%
  group_by(year, week) %>% 
  summarize(
    csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb_lisa         = n(),
    .groups                   = "drop"
  ) 

# Merge with long-term yearly dataset
lre_dat_yearly_lisa <- left_join(lre_dat_yearly, csl_year, by = "year")

# -----------------------------------------------------------------------------
# Part 2: Extract & Format 50-Year Annual CSL Baseline (1976–2024)--------
# -----------------------------------------------------------------------------
csl_annual_baseline_50yr <- lre_dat_yearly_lisa %>%
  select(year, csl_nonpup_total_emb_lisa, eulachon_ssb_est) %>%
  filter(year >= 1976, year <= 2024) %>%
  rename(csl_annual_mean = csl_nonpup_total_emb_lisa) %>%
  complete(year = 1976:2024) %>%
  arrange(year) %>%
  mutate(
    # Interpolate missing survey years (e.g., COVID gaps)
    csl_annual_mean  = zoo::na.approx(csl_annual_mean, na.rm = FALSE, rule = 2),
    log_annual_scale = log(pmax(1, csl_annual_mean))
  )

cat("--- Long-Term CSL Annual Baseline (1976–2024: Weeks 10–26) ---\n")
print(csl_annual_baseline_50yr, n = 50)

write.csv(
  csl_annual_baseline_50yr,
  file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 3: Expand Spring Grid (Weeks 10–26) & Read Weekly Base Data--------
# -----------------------------------------------------------------------------
week_all <- read.csv(
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
  row.names = NULL
)

week_all_lisa <- left_join(week_all, csl_week, by = c("year", "week"))

# Build Spring Grid strictly for Weeks 10–26
week_spring_grid <- week_all_lisa %>%
  filter(week >= 10, week <= 26) %>%
  complete(year = 2011:2024, week = 10:26) %>%
  left_join(csl_annual_baseline_50yr, by = "year") %>%
  arrange(year, week) %>%
  mutate(
    # Spring-specific radian domain across weeks 10 to 26
    week_rad         = 2 * pi * (week - 10) / (26 - 10),
    year_factor      = factor(year),
    eulachon_input   = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
    log_annual_scale = log(pmax(1, csl_annual_mean))
  )

# -----------------------------------------------------------------------------
# Part 4: Fit Stabilized Spring Models (Weeks 10–26)--------
# -----------------------------------------------------------------------------
# A. Unimodal Quadratic Spring Model (Stabilized Linear Predictor)
eulachon_spring_poly <- glm(
  eulachon_input ~ poly(week, 2, raw = TRUE) + log_annual_scale,
  data    = week_spring_grid %>% filter(!is.na(eulachon_input)),
  family  = quasipoisson(link = "log"),
  control = glm.control(maxit = 100, epsilon = 1e-6)
)

# B. Stabilized Spring GAM (fixed Tweedie p = 1.5 avoids C-level optimizer segfaults)
eulachon_spring_gam <- gam(
  eulachon_input ~ s(week, k = 6) + s(year_factor, bs = "re"),
  data    = week_spring_grid %>% filter(!is.na(eulachon_input)),
  family  = tw(p = 1.5, link = "log"),
  method  = "REML",
  control = gam.control(maxit = 200)
)

summary(eulachon_spring_gam)

# -----------------------------------------------------------------------------
# Part 5: Predict Spring Phenology Signal & Compare Models--------
# -----------------------------------------------------------------------------
spring_eval_grid <- tibble(
  week             = 10:26,
  year_factor      = factor("2011", levels = levels(week_spring_grid$year_factor)),
  log_annual_scale = mean(week_spring_grid$log_annual_scale, na.rm = TRUE)
)

# Extract predictions across Weeks 10–26
spring_eval_grid$poly_pred <- predict(eulachon_spring_poly, newdata = spring_eval_grid, type = "response")
spring_eval_grid$gam_pred  <- exp(predict(eulachon_spring_gam, newdata = spring_eval_grid, type = "link", exclude = "s(year_factor)"))

# -----------------------------------------------------------------------------
# Part 6: Plot Spring Window Model Comparison (Weeks 10–26)--------
# -----------------------------------------------------------------------------
p_spring_comp <- ggplot() +
  geom_point(
    data = week_spring_grid,
    aes(x = week, y = eulachon_input),
    alpha = 0.25, color = "darkgray"
  ) +
  geom_line(
    data = spring_eval_grid,
    aes(x = week, y = poly_pred, color = "Quadratic Polynomial GLM"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  geom_line(
    data = spring_eval_grid,
    aes(x = week, y = gam_pred, color = "Spring GAM (Weeks 10-26)"),
    linewidth = 1.2, linetype = "solid"
  ) +
  scale_x_continuous(breaks = 10:26, name = "Calendar Week Number (Spring Window: 10–26)") +
  scale_color_manual(
    values = c("Quadratic Polynomial GLM" = "firebrick", "Spring GAM (Weeks 10-26)" = "steelblue")
  ) +
  labs(
    title = "Spring Window Eulachon Phenology Fit (Weeks 10–26)",
    subtitle = "Comparing quadratic GLM vs. flexible GAM fit for spring disaggregation",
    y = "Eulachon Biomass / Run Index",
    color = "Model Fit"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p_spring_comp)

ggsave(
  filename = file.path(output_dir, "eulachon_spring_poly_vs_gam.png"),
  plot     = p_spring_comp,
  width    = 10, height = 5
)

# -----------------------------------------------------------------------------
# Part 6: Plot Spring Window Model Comparison (Weeks 10–26)--------
# -----------------------------------------------------------------------------
p_spring_comp <- ggplot() +
  geom_point(
    data = week_spring_grid,
    aes(x = week, y = eulachon_input),
    alpha = 0.25, color = "darkgray"
  ) +
  geom_line(
    data = spring_eval_grid,
    aes(x = week, y = sine_pred, color = "Harmonic Sine Model (2-Harmonics)"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  geom_line(
    data = spring_eval_grid,
    aes(x = week, y = gam_pred, color = "Spring GAM (Weeks 10-26)"),
    linewidth = 1.2, linetype = "solid"
  ) +
  scale_x_continuous(breaks = 10:26, name = "Calendar Week Number (Spring Window: 10–26)") +
  scale_color_manual(
    values = c("Harmonic Sine Model (2-Harmonics)" = "firebrick", "Spring GAM (Weeks 10-26)" = "steelblue")
  ) +
  labs(
    title = "Spring Window Eulachon Phenology Fit (Weeks 10–26)",
    subtitle = "Comparing parametric harmonic vs. flexible GAM fit for spring disaggregation",
    y = "Eulachon Biomass / Run Index",
    color = "Model Fit"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p_spring_comp)

ggsave(
  filename = file.path(output_dir, "eulachon_spring_sine_vs_gam.png"),
  plot     = p_spring_comp,
  width    = 10, height = 5
)

# -----------------------------------------------------------------------------
# Part 7: Reconstruct Weekly CSL Spring Abundances (Weeks 10–26)--------
# -----------------------------------------------------------------------------
# Fit weekly CSL GAM for Weeks 10–26
csl_train_spring <- week_spring_grid %>% 
  filter(!is.na(csl_nonpup_total_emb_lisa)) %>%
  mutate(year_factor = droplevels(year_factor))

csl_spring_gam <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, k = 6) + s(year_factor, bs = "re"),
  data    = csl_train_spring,
  family  = quasipoisson(link = "log"),
  method  = "REML"
)

# Extract relative spring phenology shape
csl_spring_eval <- tibble(
  week        = 10:26,
  year_factor = factor(levels(csl_train_spring$year_factor)[1], levels = levels(csl_train_spring$year_factor))
)

csl_spring_eval$raw_pred <- exp(predict(csl_spring_gam, newdata = csl_spring_eval, type = "link", exclude = "s(year_factor)"))

# Normalize baseline shape relative to the Spring Mean across weeks 10–26
csl_spring_eval <- csl_spring_eval %>%
  mutate(
    # Relative multiplier where average across weeks 10-26 equals 1.0
    csl_spring_relative_signal = raw_pred / mean(raw_pred)
  ) %>%
  select(week, csl_spring_relative_signal)

# Build master unscaled spring database (1993–2024 x Weeks 10–26)
eulachon_master_index <- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"))

eulachon_master_clean <- eulachon_master_index %>%
  mutate(eulachon_lbs_reconstructed = zoo::na.approx(eulachon_lbs_reconstructed, na.rm = FALSE))

# Eulachon spring signal proportions (summing to 1 across weeks 10-26)
eul_spring_signal <- spring_eval_grid %>%
  mutate(eul_prop_signal = gam_pred / sum(gam_pred)) %>%
  select(week, eul_prop_signal)

spring_master_database <- tibble(year = 1993:2024) %>%
  cross_join(tibble(week = 10:26)) %>%
  left_join(csl_annual_baseline_50yr, by = "year") %>%
  left_join(eulachon_master_clean %>% select(year, eulachon_lbs_reconstructed), by = "year") %>%
  left_join(csl_spring_signal, by = "week") %>%
  left_join(eul_spring_signal, by = "week") %>%
  left_join(
    week_spring_grid %>% select(year, week, csl_nonpup_total_emb_lisa, eulachon_ssb_4week_est),
    by = c("year", "week")
  ) %>%
  arrange(year, week) %>%
  mutate(
    # CSL modeled counts = Spring mean * relative spring signal
    csl_modeled_counts = csl_annual_mean * csl_spring_relative_signal,
    csl_final_counts   = if_else(is.na(csl_nonpup_total_emb_lisa), csl_modeled_counts, as.numeric(csl_nonpup_total_emb_lisa)),
    
    # Eulachon modeled pounds = Reconstructed annual total * spring proportion signal
    eulachon_modeled_lbs = eulachon_lbs_reconstructed * eul_prop_signal,
    eulachon_obs_interp  = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE),
    eulachon_final_lbs   = if_else(!is.na(eulachon_obs_interp), eulachon_obs_interp, eulachon_modeled_lbs)
  ) %>%
  select(
    year, week,
    csl_spring_mean_counts = csl_annual_mean,
    csl_obs_counts         = csl_nonpup_total_emb_lisa,
    csl_modeled_counts,
    csl_final_counts,
    eulachon_annual_ssb_lbs = eulachon_lbs_reconstructed,
    eulachon_ssb_obs_lbs    = eulachon_ssb_4week_est,
    eulachon_final_lbs
  )

write.csv(
  spring_master_database,
  file.path(output_dir, "weekly_unscaled_abundance_spring_10_26_1993_2024.csv"),
  row.names = FALSE
)

cat("Successfully generated and saved Spring Window (Weeks 10–26) database to:", 
    file.path(output_dir, "weekly_unscaled_abundance_spring_10_26_1993_2024.csv"), "\n")