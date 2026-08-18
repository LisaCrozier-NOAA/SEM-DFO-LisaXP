
rm(list=ls())

#Notes:
# Object Dependencies:
#   
#   Part 1: Loads csl_annual_baseline_1976_2024.csv, eulachon_master_index_lbs_1993_2024.csv, and jake_week_all.csv from copilot/outs_csl_cr. It explicitly handles the 2010 gap in eulachon_master_index_clean using zoo::na.approx() before any joins occur.
# 
# Part 2: Creates week_52_obs, fits eul_52wk_gam, and outputs eul_52wk_signal (which contains week and eul_prop_signal).
# 
# Part 3: Joining eul_52wk_signal in full_year_grid works cleanly because eul_52wk_signal was generated prior to full_year_grid.
# 
# eulachon_hybrid_lbs Calculation:
#   
#   For 2011–2024, it uses eulachon_ssb_4week_est from jake_week_all.csv (interpolated across internal missing weeks like 2020 via zoo::na.approx()), capturing real empirical timing shifts.
# 
# For 1993–2010, it falls back to eulachon_modeled_lbs (eulachon_lbs_reconstructed * eul_prop_signal).
# 
# Factor Levels & predict.gam() Warnings:
#   
#   csl_train_data uses droplevels(year_factor). fit_levels holds only observed years (2011:2018, 2021:2024).
# 
# pred_grid_clean forces unseen years (1993:2010, 2019, 2020) to anchor_level ("2011") while locking levels = fit_levels. Because newdata matches the exact levels of csl_full_year_gam, predict.gam(..., exclude = "s(year_factor)") runs 100% warning-free.
# 
# Normalization & Plotting:
#   
#   Uses scalar max checks (csl_max and eulachon_max) with ifelse() to avoid dplyr::if_else() vector length mismatches or empty-vector -Inf warnings.



# -----------------------------------------------------------------------------
# Setup & Package Loading--------
# -----------------------------------------------------------------------------
output_dir <- "copilot/outs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(tidyverse)
library(readxl)
library(mgcv)
library(zoo)
library(scales)

# -----------------------------------------------------------------------------
# Part 1: Read Baseline Data from Output Directory--------
# -----------------------------------------------------------------------------
csl_annual_baseline_1976_2024 <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"))
eulachon_master_index         <- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"))
week_all                      <- read.csv(file.path(output_dir, "jake_week_all.csv"))

# Ensure 2010 Eulachon biomass gap is interpolated before joining
eulachon_master_index_clean <- eulachon_master_index %>%
  mutate(
    eulachon_lbs_reconstructed = zoo::na.approx(eulachon_lbs_reconstructed, na.rm = FALSE)
  )

# -----------------------------------------------------------------------------
# Part 2: Extract 52-Week Baseline Eulachon Seasonal Signal--------
# -----------------------------------------------------------------------------
# A. Build 52-week observational dataset for GAM fitting
week_52_obs <- week_all %>%
  filter(year >= 2011, year <= 2024) %>%
  complete(year = 2011:2024, week = 1:52) %>%
  arrange(year, week) %>%
  mutate(
    year_factor    = factor(year),
    eulachon_input = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2)
  )

# B. Fit 52-Week Cyclic GAM for Eulachon Seasonal Baseline
eul_52wk_gam <- gam(
  eulachon_input ~ s(week, bs = "cc", k = 10) + s(year_factor, bs = "re"),
  data    = week_52_obs %>% filter(!is.na(eulachon_input)),
  family  = tw(link = "log"),
  knots   = list(week = c(1, 52)),
  method  = "REML"
)

# C. Extract population-wide 52-week baseline signal (proportions summing to 1.0)
eul_52wk_grid <- tibble(
  week        = 1:52,
  year_factor = factor("2011", levels = levels(week_52_obs$year_factor))
)

eul_52wk_grid$fit_link <- predict(eul_52wk_gam, newdata = eul_52wk_grid, type = "link", exclude = "s(year_factor)")

eul_52wk_signal <- eul_52wk_grid %>%
  mutate(
    eul_fit_raw     = exp(fit_link),
    eul_prop_signal = eul_fit_raw / sum(eul_fit_raw)
  ) %>%
  select(week, eul_prop_signal)

# -----------------------------------------------------------------------------
# Part 3: Build Complete 1993–2024 Grid with Hybrid Weekly Eulachon--------
# -----------------------------------------------------------------------------
full_year_grid <- tibble(year = 1993:2024) %>%
  cross_join(tibble(week = 1:52)) %>%
  left_join(csl_annual_baseline_1976_2024, by = "year") %>%
  left_join(eulachon_master_index_clean %>% select(year, eulachon_lbs_reconstructed), by = "year") %>%
  left_join(eul_52wk_signal, by = "week") %>%
  left_join(
    week_all %>% select(year, week, csl_nonpup_total_emb, eulachon_ssb_4week_est),
    by = c("year", "week")
  ) %>%
  arrange(year, week) %>%
  mutate(
    year_factor          = factor(year),
    log_annual_scale     = log(pmax(1, csl_annual_mean)),
    
    # Baseline historical disaggregation for pre-2011 years
    eulachon_modeled_lbs = eulachon_lbs_reconstructed * eul_prop_signal
  ) %>%
  # Interpolate missing internal weekly observational values (e.g., 2020) across modern years
  group_by(year) %>%
  mutate(
    eulachon_obs_interp = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    # HYBRID EULACHON: Use observed weekly data if available (2011-2024);
    # Fall back to reconstructed annual x baseline signal for historical years (1993-2010)
    eulachon_hybrid_lbs = if_else(
      !is.na(eulachon_obs_interp), 
      eulachon_obs_interp, 
      eulachon_modeled_lbs
    )
  )

# -----------------------------------------------------------------------------
# Part 4: Fit CSL GAM & Predict Across 1993–2024--------
# -----------------------------------------------------------------------------
# Drop unused factor levels so GAM model only sees years with observations
csl_train_data <- full_year_grid %>% 
  filter(!is.na(csl_nonpup_total_emb)) %>%
  mutate(year_factor = droplevels(year_factor))

# Fit GAM using week smooth, annual abundance scale, and year random effects
csl_full_year_gam <- gam(
  csl_nonpup_total_emb ~ s(week, bs = "cc", k = 10) + log_annual_scale + s(year_factor, bs = "re"),
  data   = csl_train_data,
  family = quasipoisson(link = "log"),
  knots  = list(week = c(1, 52)),
  method = "REML"
)

# Extract exact levels fitted by model
fit_levels   <- levels(csl_train_data$year_factor)
anchor_level <- fit_levels[1]

# Map unseen years to anchor level and enforce exact training factor levels
pred_grid_clean <- full_year_grid %>%
  mutate(
    year_factor = factor(
      if_else(as.character(year) %in% fit_levels, as.character(year), anchor_level),
      levels = fit_levels
    )
  )

# Predict population baseline across all years (warning-free)
full_year_grid$csl_pred_pop <- predict(
  csl_full_year_gam,
  newdata = pred_grid_clean,
  type    = "response",
  exclude = "s(year_factor)"
)

# Hybrid CSL Curve: Use actual observations where present; fall back to GAM prediction for gaps/historical years
full_year_grid <- full_year_grid %>%
  mutate(
    csl_hybrid = if_else(is.na(csl_nonpup_total_emb), csl_pred_pop, as.numeric(csl_nonpup_total_emb))
  )

# Normalize weekly curves to 0–1 scale per year
full_year_overlap_df <- full_year_grid %>%
  group_by(year) %>%
  mutate(
    csl_max      = max(csl_hybrid, na.rm = TRUE),
    eulachon_max = max(eulachon_hybrid_lbs, na.rm = TRUE),
    
    csl_norm      = ifelse(csl_max == 0 | is.na(csl_max), 0, csl_hybrid / csl_max),
    eulachon_norm = ifelse(eulachon_max == 0 | is.na(eulachon_max), 0, eulachon_hybrid_lbs / eulachon_max)
  ) %>%
  select(-csl_max, -eulachon_max) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Part 5: Plot 52-Week Calendar Year Phenology Overlap--------
# -----------------------------------------------------------------------------
p_52week_overlap <- ggplot(full_year_overlap_df, aes(x = week)) +
  # Raw weekly CSL observations as points (where present)
  geom_point(
    aes(y = ifelse(is.na(csl_nonpup_total_emb), NA_real_, csl_norm), color = "Observed CSL Data"),
    size = 1.8, na.rm = TRUE
  ) +
  # CSL Hybrid Curve (Solid Steelblue)
  geom_line(aes(y = csl_norm, color = "CSL Phenology (Hybrid)"), linewidth = 1) +
  # Eulachon Weekly Hybrid Signal (Dashed Gold)
  geom_line(aes(y = eulachon_norm, color = "Eulachon Hybrid Signal"), linewidth = 1, linetype = "dashed") +
  facet_wrap(~ year, ncol = 4) +
  scale_color_manual(
    values = c(
      "Observed CSL Data"      = "black",
      "CSL Phenology (Hybrid)"  = "steelblue",
      "Eulachon Hybrid Signal"  = "darkgoldenrod3"
    )
  ) +
  scale_x_continuous(breaks = seq(4, 52, by = 12), name = "Calendar Week (1 to 52)") +
  scale_y_continuous(name = "Normalized Relative Scale (0 to 1)") +
  labs(
    title = "Full Calendar Year Phenology Overlap (Weeks 1–52: 1993–2024)",
    subtitle = "Observed weekly data used for both species when available; GAM baseline fills historical years/gaps",
    color = "Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_52week_overlap)

# Save plot to output directory
ggsave(
  filename = file.path(output_dir, "eulachon_csl_52week_overlap_1993_2024.png"),
  plot     = p_52week_overlap,
  width    = 12, height = 10
)


# -----------------------------------------------------------------------------
# Setup & Package Loading--------
# -----------------------------------------------------------------------------
output_dir <- "copilot/outs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(tidyverse)
library(readxl)
library(mgcv)
library(zoo)
library(scales)

# -----------------------------------------------------------------------------
# Part 1: Read Baseline Data from Output Directory--------
# -----------------------------------------------------------------------------
csl_annual_baseline_1976_2024 <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"))
eulachon_master_index         <- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"))
week_all                      <- read.csv(file.path(output_dir, "jake_week_all.csv"))

# Ensure 2010 Eulachon biomass gap is interpolated before joining
eulachon_master_index_clean <- eulachon_master_index %>%
  mutate(
    eulachon_lbs_reconstructed = zoo::na.approx(eulachon_lbs_reconstructed, na.rm = FALSE)
  )

# -----------------------------------------------------------------------------
# Part 2: Extract 52-Week Baseline Eulachon Seasonal Signal--------
# -----------------------------------------------------------------------------
# A. Build 52-week observational dataset for GAM fitting
week_52_obs <- week_all %>%
  filter(year >= 2011, year <= 2024) %>%
  complete(year = 2011:2024, week = 1:52) %>%
  arrange(year, week) %>%
  mutate(
    year_factor    = factor(year),
    eulachon_input = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2)
  )

# B. Fit 52-Week Cyclic GAM for Eulachon Seasonal Baseline
eul_52wk_gam <- gam(
  eulachon_input ~ s(week, bs = "cc", k = 10) + s(year_factor, bs = "re"),
  data    = week_52_obs %>% filter(!is.na(eulachon_input)),
  family  = tw(link = "log"),
  knots   = list(week = c(1, 52)),
  method  = "REML"
)

# C. Extract population-wide 52-week baseline signal (proportions summing to 1.0)
eul_52wk_grid <- tibble(
  week        = 1:52,
  year_factor = factor("2011", levels = levels(week_52_obs$year_factor))
)

eul_52wk_grid$fit_link <- predict(eul_52wk_gam, newdata = eul_52wk_grid, type = "link", exclude = "s(year_factor)")

eul_52wk_signal <- eul_52wk_grid %>%
  mutate(
    eul_fit_raw     = exp(fit_link),
    eul_prop_signal = eul_fit_raw / sum(eul_fit_raw)
  ) %>%
  select(week, eul_prop_signal)

# -----------------------------------------------------------------------------
# Part 3: Build Complete 1993–2024 Grid with Hybrid Weekly Eulachon--------
# -----------------------------------------------------------------------------
full_year_grid <- tibble(year = 1993:2024) %>%
  cross_join(tibble(week = 1:52)) %>%
  left_join(csl_annual_baseline_1976_2024, by = "year") %>%
  left_join(eulachon_master_index_clean %>% select(year, eulachon_lbs_reconstructed), by = "year") %>%
  left_join(eul_52wk_signal, by = "week") %>%
  left_join(
    week_all %>% select(year, week, csl_nonpup_total_emb, eulachon_ssb_4week_est),
    by = c("year", "week")
  ) %>%
  arrange(year, week) %>%
  mutate(
    year_factor          = factor(year),
    log_annual_scale     = log(pmax(1, csl_annual_mean)),
    
    # Historical disaggregation using back-transformed SSB mapping for pre-2011 years
    eulachon_modeled_lbs = eulachon_lbs_reconstructed * eul_prop_signal
  ) %>%
  # Interpolate missing internal weekly observational values across modern years
  group_by(year) %>%
  mutate(
    eulachon_obs_interp = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    # HYBRID EULACHON: Use observed weekly data if available (2011-2024);
    # Fall back to SSB-mapped annual x baseline signal for historical years (1993-2010)
    eulachon_hybrid_lbs = if_else(
      !is.na(eulachon_obs_interp), 
      eulachon_obs_interp, 
      eulachon_modeled_lbs
    )
  )

# -----------------------------------------------------------------------------
# Part 4: Fit CSL GAM & Predict Across 1993–2024--------
# -----------------------------------------------------------------------------
# Drop unused factor levels so GAM model only sees years with observations
csl_train_data <- full_year_grid %>% 
  filter(!is.na(csl_nonpup_total_emb)) %>%
  mutate(year_factor = droplevels(year_factor))

# Fit GAM using week smooth, annual abundance scale, and year random effects
csl_full_year_gam <- gam(
  csl_nonpup_total_emb ~ s(week, bs = "cc", k = 10) + log_annual_scale + s(year_factor, bs = "re"),
  data   = csl_train_data,
  family = quasipoisson(link = "log"),
  knots  = list(week = c(1, 52)),
  method = "REML"
)

# Extract exact levels fitted by model
fit_levels   <- levels(csl_train_data$year_factor)
anchor_level <- fit_levels[1]

# Map unseen years to anchor level and enforce exact training factor levels
pred_grid_clean <- full_year_grid %>%
  mutate(
    year_factor = factor(
      if_else(as.character(year) %in% fit_levels, as.character(year), anchor_level),
      levels = fit_levels
    )
  )

# Predict population baseline across all years (warning-free)
full_year_grid$csl_pred_pop <- predict(
  csl_full_year_gam,
  newdata = pred_grid_clean,
  type    = "response",
  exclude = "s(year_factor)"
)

# Hybrid CSL Curve: Use actual observations where present; fall back to GAM prediction for gaps/historical years
full_year_grid <- full_year_grid %>%
  mutate(
    csl_hybrid = if_else(is.na(csl_nonpup_total_emb), csl_pred_pop, as.numeric(csl_nonpup_total_emb))
  )

# Normalize weekly curves to 0–1 scale per year
full_year_overlap_df <- full_year_grid %>%
  group_by(year) %>%
  mutate(
    csl_max      = max(csl_hybrid, na.rm = TRUE),
    eulachon_max = max(eulachon_hybrid_lbs, na.rm = TRUE),
    
    csl_norm      = ifelse(csl_max == 0 | is.na(csl_max), 0, csl_hybrid / csl_max),
    eulachon_norm = ifelse(eulachon_max == 0 | is.na(eulachon_max), 0, eulachon_hybrid_lbs / eulachon_max)
  ) %>%
  select(-csl_max, -eulachon_max) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Part 5: Plot 52-Week Calendar Year Phenology Overlap--------
# -----------------------------------------------------------------------------
p_52week_overlap <- ggplot(full_year_overlap_df, aes(x = week)) +
  geom_point(
    aes(y = ifelse(is.na(csl_nonpup_total_emb), NA_real_, csl_norm), color = "Observed CSL Data"),
    size = 1.8, na.rm = TRUE
  ) +
  geom_line(aes(y = csl_norm, color = "CSL Phenology (Hybrid)"), linewidth = 1) +
  geom_line(aes(y = eulachon_norm, color = "Eulachon Hybrid Signal"), linewidth = 1, linetype = "dashed") +
  facet_wrap(~ year, ncol = 4) +
  scale_color_manual(
    values = c(
      "Observed CSL Data"      = "black",
      "CSL Phenology (Hybrid)"  = "steelblue",
      "Eulachon Hybrid Signal"  = "darkgoldenrod3"
    )
  ) +
  scale_x_continuous(breaks = seq(4, 52, by = 12), name = "Calendar Week (1 to 52)") +
  scale_y_continuous(name = "Normalized Relative Scale (0 to 1)") +
  labs(
    title = "Full Calendar Year Phenology Overlap (Weeks 1–52: 1993–2024)",
    subtitle = "Observed weekly data used for both species when available; GAM baseline fills historical years/gaps",
    color = "Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_52week_overlap)

ggsave(
  filename = file.path(output_dir, "eulachon_csl_52week_overlap_1993_2024.png"),
  plot     = p_52week_overlap,
  width    = 12, height = 10
)

# -----------------------------------------------------------------------------
# Part 6: Construct and Save Unscaled Weekly Abundance Dataset--------
# -----------------------------------------------------------------------------
# Formulate full 52-week x 32-year continuous database in real units
weekly_unscaled_abundance_1993_2024 <- full_year_grid %>%
  select(
    year,
    week,
    # CSL Actual Units (Counts)
    csl_obs_counts          = csl_nonpup_total_emb,
    csl_pred_pop_counts     = csl_pred_pop,
    csl_final_counts        = csl_hybrid,
    
    # Eulachon Actual Units (Pounds)
    eulachon_annual_ssb_lbs = eulachon_lbs_reconstructed, # Back-transformed SSB mapping
    eulachon_ssb_obs_lbs    = eulachon_ssb_4week_est,    # Weekly observed 4-week SSB
    eulachon_final_lbs      = eulachon_hybrid_lbs        # Final hybrid weekly biomass in lbs
  ) %>%
  mutate(
    # Column flag for data origin
    eulachon_data_source = if_else(year <= 2010, "Historical (SSB Back-Mapped)", "Modern Observed/Interpolated"),
    csl_data_source      = if_else(is.na(csl_obs_counts), "GAM Predicted Baseline", "Observed Field Count")
  )

cat("--- UNSCALLED WEEKLY ABUNDANCE DATABASE SUMMARY (1993–2024) ---\n")
cat("Total Rows (52 weeks x 32 years):", nrow(weekly_unscaled_abundance_1993_2024), "\n")
print(head(weekly_unscaled_abundance_1993_2024, 15))

# Save unscaled database to output directory
write.csv(
  weekly_unscaled_abundance_1993_2024,
  file.path(output_dir, "weekly_unscaled_abundance_eulachon_csl_1993_2024.csv"),
  row.names = FALSE
)

cat("\nSaved unscaled weekly abundance database to:", file.path(output_dir, "weekly_unscaled_abundance_eulachon_csl_1993_2024.csv"), "\n")


# -----------------------------------------------------------------------------
# Part 7: Plot Spring Window (Apr–Jun / Weeks 14–26) Weekly Biomass & Counts--------
# -----------------------------------------------------------------------------
# Filter unscaled database strictly to Weeks 14–26 (April through June)
spring_unscaled_df <- weekly_unscaled_abundance_1993_2024 %>%
  filter(week >= 14, week <= 26)

# Dual-axis plot: Eulachon Biomass in Pounds (Left Y-Axis) & CSL Abundance (Right Y-Axis)
p_spring_unscaled <- ggplot(spring_unscaled_df, aes(x = week)) +
  # Eulachon Weekly Biomass (Pounds) - Left Axis
  geom_line(
    aes(y = eulachon_final_lbs, color = "Eulachon Biomass (lbs)"),
    linewidth = 1, linetype = "dashed"
  ) +
  # Observed CSL Points where available
  geom_point(
    aes(y = csl_obs_counts * 10000, color = "Observed CSL Counts"),
    size = 1.6, na.rm = TRUE
  ) +
  # CSL Weekly Abundance (Counts) - Scaled for secondary Y-axis overlay
  geom_line(
    aes(y = csl_final_counts * 10000, color = "CSL Abundance (Counts)"),
    linewidth = 1
  ) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = seq(14, 26, by = 4), name = "Week Number (Apr–Jun: Weeks 14–26)") +
  # Left Y-Axis for Eulachon Pounds, Right Y-Axis transformation for CSL Counts
  scale_y_continuous(
    name = "Eulachon Biomass (Pounds)",
    labels = scales::comma,
    sec.axis = sec_axis(~ . / 10000, name = "CSL Non-Pup Count", labels = scales::comma)
  ) +
  scale_color_manual(
    values = c(
      "Eulachon Biomass (lbs)" = "darkgoldenrod3",
      "CSL Abundance (Counts)" = "steelblue",
      "Observed CSL Counts"    = "black"
    )
  ) +
  labs(
    title = "Spring Window Weekly Abundance & Biomass (1993–2024: Apr–Jun)",
    subtitle = "Actual units: Eulachon in Pounds (dashed gold) vs. CSL Non-Pup Counts (blue/black points)",
    color = "Metric"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_spring_unscaled)

# Save Spring window plot
ggsave(
  filename = file.path(output_dir, "eulachon_csl_spring_unscaled_abundance_1993_2024.png"),
  plot     = p_spring_unscaled,
  width    = 12, height = 10
)