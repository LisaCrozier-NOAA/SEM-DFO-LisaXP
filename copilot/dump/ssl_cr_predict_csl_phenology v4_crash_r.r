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
# Part 1: Load Historical yearly Spring CSL Mean count (Weeks 10–26) and eulachon ssb reconstruction--------
# -----------------------------------------------------------------------------
#source(file.path(output_dir, "ssl_cr_data_sources"))
csl_annual<-csl_annual_baseline_50yr <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),row.names = NULL) %>%
              filter(year>=1993)

eulachon_annual<-eulachon_master_index <- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"))


# -----------------------------------------------------------------------------
# Part 3: Expand Grid to Full 52 Weeks (2011–2024)--------
# -----------------------------------------------------------------------------
week_all <- read.csv(file.path(output_dir,"")
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
  row.names = NULL
)

week_all_lisa <- left_join(week_all, csl_week, by = c("year", "week"))

# Expand dataset to full 52 weeks per year
week_full_year <- week_all_lisa %>%
  complete(year = 2011:2024, week = 1:52) %>%
  left_join(csl_annual_baseline_50yr, by = "year") %>%
  arrange(year, week) %>%
  mutate(
    week_rad         = 2 * pi * (week - 1) / 52,
    year_factor      = factor(year),
    eulachon_input   = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
    log_annual_scale = log(pmax(1, csl_annual_mean))
  )

# -----------------------------------------------------------------------------
# Part 4: Fit 52-Week Eulachon Cyclic GAM & Extract Baseline Signal--------
# -----------------------------------------------------------------------------
eul_train_data <- week_full_year %>% 
  filter(!is.na(eulachon_input)) %>%
  mutate(year_factor = droplevels(year_factor))

# Cyclic GAM over 52 weeks (Fixed Tweedie p = 1.5 prevents C optimizer crashes)
eulachon_52wk_gam <- gam(
  eulachon_input ~ s(week, bs = "cc", k = 10) + s(year_factor, bs = "re"),
  data    = eul_train_data,
  family  = tw(p = 1.5, link = "log"),
  knots   = list(week = c(1, 52)),
  method  = "REML"
)

# Extract 52-week population signal (proportions summing to 1.0)
eul_52wk_eval <- tibble(
  week        = 1:52,
  year_factor = factor(levels(eul_train_data$year_factor)[1], levels = levels(eul_train_data$year_factor))
)

eul_52wk_eval$raw_pred <- exp(predict(eulachon_52wk_gam, newdata = eul_52wk_eval, type = "link", exclude = "s(year_factor)"))

eul_52wk_signal <- eul_52wk_eval %>%
  mutate(eul_prop_signal = raw_pred / sum(raw_pred)) %>%
  select(week, eul_prop_signal)

# -----------------------------------------------------------------------------
# Part 5: Fit 52-Week CSL Cyclic GAM & Extract Spring-Mean Relative Signal--------
# -----------------------------------------------------------------------------
csl_train_data <- week_full_year %>% 
  filter(!is.na(csl_nonpup_total_emb_lisa)) %>%
  mutate(year_factor = droplevels(year_factor))

csl_52wk_gam <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, bs = "cc", k = 10) + s(year_factor, bs = "re"),
  data    = csl_train_data,
  family  = quasipoisson(link = "log"),
  knots   = list(week = c(1, 52)),
  method  = "REML"
)

csl_52wk_eval <- tibble(
  week        = 1:52,
  year_factor = factor(levels(csl_train_data$year_factor)[1], levels = levels(csl_train_data$year_factor))
)

csl_52wk_eval$raw_pred <- exp(predict(csl_52wk_gam, newdata = csl_52wk_eval, type = "link", exclude = "s(year_factor)"))

# Compute Spring Mean (Weeks 10–26 average) and standardize multiplier
spring_mean_baseline <- mean(csl_52wk_eval$raw_pred[10:26])

csl_52wk_signal <- csl_52wk_eval %>%
  mutate(
    # Relative multiplier where average across spring window (Weeks 10-26) equals 1.0
    csl_spring_relative_signal = raw_pred / spring_mean_baseline
  ) %>%
  select(week, csl_spring_relative_signal)

# -----------------------------------------------------------------------------
# Part 6: Build Full 1993–2024 Master Database & Crop to Weeks 10–26--------
# -----------------------------------------------------------------------------
eulachon_master_index <- read.csv(file.path(output_dir, "eulachon_master_index_lbs_1993_2024.csv"))

eulachon_master_clean <- eulachon_master_index %>%
  mutate(eulachon_lbs_reconstructed = zoo::na.approx(eulachon_lbs_reconstructed, na.rm = FALSE))

# Build full 52-week x 32-year database
full_master_database <- tibble(year = 1993:2024) %>%
  cross_join(tibble(week = 1:52)) %>%
  left_join(csl_annual_baseline_50yr, by = "year") %>%
  left_join(eulachon_master_clean %>% select(year, eulachon_lbs_reconstructed), by = "year") %>%
  left_join(csl_52wk_signal, by = "week") %>%
  left_join(eul_52wk_signal, by = "week") %>%
  left_join(
    week_full_year %>% select(year, week, csl_nonpup_total_emb_lisa, eulachon_ssb_4week_est),
    by = c("year", "week")
  ) %>%
  arrange(year, week) %>%
  mutate(
    # CSL modeled counts = Spring mean (weeks 10-26) * relative spring signal
    csl_modeled_counts = csl_annual_mean * csl_spring_relative_signal,
    csl_final_counts   = if_else(is.na(csl_nonpup_total_emb_lisa), csl_modeled_counts, as.numeric(csl_nonpup_total_emb_lisa)),
    
    # Eulachon modeled pounds = Reconstructed annual total * 52-week proportion signal
    eulachon_modeled_lbs = eulachon_lbs_reconstructed * eul_prop_signal,
    eulachon_obs_interp  = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE),
    eulachon_final_lbs   = if_else(!is.na(eulachon_obs_interp), eulachon_obs_interp, eulachon_modeled_lbs)
  )

# CROP DOWN TO WEEKS 10–26 FOR FINAL OUTPUT
spring_master_database <- full_master_database %>%
  filter(week >= 10, week <= 26) %>%
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

cat("Successfully generated 52-week cyclic models and cropped Spring Window (Weeks 10–26) database to:", 
    file.path(output_dir, "weekly_unscaled_abundance_spring_10_26_1993_2024.csv"), "\n")

# -----------------------------------------------------------------------------
# Part 7: Plot Cropped Spring Window (Weeks 10–26: 1993–2024)--------
# -----------------------------------------------------------------------------
p_spring_cropped <- ggplot(spring_master_database, aes(x = week)) +
  geom_line(
    aes(y = eulachon_final_lbs, color = "Eulachon Biomass (lbs)"),
    linewidth = 1, linetype = "dashed"
  ) +
  geom_point(
    aes(y = csl_obs_counts * 10000, color = "Observed CSL Counts"),
    size = 1.6, na.rm = TRUE
  ) +
  geom_line(
    aes(y = csl_final_counts * 10000, color = "CSL Abundance (Counts)"),
    linewidth = 1
  ) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = seq(10, 26, by = 4), name = "Week Number (Weeks 10–26)") +
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
    title = "Cropped Spring Window Weekly Abundance & Biomass (1993–2024)",
    subtitle = "Fit across 52-week cyclic domain, cropped to Weeks 10–26; CSL back-predictions match spring means",
    color = "Metric"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_spring_cropped)

ggsave(
  filename = file.path(output_dir, "eulachon_csl_spring_cropped_10_26_1993_2024.png"),
  plot     = p_spring_cropped,
  width    = 12, height = 10
)