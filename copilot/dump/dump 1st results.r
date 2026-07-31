
suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# Requires in environment:
#   results_ranked
# Optional (for functional plots + refit):
#   fit_one_combo, ak_yr, sem_yr, shark_yr, clim_tr, forage_ak_use, forage_sem_use,
#   col_salmon_juv, col_salmon_adult, safe_scale, out_dir

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("results_ranked")) stop("results_ranked not found. Run fitting script first.")

# ---- 1) Define groups ----
res <- results_ranked %>%
  mutate(
    model_group = case_when(
      n == 19 ~ "short",
      n == 24 ~ "long",
      TRUE    ~ NA_character_
    ),
    sig_cpue  = !is.na(p_cpue)  & p_cpue  < 0.05,
    sig_ssl   = !is.na(p_ssl)   & p_ssl   < 0.05,
    sig_shark = !is.na(p_shark) & p_shark < 0.05,
    sig_all3  = sig_cpue & sig_ssl & sig_shark
  ) %>%
  filter(!is.na(model_group))

if (nrow(res) == 0) stop("No models with n==19 or n==24 found in results_ranked.")

# ---- 2) Significance tally table ----
sig_tally <- res %>%
  group_by(model_group, n) %>%
  summarise(
    n_models = n(),
    n_sig_cpue = sum(sig_cpue, na.rm = TRUE),
    n_sig_ssl = sum(sig_ssl, na.rm = TRUE),
    n_sig_shark = sum(sig_shark, na.rm = TRUE),
    n_sig_all3 = sum(sig_all3, na.rm = TRUE),
    pct_sig_cpue = 100 * n_sig_cpue / n_models,
    pct_sig_ssl = 100 * n_sig_ssl / n_models,
    pct_sig_shark = 100 * n_sig_shark / n_models,
    pct_sig_all3 = 100 * n_sig_all3 / n_models,
    .groups = "drop"
  )

write_csv(sig_tally, file.path(out_dir, "phase1_sig_tally_short_vs_long.csv"))

# model_group     n n_models n_sig_cpue n_sig_ssl n_sig_shark n_sig_all3 pct_sig_cpue pct_sig_ssl pct_sig_shark pct_sig_all3
# 1 long           24      288        264       105          57         40         91.7        36.5          19.8        13.9 
# 2 short          19      288         88       190          38         13         30.6        66.0          13.2         4.51

#DUMP BELOW???-----------------------
# ---- 3) Burnham-Anderson AIC weights within each group ----
# Relative to best AIC in that group
weights_tbl <- res %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

write_csv(weights_tbl, file.path(out_dir, "phase1_results_with_aic_weights.csv"))

# ---- 4) Variable importance (sum of AIC weights where term is significant) ----
var_importance <- weights_tbl %>%
  group_by(model_group) %>%
  summarise(
    vi_cpue_sig  = sum(aic_weight[sig_cpue], na.rm = TRUE),
    vi_ssl_sig   = sum(aic_weight[sig_ssl], na.rm = TRUE),
    vi_shark_sig = sum(aic_weight[sig_shark], na.rm = TRUE),
    vi_all3_sig  = sum(aic_weight[sig_all3], na.rm = TRUE),
    .groups = "drop"
  )

write_csv(var_importance, file.path(out_dir, "phase1_variable_importance_aicw_short_vs_long.csv"))

# ---- 5) Additional "presence-only" VI (all models include all 3 terms, so should be ~1) ----
# Kept for completeness/transparency.
var_presence_importance <- weights_tbl %>%
  group_by(model_group) %>%
  summarise(
    vi_cpue_present  = sum(aic_weight, na.rm = TRUE),
    vi_ssl_present   = sum(aic_weight, na.rm = TRUE),
    vi_shark_present = sum(aic_weight, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(var_presence_importance, file.path(out_dir, "phase1_variable_presence_importance.csv"))

# ---- 6) Best significant models (all 3 significant) in each group ----
best_sig_models <- weights_tbl %>%
  filter(sig_all3) %>%
  group_by(model_group) %>%
  arrange(aic, desc(r2), .by_group = TRUE) %>%*
  