# parameters
# k0 <- 0.0
# kW <- 1.0   # warm effect on switching (+)
# kF <- 1.0   # forage effect reducing switching (+ in formula as -kF*F_t)
# Q10 <- 2.0
# p_switch <- plogis(k0 + kW * design$W_use - kF * design$F_z)

#TALLY SIGNIFICANCE--------
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


#VARIABLE IMPORTANCE--------
suppressPackageStartupMessages({
  library(tidyverse)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("results_ranked")) stop("results_ranked not found. Run fitting first.")

# 1) Keep target groups and compute within-group AIC weights
res_w <- results_ranked %>%
  mutate(
    model_group = case_when(
      n == 19 ~ "short",
      n == 24 ~ "long",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(model_group), !is.na(aic)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

write_csv(res_w, file.path(out_dir, "phase1_results_within_group_aic_weights.csv"))

# helper to compute selection importance for one selector column
selection_importance <- function(df, selector_col) {
  selector_sym <- rlang::sym(selector_col)
  
  # summed Akaike weight by selector level (Burnham-Anderson)
  imp <- df %>%
    group_by(model_group, !!selector_sym) %>%
    summarise(
      importance_aicw = sum(aic_weight, na.rm = TRUE),
      n_models = n(),
      mean_aic = mean(aic, na.rm = TRUE),
      min_aic = min(aic, na.rm = TRUE),
      mean_r2 = mean(r2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(level = !!selector_sym) %>%
    arrange(model_group, desc(importance_aicw), min_aic)
  
  # normalize within each group so top level = 1
  imp <- imp %>%
    group_by(model_group) %>%
    mutate(
      importance_rel_top1 = importance_aicw / max(importance_aicw, na.rm = TRUE),
      rank = rank(-importance_aicw, ties.method = "min")
    ) %>%
    ungroup() %>%
    mutate(selector = selector_col) %>%
    select(selector, model_group, rank, level, importance_aicw, importance_rel_top1,
           n_models, min_aic, mean_aic, mean_r2)
  
  imp
}

# 2) Compute importance tables for each selector set------
imp_sst <- selection_importance(res_w, "sst_col")
imp_ssl <- selection_importance(res_w, "ssl_col")
imp_shk <- selection_importance(res_w, "shark_col")
imp_for <- selection_importance(res_w, "forcing")

# 3) Save separate and combined
write_csv(imp_sst, file.path(out_dir, "phase1_importance_sst_col_short_long.csv"))
write_csv(imp_ssl, file.path(out_dir, "phase1_importance_ssl_col_short_long.csv"))
write_csv(imp_shk, file.path(out_dir, "phase1_importance_shark_col_short_long.csv"))
write_csv(imp_for, file.path(out_dir, "phase1_importance_forcing_short_long.csv"))

imp_all <- bind_rows(imp_sst, imp_ssl, imp_shk, imp_for)
write_csv(imp_all, file.path(out_dir, "phase1_importance_all_selectors_short_long.csv"))

# 4) Optional: top levels by selector and group
top_levels <- imp_all %>%
  group_by(selector, model_group) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

write_csv(top_levels, file.path(out_dir, "phase1_importance_top5_each_selector_short_long.csv"))

# 5) Quick plots
plot_imp <- function(df, selector_name) {
  ggplot(df, aes(x = reorder(level, importance_aicw), y = importance_aicw, fill = model_group)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste0("Selection importance by AIC weight: ", selector_name),
      x = selector_name,
      y = "Summed AIC weight (within-group)"
    )
}

p_sst <- plot_imp(imp_sst, "sst_col")
p_ssl <- plot_imp(imp_ssl, "ssl_col")
p_shk <- plot_imp(imp_shk, "shark_col")
p_for <- plot_imp(imp_for, "forcing")

ggsave(file.path(out_dir, "phase1_importance_sst_col_short_long.png"), p_sst, width = 10, height = 7, dpi = 150)
ggsave(file.path(out_dir, "phase1_importance_ssl_col_short_long.png"), p_ssl, width = 10, height = 7, dpi = 150)
ggsave(file.path(out_dir, "phase1_importance_shark_col_short_long.png"), p_shk, width = 10, height = 7, dpi = 150)
ggsave(file.path(out_dir, "phase1_importance_forcing_short_long.png"), p_for, width = 8, height = 5, dpi = 150)

cat("\nTop selector levels by group (AIC-weighted importance):\n")
print(top_levels)
