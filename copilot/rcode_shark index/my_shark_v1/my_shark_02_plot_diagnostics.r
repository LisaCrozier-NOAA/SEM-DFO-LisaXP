suppressPackageStartupMessages({
  library(tidyverse)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Requires:
# ranked_restrict (or res_shark_restrict with aic)
# If only ranked_restrict exists, this script computes within-group AIC weights.

if (!exists("ranked_restrict")) stop("ranked_restrict not found in environment.")

# ---------------------------
# 1) Prep model table
# ---------------------------
mod_tbl <- ranked %>%
  #mod_tbl <- ranked_restrict %>%
  filter(Q10 == 2, n %in% c(19, 24)) %>%
  mutate(model_group = case_when(
    n == 19 ~ "short",
    n == 24 ~ "long",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(model_group), !is.na(aic))

if (nrow(mod_tbl) == 0) stop("No rows after filtering (Q10==2 and n in {19,24}).")

# within-group Akaike weights
mod_tbl <- mod_tbl %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

# Optional: focus on plausible support only
# (drop extreme low-weight noise)
mod_tbl_plot <- mod_tbl %>%
  filter(aic_weight >= 0.01)

# ---------------------------
# 2) Construct canonical curves for overlap O(T)
# ---------------------------
# T is standardized shark temperature driver
T_grid <- tibble(T = seq(-2.5, 2.5, length.out = 200))

curve_tbl <- mod_tbl_plot %>%
  select(model_group, shark_col, shark_temp_col, overlap_form, overlap_slope,
         shark_transform, aic_weight, delta_aic, n) %>%
  distinct() %>%
  mutate(curve_id = row_number()) %>%
  tidyr::crossing(T_grid) %>%
  mutate(
    O_T = case_when(
      overlap_form == "constant" ~ 1,
      overlap_form == "linear"   ~ pmax(0.1, 1 + overlap_slope * T),
      overlap_form == "logistic" ~ plogis(overlap_slope * T) * 2,
      TRUE ~ NA_real_
    ),
    # normalized within-curve for visual comparability
    O_T_norm = O_T / mean(O_T, na.rm = TRUE)
  )

# ---------------------------
# 3) Weighted mean + band by group/form
# ---------------------------
summary_curve <- curve_tbl %>%
  group_by(model_group, overlap_form, T) %>%
  summarise(
    w_mean = weighted.mean(O_T_norm, w = aic_weight, na.rm = TRUE),
    w_q10  = {
      o <- order(O_T_norm)
      x <- O_T_norm[o]; w <- aic_weight[o] / sum(aic_weight[o], na.rm = TRUE)
      x[which(cumsum(w) >= 0.10)[1]]
    },
    w_q90  = {
      o <- order(O_T_norm)
      x <- O_T_norm[o]; w <- aic_weight[o] / sum(aic_weight[o], na.rm = TRUE)
      x[which(cumsum(w) >= 0.90)[1]]
    },
    .groups = "drop"
  )

# ---------------------------
# 4) Plot: individual weighted curves + weighted summary
# ---------------------------
p_curves <- ggplot() +
  geom_line(
    data = curve_tbl,
    aes(
      x = T, y = O_T_norm, group = curve_id,
      color = overlap_form,
      alpha = aic_weight,
      linewidth = aic_weight
    )
  ) +
  facet_wrap(~model_group, ncol = 1) +
  scale_alpha_continuous(range = c(0.08, 0.9), guide = "none") +
  scale_linewidth_continuous(range = c(0.2, 1.8), guide = "none") +
  theme_bw() +
  labs(
    title = "Shark overlap functional relationship O(T): model family",
    subtitle = "Each line = one model; alpha/width weighted by within-group AIC weight",
    x = "Standardized shark temperature driver (T)",
    y = "Normalized overlap O(T) / mean(O(T))",
    color = "Overlap form"
  )

p_summary <- ggplot(summary_curve, aes(T, w_mean, color = overlap_form, fill = overlap_form)) +
  geom_ribbon(aes(ymin = w_q10, ymax = w_q90), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~model_group, ncol = 1) +
  theme_bw() +
  labs(
    title = "Weighted mean shark overlap curves with uncertainty envelope",
    subtitle = "10th–90th weighted quantile band across supported models",
    x = "Standardized shark temperature driver (T)",
    y = "Weighted normalized overlap"
  )

# ---------------------------
# 5) Supporting barplots for selected discrete choices
# ---------------------------
choice_imp <- mod_tbl %>%
  group_by(model_group, overlap_form) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  group_by(model_group) %>%
  mutate(w_rel = w / max(w)) %>%
  ungroup()

p_choice <- ggplot(choice_imp, aes(overlap_form, w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "AIC-weighted support for overlap form",
    x = "Overlap form",
    y = "Summed AIC weight"
  )

# ---------------------------
# 6) Save outputs
# ---------------------------
write_csv(mod_tbl, file.path(out_dir, "shark_weighted_models_Q10_2_n19_n24.csv"))
write_csv(curve_tbl, file.path(out_dir, "shark_overlap_curves_weighted_points.csv"))
write_csv(summary_curve, file.path(out_dir, "shark_overlap_curves_weighted_summary.csv"))
write_csv(choice_imp, file.path(out_dir, "shark_overlap_form_importance_Q10_2.csv"))

ggsave(file.path(out_dir, "shark_overlap_curves_weighted_family.png"), p_curves, width = 10, height = 8, dpi = 150)
ggsave(file.path(out_dir, "shark_overlap_curves_weighted_summary.png"), p_summary, width = 10, height = 8, dpi = 150)
ggsave(file.path(out_dir, "shark_overlap_form_importance_Q10_2.png"), p_choice, width = 7, height = 5, dpi = 150)


print(p_curves)
print(p_summary)
print(p_choice)

message("Saved weighted functional-relationship plots and tables to: ", out_dir)

best_mod<- mod_tbl %>% group_by(model_group) %>% slice(1)
t(best_mod)
