suppressPackageStartupMessages({
  library(tidyverse)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs_5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Requires:
# ranked

if (!exists("ranked")) stop("ranked_restrict not found in environment.")

# ---------------------------
# 1) Prep model table
# ---------------------------
mod_tbl <- ranked %>%
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

curve_tbl <- mod_tbl %>%
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
  group_by( overlap_form, T) %>%
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

library(tidyverse)
# library(ggrepel) # Optional: uncomment if labels crowd/overlap each other

# 1. Isolate the rightmost point of each curve to anchor the text labels
right_labels <- curve_tbl %>%
  group_by(curve_id) %>%
  filter(T == max(T)) %>%
  ungroup()

p_curves <- ggplot() +
  # Main curves
  geom_line(
    data = curve_tbl,
    aes(
      x = T, y = O_T_norm, group = curve_id,
      color = overlap_form,
      alpha = aic_weight,
      linewidth = aic_weight
    )
  ) +
  # Right-side text labels
  geom_text(
    data = right_labels,
    aes(
      x = T, 
      y = O_T_norm, 
      label = round(overlap_slope, 2), # Displays the slope value
      color = overlap_form
    ),
    hjust = -0.2, # Pushes text slightly past the end of the line
    size = 3,
    show.legend = FALSE
  ) +
  scale_alpha_continuous(range = c(0.08, 0.9), guide = "none") +
  scale_linewidth_continuous(range = c(0.2, 1.8), guide = "none") +
  # Allow plot elements to spill into the right margin without getting clipped
  coord_cartesian(clip = "off") + 
  theme_bw() +
  theme(
    # Expand right plot margin to give the text room to breathe
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10) 
  ) +
  labs(
    title = "Shark overlap functional relationship \n O(T): model family",
    subtitle = "Each line = one model \nalpha/width weighted by within-group AIC weight",
    x = "Standardized shark temperature driver (T)",
    y = "Normalized overlap O(T) / mean(O(T))",
    color = "Overlap form"
  )

slope_and_shape_curves<-p_curves

p_summary <- ggplot(summary_curve, aes(T, w_mean, color = overlap_form, fill = overlap_form)) +
  geom_ribbon(aes(ymin = w_q10, ymax = w_q90), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
 # facet_wrap(~model_group, ncol = 1) +
  theme_bw() +
  labs(
    title = "Weighted mean shark overlap curves with uncertainty envelope",
    subtitle = "10th–90th weighted quantile band across supported models",
    x = "Standardized shark temperature driver (T)",
    y = "Weighted normalized overlap"
  )

# ---------------------------
# 5) Supporting barplots for selected discrete choices
# overlap shape---------------------------
overlap_slope_imp <- mod_tbl %>%
  group_by(model_group, overlap_slope) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  group_by(model_group) %>%
  mutate(w_rel = w / max(w)) %>%
  ungroup()

overlap_slope_choice <- ggplot(overlap_slope_imp, aes(overlap_slope, w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Overlap_slope",
    x = "overlap_slope form",
    y = "Summed AIC weight"
  )

overlap_slope_choice

overlap_imp <- mod_tbl %>%
  group_by(model_group, overlap_form) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  group_by(model_group) %>%
  mutate(w_rel = w / max(w)) %>%
  ungroup()

overlap_choice <- ggplot(overlap_imp, aes(overlap_form, w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Overlap form",
    x = "Overlap form",
    y = "Summed AIC weight"
  )

overlap_choice

#shark choice----
shark_imp <- mod_tbl %>%
  group_by(model_group, shark_col) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  mutate(w_rel = w / max(w)) %>%
  ungroup()

shark_choice <- ggplot(shark_imp, aes(shark_col, w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Shark index form",
    x = "shark index",
    y = "Summed AIC weight"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,       # Rotates the text 45 degrees
      hjust = 1,        # Align text right-edge with axis ticks
      vjust = 1         # Prevents vertical overlap
    ),
    plot.margin = margin(b = 20, r = 10, l = 10, t = 10) # Gives room so labels aren't clipped off
  )

shark_choice

#shark transformation choice----
shark_rollmean_imp <- mod_tbl %>%
  group_by(model_group, shark_transform) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  mutate(w_rel = w / max(w)) %>%
  ungroup()

shark_rollmean_choice <- ggplot(shark_rollmean_imp, aes(shark_transform, w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Shark transformation \n(scaling and smoothing)",
    x = "shark_transform",
    y = "Summed AIC weight"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,       # Rotates the text 45 degrees
      hjust = 1,        # Align text right-edge with axis ticks
      vjust = 1         # Prevents vertical overlap
    ),
    plot.margin = margin(b = 20, r = 10, l = 10, t = 10) # Gives room so labels aren't clipped off
  )

shark_rollmean_choice

#temperature choice----
library(tidyverse)

# 1. Calculate weights and reorder shark_temp_col by total weight
sst_imp <- mod_tbl %>%
  group_by(model_group, shark_temp_col) %>%
  summarise(w = sum(aic_weight), .groups = "drop") %>%
  # Sum weights across model_groups to establish overall rank for ordering
  group_by(shark_temp_col) %>%
  mutate(total_w = sum(w)) %>%
  ungroup() %>%
  # Reorder factor levels descending (.desc = TRUE for highest first)
  mutate(
    shark_temp_col = fct_reorder(shark_temp_col, total_w, .desc = TRUE),
    w_rel = w / max(w)
  )

# 2. Plot with ordered x-axis
sst_choice <- ggplot(sst_imp, aes(x = shark_temp_col, y = w, fill = model_group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "SST index",
    x = "sst index",
    y = "Summed AIC weight"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,       # Rotates the text 45 degrees
      hjust = 1,        # Align text right-edge with axis ticks
      vjust = 1         # Prevents vertical overlap
    ),
    plot.margin = margin(b = 20, r = 10, l = 10, t = 10) # Gives room so labels aren't clipped off
  )

sst_choice

#all plots---------
sst_choice
shark_choice
shark_rollmean_choice
overlap_slope_choice
overlap_choice
slope_and_shape_curves


#combined plot----------        
library(patchwork)
library(gridExtra)
# #make a table to fill in a plot slot
# myoutput<-round(c(best_par,r2=fit_metrics%>% filter(metric=="r2") %>% pull()),3)
# myoutput_df<-data.frame(
#   Parameter =names(myoutput),
#   Value     = sprintf("%.3f", myoutput),  stringsAsFactors = FALSE )
# 
# p_table <- wrap_elements(
#   tableGrob(
#     myoutput_df, 
#     rows = NULL, # removes row numbers
#     theme = ttheme_minimal(
#       core = list(fg_params = list(fontsize = 9)),
#       colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
#     )
#   )
# )

myplots <- (sst_choice +  shark_choice + shark_rollmean_choice + 
              overlap_slope_choice + overlap_choice + slope_and_shape_curves) +
  plot_layout(ncol = 3, guides = "collect") + 
  plot_annotation(title = "Shark index parameter selection by sem x16_sar") & 
  theme(
    plot.title = element_text(size = 9),      # Shrinks panel titles (e.g., "Linear predictor components...")
    plot.subtitle = element_text(size = 7),   # Shrinks panel subtitles if present
    axis.title = element_text(size = 8),      # Optional: shrinks axis labels
    axis.text = element_text(size = 7)        # Optional: shrinks axis tick text
  )

print(myplots)

ggsave(file.path(out_dir, "Shark_parameter_selection.png"), myplots, width = 8, height = 6, dpi = 170)



message("Saved weighted functional-relationship plots and tables to: ", out_dir)

best_mod<- mod_tbl %>% group_by(model_group) %>% slice(1:5);best_mod
t(best_mod)
