
#Results: selected alternative temperature variables for sharks

# Suggested shark-temp candidates for next sweep
# Based on your correlations, choose a diverse set rather than top-correlated only:
#   
# sst_wgoa_coastwatch_junjulaug (reference)
# bts_wgoa_temp_sum_195to205m (deep, moderate corr)
# swln_temp_spr_176to226m (deep spring)
# ll_wgoa_temp_sum_246to255m (deepest, weak/negative corr)
# pdo_djf
# enso_dj

#write_csv(corr_tbl, file.path(out_dir, "akyr_temp_pdo_enso_cor_with_selected_sst.csv"))

#> corr_tbl
# A tibble: 14 × 3
# variable                    cor_with_selected_sst abs_cor
# 1 bts_wgoa_temp_sum_1to5m                     0.825   0.825
# 2 secm_temp_sum_10m                           0.712   0.712
# 3 adfg_lmesh_temp_summ_bttm                   0.688   0.688
# 4 bts_egoa_temp_sum_1to5m                     0.607   0.607
# 5 swln_temp_spr_0to10m                        0.592   0.592
# 6 swln_temp_spr_176to226m                     0.515   0.515
# 7 bts_wgoa_temp_sum_195to205m                 0.514   0.514
# 8 swln_temp_fall_176to226m                    0.432   0.432
# 9 pdo_djf                                     0.387   0.387
# 10 enso_dj                                     0.338   0.338
# 11 bts_egoa_temp_sum_195to205m                 0.289   0.289
# 12 ll_egoa_temp_sum_246to255m                  0.279   0.279
# 13 swln_temp_fall_0to10m                       0.206   0.206
# 14 ll_wgoa_temp_sum_246to255m                 -0.205   0.205


suppressPackageStartupMessages({
  library(tidyverse)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("ak_yr")) stop("ak_yr not found in environment. Run annual-table build first.")

selected_sst <- "sst_wgoa_coastwatch_junjulaug"
if (!selected_sst %in% names(ak_yr)) {
  stop(paste0("Selected SST column not found in ak_yr: ", selected_sst))
}

# find candidate climate columns
temp_cols <- names(ak_yr)[stringr::str_detect(names(ak_yr), regex("temp", ignore_case = TRUE))]
pdo_cols  <- names(ak_yr)[stringr::str_detect(names(ak_yr), regex("pdo",  ignore_case = TRUE))]
enso_cols <- names(ak_yr)[stringr::str_detect(names(ak_yr), regex("enso", ignore_case = TRUE))]

# include selected SST explicitly even if it doesn't include "temp"
clim_cols <- unique(c(selected_sst, temp_cols, pdo_cols, enso_cols))
clim_cols <- clim_cols[clim_cols %in% names(ak_yr)]

if (length(clim_cols) == 0) stop("No matching temp/pdo/enso columns found in ak_yr.")

# keep only numeric columns
clim_cols <- clim_cols[sapply(ak_yr[clim_cols], is.numeric)]
if (length(clim_cols) == 0) stop("Matching columns exist but none are numeric.")

# long table
plot_dat <- ak_yr %>%
  select(year, all_of(clim_cols)) %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  mutate(group = case_when(
    variable == selected_sst ~ "selected_sst",
    str_detect(variable, regex("pdo", ignore_case = TRUE)) ~ "pdo",
    str_detect(variable, regex("enso", ignore_case = TRUE)) ~ "enso",
    str_detect(variable, regex("temp", ignore_case = TRUE)) ~ "temp",
    TRUE ~ "other"
  ))

# z-scored version for comparability
plot_dat_z <- plot_dat %>%
  group_by(variable) %>%
  mutate(value_z = {
    s <- sd(value, na.rm = TRUE)
    ifelse(is.na(s) | s == 0, 0, as.numeric(scale(value)))
  }) %>%
  ungroup()

# 1) Faceted raw values
p_raw <- ggplot(plot_dat, aes(year, value, color = variable)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "AK annual climate candidate series (raw scale)",
    subtitle = paste0("Includes selected SST: ", selected_sst),
    x = "Year", y = "Raw value"
  )

# 2) Faceted standardized values
p_z <- ggplot(plot_dat_z, aes(year, value_z, color = variable, linewidth = (variable == selected_sst))) +
  geom_line(alpha = 0.9, show.legend = FALSE) +
  scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.6)) +
  facet_wrap(~group, scales = "fixed", ncol = 1) +
  theme_bw() +
  labs(
    title = "AK annual climate candidate series (standardized)",
    subtitle = "Z-score by series; selected SST drawn thicker",
    x = "Year", y = "Z-score"
  )

# 3) Overlay of standardized series (all together)
p_overlay <- ggplot(plot_dat_z, aes(year, value_z, group = variable)) +
  geom_line(color = "grey70", linewidth = 0.45, alpha = 0.8) +
  geom_line(
    data = plot_dat_z %>% filter(variable == selected_sst),
    aes(year, value_z),
    color = "red3", linewidth = 1.2
  ) +
  theme_bw() +
  labs(
    title = "All standardized climate candidates overlaid",
    subtitle = paste0("Red = selected SST (", selected_sst, ")"),
    x = "Year", y = "Z-score"
  )

# correlation table against selected SST (full overlap)
corr_tbl <- ak_yr %>%
  select(year, all_of(clim_cols)) %>%
  summarise(across(
    all_of(setdiff(clim_cols, selected_sst)),
    ~ cor(.x, .data[[selected_sst]], use = "pairwise.complete.obs")
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "cor_with_selected_sst") %>%
  mutate(abs_cor = abs(cor_with_selected_sst)) %>%
  arrange(desc(abs_cor))

write_csv(corr_tbl, file.path(out_dir, "akyr_temp_pdo_enso_cor_with_selected_sst.csv"))
write_csv(
  tibble(variable = clim_cols) %>%
    mutate(group = case_when(
      variable == selected_sst ~ "selected_sst",
      str_detect(variable, regex("pdo", ignore_case = TRUE)) ~ "pdo",
      str_detect(variable, regex("enso", ignore_case = TRUE)) ~ "enso",
      str_detect(variable, regex("temp", ignore_case = TRUE)) ~ "temp",
      TRUE ~ "other"
    )),
  file.path(out_dir, "akyr_temp_pdo_enso_variables_used.csv")
)

ggsave(file.path(out_dir, "akyr_temp_pdo_enso_raw_facet.png"), p_raw, width = 11, height = 10, dpi = 150)
ggsave(file.path(out_dir, "akyr_temp_pdo_enso_zscore_facet.png"), p_z, width = 11, height = 10, dpi = 150)
ggsave(file.path(out_dir, "akyr_temp_pdo_enso_zscore_overlay.png"), p_overlay, width = 11, height = 6, dpi = 150)


print(p_raw)
print(p_z)
print(p_overlay)

#2 plots didn't work, probably missing data

message("Saved:\n",
        "- akyr_temp_pdo_enso_raw_facet.png\n",
        "- akyr_temp_pdo_enso_zscore_facet.png\n",
        "- akyr_temp_pdo_enso_zscore_overlay.png\n",
        "- akyr_temp_pdo_enso_cor_with_selected_sst.csv\n",
        "- akyr_temp_pdo_enso_variables_used.csv")
       