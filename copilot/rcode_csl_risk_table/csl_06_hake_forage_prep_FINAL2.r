# ==============================================================================
# Script: Ocean Forage Selection Pipeline with Robust Kable Accounting Table
# Threshold: |r| > 0.58 | Excludes Age-0 recruits & Age-1 Hake
# Routes Mackerel metrics to Sitka Herring Anchor
# ==============================================================================

library(tidyverse)
library(janitor)
library(knitr)
library(patchwork)

output_dir <- "copilot/outputs_csl_cr"

# 1. Safe Directory Check
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Data Ingestion, Merging & Explicit Non-Prey Exclusions
# -----------------------------------------------------------------------------

jake_sub <- read.csv(file.path(output_dir, "jake_lre_dat_yearly.csv")) %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_jake")))

doug_sub <- read.csv("data_Lisa/sem_master_data.csv") %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_doug")))

merged_prey <- jake_sub %>%
  inner_join(doug_sub, by = "year") %>%
  filter(year >= 1998, year <= 2024)

all_candidate_cols <- colnames(merged_prey)[
  grepl("hake|sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)
]

# Explicitly exclude non-prey juvenile life stages
non_prey_life_stages <- c("est_abund_herring_recruits_doug", "x03_hake_age1_doug")
target_cols          <- setdiff(all_candidate_cols, non_prey_life_stages)

# -----------------------------------------------------------------------------
# 2. Define Primary Anchors & Targeted Grouping (|r| > 0.58)
# -----------------------------------------------------------------------------

hake_anchor    <- "x09_dfa_hake_age5plus_doug"
sardine_anchor <- "x05_dfa_abund_sardine_doug"
sitka_anchor   <- "sitka_herring_e_go_a_doug"
herring_anchor <- "x05_herring_ncc_doug"
anchovy_anchor <- "x05_anchovy_gam_doug"

final_5 <- c(hake_anchor, sardine_anchor, sitka_anchor, herring_anchor, anchovy_anchor)

cor_mat <- cor(merged_prey %>% select(all_of(target_cols)), use = "pairwise.complete.obs")

get_correlated_series <- function(anchor, cols, mat, threshold = 0.58) {
  cors <- mat[anchor, cols]
  names(cors[abs(cors) >= threshold & !names(cors) %in% final_5])
}

# Identify all mackerel columns upfront so they group with Sitka Herring
mackerel_cols <- target_cols[grepl("mackerel", target_cols, ignore.case = TRUE)]
non_mackerel_eval <- setdiff(target_cols, c(final_5, mackerel_cols))

# Sequential Clustering: Mackerel forced to Sitka
hake_cluster    <- get_correlated_series(hake_anchor, non_mackerel_eval, cor_mat)
sardine_cluster <- get_correlated_series(sardine_anchor, setdiff(non_mackerel_eval, hake_cluster), cor_mat)

# Sitka receives all Mackerel metrics plus any other meeting the threshold
sitka_eval_pool <- c(mackerel_cols, setdiff(non_mackerel_eval, c(hake_cluster, sardine_cluster)))
sitka_cluster   <- get_correlated_series(sitka_anchor, sitka_eval_pool, cor_mat)

# Remaining pools
remaining_pool  <- setdiff(target_cols, c(final_5, hake_cluster, sardine_cluster, sitka_cluster))
herring_cluster <- get_correlated_series(herring_anchor, remaining_pool, cor_mat)
anchovy_cluster <- get_correlated_series(anchovy_anchor, setdiff(remaining_pool, herring_cluster), cor_mat)

# -----------------------------------------------------------------------------
# 3. Formatted Accounting Table using Base knitr::kable()
# -----------------------------------------------------------------------------

anchor_cor_matrix <- abs(cor_mat[target_cols, final_5])
max_r_values      <- apply(anchor_cor_matrix, 1, max, na.rm = TRUE)

accounting_df <- data.frame(
  Candidate_Metric  = target_cols,
  Assigned_Category = case_when(
    target_cols == hake_anchor ~ "FINAL 1: Hake Anchor",
    target_cols == sardine_anchor ~ "FINAL 2: Sardine Anchor",
    target_cols == sitka_anchor ~ "FINAL 3: Sitka Herring Anchor",
    target_cols == herring_anchor ~ "FINAL 4: NCC Herring Anchor",
    target_cols == anchovy_anchor ~ "FINAL 5: Executive Anchovy Choice",
    target_cols %in% hake_cluster ~ "Cluster 1: Correlated with Hake (|r| > 0.58)",
    target_cols %in% sardine_cluster ~ "Cluster 2: Correlated with Sardine (|r| > 0.58)",
    target_cols %in% sitka_cluster ~ "Cluster 3: Correlated with Sitka Herring (|r| > 0.58)",
    target_cols %in% herring_cluster ~ "Cluster 4: Correlated with NCC Herring (|r| > 0.58)",
    target_cols %in% anchovy_cluster ~ "Cluster 5: Correlated with NCC Anchovy (|r| > 0.58)",
    TRUE ~ "Uncorrelated / Low Correlation (< 0.58)"
  ),
  Max_R = round(max_r_values, 3)
) %>%
  arrange(Assigned_Category, desc(Max_R))

# Generate clean Kable table without kableExtra dependency errors
accounting_kable <- kable(
  accounting_df,
  caption = "Complete Candidate Metric Accounting (|r| > 0.58 Threshold)",
  col.names = c("Candidate Metric", "Assigned Cluster Category", "Max Absolute R"),
  align = c("l", "l", "c")
)

print(accounting_kable)

# -----------------------------------------------------------------------------
# 4. Prepare Plot Structure
# -----------------------------------------------------------------------------

long_prey <- merged_prey %>%
  select(year, all_of(target_cols)) %>%
  pivot_longer(-year, names_to = "metric_name", values_to = "raw_val") %>%
  group_by(metric_name) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max)
  ) %>%
  ungroup()

build_column_structure <- function(anchor_name, cluster_vars, group_label) {
  anchor_df <- data.frame(
    metric_name = anchor_name,
    col_group   = group_label,
    row_rank    = 1,
    is_top      = TRUE,
    panel_title = paste0("ANCHOR: ", anchor_name)
  )
  
  if (length(cluster_vars) > 0) {
    cluster_df <- data.frame(
      metric_name = cluster_vars,
      col_group   = group_label,
      row_rank    = seq_along(cluster_vars) + 1,
      is_top      = FALSE,
      panel_title = paste0("Sub ", seq_along(cluster_vars), ": ", cluster_vars)
    )
    return(bind_rows(anchor_df, cluster_df))
  } else {
    return(anchor_df)
  }
}

plot_structure <- bind_rows(
  build_column_structure(hake_anchor, hake_cluster, "1. Hake (Age 5+)"),
  build_column_structure(sardine_anchor, sardine_cluster, "2. Sardine Anchor"),
  build_column_structure(sitka_anchor, sitka_cluster, "3. Sitka Herring Anchor"),
  build_column_structure(herring_anchor, herring_cluster, "4. NCC Herring"),
  build_column_structure(anchovy_anchor, anchovy_cluster, "5. NCC Anchovy (Exec)")
) %>%
  distinct(metric_name, .keep_all = TRUE)

long_plotted <- long_prey %>%
  inner_join(plot_structure, by = "metric_name") %>%
  arrange(col_group, row_rank) %>%
  mutate(
    col_group   = factor(col_group),
    panel_title = factor(panel_title, levels = unique(plot_structure$panel_title))
  )

# -----------------------------------------------------------------------------
# 5. Plot 1: 5-Column Grid View
# -----------------------------------------------------------------------------

p_grid <- ggplot(long_plotted, aes(x = year, y = std_val, color = is_top)) +
  geom_line(aes(linewidth = is_top)) +
  geom_point(aes(size = is_top)) +
  scale_linewidth_manual(values = c("FALSE" = 0.6, "TRUE" = 1.2)) +
  scale_size_manual(values = c("FALSE" = 1.0, "TRUE" = 2.0)) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "firebrick")) +
  facet_wrap(~ col_group + panel_title, ncol = 5, dir = "v", scales = "free_y") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 5)) +
  labs(
    title = "Structured Feature Selection Grid: Final 5 Selections (Top Row) & Correlated Candidate Series",
    subtitle = "Top Row (Red) = Final 5 Selected Series; Lower Rows (Grey) = Candidate Series Correlated (|r| > 0.58)",
    x = "Year",
    y = "Standardized Value (0–1)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )

print(p_grid)

# -----------------------------------------------------------------------------
# 6. Plot 2: 5-Panel Overlay View (Bold Dark Red Anchor Over Thin Grey Correlated Series)
# -----------------------------------------------------------------------------

make_overlay_panel <- function(group_id) {
  df_sub <- long_plotted %>% filter(col_group == group_id)
  
  bg_series     <- df_sub %>% filter(!is_top)
  anchor_series <- df_sub %>% filter(is_top)
  
  ggplot() +
    geom_line(
      data = bg_series,
      aes(x = year, y = std_val, group = metric_name),
      color = "grey60", linewidth = 0.7, alpha = 0.7
    ) +
    geom_line(
      data = anchor_series,
      aes(x = year, y = std_val),
      color = "firebrick", linewidth = 1.5
    ) +
    geom_point(
      data = anchor_series,
      aes(x = year, y = std_val),
      color = "firebrick", size = 2.2
    ) +
    facet_wrap(~ col_group) +
    scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
    labs(x = NULL, y = "Std Val (0–1)") +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
}

groups_list <- levels(long_plotted$col_group)
panel_plots <- lapply(groups_list, make_overlay_panel)

p_5panel_overlay <- wrap_plots(panel_plots, ncol = 1) +
  plot_annotation(
    title = "5-Panel Ocean Prey Overlay Comparison",
    subtitle = "Dark Red Line = Primary Selected Anchor; Thin Grey Lines = Candidate Series Correlated (|r| > 0.58)"
  )

print(p_5panel_overlay)

# -----------------------------------------------------------------------------
# 7. Save Final Outputs
# -----------------------------------------------------------------------------

final_prey_df <- merged_prey %>% select(year, all_of(final_5))

write.csv(final_prey_df, file.path(output_dir, "csl_ocean_prey_5series_1998_2024.csv"), row.names = FALSE)
write.csv(accounting_df, file.path(output_dir, "prey_accounting_with_max_r.csv"), row.names = FALSE)

ggsave(file.path(output_dir, "forage_feature_selection_grid_0.58.png"), p_grid, width = 16, height = 11)
ggsave(file.path(output_dir, "forage_5panel_overlay_0.58.png"), p_5panel_overlay, width = 12, height = 14)

cat("\nPipeline Complete: Saved accounting table and both multi-panel plots to", output_dir, "\n")