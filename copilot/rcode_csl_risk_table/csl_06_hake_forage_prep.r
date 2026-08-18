
#RESULTS:
# Table: Complete Accounting of All Evaluated Candidate Prey Time Series
# 
# |Candidate_Metric              |Assigned_Category                                         |
# |:-----------------------------|:---------------------------------------------------------|
# |sardine_northern_bmass_jake   |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |herring_bmass_jake            |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |abund_herring_doug            |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |abund_sardine_doug            |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |hake_age5plus_2025_doug       |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |jackmackerel_gam_2025_doug    |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |pacificmackerel_gam_2025_doug |Cluster 1: Correlated with Hake (&#124;r&#124; > 0.69)    |
# |anchovy_northern_bmass_jake   |Cluster 2: Correlated with Sardine (&#124;r&#124; > 0.69) |
# |x05_herring_gam_doug          |Cluster 2: Correlated with Sardine (&#124;r&#124; > 0.69) |
# |sardine_ncc_doug              |Cluster 2: Correlated with Sardine (&#124;r&#124; > 0.69) |
# |sardine_gam_2025_doug         |Cluster 2: Correlated with Sardine (&#124;r&#124; > 0.69) |
# |jack_mackerel_bmass_jake      |Cluster 3: Correlated with Sitka / Jack Mackerel          |
# |hake_cpue_jake                |Cluster 5: Executive Selection Group                      |
# |x03_hake_age1_doug            |Cluster 5: Executive Selection Group                      |
# |hake_fishery_wa_2025_doug     |Cluster 5: Executive Selection Group                      |
# |x09_dfa_hake_age5plus_doug    |FINAL 1: Hake Anchor                                      |
# |x05_dfa_abund_sardine_doug    |FINAL 2: Sardine Anchor                                   |
# |sitka_herring_e_go_a_doug     |FINAL 3: Sitka Herring Stand-in                           |
# |x05_herring_ncc_doug          |FINAL 4: NCC Herring                                      |
# |x05_anchovy_gam_doug          |FINAL 5: Executive Anchovy Choice                         |




---
  title: "NCC Alternative Ocean Prey Feature Selection Analysis"
subtitle: "Evaluating Hake, Sardine, Herring, Anchovy, and Mackerel Candidate Series"
author: "Columbia River & NCC Ecology Working Group"
date: "`r Sys.Date()`"
output:
  html_document:
  toc: true
toc_float: true
theme: flatly
highlight: tango
code_folding: show
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk_set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 16, fig.height = 11)
library(tidyverse)
library(janitor)
library(knitr)

output_dir <- "copilot/outputs_csl_cr"
dir.create(output_dir, show_warnings = FALSE, recursive = TRUE)

#1. Executive Summary & Selection Rationale-----
#To incorporate alternative ocean prey into the California Sea Lion (CSL) predation risk framework, we evaluated candidate time series for Hake and key forage fish (Sardine, Anchovy, Herring, Mackerel) from sem_master_data.csv ("Doug") and jake_lre_dat_yearly.csv ("Jake").Selection Criteria & Decisions:Explicit Exclusion of Age-0 Recruits: est_abund_herring_recruits_doug was explicitly removed because Age-0 recruits are too small to serve as adult CSL forage.Missing Data Handling: sitka_herring_e_go_a_doug was chosen as a stand-in for jack_mackerel_bmass_jake ($r > 0.80$) due to multi-year data gaps in the mackerel series.Correlation Thresholding: Candidate series with high correlation ($|r| > 0.69$) were grouped under their primary anchor metric.Regional & Ecological Retentions: x05_herring_ncc_doug was retained as an orthogonal regional signal, and x05_anchovy_gam_doug was retained by executive decision due to its known role as a primary CSL prey driver in the NCC estuary.

# 1. Load and clean "Jake" dataset
jake_sub <- read.csv(file.path(output_dir, "jake_lre_dat_yearly.csv")) %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_jake")))

# 2. Load and clean "Doug" dataset
doug_sub <- read.csv("data_Lisa/sem_master_data.csv") %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_doug")))

# 3. Inner join on Year (1998–2024)
merged_prey <- jake_sub %>%
  inner_join(doug_sub, by = "year") %>%
  filter(year >= 1998, year <= 2024)

# 4. Identify all candidate columns and explicitly exclude age-0 recruits
all_candidate_cols <- colnames(merged_prey)[
  grepl("hake|sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)
]

target_cols <- setdiff(all_candidate_cols, "est_abund_herring_recruits_doug")

cat("Total candidate forage/hake time series evaluated (excluding age-0 recruits):", length(target_cols), "\n")


#2. Data Ingestion & Merging--------
# 1. Load and clean "Jake" dataset
jake_sub <- read.csv(file.path(output_dir, "jake_lre_dat_yearly.csv")) %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_jake")))

# 2. Load and clean "Doug" dataset
doug_sub <- read.csv("data_Lisa/sem_master_data.csv") %>%
  clean_names() %>%
  select(matches("year|hake|sardine|herring|anchovy|mackerel")) %>%
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_doug")))

# 3. Inner join on Year (1998–2024)
merged_prey <- jake_sub %>%
  inner_join(doug_sub, by = "year") %>%
  filter(year >= 1998, year <= 2024)

# 4. Identify all candidate columns and explicitly exclude age-0 recruits
all_candidate_cols <- colnames(merged_prey)[
  grepl("hake|sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)
]

target_cols <- setdiff(all_candidate_cols, "est_abund_herring_recruits_doug")

cat("Total candidate forage/hake time series evaluated (excluding age-0 recruits):", length(target_cols), "\n")

#3. Correlation Audit & Accounting----------
#Every candidate metric is assigned to one of our 5 Final Selected Series based on correlation strength ($|r| > 0.69$) or ecological logic.

# Define the 5 Final Selected Series
hake_anchor    <- "x09_dfa_hake_age5plus_doug"
sardine_anchor <- "x05_dfa_abund_sardine_doug"
sitka_anchor   <- "sitka_herring_e_go_a_doug"
herring_anchor <- "x05_herring_ncc_doug"
anchovy_anchor <- "x05_anchovy_gam_doug"

final_5 <- c(hake_anchor, sardine_anchor, sitka_anchor, herring_anchor, anchovy_anchor)

# Compute complete pairwise Pearson correlation matrix
cor_mat <- cor(merged_prey %>% select(all_of(target_cols)), use = "pairwise.complete.obs")

# Helper function to extract series meeting threshold (|r| > 0.69)
get_correlated_series <- function(anchor, cols, mat, threshold = 0.69) {
  cors <- mat[anchor, cols]
  names(cors[abs(cors) >= threshold & names(cors) != anchor])
}

# Group candidate metrics under top row anchors
hake_cluster    <- get_correlated_series(hake_anchor, target_cols, cor_mat)
sardine_cluster <- get_correlated_series(sardine_anchor, setdiff(target_cols, hake_cluster), cor_mat)
sitka_cluster   <- get_correlated_series(sitka_anchor, setdiff(target_cols, c(hake_cluster, sardine_cluster)), cor_mat)
herring_cluster <- get_correlated_series(herring_anchor, setdiff(target_cols, c(hake_cluster, sardine_cluster, sitka_cluster)), cor_mat)
anchovy_cluster <- setdiff(target_cols, c(final_5, hake_cluster, sardine_cluster, sitka_cluster, herring_cluster))

# Verify that 100% of candidate series are accounted for
accounted_for <- c(final_5, hake_cluster, sardine_cluster, sitka_cluster, herring_cluster, anchovy_cluster)
cat("Audit Complete: All", length(target_cols), "candidate series accounted for?", length(setdiff(target_cols, accounted_for)) == 0, "\n")

#4. Correlation Threshold Table (|r| > 0.69)-----
#This table lists every candidate metric alongside its exact Pearson correlation ($r$) with its assigned primary anchor.

build_cor_df <- function(anchor_name, cluster_vars, mat) {
  if (length(cluster_vars) == 0) return(data.frame())
  data.frame(
    Anchor_Metric = anchor_name,
    Candidate_Metric = cluster_vars,
    Pearson_R = round(mat[anchor_name, cluster_vars], 3),
    Abs_R = round(abs(mat[anchor_name, cluster_vars]), 3)
  )
}

high_cor_table <- bind_rows(
  build_cor_df(hake_anchor, hake_cluster, cor_mat),
  build_cor_df(sardine_anchor, sardine_cluster, cor_mat),
  build_cor_df(sitka_anchor, sitka_cluster, cor_mat),
  build_cor_df(herring_anchor, herring_cluster, cor_mat),
  build_cor_df(anchovy_anchor, anchovy_cluster, cor_mat)
) %>%
  arrange(Anchor_Metric, desc(Abs_R))

kable(high_cor_table, caption = "Candidate Metrics Exceeding Correlation Threshold (|r| > 0.69) with Primary Anchors")


#5. Complete Accounting Table-------
accounting_df <- data.frame(
  Candidate_Metric = target_cols,
  Assigned_Category = case_when(
    target_cols == hake_anchor ~ "FINAL 1: Hake Anchor",
    target_cols == sardine_anchor ~ "FINAL 2: Sardine Anchor",
    target_cols == sitka_anchor ~ "FINAL 3: Sitka Herring Stand-in",
    target_cols == herring_anchor ~ "FINAL 4: NCC Herring",
    target_cols == anchovy_anchor ~ "FINAL 5: Executive Anchovy Choice",
    target_cols %in% hake_cluster ~ "Cluster 1: Correlated with Hake (|r| > 0.69)",
    target_cols %in% sardine_cluster ~ "Cluster 2: Correlated with Sardine (|r| > 0.69)",
    target_cols %in% sitka_cluster ~ "Cluster 3: Correlated with Sitka / Jack Mackerel",
    target_cols %in% herring_cluster ~ "Cluster 4: Correlated with NCC Herring",
    TRUE ~ "Cluster 5: Executive Selection Group"
  )
) %>%
  arrange(Assigned_Category)

kable(accounting_df, caption = "Complete Accounting of All Evaluated Candidate Prey Time Series")

#6. Final 5 Structured Feature Selection Grid Plot--------
#The grid below displays our Final 5 Selected Time Series along the top row (highlighted in red). The panels underneath show all candidate series that met the correlation threshold ($|r| > 0.69$).


```{r grid_plot, fig.width=16, fig.height=12}
# Normalize series 0–1 for visual comparison
long_prey <- merged_prey %>%
  select(year, all_of(target_cols)) %>%
  pivot_longer(-year, names_to = "metric_name", values_to = "raw_val") %>%
  group_by(metric_name) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max)
  ) %>%
  ungroup()

# Safe helper to construct cluster rows without differing row-length errors
build_cluster_rows <- function(anchor_name, cluster_vars, group_label) {
  # Anchor row (Top Row = 1)
  anchor_row <- data.frame(
    metric_name = anchor_name,
    col_group   = group_label,
    row_rank    = 1,
    is_top      = TRUE
  )
  
  # Sub-cluster rows (only if cluster_vars is not empty)
  if (length(cluster_vars) > 0) {
    cluster_rows <- data.frame(
      metric_name = cluster_vars,
      col_group   = group_label,
      row_rank    = seq_along(cluster_vars) + 1,
      is_top      = FALSE
    )
    return(bind_rows(anchor_row, cluster_rows))
  } else {
    return(anchor_row)
  }
}

# Map metrics safely into 5 grid columns
plot_structure <- bind_rows(
  build_cluster_rows(hake_anchor, hake_cluster, "1. Hake (Age 5+)"),
  build_cluster_rows(sardine_anchor, sardine_cluster, "2. Sardine Anchor"),
  build_cluster_rows(sitka_anchor, sitka_cluster, "3. Sitka Herring (Stand-in)"),
  build_cluster_rows(herring_anchor, herring_cluster, "4. NCC Herring"),
  build_cluster_rows(anchovy_anchor, anchovy_cluster, "5. NCC Anchovy (Exec)")
)

long_plotted <- long_prey %>%
  inner_join(plot_structure, by = "metric_name", relationship = "many-to-many")

p_grid <- ggplot(long_plotted, aes(x = year, y = std_val, color = is_top)) +
  geom_line(aes(linewidth = is_top)) +
  geom_point(aes(size = is_top)) +
  scale_linewidth_manual(values = c("FALSE" = 0.6, "TRUE" = 1.2)) +
  scale_size_manual(values = c("FALSE" = 1.0, "TRUE" = 2.0)) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "firebrick")) +
  facet_grid(row_rank ~ col_group, scales = "free_y") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 5)) +
  labs(
    title = "Structured Feature Selection Grid: Final 5 Anchors (Top Row) & Correlated Candidate Series",
    subtitle = "Top Row (Red) = Final 5 Selected Series; Lower Rows (Grey) = Candidate Series Correlated (|r| > 0.69)",
    x = "Year",
    y = "Standardized Value (0–1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

print(p_grid)

#updated p-grid-------
```{r grid_plot, fig.width=16, fig.height=12}
# 1. Normalize series 0–1 for visual comparison
long_prey <- merged_prey %>%
  select(year, all_of(target_cols)) %>%
  pivot_longer(-year, names_to = "metric_name", values_to = "raw_val") %>%
  group_by(metric_name) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max)
  ) %>%
  ungroup()

# 2. Build cluster structure with clean inline titles
build_cluster_rows <- function(anchor_name, cluster_vars, group_label) {
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
  build_cluster_rows(hake_anchor, hake_cluster, "1. Hake (Age 5+)"),
  build_cluster_rows(sardine_anchor, sardine_cluster, "2. Sardine Anchor"),
  build_cluster_rows(sitka_anchor, sitka_cluster, "3. Sitka Herring (Stand-in)"),
  build_cluster_rows(herring_anchor, herring_cluster, "4. NCC Herring"),
  build_cluster_rows(anchovy_anchor, anchovy_cluster, "5. NCC Anchovy (Exec)")
) %>%
  distinct(metric_name, .keep_all = TRUE) # Prevents many-to-many join warnings

# 3. Join and order factors so Top Row stays at the top of each column
long_plotted <- long_prey %>%
  inner_join(plot_structure, by = "metric_name") %>%
  arrange(col_group, row_rank) %>%
  mutate(
    col_group   = factor(col_group),
    panel_title = factor(panel_title, levels = unique(plot_structure$panel_title))
  )

# 4. Generate Plot with Direct Header Labels Above Every Single Panel
p_grid <- ggplot(long_plotted, aes(x = year, y = std_val, color = is_top)) +
  geom_line(aes(linewidth = is_top)) +
  geom_point(aes(size = is_top)) +
  scale_linewidth_manual(values = c("FALSE" = 0.6, "TRUE" = 1.2)) +
  scale_size_manual(values = c("FALSE" = 1.0, "TRUE" = 2.0)) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "firebrick")) +
  # Use facet_wrap with dir = "v" so panels fill down Column 1, then Column 2, etc.
  facet_wrap(~ col_group + panel_title, ncol = 5, dir = "v", scales = "free_y") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 5)) +
  labs(
    title = "Structured Feature Selection Grid: Final 5 Anchors (Top Row) & Correlated Candidate Series",
    subtitle = "Red = Final 5 Selected Series; Grey = Correlated Candidate Series (|r| > 0.69)",
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


#7. Dataset Export-------
# Export the clean 5-series ocean prey dataset for CSL Risk Table merging
final_prey_df <- merged_prey %>% select(year, all_of(final_5))

write.csv(final_prey_df, file.path(output_dir, "csl_ocean_prey_5series_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "forage_feature_selection_grid.png"), p_grid, width = 16, height = 11)

cat("Export Complete: csl_ocean_prey_5series_1998_2024.csv successfully created.\n")



