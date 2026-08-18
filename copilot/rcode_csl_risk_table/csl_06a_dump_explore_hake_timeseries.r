# ==============================================================================
# Script: Compare Hake & Forage Fish Metrics ("Jake" vs "Doug" Datasets)
# ==============================================================================

library(tidyverse)
library(janitor)
library(corrplot)

output_dir <- "copilot/outputs_csl_cr"

# Target column pattern
target_pattern <- "year|hake|sardine|herring|anchovy|mackerel"

# -----------------------------------------------------------------------------
# 1. Load and Filter Both Datasets
# -----------------------------------------------------------------------------

# Dataset 1: "Jake" (jake_lre_dat_yearly.csv)
jake_raw <- read.csv(file.path(output_dir, "jake_lre_dat_yearly.csv")) %>%
  clean_names()

jake_sub <- jake_raw %>%
  select(matches(target_pattern)) %>%
  # Append "_jake" suffix to all metric columns except year
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_jake")))

# Dataset 2: "Doug" (sem_master_data.csv)
doug_raw <- read.csv("data_Lisa/sem_master_data.csv") %>%
  clean_names()

doug_sub <- doug_raw %>%
  select(matches(target_pattern)) %>%
  # Append "_doug" suffix to all metric columns except year
  rename_with(~ ifelse(.x == "year", "year", paste0(.x, "_doug")))

# -----------------------------------------------------------------------------
# 2. Merge Datasets Across Matching Years (1998–2024)
# -----------------------------------------------------------------------------

merged_prey <- jake_sub %>%
  inner_join(doug_sub, by = "year") %>%
  filter(year >= 1998, year <= 2024) %>%
  arrange(year)

cat("--- Column Comparison in Merged Dataset ---\n")
print(colnames(merged_prey))

# -----------------------------------------------------------------------------
# 3. Correlation Matrix Analysis ("Jake" vs "Doug")
# -----------------------------------------------------------------------------

# Calculate Pearson correlation matrix for all numeric metric columns
cor_data <- merged_prey %>%
  select(-year) %>%
  # Drop columns that are completely NA or zero variance
  select(where(~ sd(.x, na.rm = TRUE) > 0))

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

cat("\n--- Correlation Matrix Summary ---\n")
print(round(cor_matrix, 2))

# Optional Visual Correlation Heatmap
png(file.path(output_dir, "jake_vs_doug_forage_correlations.png"), width = 800, height = 800)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         title = "Correlation: Jake vs Doug Forage & Hake Metrics", mar = c(0,0,2,0))
dev.off()

# -----------------------------------------------------------------------------
# 4. Standardized Comparative Time-Series Plots
# -----------------------------------------------------------------------------

# Pivot to long format for standardized comparison plotting
long_prey <- merged_prey %>%
  pivot_longer(
    cols = -year,
    names_to = "metric_source",
    values_to = "raw_value"
  ) %>%
  separate(metric_source, into = c("species_metric", "source"), sep = "_(?=[^_]+$)") %>%
  group_by(species_metric) %>%
  mutate(
    # Standardize 0 to 1 for shape/timing comparison
    val_max = max(raw_value, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_value / val_max)
  ) %>%
  ungroup()

# Plot Comparative Overlay Curves
p_prey_comp <- ggplot(long_prey, aes(x = year, y = std_val, color = source, linetype = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~ species_metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("jake" = "steelblue", "doug" = "firebrick")) +
  labs(
    title = "Comparison of Hake & Forage Fish Time Series (Jake vs Doug)",
    subtitle = "Standardized 0–1 annual time series (1998–2024)",
    x = "Year",
    y = "Standardized Index (0 to 1)",
    color = "Data Source",
    linetype = "Data Source"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_prey_comp)

# Save Outputs
write.csv(merged_prey, file.path(output_dir, "jake_vs_doug_prey_merged.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "jake_vs_doug_time_series_comparison.png"), p_prey_comp, width = 12, height = 8)


#hake focus----------------------------
# ==============================================================================
# Script: Target Analysis — Hake Indices & Top Correlated Forage Species
# ==============================================================================

library(tidyverse)
library(janitor)

output_dir <- "copilot/outputs_csl_cr"

# Load merged data created in previous step (or re-merge)
merged_prey <- read.csv(file.path(output_dir, "jake_vs_doug_prey_merged.csv"))

# -----------------------------------------------------------------------------
# 1. Isolate Hake Columns & Forage Columns--------
# -----------------------------------------------------------------------------

hake_cols   <- colnames(merged_prey)[grepl("hake", colnames(merged_prey), ignore.case = TRUE)]
forage_cols <- colnames(merged_prey)[grepl("sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)]

cat("--- Identified Hake Metrics ---\n")
print(hake_cols)

cat("\n--- Identified Forage Fish Metrics ---\n")
print(forage_cols)

# -----------------------------------------------------------------------------
# 2. Hake-to-Hake Internal Consistency Matrix
# -----------------------------------------------------------------------------

hake_data <- merged_prey %>% select(all_of(hake_cols))
hake_cor_matrix <- cor(hake_data, use = "pairwise.complete.obs")

cat("\n========================================================\n")
cat("   HAKE INDEX INTERNAL CORRELATIONS (JAKE vs DOUG)      \n")
cat("========================================================\n")
print(round(hake_cor_matrix, 2))

# -----------------------------------------------------------------------------
# 3. Pairwise Correlations: Every Forage Metric vs Every Hake Metric
# -----------------------------------------------------------------------------

# Calculate pairwise correlation between all Hake cols and Forage cols
hake_forage_cors <- expand_grid(
  hake_metric   = hake_cols,
  forage_metric = forage_cols
) %>%
  rowwise() %>%
  mutate(
    pearson_r = cor(merged_prey[[hake_metric]], merged_prey[[forage_metric]], use = "complete.obs"),
    abs_r     = abs(pearson_r)
  ) %>%
  ungroup() %>%
  arrange(desc(abs_r))

cat("\n========================================================\n")
cat("   TOP FORAGE FISH METRICS CORRELATED WITH HAKE         \n")
cat("========================================================\n")
print(head(hake_forage_cors, 25))
print(hake_forage_cors, n=25)

# -----------------------------------------------------------------------------
# 4. Clean Focused Plot: Hake Indices Over Time
# -----------------------------------------------------------------------------

hake_long <- merged_prey %>%
  select(year, all_of(hake_cols)) %>%
  pivot_longer(-year, names_to = "hake_index", values_to = "raw_val") %>%
  group_by(hake_index) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max)
  ) %>%
  ungroup()

p_hake_only <- ggplot(hake_long, aes(x = year, y = std_val, color = hake_index)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Hake Index Comparison Across Time (1998–2024)",
    subtitle = "Standardized 0–1 trajectories for all available Hake metrics",
    x = "Year",
    y = "Standardized Value (0 to 1)",
    color = "Hake Metric"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_hake_only)

# -----------------------------------------------------------------------------
# 5. Clean Focused Plot: Top 4 Hake-Correlated Forage Fish
# -----------------------------------------------------------------------------

# Select top forage metrics based on absolute correlation
top_forage_names <- head(unique(hake_forage_cors$forage_metric), 4)
top_hake_name    <- hake_forage_cors$hake_metric[1]

top_compare_long <- merged_prey %>%
  select(year, all_of(c(top_hake_name, top_forage_names))) %>%
  pivot_longer(-year, names_to = "series_name", values_to = "raw_val") %>%
  group_by(series_name) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max),
    is_hake = grepl("hake", series_name, ignore.case = TRUE)
  ) %>%
  ungroup()

p_top_correlated <- ggplot(top_compare_long, aes(x = year, y = std_val, color = series_name, linetype = is_hake)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~ series_name, ncol = 2) +
  labs(
    title = "Top Forage Fish Metrics Correlated with Primary Hake Index",
    subtitle = "Standardized 0–1 trajectories comparing Hake with key forage species",
    x = "Year",
    y = "Standardized Value (0 to 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_top_correlated)

# -----------------------------------------------------------------------------
# 6. Save Outputs
# -----------------------------------------------------------------------------

write.csv(hake_forage_cors, file.path(output_dir, "hake_forage_correlation_rankings.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "hake_indices_comparison.png"), p_hake_only, width = 10, height = 6)
ggsave(file.path(output_dir, "hake_vs_top_forage_series.png"), p_top_correlated, width = 10, height = 7)


#hake, sadine and herring time series selected.--------
#now look at other forage fish-------
# ==============================================================================
# Script: Forage Fish Inter-Correlation Analysis & Anchor Comparison
# Anchor Hake Index:    x09_dfa_hake_age5plus_doug
# Anchor Sardine Index: x05_dfa_abund_sardine_doug
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# Load merged prey dataset
merged_prey <- read.csv(file.path(output_dir, "jake_vs_doug_prey_merged.csv"))

# Define key anchors
hake_anchor    <- "x09_dfa_hake_age5plus_doug"
sardine_anchor <- "x05_dfa_abund_sardine_doug"

# -----------------------------------------------------------------------------
# 1. Select Anchor Metrics & Remaining Forage Species Columns
# -----------------------------------------------------------------------------

# Find all remaining forage fish columns (anchovy, herring, mackerel, sardine)
all_forage_cols <- colnames(merged_prey)[
  grepl("sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)
]

# Create targeted dataframe containing the anchors + all forage metrics
target_cols <- unique(c("year", hake_anchor, sardine_anchor, all_forage_cols))

# Subset using a single all_of() wrapper
forage_focus <- merged_prey %>%
  select(all_of(target_cols))
cat("--- Selected Columns for Forage Inter-Correlation Analysis ---\n")
print(colnames(forage_focus))

# -----------------------------------------------------------------------------
# 2. Pairwise Correlation Analysis Among Forage Metrics & Anchors
# -----------------------------------------------------------------------------

forage_data_only <- forage_focus %>% select(-year)

# Compute complete pairwise Pearson correlation matrix
cor_matrix_forage <- cor(forage_data_only, use = "pairwise.complete.obs")

# Format into a clean, rankable table of correlations relative to the anchors
forage_cor_table <- expand_grid(
  Metric_A = colnames(forage_data_only),
  Metric_B = colnames(forage_data_only)
) %>%
  filter(Metric_A < Metric_B) %>% # Remove duplicate pairs & self-correlations
  rowwise() %>%
  mutate(
    pearson_r = cor(forage_focus[[Metric_A]], forage_focus[[Metric_B]], use = "complete.obs"),
    abs_r     = abs(pearson_r)
  ) %>%
  ungroup() %>%
  arrange(desc(abs_r))

cat("\n========================================================\n")
cat("   TOP PAIRWISE CORRELATIONS AMONG FORAGE SPECIES       \n")
cat("========================================================\n")
print(head(forage_cor_table, 20))

# Specifically highlight correlations against Anchor Sardine & Anchor Hake
anchor_cors <- forage_cor_table %>%
  filter(Metric_A %in% c(hake_anchor, sardine_anchor) | Metric_B %in% c(hake_anchor, sardine_anchor)) %>%
  arrange(desc(abs_r))

cat("\n========================================================\n")
cat("   CORRELATIONS AGAINST ANCHOR SARDINE & ANCHOR HAKE    \n")
cat("========================================================\n")
print(anchor_cors)

# -----------------------------------------------------------------------------
# 3. Standardized Time-Series Plot: Forage Fish vs Anchor Sardine & Hake
# -----------------------------------------------------------------------------

long_forage_focus <- forage_focus %>%
  pivot_longer(-year, names_to = "series_name", values_to = "raw_val") %>%
  group_by(series_name) %>%
  mutate(
    val_max = max(raw_val, na.rm = TRUE),
    std_val = if_else(val_max == 0 | is.na(val_max), 0, raw_val / val_max),
    is_anchor = series_name %in% c(hake_anchor, sardine_anchor)
  ) %>%
  ungroup()

p_forage_comparison <- ggplot(long_forage_focus, aes(x = year, y = std_val, color = series_name)) +
  geom_line(aes(linewidth = is_anchor, linetype = is_anchor)) +
  scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.3)) +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  facet_wrap(~ series_name, ncol = 3) +
  labs(
    title = "Forage Fish Time Series Relative to Anchor Hake & Sardine Metrics",
    subtitle = "Solid thick lines = Anchors (x09_dfa_hake_age5plus_doug & x05_dfa_abund_sardine_doug)",
    x = "Year",
    y = "Standardized Value (0 to 1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_forage_comparison)

# -----------------------------------------------------------------------------
# 4. Save Outputs
# -----------------------------------------------------------------------------

write.csv(forage_cor_table, file.path(output_dir, "forage_fish_pairwise_correlations.csv"), row.names = FALSE)
write.csv(anchor_cors, file.path(output_dir, "forage_correlations_with_anchors.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "forage_fish_intercomparison_faceted.png"), p_forage_comparison, width = 12, height = 8)