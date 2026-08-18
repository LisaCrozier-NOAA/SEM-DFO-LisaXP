# ==============================================================================
# Script: Structured Anchor-First Forage Correlation Analysis
# Priority 1: x09_dfa_hake_age5plus_doug vs All
# Priority 2: x05_dfa_abund_sardine_doug vs All
# Priority 3: All Remaining Pairwise Forage Correlations
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# Load merged prey dataset
merged_prey <- read.csv(file.path(output_dir, "jake_vs_doug_prey_merged.csv"))

# Define key anchors
hake_anchor    <- "x09_dfa_hake_age5plus_doug"
sardine_anchor <- "x05_dfa_abund_sardine_doug"

# -----------------------------------------------------------------------------
# 1. Select Targeted Columns
# -----------------------------------------------------------------------------

all_forage_cols <- colnames(merged_prey)[
  grepl("sardine|herring|anchovy|mackerel", colnames(merged_prey), ignore.case = TRUE)
]

target_cols <- unique(c("year", hake_anchor, sardine_anchor, all_forage_cols))

forage_focus <- merged_prey %>%
  select(all_of(target_cols))

# -----------------------------------------------------------------------------
# 2. Build Structured Priority Correlation Table
# -----------------------------------------------------------------------------

metric_cols <- setdiff(target_cols, "year")

# Generate all unique pairwise combinations
all_pairs <- expand_grid(
  var1 = metric_cols,
  var2 = metric_cols
) %>%
  filter(var1 != var2) %>%
  # Ensure unique pairs without directional duplicates
  mutate(pair_id = ifelse(var1 < var2, paste(var1, var2, sep = "___"), paste(var2, var1, sep = "___"))) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

# Calculate Pearson correlations and force anchor column ordering
forage_cor_table <- all_pairs %>%
  rowwise() %>%
  mutate(
    # Force Metric_A to be Hake anchor first, Sardine anchor second
    Metric_A = case_when(
      var1 == hake_anchor | var2 == hake_anchor ~ hake_anchor,
      var1 == sardine_anchor | var2 == sardine_anchor ~ sardine_anchor,
      TRUE ~ var1
    ),
    Metric_B = case_when(
      Metric_A == hake_anchor ~ ifelse(var1 == hake_anchor, var2, var1),
      Metric_A == sardine_anchor ~ ifelse(var1 == sardine_anchor, var2, var1),
      TRUE ~ var2
    ),
    
    # Pearson Correlation Calculation
    pearson_r = cor(forage_focus[[Metric_A]], forage_focus[[Metric_B]], use = "complete.obs"),
    abs_r     = abs(pearson_r)
  ) %>%
  ungroup() %>%
  select(Metric_A, Metric_B, pearson_r, abs_r) %>%
  # Define Priority Grouping for Sorting
  mutate(
    priority_group = case_when(
      Metric_A == hake_anchor ~ 1,
      Metric_A == sardine_anchor ~ 2,
      TRUE ~ 3
    )
  ) %>%
  # Sort by Priority Group first, then by strongest absolute correlation within each group
  arrange(priority_group, desc(abs_r)) %>%
  select(-priority_group)

# -----------------------------------------------------------------------------
# 3. Identify Metrics Strongly Correlated (|r| > 0.69) with Either Anchor
# -----------------------------------------------------------------------------

strongly_correlated_with_anchors <- forage_cor_table %>%
  filter(Metric_A %in% c(hake_anchor, sardine_anchor) & abs_r > 0.69) %>%
  pull(Metric_B) %>%
  unique()

cat("\nMetrics removed because they had |r| > 0.69 with Hake or Sardine anchor:\n")
print(strongly_correlated_with_anchors)

# -----------------------------------------------------------------------------
# 4. Print & Save Structured Output
# -----------------------------------------------------------------------------

cat("\n========================================================\n")
cat(" 1. CORRELATIONS WITH HAKE ANCHOR (x09_dfa_hake_age5plus_doug)\n")
cat("========================================================\n")
print(forage_cor_table %>% filter(Metric_A == hake_anchor), n = 30)

cat("\n========================================================\n")
cat(" 2. CORRELATIONS WITH SARDINE ANCHOR (x05_dfa_abund_sardine_doug)\n")
cat("========================================================\n")
print(forage_cor_table %>% filter(Metric_A == sardine_anchor), n = 30)

cat("\n========================================================\n")
cat(" 3. REMAINING FORAGE-TO-FORAGE CORRELATIONS\n")
cat("    (Excluding metrics with |r| > 0.69 with Hake or Sardine anchors)\n")
cat("========================================================\n")

# Filter out pairs involving any metric that already showed high correlation with anchors
forage_remaining <- forage_cor_table %>%
  filter(Metric_A != hake_anchor & Metric_A != sardine_anchor) %>%
  filter(!Metric_A %in% strongly_correlated_with_anchors & !Metric_B %in% strongly_correlated_with_anchors)

print(forage_remaining, n = 30)

# Save filtered output
write.csv(forage_cor_table, file.path(output_dir, "structured_anchor_forage_correlations.csv"), row.names = FALSE)
write.csv(forage_remaining, file.path(output_dir, "remaining_uncorrelated_forage_pairings.csv"), row.names = FALSE)