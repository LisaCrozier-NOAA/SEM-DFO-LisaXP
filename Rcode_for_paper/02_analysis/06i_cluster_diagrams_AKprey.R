# ==============================================================================
# Script: 02_analysis/06i_cluster_diagrams_AKprey.R

# Step 9: Pull A Priori preyAK Candidates (X12, X13, X14 Prefixes Only)
# ==============================================================================
# Cluster Assignments (k = 2, 3, 4) - smoltyr only:
#                               variable k2_group k3_group k4_group
#  X12_DFA_EGOA_DFA_lowerTrophic_smoltyr        1        1        1
#       X12_biomassMysidShelfSum_smoltyr        1        1        1
#                 X12_egoa_krill_smoltyr        1        1        1
#        X12_biomassDecaShelfSum_smoltyr        2        2        2
#    X13_DFA_WGOA_DFA_midTrophic_smoltyr        2        2        2
#          X13_EGOA_DFA_seabirds_smoltyr        2        2        2
#         X14_DFA_pinkSalmonAsia_smoltyr        2        2        2
#      X12_WGOA_DFA_lowerTrophic_smoltyr        1        3        3
#     X14_pinkSalmonNorthAmerica_smoltyr        1        3        3
#             X13_mid_il_capelin_smoltyr        1        3        4
#        X13_pollockGOA_age1and2_smoltyr        1        3        4
# 
# Sequential Merge History (smoltyr only):
#   > print(as.data.frame(preyAK_merges), row.names = FALSE)
# step height                                node_1                             node_2
# 1  0.413   X13_DFA_WGOA_DFA_midTrophic_smoltyr     X14_DFA_pinkSalmonAsia_smoltyr
# 2  0.428 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr   X12_biomassMysidShelfSum_smoltyr
# 3  0.489       X13_pollockGOA_age1and2_smoltyr         X13_mid_il_capelin_smoltyr
# 4  0.568     X12_WGOA_DFA_lowerTrophic_smoltyr X14_pinkSalmonNorthAmerica_smoltyr
# 5  0.605                X12_egoa_krill_smoltyr                          Cluster_2
# 6  0.617         X13_EGOA_DFA_seabirds_smoltyr                          Cluster_1
# 7  0.751       X12_biomassDecaShelfSum_smoltyr                          Cluster_6
# 8  1.070                             Cluster_3                          Cluster_4
# 9  1.793                             Cluster_5                          Cluster_8
# 10  2.351                             Cluster_7                          Cluster_9
# ==============================================================================

# ==============================================================================
# Step 9c: Re-Cluster smoltyr-Only preyAK Candidates (X12-X14) & Save Dendrogram
# ==============================================================================
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
plot_dir <- file.path(proj_dir, "Outputs/plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

# 1. Select X12, X13, and X14 columns EXCLUDING any _adultyr variables
preyAK_smoltyr_cols <- sem_data %>%
  select(matches("^X1(2|3|4)_"), -matches("_adultyr$")) %>%
  names()

cat("Total smoltyr-only X12/X13/X14 preyAK candidates:", length(preyAK_smoltyr_cols), "\n")
print(preyAK_smoltyr_cols)

# 2. Compute pairwise correlation matrix and correlation distance (1 - r)
preyAK_mat <- sem_data %>%
  select(all_of(preyAK_smoltyr_cols)) %>%
  as.data.frame()

cor_preyAK  <- cor(preyAK_mat, use = "complete.obs")
dist_preyAK <- as.dist(1 - cor_preyAK)
hc_preyAK   <- hclust(dist_preyAK, method = "ward.D2")

# Helper function to draw the horizontal dendrogram cleanly
plot_preyAK_dendrogram <- function() {
  par(mar = c(4.5, 2, 3.5, 22)) # Generous right margin so long names never clip
  plot(
    as.dendrogram(hc_preyAK),
    horiz = TRUE,
    main  = "Hierarchical Clustering of smoltyr-Only preyAK Candidates (X12-X14)\nDistance = 1 - r (Ward.D2)",
    xlab  = "Correlation Distance (1 - r)",
    cex   = 0.85
  )
  abline(v = 0.6, col = "firebrick", lty = 2, lwd = 1.5)
}

# 3a. Save high-resolution PNG to disk for manuscript/documentation
png(file.path(plot_dir,"preyAK_smoltyr_dendrogram.png"), width = 2600, height = 1600, res = 300)
plot_preyAK_dendrogram()
dev.off()
cat("\nSaved high-resolution dendrogram to:", file.path(getwd(), "preyAK_smoltyr_dendrogram.png"), "\n")

# 3b. Render in active RStudio / R plot window
plot_preyAK_dendrogram()

# 4. Print cluster groupings (k = 2, 3, 4) and sequential merge history
preyAK_clusters <- tibble(
  variable = hc_preyAK$labels,
  k2_group = cutree(hc_preyAK, k = 2),
  k3_group = cutree(hc_preyAK, k = 3),
  k4_group = cutree(hc_preyAK, k = 4)
) %>%
  arrange(k3_group, k4_group, variable)

cat("\nCluster Assignments (k = 2, 3, 4) - smoltyr only:\n")
print(as.data.frame(preyAK_clusters), row.names = FALSE)

format_node <- function(idx, labels) {
  if (idx < 0) labels[-idx] else paste0("Cluster_", idx)
}

preyAK_merges <- tibble(
  step   = seq_along(hc_preyAK$height),
  height = round(hc_preyAK$height, 3),
  node_1 = map_chr(hc_preyAK$merge[, 1], format_node, labels = hc_preyAK$labels),
  node_2 = map_chr(hc_preyAK$merge[, 2], format_node, labels = hc_preyAK$labels)
)

cat("\nSequential Merge History (smoltyr only):\n")
print(as.data.frame(preyAK_merges), row.names = FALSE)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# 1. Extract all columns starting with X12_, X13_, or X14_ and check missingness
preyAK_candidates <- sem_data %>%
  select(matches("^X1(2|3|4)_")) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing")

cat("A priori preyAK candidates (X12, X13, X14):\n")
print(as.data.frame(preyAK_candidates), row.names = FALSE)

# 2. Filter to 0-missing-year indicators first (preserves fixed.x = TRUE)
preyAK_0na <- preyAK_candidates %>%
  filter(n_missing == 0) %>%
  pull(variable)

cat("\nComplete (0-NA) preyAK indicators to test:", length(preyAK_0na), "\n")
print(preyAK_0na)

# ==============================================================================
# Step 9b: Hierarchical Clustering (hclust) of A Priori preyAK Candidates (X12-X14)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
})

# 1. Select all X12, X13, and X14 columns (all have 0 missing years)
preyAK_cols <- sem_data %>%
  select(matches("^X1(2|3|4)_")) %>%
  names()

cat("Total X12/X13/X14 preyAK candidates:", length(preyAK_cols), "\n")

# 2. Compute pairwise correlation matrix and correlation distance (1 - r)
preyAK_mat <- sem_data %>%
  select(all_of(preyAK_cols)) %>%
  as.data.frame()

cor_preyAK  <- cor(preyAK_mat, use = "complete.obs")
dist_preyAK <- as.dist(1 - cor_preyAK)
hc_preyAK   <- hclust(dist_preyAK, method = "ward.D2")

# 3. Plot horizontal dendrogram
par(mar = c(3, 4, 3, 22)) # Extra right margin for variable names
plot(
  as.dendrogram(hc_preyAK),
  horiz = TRUE,
  main  = "Hierarchical Clustering of preyAK Candidates (X12, X13, X14)\nDistance = 1 - r (Ward.D2)",
  xlab  = "Correlation Distance (1 - r)",
  cex   = 0.8
)
abline(v = 0.5, col = "firebrick", lty = 2)

# 4. Print cluster groupings (k = 2, 3, 4) and sequential merge history
preyAK_clusters <- tibble(
  variable = hc_preyAK$labels,
  k2_group = cutree(hc_preyAK, k = 2),
  k3_group = cutree(hc_preyAK, k = 3),
  k4_group = cutree(hc_preyAK, k = 4)
) %>%
  arrange(k3_group, k4_group, variable)

cat("\nCluster Assignments (k = 2, 3, 4):\n")
print(as.data.frame(preyAK_clusters), row.names = FALSE)

format_node <- function(idx, labels) {
  if (idx < 0) labels[-idx] else paste0("Cluster_", idx)
}

preyAK_merges <- tibble(
  step   = seq_along(hc_preyAK$height),
  height = round(hc_preyAK$height, 3),
  node_1 = map_chr(hc_preyAK$merge[, 1], format_node, labels = hc_preyAK$labels),
  node_2 = map_chr(hc_preyAK$merge[, 2], format_node, labels = hc_preyAK$labels)
)

cat("\nSequential Merge History:\n")
print(as.data.frame(preyAK_merges), row.names = FALSE)

#cluster results--------

# Cluster Assignments (k = 2, 3, 4):
#    print(as.data.frame(preyAK_clusters), row.names = FALSE)
#                               variable k2_group k3_group k4_group
#  X12_DFA_EGOA_DFA_lowerTrophic_smoltyr        1        1        1
#       X12_biomassMysidShelfSum_smoltyr        1        1        1
#                 X12_egoa_krill_adultyr        1        1        1
#                 X12_egoa_krill_smoltyr        1        1        1
#             X13_mid_il_capelin_adultyr        1        1        1
#        X12_biomassDecaShelfSum_smoltyr        2        2        2
#    X13_DFA_WGOA_DFA_midTrophic_smoltyr        2        2        2
#          X13_EGOA_DFA_seabirds_smoltyr        2        2        2
#         X14_DFA_pinkSalmonAsia_smoltyr        2        2        2
#      X12_WGOA_DFA_lowerTrophic_smoltyr        1        3        3
#     X14_pinkSalmonNorthAmerica_smoltyr        1        3        3
#             X13_mid_il_capelin_smoltyr        1        3        4
#        X13_pollockGOA_age1and2_smoltyr        1        3        4
#   
#        
# Sequential Merge History:
#   > print(as.data.frame(preyAK_merges), row.names = FALSE)
# step height                                node_1                             node_2
# 1  0.413   X13_DFA_WGOA_DFA_midTrophic_smoltyr     X14_DFA_pinkSalmonAsia_smoltyr
# 2  0.428 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr   X12_biomassMysidShelfSum_smoltyr
# 3  0.489       X13_pollockGOA_age1and2_smoltyr         X13_mid_il_capelin_smoltyr
# 4  0.528                X12_egoa_krill_smoltyr         X13_mid_il_capelin_adultyr
# 5  0.568     X12_WGOA_DFA_lowerTrophic_smoltyr X14_pinkSalmonNorthAmerica_smoltyr
# 6  0.617         X13_EGOA_DFA_seabirds_smoltyr                          Cluster_1
# 7  0.689                             Cluster_2                          Cluster_4
# 8  0.751       X12_biomassDecaShelfSum_smoltyr                          Cluster_6
# 9  0.950                X12_egoa_krill_adultyr                          Cluster_7
# 10  1.070                             Cluster_3                          Cluster_5
# 11  1.812                             Cluster_9                         Cluster_10
# 12  2.454                             Cluster_8                         Cluster_11
