#

# ==============================================================================
# Script: 02_analysis/06a_cluster_diagrams_semnodes_v2.R
# Step 4: Hierarchical Clustering (hclust) of Complete Peff Predators
# ==============================================================================
# ==============================================================================
# Step 4b: Hierarchical Clustering of Complete Peff Predators (Filtered)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
})

# Exclude non-salmon predator (pollock) and wrong-sign predator (N. fur seal)
exclude_preds <- c(
  "pollockGOA_SSB_smoltyr__peff_TempPrey2_z",
  "Northern_f_s_WS_adultyr__peff_Base_z"
)

peff_all_cols <- names(sem_data)[grepl("__peff_.*_z$", names(sem_data))]

peff_filtered_cols <- sem_data %>%
  select(all_of(setdiff(peff_all_cols, exclude_preds))) %>%
  keep(~ sum(is.na(.x)) == 0) %>%
  names()

cat("Remaining complete Peff predators:", length(peff_filtered_cols), "\n")
print(peff_filtered_cols)

clean_pred_labels <- peff_filtered_cols %>%
  str_replace("__peff_(.*)_z$", " [\\1]")

pred_mat <- sem_data %>%
  select(all_of(peff_filtered_cols)) %>%
  as.data.frame()
colnames(pred_mat) <- clean_pred_labels

cor_mat <- cor(pred_mat, use = "complete.obs")
cat("\nPairwise Correlation Matrix:\n")
print(round(cor_mat, 2))

dist_cor <- as.dist(1 - cor_mat)
hc_pred  <- hclust(dist_cor, method = "ward.D2")

par(mar = c(3, 4, 3, 20))
plot(
  as.dendrogram(hc_pred),
  horiz = TRUE,
  main  = "Complete Peff Salmon Predators (1998-2021)\nDistance = 1 - r (Ward.D2)",
  xlab  = "Correlation Distance (1 - r)",
  cex   = 0.85
)
abline(v = 0.5, col = "firebrick", lty = 2)

# 1. Print cluster assignments for k = 2, 3, and 4 clusters--------
cluster_table <- tibble(
  predator  = hc_pred$labels,
  k2_group  = cutree(hc_pred, k = 2),
  k3_group  = cutree(hc_pred, k = 3),
  k4_group  = cutree(hc_pred, k = 4)
) %>%
  arrange(k3_group, k4_group, predator)

print(as.data.frame(cluster_table), row.names = FALSE)

# 2. Print sequential merge steps safely
format_node <- function(idx, labels) {
  if (idx < 0) labels[-idx] else paste0("Cluster_", idx)
}

merge_steps <- tibble(
  step   = seq_along(hc_pred$height),
  height = round(hc_pred$height, 3),
  node_1 = map_chr(hc_pred$merge[, 1], format_node, labels = hc_pred$labels),
  node_2 = map_chr(hc_pred$merge[, 2], format_node, labels = hc_pred$labels)
)

cat("\nSequential Merge History (lowest height = most tightly correlated):\n")
print(as.data.frame(merge_steps), row.names = FALSE)
