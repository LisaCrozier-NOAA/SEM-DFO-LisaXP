
#RESULTS- mammal group and fish group
# === GROUP 1 CORRELATION SUB-MATRIX (MARINE MAMMALS) ===
#   > print(round(cor_mat[g1_cols, g1_cols], 3))
#                                   DFA_Steller_s_l_WS_adultyr harbour_seal_DFO_WCVI_smoltyr Killer_whales_NR_BC_adultyr ssl_seak_pup_pred_smoltyr
# DFA_Steller_s_l_WS_adultyr                         1.000                         0.772                       0.864                     0.835
# harbour_seal_DFO_WCVI_smoltyr                      0.772                         1.000                       0.532                     0.818
# Killer_whales_NR_BC_adultyr                        0.864                         0.532                       1.000                     0.649
# ssl_seak_pup_pred_smoltyr                          0.835                         0.818                       0.649                     1.000

# 
# === GROUP 2 CORRELATION SUB-MATRIX (DEMERSAL / SHARK) ===
#   > print(round(cor_mat[g2_cols, g2_cols], 3))
#                                                   DFA_sleeperSharkBSAI_predAK_smoltyr sablefishAK_biomass3plus_smoltyr
# DFA_sleeperSharkBSAI_predAK_smoltyr                               1.000                            0.787
# sablefishAK_biomass3plus_smoltyr                                  0.787                            1.000


# 
# === CROSS-GROUP CORRELATIONS ===
#   > print(round(cor_mat[g1_cols, g2_cols], 3))
#                                 DFA_sleeperSharkBSAI_predAK_smoltyr sablefishAK_biomass3plus_smoltyr
# DFA_Steller_s_l_WS_adultyr                                 -0.631                           -0.466
# harbour_seal_DFO_WCVI_smoltyr                              -0.623                           -0.540
# Killer_whales_NR_BC_adultyr                                -0.653                           -0.327
# ssl_seak_pup_pred_smoltyr                                  -0.420                           -0.289


# ==============================================================================
# Script: 02_analysis/13_group_significant_peff_predators.R
# Purpose: Identify objective groupings among the 6 significant Peff predators
#          (predator-only correlation matrix, clustering, and dimensional checks).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# 0. PATHS & INGESTION
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

peff_ts_path <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")
sweep_file   <- file.path(data_dir, "obj3_data_peff_sweep_all_models_24yr.csv")

peff_ts_df <- read_csv(peff_ts_path, show_col_types = FALSE)
sweep_df   <- read_csv(sweep_file, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 1. PULL THE 6 SIGNIFICANT PREDATORS PROGRAMMATICALLY
# ------------------------------------------------------------------------------
p_col    <- grep("^mod_p$|^pvalue$", names(sweep_df), value = TRUE, ignore.case = TRUE)[1]
name_col <- grep("^predator$|^pred$|^LisaName$", names(sweep_df), value = TRUE, ignore.case = TRUE)[1]

sig_preds_raw <- sweep_df %>%
  filter(!is.na(.data[[p_col]]), .data[[p_col]] <= 0.05) %>%
  pull(.data[[name_col]]) %>%
  unique()

# Match to peff_ts_df columns
ts_cols <- names(peff_ts_df)
pred_cols <- character()

for (p in sig_preds_raw) {
  clean <- gsub("^X[0-9]+_|^Peff_|^Pbase_", "", p)
  hit <- grep(paste0("^Peff_.*", clean), ts_cols, value = TRUE, ignore.case = TRUE)
  if (length(hit) > 0) pred_cols <- union(pred_cols, hit[1])
}

# Simplify column names for display and clustering
clean_labels <- gsub("^Peff_", "", pred_cols)

# ------------------------------------------------------------------------------
# 2. ISOLATE PREDATOR-ONLY TIME SERIES & CORRELATION MATRIX
# ------------------------------------------------------------------------------
df_preds <- peff_ts_df %>%
  select(all_of(pred_cols)) %>%
  set_names(clean_labels)

# 6x6 correlation matrix
cor_mat <- cor(df_preds, use = "pairwise.complete.obs")

cat("\n====================================================================\n")
cat("6x6 PREDATOR-ONLY CORRELATION MATRIX\n")
cat("====================================================================\n")
print(round(cor_mat, 3))

# ------------------------------------------------------------------------------
# 3. INTER-CORRELATION STRENGTH (FINDING CENTRAL/CONNECTED PREDATORS)
# ------------------------------------------------------------------------------
# For each predator, compute mean and sum of absolute correlations with all others
diag(cor_mat) <- NA  # Exclude self-correlation

conn_summary <- tibble(
  Predator       = clean_labels,
  Mean_Abs_r     = round(colMeans(abs(cor_mat), na.rm = TRUE), 3),
  Max_r          = round(apply(cor_mat, 2, function(x) max(x, na.rm = TRUE)), 3),
  Most_Correlated_With = apply(cor_mat, 2, function(x) names(which.max(x)))
) %>%
  arrange(desc(Mean_Abs_r))

cat("\n====================================================================\n")
cat("PREDATOR CONNECTIVITY (RANKED BY AVERAGE ABSOLUTE INTER-CORRELATION)\n")
cat("====================================================================\n")
print(as.data.frame(conn_summary))

# Restore diagonal
diag(cor_mat) <- 1.0

# ------------------------------------------------------------------------------
# 4. DIMENSIONALITY & HIERARCHICAL CLUSTER ANALYSIS
# ------------------------------------------------------------------------------
ev <- eigen(cor_mat)$values
cat("\nEigenvalues of 6x6 Matrix:\n")
print(round(ev, 3))
cat(sprintf("Number of dimensions with eigenvalue > 1.0: %d\n", sum(ev > 1.0)))

# Distance metric: d = 1 - r (variables that vary in the same direction are close)
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "ward.D2")

cat("\n====================================================================\n")
cat("HIERARCHICAL CLUSTERING (Ward's D2, 2-Cluster Solution):\n")
cat("====================================================================\n")
k2_groups <- cutree(hc, k = 2)
for (g in unique(k2_groups)) {
  cat(sprintf("Group %d: %s\n", g, paste(names(k2_groups[k2_groups == g]), collapse = ", ")))
}

if (sum(ev > 1.0) >= 3) {
  cat("\nHIERARCHICAL CLUSTERING (3-Cluster Solution):\n")
  k3_groups <- cutree(hc, k = 3)
  for (g in unique(k3_groups)) {
    cat(sprintf("Group %d: %s\n", g, paste(names(k3_groups[k3_groups == g]), collapse = ", ")))
  }
}

# ------------------------------------------------------------------------------
# 5. PLOT DENDROGRAM
# ------------------------------------------------------------------------------
fig_path <- file.path(proj_dir, "Outputs/figures/peff_predator_dendrogram.png")
png(fig_path, width = 8, height = 5, units = "in", res = 300)
par(mar = c(3, 4, 3, 1))
plot(hc, main = "Dendrogram of Significant Peff Predators (d = 1 - r)",
     xlab = "", sub = "", ylab = "Distance (1 - r)", las = 1)
dev.off()
cat(sprintf("\nSaved cluster dendrogram to: %s\n", fig_path))


#=================================================
#Step 1: Inspect the Within-Group Sign Structures
#=================================
# Group 1: Marine Mammals
g1_cols <- c("DFA_Steller_s_l_WS_adultyr", "harbour_seal_DFO_WCVI_smoltyr", 
             "Killer_whales_NR_BC_adultyr", "ssl_seak_pup_pred_smoltyr")
cat("=== GROUP 1 CORRELATION SUB-MATRIX (MARINE MAMMALS) ===\n")
print(round(cor_mat[g1_cols, g1_cols], 3))

# Group 2: Demersal Fish / Shark
g2_cols <- c("DFA_sleeperSharkBSAI_predAK_smoltyr", "sablefishAK_biomass3plus_smoltyr")
cat("\n=== GROUP 2 CORRELATION SUB-MATRIX (DEMERSAL / SHARK) ===\n")
print(round(cor_mat[g2_cols, g2_cols], 3))

# Cross-group correlation
cat("\n=== CROSS-GROUP CORRELATIONS ===\n")
print(round(cor_mat[g1_cols, g2_cols], 3))
