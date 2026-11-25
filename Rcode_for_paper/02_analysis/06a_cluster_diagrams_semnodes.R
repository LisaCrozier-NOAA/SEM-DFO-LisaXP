
# RESULTS------
# --- PreyNCC ---
#   Cluster 1: DFA_SeaNettle, marketsquid_GAM, DFA_abundSardine, DFA_sumPreyOfPrey_plankt..
# 
# --- PreyAK ---
#   Cluster 1: biomassDecaShelfSum, EGOA_DFA_seabirds, DFA_WGOA_DFA_midTrophic
#   Cluster 2: DFA_EGOA_DFA_lowerTrophic, biomassMysidShelfSum
#   Cluster 3: pinkSalmonNorthAmerica, WGOA_DFA_lowerTrophic, pollockGOA_age1and2
# 
# --- PredAK ---
#   Cluster 1: sablefishAK_biomass3plus, DFA_sleeperSharkBSAI_pre.., Harbour_s_WS_adultyr, DFA_halibutBiomassAge8pl..
#   Cluster 2: Northern_f_s_WS_adultyr, DFA_AllSeaLionsEMB_adultyr
# 


# ==============================================================================
# Script: 02_analysis/06a_cluster_diagrams_semnodes.R
# Purpose: Restore all significant subguilds to each SEMnode:
#          1. Select subguilds from variable importance table
#          2. Examine empirical pairwise correlation within each SEMnode, 
#          3. If there are >4, use cluster analysis split into multiple latents if they seem to be grouped
#          4. Use backward regression to prune latents and main regression to most significant, using p<0.1 cutoff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}


val_ranks_wide <- read_csv("Rcode_for_paper/Outputs/data/obj2_data_variable_importance_wide.csv",show_col_types = FALSE) %>%
  arrange(LisaName)

print(val_ranks_wide,n=Inf) 

library(dplyr)

val <- val_ranks_wide %>%
  # Row-wise maximum across the 4 DAG columns
  rowwise() %>%
  mutate(max = max(c_across(c(DAG1A_long, DAG1A_short, DAG1B_long, DAG1B_short)), na.rm = TRUE)) %>%
  ungroup() %>%
  # Assign SEM node based on Region and LisaName prefixes
  mutate(semnode = case_when(
    Region == "NCC" & grepl("X0[1-5]", LisaName)           ~ "PreyNCC",
    grepl("X06", LisaName)                                  ~ "Growth",
    grepl("X07", LisaName)                                  ~ "Abundance",
    Region == "NCC" & grepl("X08|X09|X11", LisaName)        ~ "PredNCC",
    Region == "AK"  & grepl("X1[2-4]", LisaName)            ~ "PreyAK",
    Region == "AK"  & grepl("X10|X11|X15", LisaName)            ~ "PredAK",
    TRUE                                                    ~ NA_character_
  )) %>%
  filter(max>0.05) %>%
  arrange(LisaName)

  print(val,n=Inf) 
  
print(val_ranks_wide %>% filter(grepl("X14|X13|X12",LisaName)) %>% arrange(desc(DAG1B_long)),n=Inf)

akprey_cor<-read.csv(file.path(data_out_dir,"AKPrey_cor_matrix.csv"),row.names=1);akprey_cor
print(round(akprey_cor, 3))


#Plot cluster diagram---------
# ==============================================================================
# Script: 02_analysis/06c_cluster_diagrams_max2_missing_years.R
# Purpose: Cluster candidate variables within the 4 SEM nodes (PreyNCC, PredNCC, 
#          PreyAK, PredAK), strictly filtering out any series with > 2 missing years
#          across the 1998-2021 period (n = 24 years).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
fig_dir  <- file.path(proj_dir, "Outputs/figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

# Ensure peff_df is strictly 1998-2021 as well
yr_peff_col <- grep("^year$", names(peff_df), value = TRUE, ignore.case = TRUE)
if (length(yr_peff_col) > 0) {
  peff_df <- peff_df %>% filter(.data[[yr_peff_col[1]]] >= 1998 & .data[[yr_peff_col[1]]] <= 2021)
}

val_ranks_wide <- read_csv(file.path(data_dir, "obj2_data_variable_importance_wide.csv"), show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 1. SCREEN FOR IMPORTANCE (max > 0.05) & ASSIGN NODES
# ------------------------------------------------------------------------------
val <- val_ranks_wide %>%
  mutate(
    max = pmax(DAG1A_long, DAG1A_short, DAG1B_long, DAG1B_short, na.rm = TRUE),
    semnode = case_when(
      Region == "NCC" & grepl("X0[1-5]", LisaName)           ~ "PreyNCC",
      grepl("X06", LisaName)                                  ~ "Growth",
      grepl("X07", LisaName)                                  ~ "Abundance",
      Region == "NCC" & grepl("X08|X09|X11", LisaName)        ~ "PredNCC",
      Region == "AK"  & grepl("X1[2-4]", LisaName)            ~ "PreyAK",
      Region == "AK"  & grepl("X10|X11|X15", LisaName)        ~ "PredAK",
      TRUE                                                    ~ NA_character_
    )
  ) %>%
  filter(!is.na(semnode), max > 0.05) %>%
  arrange(semnode, desc(max))

# Helper to fetch and match time series vector
get_timeseries_vector <- function(lisa_name, df_raw_sub, peff_ts_sub) {
  clean_name <- gsub("^X[0-9]+_", "", lisa_name)
  
  # 1. Try Peff variants
  peff_match <- grep(paste0("^Peff_.*", clean_name), names(peff_ts_sub), value = TRUE, ignore.case = TRUE)
  if (length(peff_match) > 0) return(peff_ts_sub[[peff_match[1]]])
  
  peff_match2 <- grep(clean_name, names(peff_ts_sub), value = TRUE, ignore.case = TRUE)
  if (length(peff_match2) > 0) return(peff_ts_sub[[peff_match2[1]]])
  
  # 2. Try raw alt data variants
  raw_match <- grep(paste0("^", clean_name, "$"), names(df_raw_sub), value = TRUE, ignore.case = TRUE)
  if (length(raw_match) > 0) return(df_raw_sub[[raw_match[1]]])
  
  raw_match2 <- grep(clean_name, names(df_raw_sub), value = TRUE, ignore.case = TRUE)
  if (length(raw_match2) > 0) return(df_raw_sub[[raw_match2[1]]])
  
  if (lisa_name %in% names(df_raw_sub)) return(df_raw_sub[[lisa_name]])
  
  return(NULL)
}

# ------------------------------------------------------------------------------
# 2. RUN CLUSTERING PER NODE WITH STRICT MISSINGNESS FILTER (NA <= 2)
# ------------------------------------------------------------------------------
target_nodes <- c("PreyNCC", "PredNCC", "PreyAK", "PredAK")

png_path <- file.path(fig_dir, "fig_cluster_dendrograms_max2_missing.png")
png(png_path, width = 14, height = 11, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))

node_clusters <- list()
filter_audit  <- list()

for (node in target_nodes) {
  vars_in_node <- val %>% filter(semnode == node) %>% pull(LisaName)
  
  cat(sprintf("\n==================================================\nProcessing Node: %s (Screened Candidates = %d)\n==================================================\n", 
              node, length(vars_in_node)))
  
  series_list <- list()
  
  for (v in vars_in_node) {
    ts_vec <- get_timeseries_vector(v, df_24, peff_df)
    
    if (is.null(ts_vec)) {
      filter_audit[[v]] <- tibble(LisaName = v, Node = node, Status = "DROPPED: Not Found in Timeseries", N_Missing = NA_integer_)
      next
    }
    
    n_miss <- sum(is.na(ts_vec))
    
    # Strict missingness filter: keep only if missing <= 2 years out of 24
    if (n_miss <= 2) {
      clean_lbl <- gsub("^X[0-9]+_", "", v)
      short_lbl <- ifelse(nchar(clean_lbl) > 26, paste0(substr(clean_lbl, 1, 24), ".."), clean_lbl)
      series_list[[short_lbl]] <- ts_vec
      filter_audit[[v]] <- tibble(LisaName = v, Node = node, Status = "RETAINED", N_Missing = n_miss)
    } else {
      cat(sprintf("  [-] DROPPED %s: %d missing years (limit is <= 2)\n", v, n_miss))
      filter_audit[[v]] <- tibble(LisaName = v, Node = node, Status = "DROPPED: >2 Missing Years", N_Missing = n_miss)
    }
  }
  
  if (length(series_list) < 2) {
    plot.new()
    title(main = sprintf("%s (n = %d: Insufficient for Clustering)", node, length(series_list)), cex.main = 1.1)
    next
  }
  
  # Aligned dataframe for current node
  df_node <- as.data.frame(series_list)
  
  # Pairwise correlation matrix
  cor_mat <- cor(df_node, use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  
  # Distance matrix based on shared direction: d = 1 - r
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Number of dimensions with eigenvalues > 1.0
  ev <- eigen(cor_mat)$values
  n_dims <- max(1, sum(ev > 1.0))
  
  # Dendrogram plot
  plot(hc, main = sprintf("%s (Retained = %d | Dims [λ>1] = %d)", node, length(series_list), n_dims),
       xlab = "", ylab = "Distance (1 - r)", sub = "", las = 2, cex = 0.85, font.lab = 2)
  
  # Outline identified sub-groups if dimensions >= 2 and at least 4 indicators exist
  k_cut <- min(n_dims, length(series_list) - 1)
  if (k_cut >= 2 && length(series_list) >= 4) {
    rect.hclust(hc, k = k_cut, border = c("red", "blue", "darkgreen", "purple")[1:k_cut])
    node_clusters[[node]] <- cutree(hc, k = k_cut)
  } else {
    node_clusters[[node]] <- cutree(hc, k = 1)
  }
}

mtext("Cluster Dendrograms: SEM Nodes Filtered to Missing <= 2 Years (Ward's D2, 1998-2021)", 
      outer = TRUE, cex = 1.3, font = 2)
dev.off()

cat(sprintf("\nPlot successfully saved to: %s\n", png_path))

# ------------------------------------------------------------------------------
# 3. PRINT CLUSTER SUMMARY
# ------------------------------------------------------------------------------
cat("\n==================================================\nSUMMARY: IDENTIFIED SUB-GROUPS BY NODE\n==================================================\n")
for (node in names(node_clusters)) {
  cat(sprintf("\n--- %s ---\n", node))
  grps <- node_clusters[[node]]
  for (g in sort(unique(grps))) {
    cat(sprintf("  Cluster %d: %s\n", g, paste(names(grps[grps == g]), collapse = ", ")))
  }
}

# Save missingness audit table
audit_df <- bind_rows(filter_audit)
write_csv(audit_df, file.path(data_dir, "obj4_semnode_data_filtering_audit.csv"))
cat(sprintf("\nAudit table saved to: %s\n", file.path(data_dir, "obj4_semnode_data_filtering_audit.csv")))