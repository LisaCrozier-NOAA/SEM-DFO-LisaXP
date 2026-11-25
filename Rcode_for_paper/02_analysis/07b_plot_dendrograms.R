

# ==============================================================================
# Script: 02_analysis/07b_plot_dendrograms.R
# Step 2: Rationale for latent groupings 
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(tidySEM)
  library(ggplot2)
})

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
table_out_dir <- file.path(master_out, "tables")
fig_out_dir <- file.path(master_out, "figures")
plot_out_dir <- file.path(master_out, "plots")
data_in_path <- file.path(data_out_dir, "sem_master_peff_zscored.csv")

source("functions/sem_utils.R")
source("functions/Lisa_sem_graph.v3_fxn.r")
source("functions/sem_output_summary_fxn.r")

sem_data <- read_csv(data_in_path, show_col_types = FALSE)
names(sem_data)

# ==============================================================================
# Step 4b: Hierarchical Clustering of Exogenous node indicators 
# ==============================================================================

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
          
          
          
#HELPER FUNCTION FOR PRINTING PLOTS -----------
          # Helper function to draw the horizontal dendrogram cleanly
          plot_dendrogram <- function(hclust_output,mytitle="Hierarchical Clustering of SEMnodes",crop_height=0.6) {
            par(mar = c(4.5, 2, 3.5, 22)) # Generous right margin so long names never clip
            plot(
              as.dendrogram(hclust_output),
              horiz = TRUE,
              main  =  mytitle,
              xlab  = "Correlation Distance (1 - r), (Ward.D2)",
              cex   = 0.85
            )
            abline(v = crop_height, col = "firebrick", lty = 2, lwd = 1.5)
          }
          
          format_node <- function(idx, labels) {
            if (idx < 0) labels[-idx] else paste0("Cluster_", idx)
          }
          
          
#PLOT PREDAK CLUSTER DENDROGRAM---------
          png(file.path(plot_out_dir,"predAK_dendrogram.png"), width = 2600, height = 1600, res = 300)
              plot_dendrogram(hclust_output=as.dendrogram(hc_pred),mytitle = "Peff Salmon Predators",crop_height=0.6)
          dev.off()
          
      # 1. Print cluster assignments for k = 2, 3, and 4 clusters
            cluster_table <- tibble(
              predator  = hc_pred$labels,
              k2_group  = cutree(hc_pred, k = 2),
              k3_group  = cutree(hc_pred, k = 3),
              k4_group  = cutree(hc_pred, k = 4)
            ) %>%
              arrange(k3_group, k4_group, predator)
            
            print(as.data.frame(cluster_table), row.names = FALSE)
            
            # 2. Print sequential merge steps safely

            merge_steps <- tibble(
              step   = seq_along(hc_pred$height),
              height = round(hc_pred$height, 3),
              node_1 = map_chr(hc_pred$merge[, 1], format_node, labels = hc_pred$labels),
              node_2 = map_chr(hc_pred$merge[, 2], format_node, labels = hc_pred$labels)
            )
    # PredAK cluster results----        
            cat("\nSequential Merge History (lowest height = most tightly correlated):\n")
            print(as.data.frame(merge_steps), row.names = FALSE)
            
            #Results summary: 5 mammals cluster together, incl all sea lions, kw, and hs WCVI
            #group 2 has an inner group w/ just sleepers and sablefish, and an outer group w/ halibut, hs coastwide and hs WS
            #decision: test inner group, bec outer group was NS in peff table
            
 #                                              predator k2_group k3_group k4_group
 #                 DFA_AllSeaLionsEMB_adultyr [AltPrey1]        1        1        1
 #                DFA_Steller_s_l_WS_adultyr [TempPrey1]        1        1        1
 #                    Killer_whales_NR_BC_adultyr [Temp]        1        1        1
 #              harbour_seal_DFO_WCVI_smoltyr [AltPrey1]        1        1        1
 #                  ssl_seak_pup_pred_smoltyr [AltPrey2]        1        1        1
 #       DFA_sleeperSharkBSAI_predAK_smoltyr [TempPrey1]        2        2        2
 #        harbour_seal_DFO_Coastwide_smoltyr [TempPrey1]        2        2        2
 #          sablefishAK_biomass3plus_smoltyr [TempPrey2]        2        2        2
 # DFA_halibutBiomassAge8plus_predAK_adultyr [TempPrey1]        2        2        3
 #                           Harbour_s_WS_adultyr [Base]        2        2        3
 #             salmonSharkGoA_predAK_smoltyr [TempPrey1]        1        3        4
            
            
#PLOT PREYAK CLUSTER DENDROGRAM---------
            
            
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

      png(file.path(plot_out_dir,"preyAK_dendrogram.png"), width = 2600, height = 1600, res = 300)
            plot_dendrogram(hclust_output=as.dendrogram(hc_preyAK),mytitle = "PreyAK",crop_height=0.8)
      dev.off()
            
  
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
            
            preyAK_merges <- tibble(
              step   = seq_along(hc_preyAK$height),
              height = round(hc_preyAK$height, 3),
              node_1 = map_chr(hc_preyAK$merge[, 1], format_node, labels = hc_preyAK$labels),
              node_2 = map_chr(hc_preyAK$merge[, 2], format_node, labels = hc_preyAK$labels)
            )
            
            cat("\nSequential Merge History:\n")
            print(as.data.frame(preyAK_merges), row.names = FALSE)
            
            #PreyAK cluster results--------
            
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
            
            
            
                        
            