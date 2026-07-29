#plot short time series in PCA map--------------

guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1);head(guild.dfas1)
sem_master_data<-read.csv("LisaXP/outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
guild.dfas1$sarSR<-sem_master_data$sarSR
importance_topvar<-read.csv("LisaXP/outputs_4/importance_topvar_daic3.csv",row.names=1);head(importance_topvar)
dredge_coef<-read.csv(file = "LisaXP/outputs_4/Top_100_dredgemodels_r2.csv", row.names = NULL);head(dredge_coef)
dredge_coef<-dredge_coef %>%
  arrange(Adj_AIC2) %>%
  mutate(rank_aic2=1:nrow(dredge_coef))


      library(tidyverse)
      library(ggrepel)

library(tidyverse)
#Rank variables by their correlation with SAR
          # 1. Define your response variable
          response_var <- "sarSR" # Or "X16_SAR" depending on which response you are targeting
          
          # 2. Extract all candidates from your lookup metadata
          candidates <- var_lookup_NCC_AK %>% 
            filter(!node_id %in% c("X06", "X07", "X16")) %>% 
            pull(Lisaname)
          
          # 3. Calculate each candidate's raw pairwise correlation with your response variable
          sar_correlations <- cor(guild.dfas1[, c(response_var, candidates)], use = "pairwise.complete.obs")[, response_var] %>%
            enframe(name = "Lisaname", value = "SAR_Correlation") %>%
            filter(Lisaname != response_var) # Remove the self-correlation line
          
          # 4. Bind the metadata (SEMnode) to create the ranked_variables object
          ranked_variables <- var_lookup_NCC_AK %>%
            select(Lisaname, SEMnode) %>%
            left_join(sar_correlations, by = "Lisaname") %>%
            # Sort them from strongest correlation to weakest correlation
            arrange(desc(abs(SAR_Correlation)))
          
          # Quick check to ensure it worked
          head(ranked_variables)
      
      # 1. Slice to post-2002 and drop the final column of the dataframe
      raw_19yr_subset <- guild.dfas1 %>% 
        filter(year > 2002) %>%
        # Drop the very last column
        select(-ncol(.)) %>%
        # Make sure we only keep the numeric variables (excluding 'year')
        select(-year)
      
      # 2. Scale and transpose the clean variables
      cleaned_19yr_data <- raw_19yr_subset %>%
        scale() %>%
        t() # Transpose: Variables are now rows, years are columns
      
      # 3. Run the PCA on the clean shapes
      pca_shapes <- prcomp(cleaned_19yr_data, scale. = FALSE)
      
      # 4. Bind PCA coordinates back to your metadata for plotting
      pca_data <- as.data.frame(pca_shapes$x) %>%
        rownames_to_column("Lisaname") %>%
        left_join(ranked_variables %>% select(Lisaname, SEMnode, SAR_Correlation), by = "Lisaname")
      
      # 5. Plot the multidimensional space
      pca.plot<-
        ggplot(pca_data, aes(x = PC1, y = PC2, color = SEMnode, label = Lisaname)) +
        geom_point(aes(size = abs(SAR_Correlation)), alpha = 0.7) +
        geom_text_repel(size = 3, max.overlaps = 15) +
        scale_size_continuous(name = "Corr with SAR") +
        theme_minimal() +
        labs(
          title = "Multidimensional Space: Trajectory Clustering",
          subtitle = "Variables grouped close together share nearly identical 19-year patterns",
          x = "PC1 (Main Trend / Regime Shape)",
          y = "PC2 (Secondary Trend / Interannual Variance)"
        )
      
      ggsave("LisaXP/outputs_4/pca.plot.png",pca.plot,units="in",width=14,height=8)
      write.csv(pca_data,"LisaXP/outputs_4/pca__guild.dfas1_data.csv",row.names = FALSE)
      
      pca_data
      pca.plot

#plot long time series in MDS Map----------------

#all years

          library(tidyverse)
          library(ggrepel)
          
          # 1. Prepare the full 24-year dataset (dropping only 'year' and the final column)
          full_data <- guild.dfas1 %>%
            select(-ncol(.)) %>% # Drop the final column
            select(-year)        # Drop the timeline tracker
          
          # 2. Compute a pairwise complete correlation matrix across all 24 years
          # This handles the NAs by pairing whatever years are available for each combo
          cor_matrix_full_years <- cor(full_data, use = "pairwise.complete.obs")
          
          # 3. Convert correlations to a physical distance metric (1 - |r|)
          # Absolute correlation ensures both strong positive and negative drivers cluster together
          dist_matrix <- as.dist(1 - abs(cor_matrix_full_years))
          
          # 4. Run Classical Multidimensional Scaling (MDS / PCoA)
          # This finds the best 2D coordinates to represent these multi-year distances
          mds_fit <- cmdscale(dist_matrix, k = 2)
          
          # 5. Extract coordinates and join your guild metadata & SAR correlations
          mds_coordinates <- as.data.frame(mds_fit) %>%
            rownames_to_column("Lisaname") %>%
            rename(Dim1 = V1, Dim2 = V2)
          
          # Calculate each variable's correlation with X16_SAR over the whole time series
          sar_full_correlations <- cor(full_data, use = "pairwise.complete.obs")[, "X16_SAR"] %>%
            enframe(name = "Lisaname", value = "SAR_Correlation")
          
          # Combine everything for plotting
          plot_data_mds <- mds_coordinates %>%
            left_join(ranked_variables %>% select(Lisaname, SEMnode), by = "Lisaname") %>%
            left_join(sar_full_correlations, by = "Lisaname") %>%
            # Label variables without a guild as Salmon Diagnostics
            mutate(SEMnode = if_else(is.na(SEMnode), "Salmon / Core Diagnostics", SEMnode))
          
          # 6. Plot the 24-Year Ecosystem Space
          mds.plot<-
            ggplot(plot_data_mds, aes(x = Dim1, y = Dim2, color = SEMnode, label = Lisaname)) +
            geom_point(aes(size = abs(SAR_Correlation)), alpha = 0.75) +
            geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5) +
            scale_size_continuous(name = "Corr with SAR (r)", range = c(2, 8)) +
            theme_minimal(base_size = 11) +
            labs(
              title = "Ecosystem Trajectory Map (Full 24-Year Time Series)",
              subtitle = "MDS Ordination: Variables clustered close together share highly matching 24-year profiles",
              x = "Dimension 1 (Main multi-year trend)",
              y = "Dimension 2 (Secondary trajectory changes)"
            ) +
            theme(
              legend.position = "right",
              plot.title = element_text(face = "bold")
            )
          
          plot_data_mds_guild.dfas1_X16SARcor<-plot_data_mds
          mds.plot_guild.dfas1_X16SARcor<-mds.plot
          
          write.csv(plot_data_mds_guild.dfas1_X16SARcor,"LisaXP/outputs_4/mds_coordinates__guild.dfas1_data.csv",row.names = FALSE)

          ggsave("LisaXP/outputs_4/mds.plot.png",mds.plot_guild.dfas1_X16SARcor,units="in",width=14,height=8)


#repeat MDS but with raw data, not just DFAs---------
          library(tidyverse)
          library(ggrepel)
          
          dat<-read.csv("LisaXP/outputs_4/sem_master_data.csv")
          dat$sarSR<-guild.dfas1$sarSR
          
          # 1. Prepare the full 24-year dataset (dropping only 'year' and the final column)
          full_data <- dat %>%
            select(-year)        # Drop the timeline tracker
          
          # 2. Compute a pairwise complete correlation matrix across all 24 years
          # This handles the NAs by pairing whatever years are available for each combo
          cor_matrix_full_years <- cor(full_data, use = "pairwise.complete.obs")
          
          # 3. Convert correlations to a physical distance metric (1 - |r|)
          # Absolute correlation ensures both strong positive and negative drivers cluster together
          dist_matrix <- as.dist(1 - abs(cor_matrix_full_years))
          
          # 4. Run Classical Multidimensional Scaling (MDS / PCoA)
          # This finds the best 2D coordinates to represent these multi-year distances
          mds_fit <- cmdscale(dist_matrix, k = 2)
          
          # 5. Extract coordinates and join your guild metadata & SAR correlations
          mds_coordinates <- as.data.frame(mds_fit) %>%
            rownames_to_column("Lisaname") %>%
            rename(Dim1 = V1, Dim2 = V2)
          
          # Calculate each variable's correlation with X16_SAR over the whole time series
          # sar_full_correlations <- cor(full_data, use = "pairwise.complete.obs")[, "X16_SAR"] %>%
          #   enframe(name = "Lisaname", value = "SAR_Correlation")
          sar_full_correlations <- cor(full_data, use = "pairwise.complete.obs")[, "sarSR"] %>%
            enframe(name = "Lisaname", value = "SR_SAR_Correlation")
          
          # Combine everything for plotting
          plot_data_mds <- mds_coordinates %>%
            left_join(ranked_variables %>% select(Lisaname, SEMnode), by = "Lisaname") %>%
            left_join(sar_full_correlations, by = "Lisaname") %>%
            # Label variables without a guild as Salmon Diagnostics
            mutate(SEMnode = if_else(is.na(SEMnode), "Salmon / Core Diagnostics", SEMnode))
          
          
          # 6. Plot the 24-Year Ecosystem Space
          
          mds.raw.SAR_SR.plot<-
            ggplot(plot_data_mds, aes(x = Dim1, y = Dim2, color = SEMnode, label = Lisaname)) +
            geom_point(aes(size = abs(SR_SAR_Correlation)), alpha = 0.75) +
            geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5) +
            scale_size_continuous(name = "Corr with SAR (r)", range = c(2, 8)) +
            theme_minimal(base_size = 11) +
            labs(
              title = "Ecosystem Trajectory Map (Full 24-Year Time Series)",
              subtitle = "MDS Ordination: Variables clustered close together share highly matching 24-year profiles",
              x = "Dimension 1 (Main multi-year trend)",
              y = "Dimension 2 (Secondary trajectory changes)"
            ) +
            theme(
              legend.position = "right",
              plot.title = element_text(face = "bold")
            )
          
          print(mds.raw.SAR_SR.plot)
          plot_data_mds_sem_master_data<-plot_data_mds
          write.csv(plot_data_mds,"LisaXP/outputs_4/mds_coordinates__sem_master_data.csv",row.names = FALSE)
          
          ggsave(filename="LisaXP/outputs_4/mds.rawSAR_SR.plot.png",plot=mds.raw.SAR_SR.plot,units="in",width=14,height=8)
          


#expore results------------------
      pca.plot; head(pca_data)
      mds.plot_guild.dfas1_X16SARcor;head(plot_data_mds_guild.dfas1_X16SARcor)
      mds.raw.SAR_SR.plot; head(plot_data_mds_sem_master_data)
      
      head(plot_data_mds)
      plot_data_mds %>% arrange(desc(SR_SAR_Correlation)) %>% head()
      plot_data_mds %>% filter(SEMnode=="Salmon / Core Diagnostics")
      
      plot_data_mds %>%
        mutate()
      
      pca_data %>% filter(Lisaname=="X13_ammod_WAI")
      
      
#ID nearest neighbors to selected preyAK to explain neg effects--------
      library(tidyverse)
      
      importance_topvar %>% filter(node=="PreyAKindNames") %>% arrange(desc(maxImport)) %>% select(node,Lisaname,total_count,maxImport)
      #               node                             Lisaname total_count maxImport
      # 1  PreyAKindNames          X12_DFA_biomassEuphShelfSum          17      0.93
      # 2  PreyAKindNames                  X12_copepodCom_EGoA           8      0.69
      # 3  PreyAKindNames                       X14_pinkSalmon           1      0.52
      # 4  PreyAKindNames              X12_copepodBiomass_WGoA           7      0.32
      # 5  PreyAKindNames                   X14_pinkSalmonAsia           1      0.17
      # 6  PreyAKindNames                     X13_hexagram_EAI           5      0.16
      # 7  PreyAKindNames X13_pollockBiomassGoAage3plus_predAK           6      0.15
      # 8  PreyAKindNames          X13_DFA_WGOA_DFA_midTrophic           1      0.14
      # 9  PreyAKindNames            X13_DFA_WGOA_DFA_seabirds           3      0.09
      # 10 PreyAKindNames                        X13_gadid_WAI           2      0.06
      # 11 PreyAKindNames                     X13_hexagram_WAI           1      0.04
      
      dredge_coef %>% filter(grepl("X13_ammod_WAI",Predictors)) %>% 
        select(rank_aic2,Formula,All_Coefficients_Sig,Coefficients) %>%
        head()
      #always negative
      dredge_coef %>% filter(grepl("X12_DFA_biomassEuphShelfSum",Predictors)) %>% 
        select(rank_aic2,Formula,All_Coefficients_Sig,Coefficients) %>%
        head()
      #always negative
      
    #use guild.dfa1s, not raw data
      
      source("LisaXP/functions/find_mds_neighbors.r")
      dat=plot_data_mds_guild.dfas1_X16SARcor
      
      target_prey_vars <- dat %>% 
        filter(SEMnode == "PreyAK") %>% 
        pull(Lisaname)
      
      # 3. Run this for all PreyAK variables and bind the results
      preyAK_neighborhoods <- map_dfr(target_prey_vars, ~find_mds_neighbors(.x, dat))
      
      # 4. View the table
      head(preyAK_neighborhoods)
      preyAK_neighborhoods %>% filter(Target_Variable=="X12_DFA_biomassEuphShelfSum")
#               Target_Variable                                 Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
# 1 X12_DFA_biomassEuphShelfSum                     X11_DFA_Harbour_p_WS PredNCC           Negative (-)                    -0.38    0.084      -0.2911234
# 2 X12_DFA_biomassEuphShelfSum X03_DFA_comMurreDietHerrSard_Yaquina_NCC PreyNCC           Negative (-)                    -0.11    0.133      -0.2211710
# 3 X12_DFA_biomassEuphShelfSum                       X02_DFA_SeaNettle_ PreyNCC           Negative (-)                    -0.70    0.139      -0.0565221
# 4 X12_DFA_biomassEuphShelfSum                   X01_Red_necked_phal_WS PreyNCC           Positive (+)                     0.57    0.171       0.1701607
# 5 X12_DFA_biomassEuphShelfSum                X15_salmonSharkGoA_predAK  PredAK           Positive (+)                     0.55    0.187      -0.1808523      preyAK_neighborhoods %>% filter(abs(Pairwise_Cor_With_Target) >0.7)
      
      
      #                         Target_Variable                                  Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
      # 1                   X12_copepodCom_WGoA                X04_CaliforniaMarketSquid_ PreyNCC           Positive (+)                     0.80    0.032     -0.16338131
      # 2                   X12_copepodCom_WGoA           X15_sablefishRecruitment_predAK  PredAK           Positive (+)                     0.82    0.049     -0.09123130
      # 3                   X12_copepodCom_WGoA               X13_DFA_WGOA_DFA_midTrophic  PreyAK           Positive (+)                     0.89    0.166     -0.13415185
      # 4                         X13_ammod_WAI           X15_DFA_sleeperSharkBSAI_predAK  PredAK           Positive (+)                     0.76    0.070     -0.28817694
      # 5                         X13_ammod_WAI                        X06_Lmu_IntSprJunH  Growth           Negative (-)                    -0.93    0.117      0.36019876
      # 6                         X13_ammod_WAI                  X10_Harbour_s_2yrLead_WS  PredAK           Positive (+)                     0.72    0.145     -0.49164283
      # 7  X13_pollockBiomassGoAage3plus_predAK           X15_DFA_sleeperSharkBSAI_predAK  PredAK           Negative (-)                    -0.76    0.058     -0.28817694
      # 8  X13_pollockBiomassGoAage3plus_predAK            X10_Californian_s_l_2yrLead_WS  PredAK           Positive (+)                     0.75    0.083      0.20664192
      # 9             X13_DFA_WGOA_DFA_seabirds                      X01_NHLlogSum_win_25 PreyNCC           Negative (-)                    -0.74    0.070      0.37382966
      # 10          X13_DFA_WGOA_DFA_midTrophic                X04_CaliforniaMarketSquid_ PreyNCC           Positive (+)                     0.88    0.138     -0.16338131
      # 11          X13_DFA_WGOA_DFA_midTrophic           X15_sablefishRecruitment_predAK  PredAK           Positive (+)                     0.91    0.162     -0.09123130
      # 12          X13_DFA_WGOA_DFA_midTrophic                       X12_copepodCom_WGoA  PreyAK           Positive (+)                     0.89    0.166     -0.27255596
      # 13          X13_DFA_WGOA_DFA_midTrophic                        X14_pinkSalmonAsia  PreyAK           Positive (+)                     0.74    0.226      0.18685792
      # 14                       X14_pinkSalmon X15_halibutBiomassAge8plus_2yrLead_predAK  PredAK           Negative (-)                    -0.74    0.061      0.02418219
      # 15                       X14_pinkSalmon                        X14_pinkSalmonAsia  PreyAK           Positive (+)                     0.94    0.088      0.18685792
      # 16                   X14_pinkSalmonAsia                            X14_pinkSalmon  PreyAK           Positive (+)                     0.94    0.088      0.11622371
      # 17                   X14_pinkSalmonAsia X15_halibutBiomassAge8plus_2yrLead_predAK  PredAK           Negative (-)                    -0.72    0.133      0.02418219
      
      #so preyAK could be neg bec----------
      
      #1. it is positively correlated with a predator or competitor
            preyAK_neighborhoods %>% filter(abs(Pairwise_Cor_With_Target) >0.7, Relationship_To_Target=="Positive (+)" ) %>% filter(Lisaname!="X12_copepodCom_WGoA") 
                  # Target_Variable                        Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
                  # 1                   X12_copepodCom_WGoA      X04_CaliforniaMarketSquid_ PreyNCC           Positive (+)                     0.80    0.032      -0.1633813
                  # 2                   X12_copepodCom_WGoA X15_sablefishRecruitment_predAK  PredAK           Positive (+)                     0.82    0.049      -0.0912313
                  # 3                   X12_copepodCom_WGoA     X13_DFA_WGOA_DFA_midTrophic  PreyAK           Positive (+)                     0.89    0.166      -0.1341518
                  # 4                         X13_ammod_WAI X15_DFA_sleeperSharkBSAI_predAK  PredAK           Positive (+)                     0.76    0.070      -0.2881769
                  # 5                         X13_ammod_WAI        X10_Harbour_s_2yrLead_WS  PredAK           Positive (+)                     0.72    0.145      -0.4916428
                  # 6  X13_pollockBiomassGoAage3plus_predAK  X10_Californian_s_l_2yrLead_WS  PredAK           Positive (+)                     0.75    0.083       0.2066419
                  # 7           X13_DFA_WGOA_DFA_midTrophic      X04_CaliforniaMarketSquid_ PreyNCC           Positive (+)                     0.88    0.138      -0.1633813
                  # 8           X13_DFA_WGOA_DFA_midTrophic X15_sablefishRecruitment_predAK  PredAK           Positive (+)                     0.91    0.162      -0.0912313
                  # 9           X13_DFA_WGOA_DFA_midTrophic              X14_pinkSalmonAsia  PreyAK           Positive (+)                     0.74    0.226       0.1868579
                  # 10                       X14_pinkSalmon              X14_pinkSalmonAsia  PreyAK           Positive (+)                     0.94    0.088       0.1868579
                  # 11                   X14_pinkSalmonAsia                  X14_pinkSalmon  PreyAK           Positive (+)                     0.94    0.088       0.1162237
      
      
      #2. it is negatively correlated with other prey or growth
            preyAK_neighborhoods %>% filter(abs(Pairwise_Cor_With_Target) >0.7, Relationship_To_Target=="Negative (-)" ) %>% 
                filter(Lisaname!="X15_DFA_sleeperSharkBSAI_predAK") %>%
              filter(Target_Variable !="X14_pinkSalmon")  %>%
              filter(Target_Variable !="X14_pinkSalmonAsia")
            
            #             Target_Variable             Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
            # 1             X13_ammod_WAI   X06_Lmu_IntSprJunH  Growth           Negative (-)                    -0.93    0.117       0.3601988
            # 2 X13_DFA_WGOA_DFA_seabirds X01_NHLlogSum_win_25 PreyNCC           Negative (-)                    -0.74    0.070       0.3738297
             
            
#so predNCC and predAK could be pos bec----------
            
            dat=plot_data_mds_guild.dfas1_X16SARcor
            
            target_prey_vars <- dat %>% 
              filter(SEMnode == "PredNCC") %>% 
              pull(Lisaname)
            
            predNCC_neighborhoods <- map_dfr(target_prey_vars, ~find_mds_neighbors(.x, dat))
            print(predNCC_neighborhoods)
            
            #1. it is positively correlated with prey
            predNCC_neighborhoods %>% filter(abs(Pairwise_Cor_With_Target) >0.7, Relationship_To_Target=="Positive (+)" ) # %>% filter(Lisaname!="X12_copepodCom_WGoA") 
            