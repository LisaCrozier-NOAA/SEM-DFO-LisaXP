
#DUMP -- first try

# 1. Isolate the complex models using your preferred flat penalty (e.g., Mult_0.4)
dredge_results_noHCI<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results_combined_noHCI.csv",row.names = NULL);head(dredge_results)
head(dredge_results_noHCI)

complex_models_only <- dredge_results_noHCI %>%
  filter(Num_Params %in% c(3, 4)) %>%
  mutate(
    # Applying the stable flat penalty at your sensitivity sweet spot
    Adjusted_AIC = Raw_AIC + (1.0  * Missing_Yrs),
    Adj_AIC2 = Raw_AIC + (2.0 *  Missing_Yrs)
  ) %>%
  # 2. Add 4 clean binary columns to check for specific guild presence
  mutate(
    PreyNCC = if_else(grepl("PreyNCC", Guild_Sources), 1, 0),
    PredNCC = if_else(grepl("PredNCC", Guild_Sources), 1, 0),
    PreyAK  = if_else(grepl("PreyAK",  Guild_Sources), 1, 0),
    PredAK  = if_else(grepl("PredAK",  Guild_Sources), 1, 0)
  ) %>%
  rowwise() %>%
  mutate(Nguilds=sum(c_across(PreyNCC:PredAK))) %>%
  
  arrange(Adj_AIC2)


#Define proxies======================
        # --- 1. Calculate the Pairwise Correlation Matrix ---
        cor_matrix_full <- cor(guild.dfas1[, all_calc_vars], use = "pairwise.complete.obs")
        
        # Define stand-in correlation threshold (e.g., |r| >= 0.5)
        # You can change this threshold to make the test more or less conservative
        cor_threshold <- 0.50
        
        growth_proxies <- names(which(abs(cor_matrix_full["X06_Lmu_IntSprJunH", ]) >= cor_threshold))
        abundance_proxies <- names(which(abs(cor_matrix_full["X07_DFA_cpue_IntSprJunHW", ]) >= cor_threshold))

        
        #check raw to break up DFAs
        cor_matrix_full <- cor(sem_master_data, use = "pairwise.complete.obs")
        gr<-which(abs(cor_matrix_full["X06_Lmu_IntSprJunH", ]) >= cor_threshold)
        surv<-which(abs(cor_matrix_full["X07_DFA_cpue_IntSprJunHW", ]) >= cor_threshold)
        keycor<-as.data.frame(round(cor_matrix_full[c(gr,surv),c("X06_Lmu_IntSprJunH","X06_DFA_IGF_mu" ,"X07_DFA_cpue_IntSprJunHW")],2)) %>% 
          arrange(desc(abs(X07_DFA_cpue_IntSprJunHW)),desc(abs(X06_Lmu_IntSprJunH)))
        keycor
        
        keycor %>% filter(abs(X06_Lmu_IntSprJunH)>0.8)
        
        round(cor_matrix_full[which(abs(cor_matrix_full["X06_Lmu_IntSprJunH", ]) >= cor_threshold),c("X06_Lmu_IntSprJunH","X07_DFA_cpue_IntSprJunHW")],2)
        growth_raw_proxies <- names(which(abs(cor_matrix_full["X06_Lmu_IntSprJunH", ]) >= cor_threshold))
        abundance_raw_proxies <- names(which(abs(cor_matrix_full["X07_DFA_cpue_IntSprJunHW", ]) >= cor_threshold))
        
#Define model group to analyze=============
        # --- 2. Isolate Your Target Group of 100 Models ---
  #p4my0------------
        p4my0 <- complex_models_only %>%
          filter(Num_Params == 4) %>%
          filter(Missing_Yrs == 0) %>%
          arrange(Adj_AIC2) %>%
          head(100)
        
        # Extract all unique predictors appearing in these 100 models
        unique_predictors <- p4my0 %>%
          separate_rows(Predictors, sep = " \\+ ") %>%
          distinct(Predictors) %>%
          pull(Predictors)
        
        # Ensure the core reference mediators are included in the pool for comparison
        reference_mediators <- c("X06_Lmu_IntSprJunH", "X07_DFA_cpue_IntSprJunHW")
        all_calc_vars <- unique(c(unique_predictors, reference_mediators))
        all_calc_vars <- all_calc_vars[all_calc_vars %in% colnames(guild.dfas1)]
        
        
        # --- 3. Test the "Stand-in" Hypothesis ---
        hypothesis_testing <- p4my0 %>%
          mutate(
            # Split the formula's predictors into a list of individual variables
            Pred_List = str_split(Predictors, " \\+ "),
            
            # Check if this model contains a Growth proxy (or Growth itself)
            Has_Growth_Proxy = map_lgl(Pred_List, ~ any(.x %in% growth_proxies)),
            
            # Check if this model contains an Abundance proxy (or Abundance itself)
            Has_Abundance_Proxy = map_lgl(Pred_List, ~ any(.x %in% abundance_proxies)),
            
            # Does it contain both?
            Meets_Both = Has_Growth_Proxy & Has_Abundance_Proxy
          )
        
        # Calculate support percentage
        percent_supported <- mean(hypothesis_testing$Meets_Both) * 100
        summary(hypothesis_testing[,c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")])
        
        cat("\n=========================================================\n")
        cat(sprintf("HYPOTHESIS TEST: %.1f%% of the Top 100 models contain both\n", percent_supported))
        cat("a Growth proxy and an Abundance proxy standing in for them.\n")
        cat("=========================================================\n\n")
        
        hypothesis_testing_p4my0<-hypothesis_testing
        write_csv(hypothesis_testing_p4my0, "LisaXP/outputs_4/dredge_p4my0_chose_proxies_for_growth_abund.csv")
        
        
#p4my5: Try another model group================  
        model_set="p4my5"
        dat <- complex_models_only %>%
          filter(Num_Params == 4) %>%
          filter(Missing_Yrs == 5) %>%
          arrange(adjaic2) %>%
          head(100)
        
        # Extract all unique predictors appearing in these 100 models
        unique_predictors <- dat %>%
          separate_rows(Predictors, sep = " \\+ ") %>%
          distinct(Predictors) %>%
          pull(Predictors)
        
        # Ensure the core reference mediators are included in the pool for comparison
        reference_mediators <- c("X06_Lmu_IntSprJunH", "X07_DFA_cpue_IntSprJunHW")
        all_calc_vars <- unique(c(unique_predictors, reference_mediators))
        all_calc_vars <- all_calc_vars[all_calc_vars %in% colnames(guild.dfas1)]
        
        
        # --- 3. Test the "Stand-in" Hypothesis ---
        hypothesis_testing <- dat %>%
          mutate(
            # Split the formula's predictors into a list of individual variables
            Pred_List = str_split(Predictors, " \\+ "),
            
            # Check if this model contains a Growth proxy (or Growth itself)
            Has_Growth_Proxy = map_lgl(Pred_List, ~ any(.x %in% growth_proxies)),
            
            # Check if this model contains an Abundance proxy (or Abundance itself)
            Has_Abundance_Proxy = map_lgl(Pred_List, ~ any(.x %in% abundance_proxies)),
            
            # Does it contain both?
            Meets_Both = Has_Growth_Proxy & Has_Abundance_Proxy
          )
        
        # Calculate support percentage
        percent_supported <- mean(hypothesis_testing$Meets_Both) * 100
        summary(hypothesis_testing[,c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")])
        
        cat("\n=========================================================\n")
        cat(sprintf("HYPOTHESIS TEST: %.1f%% of the Top 100 models contain both\n", percent_supported))
        cat("a Growth proxy and an Abundance proxy standing in for them.\n")
        cat("=========================================================\n\n")
        
        hypothesis_testing_p4my5<-hypothesis_testing
        write_csv(hypothesis_testing_p4my0, "LisaXP/outputs_4/dredge_p4my5_chose_proxies_for_growth_abund.csv")
        
        # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
        # Mode :logical    Mode :logical       Mode :logical  
        # FALSE:65         FALSE:56            FALSE:92       
        # TRUE :35         TRUE :44            TRUE :8        
        # 
        
#p3my0: Try another model group================  
        model_set="p3my0"
        dat <- complex_models_only %>%
          filter(Num_Params == 3) %>%
          filter(Missing_Yrs == 0) %>%
          arrange(adjaic2) %>%
          head(100)
        
        # Extract all unique predictors appearing in these 100 models
        unique_predictors <- dat %>%
          separate_rows(Predictors, sep = " \\+ ") %>%
          distinct(Predictors) %>%
          pull(Predictors)
        
        # Ensure the core reference mediators are included in the pool for comparison
        reference_mediators <- c("X06_Lmu_IntSprJunH", "X07_DFA_cpue_IntSprJunHW")
        all_calc_vars <- unique(c(unique_predictors, reference_mediators))
        all_calc_vars <- all_calc_vars[all_calc_vars %in% colnames(guild.dfas1)]
        
        
        # --- 3. Test the "Stand-in" Hypothesis ---
        hypothesis_testing <- dat %>%
          mutate(
            # Split the formula's predictors into a list of individual variables
            Pred_List = str_split(Predictors, " \\+ "),
            
            # Check if this model contains a Growth proxy (or Growth itself)
            Has_Growth_Proxy = map_lgl(Pred_List, ~ any(.x %in% growth_proxies)),
            
            # Check if this model contains an Abundance proxy (or Abundance itself)
            Has_Abundance_Proxy = map_lgl(Pred_List, ~ any(.x %in% abundance_proxies)),
            
            # Does it contain both?
            Meets_Both = Has_Growth_Proxy & Has_Abundance_Proxy
          )
        
        # Calculate support percentage
        percent_supported <- mean(hypothesis_testing$Meets_Both) * 100
        summary(hypothesis_testing[,c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")])
        
        cat("\n=========================================================\n")
        cat(sprintf("HYPOTHESIS TEST: %.1f%% of the Top 100 models contain both\n", percent_supported))
        cat("a Growth proxy and an Abundance proxy standing in for them.\n")
        cat("=========================================================\n\n")
        
        hypothesis_testing_p4my5<-hypothesis_testing
        write_csv(hypothesis_testing_p4my0, "LisaXP/outputs_4/dredge_p4my5_chose_proxies_for_growth_abund.csv")
        
 