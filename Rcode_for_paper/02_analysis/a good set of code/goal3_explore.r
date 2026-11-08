
data_out_dir <- "C:/Users/Lisa.Crozier C drive from work/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/data"
tbl_out_dir <- "C:/Users/Lisa.Crozier C drive from work/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/tables"


#We found a set of predators with a neg sign in the basic dag1A and dag1B models
table4_compact<- read_csv(file.path(tbl_out_dir, "Table4_Compact_Indicator_Diagnostics.csv"))


# Goal 1: topvar (daic<3)-------
importance_topvar<-read_csv(file.path(data_out_dir, "importance_topvar_daic3.csv"))
      nccpred<-importance_topvar %>%
        filter(generic_node=="PredNCC") %>%
        pull(LisaName)
      # "X08_Loons_8_WS"        "X09_DFA_HakeAge5Plus"  "X08_commonMurre_JSOES"
      
      akpred<-importance_topvar %>%
        filter(generic_node=="PredAK") %>%
        pull(LisaName)
      
      # [1] "X10_Harbour_s_2yrLead_WS"        "X15_DFA_sleeperSharks"          
      # [3] "X15_sablefishBiomass_predAK"     "X15_PacificCodBiomass_predAK"   
      # [5] "X15_sablefishRecruitment_predAK" "X15_spinyDogfishBSAI_predAK" 

# Goal 2: What were the signs for those indicators?------ 
      write_csv(indicator_analysis, file.path(data_out_dir, "indicator_analysis.csv"))
      
      indicator_signs<- indicator_analysis %>%
        filter(LisaName %in% akpred) %>%
        filter(pvalue <= 0.05) %>%
        select(LisaName,est,pvalue,AIC,Expected_Sign,Match_Status )
      
      indicator_signs %>% group_by(LisaName) %>%
        summarize(param_mean=mean(est),param_min=min(est),param_max=max(est), aic_min=min(AIC)) %>%
        arrange(param_mean)
      
      
      # LisaName                        param_mean param_min param_max aic_min
      # 1 X15_DFA_sleeperSharks               -0.842    -1.39     -0.310    85.8
      # 2 X10_Harbour_s_2yrLead_WS            -0.382    -0.475    -0.311   132. 
      # 3 X15_sablefishRecruitment_predAK     -0.334    -0.466     0.457    87.6
      # 4 X15_PacificCodBiomass_predAK         0.426     0.385     0.438    86.0
      # 5 X15_spinyDogfishBSAI_predAK          0.458     0.405     0.510    87.8
      # 6 X15_sablefishBiomass_predAK          0.744     0.426     0.930    88.6

      
      indicator_signs<- indicator_analysis %>%
        filter(LisaName %in% nccpred) %>%
        filter(pvalue <= 0.05) %>%
        select(LisaName,est,pvalue,AIC,Expected_Sign,Match_Status )
      
      indicator_signs %>% group_by(LisaName) %>%
        summarize(param_mean=mean(est),param_min=min(est),param_max=max(est), aic_min=min(AIC)) %>%
        arrange(param_mean)
      
      
      # LisaName              param_mean param_min param_max aic_min
      # 1 X09_DFA_HakeAge5Plus      -0.628    -0.671    -0.606    85.8
      # 2 X08_Loons_8_WS             0.569     0.529     0.626    88.7
      # 3 X08_commonMurre_JSOES      0.710     0.670     0.758    85.8      
      # 
      
# Goal 3: did the predators w/ positive signgs change sign when you accounted for altprey or sst?-----------
      #Extract SST models------
      
      # 1. Load the Goal 3 Comparison Summaries
      summaries <- read.csv("Rcode_for_paper/Routput_for_paper/data/Goal3_Model_Summaries_Comparison.csv")
      
      unique(summaries$ncc_predator)
      #[1] "x08_dfa_dc_corm_3_ws" "x08_large_gulls_7_ws" "x09_dfa_hakeage5plus" "x11_dfa_harbour_p_ws"
      #I think that means that murre and loons did not meet our criteria -- maybe I didn't want to look at them, not in alt prey list
      
    ncckeep<-  summaries %>%
        filter(ncc_predator %in% tolower(nccpred),delta_aic_in_pair==0) %>%
        arrange(aic)
      #I think that means the pred*sst model almost always (all but 1 altprey pair_id=31) did better than the baseline
      #but I don't know which pred got the SST -- both?
     
    top_thermal_params <- params %>%
      filter(
        model_type == "Thermal Moderation (Pred*SST)",
        pair_id == ncckeep$pair_id[1]
      )
    
    #these need to filtered to remove models where the SST term was not significant
     
      # 2. Filter exclusively for Thermal Moderation models
      thermal_models <- summaries %>%
        filter(model_type == "Thermal Moderation (Pred*SST)") %>%
        arrange(aic)
      
      # Display Top 10 Thermal Models
      print(head(thermal_models, 10))
      
      # 3. Load Parameter Estimates for the Top Thermal Model
      params <- read.csv("Rcode_for_paper/Routput_for_paper/data/Goal3_All_Parameter_Estimates.csv")
      
      top_thermal_params <- params %>%
        filter(
          model_type == "Thermal Moderation (Pred*SST)",
          pair_id == thermal_models$pair_id[1]
        )
      
      print(top_thermal_params)
      