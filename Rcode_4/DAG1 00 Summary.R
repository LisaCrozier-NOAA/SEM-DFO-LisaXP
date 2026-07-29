#General summary and progress results to date

#1. Resolved main salmon regression pathways----
    #1a. DAG1A & B resolved that growth does not feed into cpue, and mixed into SAR in DAG1B. DAG1C shows different growth metrics are better for predicting SAR
        source("LisaXP/Rcode_4/DAG1 100 RegressionPathSignif.r")
        signif_summary<-read.csv("LisaXP/outputs_4/MainRegression_Significance_Summary.csv");signif_summary
       #  model_id    path_label           Total_Models_Evaluated Signif_Overall Pct_Signif_Overall Top_100_Evaluated Signif_In_Top_100 Pct_Signif_In_Top_100
       #  1 DAG1A_long  Abundance -> SAR                      46648          46648              100                 100               100                   100
       #  2 DAG1A_short Abundance -> SAR                     188160         177159               94.2               100               100                   100
       #  3 DAG1B_long  Abundance -> SAR                      46648          46631              100                 100               100                   100
       #  4 DAG1B_short Abundance -> SAR                     188160         166435               88.5               100                85                    85
       #  5 DAG1C_long  Abundance -> SAR                      46648          36736               78.8               100                 6                     6
       #  6 DAG1C_short Abundance -> SAR                     188160         117530               62.5               100               100                   100
       #  7 DAG1A_long  Growth -> Abundance                   46648           3332                7.1               100                 1                     1
       #  8 DAG1A_short Growth -> Abundance                  188160          10976                5.8               100                 6                     6
       #  9 DAG1B_long  Growth -> SAR                         46648          21633               46.4               100                43                    43
       #  10 DAG1B_short Growth -> SAR                        188160          52940               28.1               100                54                    54
       #  11 DAG1C_long  Growth -> SAR                         46648          17712               38                 100                97                    97
       #  12 DAG1C_short Growth -> SAR                        188160          49531               26.3               100               100                   100
       #  13 DAG1C_long  PreyNCC -> Abundance                  46648           5880               12.6               100                 0                     0
       #  14 DAG1C_short PreyNCC -> Abundance                 188160          22400               11.9               100                 0                     0
        
        # 
        #Growth -> Abundance sign in only 1% and 6% of models
        #Abundance -> SAR almost always sign, but not in DAG1C_long because predNCC->SAR overwhelmed the signal in DAG1C_long (hake/mackerel)
        #             DAG1B_short was only 85% because something about the new AKpred seems to have shifted the bulk of the prediction over to IGF, leaving Abund small to insignif
        #Growth -> SAR always signif in DAG1C, but weak in DAG1B (43%long, 54%short)
        
    #1b. DAG1C resolved that we do not need to worry about indirect effects of preyNCC and predNCC 
        #summary_table_noHCI<-read.csv("LisaXP/outputs_4/MainRegression_Significance_Summary_noHCI_shortnames.csv",row.names = NULL)
        print(as.tibble(summary_table_noHCI),n=InF)
        
            model_id    path_label           Total_Models_Evaluated Signif_Overall Pct_Signif_Overall Top_N_Evaluated Signif_In_Top_25 Pct_Signif_In_Top_25
        13 DAG1C_long  PreyNCC -> Abundance                  46648           5880               12.6                25                 3                    12
        14 DAG1C_short PreyNCC -> Abundance                 188160          22400               11.9                25                 0                     0

        13 DAG1C_long  PreyNCC -> Abundance                  46648           5880               12.6               100                38                    38
        14 DAG1C_short PreyNCC -> Abundance                 188160          22400               11.9               100                 1                     1
        
        
              Model       Path                  N_All Sig_All Pos_All Neg_All Sig_100 Pos_100 Neg_100 Hypothesis                                      Support    
        
        5 DAG1C_long  Abundance -> SAR      46648    78.8    99.7   0.300       6      99       1 Early marine survival predicts SAR              No (Not Si…
        6 DAG1C_short Abundance -> SAR     188160    62.5    99     1         100     100       0 Early marine survival predicts SAR              Yes        
                                                                                                                                                      
        13 DAG1C_long  PredNCC -> SAR        46648    31.2    42.3  57.7        99       1      99 NCC predators reduce survival beyond JSOES cpue Yes (but r…
        14 DAG1C_short PredNCC -> SAR       188160    30.3    49.3  50.7        99       1      99 NCC predators reduce survival beyond JSOES cpue Yes (but r…
                                                                                                                                                                                                                                                                                                     15        

                
        #why all coef went screwy:
        #When you try to estimate both the direct and indirect pathways simultaneously in a single model (like in DAG1C), 
        #you are forcing your estimators to partition the shared variance between highly collinear steps. 
        #This is exactly why the coefficients became "screwy" — mathematically, 
        #this is called coefficient suppression or multicollinearity-driven sign flipping, and it is a structural artifact, not a fair biological test.
        
        
#2. Resolved best predictors of salmon indicators--------
    #We have good predictors of IGF (preyNCC), but other nodes are mixed pos and neg coef, 
    # As expected: preyNCC
    # Proxy for larger ecosystem condition:
    # Statistical artefact bc of close correlation or noise (stragglers)
    
    #We have good predictors of IGF and cpue, but for predators DFAs were as expected while stragglers were probably spurious (noise) or cross correlated with something else (pos)
    #PreyNCC: only sardine was negative, and that just means it was tracking herring or that fewer sardines means more food for salmon
    #PredNCC: All neg signs on predators were DFAs, all positive signs on predators were stragglers 
    #PredAK:  CR seals was strongly selected and always positive. Pcod and halibut also positive. All of these were stragglers
    #         WS seals, ssl, ATF and sleeper sharks were negative. Those are 2 DFAs and 2 stragglers    
    #PreyAK: X12_copepodCom_EGoA was positive and selected 51% in DAG1C_long, but all other preyAK were negative!
    
    
    
    source("LisaXP/Rcode_4/DAG1 101 AllPathSignif.r")
          #top100_indicator_behavior<-read_csv("LisaXP/outputs_4/Top_100_Indicator_Sign_and_Significance_ABC.csv")
    source("LisaXP/Rcode_4/DAG1 102 AddExplanTable to AllPathSignif.r")
          eco <- ecological_synthesis <-read_csv("LisaXP/outputs_4/Ecological_Synthesis_Summary_Table.csv")
          eco %>% head()
          
          #Successes
          eco %>% select(Node,Node2,Indicator,Type,Hyp,OK,Explanation) %>%
            group_by(Indicator) %>%
            slice(1) %>%
            ungroup() %>%
            filter(OK==1) %>%
            arrange(Node,Hyp)
          
          # Node      Node2      Indicator                            Type      Hyp         OK Explanation                                                          
          # 1 Abundance Abundance  X07_DFA_cpue_IntSprJunHW             DFA       Positive     1 Better NCC survival carries through AK stage                         
          # 2 Growth    Growth     X06_DFA_IGF_mu                       DFA       Positive     1 Bigger is better                                                     
          # 3 Growth    Growth     X06_Lmu_IntSprJunH                   Straggler Positive     1 Bigger is better                                                     
          # 4 PredAK    PredAK     X10_DFA_ssl.est.wholerange_2yrLead   DFA       Negative     1 Direct predator (deepwater top-down consumer)                        
          # 5 PredAK    PredAK     X10_Harbour_s_2yrLead_WS             Straggler Negative     1 Direct predator                                                      
          # 6 PredAK    PredAK     X15_ArrowtoothFlounderBiomass_predAK Straggler Negative     1 Direct predator (deepwater top-down consumer)                        
          # 7 PredAK    PredAK     X15_DFA_sleeperSharkBSAI_predAK      DFA       Negative     1 Direct predator (deepwater top-down consumer)                        
          # 8 PredNCC   PredNCC    X09_DFA_HakeAge5Plus                 DFA       Negative     1 Direct predator                                                      
          # 9 PreyAK    Competitor X14_pinkSalmon                       Straggler Negative     1 Competitor (more pink salmon = less food for Chinook salmon)         
          # 10 PreyAK    PreyAK     X12_copepodCom_EGoA                  Straggler Positive     1 Direct trophic support (cold-water lipid-rich copepod signal)        
          # 11 PreyNCC   Competitor X05_DFA_abundSardine                 DFA       Negative     1 Competitor (tracking herring / fewer sardines = more food for salmon)
          # 12 PreyNCC   PreyNCC    X01_DFA_sumPreyOfPrey_planktonJun    DFA       Positive     1 Direct trophic support (more food = more growth)                     

          #Discrepancies
          eco %>% select(Node,Node2,Indicator,Type,Hyp,OK,Explanation) %>%
            group_by(Indicator) %>%
            slice(1) %>%
            ungroup() %>%
            filter(OK==0) %>%
            arrange(Node,Hyp)
          
          # Node    Node2   Indicator                                 Type      Hyp         OK Explanation                                                                       
          # 1 PredAK  PredAK  X10_Harbor_seal_CR_2yrLead                Straggler Negative     0 Spurious (Noisy / high missing data; trend-rider on short series)                 
          # 2 PredAK  PredAK  X15_PacificCodBiomass_predAK              Straggler Negative     0 Ecosystem proxy (Pcod had negative response to heat wave, favored other predators)
          # 3 PredAK  PredAK  X15_halibutBiomassAge8plus_2yrLead_predAK Straggler Negative     0 Ecosystem proxy (Halibut declines in warm water favor competitors)                
          # 4 PredNCC PredNCC X08_Loons_8_WS                            Straggler Negative     0 Coincidence (Loon population co-trending with CPUE)                               
          # 5 PredNCC PredNCC X08_commonMurre_JSOES                     Straggler Negative     0 Ecosystem proxy (Murre ~ Habitat Compression Index)                               
          # 6 PreyAK  PreyAK  X12_DFA_biomassEuphShelfSum               DFA       Positive     0 Investigate neighborhood relationships further                                    
          # 7 PreyAK  PreyAK  X12_copepodBiomass_WGoA                   Straggler Positive     0 Investigate neighborhood relationships further                                    
          # 8 PreyAK  PreyAK  X12_copepodCom_WGoA                       Straggler Positive     0 Investigate neighborhood relationships further                                    
          # 9 PreyAK  PreyAK  X13_DFA_WGOA_DFA_midTrophic               DFA       Positive     0 Investigate neighborhood relationships further                                    
          # 10 PreyAK  PreyAK  X13_DFA_WGOA_DFA_seabirds                 DFA       Positive     0 Investigate neighborhood relationships further                                    
          # 11 PreyAK  PreyAK  X13_hexagram_EAI                          Straggler Positive     0 Investigate neighborhood relationships further                                    
      
          
          #Investigate further--------
          eco %>% select(Node,Indicator,Type,Pct_Pos,Pct_Neg,Hyp,OK,Explanation) %>%
            group_by(Indicator) %>%
            slice(1) %>%
            ungroup() %>%
            filter(OK==0) %>%
            arrange(Node) %>%
            filter(grepl("proxy|neighborhood",Explanation))
          