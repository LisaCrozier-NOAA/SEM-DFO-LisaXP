#Lisa slow processing of coef before generating DAG1 102 formal script

#Sort expected vs unexpected behavior-------
#PreyNCC----
#only sardine was negative, and that just means it was tracking herring or that fewer sardines means more food for salmon

x<-top100_indicator_behavior %>% #filter(Times_Selected_In_Top100>20) %>% 
  filter(rhs=="PreyNCC") %>% 
  filter(model_id != "DAG1C_long") %>%
  filter(model_id !="DAG1C_short")
x    

x %>% filter(Pct_Positive>1)
# model_id    rhs     Lisaname                          Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PreyNCC X04_marketsquid_GAM                                      1              1          100            0                 1             100         0.91 
# 2 DAG1A_short PreyNCC X01_DFA_sumPreyOfPrey_planktonJun                       36             36          100            0                36             100         0.931
# 3 DAG1B_long  PreyNCC X04_marketsquid_GAM                                      7              7          100            0                 7             100         0.91 
# 4 DAG1B_short PreyNCC X01_DFA_sumPreyOfPrey_planktonJun                       29             29          100            0                29             100         0.931

#only sardine was negative, and that just means it was tracking herring or that fewer sardines means more food for salmon
x %>% filter(Pct_Negative >1)

# model_id    rhs     Lisaname             Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PreyNCC X05_DFA_abundSardine                       99              0            0          100                99             100        -0.934
# 2 DAG1A_short PreyNCC X05_DFA_abundSardine                       64              0            0          100                64             100        -0.944
# 3 DAG1B_long  PreyNCC X05_DFA_abundSardine                       93              0            0          100                93             100        -0.934
# 4 DAG1B_short PreyNCC X05_DFA_abundSardine                       71              0            0          100                71             100        -0.945



#PredNCC------
#All neg signs on predators were DFAs
#All positive signs on predators were stragglers 

#I ignored 1C bec the double paths to cpue and sar make this summary confusing
x<-top100_indicator_behavior %>% filter(Times_Selected_In_Top100>10) %>% 
  filter(rhs=="PredNCC") %>% 
  filter(model_id != "DAG1C_long") %>%
  filter(model_id !="DAG1C_short")

#All positive signs on predators were stragglers 
x %>% filter(Pct_Positive>50) %>% arrange(Lisaname,desc(Times_Selected_In_Top100))
#         model_id    rhs     Lisaname              Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PredNCC X08_Loons_8_WS                              27             27          100            0                27             100         0.529
# 2 DAG1A_long  PredNCC X11_Harbor_seal_CR                           1              1          100            0                 1             100         0.439
# 3 DAG1A_short PredNCC X08_commonMurre_JSOES                       37             37          100            0                37             100         0.755
# 4 DAG1A_short PredNCC X08_Loons_8_WS                              24             24          100            0                24             100         0.611
# 5 DAG1A_short PredNCC X11_Harbor_seal_CR                           6              6          100            0                 6             100         0.643
# 6 DAG1B_long  PredNCC X08_Loons_8_WS                              27             27          100            0                27             100         0.545
# 7 DAG1B_long  PredNCC X08_Large_gulls_7_WS                         2              2          100            0                 0               0         0.196
# 8 DAG1B_long  PredNCC X11_DFA_Harbour_p_WS                         2              2          100            0                 0               0         0.011
# 9 DAG1B_long  PredNCC X11_Harbor_seal_CR                           2              2          100            0                 0               0         0.094
# 10 DAG1B_short PredNCC X08_commonMurre_JSOES                       43             43          100            0                43             100         0.67 
# 11 DAG1B_short PredNCC X08_Loons_8_WS                              16             16          100            0                16             100         0.599



#All neg signs on predators were DFAs
x %>% filter(Pct_Negative >1)
# model_id    rhs     Lisaname                   Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PredNCC X09_DFA_HakeAge5Plus                             72              0            0          100                72             100        -0.606
# 2 DAG1A_short PredNCC X09_DFA_HakeAge5Plus                             33              0            0          100                33             100        -0.661
# 3 DAG1B_long  PredNCC X09_DFA_HakeAge5Plus                             63              0            0          100                63             100        -0.613
# 4 DAG1B_long  PredNCC X08_DFA_DC_corm_3_WS                              2              0            0          100                 0               0        -0.228
# 5 DAG1B_long  PredNCC X09_DFA_ChinAbundSnakeFall                        2              0            0          100                 0               0        -0.144
# 6 DAG1B_short PredNCC X09_DFA_HakeAge5Plus                             41              0            0          100                41             100        -0.67 


#PreyAK------
#X12_copepodCom_EGoA was positive and selected 51% in DAG1C_long, but all other preyAK were negative!

x<-top100_indicator_behavior %>% #filter(Times_Selected_In_Top100>20) %>% 
  filter(rhs=="PreyAK")
x
x %>% filter(Pct_Positive>50) %>% filter(Times_Selected_In_Top100>10)
x %>% filter(Pct_Negative>50) %>% filter(Times_Selected_In_Top100>20)
# model_id    rhs    Lisaname                    Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PreyAK X12_copepodBiomass_WGoA                           28              0            0          100                28             100        -0.418
# 2 DAG1B_long  PreyAK X12_DFA_biomassEuphShelfSum                       45              0            0          100                45             100        -0.617
# 3 DAG1B_short PreyAK X12_DFA_biomassEuphShelfSum                       36              0            0          100                36             100        -0.689
# 4 DAG1C_short PreyAK X13_DFA_WGOA_DFA_midTrophic                       21              0            0          100                21             100        -0.315
# >

x %>% filter(Pct_Negative >1)

# x %>% filter(Pct_Positive>50) %>% filter(Times_Selected_In_Top100>10)
# # A tibble: 1 × 10
# model_id   rhs    Lisaname            Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
#   1 DAG1C_long PreyAK X12_copepodCom_EGoA                       51             51          100            0                49            96.1         0.351


# >       x %>% filter(Pct_Negative>50) %>% filter(Times_Selected_In_Top100>10)
# # A tibble: 14 × 10
# model_id    rhs    Lisaname                    Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# 1 DAG1A_long  PreyAK X12_copepodBiomass_WGoA                           28              0          0          100                  28           100          -0.418
# 2 DAG1A_long  PreyAK X12_DFA_biomassEuphShelfSum                       14              0          0          100                  13            92.9        -0.376
# 3 DAG1A_short PreyAK X12_DFA_biomassEuphShelfSum                       15              0          0          100                  12            80          -0.442
# 4 DAG1A_short PreyAK X12_copepodBiomass_WGoA                           12              0          0          100                  12           100          -0.342
# 5 DAG1A_short PreyAK X13_hexagram_EAI                                  12              0          0          100                  12           100          -0.247
# 6 DAG1A_short PreyAK X13_DFA_WGOA_DFA_seabirds                         11              3         27.3         72.7                 8            72.7        -0.62 
# 7 DAG1B_long  PreyAK X12_DFA_biomassEuphShelfSum                       45              0          0          100                  45           100          -0.617
# 8 DAG1B_long  PreyAK X12_copepodBiomass_WGoA                           19              0          0          100                  19           100          -0.441
# 9 DAG1B_short PreyAK X12_DFA_biomassEuphShelfSum                       36              0          0          100                  36           100          -0.689
# 10 DAG1B_short PreyAK X12_copepodBiomass_WGoA                           11              0          0          100                  11           100          -0.318
# 11 DAG1C_long  PreyAK X14_pinkSalmon                                    17              1          5.9         94.1                15            88.2        -0.354
# 12 DAG1C_short PreyAK X13_DFA_WGOA_DFA_midTrophic                       21              0          0          100                  21           100          -0.315
# 13 DAG1C_short PreyAK X12_copepodCom_WGoA                               14              0          0          100                  12            85.7        -0.242
# 14 DAG1C_short PreyAK X12_copepodBiomass_WGoA                           13              2         15.4         84.6                12            92.3        -0.085


#PredAK---------
#CR seals was strongly selected and always positive. Pcod and halibut also positive
#WS seals, ssl, ATF and sleeper sharks were negative    

x<-top100_indicator_behavior %>% filter(Times_Selected_In_Top100>10) %>% 
  filter(rhs=="PredAK") 
x    

x %>% filter(Pct_Positive>50) %>% arrange(Lisaname,desc(Times_Selected_In_Top100))
x %>% filter(Pct_Negative >50) %>% arrange(Lisaname,desc(Times_Selected_In_Top100))

#CR seals was strongly selected and always positive
# model_id    rhs    Lisaname                                  Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# <chr>       <chr>  <chr>                                                        <int>          <dbl>        <dbl>        <dbl>             <dbl>           <dbl>         <dbl>
# 1 DAG1C_short PredAK X10_Harbor_seal_CR_2yrLead                                      96             96        100            0                  96           100           0.956
# 2 DAG1A_short PredAK X10_Harbor_seal_CR_2yrLead                                      77             77        100            0                  77           100           0.563
# 3 DAG1B_short PredAK X10_Harbor_seal_CR_2yrLead                                      51             51        100            0                  51           100           0.535
# 4 DAG1A_long  PredAK X10_Harbor_seal_CR_2yrLead                                      19             19        100            0                  13            68.4         0.325
# 5 DAG1B_long  PredAK X10_Harbor_seal_CR_2yrLead                                      11             11        100            0                   7            63.6         0.297
# 6 DAG1B_long  PredAK X15_PacificCodBiomass_predAK                                    12             12        100            0                  10            83.3         0.392
# 7 DAG1C_long  PredAK X15_halibutBiomassAge8plus_2yrLead_predAK                       13             11         84.6         15.4                11            84.6         0.26 
# 8 DAG1B_long  PredAK X15_halibutBiomassAge8plus_2yrLead_predAK                       11             11        100            0                   9            81.8         0.477


x %>% filter(Pct_Negative >50) %>% arrange(Lisaname,desc(Times_Selected_In_Top100))
# # A tibble: 8 × 10
# model_id    rhs    Lisaname                             Times_Selected_In_Top100 Times_Positive Pct_Positive Pct_Negative Times_Significant Pct_Significant Mean_Estimate
# <chr>       <chr>  <chr>                                                   <int>          <dbl>        <dbl>        <dbl>             <dbl>           <dbl>         <dbl>
# 1 DAG1C_long  PredAK X10_DFA_ssl.est.wholerange_2yrLead                         12              4         33.3         66.7                 7            58.3        -0.187
# 2 DAG1A_long  PredAK X10_Harbour_s_2yrLead_WS                                   18              0          0          100                  13            72.2        -0.353
# 3 DAG1B_long  PredAK X10_Harbour_s_2yrLead_WS                                   14              0          0          100                   9            64.3        -0.338
# 4 DAG1C_long  PredAK X15_ArrowtoothFlounderBiomass_predAK                       11              1          9.1         90.9                 9            81.8        -0.727
# 5 DAG1B_short PredAK X15_DFA_sleeperSharkBSAI_predAK                            22              0          0          100                  22           100          -1.14 
# 6 DAG1A_long  PredAK X15_DFA_sleeperSharkBSAI_predAK                            17              0          0          100                  11            64.7        -0.356
# 7 DAG1A_short PredAK X15_DFA_sleeperSharkBSAI_predAK                            15              0          0          100                  15           100          -0.667
# 8 DAG1B_long  PredAK X15_DFA_sleeperSharkBSAI_predAK                            14              0          0          100                  12            85.7        -0.6  

#WS seals, ssl, ATF and sleeper sharks were negative


#Test the cross-correlation hypothesis---------
#preyNCC: as expected, no explanation necessary
#predNCC:              
#murre ~ HCI, joint environmental variable
#loon ~ cpue, coincidence
#predAK:
# X10_Harbor_seal_CR_2yrLead -- spurious crappy time series, too much missing data
# Pcod is positively corr w/ cpue and negatively with our strong hake predictor of cpue. It is closest to sablefish in AK and anchovy in NCC, negatively corr. Probably spurious
# halibut is negatively corr w/ pink salmon, ssl, and pompano, a warm water fish. So this could represent an ecosystem effect with warm water favoring competitors and the cold water predator is just along for the ride
#preyAK -- way too complicated for me to figure out  
#preyAK have neg sign when we thought they should have a positive sign
#they are positively correlated with some potential predators and competitors (X15_salmonSharkGoA_predAK,X15_sablefishRecruitment_predAK,X14_pinkSalmonAsia)
#they are negatively correlated with other predators, so it isn't a slam dunk
#different variables turned up in the dredge, specifically X13_ammod_WAI  which has a -0.93 corr w/ Growth (X06_Lmu_IntSprJunH), so probably that is just filling in for growth                    




# Define your target Alaskan variables to investigate
wrongsign_targets_predNCC <- c(
  "X08_Loons_8_WS",      
  "X08_commonMurre_JSOES" # The strongly selected, always positive straggler predator
)
wrongsign_targets_predAK <- c(
  "X10_Harbor_seal_CR_2yrLead",
  "X15_PacificCodBiomass_predAK" 
  "X15_halibutBiomassAge8plus_2yrLead_predAK",     
)

wrongsign_targets_preyAK <- c(
  "X12_DFA_biomassEuphShelfSum",
  "X12_copepodBiomass_WGoA", 
  "X13_DFA_WGOA_DFA_midTrophic"     
)

# Run the neighbor function for each target
alaskan_neighborhoods <- map_dfr(wrongsign_targets_preyAK, function(target) {
  find_mds_neighbors(target, plot_data_mds, n_neighbors = 5)
})
print(alaskan_neighborhoods %>% filter(Pairwise_Cor_With_Target<=-0.1))

alaskan_neighborhoods <- map_dfr(wrongsign_targets_predAK, function(target) {
  find_mds_neighbors(target, plot_data_mds, n_neighbors = 5)
})

# Print the results to look at the correlations and SEMnodes
print(as_tibble(alaskan_neighborhoods), n = Inf)

#X10_Harbor_seal_CR_2yrLead -- spurious crappy time series, too much missing data


#  Target_Variable                           Lisaname                           SEMnode                   Relationship_To_Target Pairwise_Cor_With_Target Distance SR_SAR_Correlation
# 6 X10_Harbor_seal_CR_2yrLead                gadid_EAI                          Salmon / Core Diagnostics Positive (+)                               0.28    0.03             0.312  
# 7 X10_Harbor_seal_CR_2yrLead                X08_DFA_DC_corm_3_WS               PredNCC                   Positive (+)                               0.55    0.076           -0.168  
# 8 X10_Harbor_seal_CR_2yrLead                X06_Lmu_IntSprMayW                 Growth                    Negative (-)                              -0.73    0.078            0.128  
# 9 X10_Harbor_seal_CR_2yrLead                X09_canaryRockfish                 PredNCC                   Positive (+)                               0.13    0.115           -0.0729 
# 10 X10_Harbor_seal_CR_2yrLead                ammod_EAI                          Salmon / Core Diagnostics Negative (-)                              -0.68    0.124           -0.222  
# 1 X15_halibutBiomassAge8plus_2yrLead_predAK X01_NHLlogSum_win_05               PreyNCC                   Positive (+)                               0.72    0.038           -0.0262 
# 2 X15_halibutBiomassAge8plus_2yrLead_predAK X14_pinkSalmon                     PreyAK                    Negative (-)                              -0.74    0.046           -0.0208 
# 3 X15_halibutBiomassAge8plus_2yrLead_predAK X04_PacificPompano_                PreyNCC                   Negative (-)                              -0.49    0.06            -0.319  
# 4 X15_halibutBiomassAge8plus_2yrLead_predAK X14_pinkSalmonAsia                 PreyAK                    Negative (-)                              -0.72    0.095            0.0739 
# 5 X15_halibutBiomassAge8plus_2yrLead_predAK X10_DFA_ssl.est.wholerange_2yrLead PredAK                    Negative (-)                              -0.9     0.131           -0.0535 
# 11 X15_PacificCodBiomass_predAK              X15_sablefishBiomass_predAK        PredAK                    Negative (-)                              -0.79    0.029           -0.00953
# 12 X15_PacificCodBiomass_predAK              X05_anchovy_GAM                    PreyNCC                   Negative (-)                              -0.85    0.032           -0.211  
# 13 X15_PacificCodBiomass_predAK              X09_DFA_HakeAge5Plus               PredNCC                   Negative (-)                              -0.53    0.083           -0.569  
# 14 X15_PacificCodBiomass_predAK              X07_DFA_cpue_IntSprJunHW           Abundance                 Positive (+)                               0.29    0.135            0.626  
# 15 X15_PacificCodBiomass_predAK              X01_NHLlogSum_win_15               PreyNCC                   Positive (+)                               0.44    0.144            0.145  
# 

alaskan_neighborhoods %>% filter(Target_Variable=="X15_halibutBiomassAge8plus_2yrLead_predAK")
alaskan_neighborhoods %>% filter(Target_Variable=="X15_PacificCodBiomass_predAK")
# X10_Harbor_seal_CR_2yrLead -- spurious crappy time series, too much missing data
# Pcod is positively corr w/ cpue and negatively with our strong hake predictor of cpue. It is closest to sablefish in AK and anchovy in NCC, negatively corr. Probably spurious
# halibut is negatively corr w/ pink salmon, ssl, and pompano, a warm water fish. So this could represent an ecosystem effect with warm water favoring competitors and the cold water predator is just along for the ride


predNCC_neighborhoods <- map_dfr(wrongsign_targets_predNCC, function(target) {
  find_mds_neighbors(target, plot_data_mds_guild.dfas1_X16SARcor, n_neighbors = 5)
})

# Print the results to look at the correlations and SEMnodes
print(as_tibble(predNCC_neighborhoods), n = Inf)
predNCC_neighborhoods %>% filter(Pairwise_Cor_With_Target>0.5)
#murre ~ HCI, joint environmental variable
#loon ~ cpue, coincidence
# Target_Variable                 Lisaname   SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
# 1        X08_Loons_8_WS         X13_hexagram_WAI    PreyAK           Positive (+)                     0.58    0.060       0.3908780
# 2        X08_Loons_8_WS X07_DFA_cpue_IntSprJunHW Abundance           Positive (+)                     0.54    0.135       0.5486211
# 3 X08_commonMurre_JSOES           X01_habCompInd   PreyNCC           Positive (+)                     0.70    0.122       0.6020082

# Target_Variable       Lisaname                   SEMnode   Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
# 1 X08_Loons_8_WS        X01_Sabine_gull_WS         PreyNCC   Positive (+)                               0.36    0.058          0.0361
# 2 X08_Loons_8_WS        X13_hexagram_WAI           PreyAK    Positive (+)                               0.58    0.06           0.391 
# 3 X08_Loons_8_WS        X01_Red_phal_WS            PreyNCC   Negative (-)                              -0.08    0.12           0.0621
# 4 X08_Loons_8_WS        X06_StomFull_May           Growth    Negative (-)                              -0.16    0.128         -0.120 
# 5 X08_Loons_8_WS        X07_DFA_cpue_IntSprJunHW   Abundance Positive (+)                               0.54    0.135          0.549 
# 6 X08_commonMurre_JSOES X04_Sablefish__JSOES       PreyNCC   Negative (-)                              -0.37    0.045          0.188 
# 7 X08_commonMurre_JSOES X12_copepodCom_EGoA        PreyAK    Negative (-)                              -0.57    0.064          0.0764
# 8 X08_commonMurre_JSOES X14_pinkSalmonNorthAmerica PreyAK    Negative (-)                              -0.19    0.091         -0.254 
# 9 X08_commonMurre_JSOES X09_chilipepper            PredNCC   Positive (+)                               0       0.115         -0.124 
# 10 X08_commonMurre_JSOES X01_habCompInd             PreyNCC   Positive (+)                               0.7     0.122          0.602 
# 

# Target_Variable       Lisaname                  SEMnode                   Relationship_To_Target Pairwise_Cor_With_Target Distance SR_SAR_Correlation
# 1 X08_Loons_8_WS        X06_StomFull_May          Growth                    Negative (-)                              -0.16    0.03             -0.0962
# 2 X08_Loons_8_WS        estAbundHerringRecruits   Salmon / Core Diagnostics Negative (-)                              -0.21    0.075            -0.115 
# 3 X08_Loons_8_WS        X01_Sabine_gull_WS        PreyNCC                   Positive (+)                               0.36    0.082             0.0415
# 4 X08_Loons_8_WS        HakeAge5Plus_2025         Salmon / Core Diagnostics Negative (-)                              -0.34    0.088            -0.470 
# 5 X08_Loons_8_WS        jackmackerel_GAM_2025     Salmon / Core Diagnostics Negative (-)                              -0.37    0.095            -0.285 
# 6 X08_commonMurre_JSOES pacificmackerel_GAM_2025  Salmon / Core Diagnostics Negative (-)                              -0.79    0.033            -0.308 
# 7 X08_commonMurre_JSOES iphc20bigSkate            Salmon / Core Diagnostics Negative (-)                              -0.34    0.068            -0.485 
# 8 X08_commonMurre_JSOES darkblotchedRockfish_2025 Salmon / Core Diagnostics Positive (+)                               0.49    0.072             0.307 
# 9 X08_commonMurre_JSOES X04_Sablefish__JSOES      PreyNCC                   Negative (-)                              -0.37    0.073             0.129 
# 10 X08_commonMurre_JSOES X01_NHLlogSum_win_15      PreyNCC                   Positive (+)                               0.41    0.089             0.145 
# 
# murre corr w/ hake DFA (why are these all sep?)

#preyAK neighborhoods-----
print(alaskan_neighborhoods %>% filter(Pairwise_Cor_With_Target<=-0.2))
#Negative
#               Target_Variable                                 Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
# 1 X12_DFA_biomassEuphShelfSum                     X11_DFA_Harbour_p_WS PredNCC           Negative (-)                    -0.38    0.084      -0.2911234
# 3 X12_DFA_biomassEuphShelfSum                       X02_DFA_SeaNettle_ PreyNCC           Negative (-)                    -0.70    0.139      -0.0565221
# 5     X12_copepodBiomass_WGoA                     X01_NHLlogSum_win_25 PreyNCC           Negative (-)                    -0.37    0.116       0.3738297
# 6     X12_copepodBiomass_WGoA              X10_Northern_f_s_2yrLead_WS  PredAK           Negative (-)                    -0.32    0.127       0.3731474
# 7     X12_copepodBiomass_WGoA               X10_Harbor_seal_CR_2yrLead  PredAK           Negative (-)                    -0.52    0.145       0.3265008

# 2 X12_DFA_biomassEuphShelfSum X03_DFA_comMurreDietHerrSard_Yaquina_NCC PreyNCC           Negative (-)                    -0.11    0.133      -0.2211710
# 4     X12_copepodBiomass_WGoA                     X08_DFA_DC_corm_3_WS PredNCC           Negative (-)                    -0.17    0.098      -0.2326698

#Positive
#                Target_Variable                        Lisaname SEMnode Relationship_To_Target Pairwise_Cor_With_Target Distance SAR_Correlation
# 1 X12_DFA_biomassEuphShelfSum          X01_Red_necked_phal_WS PreyNCC           Positive (+)                     0.57    0.171      0.17016070
# 2 X12_DFA_biomassEuphShelfSum       X15_salmonSharkGoA_predAK  PredAK           Positive (+)                     0.55    0.187     -0.18085232
# 3     X12_copepodBiomass_WGoA              X09_canaryRockfish PredNCC           Positive (+)                     0.50    0.115     -0.17556108
# 4 X13_DFA_WGOA_DFA_midTrophic        X15_sharkCatchGoA_predAK  PredAK           Positive (+)                     0.69    0.131      0.02445275
# 5 X13_DFA_WGOA_DFA_midTrophic      X04_CaliforniaMarketSquid_ PreyNCC           Positive (+)                     0.88    0.142     -0.16338131
# 6 X13_DFA_WGOA_DFA_midTrophic             X12_copepodCom_WGoA  PreyAK           Positive (+)                     0.89    0.161     -0.27255596
# 7 X13_DFA_WGOA_DFA_midTrophic X15_sablefishRecruitment_predAK  PredAK           Positive (+)                     0.91    0.169     -0.09123130
# 8 X13_DFA_WGOA_DFA_midTrophic              X14_pinkSalmonAsia  PreyAK           Positive (+)                     0.74    0.230      0.18685792

#X12_DFA_biomassEuphShelfSum ~ X15_salmonSharkGoA_predAK (r=0.55) -- but relati
#preyAK have neg sign when we thought they should have a positive sign
#they are positively correlated with some potential predators and competitors (X15_salmonSharkGoA_predAK,X15_sablefishRecruitment_predAK,X14_pinkSalmonAsia)
#they are negatively correlated with other predators, so it isn't a slam dunk

#Cross-correlation test results-------------
#preyNCC: as expected, no explanation necessary
#predNCC:              
#murre ~ HCI, joint environmental variable
#loon ~ cpue, coincidence
#predAK:
# X10_Harbor_seal_CR_2yrLead -- spurious crappy time series, too much missing data
# Pcod is positively corr w/ cpue and negatively with our strong hake predictor of cpue. It is closest to sablefish in AK and anchovy in NCC, negatively corr. Probably spurious
# halibut is negatively corr w/ pink salmon, ssl, and pompano, a warm water fish. So this could represent an ecosystem effect with warm water favoring competitors and the cold water predator is just along for the ride
#preyAK --  too complicated for me to figure out  
#preyAK have neg sign when we thought they should have a positive sign
#they are positively correlated with some potential predators and competitors (X15_salmonSharkGoA_predAK,X15_sablefishRecruitment_predAK,X14_pinkSalmonAsia)
#they are negatively correlated with other predators, so it isn't a slam dunk
#different variables turned up in the dredge, specifically X13_ammod_WAI  which has a -0.93 corr w/ Growth (X06_Lmu_IntSprJunH), so probably that is just filling in for growth                    


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

