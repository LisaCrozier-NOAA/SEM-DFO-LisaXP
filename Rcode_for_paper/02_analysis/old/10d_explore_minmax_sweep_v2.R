#Explore results from Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA.csv
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
plot_dir     <- file.path(master_out, "plots")
sem_data_all_data  <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_24yr_noNA <- read.csv(file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_data_short_noNA <- sem_data_all_data %>% 
  filter(complete.cases(.))

#final results-------
altprey_results_24yr
altprey_results_19yr

altprey_results_24yr
# A tibble: 6 × 8
# predator                                     selected_prey                                        fitted_weights base_beta beta_peff delta_aic base_r2 final_r2
# <chr>                                        <chr>                                                <chr>              <dbl>     <dbl>     <dbl>   <dbl>    <dbl>
# 1 X11_DFA_Harbour_p_WS                         X01_habCompInd_smoltyr + X05_anchovy_GAM_smoltyr     1; 0.7            -0.349    -1.20     -10.8    0.385    0.608
# 2 X08_DFA_DC_corm_3_WS                         X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr 1; 0.5            -0.118    -0.978    -13.5    0.309    0.607
# 3 X09_DFA_HakeAge5Plus_smoltyr                 X04_marketsquid_GAM_smoltyr + X05_anchovy_GAM_smolt… 0.3; 0.5          -0.605    -1.21      -6.32   0.464    0.588
# 4 X10_Harbour_s_2yrLead_WS_adultyr             X01_habCompInd_adultyr + X09_DFA_HakeAge5Plus_adult… 0.9; 0.8          -0.546    -2.79      -5.86   0.464    0.580
# 5 X09_DFA_ChinAbundSnakeFall_smoltyr           X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr 1; 0.8             0.106    -1.75     -10.1    0.308    0.545
# 6 X15_ArrowtoothFlounderBiomass_predAK_adultyr X13_mid_il_capelin_adultyr                           1                 -0.187    -0.479     -8.75   0.343    0.544

 altprey_results_19yr
# # A tibble: 15 × 8
# predator                                          selected_prey                                  fitted_weights base_beta beta_peff delta_aic base_r2 final_r2
# 1 X08_DFA_DC_corm_3_WS                              X01_habCompInd_smoltyr + X05_herring_GAM_smol… 1; 0.7            -0.179    -1.56    -17.8     0.299    0.726
# 2 X09_DFA_HakeAge5Plus_smoltyr                      X04_marketsquid_GAM_smoltyr + X05_anchovy_GAM… 0.2; 0.5          -0.817    -1.43     -8.83    0.562    0.725
# 3 X15_sablefishBiomass_predAK_smoltyr               X13_pollock_age1plus_smoltyr + X13_sitkaHerri… 0.9; 1             0.150    -7.07    -17.3     0.298    0.718
# 4 X15_salmonSharkBSAI_predAK_smoltyr                X13_pollock_age1plus_smoltyr + X13_sitkaHerri… 0.7; 0.9          -0.229    -2.47    -15.2     0.311    0.690
# 5 X10_Harbour_s_2yrLead_WS_smoltyr                  X01_habCompInd_smoltyr + X04_marketsquid_GAM_… 1; 0.2            -0.513    -1.11    -11.0     0.375    0.650
# 6 X11_DFA_Harbour_p_WS                              X01_habCompInd_smoltyr + X05_DFA_abundSardine… 1; 0.5            -0.373    -1.71    -10.4     0.384    0.643
# 7 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr X21_sst_egoa_junjulaug_adultyr + X13_capelin_… 0.4; 1            -0.377    -0.962    -9.74    0.398    0.640
# 8 X15_ArrowtoothFlounderBiomass_predAK_adultyr      X21_sst_egoa_junjulaug_adultyr + X13_capelin_… 0.2; 1            -0.188    -0.731   -11.5     0.321    0.629
# 9 X15_DFA_sleeperSharks_adultyr                     X14_pinkSalmonNorthAmerica_adultyr + X13_cape… 0.4; 1            -0.412    -0.735    -4.85    0.508    0.619
# 10 X15_DFA_sleeperSharks_smoltyr                     X13_sitkaHerring_EGoA_smoltyr + X13_capelin_W… 1; 1              -0.613    -0.810    -0.438   0.592    0.601
# 11 X10_Harbour_s_2yrLead_WS_adultyr                  X01_habCompInd_adultyr + X05_eulachon_during_… 0.6; 1            -0.578    -0.862    -5.13    0.468    0.594
# 12 X10_Californian_s_l_2yrLead_WS_smoltyr            None (Unattenuated)                            None               0.608     0.608     0       0.580    0.580
# 13 X15_salmonSharkGoA_predAK_adultyr                 X13_pollock_age1plus_adultyr + X13_capelin_WG… 0.9; 1            -0.294    -3.68     -8.30    0.328    0.566
# 14 X09_chilipepper_adultyr                           X05_anchovy_GAM_adultyr + X05_herring_GAM_adu… 1; 0.9            -0.236    -3.23     -8.53    0.307    0.558
# 15 X15_spinyDogfishBSAI_predAK_smoltyr               X13_capelin_WGoA_smoltyr + X12_egoa_krill_smo… 1; 1               0.188    -0.962    -8.10    0.304    0.545


#19yr test---------------
altprey_results<-read_csv(file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_19yr_noNA_dougsmoothed.csv"));head(altprey_results)
altprey_results<-altprey_results %>% filter(!grepl("X09_DFA_ChinAbundSnakeFall_adultyr|X09_DFA_HakeAge5Plus_adultyr|sharkCatchGoA|X15_spinyDogfishGoA_predAK_adultyr",predator)) %>% distinct()
altprey_results_19yr<-altprey_results %>% filter(final_r2 >0.5) %>% select(predator,selected_prey,fitted_weights,base_beta,beta_peff,delta_aic, base_r2, final_r2)
head(altprey_results_19yr)

altprey_results %>% filter(beta_peff>0) #all NONE
altprey_results %>% filter(beta_peff<0,base_beta>0) %>% select(predator,selected_prey,fitted_weights,base_beta,beta_peff,delta_aic, base_r2, final_r2)
altprey_results %>% filter(pvalue >0.05) %>% select(predator,selected_prey,fitted_weights,pvalue,base_beta,beta_peff,delta_aic, base_r2, final_r2)
altprey_results %>% filter(final_r2 >0.5) %>% select(predator,selected_prey,fitted_weights,pvalue,base_beta,beta_peff,delta_aic, base_r2, final_r2)

#top models (final_r2 >0.5)
# predator                                          selected_prey                                                 fitted_weights     pvalue base_beta beta_peff delta_aic base_r2 final_r2
# 1 X08_DFA_DC_corm_3_WS                              X01_habCompInd_smoltyr + X05_herring_GAM_smoltyr              1; 0.7            2.32e-8    -0.179    -1.56    -17.8     0.299    0.726
# 2 X09_DFA_HakeAge5Plus_smoltyr                      X04_marketsquid_GAM_smoltyr + X05_anchovy_GAM_smoltyr         0.2; 0.5          2.59e-8    -0.817    -1.43     -8.83    0.562    0.725
# 3 X15_sablefishBiomass_predAK_smoltyr               X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr  0.9; 1            4.97e-8     0.150    -7.07    -17.3     0.298    0.718
# 4 X15_salmonSharkBSAI_predAK_smoltyr                X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr  0.7; 0.9          4.57e-7    -0.229    -2.47    -15.2     0.311    0.690
# 5 X10_Harbour_s_2yrLead_WS_smoltyr                  X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr          1; 0.2            6.72e-6    -0.513    -1.11    -11.0     0.375    0.650
# 6 X11_DFA_Harbour_p_WS                              X01_habCompInd_smoltyr + X05_DFA_abundSardine_smoltyr         1; 0.5            9.77e-6    -0.373    -1.71    -10.4     0.384    0.643
# 7 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr X21_sst_egoa_junjulaug_adultyr + X13_capelin_WGoA_adultyr     0.4; 1            1.18e-5    -0.377    -0.962    -9.74    0.398    0.640
# 8 X15_ArrowtoothFlounderBiomass_predAK_adultyr      X21_sst_egoa_junjulaug_adultyr + X13_capelin_WGoA_adultyr     0.2; 1            2.13e-5    -0.188    -0.731   -11.5     0.321    0.629
# 9 X15_DFA_sleeperSharks_adultyr                     X14_pinkSalmonNorthAmerica_adultyr + X13_capelin_WGoA_adultyr 0.4; 1            3.45e-5    -0.412    -0.735    -4.85    0.508    0.619
# 10 X15_DFA_sleeperSharks_smoltyr                     X13_sitkaHerring_EGoA_smoltyr + X13_capelin_WGoA_smoltyr      1; 1              8.33e-5    -0.613    -0.810    -0.438   0.592    0.601
# 11 X10_Harbour_s_2yrLead_WS_adultyr                  X01_habCompInd_adultyr + X05_eulachon_during_chinook_adultyr  0.6; 1            1.15e-4    -0.578    -0.862    -5.13    0.468    0.594
# 12 X10_Californian_s_l_2yrLead_WS_smoltyr            None (Unattenuated)                                           None              2.05e-4     0.608     0.608     0       0.580    0.580
# 13 X15_salmonSharkGoA_predAK_adultyr                 X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       0.9; 1            3.71e-4    -0.294    -3.68     -8.30    0.328    0.566
# 14 X09_chilipepper_adultyr                           X05_anchovy_GAM_adultyr + X05_herring_GAM_adultyr             1; 0.9            4.98e-4    -0.236    -3.23     -8.53    0.307    0.558
# 15 X15_spinyDogfishBSAI_predAK_smoltyr               X13_capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr             1; 1              7.87e-4     0.188    -0.962    -8.10    0.304    0.545

#can eliminate insign pred
#24yr test------------

altprey_results<-read_csv(file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA_akferris.csv"));head(altprey_results)
altprey_results<-altprey_results %>% filter(!grepl("X09_DFA_ChinAbundSnakeFall_adultyr|X09_DFA_HakeAge5Plus_adultyr|sharkCatchGoA|X15_spinyDogfishGoA_predAK_adultyr",predator)) %>% distinct()
altprey_results_24yr<-altprey_results %>% filter(final_r2 >0.5) %>% select(predator,selected_prey,fitted_weights,base_beta,beta_peff,delta_aic, base_r2, final_r2)
nrow(altprey_results_24yr)

#should I reduce the predator list?
print(altprey_results %>% select(predator) %>% arrange(predator),n=Inf)
# X08_DFA_DC_corm_3_WS X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr 1; 0.5            -0.118    -0.978 0.0000152       -0.628          -0.268                 62.5 -3.88     9.65     -13.5   0.309    0.607
# X08_DFA_DC_corm_3_WS X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr 1; 0.5            -0.118    -0.978 0.0000152       -0.628          -0.268                 62.5 -3.88     9.65     -13.5   0.309    0.607

altprey_results<-altprey_results %>% filter(!grepl("X09_DFA_ChinAbundSnakeFall_adultyr|X09_DFA_HakeAge5Plus_adultyr|sharkCatchGoA|X15_spinyDogfishGoA_predAK_adultyr",predator)) %>% distinct()
nrow(altprey_results) #33
head(altprey_results)
altprey_results[2:3,]  

altprey_results %>% filter(beta_peff>0)
altprey_results %>% filter(beta_peff<0,base_beta>0) %>% select(predator,selected_prey,fitted_weights,base_beta,beta_peff,delta_aic, base_r2, final_r2)
altprey_results %>% filter(pvalue >0.05) %>% select(predator,selected_prey,fitted_weights,pvalue,base_beta,beta_peff,delta_aic, base_r2, final_r2)
altprey_results %>% filter(final_r2 >0.4) %>% select(predator,selected_prey,fitted_weights,pvalue,base_beta,beta_peff,delta_aic, base_r2, final_r2)

#still true that all positive final coef are unattenuated, so they were positive before also -- just proxies
# predator                               selected_prey       fitted_weights base_beta beta_peff  pvalue std_eff_beta mean_net_impact realized_damping_pct   aic base_aic delta_aic base_r2 final_r2
# 1 X10_Californian_s_l_2yrLead_WS_smoltyr None (Unattenuated) None               0.485     0.485 0.00250        0.442           0.284                    0  2.21     2.21         0   0.493    0.493
# 2 X15_PacificCodBiomass_predAK_smoltyr   None (Unattenuated) None               0.270     0.270 0.108          0.288           0.136                    0  7.50     7.50         0   0.368    0.368
# 3 X10_Northern_f_s_2yrLead_WS_adultyr    None (Unattenuated) None               0.279     0.279 0.120          0.258           0.135                    0  7.65     7.65         0   0.364    0.364
# 4 X10_Californian_s_l_2yrLead_WS_adultyr None (Unattenuated) None               0.260     0.260 0.148          0.237           0.138                    0  7.95     7.95         0   0.356    0.356


#many that were positive are now negative (altho you really should only keep the daic>4 for 2 extra param)
altprey_results %>% filter(beta_peff<0,base_beta>0) %>% select(predator,selected_prey,fitted_weights,base_beta,beta_peff,delta_aic, base_r2, final_r2)
# A tibble: 14 × 8
# predator                                   selected_prey                                                fitted_weights base_beta beta_peff delta_aic base_r2 final_r2
# 1 X09_DFA_ChinAbundSnakeFall_smoltyr         X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr         1; 0.8            0.106     -1.75    -10.1     0.308    0.545
# 2 X15_sablefishBiomass_predAK_smoltyr        X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr 0.8; 1            0.131     -2.53     -6.62    0.315    0.480
# 3 X10_Northern_f_s_2yrLead_WS_smoltyr        X01_habCompInd_smoltyr + X05_anchovy_GAM_smoltyr             1; 0.9            0.0256    -0.728    -5.73    0.301    0.449
# 4 X15_sablefishBiomass_predAK_adultyr        X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr    1; 1              0.228     -1.70     -3.45    0.332    0.422
# 5 X11_ssl_seak_pup_pred                      X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr   1; 1              0.183     -0.366    -3.05    0.343    0.421
# 6 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr         1; 0.8            0.238     -1.61     -2.28    0.346    0.406
# 7 X15_salmonSharkGoA_predAK_smoltyr          X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr   1; 1              0.0565    -1.09     -3.81    0.302    0.404
# 8 X15_spinyDogfishGoA_predAK_smoltyr         X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr   1; 1              0.0226    -0.368    -3.42    0.301    0.393
# 9 X15_PacificCodBiomass_predAK_adultyr       X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr    0.3; 1            0.0389    -0.293    -2.36    0.302    0.367
# 10 X08_Large_gulls_7_WS                       X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr         1; 0.2            0.212     -0.427    -0.955   0.323    0.350
# 11 X10_DFA_ssl.est.wholerange_2yrLead_adultyr X01_habCompInd_adultyr + X04_marketsquid_GAM_adultyr         1; 1              0.132     -1.63     -0.504   0.317    0.331
# 12 X08_Loons_8_WS                             X01_habCompInd_smoltyr + X05_DFA_abundSardine_smoltyr        1; 0.5            0.227     -0.296    -0.158   0.322    0.326
# 13 X15_sablefishRecruitment_predAK_adultyr    X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr    1; 1              0.0774    -0.585    -0.236   0.306    0.313
# 14 X10_AllSeaLionsEMB_adultyr                 X04_marketsquid_GAM_adultyr + X05_herring_GAM_adultyr        1; 1              0.0301    -2.71     -0.334   0.301    0.310
# > 

#if we use pvalue as a constraint, the final r2 totally suck and daic are low -- no advantage to these models
altprey_results %>% filter(pvalue >0.05) %>% select(predator,selected_prey,fitted_weights,pvalue,base_beta,beta_peff,delta_aic, base_r2, final_r2)
# A tibble: 14 × 9
# predator                                          selected_prey                                              fitted_weights pvalue base_beta beta_peff delta_aic base_r2 final_r2
# 1 X15_spinyDogfishGoA_predAK_smoltyr                X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr 1; 1           0.0549    0.0226    -0.368    -3.42    0.301    0.393
# 2 X15_PacificCodBiomass_predAK_smoltyr              None (Unattenuated)                                        None           0.108     0.270      0.270     0       0.368    0.368
# 3 X15_PacificCodBiomass_predAK_adultyr              X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr  0.3; 1         0.111     0.0389    -0.293    -2.36    0.302    0.367
# 4 X15_salmonSharkGoA_predAK_adultyr                 X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr  0.6; 1         0.118    -0.222     -0.717    -1.42    0.326    0.365
# 5 X10_Northern_f_s_2yrLead_WS_adultyr               None (Unattenuated)                                        None           0.120     0.279      0.279     0       0.364    0.364
# 6 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr X13_mid_il_capelin_adultyr                                 1              0.123    -0.138     -0.280    -1.72    0.316    0.363
# 7 X10_Californian_s_l_2yrLead_WS_adultyr            None (Unattenuated)                                        None           0.148     0.260      0.260     0       0.356    0.356
# 8 X08_Large_gulls_7_WS                              X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr       1; 0.2         0.176     0.212     -0.427    -0.955   0.323    0.350
# 9 X15_sablefishRecruitment_predAK_smoltyr           X13_sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr     1; 1           0.197    -0.0674    -0.249    -1.46    0.305    0.346
# 10 X10_DFA_ssl.est.wholerange_2yrLead_adultyr        X01_habCompInd_adultyr + X04_marketsquid_GAM_adultyr       1; 1           0.294     0.132     -1.63     -0.504   0.317    0.331
# 11 X08_Loons_8_WS                                    X01_habCompInd_smoltyr + X05_DFA_abundSardine_smoltyr      1; 0.5         0.336     0.227     -0.296    -0.158   0.322    0.326
# 12 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr X13_pollock_age1plus_smoltyr + X13_mid_il_capelin_smoltyr  0.5; 1         0.347    -0.0802    -0.181    -0.689   0.305    0.325
# 13 X15_sablefishRecruitment_predAK_adultyr           X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr  1; 1           0.510     0.0774    -0.585    -0.236   0.306    0.313
# 14 X10_AllSeaLionsEMB_adultyr                        X04_marketsquid_GAM_adultyr + X05_herring_GAM_adultyr      1; 1           0.552     0.0301    -2.71     -0.334   0.301    0.310

