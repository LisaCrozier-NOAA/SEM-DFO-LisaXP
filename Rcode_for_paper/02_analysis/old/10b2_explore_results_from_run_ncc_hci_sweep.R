print(
  ncc_sweep_table %>% 
    select(predator, base_beta, abs_minmax_a, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg),
  n = 50
)


#HAKE AND LOONS BEST PREDICTORS OVERALL, AND NEITHER OF THEM WERE IMPROVED BY ADDING TEMPERATURE. SO WE CAN IGNORE HCI FOR NCC SMOLTS IN 24 YR RUN.
  
  print(
    ncc_sweep_table %>% 
 #     filter(Q1_aic_improved=="TRUE") %>% 
      select(predator,base_aic,abs_final_aic, base_beta, abs_minmax_a, abs_weight_w, abs_beta_peff, abs_delta_aic, Q1_aic_improved, Q3_sign_flipped_pos_to_neg),
    n = 50
  )
  predator                                   base_aic abs_final_aic base_beta abs_minmax_a abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved Q3_sign_flipped_pos_to_neg
  <chr>                                         <dbl>         <dbl>     <dbl>        <dbl>        <dbl>         <dbl>         <dbl> <lgl>           <lgl>                     
    1 X09_DFA_HakeAge5Plus_smoltyr                  -7.44        -7.29   -0.555            0.5          0.1        -0.562       0.156   FALSE           FALSE                     
  2 X08_Loons_8_WS                                -4.86        -4.07    0.544            0.5          0.1         0.541       0.794   FALSE           FALSE                     
  3 X08_DFA_DC_corm_3_WS                           2.47        -2.73   -0.216            0            1          -0.400      -5.19    TRUE            FALSE                     
  4 X10_Californian_s_l_2yrLead_WS_smoltyr         3.47        -1.75   -0.0935           0            1          -0.436      -5.22    TRUE            FALSE                     
  5 X09_DFA_ChinAbundSnakeFall_smoltyr             3.31        -0.759  -0.127            0            1          -0.442      -4.07    TRUE            FALSE                     
  6 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr     3.39        -0.461  -0.106            0            1          -0.441      -3.85    TRUE            FALSE                     
  7 X10_Northern_f_s_2yrLead_WS_smoltyr            3.76        -0.365  -0.0180           0            1          -0.397      -4.13    TRUE            FALSE                     
  8 X11_DFA_Harbour_p_WS                           3.77         0.167   0.00190          0            1          -0.420      -3.60    TRUE            TRUE                      
  9 X08_Large_gulls_7_WS                           2.85         0.340   0.203            0            1          -0.406      -2.51    TRUE            TRUE                      
  10 X10_Harbour_s_2yrLead_WS_smoltyr               1.69         1.83    0.307            0.5          0.1         0.307       0.143   FALSE           FALSE                     
  11 X10_AllSeaLionsEMB_smoltyr                     2.11         2.11   -0.266            0.5          0.1        -0.265       0.00497 FALSE           FALSE                     
  12 X10_AllSeaLionsEMB_2025_raw_sc_smoltyr         2.11         2.11   -0.266            0.5          0.1        -0.265       0.00497 FALSE           FALSE                     
  13 X10_DFA_ssl.est.wholerange_2yrLead_adultyr     2.70         2.71   -0.162            0.5          0.1        -0.162       0.0138  FALSE           FALSE                     
  14 X10_Northern_f_s_2yrLead_WS_adultyr            2.73         2.80    0.170            0.5          0.1         0.170       0.0675  FALSE           FALSE                     
  15 X10_Harbour_s_2yrLead_WS_adultyr               3.15         3.10   -0.164            0.3          0.7        -0.189      -0.0478  TRUE            FALSE                     
  16 X10_Californian_s_l_2yrLead_WS_adultyr         3.73         3.25   -0.0364           0            1          -0.119      -0.482   TRUE            FALSE                     
  17 X10_AllSeaLionsEMB_2025_raw_sc_adultyr         3.48         3.48   -0.113            0.5          0.1        -0.112       0.00705 FALSE           FALSE                     
  18 X10_AllSeaLionsEMB_adultyr                     3.48         3.48   -0.113            0.5          0.1        -0.112       0.00705 FALSE           FALSE                     
  > 
    
    #19yrs -- still murre (pos) then hake (neg)
    
  #   predator                                   base_aic abs_final_aic base_beta abs_minmax_a abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved Q3_sign_flipped_pos_to_neg
  # <chr>                                         <dbl>         <dbl>     <dbl>        <dbl>        <dbl>         <dbl>         <dbl> <lgl>           <lgl>                     
  #   1 X08_commonMurre_JSOES                        -9.37         -9.12    0.559            0.5          0.1        0.583        0.254   FALSE           FALSE                     
  # 2 X09_DFA_HakeAge5Plus_smoltyr                 -9.14         -8.94   -0.535            0.5          0.1       -0.542        0.199   FALSE           FALSE                     
  # 3 X08_Loons_8_WS                               -6.40         -5.67    0.567            0.5          0.1        0.571        0.722   FALSE           FALSE                     
  # 4 X10_Northern_f_s_2yrLead_WS_adultyr          -3.08         -2.93    0.370            0.5          0.1        0.377        0.155   FALSE           FALSE                     
  # 5 X10_Californian_s_l_2yrLead_WS_smoltyr        1.81         -1.60   -0.0874           0            1         -0.386       -3.42    TRUE            FALSE                     
  # 6 X08_Large_gulls_7_WS                          0.382        -0.645   0.260            0            1         -0.357       -1.03    TRUE            TRUE                      
  # 7 X10_Northern_f_s_2yrLead_WS_smoltyr           1.43         -0.623   0.153            0            1         -0.342       -2.06    TRUE            TRUE                      
  # 8 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr    2.06          0.125  -0.00117          0            1         -0.353       -1.93    TRUE            FALSE                     
  # 9 X10_Harbour_s_2yrLead_WS_smoltyr              1.73          0.307   0.150            0            1         -0.345       -1.42    TRUE            TRUE                      
  # 10 X08_DFA_DC_corm_3_WS                          1.74          0.491   0.105            0            1         -0.263       -1.25    TRUE            TRUE                      
  # 11 X11_DFA_Harbour_p_WS                          1.29          0.549   0.156            0            1         -0.275       -0.741   TRUE            TRUE                      
  # 12 X10_AllSeaLionsEMB_2025_raw_sc_smoltyr        0.886         0.889  -0.217            0.5          0.1       -0.217        0.00295 FALSE           FALSE                     
  # 13 X10_AllSeaLionsEMB_smoltyr                    0.886         0.889  -0.217            0.5          0.1       -0.217        0.00295 FALSE           FALSE                     
  # 14 X09_DFA_ChinAbundSnakeFall_smoltyr            1.99          1.32    0.0405           0            1         -0.191       -0.666   TRUE            TRUE                      
  # 15 X10_Harbour_s_2yrLead_WS_adultyr              1.59          1.54   -0.144            0.4          0.6       -0.177       -0.0477  TRUE            FALSE                     
  # 16 X10_DFA_ssl.est.wholerange_2yrLead_adultyr    1.63          1.64   -0.116            0.5          0.1       -0.115        0.00995 FALSE           FALSE                     
  # 17 X10_Californian_s_l_2yrLead_WS_adultyr        2.00          1.94    0.0440           0            1         -0.0570      -0.0637  TRUE            TRUE                      
  # 18 X10_AllSeaLionsEMB_2025_raw_sc_adultyr        1.96          1.97   -0.0635           0.5          0.1       -0.0622       0.00409 FALSE           FALSE                     
  # 19 X10_AllSeaLionsEMB_adultyr                    1.96          1.97   -0.0635           0.5          0.1       -0.0622       0.00409 FALSE           FALSE                     
  # > 
    
    
    #ALTPREY 19YRS
  #   predator                                   base_aic abs_final_aic    base_r2 abs_final_r2 base_beta abs_prey           abs_weight_w abs_beta_peff abs_delta_aic Q1_aic_improved Q3_sign_flipped_pos_…¹
  # 1 X08_Loons_8_WS                               -6.40       -13.1    0.359            0.549    0.567   X04_marketsquid_G…          0.9         0.627        -6.66  TRUE            FALSE                 
  # 2 X10_Northern_f_s_2yrLead_WS_adultyr          -3.08       -12.5    0.237            0.535    0.370   X05_herring_GAM_a…          1           1.04         -9.42  TRUE            FALSE                 
  # 3 X08_commonMurre_JSOES                        -9.37       -11.5    0.452            0.510    0.559   X05_anchovy_GAM_s…          1           0.622        -2.13  TRUE            FALSE                 
  # 4 X10_Harbour_s_2yrLead_WS_smoltyr              1.73       -11.3    0.0173           0.504    0.150   X09_DFA_HakeAge5P…          1           1.33        -13.0   TRUE            FALSE                 
  # 5 X09_DFA_HakeAge5Plus_smoltyr                 -9.14        -9.60   0.445            0.459   -0.535   X05_eulachon_duri…          0.3        -0.551        -0.464 TRUE            FALSE                 
  # 6 X10_Northern_f_s_2yrLead_WS_smoltyr           1.43        -8.42   0.0325           0.424    0.153   X09_DFA_HakeAge5P…          1           0.686        -9.85  TRUE            FALSE                 
  # 7 X10_Californian_s_l_2yrLead_WS_adultyr        2.00        -5.71   0.00297          0.336    0.0440  X04_marketsquid_G…          1           0.711        -7.71  TRUE            FALSE                 
  # 8 X08_Large_gulls_7_WS                          0.382       -3.23   0.0845           0.243    0.260   X04_marketsquid_G…          1           0.391        -3.61  TRUE            FALSE                 
  # 9 X11_DFA_Harbour_p_WS                          1.29        -2.71   0.0397           0.222    0.156   X05_anchovy_GAM_s…          1           0.959        -4.00  TRUE            FALSE                 
  # 10 X09_DFA_ChinAbundSnakeFall_smoltyr            1.99        -1.14   0.00369          0.155    0.0405  X09_DFA_HakeAge5P…          1           0.497        -3.13  TRUE            FALSE                 
  # 11 X10_Californian_s_l_2yrLead_WS_smoltyr        1.81        -0.828  0.0128           0.141   -0.0874  X09_DFA_HakeAge5P…          1           0.391        -2.64  TRUE            FALSE                 
  # 12 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr    2.06        -0.151  0.00000230       0.110   -0.00117 X09_DFA_HakeAge5P…          1           0.386        -2.21  TRUE            FALSE                 
  # 13 X10_DFA_ssl.est.wholerange_2yrLead_adultyr    1.63        -0.0244 0.0226           0.104   -0.116   X04_marketsquid_G…          1           0.608        -1.65  TRUE            FALSE                 
  # 14 X08_DFA_DC_corm_3_WS                          1.74         0.260  0.0167           0.0904   0.105   X04_marketsquid_G…          1           0.262        -1.48  TRUE            FALSE                 
  # 15 X10_AllSeaLionsEMB_smoltyr                    0.886        0.702  0.0599           0.0690  -0.217   X04_marketsquid_G…          0.6        -0.471        -0.185 TRUE            FALSE                 
  # 16 X10_AllSeaLionsEMB_2025_raw_sc_smoltyr        0.886        0.702  0.0599           0.0690  -0.217   X04_marketsquid_G…          0.6        -0.471        -0.185 TRUE            FALSE                 
  # 17 X10_Harbour_s_2yrLead_WS_adultyr              1.59         0.764  0.0245           0.0659  -0.144   X05_DFA_abundSard…          0.6        -0.481        -0.824 TRUE            FALSE                 
  # 18 X10_AllSeaLionsEMB_adultyr                    1.96         1.76   0.00502          0.0156  -0.0635  X04_marketsquid_G…          0.9        -0.331        -0.204 TRUE            FALSE                 
  # 19 X10_AllSeaLionsEMB_2025_raw_sc_adultyr        1.96         1.76   0.00502          0.0156  -0.0635  X04_marketsquid_G…          0.9        -0.331        -0.204 TRUE            FALSE                 
  # #