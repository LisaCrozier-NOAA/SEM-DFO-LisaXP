#Results
write.csv(all_full_pruned_models, file.path(out_dir, "two_stage_fully_pruned_sem_sweep.csv"), row.names = FALSE)
pred_lookup <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv"), row.names = NULL) %>%
  clean_names() %>%
  mutate(across(everything(), tolower))


my_pred_lookup<-pred_lookup %>% 
  select(short_name_lower,region,pred_data_col,altprey1_data_col,altprey2_data_col,altprey3_data_col)

print(my_pred_lookup) 


results<-  top_models_summary %>%
  select(ncc_predator, ak_predator, active_ncc_terms, active_ak_terms, aic, delta_aic, sar_r2, cfi, pvalue) %>%
  filter(delta_aic<10)%>%
  as_tibble()

results$active_ak_terms[1]
print(results %>% arrange(ncc_predator,aic),n=Inf)
print(results %>% arrange(aic),n=Inf) %>% head(5)

print(
  top_params %>%
    select(pair_id, delta_aic, stage = lhs, term, resolved_species, term_type, est, pvalue, effect_direction),
  n = Inf
)

#ncc pred: hake, gulls and harbour porpoise had sign main effects. cormorants etc had 1 model with sign main effects, but the aic on that was really bad (daic=9.45)
#I'm eliminating all models with a positive predator effect (gulls and cormorants) and an interaction effect without a main effect (arrowtooth and x15_dfa_sleepersharks)
#I'm eliminating all models with a positive predator effect (gulls and cormorants) and an interaction effect without a main effect (arrowtooth and x15_dfa_sleepersharks)

#I'm eliminating all models with daic>10
print(my_pred_lookup %>% filter(pred_data_col=="x08_large_gulls_7_ws")) 
print(my_pred_lookup %>% filter(pred_data_col=="x10_dfa_ssl.est.wholerange_2yrlead")) 

results<-  top_models_summary %>%
  filter(delta_aic<10)%>%
  select(ncc_predator, ak_predator, active_ncc_terms, active_ak_terms, aic, delta_aic, sar_r2, cfi, pvalue) %>%
  as_tibble()

print(top_params %>%
  filter(ncc_predator == "x08_large_gulls_7_ws",ak_predator=="x10_dfa_ssl.est.wholerange_2yrlead") ,
  n=Inf)
print(top_params %>%
        filter(ncc_predator == "x08_large_gulls_7_ws",ak_predator=="x10_californian_s_l_2yrlead_ws") ,
      n=Inf)

results %>% head()
    #gulls*hake + ssl.wholerange*hake + ssl.wholerange*herring had sign int

#top  models (aic<=115, daic<2)
    # ncc_predator         ak_predator                          active_ncc_terms                    active_ak_terms                                          aic delta_aic sar_r2   cfi pvalue
    # 1 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yrlead   x08_large_gulls_7_ws, ncc_int_1   x10_dfa_ssl.est.wholerange_2yrlead, ak_int_1, ak_int_2  113.      0     0.636 0.901  0.164
    # 2 x09_dfa_hakeage5plus x15_arrowtoothflounderbiomass_predak x09_dfa_hakeage5plus              ak_int_1, ak_int_2                                      114.      1.30  0.660 1      0.479
    # 3 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws       x08_large_gulls_7_ws, ncc_int_1   x10_californian_s_l_2yrlead_ws, ak_int_1                115.      1.94  0.565 0.942  0.247

#altprey  
    #   short_name_lower              region        pred_data_col                       altprey1_data_col             altprey2_data_col             altprey3_data_col
    # 1 large_gulls_7_ws_2026         ncc         x08_large_gulls_7_ws                  x09_dfa_hakeage5plus          x05_dfa_abundsardine              <NA>
    # 1 ssl.est.wholerange_2yrlead    ncc         x10_dfa_ssl.est.wholerange_2yrlead    x09_dfa_hakeage5plus_2yrlead  x05_dfa_abundsardine_2yrlead    eulachon_during_chinook_2yrlead

#coef
    # lhs                      term                                   est     se     z   pvalue pair_id ncc_predator         ak_predator                       aic   cfi global_p resolved_species term_type effect_direction delta_aic
    # 1 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws                0.472  0.205   2.31 0.0210        14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 x08_large_gulls… Main Pre… Co-variance / B…         0
    # 2 x07_dfa_cpue_intsprjunhw ncc_int_1                          -0.121  0.0355 -3.40 0.000679      14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 x09_dfa_hakeage… Alternat… Apparent Compet…         0
    # 3 x16_sar                  x07_dfa_cpue_intsprjunhw            0.490  0.134   3.67 0.000245      14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 Early CPUE Link  Stage Li… Positive Link            0
    # 4 x16_sar                  x10_dfa_ssl.est.wholerange_2yrlead  0.608  0.261   2.33 0.0197        14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 x10_dfa_ssl.est… Main Pre… Co-variance / B…         0
    # 5 x16_sar                  ak_int_1                           -0.0653 0.0316 -2.07 0.0388        14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 x09_dfa_hakeage… Alternat… Apparent Compet…         0
    # 6 x16_sar                  ak_int_2                           -0.190  0.0651 -2.92 0.00355       14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yr…  113. 0.901    0.164 x05_dfa_abundsa… Alternat… Apparent Compet…         0

#gulls (and cormorants) have a positive main effect in all models. fuck it. I'm not leaving that.
    # lhs                      term                               est     se     z   pvalue pair_id ncc_predator         ak_predator                      aic   cfi global_p resolved_species      term_type effect_direction delta_aic
    # 1 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws            0.472  0.205   2.31 0.0210        15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws  115. 0.942    0.247 x08_large_gulls_7_ws  Main Pre… Co-variance / B…      1.94
    # 2 x07_dfa_cpue_intsprjunhw ncc_int_1                      -0.121  0.0355 -3.40 0.000679      15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws  115. 0.942    0.247 x09_dfa_hakeage5plus  Alternat… Apparent Compet…      1.94
    # 3 x16_sar                  x07_dfa_cpue_intsprjunhw        0.466  0.148   3.14 0.00167       15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws  115. 0.942    0.247 Early CPUE Link       Stage Li… Positive Link         1.94
    # 4 x16_sar                  x10_californian_s_l_2yrlead_ws  0.843  0.278   3.03 0.00242       15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws  115. 0.942    0.247 x10_californian_s_l_… Main Pre… Co-variance / B…      1.94
    # 5 x16_sar                  ak_int_1                       -0.0949 0.0348 -2.73 0.00635       15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws  115. 0.942    0.247 x09_dfa_hakeage5plus… Alternat… Apparent Compet…      1.94

print(top_params %>%
        filter(ncc_predator == "x08_large_gulls_7_ws",term=="x08_large_gulls_7_ws") ,
      n=Inf)
# lhs                      term                   est    se     z  pvalue pair_id ncc_predator         ak_predator                            aic   cfi global_p resolved_species     term_type          effect_direction delta_aic
# <chr>                    <chr>                <dbl> <dbl> <dbl>   <dbl>   <int> <chr>                <chr>                                <dbl> <dbl>    <dbl> <chr>                <chr>              <chr>                <dbl>
#   1 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.472 0.205  2.31 0.0210       14 x08_large_gulls_7_ws x10_dfa_ssl.est.wholerange_2yrlead    113. 0.901    0.164 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      0   
# 2 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.472 0.205  2.31 0.0210       15 x08_large_gulls_7_ws x10_californian_s_l_2yrlead_ws        115. 0.942    0.247 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      1.94
# 3 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      18 x08_large_gulls_7_ws x15_arrowtoothflounderbiomass_predak  116. 1        0.538 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      3.02
# 4 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      26 x08_large_gulls_7_ws x11_ssl_seak_pup_pred                 120. 0.980    0.359 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      7.20
# 5 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      22 x08_large_gulls_7_ws x15_spinydogfishgoa_predak            121. 0.930    0.230 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      8.19
# 6 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      19 x08_large_gulls_7_ws x15_sablefishbiomass_predak           123. 0.913    0.189 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…      9.89
# 7 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      24 x08_large_gulls_7_ws x15_dfa_sleepersharks                 123. 1        0.517 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…     10.2 
# 8 x07_dfa_cpue_intsprjunhw x08_large_gulls_7_ws 0.516 0.185  2.79 0.00528      21 x08_large_gulls_7_ws x15_pacificcodbiomass_predak          126. 0.917    0.187 x08_large_gulls_7_ws Main Predator Eff… Co-variance / B…     13.6 

print(results %>%
        filter(ncc_predator == "x08_dfa_dc_corm_3_ws") ,
      n=Inf)

print(top_params %>%
        filter(ncc_predator == "x08_dfa_dc_corm_3_ws",term=="x08_dfa_dc_corm_3_ws") ,
      n=Inf)

#in the next best model, the main predator effect was eliminated, despite maintaining interaction effects (arrowtooth*pollock + arrowtooth*herring), 
print(top_params %>%
        filter(ncc_predator == "x09_dfa_hakeage5plus",ak_predator=="x15_arrowtoothflounderbiomass_predak") ,
      n=Inf)

    # lhs                      term                         est     se     z      pvalue pair_id ncc_predator         ak_predator                            aic   cfi global_p resolved_species   term_type effect_direction delta_aic
    # 1 x07_dfa_cpue_intsprjunhw x09_dfa_hakeage5plus     -0.603  0.163  -3.71 0.000211         31 x09_dfa_hakeage5plus x15_arrowtoothflounderbiomass_predak  114.     1    0.479 x09_dfa_hakeage5p… Main Pre… Top-Down Predat…      1.30
    # 2 x16_sar                  x07_dfa_cpue_intsprjunhw  0.373  0.126   2.96 0.00307          31 x09_dfa_hakeage5plus x15_arrowtoothflounderbiomass_predak  114.     1    0.479 Early CPUE Link    Stage Li… Positive Link         1.30
    # 3 x16_sar                  ak_int_1                 -0.265  0.0535 -4.96 0.000000713      31 x09_dfa_hakeage5plus x15_arrowtoothflounderbiomass_predak  114.     1    0.479 x13_pollock_age1p… Alternat… Apparent Compet…      1.30
    # 4 x16_sar                  ak_int_2                  0.0564 0.0263  2.14 0.0322           31 x09_dfa_hakeage5plus x15_arrowtoothflounderbiomass_predak  114.     1    0.479 x13_stka_herr_mat… Alternat… Buffering / Pre…      1.30
x<-unique(top_params %>%
        filter(term_type=="Main Predator Effect",lhs=="x16_sar") %>%
        pull(term));sort(x)
#these terms have main effects in ak stage:
# "x10_californian_s_l_2yrlead_ws"     "x10_dfa_ssl.est.wholerange_2yrlead" "x10_northern_f_s_2yrlead_ws"        "x11_ssl_seak_pup_pred"              "x15_sablefishbiomass_predak"        "x15_spinydogfishgoa_predak"        
#no sleeper sharks!!!!

results[1,]
results %>% filter(ak_predator %in% x)
results %>% filter(!ak_predator %in% x)

#we are down to model 5, but that is my hake + ssl model, so happy about that. 
print(top_params %>%
        filter(ncc_predator == "x09_dfa_hakeage5plus",ak_predator=="x11_ssl_seak_pup_pred") ,
      n=Inf)

# lhs                      term                         est     se     z   pvalue pair_id ncc_predator         ak_predator             aic   cfi global_p resolved_species      term_type                effect_direction delta_aic
# 1 x07_dfa_cpue_intsprjunhw x09_dfa_hakeage5plus     -0.603  0.163  -3.71 0.000211      39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 x09_dfa_hakeage5plus  Main Predator Effect     Top-Down Predat…      5.49
# 2 x16_sar                  x07_dfa_cpue_intsprjunhw  0.390  0.147   2.66 0.00786       39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 Early CPUE Link       Stage Link (CPUE -> SAR) Positive Link         5.49
# 3 x16_sar                  x11_ssl_seak_pup_pred    -0.991  0.298  -3.32 0.000899      39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 x11_ssl_seak_pup_pred Main Predator Effect     Top-Down Predat…      5.49
# 4 x16_sar                  ak_int_1                  0.0625 0.0235  2.65 0.00795       39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 x13_stka_herr_matbiom Alternate Prey Interact… Buffering / Pre…      5.49
# 5 x16_sar                  ak_int_2                  0.0744 0.0280  2.66 0.00786       39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 x13_mid_il_capelin    Alternate Prey Interact… Buffering / Pre…      5.49
# 6 x16_sar                  ak_int_3                  0.0593 0.0251  2.36 0.0182        39 x09_dfa_hakeage5plus x11_ssl_seak_pup_pred  118. 0.919    0.176 x12_egoa.krill        Alternate Prey Interact… Buffering / Pre…      5.49



