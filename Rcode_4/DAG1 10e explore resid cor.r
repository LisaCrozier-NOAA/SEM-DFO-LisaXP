#Let's focus on the residual from cpue
#Input data----------
          loadings_df_1998_topbio_allwgoa<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/mybio_DFAloadings_topbio_allwgoa.csv")
          loadings_df_1998_topbio_allwgoa
          
          trend_df_1998_topbio_allwgoa<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/mybio_DFAtrends_topbio_allwgoa.csv") %>%
            select(year,trend) %>%
            rename(wgoa_cap.pcod=trend)
          trend_df_1998_herring_egoa<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/mybio_DFAtrends_herring_egoa.csv") %>%
            select(year,trend) %>%
            rename(egoa_herring=trend)
          
          #Data for modeling
          sem_master_data<-read.csv("LisaXP/outputs_4/sem_master_data.csv",row.names=1) %>%
            mutate(year=1998:2021) %>%
            inner_join(trend_df_1998_topbio_allwgoa,by="year") %>%
            inner_join(trend_df_1998_herring_egoa,by="year") %>%
            mutate(sarSR.sc = scale(sarSR),
                   sarSR.sc = scale(sarSR),
                   cpue.resid.SR = X07_DFA_cpue_IntSprJunHW/scale(sarSR)
          names(sem_master_data)
          
          sem_master_data$resid.cpue.SR<- scale(sem_master_data$sarSR) - sem_master_data$X07_DFA_cpue_IntSprJunHW
          #sem_master_data$resid.cpue.herring.SAR<- sem_master_data$sarSR) - sem_master_data$
          
          sem_master_data<-sem_master_data %>%
            mutate(resid.cpue.herring.SAR = X16_SAR - X07_DFA_cpue_IntSprJunHW - egoa_herring)%>%
            mutate(resid.cpue.herring.wgoacap.SAR = X16_SAR - X07_DFA_cpue_IntSprJunHW - egoa_herring - wgoa_cap.pcod) %>%
            mutate(ssl.altprey = )
          
          # The Math of Linear Differences Matches SEM Exactly
          #For two standardized variables with standard deviation $\sigma = 1$:$$\text{Difference} = \text{sarSR} - X07$$This difference directly isolates Alaska-stage net survival performance:
          #Positive value ($\text{sarSR} - X07 > 0$): Adult returns outperformed early Alaska CPUE (higher survival in Alaska/ocean).
          #Negative value ($\text{sarSR} - X07 < 0$): Adult returns underperformed early Alaska CPUE (lower survival in Alaska/ocean).
          #Because SEM and linear models operate on additive linear combinations ($\text{sarSR} = \beta_1 X07 + \beta_2 \text{Indicator} + \epsilon$), 
          #the linear difference ($\text{sarSR} - X07$) preserves the exact covariance structure of the regression model.
          
          plot(resid.cpue.SR~year,data=sem_master_data,type='b')
          lines(egoa_herring~year,data=sem_master_data,col=2)
          lines(wgoa_cap.pcod~year,data=sem_master_data,col=3)
          abline(h=0)

          "BTS.ATF","BTS.SBLF","BTS.ATKA","BTS.APEX.PRED.BIOM","BTS.MOTILE.EPI.BIOM"
          cor_results1<-round(cor(sem_master_data$resid.cpue.SR,sem_master_data,use="pair"),2)
          sorted_cor1 <- sort(round(cor_results1[1, ], 2), decreasing = TRUE)
          sorted_cor1
          
          cor_results2<-round(cor(sem_master_data$resid.cpue.herring.SAR,sem_master_data,use="pair"),2)
          sorted_cor2 <- sort(round(cor_results2[1, ], 2), decreasing = TRUE)
          sorted_cor2[5:10]
          tail(sorted_cor2)
          sorted_cor2

          cor_results3<-round(cor(sem_master_data$resid.cpue.herring.wgoacap.SAR,sem_master_data,use="pair"),2)
          sorted_cor3 <- sort(round(cor_results3[1, ], 2), decreasing = TRUE)
          sorted_cor3[5:10]
          tail(sorted_cor3)
          sorted_cor2
          
          #herring + sablefish aic=55.008
          #herring + wgoa capelin aic=51.543
          
  X12_DFA_biomassEuphShelfSum
                        ##+ #aic 58.751
                        #X15_salmonSharkGoA_predAK
                  
          single_model_text <- glue('
             # Structural Model
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +
              egoa_herring + 
              X15_sablefishBiomass_predAK
          
          ')
          
          single_model_text <- glue('
             # Structural Model
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +
              egoa_herring + 
              wgoa_cap.pcod
          
          ')
          
          single_model_text <- glue('
             # Structural Model
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +
              egoa_herring + 
              wgoa_cap.pcod +
              X10_DFA_ssl.est.wholerange_2yrLead
          
          ')
          
          fit_single <- sem(single_model_text, data = sem_master_data,std.lv = TRUE,missing = "ML") # 'ML' helps if you have some missing years
          sem_output_summary_fxn(fit_single)
          
          DAG1C_long_topmodel<-fit_single
          
#