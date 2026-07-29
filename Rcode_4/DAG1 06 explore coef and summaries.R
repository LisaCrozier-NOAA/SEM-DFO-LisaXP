
#1. when was IGF a signif predictor of cpue or sar?#   IGF a signif predictor of cpue (A) or sar (B)?

#2. how did seabirds do? were they over important?


#find names in top models
        topmodels<-read.csv("LisaXP/outputs_4/Doug.topmodels100.csv",row.names=1);
        head(topmodels)
        x<-   topmodels %>%
          group_by(model_id) %>%
          slice(1) %>% 
          ungroup() 
        
        x4<-x %>% select(model_id,PreyNCCindNames,PredNCCindNames,PreyAKindNames,PredAKindNames);x4
        # model_id    PreyNCCindNames       PredNCCindNames             PreyAKindNames          PredAKindNames             
        # 1 DAG1A_long  05.ForageFishNCC_DFA1 09.PredFishNCC_b_DFA1       12.ZooPreyAK_0_smoothed 10.PredMammalNCC_1_smoothed
        # 2 DAG1A_short 05.ForageFishNCC_DFA1 08.PredBirdNCC_b_0_smoothed 12.ZooPreyAK_0_smoothed 10.PredMammalNCC_0_smoothed
        # 3 DAG1B_long  05.ForageFishNCC_DFA1 09.PredFishNCC_b_DFA1       12.ZooPreyAK_DFA1       10.PredMammalNCC_1_smoothed
        # 4 DAG1B_short 05.ForageFishNCC_DFA1 08.PredBirdNCC_b_0_smoothed 12.ZooPreyAK_0_smoothed 10.PredMammalNCC_0_smoothed


#1. When was preyAK negative?---------

        model="DAG1A_long"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
         head(modNum_indicators)
         
         modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
         
        #when was PreyAK neg?
        p1Along_top<-      modelcoefs%>% filter(modNum == modNum_top[1])
          p1Along_top
        p1Along<-    modelcoefs%>% 
          filter(est < 0)  
        100*nrow(p1Along)/nrow(modelcoefs)  #3.1%
        p1Along$modNum
        
        unique(modelcoefs$PreyAKindNames)
        
        tapply(modelcoefs$est,list(modelcoefs$PreyAKindNames),mean,na.rm=T)
        tapply(modelcoefs$est,list(modelcoefs$PreyAKindNames),min,na.rm=T)
        tapply(modelcoefs$est,list(modelcoefs$PreyAKindNames),max,na.rm=T)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
    #just top models (only 8 models daic<3)
        #   PreyAKindNames             est_mu est_min est_max
        # 1 12.ZooPreyAK_0_smoothed    -0.341  -0.341  -0.341
        # 2 12.ZooPreyAK_DFA1          -0.370  -0.407  -0.297
        # 3 13.FishPreyAK_b_2_smoothed -0.293  -0.293  -0.293
        # 4 13.FishPreyAK_b_3_smoothed  0.340   0.340   0.340
        # 5 13.FishPreyAK_b_4_smoothed -0.542  -0.542  -0.542
        # 6 13.FishPreyAK_b_DFA1       -0.517  -0.517  -0.517  
        
    #all models
        #   PreyAKindNames               est_mu  est_min est_max
        # 1 12.ZooPreyAK_0_smoothed    -0.418   -0.494   -0.341 
        # 2 12.ZooPreyAK_1_smoothed     0.147    0.0821   0.261 
        # 3 12.ZooPreyAK_2_smoothed    -0.129   -0.348    0.0459
        # 4 12.ZooPreyAK_DFA1          -0.328   -0.453   -0.198 
        # 5 13.FishPreyAK_DFA1         -0.00909 -0.441    0.268 
        # 6 13.FishPreyAK_b_0_smoothed -0.311   -0.452   -0.117 
        # 7 13.FishPreyAK_b_1_smoothed  0.133    0.0252   0.210 
        # 8 13.FishPreyAK_b_2_smoothed -0.287   -0.366   -0.193 
        # 9 13.FishPreyAK_b_3_smoothed  0.191    0.0364   0.340 
        # 10 13.FishPreyAK_b_4_smoothed -0.120   -0.542    0.151 
        # 11 13.FishPreyAK_b_DFA1       -0.300   -0.517   -0.122 
        # 12 14.CompAK_0_smoothed        0.203   -0.00665  0.347 
        # 13 14.CompAK_1_smoothed       -0.247   -0.456   -0.147 
        # 14 14.CompAK_smoothed          0.303    0.0815   0.505        
        
        model="DAG1B_long"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        
        PreyAKindNames    est_mu est_min est_max
        <chr>              <dbl>   <dbl>   <dbl>
          1 12.ZooPreyAK_DFA1 -0.773  -0.773  -0.773

        
        
        model="DAG1C_long"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        
        PreyAKindNames    est_mu est_min est_max
        1 12.ZooPreyAK_1_smoothed  0.392   0.309   0.504
        
        
#short=========
        
        model="DAG1A_short"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
 
        # PreyAKindNames             est_mu est_min est_max
        # 1 12.ZooPreyAK_0_smoothed    -0.342  -0.342  -0.342
        # 2 12.ZooPreyAK_1_smoothed     0.221   0.221   0.221
        # 3 12.ZooPreyAK_DFA1          -0.548  -0.548  -0.548
        # 4 13.FishPreyAK_b_2_smoothed -0.247  -0.247  -0.247
        
        model="DAG1B_short"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        
            #        PreyAKindNames    est_mu est_min est_max
        # 1 12.ZooPreyAK_0_smoothed    -0.316  -0.316  -0.316
        # 2 12.ZooPreyAK_1_smoothed     0.203   0.203   0.203
        # 3 12.ZooPreyAK_DFA1          -0.677  -0.908  -0.291
        # 4 13.FishPreyAK_b_1_smoothed  0.365   0.365   0.365
        # 5 13.FishPreyAK_b_2_smoothed -0.224  -0.224  -0.224
        # 6 13.FishPreyAK_b_4_smoothed -0.395  -0.501  -0.236
        # 7 13.FishPreyAK_b_DFA1       -0.927  -0.927  -0.927

        
        
        model="DAG1C_short"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv")) %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyAK")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>% group_by(PreyAKindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        
        #       PreyAKindNames    est_mu est_min est_max
        # 1 13.FishPreyAK_DFA1   -0.189  -0.189  -0.189
        # 2 14.CompAK_0_smoothed -0.133  -0.133  -0.133
        # 3 14.CompAK_smoothed   -0.131  -0.131  -0.131
        
#Were there indirect effects of NCCpredators and prey on SAR or cpue?-------
        #just 5 models in the daic<3, excluding HCI 
        model="DAG1C_long"  
        #modNum_indicators <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))# %>% mutate(daic=AIC-AIC[1]) %>% filter(daic<3) %>% select(modNum,PreyAKindNames) 
        modNum_indicators0 <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))
        #remove HCI
        head(guilds) #habCompInd 01.ZooPreyNCC_JSOES_1
        modNum_indicators<-modNum_indicators0 %>% filter(PreyNCCindNames!="01.ZooPreyNCC_JSOES_1_smoothed")
        table(modNum_indicators$PreyNCCindNames)
        
        modNum_indicators<-modNum_indicators0 %>% 
              filter(PreyNCCindNames!="01.ZooPreyNCC_JSOES_1_smoothed")%>% 
              mutate(daic=AIC-AIC[1]) %>% 
              filter(daic<3) %>% 
              select(modNum,PreyNCCindNames,PredNCCindNames) 
        
        
#SAR~predNCC        
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PredNCC")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        PredNCC_indirect_mu<-modelcoefs %>% group_by(PredNCCindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        PredNCC_indirect_p<-modelcoefs %>% group_by(PredNCCindNames) %>% summarise(pvalue_mu=mean(pvalue),pvalue_min=min(pvalue),pvalue_max=max(pvalue))
        cbind(PredNCC_indirect_mu,PredNCC_indirect_p[,-1])
#only 1 pred in dAIC<3?   
          #yes, in DAG1C_long, hake important=0.788  
        importance_summary %>% arrange(desc(DAG1C_long)) %>% filter(SEMnode=="PredNCC")
        
        # PredNCCindNames           est_mu    est_min    est_max    pvalue_mu   pvalue_min  pvalue_max
        # 1 09.PredFishNCC_b_DFA1 -0.5650225 -0.6947849 -0.3451518 0.0004908864 1.332268e-15 0.002454428
        
        #SAR ~ Hake always neg, always signif---------
        #indirect effects were better than direct effects, since cpue->sar was NS
        #would that still be true for SR SAR only, since JSOES is esp weak for UC?
        #TBD
        
#SAR ~ preyNCC        
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="logSAR", rhs=="PreyNCC")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        PreyNCC_indirect_mu<-modelcoefs %>% group_by(PreyNCCindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        PreyNCC_indirect_p<-modelcoefs %>% group_by(PreyNCCindNames) %>% summarise(pvalue_mu=mean(pvalue),pvalue_min=min(pvalue),pvalue_max=max(pvalue))
        x<-cbind(PreyNCC_indirect_mu,PreyNCC_indirect_p[,-1]);x[,-1]<-round(x[,-1],3);x

        #SAR~PreyNCC -- always neg, always highly sign, marketsquid_GAM_2025 or  herring_GAM ------     
        #             PreyNCCindNames     est_mu    est_min    est_max    pvalue_mu   pvalue_min   pvalue_max
        # 1     04.CompFishNCC_smoothed -0.8866272 -1.0233334 -0.7018768 3.882453e-06 2.175054e-07 8.990240e-06
        # 2 05.ForageFishNCC_1_smoothed -0.4615349 -0.5652372 -0.3578326 7.767271e-06 2.151571e-06 1.338297e-05
        
        #               PreyNCCindNames est_mu est_min est_max pvalue_mu pvalue_min pvalue_max
        # 1     04.CompFishNCC_smoothed -0.887  -1.023  -0.702         0          0          0
        # 2 05.ForageFishNCC_1_smoothed -0.462  -0.565  -0.358         0          0          0
        
        
        guilds %>% filter(guild=="04.CompFishNCC")
        guilds %>% filter(guild=="05.ForageFishNCC_1")
        
        
#cpue ~ preyNCC -- never signif------
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
          filter(lhs=="Abundance", rhs=="PreyNCC")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        nrow(modNum_indicators)
        head(modelcoefs)
        head(modNum_indicators)
        
        PreyNCC_altprey_mu<-modelcoefs %>% group_by(PreyNCCindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        PreyNCC_altprey_p<-modelcoefs %>% group_by(PreyNCCindNames) %>% summarise(pvalue_mu=mean(pvalue),pvalue_min=min(pvalue),pvalue_max=max(pvalue))
        x<-cbind(PreyNCC_altprey_mu,PreyNCC_altprey_p[,-1]);x[,-1]<-round(x[,-1],3);x

        #               PreyNCCindNames est_mu est_min est_max pvalue_mu pvalue_min pvalue_max
        # 1     04.CompFishNCC_smoothed -0.012  -0.012  -0.012     0.949      0.949      0.949
        # 2 05.ForageFishNCC_1_smoothed  0.038   0.038   0.038     0.814      0.814      0.814
                 
#short===========    
        model="DAG1C_short"  #only 1 model
        modNum_indicators0 <- read.csv(paste0(path, model, "/SEMresultsByClus.csv"))
        modNum_indicators<-modNum_indicators0 %>% 
          filter(PreyNCCindNames!="01.ZooPreyNCC_JSOES_1_smoothed")%>% 
          mutate(daic=AIC-AIC[1]) %>% 
          filter(daic<3) %>% 
          select(modNum,PreyNCCindNames,PredNCCindNames) 
 
        
        #SAR~predNCC        
        modelcoefs <- read.csv(paste0(path, model, "/parameterEstimates.csv")) %>%
       #   filter(lhs=="logSAR", rhs=="PredNCC")  %>%
          filter(modNum %in% modNum_indicators$modNum) %>%
          left_join(modNum_indicators,by="modNum")
        head(modelcoefs)
        head(modNum_indicators)
        
        modelcoefs %>%
             filter(lhs=="logSAR"|lhs=="Abundance", rhs=="PredNCC"|rhs=="PreyNCC")
        # modNum       lhs op     rhs        est         se          z       pvalue   ci.lower   ci.upper         PreyNCCindNames             PredNCCindNames
        # 1 120396 Abundance  ~ PreyNCC  0.1463142 0.21107611   0.693182 0.4881953332 -0.2673874  0.5600157 04.CompFishNCC_smoothed 08.PredBirdNCC_b_0_smoothed
        # 2 120396 Abundance  ~ PredNCC  0.7589601 0.21107611   3.595670 0.0003235572  0.3452586  1.1726617 04.CompFishNCC_smoothed 08.PredBirdNCC_b_0_smoothed
        # 3 120396    logSAR  ~ PreyNCC -0.5938803 0.05868810 -10.119262 0.0000000000 -0.7089069 -0.4788537 04.CompFishNCC_smoothed 08.PredBirdNCC_b_0_smoothed
        # 4 120396    logSAR  ~ PredNCC -0.9899608 0.05251949 -18.849398 0.0000000000 -1.0928971 -0.8870245 04.CompFishNCC_smoothed 08.PredBirdNCC_b_0_smoothed  
        
        #SAR~pred murre***, neg
        #SAR~pred market squid***, neg
        #cpue ~ prey market squid NS
        
        PredNCC_indirect_mu<-modelcoefs %>% group_by(PredNCCindNames) %>% summarise(est_mu=mean(est),est_min=min(est),est_max=max(est))
        PredNCC_indirect_p<-modelcoefs %>% group_by(PredNCCindNames) %>% summarise(pvalue_mu=mean(pvalue),pvalue_min=min(pvalue),pvalue_max=max(pvalue))
        x<-cbind(PredNCC_indirect_mu,PredNCC_indirect_p[,-1]);x[,-1]<-round(x[,-1],3);x
        
        