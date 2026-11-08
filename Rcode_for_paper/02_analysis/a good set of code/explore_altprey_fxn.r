#let's figure out how this really works


pred_lookup <- read.csv(file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv"), stringsAsFactors = FALSE)
sem_data_noNA     <- read.csv(file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_orig     <- read.csv("data_Lisa/sem_master_data.csv", stringsAsFactors = FALSE)

plot(sem_data_noNA$year,sem_data_noNA$sitkaHerring_EGoA)

sem_complete_data <- sem_data_noNA %>% filter(complete.cases(.))

# Apply +5.0 domain shift across numeric columns (excluding Year)
sem_plus5 <- sem_complete_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ .x + 5.0))


dat<-sem_plus5 %>% select(Year,X07_DFA_cpue_IntSprJunHW,X16_SAR,X21_sst_egoa_junjulaug_smoltyr, 
                                  X15_spinyDogfishGoA_predAK_smoltyr,X11_ssl_seak_pup_pred,X10_Harbour_s_2yrLead_WS_adultyr,
                                  X13_mid_il_capelin_smoltyr, X12_egoa_krill_smoltyr, X13_sitkaHerring_EGoA_smoltyr,
                                  X09_DFA_HakeAge5Plus , X05_DFA_abundSardine_adultyr)
dat<- dat %>% mutate(SAR_residual=res_df$SAR_residual)
dat
n=ncol(dat)-1
matplot(dat$Year,dat[,-1],type='b',col=1:n,lty=1:n)

matplot(dat$Year,dat[,c("SAR_residual","X11_ssl_seak_pup_pred","X13_mid_il_capelin_smoltyr","X12_egoa_krill_smoltyr")],type='b',col=1:n,lty=1:n)
#==================
# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

logistic_gate <- function(a_norm, gamma) {
  1 / (1 + exp(-gamma * (a_norm - 0.5)))
}

evaluate_logistic_switching <- function(pred_col, candidate_col, df_model) {
  
  P <- df_model[[pred_col]]
  A <- df_model[[candidate_col]]
  A_norm <- min_max_scale(A)
  
  # Grid bounds: w in [0.0, 1.0], gamma in [1, 10]
  w_grid     <- seq(0.1, 1.0, by = 0.1)
  gamma_grid <- c(1, 2, 4, 6, 8, 10)
  
  best_fit <- NULL
  best_aic <- Inf
  best_params <- list(w = 0, gamma = 0, est = NA_real_, pval = NA_real_)
  
  # Fit baseline un-moderated model first
  base_syntax <- paste0(target_sar, " ~ ", target_cpue, " + ", pred_col)
  base_sem    <- tryCatch(sem(base_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  
  if (!is.null(base_sem) && lavInspect(base_sem, "converged")) {
    best_aic <- fitMeasures(base_sem, "aic")
    pe <- parameterEstimates(base_sem) %>% filter(op == "~" & rhs == pred_col)
    if (nrow(pe) > 0) {
      best_params <- list(w = 0, gamma = 0, est = pe$est[1], pval = pe$pvalue[1])
    }
  }
  
  # Grid search over w and gamma
  for (w_val in w_grid) {
    for (g_val in gamma_grid) {
      
      # Construct Effective Predator Impact P_eff
      weight_mod <- 1 - (w_val * logistic_gate(A_norm, g_val))
      df_model$P_eff <- P * weight_mod
      
      sem_syntax <- paste0(target_sar, " ~ ", target_cpue, " + P_eff")
      fit <- tryCatch(sem(sem_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
      
      if (!is.null(fit) && lavInspect(fit, "converged")) {
        pe <- parameterEstimates(fit) %>% filter(op == "~" & rhs == "P_eff")
        
        if (nrow(pe) > 0 && !is.na(pe$pvalue[1]) && pe$est[1] < 0) { # Must be negative predation
          current_aic <- fitMeasures(fit, "aic")
          
          if (current_aic < best_aic) {
            best_aic <- current_aic
            best_fit <- fit
            best_params <- list(w = w_val, gamma = g_val, est = pe$est[1], pval = pe$pvalue[1])
          }
        }
      }
    }
  }
  
  return(list(aic = best_aic, params = best_params))
}
#--------------------------------------
sem_dat<-dat
ssl_capelin_res <- evaluate_logistic_switching("X11_ssl_seak_pup_pred", "X13_mid_il_capelin_smoltyr", sem_dat)
ssl_krill <- evaluate_logistic_switching("X11_ssl_seak_pup_pred", "X12_egoa_krill_smoltyr", sem_dat)

#plot krill raw, max/min scaled, logistic transformed, and final weighting
  prey<-"X12_egoa_krill_smoltyr"
      gg=ssl_krill$params$gamma;gg #10
      ww=ssl_krill$params$w;ww #1
      
      A_norm <- min_max_scale(sem_dat[,prey])
      weight_mod_krill <- 1 - (1 * logistic_gate( min_max_scale(sem_dat[,prey]), 6))
      P_eff_krill <-Peff<- sem_dat$X11_ssl_seak_pup_pred * weight_mod

  plot(sem_dat$Year,sem_dat[,prey],col=1,type='l',ylim=c(0,8),main=prey)
  #    lines(sem_dat$Year,A_norm,col=2,type='l')
  #    lines(sem_dat$Year,logistic_gate(A_norm,gamma=gg),col=4,type='l')
      lines(sem_dat$Year,1 - (ww * logistic_gate(A_norm, gg)),col=3,type='l')

    lines(X11_ssl_seak_pup_pred~Year,data=sem_dat,col=1,lwd=3)
    lines(sem_dat$Year,P_eff,col=4,type='l',lwd=3)

    
    
#plot capelin raw, max/min scaled, logistic transformed, and final weighting
    prey<-"X13_mid_il_capelin_smoltyr"
    gg=ssl_capelin_res$params$gamma;gg #10
    ww=ssl_capelin_res$params$w;ww #1
    
    A_norm <- min_max_scale(sem_dat[,prey])
    weight_mod_capelin<-weight_mod <- 1 - (1 * logistic_gate( min_max_scale(sem_dat[,prey]), 6))
    Peff<- sem_dat$X11_ssl_seak_pup_pred * weight_mod
    
    plot(sem_dat$Year,sem_dat[,prey],col=1,type='l',ylim=c(0,8),main=prey)
    lines(sem_dat$Year,A_norm,col=2,type='l')
    lines(sem_dat$Year,logistic_gate(A_norm,gamma=gg),col=4,type='l')
    lines(sem_dat$Year,1 - (ww * logistic_gate(A_norm, gg)),col=3,type='l')
    
    lines(X11_ssl_seak_pup_pred~Year,data=sem_dat,col=1,lwd=3)
    lines(sem_dat$Year,P_eff,col=4,type='l',lwd=3,lty=2)

    #Both-------
    pred="X11_ssl_seak_pup_pred"
    prey1 = "X13_mid_il_capelin_smoltyr"
    prey2 = "X12_egoa_krill_smoltyr"
    prey3 = "X13_sitkaHerring_EGoA_smoltyr"
    mod_prey1 <- evaluate_logistic_switching(pred, prey1, sem_dat)
    mod_prey2 <- evaluate_logistic_switching(pred, prey2, sem_dat)
    mod_prey3 <- evaluate_logistic_switching(pred, prey3, sem_dat)
    wt_prey1 = 1 - (mod_prey1$params$w * logistic_gate(A_norm, mod_prey1$params$gamma))
    wt_prey2 = 1 - (mod_prey2$params$w * logistic_gate(A_norm, mod_prey2$params$gamma))
    wt_prey3 = 1 - (mod_prey3$params$w * logistic_gate(A_norm, mod_prey3$params$gamma))
    
    plot(sem_dat$Year,sem_dat[,pred],col=1,type='l',ylim=c(0,8),main=paste(pred,"\n","prey1=",prey1,"\n","prey2=",prey2))
    lines(sem_dat$Year,1-Anorm_prey1,col=3)
    lines(sem_dat$Year,wt_prey1,col=3,lty=2)

    
    lines(sem_dat$Year,1-Anorm_prey2,col=4)
    lines(sem_dat$Year,wt_prey2,col=4,lty=2)
    
    lines(sem_dat$Year,sem_dat[,pred]*wt_prey1,col="pink")
    lines(sem_dat$Year,sem_dat[,pred]*wt_prey1*wt_prey2,col=2)
    lines(sem_dat$Year,-sem_dat[,"SAR_residual"]+5,col="purple")

        legend("topleft",legend=c(pred,"Peff",prey1,prey2),col=c(1,2,3,4),lty=1,bty='n')
    
    
#dump prey3 
        plot(sem_dat$Year,sem_dat[,prey3],col=1,type='l',ylim=c(0,8),main=prey3)
        
    lines(sem_dat$Year,1-Anorm_prey2,col=4,type='l')
    lines(sem_dat$Year,1-Anorm_prey3,col=6,type='l')
    lines(sem_dat$Year,wt_prey2,col=4,type='l')
    lines(sem_dat$Year,wt_prey3,col=6,type='l')

    #compare max/min and logistic for all prey
    pred="X11_ssl_seak_pup_pred"
    prey1 = "X13_mid_il_capelin_smoltyr"
    prey2 = "X12_egoa_krill_smoltyr"
    prey3 = "X13_sitkaHerring_EGoA_smoltyr"
    Anorm_prey1=min_max_scale(sem_dat[,prey1])
    Anorm_prey2=min_max_scale(sem_dat[,prey2])
    Anorm_prey3=min_max_scale(sem_dat[,prey3])
    mod_prey1 <- evaluate_logistic_switching(pred, prey1, sem_dat)
    mod_prey2 <- evaluate_logistic_switching(pred, prey2, sem_dat)
    mod_prey3 <- evaluate_logistic_switching(pred, prey3, sem_dat)
    wt_prey1 = 1 - (mod_prey1$params$w * logistic_gate(min_max_scale(sem_dat[,prey1]), mod_prey1$params$gamma))
    wt_prey2 = 1 - (mod_prey2$params$w * logistic_gate(min_max_scale(sem_dat[,prey2]), mod_prey2$params$gamma))
    wt_prey3 = 1 - (mod_prey3$params$w * logistic_gate(min_max_scale(sem_dat[,prey3]), mod_prey3$params$gamma))
    wt_prey1 = 1 - Anorm_prey1
    wt_prey2 = 1 - (mod_prey2$params$w * logistic_gate(min_max_scale(sem_dat[,prey2]), mod_prey2$params$gamma))
    wt_prey3 = 1 - (mod_prey3$params$w * logistic_gate(min_max_scale(sem_dat[,prey3]), mod_prey3$params$gamma))
    
    
    
    
ssl_capelin_res
      # $aic
      # aic 
      # 63.785 
      # 
      # $params
      # $params$w
      # [1] 1
      # 
      # $params$gamma
      # [1] 6
      # 
      # $params$est
      # [1] -0.1875094
      # 
      # $params$pval
      # [1] 0.08883102


ssl_krill
      # $aic
      # aic 
      # 64.938 
      # 
      # $params
      # $params$w
      # [1] 1
      # 
      # $params$gamma
      # [1] 10
      # 
      # $params$est
      # [1] -0.1149279
      # 
      # $params$pval
      # [1] 0.2011828

