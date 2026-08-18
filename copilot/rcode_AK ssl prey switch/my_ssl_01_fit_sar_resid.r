#constrain optimizer so drop in forage fish has to be bad for salmon

#results-- too extreme 0/1 prey switching
# Parameter   Value
# 1        k0 -99.255
# 2        kT  19.844
# 3        kF 152.218
# 4 B_capelin   0.311
# 5 B_herring   0.689
# 6        r2   0.332


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

# =========================================================
# Constrained SSL switch fit
# =========================================================

# ---------------------------
# 0) User paths---------
# ---------------------------
      path_ssl_dat <- "copilot/outputs_2/ssl.dat.csv"   # <- update if needed
      out_dir      <- "copilot/outputs_4"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      ssl.dat <- read_csv(path_ssl_dat, show_col_types = FALSE) %>% clean_names() %>%
        select(-salmon_resid)
      df_sharkonly <- read.csv(file.path("copilot/outputs_4/resid_x16_ishark_cpue.csv"))
      guild.dfasAK <- read.csv(file.path("data_Lisa/guild.dfasAK.csv")) %>% clean_names()
      ssl.dat
      
      ssl.dat<-left_join(ssl.dat, guild.dfasAK %>% select(year,
                                                          x13_stka_herr_matbiom,
                                                          x13_egoa_herring),
                join_by(year)
      )
      ssl.dat<-left_join(ssl.dat, df_sharkonly %>% select(year,sar_residuals),join_by(year)) %>%
        rename(sar_ishark_resid=sar_residuals)
                              names(ssl.dat)                            
      # ---------------------------
      # 1) Read + check
      # ---------------------------
    #  ssl_var<-"ssl_model_eric"; ssl.predict<-ssl.dat
    #  ssl_var<-"ssl_west_pup_pred"
    #  ssl_var<-"ssl_seak_pup_pred"
    #  ssl_var<-"ssl_west_np_pred"
    #  ssl_var<-"ssl_east_pup_pred"
      
       ssl_var<-"ssl_east_np_pred" #good, r2=0.332
    #   ssl_var<-"ssl_cent_pup_pred"
      
    #   ssl_var<-"ssl_cent_np_pred"#good for sar -- r2=0.26

          ssl.dat <- ssl.dat  %>% mutate(ssl = .data[[ssl_var]])
      
names(ssl.dat)



      req <- c(
        "year",
#        "x16_sar",
        "sar_ishark_resid",
        "ssl",
        "sst_wgoa_coastwatch_junjulaug",
        "capelin_avg",
        "herr_avg",
        "x13_stka_herr_matbiom",
        "x13_egoa_herring"
      )
      miss <- setdiff(req, names(ssl.dat))
      if (length(miss) > 0) {
        stop("ssl.dat missing required columns: ", paste(miss, collapse = ", "))
      }
      
      d0 <- ssl.dat %>%
        select(all_of(req)) %>%
        filter(if_all(everything(), ~is.finite(.x))) %>%
        arrange(year)
      
      if (nrow(d0) < 10) stop("Too few complete rows to fit model.")

# ---------------------------
# 2) Scale inputs/target for stable optimization-------
# ---------------------------
      z <- function(x) as.numeric(scale(x))
      
      d <- d0 %>%
        mutate(
          y   = z(sar_ishark_resid),
          ssl = z(ssl),
          sst = z(sst_wgoa_coastwatch_junjulaug),
          f1  = z(capelin_avg),
          f2  = z(x13_stka_herr_matbiom)
        )

# ---------------------------
# 3) Constrained parameterization---------
# ---------------------------
            # unconstrained theta = (k0, kT, log_kF, w1, w2)
            # transforms:
            #   kF = exp(log_kF) > 0
            #   a,b from softmax => a>=0,b>=0,a+b=1
            
            theta_to_par <- function(theta) {
              k0 <- theta[1]
              kT <- theta[2]
              kF <- exp(theta[3])
              
              w <- theta[4:5]
              w <- w - max(w)                  # stability
              ew <- exp(w)
              ab <- ew / sum(ew)
              B_capelin <- ab[1]
              B_herring <- ab[2]
              
              c(k0 = k0, kT = kT, kF = kF, B_capelin = B_capelin, B_herring = B_herring)
            }
            
            predict_sslw <- function(theta, d) {
              p <- theta_to_par(theta)
              forage <- p["B_capelin"] * d$f1 + p["B_herring"] * d$f2
              eta <- p["k0"] + p["kT"] * d$sst - p["kF"] * forage
              p_switch <- plogis(eta)
              yhat <- -(d$ssl * p_switch) #when ssl are high, they reduce salmon survival
              list(
                yhat = yhat,
                eta = eta,
                p_switch = p_switch,
                forage = forage,
                par = p
              )
            }
            
            obj_sse <- function(theta, d) {
              pred <- predict_sslw(theta, d)$yhat
              sum((d$y - pred)^2)
            }

# ---------------------------
# 4) Multi-start optimization--------
# ---------------------------
        set.seed(20260801)
        
        # sensible deterministic starts + random starts
        base_starts <- rbind(
          c(0,  1.5, log(2.0),  0,  0),   # a=b=0.5
          c(0,  1.0, log(1.0),  1, -1),   # a>b
          c(0,  1.0, log(1.0), -1,  1),   # b>a
          c(-1, 2.0, log(3.0),  0,  0)
        )
        
        rand_starts <- cbind(
          runif(40, -3, 3),      # k0
          runif(40, -5, 5),      # kT
          runif(40, log(0.05), log(20)),  # log_kF
          runif(40, -3, 3),      # w1
          runif(40, -3, 3)       # w2
        )
        
        starts <- rbind(base_starts, rand_starts)
        
        fits <- vector("list", nrow(starts))
        vals <- rep(Inf, nrow(starts))
        
        for (i in seq_len(nrow(starts))) {
          fit_i <- try(
            optim(
              par = starts[i, ],
              fn = obj_sse,
              d = d,
              method = "BFGS",
              control = list(maxit = 5000, reltol = 1e-10)
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_i, "try-error")) {
            fits[[i]] <- fit_i
            vals[i] <- fit_i$value
          }
        }
        
        if (all(!is.finite(vals))) stop("All optimization starts failed.")
        
        best_i <- which.min(vals)
        best_fit <- fits[[best_i]]
        best_theta <- best_fit$par
        best_par <- theta_to_par(best_theta)
        best_pred <- predict_sslw(best_theta, d)

# ---------------------------
# 5) Fit metrics----------
# ---------------------------
        resid <- d$y - best_pred$yhat
        sse <- sum(resid^2)
        rmse <- sqrt(mean(resid^2))
        r2 <- 1 - sse / sum((d$y - mean(d$y))^2)
        cor_yy <- cor(d$y, best_pred$yhat)
        
        fit_metrics <- tibble(
          metric = c("n", "sse", "rmse", "r2", "cor_y_yhat", "optim_convergence", "optim_value"),
          value = c(nrow(d), sse, rmse, r2, cor_yy, best_fit$convergence, best_fit$value)
        )

# ---------------------------
# 6) Output tables---------
# ---------------------------
          par_tbl <- tibble(
            parameter = names(best_par),
            estimate = as.numeric(best_par)
          )
          
          pred_tbl <- d %>%
            transmute(
              year,
              sar_ishark_resid_raw = d0$sar_ishark_resid,
              sar_ishark_resid_z = y,
              ssl_raw = d0$ssl,
              ssl_z = ssl,
              sst_z = sst,
              forage1_z = f1,
              forage2_z = f2,
              forage_weighted_z = best_pred$forage,
              eta = best_pred$eta,
              p_switch = best_pred$p_switch,
              ssl_weighted_pred_z = best_pred$yhat
            )
          print(par_tbl)
          print(fit_metrics)
          print(pred_tbl)
          
          ssl.predict[,ssl_var]<-pred_tbl %>% pull(ssl_weighted_pred_z)
          
           write_csv(par_tbl, file.path(out_dir, paste0("ssl_isharkresid_",ssl_var,"_best_parameters.csv")))
           write_csv(fit_metrics, file.path(out_dir, paste0("ssl_isharkresid_",ssl_var,"_fit_metrics.csv")))
           write_csv(pred_tbl, file.path(out_dir, paste0("ssl_isharkresid_",ssl_var,"_predictions.csv")))
           write_csv(ssl.predict, file.path(out_dir, paste0("ssl_isharkresid_",ssl_var,"_all_predictions.csv")))
           
# ---------------------------
# 7) Plots-------
# ---------------------------

    # A) observed SSL, weighted SSL, salmon residual over year (z scale)
          pA_df <- pred_tbl %>%
            select(year, ssl_observed_z = ssl_z, ssl_weighted_z = ssl_weighted_pred_z, sar_ishark_resid_z) %>%
            pivot_longer(-year, names_to = "series", values_to = "value")
          
          pA <- ggplot(pA_df, aes(year, value, color = series)) +
            geom_hline(yintercept = 0, linetype = 2, color = "gray55") +
            geom_line(linewidth = 1) +
            geom_point(size = 1.5) +
            theme_bw() +
            labs(
              title = "SSL, Observed/Weighted, Salmon Resid",
              subtitle = sprintf("k0=%.3f, kT=%.3f, kF=%.3f, a=%.3f, b=%.3f (a+b=1)",
                                 best_par["k0"], best_par["kT"], best_par["kF"], best_par["B_capelin"], best_par["B_herring"]),
              x = "Year", y = "z-score", color = NULL
            )
          

# B) forage + temp diagnostics
          diag_df <- pred_tbl %>%
            mutate(
              temp_unweighted = sst_z,
              temp_weighted = best_par["kT"] * sst_z,
              forage1_contrib = best_par["B_capelin"] * forage1_z,
              forage2_contrib = best_par["B_herring"] * forage2_z,
              forage_total = forage_weighted_z,
              forage_kf_term = -best_par["kF"] * forage_weighted_z
            )
          
          pB1_df <- diag_df %>%
            select(year, forage1_z, forage2_z, forage1_contrib, forage2_contrib, forage_total) %>%
            pivot_longer(-year, names_to = "series", values_to = "value")
          
          pB1 <- ggplot(pB1_df, aes(year, value, color = series)) +
            geom_hline(yintercept = 0, linetype = 2, color = "gray55") +
            geom_line(linewidth = 1) +
            geom_point(size = 1.3) +
            theme_bw() +
            labs(title = "Forage inputs and weights (z-scale)",
                 x = "Year", y = "value", color = NULL)
          
          pB2_df <- diag_df %>%
            select(year, temp_unweighted, temp_weighted, forage_kf_term, eta, p_switch) %>%
            pivot_longer(-c(year, p_switch), names_to = "series", values_to = "value")
          
          pB2 <- ggplot(pB2_df, aes(year, value, color = series)) +
            geom_hline(yintercept = 0, linetype = 2, color = "gray55") +
            geom_line(linewidth = 1) +
            geom_point(size = 1.3) +
            theme_bw() +
            labs(title = "Linear predictor components (z-scale inputs)",
                 subtitle = "eta = k0 + kT*sst - kF*forage",
                 x = "Year", y = "component value", color = NULL)
          
          pB3 <- ggplot(diag_df, aes(year, p_switch)) +
            geom_line(linewidth = 1.1, color = "blue3") +
            geom_point(size = 1.5, color = "blue4") +
            theme_bw() +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = "Switch probability over time", x = "Year", y = "p_switch")
          

# C) observed vs predicted scatter
        pC <- ggplot(pred_tbl, aes(sar_ishark_resid_z, ssl_weighted_pred_z)) +
          geom_point(size = 2) +
          geom_abline(slope = 1, intercept = 0, linetype = 2) +
          theme_bw() +
          labs(title = "Observed vs predicted (z-scale)",
               x = "Observed salmon residual (z)",
               y = "Predicted weighted SSL (z)")
        
        #ggsave(file.path(out_dir, "ssl_switch_constrained_obs_vs_pred.png")), pC, width = 6.2, height = 5.2, dpi = 170)

#combined plot----------        
library(patchwork)
        library(gridExtra)
        #make a table to fill in a plot slot
                myoutput<-round(c(best_par,r2=fit_metrics%>% filter(metric=="r2") %>% pull()),3)
                myoutput_df<-data.frame(
                   Parameter =names(myoutput),
                   Value     = sprintf("%.3f", myoutput),  stringsAsFactors = FALSE )
                
                p_table <- wrap_elements(
                  tableGrob(
                    myoutput_df, 
                    rows = NULL, # removes row numbers
                    theme = ttheme_minimal(
                      core = list(fg_params = list(fontsize = 9)),
                      colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
                    )
                  )
                )

      myplots <- (pA + pB1 + pB2 + pB3 + pC + p_table) + 
        plot_layout(ncol = 2, guides = "collect") + 
        plot_annotation(title = paste(ssl_var, "constrained prey switch")) & 
        theme(
          plot.title = element_text(size = 9),      # Shrinks panel titles (e.g., "Linear predictor components...")
          plot.subtitle = element_text(size = 7),   # Shrinks panel subtitles if present
          axis.title = element_text(size = 8),      # Optional: shrinks axis labels
          axis.text = element_text(size = 7)        # Optional: shrinks axis tick text
        )
      
      print(myplots)
      
ggsave(file.path(out_dir, paste0("ssl_isharkresid_",ssl_var,".png")), myplots, width = 6.2, height = 5.2, dpi = 170)

# ---------------------------
# 8) Console summary
# ---------------------------
cat("\nBest constrained parameters:\n")
print(ssl_var)
print(myoutput_df)




