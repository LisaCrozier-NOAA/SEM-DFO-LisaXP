# See if SSL index works with the indices I already had
 ssl.weight<-function(k0=0, kT=1.5, kF=2,a,b,
                      ssl=ssl.dat$ssl.model.eric,
                      sst=ssl.dat$sst_wgoa_coastwatch_junjulaug,
                      forage1=ssl.dat$capelin.avg,
                      forage2=ssl.dat$herr.avg)
 {
                forage=a*forage1 + b*forage2
                p_switch <-plogis(k0 + kT * sst - kF * forage)
                ssl.w <- ssl * p_switch 
                return(ssl.w)
               }
# For mainline model going forward:
#   
  # k0 = 0
  # kT = 1.5
  # kF = 2.0


#CLIMATE DATA----
  sst<-ak_yr %>% select(year,contains("sst"),contains("pdo"),contains("enso")) %>%
    filter(year>=1998,year<=2021) %>%
    mutate(across(-year,scale))  %>%
  relocate(sst_wgoa_coastwatch_junjulaug, .after=year)
  head(sst)
  cand_sst=names(sst)[-1];cand_sst
        
  # clim<-clim_tr %>% 
  #         select(year,contains("1985")) %>%
  #         drop_na() %>%
  #         left_join(ak_yr %>% select(year,sst_wgoa_coastwatch_junjulaug),join_by(year)) %>%
  #         mutate(sst.sc=scale(sst_wgoa_coastwatch_junjulaug)) %>%
  #         select(-sst_wgoa_coastwatch_junjulaug)
  #       head(clim)
  #       
  #       dat<-clim
  #       head(dat)
  #       matplot(dat$year,dat[,-1],type='l',col=1:(ncol(dat)-1))
  #       legend("topleft",legend=names(dat)[-1],col=1:(ncol(dat)-1),lty=1,bty='n')
  #       matlines(dat$year,dat$sst.sc,lwd=3,col=5)
  #       matlines(dat$year,dat$egoa_clim_1985_t1,lwd=3)
  #       matlines(dat$year,dat$wgoa_clim_1985_t1,lwd=3,col=3)
  #       
  #       round(cor(dat$sst.sc,dat[,-1]),2)
  #       # egoa_clim_1985_t1 egoa_clim_1985_t2 wgoa_clim_1985_t1 wgoa_clim_1985_t2 sst.sc
  #       # [1,]               0.8             -0.18              0.68              0.55      1
  #       #NOT THAT CLOSE! IS IT A DFA/MARSS PROBLEM?

#SSL DATA------
        path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-NMFS/SEM-DFO-NMFS/Data/NCC/ssl_counts.csv"
        ssl.orig<-read.csv(path) %>% select(-est_se) %>%
                  rename(ssl.count.eric=count,ssl.model.eric=est);head(ssl.orig)
        
        # #DON'T USE (2 yr lag)
        #     ssl.dfa<-sem_yr %>% select(year,x10_dfa_ssl_est_wholerange_2yr_lead) %>%
        #       mutate(ssl.yr=year+2);head(ssl.dfa)
        #     sem2<-sem %>%
        #           left_join(sem,ssl.orig,join_by(year));names(sem2) #doesn't start until 2000?
        
        # path2<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/2026_06_29_SEM_AKPred/allData.rds"
        # allData.sem<-readRDS(path2) %>%
        #   filter(grepl("ssl",shortName));head(allData.sem)
        
   #Look at what Doug is using, pretty steep     
          matplot(ssl.orig$year,scale(ssl.orig[,-1]),type='l',col=1:(ncol(ssl.orig)-1),ylim=c(-2,3))
            #  matlines(ssl.dfa$ssl.yr,scale(ssl.dfa$x10_dfa_ssl_est_wholerange_2yr_lead)+1,lwd=3,col=5)
              mylegend=c("raw data from Eric","chasco model","DFA all mammals")
              legend("topleft",legend=mylegend,col=1:length(mylegend),lty=1,bty='n')
              
        #new ak data      
              ssl.ak<-ak_yr %>% select(year,contains("ssl")) %>%
                relocate(ssl_west_pup_pred,.after=year)
              head(ssl.ak)
                # drop_na() %>%
                # mutate(ak_sc=scale(ssl_west_pup_pred)) %>%
                # mutate(across(-year,scale)) 
                # 
        #merge w/ eric modeled sst whole range
               
             ssl.all<-   left_join(ssl.orig, ssl.ak,join_by(year)) %>%
               filter(year>=1998,year<=2021) %>%
                mutate(across(-year,scale)) 
                head(ssl.all)
             
                dat<-ssl.all 
        matplot(dat$year,dat[,-1],type='l',col=1:(ncol(dat)-1))
              legend("topleft",legend=names(dat)[-1],col=1:(ncol(dat)-1),lty=1,bty='n')
              matlines(dat$year,dat$ssl.model.eric,col=2,lwd=3)
              matlines(dat$year,dat$ssl_west_pup_pred,lwd=3,col=3)
              
              cand_ssl=names(ssl.all)[-1];cand_ssl
              
        
 #forage fish ----
              library(dplyr)
              
              forage <- ak_yr %>% 
                select(year, contains("capelin"), contains("herr")) %>%
                filter(year >= 1998, year <= 2021) %>%
                mutate(across(-year, ~as.numeric(scale(.)))) %>%
                rowwise() %>%
                mutate(
                  herr.avg    = mean(c_across(contains("herr")), na.rm = TRUE),
                  capelin.avg = mean(c_across(contains("capelin")), na.rm = TRUE)
                ) %>%
                ungroup()
              
              head(forage)
              cand_forage=c("capelin.avg","herr.avg");cand_forage
              
              
              dat<-forage 
              matplot(dat$year,dat[,-1],type='l',col=1:(ncol(dat)-1))
              legend("topleft",legend=names(dat)[-1],col=1:(ncol(dat)-1),lty=1,bty='n')
              matlines(dat$year,dat$herr.avg,col=3,lwd=3)
              matlines(dat$year,dat$capelin.avg,lwd=3,col=1)
#get salmon data--------
              salmon<-sem_yr %>% select(year,contains("x07"),contains("x16"));head(salmon)
#merge all data -------              
              ssl.dat<-  left_join(ssl.all,forage %>% select(year,herr.avg,capelin.avg),join_by(year))             
              ssl.dat<-  left_join(ssl.dat,sst ,join_by(year))             
              ssl.dat<-  left_join(ssl.dat,salmon ,join_by(year))  
              ssl.dat<-  ssl.dat %>% 
                  mutate(across(-year, ~as.numeric(scale(.))))   %>%
                mutate(salmon.resid=x16_sar-x07_dfa_cpue_int_spr_jun_hw)
              names(ssl.dat)
              summary(ssl.dat)
              
              write.csv(ssl.dat,"copilot/outputs_2/ssl.dat.csv",row.names = FALSE)
              
              cand_ssl=names(ssl.all)[-1];cand_ssl
              cand_forage=c("capelin.avg","herr.avg");cand_forage
              cand_sst=names(sst)[-1];cand_sst
              
#plot raw assumptions------              
              
              ssl.dat$ssl.w1<-    ssl.weight(ssl = ssl.dat$ssl.count.eric, sst=ssl.dat$sst_wgoa_coastwatch_junjulaug,forage=ssl.dat$capelin.avg) 
              ssl.dat$ssl.w2<-    ssl.weight(ssl = ssl.dat$ssl.model.eric, sst=ssl.dat$sst_wgoa_coastwatch_junjulaug,forage=ssl.dat$capelin.avg) 
              
              plot(ssl.w1~year,data=ssl.dat,type='l')                
                  lines(ssl.count.eric~year,data=ssl.dat,lty=2)                
                  lines(ssl.w2~year,data=ssl.dat,lty=2,col=3)                
                  lines(salmon.resid~year,data=ssl.dat,col=2)                
  #now model--------              
                  suppressPackageStartupMessages({
                    library(tidyverse)
                  })
                  
                  # ------------------------------------------------------------
                  # Requires: ssl.dat in environment with columns:
                  # salmon.resid, ssl.model.eric, sst_wgoa_coastwatch_junjulaug,
                  # capelin.avg, herr.avg
                  # ------------------------------------------------------------
                  
                  req <- c("salmon.resid", "ssl.model.eric", "sst_wgoa_coastwatch_junjulaug", "capelin.avg", "herr.avg")
                  miss <- setdiff(req, names(ssl.dat))
                  if (length(miss) > 0) stop("ssl.dat missing columns: ", paste(miss, collapse = ", "))
                  
                  d <- ssl.dat %>%
                    select(all_of(req)) %>%
                    filter(if_all(everything(), ~is.finite(.x)))
                  
                  if (nrow(d) < 10) stop("Too few complete rows after filtering.")
                  
                  # Optional scaling for stable optimization (recommended)
                  d <- d %>%
                    mutate(
                      y = as.numeric(scale(salmon.resid)),
                      ssl = as.numeric(scale(ssl.model.eric)),
                      sst = as.numeric(scale(sst_wgoa_coastwatch_junjulaug)),
                      f1  = as.numeric(scale(capelin.avg)),
                      f2  = as.numeric(scale(herr.avg))
                    )
                  
                  # Model prediction
                  pred_sslw <- function(par, ssl, sst, f1, f2) {
                    k0 <- par[1]; kT <- par[2]; kF <- par[3]; a <- par[4]; b <- par[5]
                    forage <- a*f1 + b*f2
                    p_switch <- plogis(k0 + kT*sst - kF*forage)
                    ssl * p_switch
                  }
                  
                  # Objective: SSE to salmon residuals
                  obj_sse <- function(par, d) {
                    yhat <- pred_sslw(par, d$ssl, d$sst, d$f1, d$f2)
                    sum((d$y - yhat)^2)
                  }
                  
                  # Multi-start to reduce local-minimum risk
                  set.seed(42)
                  starts <- rbind(
                    c(0, 1.5, 2.0, 0.5, 0.5),              # your current defaults-ish
                    matrix(runif(5*30, -2, 2), ncol = 5)   # random starts
                  )
                  
                  fit_list <- vector("list", nrow(starts))
                  
                  for (i in seq_len(nrow(starts))) {
                    fit_list[[i]] <- optim(
                      par = starts[i, ],
                      fn = obj_sse,
                      d = d,
                      method = "L-BFGS-B",
                      lower = c(-5, -5,  0, -5, -5),   # kF constrained nonnegative
                      upper = c( 5,  5, 10,  5,  5)
                    )
                  }
                  
                  best_idx <- which.min(vapply(fit_list, `[[`, numeric(1), "value"))
                  best <- fit_list[[best_idx]]
                  par_hat <- best$par
                  names(par_hat) <- c("k0", "kT", "kF", "a", "b")
                  
                  # Predictions and fit stats
                  d$yhat <- pred_sslw(par_hat, d$ssl, d$sst, d$f1, d$f2)
                  sse <- sum((d$y - d$yhat)^2)
                  sstot <- sum((d$y - mean(d$y))^2)
                  r2 <- 1 - sse/sstot
                  rmse <- sqrt(mean((d$y - d$yhat)^2))
                  cor_y <- cor(d$y, d$yhat)
                  
                  cat("\nBest parameters:\n")
                  print(par_hat)
                  cat("\nFit metrics:\n")
                  cat("SSE :", round(sse, 4), "\n")
                  cat("RMSE:", round(rmse, 4), "\n")
                  cat("R2  :", round(r2, 4), "\n")
                  cat("Cor :", round(cor_y, 4), "\n")
                  
                  # Save results
                  out_dir <- "copilot/outputs_2"
                  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                  
                  readr::write_csv(
                    tibble(parameter = names(par_hat), estimate = as.numeric(par_hat)),
                    file.path(out_dir, "ssl_weight_optimized_parameters.csv")
                  )
                  
                  readr::write_csv(
                    d %>% transmute(y, yhat),
                    file.path(out_dir, "ssl_weight_observed_vs_predicted.csv")
                  )
                  
                  # Optional quick plot
                  p <- ggplot(d, aes(y, yhat)) +
                    geom_point() +
                    geom_abline(slope = 1, intercept = 0, linetype = 2) +
                    theme_bw() +
                    labs(title = "Optimized SSL weight vs salmon.resid (scaled)",
                         x = "Observed salmon.resid (z)", y = "Predicted ssl.w (z-scale target)")
                  
                  
                  ggsave(file.path(out_dir, "ssl_weight_fit_scatter.png"), p, width = 6, height = 5, dpi = 160)

                  print(p)                  
                  print(par_hat)
#plot predicted ssl weights------
                  suppressPackageStartupMessages({
                    library(tidyverse)
                    library(readr)
                    library(janitor)
                  })
                  
                  # ---- paths ----
                  path_ssl_dat <- "copilot/outputs_2/ssl.dat.csv"  # change if your ssl.dat source is elsewhere
                  path_par     <- "copilot/outputs_2/ssl_weight_optimized_parameters.csv"
                  out_dir      <- "copilot/outputs_2"
                  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                  
                  # ---- read data ----
                  ssl.dat <- read_csv(path_ssl_dat, show_col_types = FALSE) %>% clean_names()
                  par_tbl <- read_csv(path_par, show_col_types = FALSE)
                  
                  # expected columns in ssl.dat (cleaned names)
                  req <- c("year", "ssl_model_eric", "sst_wgoa_coastwatch_junjulaug", "capelin_avg", "herr_avg", "salmon_resid")
                  miss <- setdiff(req, names(ssl.dat))
                  if (length(miss) > 0) stop("ssl.dat missing columns: ", paste(miss, collapse = ", "))
                  
                  # ---- parameters ----
                  # expects columns: parameter, estimate
                  p <- setNames(par_tbl$estimate, par_tbl$parameter)
                  need_par <- c("k0","kT","kF","a","b")
                  if (!all(need_par %in% names(p))) stop("Parameter file missing one of: ", paste(need_par, collapse=", "))
                  
                  k0 <- p[["k0"]]; kT <- p[["kT"]]; kF <- p[["kF"]]; a <- p[["a"]]; b <- p[["b"]]
                  
                  # ---- compute weighted SSL ----
                  plot_df <- ssl.dat %>%
                    mutate(
                      forage = a * capelin_avg + b * herr_avg,
                      p_switch = plogis(k0 + kT * sst_wgoa_coastwatch_junjulaug - kF * forage),
                      ssl_weighted = ssl_model_eric * p_switch
                    ) %>%
                    select(year, ssl_observed = ssl_model_eric, ssl_weighted, salmon_resid) %>%
                    arrange(year)
                  
                  # long format for plotting raw units
                  plot_long <- plot_df %>%
                    pivot_longer(-year, names_to = "series", values_to = "value")
                  
                  # ---- Plot 1: raw units (single panel) ----
                  p1 <- ggplot(plot_long, aes(year, value, color = series)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.6) +
                    theme_bw() +
                    labs(
                      title = "Observed SSL, Weighted SSL, and Salmon Residual by Year",
                      subtitle = sprintf("weighted SSL = ssl * logit^-1(k0 + kT*sst - kF*(a*capelin + b*herring)); k0=%.3f, kT=%.3f, kF=%.3f, a=%.3f, b=%.3f",
                                         k0, kT, kF, a, b),
                      x = "Year", y = "Value", color = NULL
                    )
                  
                  ggsave(file.path(out_dir, "timeseries_ssl_observed_weighted_salmonresid_raw.png"), p1, width = 10, height = 5.8, dpi = 170)
                  
                  # ---- Plot 2: standardized overlay (shape comparison) ----
                  z <- function(x) as.numeric(scale(x))
                  plot_z <- plot_df %>%
                    mutate(
                      ssl_observed_z = z(ssl_observed),
                      ssl_weighted_z = z(ssl_weighted),
                      salmon_resid_z = z(salmon_resid)
                    ) %>%
                    select(year, ssl_observed_z, ssl_weighted_z, salmon_resid_z) %>%
                    pivot_longer(-year, names_to = "series", values_to = "zvalue")
                  
                  p2 <- ggplot(plot_z, aes(year, zvalue, color = series)) +
                    geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.5) +
                    theme_bw() +
                    labs(
                      title = "Standardized Time Series (Shape Comparison)",
                      subtitle = "All series z-scored",
                      x = "Year", y = "Standardized value (z)", color = NULL
                    )
                  
                  ggsave(file.path(out_dir, "timeseries_ssl_observed_weighted_salmonresid_z.png"), p2, width = 10, height = 5.8, dpi = 170)
                  
                  # optional csv output
                  write_csv(plot_df, file.path(out_dir, "timeseries_ssl_observed_weighted_salmonresid.csv"))
                  
                  message("Wrote plots and data to: ", out_dir)
                  
                  print(p1)
                  print(p2)
#plot forage fish and temperature inputs to model---------
                  suppressPackageStartupMessages({
                    library(tidyverse)
                    library(readr)
                    library(janitor)
                  })
                  
                  # ---------------------------
                  # Paths
                  # ---------------------------
                  path_ssl_dat <- "copilot/outputs_2/ssl.dat.csv"   # change if needed
                  path_par     <- "copilot/outputs_2/ssl_weight_optimized_parameters.csv"
                  out_dir      <- "copilot/outputs_2"
                  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                  
                  # ---------------------------
                  # Read inputs
                  # ---------------------------
                  ssl.dat <- read_csv(path_ssl_dat, show_col_types = FALSE) %>% clean_names()
                  par_tbl <- read_csv(path_par, show_col_types = FALSE)
                  
                  req_cols <- c("year", "sst_wgoa_coastwatch_junjulaug", "capelin_avg", "herr_avg")
                  miss <- setdiff(req_cols, names(ssl.dat))
                  if (length(miss) > 0) stop("ssl.dat missing columns: ", paste(miss, collapse = ", "))
                  
                  p <- setNames(par_tbl$estimate, par_tbl$parameter)
                  need_par <- c("k0","kT","kF","a","b")
                  if (!all(need_par %in% names(p))) {
                    stop("Parameter file must contain: ", paste(need_par, collapse = ", "))
                  }
                  
                  k0 <- p[["k0"]]; kT <- p[["kT"]]; kF <- p[["kF"]]; a <- p[["a"]]; b <- p[["b"]]
                  
                  # ---------------------------
                  # Build diagnostic table
                  # ---------------------------
                  d <- ssl.dat %>%
                    transmute(
                      year = as.integer(year),
                      sst = sst_wgoa_coastwatch_junjulaug,
                      capelin = capelin_avg,
                      herring = herr_avg
                    ) %>%
                    arrange(year) %>%
                    mutate(
                      forage_unweighted = capelin + herring,
                      forage_weighted   = a*capelin + b*herring,
                      
                      temp_unweighted   = sst,
                      temp_weighted     = kT*sst,
                      
                      term_intercept    = k0,
                      term_temp         = kT*sst,
                      term_forage       = -kF*forage_weighted,
                      eta               = term_intercept + term_temp + term_forage,
                      p_switch          = plogis(eta)
                    )
                  
                  # ---------------------------
                  # Plot A: forage raw components + unweighted/weighted totals
                  # ---------------------------
                  forage_long <- d %>%
                    select(year, capelin, herring, forage_unweighted, forage_weighted) %>%
                    pivot_longer(-year, names_to = "series", values_to = "value")
                  
                  p_forage <- ggplot(forage_long, aes(year, value, color = series)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.4) +
                    theme_bw() +
                    labs(
                      title = "Forage Time Series: Raw Inputs and Combined Terms",
                      subtitle = sprintf("forage_weighted = a*capelin + b*herring   (a=%.3f, b=%.3f)", a, b),
                      x = "Year", y = "Forage value", color = NULL
                    )
                  
                  ggsave(file.path(out_dir, "diag_forage_raw_vs_weighted.png"), p_forage, width = 10, height = 5.8, dpi = 170)
                  
                  # ---------------------------
                  # Plot B: temperature raw vs weighted
                  # ---------------------------
                  temp_long <- d %>%
                    select(year, temp_unweighted, temp_weighted) %>%
                    pivot_longer(-year, names_to = "series", values_to = "value")
                  
                  p_temp <- ggplot(temp_long, aes(year, value, color = series)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.4) +
                    theme_bw() +
                    labs(
                      title = "Temperature Time Series: Raw vs Weighted",
                      subtitle = sprintf("temp_weighted = kT * SST   (kT=%.3f)", kT),
                      x = "Year", y = "Temperature term value", color = NULL
                    )
                  
                  ggsave(file.path(out_dir, "diag_temperature_raw_vs_weighted.png"), p_temp, width = 10, height = 5.8, dpi = 170)
                  
                  # ---------------------------
                  # Plot C: logistic predictor components + p_switch
                  # ---------------------------
                  comp_long <- d %>%
                    select(year, term_intercept, term_temp, term_forage, eta) %>%
                    pivot_longer(-year, names_to = "component", values_to = "value")
                  
                  p_eta <- ggplot(comp_long, aes(year, value, color = component)) +
                    geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.3) +
                    theme_bw() +
                    labs(
                      title = "Switch Linear Predictor Components",
                      subtitle = sprintf("eta = k0 + kT*SST - kF*forage_weighted   (k0=%.3f, kF=%.3f)", k0, kF),
                      x = "Year", y = "Linear predictor scale", color = NULL
                    )
                  
                  p_ps <- ggplot(d, aes(year, p_switch)) +
                    geom_line(linewidth = 1.1, color = "blue3") +
                    geom_point(size = 1.5, color = "blue4") +
                    theme_bw() +
                    labs(
                      title = "Switch Probability Over Time",
                      subtitle = "p_switch = logistic(eta)",
                      x = "Year", y = "p_switch"
                    ) +
                    ylim(0, 1)
                  
                  ggsave(file.path(out_dir, "diag_switch_eta_components.png"), p_eta, width = 10, height = 5.8, dpi = 170)
                  ggsave(file.path(out_dir, "diag_switch_probability.png"), p_ps, width = 10, height = 5.2, dpi = 170)
                  
                  # ---------------------------
                  # Optional standardized comparison for shape only
                  # ---------------------------
                  z <- function(x) as.numeric(scale(x))
                  d_z <- d %>%
                    transmute(
                      year,
                      sst_z = z(sst),
                      temp_weighted_z = z(temp_weighted),
                      forage_unweighted_z = z(forage_unweighted),
                      forage_weighted_z = z(forage_weighted)
                    ) %>%
                    pivot_longer(-year, names_to = "series", values_to = "zvalue")
                  
                  p_z <- ggplot(d_z, aes(year, zvalue, color = series)) +
                    geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
                    geom_line(linewidth = 1) +
                    geom_point(size = 1.3) +
                    theme_bw() +
                    labs(
                      title = "Standardized Inputs: Unweighted vs Weighted Shapes",
                      x = "Year", y = "z-score", color = NULL
                    )
                  
                  ggsave(file.path(out_dir, "diag_inputs_weighted_vs_unweighted_z.png"), p_z, width = 10, height = 5.8, dpi = 170)
                  
                  # Save underlying table
                  write_csv(d, file.path(out_dir, "diag_forage_temp_weighting_timeseries.csv"))
                  
                  message("Done. Wrote diagnostics to: ", out_dir)
                  
                  
                  print(p_z)
                  
                  # ---- Show all generated plots in RStudio and save a contact sheet ----
                  suppressPackageStartupMessages({
                    library(png)
                    library(grid)
                    library(gridExtra)
                  })
                  
                  plot_dir <- "copilot/outputs_2"
                  png_files <- list.files(plot_dir, pattern = "\\.png$", full.names = TRUE)
                  
                  if (length(png_files) == 0) {
                    message("No PNG files found in: ", plot_dir)
                  } else {
                    message("Found ", length(png_files), " plot files:")
                    print(basename(png_files))
                    
                    # 1) Print each figure one-by-one in the Plots pane
                    for (f in png_files) {
                      img <- png::readPNG(f)
                      grid::grid.newpage()
                      grid::grid.raster(img)
                      grid::grid.text(
                        label = basename(f), x = unit(0.02, "npc"), y = unit(0.98, "npc"),
                        just = c("left", "top"), gp = gpar(col = "white", fontsize = 11, fontface = "bold")
                      )
                      Sys.sleep(0.4)  # small pause so you can see each one
                    }
                    
                    # 2) Build a contact sheet image (all plots in one file)
                    grobs <- lapply(png_files, function(f) {
                      img <- png::readPNG(f)
                      grobTree(
                        rasterGrob(img, interpolate = TRUE),
                        textGrob(basename(f), x = 0.02, y = 0.98, just = c("left","top"),
                                 gp = gpar(col = "white", fontsize = 10, fontface = "bold"))
                      )
                    })
                    
                    n <- length(grobs)
                    ncol <- if (n <= 2) 1 else if (n <= 6) 2 else 3
                    nrow <- ceiling(n / ncol)
                    
                    contact <- gridExtra::arrangeGrob(grobs = grobs, ncol = ncol)
                    out_file <- file.path(plot_dir, "ALL_PLOTS_CONTACT_SHEET.png")
                    ggsave(out_file, contact, width = 4.8 * ncol, height = 3.4 * nrow, dpi = 160)
                    
                    message("Saved contact sheet: ", out_file)
                  }
                  