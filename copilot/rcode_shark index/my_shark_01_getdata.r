#Aug 6, 2026, preparing for Eve handoff

out_dir    <- "copilot/outputs_5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#Load data
        suppressPackageStartupMessages({
          library(tidyverse)
          library(janitor)
          library(lavaan)
          library(glue)
        })
        
        
        #sem_master_data.csv for doug-created DFAs -- MARSS or DFA smoothed, predators lagged 2 years
        
        # data_Lisa/goa_prey_clim_raw_trends_avg.csv for AK data already combined see Lisa.coalate.DFAtrends.r for source file
        # data_Lisa/ak_yr for all raw AK data
        # no trends in #Ferris-DFAIndicators-goa\mybioclim_egoa_wgoa.csv for AK data already combined
        #EGOA_EcoState_Data_Jan2023.csv & WGOA_EcoState_Data_Jan2023.csv for new AK data
        
        #AKshark.table19.3_2022assess.csv for raw shark data
        
        #SEM-DFO-NMFS/Data/NCC/ssl_counts.csv for raw Doug data to get all years
        
        #SSL DATA------
        #path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-NMFS/SEM-DFO-NMFS/Data/NCC/ssl_counts.csv"
        # ssl.orig<-read.csv(path) %>% select(-est_se) %>%
        #   rename(ssl.count.eric=count,ssl.model.eric=est);head(ssl.orig)
        
        
        sem.dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
          clean_names()
          names(sem.dat)
        
        ak_yr<-read.csv(file.path("data_Lisa/ak_yr.csv"))
        
        goa_prey<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
          clean_names() %>%
          select(1:12)
          names(goa_prey)
          
        sst_dat<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
            clean_names() %>%
            select(19:32)
          names(sst_dat)
        
        shark_raw<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
          clean_names() %>%
          select(contains("goa"))
          names(shark_raw)
        
        ssl_raw<-read.csv(file.path("copilot/outputs_2/ssl.dat.csv"))%>% 
            clean_names() %>%
          select(year,contains("ssl"))
          names(ssl_raw)

        salmon_dat<-sem.dat %>% select(contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
          names(salmon_dat)

            
