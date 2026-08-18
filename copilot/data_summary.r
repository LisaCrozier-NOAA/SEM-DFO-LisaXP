

#Master file for all candidate data analyzed


#sem_master_data.csv for doug-created DFAs -- MARSS or DFA smoothed, predators lagged 2 years

# data_Lisa/goa_prey_clim_raw_trends_avg.csv for AK data already combined see Lisa.coalate.DFAtrends.r for source file
# data_Lisa/ak_yr for all raw AK data
# no trends in #Ferris-DFAIndicators-goa\mybioclim_egoa_wgoa.csv for AK data already combined
      #EGOA_EcoState_Data_Jan2023.csv & WGOA_EcoState_Data_Jan2023.csv for new AK data

#AKshark.table19.3_2022assess.csv for raw shark data

#SEM-DFO-NMFS/Data/NCC/ssl_counts.csv for raw Doug data to get all years

#SSL DATA------
#path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-NMFS/SEM-DFO-NMFS/Data/NCC/ssl_counts.csv"
ssl.orig<-read.csv(path) %>% select(-est_se) %>%
  rename(ssl.count.eric=count,ssl.model.eric=est);head(ssl.orig)



# ---------------------------
# 0) Paths
# ---------------------------
path_sem   <- "data_Lisa/sem_master_data.csv"
path_shark <- "data_Lisa/AKshark.table19.3_2022assess.csv"

path_egoa  <- "Ferris-DFAIndicators-goa/data/EGOA_EcoState_Data_Jan2023.csv"
path_wgoa  <- "Ferris-DFAIndicators-goa/data/WGOA_EcoState_Data_Jan2023.csv"

out_dir    <- "copilot/outputs_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



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
