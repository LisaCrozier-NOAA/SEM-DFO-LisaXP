
#Script: 03_altprey_dat_10212026.r
#Purpose -- to create master data file for modeling alternative prey and SST effects on predators
#Output -- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")

library(tidyverse)
library(lubridate)
library(MARSS)

# --- Configuration & Paths ---
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM"
#doug_dir     <- "2026_10_21_SEM"

master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")
fig_out_dir  <- file.path(master_out, "figures")

crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")

master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

dfa_file  <- file.path(doug_dir, "shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds")
singletons_file  <- file.path(doug_dir, "allData.rds")


eulachon_file <- "copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv"
if (!file.exists(eulachon_file)) {
  eulachon_file <- file.path(dirname(rootdir), "copilot/outputs_8/eulachon_scaled_annual_weekly_during_chinook.csv")
}

ak_yr_file <- "data_Lisa/ak_yr.csv"
if (!file.exists(ak_yr_file)) {
  ak_yr_file <- file.path(dirname(rootdir), "data_Lisa/ak_yr.csv")
}


master_df <- readRDS(dfa_file)%>% mutate(year=year(date))
singletons_df <- readRDS(singletons_file)%>% mutate(year=year(date)) %>%
  filter(year>=1998, year<=2021)

prey_from_singletons_df<-c("biomassEuphShelfSum","sitkaHerring_EGoA","capelin_WGoA","pollockGOA_age1and2",
                           #   "biomassEuphShelfSum_adultyr", #excluded because it has NAs and I am not planning to use it right now
                           "sitkaHerring_EGoA_adultyr","capelin_WGoA_adultyr","pollockGOA_age1and2_adultyr",
                           "HCI","HCI_adultyr","SST_NCC_aprjun_DFObccm","SST_NCC_aprjun_DFObccm_adultyr","SST_EGOA_junaug","SST_EGOA_junaug_adultyr")




# get_prey_pool <- function(guild_num, stage_suffix, pred_col = "") {
#   sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
#   
#   m_squid     <- paste0("X04_marketsquid_GAM", sfx)
#   anchovy     <- paste0("X05_anchovy_GAM", sfx)
#   sardine     <- paste0("X05_DFA_abundSardine", sfx)
#   herring_ncc <- paste0("X05_herring_NCC", sfx)
#   eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
#   hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
#   
#   euphausiids <- paste0("biomassEuphShelfSum", sfx) #in DFA: X12_DFA_EGOA_DFA_lowerTrophic
#   
#   pollock     <- paste0("X13_pollockGOA_age1and2", sfx)
#   herring_ak  <- paste0("sitkaHerring_EGoA", sfx)
#   capelin     <- paste0("capelin_WGoA", sfx)
#   krill_ak    <- paste0("X12_egoa_krill", sfx)
#   pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
  
  
sort(unique(singletons_df$shortName))
unique(singletons_df$shortName)
  head(master_crosswalk)
  head(singletons_df)
  head(master_df)
  range(master_df$year)
  master_df %>% filter(year==1998)
  names(dat1)
  
  dat1<-master_df  %>%
    rename(DFAname=shortName) %>%
    right_join(master_crosswalk %>% select(DFAname,LisaName),by="DFAname") %>%
    select(year,LisaName,finalVal) %>%
    pivot_wider(id_cols = year,names_from = LisaName,values_from = finalVal) %>% 
    arrange(year)
  head(dat1)
  
  range(dat1$year)
  
  dat2<-singletons_df  %>%
    filter(shortName %in% prey_from_singletons_df) %>%
    select(year,shortName,finalVal) %>%
    pivot_wider(id_cols = year,names_from = shortName,values_from = finalVal) %>% 
    arrange(year)
  head(dat2)
  print(dat2,n=Inf)
  
  dat_prey<-full_join(dat1,dat2,by="year")
  head(dat_prey)
  print(names(dat_prey))
  

#add extra files-----
  eulachon_scaled <- read.csv(eulachon_file, stringsAsFactors = FALSE) %>%
    select(Year, eulachon_during_chinook, eulachon_annual_scaled) %>%
    rename(
      X05_eulachon_during_chinook_smoltyr = eulachon_during_chinook,
      X05_eulachon_annual_scaled_smoltyr   = eulachon_annual_scaled
    ) %>%
    mutate(
      X05_eulachon_during_chinook_adultyr = dplyr::lead(X05_eulachon_during_chinook_smoltyr, 2),
      X05_eulachon_annual_scaled_adultyr   = dplyr::lead(X05_eulachon_annual_scaled_smoltyr, 2)
    )
  
  ak_yr_selected <- read.csv(ak_yr_file, row.names = NULL, stringsAsFactors = FALSE) %>%
    rename(Year = year) %>%
    filter(Year >= 1996 & Year <= 2025) %>%
    select(
      Year,
      X11_ssl_seak_pup_pred        = ssl_seak_pup_pred,
      X12_egoa_krill_smoltyr       = secm_euph_dens,
      X13_mid_il_capelin_smoltyr   = mid_il_capelin,
      X21_sst_egoa_junjulaug_smoltyr = sst_egoa_coastwatch_junjulaug
    ) %>%
    mutate(
      X12_egoa_krill_adultyr         = dplyr::lead(X12_egoa_krill_smoltyr, 2),
      X13_mid_il_capelin_adultyr     = dplyr::lead(X13_mid_il_capelin_smoltyr, 2),
      X21_sst_egoa_junjulaug_adultyr = dplyr::lead(X21_sst_egoa_junjulaug_smoltyr, 2)
    )
  
  #X12_egoa_krill_adultyr X13_mid_il_capelin_adultyr X21_sst_egoa_junjulaug_adultyr all missing 2021 (i.e., 2023 not in Bridget's database)
  yr2020<-ak_yr_selected %>% filter(Year==2020) %>% select(X12_egoa_krill_adultyr,X13_mid_il_capelin_adultyr,X21_sst_egoa_junjulaug_adultyr);yr2020
  ak_yr_selected[26,c("X12_egoa_krill_adultyr", "X13_mid_il_capelin_adultyr", "X21_sst_egoa_junjulaug_adultyr")]<- yr2020
  
  #merge files and save alt data and predator file----------
  head(dat_prey)
  print(dat2,n=Inf)
  
  dat_all<-full_join(dat_prey,eulachon_scaled,by=c("year" = "Year"))
  dat_all<-full_join(dat_all,ak_yr_selected,by=c("year" = "Year")) %>%
  filter(year>=1998,year<=2021) %>% arrange(year)
  names(dat_all)
  print(dat_all %>% select(year,contains("X11")),n=Inf)
  
  
  # --- Standardize all non-suffixed base prey/predator columns to _smoltyr ---
  names(dat_all) <- ifelse(
    names(dat_all) != "year" & !grepl("_smoltyr$|_adultyr$", names(dat_all)),
    paste0(names(dat_all), "_smoltyr"),
    names(dat_all)
  )
  
  # Save updated master file
  write.csv(dat_all, file = file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv"))
  
  valid_cols <- names(dat_all)[colSums(!is.na(dat_all)) == 24]
  print(valid_cols)
  
  