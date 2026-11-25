


#Script: 01_data_prep/04_harvest.r


suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lavaan)
  library(gt)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS & MASTER DATA LOAD
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
meta_dir   <- file.path(proj_dir, "metadata")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")

pinn_path <- file.path(proj_dir, "metadata", "Braun et al 2025 ACOE Pinn Pred Bonn Tailrace Table5.xlsx")
      pinn <- read_xlsx(pinn_path, sheet = 1,skip=3)       %>%
        rename(pinn_pred= `Percent of\r\nRun Consumed`,
               nChin_eaten = `Chinook Consumption Estimate`,
               nChin_bonn = `Bonneville Dam Spring Chinook\r\nPassage`) %>%
        mutate(smoltyr = as.numeric(Year)-2) %>%
        select(smoltyr,nChin_bonn,nChin_eaten,pinn_pred)

        pinn_ext <- tibble(
          smoltyr     = c(1998, 1999),
          nChin_bonn  = c(NA_real_, NA_real_),
          nChin_eaten = c(NA_real_, NA_real_),
          pinn_pred   = c(0.003, 0.003)
        ) %>%
          bind_rows(pinn) %>%
          arrange(smoltyr) %>%
          filter(smoltyr >= 1998, smoltyr <= 2021) %>%
          mutate(
            # Instantaneous pinniped mortality rate: F = -ln(1 - M)
            lnSurv_pinn = -log(pmax(1e-4, 1 - pinn_pred))
          ) 
        print(pinn_ext,n=Inf)

harvest_path <- file.path(proj_dir, "metadata", "JointStaffReport 2025 SPRINGCHINOOK.table5.csv")
      harvest <- read.csv(harvest_path, stringsAsFactors = FALSE, row.names = NULL) %>%
        rename(harvest_z1_to_z5= Harvest.impact.downstream.of.Bonneville.Dam..Zones.1.5..Grand.Total,
               nmouth = Upriver.Run,
               nbonn = Bonneville.Dam.Count) %>%
        filter(!id %in% c('80-84','85-89','90-94')) %>%
        mutate(smoltyr = as.numeric(id)-2) %>%
        filter(smoltyr>=1998, smoltyr<=2021) %>%
        select(smoltyr,id,nmouth,nbonn,harvest_z1_to_z5)


bird_path <- file.path(proj_dir, "metadata", "aggregate_annual_estuary_avian_predation_on_chinook_fromBON_1998_2018.csv")
      bird <- read.csv(bird_path, stringsAsFactors = FALSE, row.names = NULL) %>% 
        group_by(year) %>%
        summarize(birdmort=mean(mean)) %>%
        select(year,birdmort)
  
      
data_path <- file.path(proj_dir, "metadata", "sem_altprey_data_1998_2021_10212026.csv")
      model_dat <- read.csv(data_path, stringsAsFactors = FALSE, row.names = 1);head(model_dat)
      

head(harvest)
head (bird)
head(model_dat %>% select(year,X16_DFA_sar_UCsprCh_smoltyr))
  
     