#Explore simple correlations with HCI and SST

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION SWITCH & PATHS
# ------------------------------------------------------------------------------
#years_included <- "19yr" # Switch: "19yr" or "24yr"
years_included <- "24yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)

names(raw_df)
dat<-raw_df %>% select(Year,X01_habCompInd_smoltyr,X21_sst_egoa_junjulaug_adultyr) %>%
  mutate(hci=-X01_habCompInd_smoltyr) %>% 
  select(Year,hci,X21_sst_egoa_junjulaug_adultyr) 

n=ncol(dat)-1
matplot(dat$Year,dat[,-1],type='b',lty=1:n,col=1:n)
legend("topright",legend=names(dat)[-1],lty=1:n,col=1:n)
