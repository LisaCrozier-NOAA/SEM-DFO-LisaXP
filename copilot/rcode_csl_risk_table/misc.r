library(readxl)


 csl_weekly_1998_2024<-read.csv(file.path(output_dir, "csl_reconstructed_weekly_1998_2024.csv"), row.names = NULL)
 eulachon_weekly_1998_2024<- read.csv(file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"), row.names = NULL)
 
 bonn_daily<- read_excel("data_Lisa/Adult_BONpassage_1976_2024.xlsx",sheet=1);head(bonn_daily)
 shad_daily<- read.csv("data_Lisa/SHAD_bonn_19382026.csv",row.names = NULL);head(shad_daily)
 
 unique(bonn_daily$parameter)
 
 "C:\Users\Lisa.Crozier\Documents\Marine survival\SEM-DFO-LisaXP\data_Lisa\Adult_BONpassage_1976_2024.xlsx"
 
 "C:\Users\Lisa.Crozier\Documents\Marine survival\SEM-DFO-LisaXP\data_Lisa\Adult_oncorhynchus_bonn_1976_2024.csv"
 
 
 
 library(janitor)
 
 #hake and forage fish -- first Doug, then Jake
 output_dir <- "copilot/outputs_csl_cr"
 
 #Doug data---------
 salmon_dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
   clean_names()  %>% 
   select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
 names(salmon_dat)
 
 
 #jake data--------
 jake_annual<-read.csv(file.path(output_dir,"jake_lre_dat_yearly.csv")) 
 jake_annual %>% select(
   grepl("year|hake|sardine|herring|anchovy|mackerel",names(jake_annual)))
 
 names(jake_annual)
 
 altprey<-c("hake_cpue")
 spp_northern_names <- c("Engraulis mordax", "Sardinops sagax")
 spp_all_names <- c("Clupea pallasii", "Trachurus symmetricus")
 
 
 #REsults of hake and forage fish corr------
 Metric_A                   Metric_B                        pearson_r  abs_r
 1 x09_dfa_hake_age5plus_doug sardine_northern_bmass_jake       -0.941  0.941 
 2 x09_dfa_hake_age5plus_doug herring_bmass_jake                 0.938  0.938 
 3 x09_dfa_hake_age5plus_doug abund_sardine_doug                -0.827  0.827 
 4 x09_dfa_hake_age5plus_doug abund_herring_doug                 0.711  0.711 
 5 x09_dfa_hake_age5plus_doug x05_anchovy_gam_doug               0.709  0.709 
 6 x09_dfa_hake_age5plus_doug pacificmackerel_gam_2025_doug      0.702  0.702 
 7 x09_dfa_hake_age5plus_doug jackmackerel_gam_2025_doug         0.699  0.699 
 8 x09_dfa_hake_age5plus_doug x05_dfa_abund_sardine_doug        -0.361  0.361 
 
 
 Metric_A                   Metric_B                        pearson_r   abs_r
 1 x05_dfa_abund_sardine_doug sardine_northern_bmass_jake       0.970   0.970  
 2 x05_dfa_abund_sardine_doug abund_sardine_doug                0.924   0.924  
 3 x05_dfa_abund_sardine_doug sardine_gam_2025_doug             0.890   0.890  
 4 x05_dfa_abund_sardine_doug abund_herring_doug               -0.866   0.866  
 5 x05_dfa_abund_sardine_doug anchovy_northern_bmass_jake       0.854   0.854  
 6 x05_dfa_abund_sardine_doug x05_herring_gam_doug             -0.822   0.822  
 7 x05_dfa_abund_sardine_doug sardine_ncc_doug                  0.805   0.805  
 8 x05_dfa_abund_sardine_doug pacificmackerel_gam_2025_doug    -0.322   0.322 
 
 
 
 print(strongly_correlated_with_anchors)
 [1] "sardine_northern_bmass_jake"   "herring_bmass_jake"           
 [3] "abund_sardine_doug"            "abund_herring_doug"           
 [5] "x05_anchovy_gam_doug"          "pacificmackerel_gam_2025_doug"
 [7] "jackmackerel_gam_2025_doug"    "sardine_gam_2025_doug"        
 [9] "anchovy_northern_bmass_jake"   "x05_herring_gam_doug"         
 [11] "sardine_ncc_doug"             