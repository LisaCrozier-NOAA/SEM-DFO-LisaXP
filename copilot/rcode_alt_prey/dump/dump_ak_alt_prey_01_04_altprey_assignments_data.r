

# This script loads the current list of indicators (combined_guildfiles), 
# reduces it to the predators of interest for alternate prey,
# and assigns manually 2 options for alternate prey
# it outputs those assignments to "copilot/outputs_8/nccpred_altprey_assigned.csv" using a lower case version of short_name_lower

library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(janitor)
###########################################################################
#1. decide on altprey1 & altprey2 for all ncc predators--------
###########################################################################

#get the current list of indicators for SEM------

outputDir <- file.path( "copilot/outputs_8")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

combined_guildfiles<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/guildsWithExclude.csv", row.names = NULL) %>% 
      clean_names() %>% 
      rename(sem_name=se_mnode) %>%
      rename(short_name_lower=short_name) %>%
      mutate(across(where(is.character), tolower)) %>% 
      filter(!is.na(latest_guild) )     
    
head(combined_guildfiles)
combined_guildfiles[which(is.na(combined_guildfiles$short_name_lower)),]
    
    

#look at all of the predators by shortname, decide on preferred prey---------
nccpred <- combined_guildfiles %>% 
  filter(grepl("08|09|10|11", latest_guild)  ) %>% 
  select(indicator, short_name_lower, sem_name,latest_guild) %>%
  mutate(
    altprey1 = case_when(
      # Focal Chinook
      grepl("killer.whale", indicator) ~ "chin",
      
      # Early-run River / Estuary Target
      grepl("allsealionsbonn", indicator) ~ "eulachon",
      
      # Seabirds (Guild 08)
      grepl("08\\.predbirdncc|peli|tern|murre|corm|alcid|shrw|shearwater|grebe", indicator) ~ "market_squid",
      
      # Invertebrate / Planktivorous Feeders
      grepl("hake|rockfish|chilipepper|mackerel|chinook abundance", indicator) ~ "krill",
      
      # Standardized Hake Primary (Marine Mammals & Demersal Predators)
      .default = "hake"
    ),
    
    altprey2 = case_when(
      # Focal Salmonid Predators
      grepl("killer.whale|allsealionsbonn|allsealionsemb", indicator) ~ "salmonids",
      
      # ANCHOVY PREFERRERS
      # - California Sea Lions
      # - Seabirds (Guild 08)
      # - Pelagic Mackerels, Rockfish, & Adult Chinook Stocks
      grepl("californian_s_l", indicator) ~ "anchovy",
      grepl("08\\.predbirdncc|peli|tern|murre|corm|alcid|shrw|shearwater|grebe", indicator) ~ "anchovy",
      grepl("rockfish|chilipepper|mackerel|chinook abundance", indicator) ~ "anchovy",
      
      # HAKE PREFERRERS (Adult Hake diet)
      grepl("hake", indicator) ~ "anchovy",
      
      # HERRING PREFERRERS (Default for remaining mammals & demersals)
      # - Steller Sea Lions, Harbor Seals, Fur Seals, Harbor Porpoise
      # - Lingcod, Halibut, Skates, Dogfish
      .default = "herring"
    )
  ) %>% 
  arrange(altprey1, altprey2)

# Display complete output
print(as_tibble(nccpred), n = Inf)

#save prey assignments----------
write.csv(nccpred, file.path(outputDir, "nccpred_altprey_assigned.csv"), row.names = FALSE)



#2. get prey datasets------------

# This script processes eulachon annual and weekly datafiles from outputs_csl_cr
# then it joins that with the extended prey datafiles from doug_01 and doug_02 (clusDataDFA_wide_1996_2025.csv) and 

library(dplyr)
library(stringr)
library(imputeTS)

# ---------------------------------------------------------------------------
# Load and Process Eulachon Datasets-------
# ---------------------------------------------------------------------------

        # Annual Run Index (Reconstructed Lbs)
        eulachon_annual <- read.csv(
          "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/eulachon_master_index_lbs_1993_2024.csv",
          row.names = NULL
        ) %>%
          select(year, eulachon_lbs_reconstructed) %>%
          arrange(year)
        
        eulachon_annual$eulachon_lbs_interp <- na_interpolation(eulachon_annual$eulachon_lbs_reconstructed)
        
        # Weekly Overlap Index (Active Chinook Window)
        eulachon_weekly <- read.csv(
          "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/copilot/outputs_csl_cr/master_estuary_2wklagweekly_all_species_1998_2024.csv",
          row.names = NULL
        ) %>%
          group_by(year) %>%
          summarise(eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE), .groups = "drop") %>%
          select(year, eulachon_during_chinook)
        
        # Combine and apply Log + Z-score scaling
        eulachon_clean <- full_join(eulachon_weekly, eulachon_annual %>% select(year, eulachon_lbs_interp), by = "year") %>%
          rename(Year = year) %>%
          mutate(
            # Primary index: Weekly overlap
            eulachon_weekly_scaled = as.vector(scale(log1p(eulachon_during_chinook))),
            # Secondary index: Annual reconstructed run size
            eulachon_annual_scaled = as.vector(scale(log1p(eulachon_lbs_interp)))
          )
        
        write.csv(eulachon_clean, file.path(outputDir, "eulachon_scaled_annual_weekly_during_chinook.csv"), row.names = FALSE)
        
# ---------------------------------------------------------------------------
# Load sem_data and select prey series -----
# ---------------------------------------------------------------------------

sem_data<-read.csv("data_Lisa/sem_master_data.csv",row.names = NULL)
names(sem_data)

clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))
names(clusDataDFA_wide)


# Construct wide prey dataset with exact raw source column names and 2-year leads
altprey_guildnames <- clusDataDFA_wide %>%
  select(
    Year,
    X09.PredFishNCC_DFA1,        # chin
    X09.PredFishNCC_b_DFA1,      # hake
    X12.ZooPreyAK_DFA1,          # krill
    X04.CompFishNCC_smoothed,    # market_squid
    X05.ForageFishNCC_0_smoothed,# anchovy
    X13.FishPreyAK_b_DFA1        # herring
  ) %>%
  left_join(eulachon_clean %>%
            select(Year, eulachon_during_chinook, eulachon_annual), 
              by = "Year") %>%
  mutate(
    # Generate matching _2yrlead columns using exact source names
    X09.PredFishNCC_DFA1_2yrlead         = lead(X09.PredFishNCC_DFA1, 2),
    X09.PredFishNCC_b_DFA1_2yrlead       = lead(X09.PredFishNCC_b_DFA1, 2),
    X12.ZooPreyAK_DFA1_2yrlead           = lead(X12.ZooPreyAK_DFA1, 2),
    X04.CompFishNCC_smoothed_2yrlead     = lead(X04.CompFishNCC_smoothed, 2),
    X05.ForageFishNCC_0_smoothed_2yrlead = lead(X05.ForageFishNCC_0_smoothed, 2),
    X13.FishPreyAK_b_DFA1_2yrlead        = lead(X13.FishPreyAK_b_DFA1, 2),
    eulachon_during_chinook_2yrlead      = lead(eulachon_during_chinook, 2),
    eulachon_annual_2yrlead              = lead(eulachon_annual, 2)
  )

#rename manually to better sem_data names
ncc_altprey_data <- altprey_guildnames  %>%
  rename(
    X09_DFA_ChinAbundSnakeFall = X09.PredFishNCC_DFA1,        # chin
    X09_DFA_HakeAge5Plus = X09.PredFishNCC_b_DFA1,      # hake
    X12_DFA_biomassEuphShelfSum = X12.ZooPreyAK_DFA1,          # krill
    X04_marketsquid_GAM = X04.CompFishNCC_smoothed,    # market_squid
    X05_anchovy_GAM = X05.ForageFishNCC_0_smoothed,# anchovy
    X05_DFA_abundSardine = X13.FishPreyAK_b_DFA1,        # herring
    
    X09_DFA_ChinAbundSnakeFall_2yrlead = X09.PredFishNCC_DFA1_2yrlead         ,
    X09_DFA_HakeAge5Plus_2yrlead = X09.PredFishNCC_b_DFA1_2yrlead      ,
    X12_DFA_biomassEuphShelfSum_2yrlead =  X12.ZooPreyAK_DFA1_2yrlead          ,
    X04_marketsquid_GAM_2yrlead = X04.CompFishNCC_smoothed_2yrlead    ,
    X05_anchovy_GAM_2yrlead = X05.ForageFishNCC_0_smoothed_2yrlead ,
    X05_DFA_abundSardine_2yrlead = X13.FishPreyAK_b_DFA1_2yrlead        
   
  )

# Save master prey dataset -----------
write.csv(ncc_altprey_data, file.path(outputDir, "ncc_altprey_data_1996_2025.csv"), row.names = FALSE)


