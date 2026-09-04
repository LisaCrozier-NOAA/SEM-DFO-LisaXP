

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
#decide on altprey1 & altprey2 for all ncc predators--------
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
