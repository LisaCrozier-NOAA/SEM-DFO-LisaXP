
library(lubridate)
library(tidyverse)

"C:\Users\Lisa.Crozier C drive from work\Marine survival\Doug results\analyzeAKindices\guilds.csv"
path_pc<-"C:/Users/Lisa.Crozier C drive from work/Marine survival/Doug results/analyzeAKindices"
path_doug<-c("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices",
            "C:/Users/Lisa.Crozier C drive from work/Marine survival/Doug results/analyzeAKindices")[2]

doug_data<-readRDS(file.path("Rcode_for_paper","metadata","Doug_10132026_adultyr_allData.rds"));names(doug_data)

doug_guilds<-read.csv(file.path(path_doug,"guilds.csv"),na.strings = c("","NA"))%>% select(-mismatch_09oct26_vs_shortName);names(doug_guilds)
guilds<-doug_guilds %>% select(node,SEMnode,region,source_short,source_raw,shortName,Lisa_21oct26,region);guilds %>% head()

guilds_current<-guilds %>% filter(!is.na(Lisa_21oct26)) ;guilds_current %>% head()
guilds_current %>% group_by(SEMnode) %>%
  summarize(source_short=n())

print(guilds_current %>% group_by(source_short) %>%
  summarize(SEMnode=n()),n=Inf)

guilds_current %>% head()
guilds_current %>% filter(is.na(source_short))
guilds_current %>% filter(source_short=="erddap")
erddap
tapply(guilds_current,list(region))

guilds_current %>% group_by(region) %>%
  summarize(SEMnode=n())


guilds_current %>%   summarize(region=n())

guilds<-doug_guilds %>% select(shortName,Lisa_21oct26,region);guilds

guilds %>% filter(region=="CC")


doug_guilds_10092026 <-doug_guilds
doug_dat<-doug_data %>% select(shortName,yearType,date,finalVal,smoothed) %>%
  left_join(guilds,by="shortName") %>%
  mutate(date=year(date)) %>%
  rename(year=date) %>%
  filter(year>=1996, year<=2023)
tail(doug_dat)

doug_guilds  %>% filter(!is.na(Lisa_21oct26)|exclusionReason) %>% head()
doug_guilds %>% select(shortName,SEMnode,exclusionReason,replacementShortName,pastGuild,Lisa_21oct26)

doug_dat %>% filter(year==1998,!is.na(latestGuild)) %>% head()
doug_dat %>% filter(year==1998,yearType=="adult",!is.na(latestGuild)) %>% head()

doug_dat %>% filter(latestGuild=="12.ZooPreyAK",year==1998)
doug_dat %>% filter(grepl("15.PredFishAK",latestGuild),year==1998)
summary(doug_data$imputed)           

doug_dat %>% filter(grepl("lead",shortName))
doug_data %>% filter(grepl("adult",yearType)) %>% head()
doug_data %>% filter(grepl("adult",yearType)) %>% head()

unique(doug_dat$shortName)
unique(doug_dat$latestGuild)

doug_dat_wid<- pivot_wider(doug_dat %>%
                             filter(grepl("15.",latestGuild)) %>% 
                             select(shortName,year,finalVal),
                           id_cols=everything(),names_from = shortName,values_from = finalVal)
head(doug_dat_wid)
names((doug_dat_wid))


#What is the current species list?
#new guilds shortname scheme
doug_guilds<-read.csv(file.path(path_doug,"guilds.csv"))%>% 
  filter(!is.na(latestGuild))%>%
  mutate(node=paste0("X",substr(latestGuild,1,2)), 
         guild_main=str_split_i(latestGuild,"_",1),
         guild2=str_split_i(latestGuild,"_",2),
         guild3=str_split_i(latestGuild,"_",3)) %>% 
  mutate(
         guild_sub=case_when(
           guild2 %in% c("WPbirds", "JSOES",   "NHL","b") ~ paste(guild_main,guild2,sep="_"),
          .default = guild_main) ) %>%
  select(-guild2,-guild3)

         
head(doug_guilds)

unique(doug_guilds$exclusionReason)
doug_guilds %>% filter(node=="X15")  %>%
  select(shortName,source_short, source_raw,node,   guild_main,       guild_sub)

doug_guilds %>% filter(node=="X13")  %>%
  select(shortName,source_short, source_raw,node,   guild_main,       guild_sub)

#Lisa guild assignments---------

doug_guilds_Lisa %>% filter(node=="X13", !is.na(Lisa_15Oct2026))  %>%
  select(shortName,source_short, source_raw,node,   guild_main,       guild_sub,Lisa_15Oct2026)


write.csv(doug_guilds_Lisa,"Lisa_guilds.csv",row.names = FALSE)



# Explain exclusions and replaced indicators-------
      library(dplyr)
      library(tibble)
      
      
      # 1. Define exclusions for existing indicators
      exclusion_map <- tribble(
        ~shortName,                   ~Lisa_15Oct2026, ~exclusionReason,                       ~replacementShortName,
        "sharkCatchGoA_predAK_2026",  NA_character_,   "selected_spp_specific_indicators",     "sleeperSharkGoA_predAK_2026",
        
        "rhinoAuk_EGoA",              NA_character_,   "included_in_EGOA_DFA_seabirds_2026",  "EGOA_DFA_seabirds_2026",
        "ammod_EAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        "gadid_EAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        "ammod_WAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        "gadid_WAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        "hexagram_EAI",               NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        "hexagram_WAI",               NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
        
        "pollockBiomassGoAage3plus_predAK_2026",NA_character_, "need_larger_preds",           "pollockGOA_SSB",
        "pollockBiomassAIage1plus_predAK_2026", NA_character_, "Aleutians_instead_of_GOA",    "pollockGOA_age1and2",
        "ArrowtoothFlounderBiomass_predAK_2026",NA_character_, "need_larger_preds",           "ArrowtoothFlounderFSB",
        "PacificCodBiomass_predAK_2026",        NA_character_, "need_larger_preds",           "PacificCodGOA_FSB",
        "sablefishBiomass_predAK_2026",         NA_character_, "need_larger_preds",           "sablefishAK_biomass3plus"
        
      )
      
      # 2. Merge changes into guild database
      doug_guilds_Lisa <- doug_guilds %>%
        left_join(exclusion_map, by = "shortName", suffix = c(".orig", ".map")) %>%
        mutate(
          # Takes the tribble value if present; otherwise defaults to original 
          Lisa_15Oct2026 = coalesce(Lisa_15Oct2026, guild_sub),
          exclusionReason = coalesce(exclusionReason.map, exclusionReason.orig),
          replacementShortName = coalesce(replacementShortName.map,replacementShortName.orig )
        ) %>%
        select(-ends_with(".orig"), -ends_with(".map")) %>% # Clean up suffix columns
        arrange(node, shortName, Lisa_15Oct2026)


#add new indicators------------
      #new indicators:
      # pollockGOA_SSB
      # pollockGOA_SSB_lag1
      # pollockGOA_age1and2
      # ArrowtoothFlounderFSB                        
      # PacificCodGOA_FSB
      # sablefishAK_biomass3plus                     
      # sablefishAK_biomass4plus                     

names(doug_guilds_Lisa)
# [1] "ID"                   "indicator"            "shortName"            "dataset"              "region"               "subregion"            "trophicGroup"         "impute"               "logTransform"         "incl2026"            
# [11] "SEMnode"              "shortName_09oct26"    "pastGuild"            "latestGuild"          "source_short"         "source_raw"           "node"                 "guild_main"           "guild_sub"            "Lisa_15Oct2026"      
# [21] "exclusionReason"      "replacementShortName"

doug_guilds_Lisa %>% filter(shortName=="ArrowtoothFlounderBiomass_predAK_2026")

doug_guilds_Lisa<-doug_guilds_Lisa_start
#pollockGOA_SSB
newrow<-doug_guilds_Lisa %>% filter(shortName=="pollockBiomassGoAage3plus_predAK_2026")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
                  c("pollockGOA_SSB","pollockGOA_SSB","X15","15.PredFishAK","15.PredFishAK")
          newrow

          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
#pollockGOA_SSB_lag1          
          doug_guilds_Lisa %>% filter(shortName=="pollockGOA_SSB")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="pollockGOA_SSB")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("pollockGOA_SSB_lag1","pollockGOA_SSB_lag1","X15","15.PredFishAK","15.PredFishAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
          
#pollockGOA_age1and2          
          doug_guilds_Lisa %>% filter(shortName=="pollockGOA_SSB")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="pollockGOA_SSB")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("pollockGOA_age1and2","pollockGOA_age1and2","X13","13.FishPreyAK","13.FishPreyAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
#ArrowtoothFlounderFSB          
          doug_guilds_Lisa %>% filter(shortName=="ArrowtoothFlounderBiomass_predAK_2026")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="ArrowtoothFlounderBiomass_predAK_2026")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("ArrowtoothFlounderFSB","ArrowtoothFlounderFSB","X15","15.PredFishAK","15.PredFishAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
          
#PacificCodGOA_FSB          
          doug_guilds_Lisa %>% filter(shortName=="ArrowtoothFlounderBiomass_predAK_2026")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="ArrowtoothFlounderBiomass_predAK_2026")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("PacificCodGOA_FSB","PacificCodGOA_FSB","X15","15.PredFishAK","15.PredFishAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
          
          
#sablefishAK_biomass3plus          
          doug_guilds_Lisa %>% filter(shortName=="sablefishBiomass_predAK_2026")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="sablefishBiomass_predAK_2026")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("sablefishAK_biomass3plus","sablefishAK_biomass3plus","X15","15.PredFishAK","15.PredFishAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
          
          
          
#sablefishAK_biomass4plus          
          doug_guilds_Lisa %>% filter(shortName=="sablefishBiomass_predAK_2026")
          
          newrow<-doug_guilds_Lisa %>% filter(shortName=="sablefishBiomass_predAK_2026")
          newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
          newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
            c("sablefishAK_biomass4plus","sablefishAK_biomass4plus","X15","15.PredFishAK","15.PredFishAK")
          newrow
          doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 

          
          tail(doug_guilds_Lisa)
          
doug_guilds_Lisa_start<-doug_guilds_Lisa


#review notes about questions on original exclude list ----------         
doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")

doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")
doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")
doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")
doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")

                    
          
      #how stronlgy correlated are birds and seanettles?------
path_doug<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/"
#"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/completeness.csv"
doug_complete<-read.csv(file.path(path_doug,"2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/completeness.csv")) %>%
    rename(year=date) %>% mutate(year=year(year));
head(doug_complete[,1:5])

head(master_crosswalk)
crosswalk<-master_crosswalk

doug_data<-readRDS(file.path("Rcode_for_paper","metadata","Doug_10132026_adultyr_allData.rds"));head(doug_data)
doug_guilds<-read.csv(file.path(path_doug,"guilds.csv"))%>% filter(!is.na(latestGuild));head(doug_guilds)


library(dplyr)
library(tibble)

# 1. Create a named vector: names are old column names (dfa_cols), values are new column names (LisaName)
name_map <- deframe(master_crosswalk %>% select(dfa_cols, LisaName))

# 2. Rename columns dynamically
doug_complete <- doug_complete %>%
  rename_with(~ ifelse(.x %in% names(name_map), name_map[.x], .x))

master_crosswalk %>% filter(grepl("X01",LisaName))

dat<-doug_complete[,c()]
n<-ncol(dat)-1
matplot(dat$year,dat[,c("X08")])

          