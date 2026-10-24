do you see any obvious errors or contradictions in this code? I'm sure it could have been written much more elegantly and less prone to errors: 


#Update to guilds.csv 16 oct, 2026

library(lubridate)
library(tidyverse)

path_doug<-c("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices","C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices")[1]
doug_guilds0<-read.csv(file.path(path_doug,"guilds.csv"));head(doug_guilds)


  #1. Parse out original main guild, which is used for Lisaguilds prior to DFA separating stragglers--------
        doug_guilds<- doug_guilds0 %>%
          mutate(node=paste0("X",substr(latestGuild,1,2)), 
                 guild_main=str_split_i(latestGuild,"_",1),
                 guild2=str_split_i(latestGuild,"_",2),
                 guild3=str_split_i(latestGuild,"_",3)) %>% 
          mutate(
            guild_sub=case_when(
              guild2 %in% c("WPbirds", "JSOES",   "NHL","b") ~ paste(guild_main,guild2,sep="_"),
              .default = guild_main) ) %>%
          select(-guild2,-guild3)


        
  #2. Explain exclusions and replaced indicators-------
        library(dplyr)
        library(tibble)
        
        
        # 1. Define exclusions for existing indicators
        exclusion_map <- tribble(
          ~shortName,                   ~Lisa_15Oct2026, ~exclusionReason,                       ~replacementShortName,
          "sharkCatchGoA_predAK_2026",  NA_character_,   "selected_spp_specific_indicators",     "sleeperSharkGoA_predAK_2026",
          
          "rhinoAuk_EGoA",              NA_character_,   "included_in_EGOA_DFA_seabirds_2026",  "EGOA_DFA_seabirds_2026",
          "ammod_EAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          "gadid_EAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          "ammod_WAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          "gadid_WAI",                  NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          "hexagram_EAI",               NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          "hexagram_WAI",               NA_character_,   "excluding_AI_stragglers",             "WGOA_DFA_lowerTrophic_2026",
          
          "pollockBiomassGoAage3plus_predAK_2026",NA_character_, "need_larger_preds",           "pollockGOA_SSB",
          "pollockBiomassAIage1plus_predAK_2026", NA_character_, "Aleutians_instead_of_GOA",    "pollockGOA_age1and2",
          "ArrowtoothFlounderBiomass_predAK_2026",NA_character_, "need_larger_preds",           "ArrowtoothFlounderFSB",
          "PacificCodBiomass_predAK_2026",        NA_character_, "need_larger_preds",           "PacificCodGOA_FSB",
          "sablefishBiomass_predAK_2026",         NA_character_, "need_larger_preds",           "sablefishAK_biomass3plus"
          
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
        
        
  #3.Add new indicators------------
        #new indicators:
        # pollockGOA_SSB
        # pollockGOA_SSB_lag1
        # pollockGOA_age1and2
        # ArrowtoothFlounderFSB                        
        # PacificCodGOA_FSB
        # sablefishAK_biomass3plus                     
        # sablefishAK_biomass4plus                     
        
        names(doug_guilds_Lisa)
        # [1] "ID"                   "indicator"            "shortName"            "dataset"              "region"               "subregion"            "trophicGroup"         "impute"               "logTransform"         "incl2026"            
        # [11] "SEMnode"              "shortName_09oct26"    "pastGuild"            "latestGuild"          "source_short"         "source_raw"           "node"                 "guild_main"           "guild_sub"            "Lisa_15Oct2026"      
        # [21] "exclusionReason"      "replacementShortName"
        

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
        
#-------------------------------------------------------        
  #4. Review notes about previously excluded indicators that should be put back in ----------    
        #              shortName SEMnode                                    exclusionReason
        # 1         chubMack_NCC PreyNCC ?_counts_too_low_selected_pacificmackerel_GAM_2025
        # 2        predSealAbund  PredAK                                     ?_ends_in_2019
        # 3       Lmu_IntSprMayH  Growth                                    ?_justFWeffects
        # 4        northAnch_NCC PreyNCC                       ?_lastused:labelDoug_30may25
        # 5 sablefish_WCGFS_2025 PredNCC                       ?_lastused:labelLisa_21aug25
        # 6     jackMackerel_BRT PreyNCC               ?_selected_gam_over_brt_for_all_sdms
        
        #Notes: 
        #Need to change source for some of the JSOES files
        #DFO files -- checking with Cam
        #Some probably should be added back in: "last used"
        
        #          shortName data_updated            ExplainExclusion2
        # 1             chubMack_NCC         <NA> ?_lastused:labelDoug_30may25
        # 2         jackMackerel_BRT         <NA> ?_lastused:labelDoug_30may25
        # 3          marketSquid_NCC         <NA> ?_lastused:labelDoug_30may25
        # 4            northAnch_NCC         <NA> ?_lastused:labelDoug_30may25
        # 5           Lmu_IntSprMayH         <NA>              ?_justFWeffects
        # 6      spottedRatfish_2025         <NA> ?_lastused:labelLisa_16sep25
        # 7     sablefish_WCGFS_2025         <NA> ?_lastused:labelLisa_21aug25
        # 8            predSealAbund         <NA>          ?_labelDoug_30may25
#ok to exclude        
        doug_guilds_Lisa %>% filter(shortName=="chubMack_NCC")
        doug_guilds_Lisa %>% filter(shortName=="jackMackerel_BRT")
        doug_guilds_Lisa %>% filter(shortName=="marketSquid_NCC")
        doug_guilds_Lisa %>% filter(shortName=="CaliforniaMarketSquid__2025")
        doug_guilds_Lisa %>% filter(shortName=="predSealAbund")
        doug_guilds_Lisa[which(doug_guilds_Lisa$shortName=="predSealAbund"),"replacementShortName"]<- "harbour_seal_DFO_WCVI"
          
        
#add back in 
doug_guilds_Lisa_start<-doug_guilds_Lisa
        doug_guilds_Lisa %>% filter(shortName=="sablefish_WCGFS_2025")
        doug_guilds_Lisa %>% filter(shortName=="predSealAbund")
        doug_guilds_Lisa %>% filter(shortName=="Lmu_IntSprMayH")
        
   
        doug_guilds_Lisa %>% filter(shortName=="sablefish_WCGFS_2025")
        doug_guilds_Lisa[which(doug_guilds_Lisa$shortName=="sablefish_WCGFS_2025"),
                         c("exclusionReason","node","guild_main","Lisa_15Oct2026")] <- 
                         c(NA,"X09","09.PredFishNCC","09.PredFishNCC")
        
        
        doug_guilds_Lisa %>% filter(shortName=="Lmu_IntSprJunH_2025")
        doug_guilds_Lisa %>% filter(shortName=="Lmu_IntSprMayH")
        doug_guilds_Lisa[which(doug_guilds_Lisa$shortName=="Lmu_IntSprMayH"),
                         c("exclusionReason","node","guild_main","Lisa_15Oct2026")] <- 
                          c(NA,"X06","06.Cond1NCC","06.Cond1NCC")
        
        
        
        doug_guilds_Lisa %>% filter(shortName=="predSealAbund")
        doug_guilds_Lisa %>% filter(shortName=="Lmu_IntSprMayH")
        doug_guilds_Lisa[which(doug_guilds_Lisa$shortName=="predSealAbund"),
                         c("exclusionReason","node","guild_main","Lisa_15Oct2026")] <- 
          c(NA,"X06","06.Cond1NCC","06.Cond1NCC")
        
        
#new data from PACEA        
        newrow<-doug_guilds_Lisa %>% filter(shortName=="predSealAbund")
        newrow[,c("pastGuild","latestGuild","guild_sub","exclusionReason","replacementShortName","shortName_09oct26")]<-NA
        newrow[,c("source_short","source_raw")]<-c("pacea_github_dfo","harbour_seals.rda")
        newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
          c("harbour_seal_DFO_WCVI","harbour_seal_DFO_WCVI","X11","11.PredMammalSmolt","11.PredMammalSmolt")
        newrow
        doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 

        newrow<-doug_guilds_Lisa %>% filter(shortName=="harbour_seal_DFO_WCVI")
         newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
          c("harbour_seal_DFO_Coastwide","harbour_seal_DFO_Coastwide","X11","11.PredMammalSmolt","11.PredMammalSmolt")
        newrow
        doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
        
        newrow<-doug_guilds_Lisa %>% filter(shortName=="sitkaHerring_EGoA")
        newrow[,c("dataset","subregion")]<-c("pacea_github_dfo","BCcoast")
        newrow[,c("source_short","source_raw")]<-c("pacea_github_dfo","herring_spawning_biomass.rda")
        newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026")]<-
          c("HerringSSB_BCcoast_wcvi_cc_prd","HerringSSB_BCcoast_wcvi_cc_prd","X13","13.FishPreyAK","13.FishPreyAK")
        newrow
        doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
        
        
        newrow<-doug_guilds_Lisa %>% filter(shortName=="sitkaHerring_EGoA")
        newrow[,c("dataset","subregion")]<-c("pacea_github_dfo","BCcoast")
        newrow[,c("source_short","source_raw")]<-c("pacea_github_dfo","bccm_surface_temperature_full_03.rds")
        newrow[,c("indicator","shortName","node","guild_main","Lisa_15Oct2026","trophicGroup")]<-
          c("sst_bccm_full_grid_mean_aprmayjun","sst_bccm_full_grid_mean_aprmayjun","X21","X21.Climate","X21.Climate","climate")
        newrow
        doug_guilds_Lisa<-rbind(doug_guilds_Lisa,newrow) 
        
        write.csv(doug_guilds_Lisa, "guilds_Lisa_10162026.csv", row.names=F)
        