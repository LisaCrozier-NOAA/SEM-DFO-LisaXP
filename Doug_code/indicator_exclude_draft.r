#this process is NOT COMPLETE:
#so far (6/25/2026) I have looked at all the indicators through X09, EXCLUDING X08!
#I made temporary files for each node, which are meant to be strung together at the end
#I did not look back through all of the notes. Need to verify sumPreyofPrey and sumPrey, and a few others marked with "?_"
#I left in the NAs in many cases even though that is not necessary to write out-- it is the default. But that way I could double check which shortNames were not mentioned in the list
#this all needs to be reviewed before final

#notes for Doug -- see if Murre are correlate w/ any of the WS productivity birds (01.ZooPreyNCC_WPbirds) -- if not, drop that group?

library(dplyr)
library(tidyr)
library(readxl)
library(stringr)

#Input files
path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/"
    guildfile<-read.csv(paste0(path,"analyzeAKindices/guilds.csv"),stringsAsFactors = T,na.strings = c("", "NA", " "))
    indicatorsfile<-read.csv(paste0(path,"analyzeAKindices/indicators.csv"))
    exclude.column<-read_xlsx(paste0(path,"analyzeAKindices/indicators LCAug25_exclude column.xlsx"))
    
#Add a consistent guild name to each indicator------------

          
          # Helper function 
          get_guild_value <- function(x) {
            # Remove NA values from the row vector
            non_nas <- na.omit(x)
            
            # If the whole row was NA, return NA
            if (length(non_nas) == 0) return(NA_character_)
            
            # If all non-NA values are the same, or if there's only one, return that value.
            # Otherwise, return the last non-NA value assigned.
            if (length(unique(non_nas)) == 1) {
              return(non_nas[1])
            } else {
              return(tail(non_nas, 1))
            }
          }
          
          # Apply to your dataframe
          guildfile_updated <- guildfile %>%
            rowwise() %>%
            mutate(
              guild = get_guild_value(c_across(labelLisa:last_col(1)))
            ) %>%
            ungroup() %>%
            # Move the new 'guild' column to position 3 (before labelLisa)
            relocate(guild, .before = labelLisa) %>% 
            # 2. Now that guild has shifted, the last two columns are back to being "last"
            relocate(last_col(1):last_col(), .after = guild) %>%
            arrange(guild,shortName)
          
          
       #    write.csv(guildfile_updated,paste0(path,"guildfile_updated.csv"),row.names=FALSE)
          
          guildfile<-guildfile_updated %>% 
            select(1:5)

#collect the current indicators file and the exclude column list-----------
          indicatorsfile2<-indicatorsfile %>% select(shortName,reasonExcluded,dataset,region,subregion)
          
          #merge exclusion column with indicator file
          indicatorsfile2<- 
            left_join(indicatorsfile2,exclude.column %>% select(shortName,ExplainExclusion),by="shortName") #%>%
          
          #merge indicators file with guild file
          indicatorsfile2<- 
            left_join(indicatorsfile2,guildfile,by=c("shortName")) 
          
          names(indicatorsfile2)

          
#Define exclusions 1 node at a time=============
#X01-----------
          tmp<- indicatorsfile2 %>% select(indicator,shortName,guild,Doug_24jun26,ExplainExclusion,reasonExcluded,SEMnode.for.AIC.analysis,dataset,region,subregion) %>%
            arrange(guild,shortName)
          # filter(is.na(Doug_24jun26))
          head(tmp)
          
          tmp1<-tmp %>% filter(guild %in% c("01.ZooPreyNCC","01.ZooPreyNCC_NHL","01.ZooPreyNCC_JSOES","01.ZooPreyNCC_WPbirds","01.ZooPreyNCCb","PreyNCC"))
          
          tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
            
          
          mutate(data_updated=case_when(
            shortName=="CanMagM_planktonJun_2025b"  ~ "sumPrey_planktonJun_2026",
            shortName=="CanMagM_planktonJuneNCC"  ~ "sumPrey_planktonJun_2026",
            shortName=="CrabMegalopae"  ~ "sumPrey_planktonJun_2026",
            shortName=="CrabMegalopae_2025"  ~ "sumPrey_planktonJun_2026",
            shortName=="DecapodAll_planktonJun_2025b"  ~ "sumPrey_planktonJun_2026",
            shortName=="DecapodAll_planktonJuneNCC"  ~ "sumPrey_planktonJun_2026",
            shortName=="Fish_planktonJuneNCC"  ~ "sumPrey_planktonJun_2026",
            shortName=="Hyperiid_planktonJuneNCC"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="Insect_planktonJuneNCC"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="LarCop_planktonJuneNCC"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="Pteropod_planktonJuneNCC"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="PlanktonJuneNCC"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="PlanktonJuneNCC_Burke_2025"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="PlanktonJuneNCC_IGF"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="PlanktonJuneNCC_Ruz1_2025"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="PlanktonJuneNCC_Ruz2_2025"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="sumAll_planktonJun_2025b"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="sumPrey_planktonJun_2025b"  ~ "sumPreyOfPrey_planktonJun_2026",
            shortName=="sumPreyOfPrey_planktonJun_2025b"  ~ "sumPreyOfPrey_planktonJun_2026",
            
            shortName=="sumPrey_planktonJun_2026"  ~ NA,
            shortName=="sumPreyOfPrey_planktonJun_2026"  ~ NA,
            

            shortName=="copeBiomassNH05north"  ~ "copeLogBiomassNorth_2026",
            shortName=="copeBiomassNH05south"  ~ "copeLogBiomassSouth_2026",
            shortName=="copeBiomassNorth_2025b"  ~ "copeLogBiomassNorth_2026",
            shortName=="copeBiomassSouth_2025b"  ~ "copeLogBiomassSouth_2026",
            shortName=="copeLogBiomassNorth_2026"  ~ NA,
            shortName=="copeLogBiomassSouth_2026"  ~ NA,
            shortName=="copeCIaxis1"  ~ NA,
            shortName=="copeCIaxis2"  ~ NA,
            shortName=="copeCIaxis1_2025b"  ~ NA,
            shortName=="copeCIaxis2_2025b"  ~ NA,
            
            shortName=="krill_NHL_2025"  ~ "Tspin_NHL_2026",
            shortName=="muCanMagM_NHLbiomass"  ~ "NHLlogSum_win_05_2026",
            shortName=="muDecapodAll_NHLbiomass"  ~ "NHLlogSum_win_05_2026",
            shortName=="muFish_NHLbiomass"  ~ "NHLlogSum_win_05_2026",
            shortName=="muLarCop_NHLbiomass"  ~ "NHLlogSum_win_05_2026",
            shortName=="muPteropod_NHLbiomass"  ~ "NHLlogSum_win_05_2026",
            shortName=="NHLsum_win_05_2025b"  ~ "NHLlogSum_win_05_2026",
            shortName=="NHLsum_win_15_2025b"  ~ "NHLlogSum_win_15_2026",
            shortName=="NHLsum_win_25_2025b"  ~ "NHLlogSum_win_25_2026",

            shortName=="muCanMagM_NHLdensity"  ~ "NHLlogSum_win_05_2026",
            shortName=="muDecapodAll_NHLdensity"  ~ "NHLlogSum_win_05_2026",
            shortName=="muFish_NHLdensity"  ~ "NHLlogSum_win_05_2026",
            shortName=="muLarCop_NHLdensity"  ~ "NHLlogSum_win_05_2026",
            shortName=="muPteropod_NHLdensity"  ~ "NHLlogSum_win_05_2026",
            
            
            shortName=="Epac_NHL_2026"  ~ NA,
            shortName=="Tspin_NHL_2026"  ~ NA,
            shortName=="NHLlogSum_win_05_2026"  ~ NA,
            shortName=="NHLlogSum_win_15_2026"  ~ NA,
            shortName=="NHLlogSum_win_25_2026"  ~ NA,

            shortName=="cassAuk_JSOES_2026"  ~ NA,
            shortName=="cassAuk_NWFSC"  ~ "cassAuk_JSOES_2026",
            shortName=="Cassin_s_aukl_WS_2026"  ~ NA,
            shortName=="habCompInd"  ~ NA,
            shortName=="Parakeet_aukl_WS_2026"  ~ NA,
            
            .default = NA)) %>%
            
           mutate(ExplainExclusion2=case_when(
            shortName=="copeCIaxis1"  ~ "selected_copeLogBiomassNorth_2026",
            shortName=="copeCIaxis2"  ~ "selected_copeLogBiomassSouth_2026",

            shortName=="Parakeet_aukl_WS_2026"  ~ "straggler_NS",
            shortName=="Red_necked_phal_WS_2026"  ~ "straggler_NS",
            shortName=="Red_phal_WS_2026"  ~ "straggler_NS",
            shortName=="Sabine_gull_WS_2026"  ~ "straggler_NS",
 
            shortName=="crabLandingWA"  ~ "switched_to_JSOES_sumPrey_planktonJun_2026",
            
            .default = NA))
            
         X01<- tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
#X02=================            

         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("02.CompJellyNCC"))
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
           
           
            mutate(data_updated=case_when(
                shortName=="EggyolkJelly__2025"  ~ NA,
                shortName=="MoonJelly__2025"  ~ NA,
                shortName=="PyrosomaAtlanticum__2025"  ~ NA,
                shortName=="SeaNettle__2025"  ~ NA,
                shortName=="WaterJelly__2025"  ~ NA,
                shortName=="eggJel_NCC"  ~ "EggyolkJelly__2025",
                shortName=="moonJel_NCC"  ~ "MoonJelly__2025",
                shortName=="seaNet_NCC"  ~ "SeaNettle__2025",
                shortName=="waterJel_NCC"  ~ "WaterJelly__2025",
             .default = NA)) %>%
                 
            mutate(ExplainExclusion2=case_when(
                shortName=="EggyolkJelly__2025"  ~ "Daly:diff diet",
                shortName=="WaterJelly__2025"  ~ "Daly:diff diet",
                
                .default = NA))
         
         X02<- 
           tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
#X03=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("03.FishPreyNCC","03.FishPreyNCC_0"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%

           mutate(data_updated=case_when(
             shortName=="HakeAge1_2025"  ~ NA,
             shortName=="IchNear_2025"  ~ "IchNear_2026",
             shortName=="IchNear_2026"  ~ NA,
             shortName=="comMurreDietFlat_Yaquina_NCC"  ~ NA,
             shortName=="comMurreDietHerrSard_Yaquina_NCC"  ~ NA,
             shortName=="comMurreDietRock_Yaquina_NCC"  ~ NA,
             shortName=="comMurreDietSLance_Yaquina_NCC"  ~ NA,
             shortName=="comMurreDietSmelt_Yaquina_NCC"  ~ NA,
             shortName=="hake_NCC"  ~ "yoyHake_JSOES_NCC",
             shortName=="lagAbundAnchovy"  ~ NA,
             shortName=="lagAbundHerring"  ~ NA,
             shortName=="lagAbundSardine"  ~ NA,
             shortName=="rhinoAukDietAnch_NCC"  ~ NA,
             shortName=="rhinoAukDietHerr_NCC"  ~ NA,
             shortName=="rhinoAukDietRock_NCC"  ~ NA,
             shortName=="rhinoAukDietSLance_NCC"  ~ NA,
             shortName=="rhinoAukDietSmelt_NCC"  ~ NA,
             shortName=="surfSmelt_NCC"  ~ NA,
             shortName=="whiteBaitSmelt_NCC"  ~ NA,
             shortName=="yoyHake_JSOES_NCC"  ~ NA,
             shortName=="yoyRock_JSOES_NCC"  ~ NA,
             shortName=="IchNear"  ~ "IchNear_2026",
             
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="hake_NCC"  ~ "Daly:toolate",
             shortName=="yoyHake_JSOES_NCC"  ~ "Daly:toolate",
             shortName=="comMurreDietFlat_Yaquina_NCC"  ~ "missingdata_shortNS",
#             shortName=="lagAbundAnchovy"  ~ "missingdata_shortNS",   #kept bec in DFA
             shortName=="lagAbundHerring"  ~ "missingdata_shortNS",
             shortName=="lagAbundSardine"  ~ "missingdata_shortNS",
             shortName=="rhinoAukDietAnch_NCC"  ~ "missingdata_shortNS",
             shortName=="rhinoAukDietHerr_NCC"  ~ "missingdata_shortNS",
#             shortName=="rhinoAukDietRock_NCC"  ~ "missingdata_shortNS",  #kept bec in DFA
             shortName=="rhinoAukDietSLance_NCC"  ~ "missingdata_shortNS",
             shortName=="rhinoAukDietSmelt_NCC"  ~ "missingdata_shortNS",
             
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X03<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)

#X04=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("04.CompFishNCC","04.CompFishNCC_0"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="CaliforniaMarketSquid__2025"  ~ NA,
             shortName=="PacificPompano__2025"  ~ NA,
             shortName=="Sablefish__JSOES_2025"  ~ NA,
             shortName=="chubMack_NCC"  ~ NA,
             shortName=="juvChum_NCC"  ~ NA,
             shortName=="marketSquid_NCC"  ~ NA,
             shortName=="marketsquid_GAM_2025"  ~ NA,
             shortName=="pompano_NCC"  ~ NA,
             shortName=="HakeAge1"  ~ NA,
            
             
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="chubMack_NCC"  ~ "?_lastused:labelDoug_30may25",
             shortName=="jackMackerel_BRT"  ~ "?_lastused:labelDoug_30may25",
             shortName=="juvChum_NCC"  ~ "Daly:notcompetitor",
             shortName=="marketSquid_NCC"  ~ "?_lastused:labelDoug_30may25",
             shortName=="HakeAge1"  ~ "straggler_NS",
             shortName=="pompano_NCC"  ~ "straggler_NS",

             
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X04<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X05=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("05.ForageFishNCC","05.ForageFishNCC_0"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="northAnch_NCC"  ~ NA,
             shortName=="abundAnchovy"  ~ NA,
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="northAnch_NCC"  ~ "?_lastused:labelDoug_30may25",
             shortName=="abundAnchovy"  ~ "missingdata_shortNS",
             
             
             
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X05<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
         
#X06=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("06.Cond1NCC"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="IGF_mu"  ~ "IGF_mu_2025",
             shortName=="Lmu_IntSprJunH"  ~ "Lmu_IntSprJunH_2025",
             shortName=="Lmu_IntSprMayH"  ~ NA,
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="Lmu_IntSprMayH"  ~ "?_justFWeffects", #double check this
             
             
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X06<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)

         
#X07=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("07.Cond2NCC"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="cpue_IntSprJunH"  ~ "cpue_IntSprJunHW_2025",
             .default = NA)) %>%
           
            mutate(ExplainExclusion2=case_when(
              shortName=="cpue_IntSprJunH"  ~ NA,
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X07<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
#X08=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(str_detect(guild, "08\\.PredBird"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PredNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="commonMurre_NWFSC"  ~ "commonMurre_JSOES_2026",
             shortName=="sootyShear_NWFSC"  ~ "sootyShear_NWFSC_JSOES_2026",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             is.na(data_updated) ~ "notsalmonpred_or_grouped",
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X08<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X09=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("09.PredFishNCC"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PredNCC") %>%
           
           mutate(data_updated=case_when(
             shortName=="HakeAge5Plus"  ~ "HakeAge5Plus_2025",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="sablefish_WCGFS_2025"  ~ "?_lastused:labelLisa_21aug25",
             shortName=="spottedRatfish_2025"  ~ "?_lastused:labelLisa_16sep25",
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X09<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
#X10=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("10.PredMammalNCC","10.PredMammalNCC_0","10.PredMammalNCC_1","BirdPred_AK","MammalPred","Pred_AK"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PredAK") %>%
           mutate(data_updated=case_when(
             shortName=="AllPinnDayBonn"  ~ "AllSeaLionsBonn_2yrLead_2025",
             shortName=="AllSeaLionsBonn_2025"  ~ "AllSeaLionsBonn_2yrLead_2025",
             shortName=="AllSeaLionsEMB_2025"  ~ "AllSeaLionsEMB_2yrLead_2025",
             shortName=="BonnPinn"  ~ "AllSeaLionsBonn_2yrLead_2025",
             shortName=="CSL_Bonn"  ~ "AllSeaLionsBonn_2yrLead_2025",
             shortName=="CSL_EMB"  ~ "AllSeaLionsEMB_2yrLead_2025",
             shortName=="CSLspring_EastMooringBasin"  ~ "AllSeaLionsEMB_2yrLead_2025",
             shortName=="Harbor_seal_CR"  ~ "Harbor_seal_CR_2yrLead_2025",
             shortName=="Killer_whales_BC"  ~ "Killer.whales.BC_2yrLead_2025",
             shortName=="Killer_whales_BC_2025"  ~ "Killer.whales.BC_2yrLead_2025",
             shortName=="Killer.whales.BC_2yrLead_2025"  ~ NA,
             shortName=="SSL_Bonn"  ~ "AllSeaLionsBonn_2yrLead_2025",
             shortName=="SSL_EMB"  ~ "AllSeaLionsEMB_2yrLead_2025",
             shortName=="Steller_s_l_WS_2026"  ~ "Steller_s_l_2yrLead_WS_2026",
             shortName=="ssl_est_wholerange"  ~ "ssl.est.wholerange_2yrLead",
             
             shortName=="SSL_CAI"  ~ "ssl_est_wholerange",
             shortName=="SSL_EAI"  ~ "ssl_est_wholerange",
             shortName=="SSL_WAI"  ~ "ssl_est_wholerange",
             shortName=="SSL_EGoA"  ~ "ssl_est_wholerange",
             shortName=="SSL_WGoA"  ~ "ssl_est_wholerange",
             shortName=="sealCount"  ~ "predSealAbund",

             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="kittiwake_WGoA"  ~ "notsalmonpredator",
             shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
             shortName=="predSealAbund"  ~ "?_labelDoug_30may25",
             shortName=="Californian_s_l_WS_2026"  ~ "?_lastused:",
             shortName=="Northern_f_s_WS_2026"  ~ "?_lastused:shiftDoug_07apr26",
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X10<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X11=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("11.PredMammalSmolt"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PredAK") %>%
           mutate(data_updated=case_when(
             shortName=="AllPinnDayBonn"  ~ "AllSeaLionsBonn_2yrLead_2025",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X11<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X12=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("12.ZooPreyAK","Prey_WCVI"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyAK") %>%
           mutate(data_updated=case_when(
             shortName=="pelag_CAI"  ~ "pollockBiomass_predAK_2026",
             shortName=="pelag_EAI"  ~ "pollockBiomass_predAK_2026",
             shortName=="pelag_WAI"  ~ "pollockBiomass_predAK_2026",
             
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="biomassAmphShelfFall"  ~ "timingmismatchforSpringrun",
             shortName=="biomassDecaShelfFall"  ~ "timingmismatchforSpringrun",
             shortName=="biomassEuphInletFall"  ~ "timingmismatchforSpringrun",
             shortName=="biomassEuphInletSum"  ~ "timingmismatchforSpringrun",
             shortName=="biomassEuphShelfFall"  ~ "timingmismatchforSpringrun",
             shortName=="biomassMysidShelfFall"  ~ "timingmismatchforSpringrun",
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X12<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X13=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("13.FishPreyAK","Prey_AK"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyAK") %>%
           mutate(data_updated=case_when(
             shortName=="pelag_CAI"  ~ "pollockBiomass_predAK_2026",
             shortName=="pelag_EAI"  ~ "pollockBiomass_predAK_2026",
             shortName=="pelag_WAI"  ~ "pollockBiomass_predAK_2026",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X13<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
         
#X14=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("14.CompAK","14.CompAK_0","14.CompAK_1"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyAK") %>%
           mutate(data_updated=case_when(
             shortName=="pinkSalmon"  ~ "pinkSalmon_2025",
             shortName=="pinkSalmonAsia"  ~ "pinkSalmonAsia_2025",
             shortName=="pinkSalmonNorthAmerica"  ~ "pinkSalmonNorthAmerica_2025",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X14<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)
         
         
#X15=================            
         
         unique(tmp$guild)
         tmp1<-tmp %>% filter(guild %in% c("15.PredFishAK","Pred_AK"))
         tmp1[,2:4]
         
         tmp1<-tmp1 %>% mutate(SEMnode="PreyAK") %>%
           mutate(data_updated=case_when(
             shortName=="ApexBio_CAI"  ~ "ArrowtoothFlounderBiomass_predAK_2026",
             shortName=="ApexBio_EAI"  ~ "ArrowtoothFlounderBiomass_predAK_2026",
             shortName=="ApexBio_WAI"  ~ "ArrowtoothFlounderBiomass_predAK_2026",
             shortName=="ApexFish_EGoA"  ~ "PacificCodBiomass_predAK_2026",
             shortName=="ApexFish_WGoA"  ~ "PacificCodBiomass_predAK_2026",
             .default = NA)) %>%
           
           mutate(ExplainExclusion2=case_when(
             shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
             .default = NA))
         
         tmp1 %>% select(shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild) %>% arrange(data_updated,ExplainExclusion2)
         
         
         X15<-tmp1 %>% select(indicator,shortName,data_updated,ExplainExclusion2,Doug_24jun26,SEMnode,guild,ExplainExclusion,reasonExcluded) %>% arrange(data_updated,ExplainExclusion2)

         
#Combine all objects --still missing X08=========
         object_names <- sprintf("X%02d", 1:15)
         
         # 2. Grab the objects from the environment and bind them
         combined_guildfiles <- bind_rows(mget(object_names))
         
         write.csv(combined_guildfiles,paste0(path,"analyzeAKindices/guildfile.excludecol.csv"),row.names = FALSE)        
         
         
         
#Still to be clarified-------------         
         tmp1<-combined_guildfiles %>% filter(str_detect(ExplainExclusion2, fixed("?_")))
         tmp1[,2:4]

         #          shortName data_updated            ExplainExclusion2
         # 1             chubMack_NCC         <NA> ?_lastused:labelDoug_30may25
         # 2         jackMackerel_BRT         <NA> ?_lastused:labelDoug_30may25
         # 3          marketSquid_NCC         <NA> ?_lastused:labelDoug_30may25
         # 4            northAnch_NCC         <NA> ?_lastused:labelDoug_30may25
         # 5           Lmu_IntSprMayH         <NA>              ?_justFWeffects
         # 6      spottedRatfish_2025         <NA> ?_lastused:labelLisa_16sep25
         # 7     sablefish_WCGFS_2025         <NA> ?_lastused:labelLisa_21aug25
         # 8            predSealAbund         <NA>          ?_labelDoug_30may25
         # 9  Californian_s_l_WS_2026         <NA>                  ?_lastused:
         # 10    Northern_f_s_WS_2026         <NA> ?_lastused:shiftDoug_07apr26         
                  
