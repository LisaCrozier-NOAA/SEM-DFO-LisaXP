

#I'm just using this script to get the current indicators
#I should have just read in "analyzeAKindices/guildsWithExclude.csv" from his directory and skipped this whole script

#Not necessary, but this script builds the final guildfile from indicators.csv, guilds.csv and our earlier exclude list.
#It explains all the indicators that were excluded, but these explanations and even decisions ARE NOT COMPLETE OR ALL CORRECT
#This script outputs combined_guildfiles in "analyzeAKindices/guildsWithExclude.csv" 





# this process is NOT COMPLETE:
# so far (6/25/2026) I have looked at all the indicators through X09, EXCLUDING X08!
# I made temporary files for each node, which are meant to be strung together at the end
# I did not look back through all of the notes. Need to verify sumPreyofPrey and sumPrey, and a few others marked with "?_"
# I left in the NAs in many cases even though that is not necessary to write out-- it is the default. But that way I could double check which shortNames were not mentioned in the list
# this all needs to be reviewed before final

# notes for Doug -- see if Murre are correlate w/ any of the WS productivity birds (01.ZooPreyNCC_WPbirds) -- if not, drop that group?
# Answer: commonMurre_JSOES_2026 is not correlated with any of the 01.ZooPreyNCC_WPbirds indicators

library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(janitor)


################Inf###########################################################################
# Input files

rootdir<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"

# Doug directory
# Change index to 1 for Lisa's path
path <- c("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results", "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/")[1]

workingDir <- file.path(path, "analyzeAKindices")
guildfile <- read.csv(file.path(workingDir, "guilds.csv"), stringsAsFactors=T, na.strings=c("", "NA", " "))
indicatorsfile <- read.csv(file.path(workingDir, "indicators.csv"))
excludeColumn <- read_xlsx(file.path(workingDir, "indicators LCAug25_exclude column.xlsx"))

outputFile <- file.path(workingDir, "guildsWithExclude.csv")

# Create output from scratch or update existing output file
createFromScratch <- TRUE

###########################################################################
# Functions
###########################################################################

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

###########################################################################
# Run
###########################################################################

setwd(workingDir)

# collect the current indicators file and the exclude column list-----------
indicatorsfile2 <- indicatorsfile %>% select(shortName, reasonExcluded, dataset, region, subregion)

# merge exclusion column with indicator file
indicatorsfile2 <- 
    left_join(indicatorsfile2, excludeColumn %>% select(shortName, ExplainExclusion), by="shortName") #%>%

# Create output from scratch or update existing output file
if(createFromScratch) {
    # Add a consistent guild name to each indicator------------
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
        arrange(guild, shortName)
    
    write.csv(guildfile_updated, file.path(path, "guildfile_updated.csv"), row.names=FALSE)
    
    guildfile <- guildfile_updated %>% 
        select(1:5)
    
    cat("Renaming", names(guildfile)[4], "to 'latestGuild'")
    names(guildfile)[4] <- "latestGuild"
    
    # merge indicators file with guild file
    indicatorsfile2 <- 
        left_join(indicatorsfile2, guildfile, by=c("shortName")) 
    names(indicatorsfile2)
} else {

    prevGuilds <- read.csv(outputFile)
    
    prevGuilds <- prevGuilds %>% rename(latestGuild=last_col()) %>% 
        rowwise() %>% mutate(guild=get_guild_value(c_across(pastGuild:last_col())))
    
    indicatorsfile2 <- left_join(indicatorsfile2, prevGuilds, by="shortName")
}

#Define exclusions 1 node at a time=============
#X01-----------
tmp <- indicatorsfile2 %>% select(indicator, shortName, guild, latestGuild) %>%
    arrange(guild, shortName)
# filter(is.na(latestGuild))
head(tmp)

tmp1 <- tmp %>% filter(str_starts(guild, "01.ZooPreyNCC"))
tmp1 <- tmp1 %>% mutate(SEMnode="PreyNCC") %>%
    
    
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

X01 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)
#X02=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "02.CompJellyNCC"))

tmp1 <- tmp1 %>% mutate(SEMnode="PreyNCC") %>%
    
    
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

X02 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X03=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "03.FishPreyNCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyNCC") %>%
    
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

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X03 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X04=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "04.CompFishNCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyNCC") %>%
    
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

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X04 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X05=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "05.ForageFishNCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyNCC") %>%
    
    mutate(data_updated=case_when(
        shortName=="northAnch_NCC"  ~ NA,
        shortName=="abundAnchovy"  ~ NA,
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="northAnch_NCC"  ~ "?_lastused:labelDoug_30may25",
        shortName=="abundAnchovy"  ~ "missingdata_shortNS",
        
        
        
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X05 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X06=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "06.Cond1NCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="Growth") %>%
    
    mutate(data_updated=case_when(
        shortName=="IGF_mu"  ~ "IGF_mu_2025",
        shortName=="Lmu_IntSprJunH"  ~ "Lmu_IntSprJunH_2025",
        shortName=="Lmu_IntSprMayH"  ~ NA,
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="Lmu_IntSprMayH"  ~ "?_justFWeffects", #double check this
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X06 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X07=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "07.Cond2NCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="Abundance") %>%
    
    mutate(data_updated=case_when(
        shortName=="cpue_IntSprJunH"  ~ "cpue_IntSprJunHW_2025",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="cpue_IntSprJunH"  ~ NA,
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X07 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X08=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "08.PredBird"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PredNCC") %>%
    
    mutate(data_updated=case_when(
        shortName=="commonMurre_NWFSC"  ~ "commonMurre_JSOES_2026",
        shortName=="sootyShear_NWFSC"  ~ "sootyShear_NWFSC_JSOES_2026",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        is.na(data_updated) ~ "notsalmonpred_or_grouped",
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X08 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X09=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "09.PredFishNCC"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PredNCC") %>%
    
    mutate(data_updated=case_when(
        shortName=="HakeAge5Plus"  ~ "HakeAge5Plus_2025",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="sablefish_WCGFS_2025"  ~ "?_lastused:labelLisa_21aug25",
        shortName=="spottedRatfish_2025"  ~ "?_lastused:labelLisa_16sep25",
        shortName=="synWCVIarrowtoothFloundMature" ~ "Failed completeness criterion",
        shortName=="synWCVIbigSkate" ~ "Failed completeness criterion",
        shortName=="synWCVIcanaryRockMature" ~ "Failed completeness criterion",
        shortName=="synWCVIcodMature" ~ "Failed completeness criterion",
        shortName=="synWCVIlingcodImmature" ~ "Failed completeness criterion",
        shortName=="synWCVIlongnoseSkate" ~ "Failed completeness criterion",
        shortName=="synWCVIspinyDogMature" ~ "Failed completeness criterion",
        shortName=="synWCVIyellowTailRockMature" ~ "Failed completeness criterion",
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X09 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X10=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "10.PredMammalNCC|BirdPred_AK|MammalPred"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PredAK") %>%
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
        shortName=="Killer.whales.BC_2yrLead_2025"  ~ "Killer.whales.NR.BC_2yrLead_2026",
        shortName=="Northern_f_s_WS_2026" ~ "Northern_f_s_2yrLead_WS_2026",
        shortName=="Californian_s_l_WS_2026" ~ "Californian_s_l_2yrLead_WS_2026",
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
        shortName=="Californian_s_l_WS_2026"  ~ "replaced with 2-year lead version: Californian_s_l_2yrLead_WS_2026",
        shortName=="Northern_f_s_WS_2026"  ~ "replaced with 2-year lead version: Northern_f_s_2yrLead_WS_2026",
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X10 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X11=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "11.PredMammalSmolt"))
tmp1[,2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PredNCC") %>%
    mutate(data_updated=case_when(
        shortName=="AllPinnDayBonn"  ~ "AllSeaLionsBonn_2yrLead_2025",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X11 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X12=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "12.ZooPreyAK|Prey_WCVI"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyAK") %>%
    mutate(data_updated=case_when(
        shortName=="pelag_CAI"  ~ "pollockBiomass_predAK_2026",
        shortName=="pelag_EAI"  ~ "pollockBiomass_predAK_2026",
        shortName=="pelag_WAI"  ~ "pollockBiomass_predAK_2026",
        
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="biomassAmphShelfFall"  ~ "timingmismatchforSpringrun",
        shortName=="biomassDecaShelfFall"  ~ "timingmismatchforSpringrun",
        shortName=="biomassEuphInletFall"  ~ "spatialmismatchforSpringrun",
        shortName=="biomassEuphInletSum"  ~ "spatialmismatchforSpringrun",
        shortName=="biomassEuphShelfFall"  ~ "timingmismatchforSpringrun",
        shortName=="biomassMysidShelfFall"  ~ "timingmismatchforSpringrun",
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X12 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X13=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "13.FishPreyAK|Prey_AK"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyAK") %>%
    mutate(data_updated=case_when(
        shortName=="pelag_CAI"  ~ "pollockBiomass_predAK_2026",
        shortName=="pelag_EAI"  ~ "pollockBiomass_predAK_2026",
        shortName=="pelag_WAI"  ~ "pollockBiomass_predAK_2026",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)


X13 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X14=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "14.CompAK"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PreyAK") %>%
    mutate(data_updated=case_when(
        shortName=="pinkSalmon"  ~ "pinkSalmon_2025",
        shortName=="pinkSalmonAsia"  ~ "pinkSalmonAsia_2025",
        shortName=="pinkSalmonNorthAmerica"  ~ "pinkSalmonNorthAmerica_2025",
        .default = NA)) %>%
    
    mutate(ExplainExclusion2=case_when(
        shortName=="Guadelupe_f_s_WS_2026"  ~ "notsalmonpred", #check if actually straggler?
        .default = NA))

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)

X14 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X15=================            

unique(tmp$guild)
tmp1 <- tmp %>% filter(str_starts(guild, "15.PredFishAK|Pred_AK"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="PredAK") %>%
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

tmp1 %>% select(shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% arrange(data_updated, ExplainExclusion2)

X15 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#X16================= 
# SAR
tmp1 <- tmp %>% filter(str_starts(guild, "SAR"))
tmp1[, 2:4]

tmp1 <- tmp1 %>% mutate(SEMnode="SAR", data_updated=NA, ExplainExclusion2=NA)

X16 <- tmp1 %>% select(indicator, shortName, data_updated, ExplainExclusion2, latestGuild, SEMnode, guild) %>% 
    arrange(data_updated, ExplainExclusion2)

#Combine all objects --still missing X08=========
object_names <- sprintf("X%02d", 1:16)

# 2. Grab the objects from the environment and bind them
combined_guildfiles <- bind_rows(mget(object_names))

combined_guildfiles <- combined_guildfiles %>% rename(replacementShortName=data_updated, exclusionReason=ExplainExclusion2, pastGuild=guild) %>% 
    select(indicator, shortName, SEMnode, exclusionReason, replacementShortName, pastGuild, latestGuild)

# Load the old file if we're updating it
if(!createFromScratch) {
    prevGuilds <- read.csv(outputFile)
    
    prevGuilds <- prevGuilds %>% select(indicator, shortName, pastGuild:last_col())
    
    prevGuilds$pastGuild <- NULL
    
    combined_guildfiles$latestGuild <- NULL

    combined_guildfiles <- left_join(combined_guildfiles, prevGuilds, by=c("indicator", "shortName"))
}
write.csv(combined_guildfiles, outputFile, row.names = F, quote=F)        

#Still to be clarified-------------         
tmp1 <- combined_guildfiles %>% filter(str_detect(exclusionReason,  fixed("?_")))
tmp1[, 2:4]

#          shortName data_updated            ExplainExclusion2
# 1             chubMack_NCC         <NA> ?_lastused:labelDoug_30may25
# 2         jackMackerel_BRT         <NA> ?_lastused:labelDoug_30may25
# 3          marketSquid_NCC         <NA> ?_lastused:labelDoug_30may25
# 4            northAnch_NCC         <NA> ?_lastused:labelDoug_30may25
# 5           Lmu_IntSprMayH         <NA>              ?_justFWeffects
# 6      spottedRatfish_2025         <NA> ?_lastused:labelLisa_16sep25
# 7     sablefish_WCGFS_2025         <NA> ?_lastused:labelLisa_21aug25
# 8            predSealAbund         <NA>          ?_labelDoug_30may25

# Notes:
# chubMack_NCC, jackMackerel_BRT, marketSquid_NCC, northAnch_NCC, predSealAbund, and Lmu_IntSprMayH: these were all excluded in 
#   "LisaDataProcessScripts 2025/guilds.Lisa.07152025.csv," which was committed to the repo by Lisa on July 15, 2025 
#   (documenting guild choices and datasets for group review and Doug)
# spottedRatfish_2025: excluded in "guilds simple LC Sep 24 2025.csv", which was committed to repo by Lisa on Sep. 25, 2025 (guilds and indicators messy)
# sablefish_WCGFS_2025: looks like this was inadvertently excluded because "LisaDataProcessScripts 2025/guilds.Lisa.07152025.csv," which
#   was committed to repo by Lisa on 08Sep25 (changed some birds from predators to zooplankton indicators because they don't eat salmon) was based on an
#   old guilds.csv file that used a different shortName => check with Lisa to see if this should be reincorporated


