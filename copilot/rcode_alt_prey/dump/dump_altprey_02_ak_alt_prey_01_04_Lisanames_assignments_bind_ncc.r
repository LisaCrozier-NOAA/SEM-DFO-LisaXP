

nccpred_dfa_prey_mapped <- read.csv(file.path(outputDir, "nccpred_dfa_prey_mapped.csv"), row.names = NULL)
print(as_tibble(nccpred_dfa_prey_mapped %>% select(constituent_predators, altprey1, altprey2,altprey1_data_col,  altprey2_data_col)), n = Inf)


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
#1. decide on altprey1 & altprey2 for all ak predators--------
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

clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))
sem_data<-read.csv("data_Lisa/sem_master_data.csv",row.names = NULL)
  
    
library(dplyr)

# Look at all Alaska fish predators (Guild 15) and assign preferred alternate prey
akpred <- combined_guildfiles %>% 
  filter(grepl("15", latest_guild)) %>% 
  select(indicator, short_name_lower, sem_name, latest_guild) %>%
  mutate(
    # --- ALTPREY 1: Pollock first priority for most, Krill for Spiny Dogfish ---
    altprey1 = case_when(
      # Spiny Dogfish: Krill primary
      grepl("spinydogfish", indicator) ~ "krill",
      
      # Groundfish, Sharks (including Salmon Shark), and Demersals: Pollock primary
      .default = "pollock"
    ),
    
    # --- ALTPREY 2: Herring for groundfish/dogfish, Pink Salmon for Salmon Sharks ---
    altprey2 = case_when(
      # Salmon Sharks: Secondary target is Pink Salmon / Salmonids
      grepl("salmonshark", indicator) ~ "pink_salmon",
      
      # Default for groundfish, sleepers, cod, halibut, sablefish, dogfish
      .default = "herring"
    ),
    
    # --- ALTPREY 3: Capelin as third choice ---
    altprey3 = case_when(
      # Sablefish: Tertiary switch to krill
      grepl("sablefish", indicator) ~ "krill",
      
      # Capelin retained only for small-pelagic/generalist foragers
      grepl("pacificcod|arrowtooth|spinydogfish|sharkcatch", indicator) ~ "capelin",
      
      # Set to NA for large apex/benthic predators (Halibut, Sleeper Shark, Salmon Shark)
      .default = NA_character_
      )
  ) %>% 
  arrange(altprey1, altprey2, altprey3)

#add ssl row
      akpred<-rbind(akpred,c("ssl_seak_pup_pred","ssl_seak_pup_pred","predak",NA_character_,NA_character_,"herring","capelin"))
      akpred
 
# Save updated Alaska prey assignments --------
      write.csv(akpred, file.path(outputDir, "akpred_altprey_assigned.csv"), row.names = FALSE)


#get akpred data-----
akpred_shortnames <- combined_guildfiles %>% filter(grepl("15", latest_guild)) %>% pull(short_name_lower)
akpred_shortnames
# grep("pollock",names(sem_data),value=T)
names(sem_data)

#ssl pup and sem_data
      ssl<-read.csv("copilot/outputs_2/ssl.dat.csv",row.names = NULL)  %>%
        select(Year=year,X10_ssl_seak_pup_pred=ssl_seak_pup_pred) 
      akpred_from_sem_data<-read.csv("data_Lisa/sem_master_data.csv",row.names = NULL) %>%
        rename(Year=year) %>% 
        filter(Year>=1996) %>% 
        select(Year,
               X15_ArrowtoothFlounderBiomass_predAK,X15_sablefishBiomass_predAK,
               X15_sablefishRecruitment_predAK, X15_PacificCodBiomass_predAK,
               X15_spinyDogfishGoA_predAK,X15_halibutBiomassAge8plus_2yrLead_predAK,
               X15_DFA_sleeperSharkBSAI_predAK,
               X15_spinyDogfishBSAI_predAK, X15_salmonSharkBSAI_predAK,
               X15_salmonSharkGoA_predAK,X15_sharkCatchGoA_predAK)

  
  akpred_data<- left_join(akpred_from_sem_data,ssl, by = "Year")

  # Save akpred_data --------
  write.csv(akpred_data, file.path(outputDir, "akpred_data.csv"), row.names = FALSE)
  
names(akpred_data)

      # Display complete output table
      print(as_tibble(akpred), n = Inf)
      

      

# 2. Get prey data----------
      #unique prey needed-------
      unique(c(akpred$altprey1,akpred$altprey2,akpred$altprey3))
      # "krill"       "pollock"   "herring"     "pink_salmon" "capelin"   
      #  X13_pollock_age1plus (sem_data: pollock), X12_egoa.krill (krill), X13_stka_herr_matbiom (herring), X14_pinkSalmonNorthAmerica (pink_salmon),X13_mid_il_capelin (capelin)
      
      
      #look at options
      rankedIndicators <- read.csv(file.path(outputDir, "rankedIndicators.csv"), row.names = NULL) %>%
        mutate(dfa_cols=paste0("X",DFAname))%>% 
        relocate(dfa_cols,.after=DFAname) %>%
        filter(grepl("12|13|14|15",DFAname)) %>%
        select(dfa_cols,rankedIndicators)
      
      rankedIndicators
      
      #haven't included EGOA into Doug's process yet
      names(sem_data)
      akprey_from_sem_data<-read.csv("data_Lisa/sem_master_data.csv",row.names = NULL)%>%
        rename(Year=year) %>% 
        filter(Year>=1996) %>% 
        select(Year,X13_pollock_age1plus=pollockBiomassAIage1plus_predAK_2026,X14_pinkSalmonNorthAmerica) ;names(akprey_from_sem_data)
      
      akprey_from_ak_yr_data<-read.csv("data_Lisa/ak_yr.csv",row.names = NULL) %>%
        rename(Year=year) %>% 
        filter(Year>=1996) %>% 
        select(Year,X12_egoa.krill=secm_euph_dens,X13_mid_il_capelin=mid_il_capelin, X13_stka_herr_matbiom=stka_herr_matbiom) ;names(akprey_from_ak_yr_data)
      
#      Construct wide prey dataset with exact raw source column names and 2-year leads
      akprey_data <- left_join(akprey_from_sem_data,akprey_from_ak_yr_data , by = "Year")

      ak_data <- left_join(akpred_data,akprey_data , by = "Year")
      
      names(ak_data)
      print(as_tibble(ak_data),n=Inf)

      
#Save ak pred and prey data------
      write.csv(ak_data, file.path(outputDir, "ak_data_1998_2021.csv"), row.names = FALSE)

      #X15_spinyDogfishBSAI_predAK & X15_salmonSharkBSAI_predAK are short time series, the rest are long
      
      
# # 3. Reduce pred list to unique dfas --------
      #I don't think we need this anymore
# # Add dfa names for predators, check for conflicts in alt prey
#       akpred<-read.csv( file.path(outputDir, "akpred_altprey_assigned.csv"), row.names = NULL)
#       nrow(akpred)
#       # 3. Match individual indicators w/ dfa names via grepl (case-insensitive substring search inside the collapsed string)
#       akpred_dfa_matched <- akpred %>%
#         # Join all combinations to test each short_name against each collapsed DFA string
#         left_join(rankedIndicators, by = character(), relationship = "many-to-many") %>%
#         # Keep rows where short_name_lower is found inside the raw collapsed string
#         # (or keep as NA if an indicator isn't part of any DFA in rankedIndicators)
#         filter(
#           is.na(rankedIndicators) | 
#             mapply(function(pattern, text) grepl(pattern, text, ignore.case = TRUE), 
#                    short_name_lower, rankedIndicators)
#         ) %>%
#         relocate(dfa_cols,.before=latest_guild)
#       
#       akpred_dfa_matched
#       nrow(akpred_dfa_matched)
#       
#       akpred_dfa_matched$short_name_lower
#       akpred$short_name_lower
#       
#       prey_conflicts <- akpred_dfa_matched %>%
#         group_by(dfa_cols) %>%
#         summarise(
#           num_predators       = n(),
#           predators_in_dfa    = paste(unique(short_name_lower), collapse = ", "),
#           altprey1_choices    = paste(unique(altprey1[!is.na(altprey1)]), collapse = " | "),
#           altprey2_choices    = paste(unique(altprey2[!is.na(altprey2)]), collapse = " | "),
#           altprey3_choices    = paste(unique(altprey3[!is.na(altprey2)]), collapse = " | "),
#           altprey1_conflict   = n_distinct(altprey1[!is.na(altprey1)]) > 1,
#           altprey2_conflict   = n_distinct(altprey2[!is.na(altprey2)]) > 1,
#           altprey3_conflict   = n_distinct(altprey3[!is.na(altprey2)]) > 1,
#           .groups = "drop"
#         ) %>%
#         filter(altprey1_conflict | altprey2_conflict | altprey3_conflict)
#       
#       prey_conflicts
#       prey_conflicts %>% select(dfa_cols,num_predators,altprey1_choices,altprey2_choices,altprey3_choices,predators_in_dfa)
#       
#       
#       # Collapse to Unique DFA Lines---------
#       akpred_dfa_unique <- akpred_dfa_matched %>%
#         filter(!is.na(dfa_cols)) %>%
#         mutate(is_2yrlead="FALSE") %>%
#         group_by(dfa_cols) %>%
#         summarise(
#           # Aggregate all constituent short names for reference
#           constituent_predators = paste(unique(short_name_lower), collapse = " | "),
#           altprey1              = unique(altprey1[!is.na(altprey1)])[1],
#           altprey2              = unique(altprey2[!is.na(altprey2)])[1],
#           altprey3              = unique(altprey3[!is.na(altprey3)])[1],
#           is_2yrlead              = unique(is_2yrlead[!is.na(is_2yrlead)])[1],
#           .groups = "drop"
#         )
#       
#       # View the final unconflicted DFA table
#       print(as_tibble(akpred_dfa_unique), n = Inf)
#       print(as_tibble(nccpred_dfa_unique), n = Inf)
#       
#       akpred_dfa_unique$dfa_cols
#       nrow(akpred_dfa_unique)
#       #save dfa prey assignments----------
#       write.csv(akpred_dfa_unique, file.path(outputDir, "akpred_dfa_altprey_assigned.csv"), row.names = FALSE)
# 
#       
      
 # ---------------------------------------------------------------------------
 #4. Map data source column names into akpred_dfa_unique-------------------------
 # ---------------------------------------------------------------------------
 
      names(nccpred_dfa_matched_sem_names)
      names(akpred_dfa_prey_mapped)
      names()
      akpred_dfa_prey_mapped %>% 
        rename()
      
      akpred <- read.csv(file.path(outputDir, "akpred_altprey_assigned.csv"), row.names = NULL)
      ak_data <- read.csv(file.path(outputDir, "ak_data_1998_2021.csv"), row.names = NULL)
      names(ak_data)
#      akpred_dfa_unique <- read.csv(file.path(outputDir, "akpred_dfa_altprey_assigned.csv"), row.names = NULL)
      
      #  X12_egoa.krill (krill), X13.FishPreyAK_DFA1 (pollock), X13_stka_herr_matbiom (herring), X14.CompAK_1_smoothed (pink_salmon),X13_mid_il_capelin (capelin)
      
      akpred_dfa_prey_mapped <- akpred %>%
        mutate(
          # Map exact source column name for altprey1
          altprey1_data_col = case_when(
            altprey1 == "pollock"  ~ "X13_pollock_age1plus",
            altprey1 == "krill"   ~ "X12_egoa.krill",
            TRUE ~ NA_character_
          ),
          
          # Map exact source column name for altprey2
          altprey2_data_col = case_when(
            altprey2 == "herring"   ~ "X13_stka_herr_matbiom",
            altprey2 == "pink_salmon"  ~ "X14_pinkSalmonNorthAmerica",
            TRUE ~ NA_character_
          ),
          
          # Map exact source column name for altprey3
          altprey3_data_col = case_when(
            altprey3 == "capelin"   ~ "X13_mid_il_capelin",
            altprey3 == "krill"     ~ "X12_egoa.krill",
            TRUE ~ NA_character_
          )
        ) %>%
        arrange(altprey1,altprey2,altprey3)
      
      akpred_dfa_prey_mapped
      
      #add pred_name_sem_data-----
      names(ak_data)
      names(nccpred_dfa_matched_sem_names)
      names(akpred_dfa_prey_mapped)
      akpred_dfa_prey_mapped %>% arrange(short_name_lower) %>% pull(short_name_lower) 
      
      library(dplyr)
      
      akpred_dfa_prey_mapped_sem_name <- akpred_dfa_prey_mapped %>%
        mutate(
          pred_name_sem_data = case_when(
            short_name_lower == "arrowtoothflounderbiomass_predak_2026"   ~ "X15_ArrowtoothFlounderBiomass_predAK",
            short_name_lower == "sablefishbiomass_predak_2026"           ~ "X15_sablefishBiomass_predAK",
            short_name_lower == "sablefishrecruitment_predak_2026"       ~ "X15_sablefishRecruitment_predAK",
            short_name_lower == "pacificcodbiomass_predak_2026"          ~ "X15_PacificCodBiomass_predAK",
            short_name_lower == "spinydogfishgoa_predak_2026"           ~ "X15_spinyDogfishGoA_predAK",
            short_name_lower == "halibutbiomassage8plus_2yrlead_predak_2026" ~ "X15_halibutBiomassAge8plus_2yrLead_predAK",
            short_name_lower == "sleepersharkbsai_predak_2026"          ~ "X15_DFA_sleeperSharkBSAI_predAK",
            short_name_lower == "spinydogfishbsai_predak_2026"          ~ "X15_spinyDogfishBSAI_predAK",
            short_name_lower == "salmonsharkbsai_predak_2026"           ~ "X15_salmonSharkBSAI_predAK",
            short_name_lower == "salmonsharkgoa_predak_2026"            ~ "X15_salmonSharkGoA_predAK",
            short_name_lower == "sharkcatchgoa_predak_2026"             ~ "X15_sharkCatchGoA_predAK",
            
            # Retain original values for any items not listed above
            TRUE ~ short_name_lower
          )
        )
      
      # akpred_dfa_prey_mapped_sem_name <- akpred_dfa_prey_mapped %>%
      # rename(X15_ArrowtoothFlounderBiomass_predAK = arrowtoothflounderbiomass_predak_2026,
      #        X15_sablefishBiomass_predAK = sablefishbiomass_predak_2026,
      #        X15_sablefishRecruitment_predAK = sablefishrecruitment_predak_2026, 
      #        X15_PacificCodBiomass_predAK = pacificcodbiomass_predak_2026,
      #        X15_spinyDogfishGoA_predAK = spinydogfishgoa_predak_2026,
      #        X15_halibutBiomassAge8plus_2yrLead_predAK = halibutbiomassage8plus_2yrlead_predak_2026,
      #        X15_DFA_sleeperSharkBSAI_predAK = sleepersharkbsai_predak_2026,
      #        X15_spinyDogfishBSAI_predAK = spinydogfishbsai_predak_2026, 
      #        X15_salmonSharkBSAI_predAK = salmonsharkbsai_predak_2026,
      #        X15_salmonSharkGoA_predAK = salmonsharkgoa_predak_2026,
      #        X15_sharkCatchGoA_predAK = salmonsharkgoa_predak_2026
      #        )
      
#Save akpred w/ prey column names-------
      write.csv(akpred_dfa_prey_mapped, file.path(outputDir, "akpred_dfa_prey_mapped_sem_data_names.csv"), row.names = FALSE)

      
#5. Bind NCC to AK prey------      
      nccpred_dfa_prey_mapped <-  read.csv( file.path(outputDir, "nccpred_dfa_prey_mapped.csv"), row.names = NULL) %>%
        mutate(altprey3=NA_character_,altprey3_data_col=NA_character_) %>%
        relocate(altprey3, .after = altprey2)
      names(nccpred_dfa_prey_mapped)
      names(akpred_dfa_prey_mapped)
      
      all_pred_dfa_altprey<-rbind(nccpred_dfa_prey_mapped,akpred_dfa_prey_mapped) %>%
          rename(pred_dfa_cols=dfa_cols)
      all_pred_dfa_altprey %>% select(-constituent_predators)
      
      all_pred_dfa_altprey$pred_dfa_cols

      
#Save all_pred_dfa_altprey lookup table-------
      names(nccpred_dfa_matched_sem_names)
      all_pred_dfa_altprey_lookup_sem_names<- all_pred_dfa_altprey %>%
        rename(altprey1_name_sem_data = altprey1_data_col,
               altprey2_name_sem_data = altprey2_data_col,
               altprey3_name_sem_data = altprey3_data_col)
      
      write.csv(all_pred_dfa_altprey, file.path(outputDir, "all_pred_dfa_altprey.csv"), row.names = FALSE)
      
      names(all_pred_dfa_altprey)
            
      #change all the column names back to Lisanames????------
var_lookup<-read.csv("C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/data_Lisa/var_lookup_NCC_AK.SAR.csv",row.names = 1)
var_lookup

names(all_pred_dfa_altprey)
head(all_pred_dfa_altprey)

# all_pred_dfa_altprey<-all_pred_dfa_altprey %>%
#   select(pred_dfa_cols,)
#   left_join(var_lookup,join_by()


#Merge prey data?
      #Save prey data------
altprey_data_ak<-read.csv( file.path(outputDir, "altprey_data_ak_1996_2025.csv"), row.names = NULL);names(altprey_data_ak)

names(altprey_data_ak)

altprey_data_ak_lisanames <- names(altprey_data_ak)
altprey_data_ncc<-read.csv( file.path(outputDir, "altprey_master_source_names_1996_2025.csv"), row.names = NULL);names(altprey_data_ncc)


names(altprey_data_ak)[2:4]<-c("")
head(var_lookup)

var_lookup %>% filter(guildname==names(altprey_data_ak)[2]) %>% pull(Lisaname)
var_lookup %>% filter(guildname==names(altprey_data_ak)[3]) %>% pull(Lisaname)
var_lookup %>% filter(guildname==names(altprey_data_ak)[2]) %>% pull(Lisaname)


var_lookup %>% filter(guildname==all_pred_dfa_altprey$pred_dfa_cols[2]) %>% pull(Lisaname)
