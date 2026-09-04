

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

clusDataDFA_wide <- read.csv(file.path(outputDir, "clusDataDFA_wide_1996_2025.csv"))

    
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

akpred 

#add ssl pup info-----
      ssl<-read.csv("copilot/outputs_2/ssl.dat.csv",row.names = NULL)  %>%
        select(Year=year,X10_ssl_seak_pup_pred=ssl_seak_pup_pred) 
      
      head(ssl)
      
      akpred<-rbind(akpred,c("ssl_seak_pup_pred","ssl_seak_pup_pred","predak",NA_character_,NA_character_,"herring","capelin"))
      akpred
      
      # Display complete output table
      print(as_tibble(akpred), n = Inf)
      
      #unique prey needed-------
      unique(c(akpred$altprey1,akpred$altprey2,akpred$altprey3))
     # "krill"       "pollock"   "herring"     "pink_salmon" "capelin"   
     #  X12_egoa.krill (krill), X13.FishPreyAK_DFA1 (pollock), X13_stka_herr_matbiom (herring), X14.CompAK_1_smoothed (pink_salmon),X13_mid_il_capelin (capelin)
      
      
# Save updated Alaska prey assignments manifest--------
      write.csv(akpred, file.path(outputDir, "akpred_altprey_assigned.csv"), row.names = FALSE)

# 2. Get prey data----------
      #look at options
      rankedIndicators <- read.csv(file.path(outputDir, "rankedIndicators.csv"), row.names = NULL) %>%
        mutate(dfa_cols=paste0("X",DFAname))%>% 
        relocate(dfa_cols,.after=DFAname) %>%
        filter(grepl("12|13|14|15",DFAname)) %>%
        select(dfa_cols,rankedIndicators)
      
      rankedIndicators
      
      #haven't included EGOA into Doug's process yet
      egoa.krill<-read.csv("data_Lisa/ak_yr.csv",row.names = NULL) %>%
        rename(Year=year) %>% 
        filter(Year>=1996) %>% 
        select(Year,X12_egoa.krill=secm_euph_dens,X13_mid_il_capelin=mid_il_capelin, X13_stka_herr_matbiom=stka_herr_matbiom) ;names(egoa.krill)
      
      # Construct wide prey dataset with exact raw source column names and 2-year leads
      altprey_data_ak <- clusDataDFA_wide %>%
        select(
          Year,
          X13.FishPreyAK_DFA1,        # WGOA_DFA_midTrophic_2026 - pollockBiomassAIage1plus_predAK_2026 - WGOA_DFA_lowerTrophic_2026 - rhinoAuk_EGoA - sitkaHerring_EGoA
          X13.FishPreyAK_b_DFA1,      # WGOA_DFA_seabirds_2026 - capelin_WGoA - ammod_EAI - estAbundHerringRecruits - gadid_EAI
          X14.CompAK_1_smoothed,      # pinkSalmonNorthAmerica_2025
        ) %>%
        left_join(egoa.krill, by = "Year")

      print(as_tibble(altprey_data_ak),n=Inf)

      
#Save prey data------
      write.csv(altprey_data_ak, file.path(outputDir, "altprey_data_ak_1996_2025.csv"), row.names = FALSE)


# 3. Reduce pred list to unique dfas --------
# Add dfa names for predators, check for conflicts in alt prey
      akpred<-read.csv( file.path(outputDir, "akpred_altprey_assigned.csv"), row.names = NULL)
      
      # 3. Match individual indicators w/ dfa names via grepl (case-insensitive substring search inside the collapsed string)
      akpred_dfa_matched <- akpred %>%
        # Join all combinations to test each short_name against each collapsed DFA string
        left_join(rankedIndicators, by = character(), relationship = "many-to-many") %>%
        # Keep rows where short_name_lower is found inside the raw collapsed string
        # (or keep as NA if an indicator isn't part of any DFA in rankedIndicators)
        filter(
          is.na(rankedIndicators) | 
            mapply(function(pattern, text) grepl(pattern, text, ignore.case = TRUE), 
                   short_name_lower, rankedIndicators)
        ) %>%
        relocate(dfa_cols,.before=latest_guild)
      
      akpred_dfa_matched
      
      prey_conflicts <- akpred_dfa_matched %>%
        group_by(dfa_cols) %>%
        summarise(
          num_predators       = n(),
          predators_in_dfa    = paste(unique(short_name_lower), collapse = ", "),
          altprey1_choices    = paste(unique(altprey1[!is.na(altprey1)]), collapse = " | "),
          altprey2_choices    = paste(unique(altprey2[!is.na(altprey2)]), collapse = " | "),
          altprey3_choices    = paste(unique(altprey3[!is.na(altprey2)]), collapse = " | "),
          altprey1_conflict   = n_distinct(altprey1[!is.na(altprey1)]) > 1,
          altprey2_conflict   = n_distinct(altprey2[!is.na(altprey2)]) > 1,
          altprey3_conflict   = n_distinct(altprey3[!is.na(altprey2)]) > 1,
          .groups = "drop"
        ) %>%
        filter(altprey1_conflict | altprey2_conflict | altprey3_conflict)
      
      prey_conflicts
      prey_conflicts %>% select(dfa_cols,num_predators,altprey1_choices,altprey2_choices,altprey3_choices,predators_in_dfa)
      
      
      # Collapse to Unique DFA Lines---------
      akpred_dfa_unique <- akpred_dfa_matched %>%
        filter(!is.na(dfa_cols)) %>%
        mutate(is_2yrlead="FALSE") %>%
        group_by(dfa_cols) %>%
        summarise(
          # Aggregate all constituent short names for reference
          constituent_predators = paste(unique(short_name_lower), collapse = " | "),
          altprey1              = unique(altprey1[!is.na(altprey1)])[1],
          altprey2              = unique(altprey2[!is.na(altprey2)])[1],
          altprey3              = unique(altprey3[!is.na(altprey3)])[1],
          is_2yrlead              = unique(is_2yrlead[!is.na(is_2yrlead)])[1],
          .groups = "drop"
        )
      
      # View the final unconflicted DFA table
      print(as_tibble(akpred_dfa_unique), n = Inf)
      print(as_tibble(nccpred_dfa_unique), n = Inf)
      
      #save dfa prey assignments----------
      write.csv(akpred_dfa_unique, file.path(outputDir, "akpred_dfa_altprey_assigned.csv"), row.names = FALSE)

      
      
 # ---------------------------------------------------------------------------
 #4. Map data source column names into akpred_dfa_unique-------------------------
 # ---------------------------------------------------------------------------
      
      akpred_dfa_unique <- read.csv(file.path(outputDir, "akpred_dfa_altprey_assigned.csv"), row.names = NULL)
      
      #  X12_egoa.krill (krill), X13.FishPreyAK_DFA1 (pollock), X13_stka_herr_matbiom (herring), X14.CompAK_1_smoothed (pink_salmon),X13_mid_il_capelin (capelin)
      
      akpred_dfa_prey_mapped <- akpred_dfa_unique %>%
        mutate(
          # Map exact source column name for altprey1
          altprey1_data_col = case_when(
            altprey1 == "pollock"  ~ "X13.FishPreyAK_DFA1",
            altprey1 == "krill"   ~ "X12_egoa.krill",
            TRUE ~ NA_character_
          ),
          
          # Map exact source column name for altprey2
          altprey2_data_col = case_when(
            altprey2 == "herring"   ~ "X13_stka_herr_matbiom",
            altprey2 == "pink_salmon"  ~ "X14.CompAK_1_smoothed",
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
#Save akpred w/ prey column names-------
      write.csv(akpred_dfa_prey_mapped, file.path(outputDir, "akpred_dfa_prey_mapped.csv"), row.names = FALSE)

      
#5. Bind NCC to AK prey------      
      nccpred_dfa_prey_mapped <-  read.csv( file.path(outputDir, "nccpred_dfa_prey_mapped.csv"), row.names = NULL) %>%
        mutate(altprey3=NA_character_,altprey3_data_col=NA_character_) %>%
        relocate(altprey3, .after = altprey2)
      names(nccpred_dfa_prey_mapped)
      names(akpred_dfa_prey_mapped)
      
      all_pred_dfa_altprey<-rbind(nccpred_dfa_prey_mapped,akpred_dfa_prey_mapped) %>%
          rename(pred_dfa_cols=dfa_cols)
      all_pred_dfa_altprey %>% select(-constituent_predators)

#Save all_pred_dfa_altprey-------
      write.csv(all_pred_dfa_altprey, file.path(outputDir, "all_pred_dfa_altprey.csv"), row.names = FALSE)
      
            