library(dplyr)
library(ggplot2)
library(Cairo)
library(lubridate)

#write.csv(importance_topvar,"LisaXP/outputs_4/importance_topvar_daic3.csv")

# ==============================================================================
# PREPARATION & DISCOVERY
# ==============================================================================

# 1. Filter the top-performing variables down to just the DFA-derived indicators
dfa_priority_list <- importance_topvar %>%
  filter(grepl("DFA", indNames)) %>%
  distinct(indNames) %>%
  pull(indNames)

message(paste("Found", length(dfa_priority_list), "important DFA indicators within dAIC <= 3."))

# Open a PDF file to store all deconstruction plots
CairoPDF("LisaXP/outputs_4/DFA_Top_Variables_Deconstruction.pdf", width = 12, height = 10)

# Initialize a container to store our lavaan measurement model syntax
latent_syntax_blocks <- c("# ==========================================\n# AUTOMATED MEASUREMENT MODEL SYNTAX\n# ==========================================\n")

# ==============================================================================
# MAIN TRACKING & RENDERING LOOP
# ==============================================================================
#dfa_ind <- dfa_priority_list[1]

for (dfa_ind in dfa_priority_list) {
  
  # Standardize naming format to ensure we can match against var_lookup_NCC_AK
  # e.g., if it comes in as "09.PredFishNCC_b_DFA1", ensure we can match it
  lookup_row <- var_lookup_NCC_AK %>%
    filter(guildname == dfa_ind | paste0("X", guildname) == dfa_ind | guildname == paste0("X", dfa_ind))
  
  # If it's a smoothed variant or slightly different prefix, use loose regex matching
  if (nrow(lookup_row) == 0) {
    clean_search <- gsub("^X|_DFA1$", "", dfa_ind)
    lookup_row <- var_lookup_NCC_AK %>% 
      filter(grepl(clean_search, guildname))
  }
  
  # Skip safely if no reference is found in your key
  if (nrow(lookup_row) == 0) {
    warning(paste("Could not map indicator:", dfa_ind, "in var_lookup_NCC_AK. Skipping."))
    next
  }
  
  target_lisa_name <- lookup_row$Lisaname[1]
  full_guildname   <- lookup_row$guildname[1]
  core_dfa_name    <- gsub("^X|\\.DFA1$|_DFA1$", "", full_guildname)
  
  # --- Step A: Extract Raw Components from guild.data ---
  raw_components_df <- guild.data %>%
    filter(DFAguild == core_dfa_name)
  
  if (nrow(raw_components_df) == 0) {
    warning(paste("No raw component rows found in guild.data for core name:", core_dfa_name))
    next
  }
  
  raw_ts_long <- raw_components_df %>%
    mutate(year = year(date)) %>%
    filter(year >= 1998, year <= 2021) %>%
    select(year, shortName, finalVal) %>%
    rename(Variable = shortName, Value = finalVal) %>%
    mutate(Source = "Raw Input")
  
  # Get a vector of the unique raw component names for the SEM syntax step later
  raw_short_names <- unique(raw_ts_long$Variable)
  
  # --- Step B: Extract Final Output from guild.dfas1 ---
  if (!target_lisa_name %in% colnames(guild.dfas1)) {
    warning(paste("Lisaname column", target_lisa_name, "not found in guild.dfas1. Skipping plot."))
    next
  }
  
  final_dfa_long <- guild.dfas1 %>%
    select(year, all_of(target_lisa_name)) %>%
    rename(Value = !!sym(target_lisa_name)) %>%
    mutate(
      Variable = target_lisa_name,
      Source = "Final DFA Output"
    )
  
  # --- Step C: Combine and Render Plot to PDF ---
  combined_plot_data <- bind_rows(raw_ts_long, final_dfa_long) %>%
    mutate(Plot_Label = case_when(
      Source == "Final DFA Output" ~ paste("★ FINAL DFA OUTPUT:", target_lisa_name),
      Variable %in% names(my_pretty_names) ~ paste(my_pretty_names[Variable]),
      TRUE ~ paste(Variable)
    ))
  
  # Adjust dynamic panel height multiplier based on how many sub-components exist
  num_panels <- length(unique(combined_plot_data$Plot_Label))
  
  deconstruction_plot <- ggplot(combined_plot_data, aes(x = year, y = Value)) +
    geom_line(aes(color = Source), size = 1.2, show.legend = FALSE) +
    geom_point(aes(color = Source), size = 1.8, show.legend = FALSE) +
    scale_color_manual(values = c("Final DFA Output" = "firebrick", "Raw Input" = "darkgray")) +
    facet_wrap(~Plot_Label, scales = "free_y", ncol = 1) +  
    theme_minimal(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    ) +
    labs(
      title = paste("DFA Deconstruction:", core_dfa_name),
      x = "Year",
      y = "Value"
    )
  
  print(deconstruction_plot)
  
  # --- Step D: Construct the lavaan Measurement Syntax ---
  # Cleans up core name to create a valid, legal latent variable name (e.g., Latent_PredFishNCC_b)
  latent_var_name <- paste0("Latent_", gsub("[^a-zA-Z0-9_]", "_", core_dfa_name))
  
  # Group indicators using standard lavaan format: Latent =~ indicator1 + indicator2 + indicator3
  syntax_line <- paste0(latent_var_name, " =~ ", paste(raw_short_names, collapse = " + "))
  latent_syntax_blocks <- c(latent_syntax_blocks, syntax_line)
}

dev.off()
message("PDF plotting complete.")


library(dplyr)
library(tidyr)
library(lubridate)

#Create dataset ready to enter the SEM-----------
        #with all of the raw time series (how is this going to work with missing data?)
        # ==============================================================================
        # 1. ISOLATE ALL TARGET RAW COMPONENT NAMES
        # ==============================================================================
        # Look up all unique core DFA guilds that made the dAIC <= 3 cut
        top_dfa_guilds <- importance_topvar %>%
          filter(grepl("DFA", indNames)) %>%
          distinct(indNames) %>%
          # Map names cleanly back to the core DFA identifiers
          mutate(clean_ind = gsub("^X|_DFA1$", "", indNames)) %>%
          left_join(
            var_lookup_NCC_AK %>% 
              mutate(clean_guild = gsub("^X|\\.DFA1$|_DFA1$", "", guildname)),
            by = c("clean_ind" = "clean_guild")
          ) %>%
          mutate(core_dfa_name = coalesce(clean_ind, Lisaname)) %>%
          distinct(core_dfa_name) %>%
          pull(core_dfa_name)
        
        # ==============================================================================
        # 2. EXTRACT AND WIDEN THE RAW DATA (guild.data)
        # ==============================================================================
        raw_components_wide <- guild.data %>%
          # Filter to include only the raw series matching our high-priority DFAs
          filter(DFAguild %in% top_dfa_guilds) %>%
          mutate(year = year(date)) %>%
          filter(year >= 1998, year <= 2021) %>%
          # Deduplicate to ensure no structural padding errors during pivot
          distinct(year, shortName, .keep_all = TRUE) %>%
          select(year, shortName, finalVal) %>%
          # Pivot the raw rows into unique standalone columns
          pivot_wider(names_from = shortName, values_from = finalVal)%>%
          mutate(across(-year, ~ as.numeric(scale(.))))
        
        # ==============================================================================
        # 3. MERGE RAW CHANNELS WITH THE BASELINE DATA (guild.dfas1)
        # ==============================================================================
        # Combine your core non-DFA structural variables with the new raw wide variables
        sem_master_data <- guild.dfas1 %>%
          select(year, everything()) %>%
          left_join(raw_components_wide, by = "year")
        
        # ==============================================================================
        # 4. SANITY CHECKS
        # ==============================================================================
        message("--- Master SEM Dataset Created Successfully! ---")
        message(paste("Total Rows (Years):", nrow(sem_master_data)))
        message(paste("Total Columns Available:", ncol(sem_master_data)))
        
        # Print a quick snapshot of the columns added
        new_raw_cols <- setdiff(colnames(raw_components_wide), "year")
        print("New individual raw indicators added as available columns:")
        print(new_raw_cols)
        
        write.csv(sem_master_data,"LisaXP/outputs_4/sem_master_data.csv",row.names = FALSE)




# ==============================================================================
# PRINT MEASUREMENT MODEL CODES FOR SEM
# ==============================================================================
cat(paste(latent_syntax_blocks, collapse = "\n\n"))



# ==========================================
# AUTOMATED MEASUREMENT MODEL SYNTAX
# ==========================================


# Latent_06_Cond1NCC =~ IGF_mu_2025 + StomFull_Jun_2025
# 
# Latent_05_ForageFishNCC =~ sardine_NCC + abundHerring + abundSardine + sardine_GAM_2025
# 
# Latent_09_PredFishNCC_b =~ iphc20bigSkate + HakeAge5Plus_2025 + darkblotchedRockfish_2025 + greestripedRockfish_2025 + hakeFisheryWA_2025 + jackmackerel_GAM_2025 + pacificmackerel_GAM_2025 + PacificSpinyDogfish_2025
# 
# Latent_12_ZooPreyAK =~ crestedAuk_WAI + leastAuk_WAI + zoop_EGoA + biomassAmphShelfSum + biomassEuphShelfSum + biomassMysidShelfSum
# 
# Latent_13_FishPreyAK_b =~ gadid_EAI + ammod_EAI + capelin_WGoA + estAbundHerringRecruits + WGOA_DFA_seabirds_2026
# 
# Latent_13_FishPreyAK =~ rhinoAuk_EGoA + sitkaHerring_EGoA + WGOA_DFA_lowerTrophic_2026 + WGOA_DFA_midTrophic_2026 + pollockBiomassAIage1plus_predAK_2026
# 
# Latent_15_PredFishAK_b =~ sleeperSharkBSAI_predAK_2026 + sleeperSharkGoA_predAK_2026