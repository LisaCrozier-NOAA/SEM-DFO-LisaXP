
 #Step 1. Find which indicators went into the important DFAs
 
  
  library(dplyr)

# 1. Target the specific Lisaname we want to investigate
target_lisa_name <- "X09_DFA_HakeAge5Plus"
target_lisa_name <- "X09_DFA_HakeAge5Plus"

# 2. Trace backwards to get the underlying 'guildname' from your lookup table
# This finds rows where Lisaname matches, and cleans the 'X' prefix off the guildname
target_guild_info <- var_lookup_NCC_AK %>%
  filter(Lisaname == target_lisa_name)

if(nrow(target_guild_info) == 0) {
  stop("Could not find the target name in var_lookup_NCC_AK!")
}

# Pull the guildname (e.g., "X09.PredFishNCC_b_DFA1")
full_guildname <- target_guild_info$guildname[1]

# Strip off the leading "X" and the trailing "_DFA1" to match your guild.data format
# This converts "X09.PredFishNCC_b_DFA1" -> "09.PredFishNCC_b"
core_dfa_name <- gsub("^X|\\.DFA1$|_DFA1$", "", full_guildname)

message(paste("Successfully traced backwards! Core DFA identifier is:", core_dfa_name))


# 3. Query guild.data to pull all original raw time series names collapsed into this DFA
raw_components <- guild.data %>%
  filter(DFAguild == core_dfa_name) %>%
  select(shortName) %>%
  distinct() %>%
  pull(shortName)

# 4. Print out your list of raw time series indicators!
print("The raw time series combined into this DFA are:")
print(raw_components)
