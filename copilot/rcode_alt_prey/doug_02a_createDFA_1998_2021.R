# ==============================================================================
# Script: doug_02_createDFA_19982021.R
# Purpose: Fit DFAs (loading > 0.2) & MARSS Kalman Smoothing across Guilds (1998–2021)
# Output: copilot/outputs_10/
# ==============================================================================

library(tidyverse)
library(MARSS)
library(lubridate)
library(openxlsx2)

# ---------------------------------------------------------------------------
# Constants & Paths
# ---------------------------------------------------------------------------
rootdir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"

path <- c(
  "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results", 
  "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/"
)[1]

workingDir <- file.path(path, "analyzeAKindices")

indicatorsFile <- file.path(workingDir, "indicators.csv")
dataDir        <- file.path(workingDir, "data")
#outputDir      <- file.path(rootdir, "copilot/outputs_10")
outputDir      <- file.path(rootdir, "copilot/outputs_altprey")
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

guildFile  <- file.path(workingDir, "guildsWithExclude.csv")
guildCol   <- "latestGuild"

numFactors <- 1
minLoading <- 0.2

# ---------------------------------------------------------------------------
# Load Input Data
# ---------------------------------------------------------------------------
datWide <- read.csv(file.path(outputDir, "datWide_1998_2021_qualified.csv"))

clusData <- datWide %>% 
  pivot_longer(cols = -Year, names_to = "shortName", values_to = "finalVal") %>% 
  rename(year = Year) %>% 
  mutate(date = ymd(paste0(year, "-01-01"))) %>% 
  filter(!is.na(finalVal))

guilds_master <- read.csv(guildFile)

guilds <- guilds_master %>% 
  mutate(guild = .data[[guildCol]]) %>% 
  filter(shortName %in% clusData$shortName, guild != "") %>% 
  rename(guildSEMnode = SEMnode) %>% 
  select(shortName, guild, guildSEMnode)

clusData <- left_join(clusData, guilds, by = "shortName") %>% 
  filter(!is.na(guild))

unique_guild_combos <- guilds %>% 
  distinct(guild, guildSEMnode) %>% 
  filter(guild != "", guildSEMnode != "") %>% 
  arrange(guildSEMnode, guild)

# ---------------------------------------------------------------------------
# Run DFAs & MARSS Smoothing across Guilds
# ---------------------------------------------------------------------------
clusDataDFAlist <- list()
loadingsList    <- list()
smoothedIndList <- list()
marssFitList    <- list() # Preserves full MARSS model objects for extension

for (i in seq_len(nrow(unique_guild_combos))) {
  thisRow          <- unique_guild_combos[i, ]
  thisGuild        <- thisRow$guild
  thisGuildSEMnode <- thisRow$guildSEMnode
  
  thisClusData <- clusData %>% 
    filter(guild == thisGuild, guildSEMnode == thisGuildSEMnode)
  
  unique_inds <- unique(thisClusData$shortName)
  
  # CASE 1: Multi-Indicator Guild -> Run Dynamic Factor Analysis (DFA)
  if (length(unique_inds) > 1) {
    cat("Fitting DFA for Guild:", thisGuild, "(", length(unique_inds), "indicators)\n")
    
    wide_mat <- thisClusData %>% 
      select(shortName, date, finalVal) %>% 
      pivot_wider(id_cols = date, names_from = shortName, values_from = finalVal) %>% 
      arrange(date)
    
    thisDate <- wide_mat$date
    
    mat_data <- wide_mat %>% 
      select(-date) %>% 
      mutate(across(everything(), ~ as.vector(scale(.)))) %>% 
      t()
    
    fit <- MARSS(mat_data, model = list(m = numFactors), form = "dfa", method = "BFGS", silent = TRUE)
    
    # Store fitted object for future state projection
    marssFitList[[paste0(thisGuild, "_DFA1")]] <- fit
    
    # Enforce sign consistency: match sign of highest absolute loading
    Z_matrix <- fit$par$Z
    maxIndex <- which.max(abs(Z_matrix[, 1]))
    if (Z_matrix[maxIndex, 1] < 0) {
      fit$par$Z[, 1]   <- -fit$par$Z[, 1]
      fit$states[1, ]  <- -fit$states[1, ]
    }
    
    processes <- fit$states
    Z_est     <- fit$par$Z
    
    loadingsDF <- data.frame(
      guild        = thisGuild,
      guildSEMnode = thisGuildSEMnode,
      indicator    = rownames(mat_data),
      Z_est        = as.numeric(Z_est)
    )
    
    loadingsDF$newGuild <- loadingsDF$guild
    strag_idx <- 0
    for (r in seq_len(nrow(loadingsDF))) {
      if (abs(loadingsDF[r, "Z_est"]) < minLoading) {
        loadingsDF[r, "newGuild"] <- paste0(loadingsDF[r, "guild"], "_", strag_idx)
        strag_idx <- strag_idx + 1
      }
    }
    loadingsList[[length(loadingsList) + 1]] <- loadingsDF
    
    thisDFA <- data.frame(
      SEMlatent    = thisGuildSEMnode,
      guild        = thisGuild,
      shortName    = paste0(thisGuild, "_DFA1"),
      date         = thisDate,
      finalVal     = processes[1, ]
    )
    clusDataDFAlist[[length(clusDataDFAlist) + 1]] <- thisDFA
    
  } else {
    # CASE 2: Single-Indicator Guild / Straggler -> Apply MARSS Kalman Smoothing
    cat("Smoothing Straggler:", thisGuild, "(1 indicator:", unique_inds, ")\n")
    
    thisDate <- thisClusData$date
    vec_data <- scale(thisClusData$finalVal) %>% as.vector()
    
    fit_single <- MARSS(vec_data, fit = FALSE)
    fit_single$par <- fit_single$start
    kfList <- MARSSkf(fit_single)
    smoothed_vals <- as.numeric(t(kfList$xtT))
    
    # Store fitted object
    marssFitList[[paste0(thisGuild, "_smoothed")]] <- fit_single
    
    thisDFA <- data.frame(
      SEMlatent    = thisGuildSEMnode,
      guild        = thisGuild,
      shortName    = paste0(thisGuild, "_smoothed"),
      date         = thisDate,
      finalVal     = smoothed_vals
    )
    clusDataDFAlist[[length(clusDataDFAlist) + 1]] <- thisDFA
    
    smoothedIndList[[length(smoothedIndList) + 1]] <- data.frame(
      guild            = thisGuild,
      DFAname          = paste0(thisGuild, "_smoothed"),
      rankedIndicators = unique_inds
    )
  }
}

# ---------------------------------------------------------------------------
# Combine and Export Results
# ---------------------------------------------------------------------------
clusDataDFA <- bind_rows(clusDataDFAlist)
saveRDS(clusDataDFA, file.path(outputDir, "clusDataDFA.rds"))

# Save fitted MARSS models for out-of-sample extension
saveRDS(marssFitList, file.path(outputDir, "marss_fits.rds"))

# Export Wide Version of DFA Indicators
clusDataDFA_wide <- clusDataDFA %>% 
  mutate(Year = year(date)) %>% 
  pivot_wider(id_cols = Year, names_from = shortName, values_from = finalVal) %>% 
  arrange(Year)

write.csv(clusDataDFA_wide, file.path(outputDir, "clusDataDFA_wide_1998_2021.csv"), row.names = FALSE)

# Export Loadings and Guild Assignments
loadings <- bind_rows(loadingsList)
write_xlsx(loadings, file.path(outputDir, "loadings.xlsx"))

newGuilds <- loadings %>% select(indicator, newGuild) %>% rename(shortName = indicator)
guilds_updated <- left_join(guilds_master, newGuilds, by = "shortName")
write.csv(guilds_updated, file.path(outputDir, "newGuilds.csv"), row.names = FALSE)

# Export Ranked Indicators Summary
DFAnames <- loadings %>% 
  mutate(DFAname = paste0(guild, "_DFA1")) %>% 
  arrange(guild, desc(abs(Z_est))) %>% 
  group_by(guild, DFAname) %>% 
  summarize(rankedIndicators = paste(indicator, collapse = " - "), .groups = "drop")

if (length(smoothedIndList) > 0) {
  DFAnames <- bind_rows(DFAnames, bind_rows(smoothedIndList))
}

write.csv(DFAnames %>% arrange(DFAname), file.path(outputDir, "rankedIndicators.csv"), row.names = FALSE)

cat("\nDFA Pipeline Complete! All outputs and marss_fits.rds saved to:", outputDir, "\n")
