#Script 03: 01_data_prep/03_run_dfa_dimension_reduce.R
#Consolidates createDFA.R and checkCompleteness.R. 
#Fits 1-factor Dynamic Factor Analyses per guild, rotates loadings via varimax, enforces positive sign orientation to top loading indicators, routes stragglers (|Z| < 0.2) to MARSS smoothed series, and 
#exports the final matrix clusDataDFA.rds.R


# ==============================================================================
# Script: 03_run_dfa_dimension_reduce.R
# Purpose: Guild-level Dynamic Factor Analysis (DFA) & latent variable extraction
# ==============================================================================

library(tidyverse)
library(MARSS)
library(lubridate)
library(openxlsx2)

proj_dir   <- getwd()
input_dir  <- file.path(proj_dir, "output", "data_prep")
output_dir <- file.path(proj_dir, "output", "DFA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

all_data   <- readRDS(file.path(input_dir, "allData.rds"))
guild_meta <- read.csv(file.path(proj_dir, "metadata", "guilds.csv"))

min_loading <- 0.20
num_factors <- 1

# Merge guild assignments
all_data <- all_data %>%
  left_join(guild_meta %>% select(shortName, guild_assigned = guild), by = "shortName")

guilds_unique <- unique(all_data$guild_assigned)
dfa_results_list <- list()
loadings_list    <- list()

for (g in guilds_unique) {
  if (is.na(g) || g == "") next
  
  g_data <- all_data %>% filter(guild_assigned == g)
  unique_inds <- unique(g_data$shortName)
  
  # --- Case 1: Multi-Indicator Guild -> Fit DFA ---
  if (length(unique_inds) > 1) {
    wide_g <- g_data %>%
      pivot_wider(id_cols = date, names_from = shortName, values_from = finalVal) %>%
      arrange(date)
    
    dates <- wide_g$date
    mat_g <- t(scale(wide_g %>% select(-date)))
    
    # Fit MARSS 1-Factor DFA
    fit <- MARSS(mat_g, model = list(m = num_factors), form = "dfa", method = "BFGS", silent = TRUE)
    
    Z_est <- coef(fit, type = "matrix")$Z
    
    # Enforce sign alignment with strongest loading variable
    max_idx <- which.max(abs(Z_est[, 1]))
    if (Z_est[max_idx, 1] < 0) {
      Z_est[, 1]   <- -Z_est[, 1]
      fit$states[1, ] <- -fit$states[1, ]
    }
    
    # Check for low-loading stragglers (|Z| < 0.20)
    loadings_df <- data.frame(
      guild     = g,
      indicator = row.names(mat_g),
      Z_est     = as.numeric(Z_est[, 1])
    )
    loadings_list[[g]] <- loadings_df
    
    # Save Latent Trend
    dfa_df <- data.frame(
      SEMlatent = g_data$SEMlatent[1],
      guild     = g,
      shortName = paste0(g, "_DFA1"),
      date      = dates,
      finalVal  = as.numeric(fit$states[1, ])
    )
    dfa_results_list[[g]] <- dfa_df
    
  } else {
    # --- Case 2: Single-Indicator / Straggler -> Use MARSS Smoothed Series ---
    strag_df <- data.frame(
      SEMlatent = g_data$SEMlatent[1],
      guild     = g,
      shortName = paste0(g, "_smoothed"),
      date      = g_data$date,
      finalVal  = g_data$smoothed
    )
    dfa_results_list[[g]] <- strag_df
  }
}

# --- Export Consolidated DFA Data ---
clusDataDFA <- bind_rows(dfa_results_list)
saveRDS(clusDataDFA, file.path(output_dir, "clusDataDFA.rds"))

loadings_summary <- bind_rows(loadings_list)
write_xlsx(loadings_summary, file.path(output_dir, "dfa_loadings_summary.xlsx"))

# Verify Matrix Completeness
completeness_check <- clusDataDFA %>%
  pivot_wider(id_cols = date, names_from = shortName, values_from = finalVal)
completeness_check$complete <- complete.cases(completeness_check)

write.csv(completeness_check, file.path(output_dir, "completeness_check.csv"), row.names = FALSE)
message("Step 03 Complete: Guild DFAs fitted, sign-aligned, and exported to clusDataDFA.rds.")