

# ==============================================================================
# Script: 02_fit_locked_dfas_and_project_shifts_1995_2025.R
# Purpose: Lock Doug's 1998-2021 DFA parameters, project across extended years 
#          without re-fitting, and apply explicit _smoltyr / _adultyr 2-year shifts.
# Output: metadata/datWide_1998_2021_qualified.csv
#         metadata/datWide_1996_2025_extended_shifted.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(MARSS)
library(openxlsx2)

# ------------------------------------------------------------------------------
# 1. SETUP PATHS
# ------------------------------------------------------------------------------
rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
path      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"
outputDir <- file.path(rootdir, "metadata")

raw_file  <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw.csv")
guildFile <- file.path(path, "guildsWithExclude.csv")
guildCol  <- "latestGuild"

if (!file.exists(raw_file)) stop("Missing raw file. Run 01 script first!")

# Load raw 1995-2025 data
raw_data <- read.csv(raw_file, stringsAsFactors = FALSE)
guilds_master <- read.csv(guildFile, stringsAsFactors = FALSE)

# Define baseline index (1998-2021)
baseline_years <- 1998:2021
full_years     <- raw_data$Year

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS FOR LOCKED SCALING & KALMAN SMOOTHING
# ------------------------------------------------------------------------------
# Scale extended vector using ONLY 1998-2021 mean and SD
scale_baseline_lock <- function(vec, yrs) {
  base_idx  <- which(yrs %in% baseline_years)
  base_mean <- mean(vec[base_idx], na.rm = TRUE)
  base_sd   <- sd(vec[base_idx], na.rm = TRUE)
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

# ------------------------------------------------------------------------------
# 3. PROCESS GUILDS & FIT LOCKED MARSS MODELS
# ------------------------------------------------------------------------------
guilds <- guilds_master %>%
  mutate(guild = .data[[guildCol]]) %>%
  filter(shortName %in% names(raw_data), guild != "") %>%
  rename(guildSEMnode = SEMnode) %>%
  select(shortName, guild, guildSEMnode)

unique_guilds <- guilds %>% distinct(guild, guildSEMnode) %>% arrange(guildSEMnode, guild)

baseline_results <- list()
extended_results <- list()
marss_fits_locked <- list()

for (i in seq_len(nrow(unique_guilds))) {
  g_row  <- unique_guilds[i, ]
  g_name <- g_row$guild
  g_node <- g_row$guildSEMnode
  
  g_inds <- guilds %>% filter(guild == g_name, guildSEMnode == g_node) %>% pull(shortName)
  
  if (length(g_inds) > 1) {
    # --- MULTI-INDICATOR GUILD (DFA) ---
    # 1. Build 1998-2021 Baseline Matrix for Fitting
    base_mat <- raw_data %>%
      filter(Year %in% baseline_years) %>%
      select(all_of(g_inds)) %>%
      mutate(across(everything(), ~ as.vector(scale(.)))) %>%
      t()
    
    # Fit MARSS DFA strictly on 1998-2021 baseline
    fit_base <- MARSS(base_mat, model = list(m = 1), form = "dfa", method = "BFGS", silent = TRUE)
    
    # Enforce sign alignment on highest absolute loading
    Z_est <- fit_base$par$Z
    max_idx <- which.max(abs(Z_est[, 1]))
    sign_flip <- if (Z_est[max_idx, 1] < 0) -1 else 1
    
    fit_base$par$Z[, 1] <- fit_base$par$Z[, 1] * sign_flip
    fit_base$states[1, ] <- fit_base$states[1, ] * sign_flip
    
    # Store locked fit object
    marss_fits_locked[[paste0(g_name, "_DFA1")]] <- fit_base
    
    # Baseline trend (1998-2021)
    base_trend <- fit_base$states[1, ]
    baseline_results[[paste0(g_name, "_DFA1")]] <- tibble(Year = baseline_years, val = base_trend)
    
    # 2. Project across FULL 1995-2025 using LOCKED Z and Covariance
    ext_mat <- raw_data %>%
      select(all_of(g_inds)) %>%
      mutate(across(everything(), ~ scale_baseline_lock(., full_years))) %>%
      t()
    
    marss_proj <- MARSS(ext_mat, model = list(m = 1), form = "dfa", fit = FALSE, silent = TRUE)
    marss_proj$par <- fit_base$par # Inject locked parameters
    
    ext_trend <- as.numeric(MARSSkf(marss_proj)$xtT[1, ])
    extended_results[[paste0(g_name, "_DFA1")]] <- tibble(Year = full_years, val = ext_trend)
    
  } else if (length(g_inds) == 1) {
    # --- SINGLE INDICATOR / STRAGGLER (MARSS Kalman Filter) ---
    ind_name <- g_inds[1]
    
    # 1. Fit baseline univariate MARSS on 1998-2021
    base_vec <- scale(raw_data %>% filter(Year %in% baseline_years) %>% pull(!!sym(ind_name))) %>% as.vector()
    fit_single <- MARSS(base_vec, fit = FALSE, silent = TRUE)
    fit_single$par <- fit_single$start
    
    base_trend <- as.numeric(t(MARSSkf(fit_single)$xtT))
    baseline_results[[paste0(g_name, "_smoothed")]] <- tibble(Year = baseline_years, val = base_trend)
    marss_fits_locked[[paste0(g_name, "_smoothed")]] <- fit_single
    
    # 2. Project full window (1995-2025) with baseline locked scaling
    ext_vec <- scale_baseline_lock(raw_data[[ind_name]], full_years)
    ext_proj <- MARSS(ext_vec, fit = FALSE, silent = TRUE)
    ext_proj$par <- fit_single$par
    
    ext_trend <- as.numeric(t(MARSSkf(ext_proj)$xtT))
    extended_results[[paste0(g_name, "_smoothed")]] <- tibble(Year = full_years, val = ext_trend)
  }
}

# ------------------------------------------------------------------------------
# 4. CONSTRUCT EXTENDED MATRIX AND APPLY EXPLICIT 2-YEAR SHIFTS
# ------------------------------------------------------------------------------
# Combine projected trend long data into wide matrix (1995-2025)
ext_wide <- bind_rows(
  lapply(names(extended_results), function(nm) {
    extended_results[[nm]] %>% mutate(shortName = nm)
  })
) %>%
  pivot_wider(id_cols = Year, names_from = shortName, values_from = val) %>%
  arrange(Year)

# Explicit Shift Application Logic:
# Adult Predators (X10, X15): Original column is _adultyr. Create _smoltyr via lag(2).
# Alternate Prey (X05, X13, X14): Original column is _smoltyr. Create _adultyr via lead(2).

shifted_matrix <- ext_wide %>%
  # --- Adult Mammal Predators (X10) ---
  rename_with(~ paste0(.x, "_adultyr"), contains("10.PredMammal")) %>%
  mutate(across(contains("10.PredMammal") & ends_with("_adultyr"), 
                ~ dplyr::lag(.x, 2), 
                .names = "{gsub('_adultyr', '_smoltyr', .col)}")) %>%
  
  # --- Adult Groundfish Predators (X15) ---
  rename_with(~ paste0(.x, "_adultyr"), contains("15.PredFish")) %>%
  mutate(across(contains("15.PredFish") & ends_with("_adultyr"), 
                ~ dplyr::lag(.x, 2), 
                .names = "{gsub('_adultyr', '_smoltyr', .col)}")) %>%
  
  # --- Alternate Prey (X05, X13, X14) ---
  rename_with(~ paste0(.x, "_smoltyr"), contains("05.Forage") | contains("13.FishPrey") | contains("14.CompAK")) %>%
  mutate(across((contains("05.Forage") | contains("13.FishPrey") | contains("14.CompAK")) & ends_with("_smoltyr"), 
                ~ dplyr::lead(.x, 2), 
                .names = "{gsub('_smoltyr', '_adultyr', .col)}"))

# ------------------------------------------------------------------------------
# 5. EXPORT FINAL MATRICES & LOCKED MARSS CONTAINER
# ------------------------------------------------------------------------------
# Save locked MARSS fits container
saveRDS(marss_fits_locked, file.path(outputDir, "marss_fits_locked_1998_2021.rds"))

# Export 1: Exact Baseline Matrix (1998-2021)
write.csv(ext_wide %>% filter(Year >= 1998 & Year <= 2021), 
          file.path(outputDir, "datWide_1998_2021_qualified.csv"), row.names = FALSE)

# Export 2: Full Extended & Shifted Target Matrix (1996-2025 window)
write.csv(shifted_matrix %>% filter(Year >= 1996 & Year <= 2025), 
          file.path(outputDir, "datWide_1996_2025_extended_shifted.csv"), row.names = FALSE)

cat("\n==============================================================================\n")
cat("FIXED-PARAMETER EXTENSION & SHIFT MATRIX GENERATED SUCCESSFULLY!\n")
cat("Outputs saved to:", outputDir, "\n")
cat(" - datWide_1998_2021_qualified.csv\n")
cat(" - datWide_1996_2025_extended_shifted.csv\n")
cat(" - marss_fits_locked_1998_2021.rds\n")
cat("==============================================================================\n")
