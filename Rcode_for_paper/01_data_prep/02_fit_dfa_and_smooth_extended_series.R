


# ==============================================================================
# Script: 02_fit_dfa_and_smooth_extended_series.R
# Purpose: 1. Fit MARSS DFA state-space models and smooth extended series (1995-2025).
#          2. Re-impose boundary NAs from completeness.csv BEFORE shifting.
#          3. Apply 2-year cohort leads/lags for smoltyr and adultyr alignment.
#          4. Export cleanly differentiated clusDataDFA_wide files with ALPHANUMERIC
#             column ordering (Year, X01_..., X02_..., ..., X16_SAR).
# Outputs: 
#   - metadata/clusDataDFA_wide_1995_2025_extended_dfaCols.csv
#   - metadata/clusDataDFA_wide_1996_2025_extended_shifted_dfaCols.csv
#   - metadata/clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(MARSS)

# --- Configuration & Paths ---
rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
path      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"
outputDir <- file.path(rootdir, "metadata")

raw_file          <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw_shortnames.csv")
guildFile         <- file.path(path, "guildsWithExclude.csv")
crosswalk_file    <- file.path(outputDir, "master_name_crosswalk.csv")
completeness_file <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/DAG1A_long/completeness.csv"

# Fallback check for raw file
if (!file.exists(raw_file)) {
  raw_file_alt <- file.path(outputDir, "datWide_1995_2025_reprocessed_raw.csv")
  if (file.exists(raw_file_alt)) raw_file <- raw_file_alt else stop("Missing raw datWide input file!")
}

raw_data      <- read.csv(raw_file, stringsAsFactors = FALSE)
guilds_master <- read.csv(guildFile, stringsAsFactors = FALSE)
complete_df   <- read.csv(completeness_file, stringsAsFactors = FALSE)
crosswalk_df  <- read.csv(crosswalk_file, stringsAsFactors = FALSE)

baseline_years <- 1998:2021
full_years     <- raw_data$Year

# Lock scaling parameters to 1998-2021 baseline
scale_baseline_lock <- function(vec, yrs) {
  base_idx  <- which(yrs %in% baseline_years)
  base_mean <- mean(vec[base_idx], na.rm = TRUE)
  base_sd   <- sd(vec[base_idx], na.rm = TRUE)
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

guilds <- guilds_master %>%
  mutate(guild = latestGuild) %>%
  filter(shortName %in% names(raw_data), guild != "") %>%
  rename(guildSEMnode = SEMnode) %>%
  select(shortName, guild, guildSEMnode)

unique_guilds <- guilds %>% distinct(guild, guildSEMnode) %>% arrange(guildSEMnode, guild)

extended_results  <- list()
marss_fits_locked <- list()

cat("Fitting DFAs and MARSS Kalman smoothing across guilds...\n")

for (i in seq_len(nrow(unique_guilds))) {
  g_row  <- unique_guilds[i, ]
  g_name <- g_row$guild
  g_node <- g_row$guildSEMnode
  g_inds <- guilds %>% filter(guild == g_name, guildSEMnode == g_node) %>% pull(shortName)
  
  if (length(g_inds) > 1) {
    # --- MULTI-INDICATOR GUILD (DFA) ---
    base_mat <- raw_data %>%
      filter(Year %in% baseline_years) %>%
      select(all_of(g_inds)) %>%
      mutate(across(everything(), ~ as.vector(scale(.)))) %>%
      t()
    
    fit_base <- MARSS(base_mat, model = list(m = 1), form = "dfa", method = "BFGS", silent = TRUE)
    
    # Enforce sign alignment
    Z_est <- fit_base$par$Z
    max_idx <- which.max(abs(Z_est[, 1]))
    if (Z_est[max_idx, 1] < 0) {
      fit_base$par$Z[, 1]  <- -fit_base$par$Z[, 1]
      fit_base$states[1, ] <- -fit_base$states[1, ]
    }
    
    marss_fits_locked[[paste0(g_name, "_DFA1")]] <- fit_base
    
    # Project full window (1995-2025)
    ext_mat <- raw_data %>%
      select(all_of(g_inds)) %>%
      mutate(across(everything(), ~ scale_baseline_lock(., full_years))) %>%
      t()
    
    marss_proj <- MARSS(ext_mat, model = list(m = 1), form = "dfa", fit = FALSE, silent = TRUE)
    marss_proj$par <- fit_base$par
    
    ext_trend <- as.numeric(MARSSkf(marss_proj)$xtT[1, ])
    extended_results[[paste0(g_name, "_DFA1")]] <- tibble(Year = full_years, val = ext_trend)
    
  } else if (length(g_inds) == 1) {
    # --- SINGLE INDICATOR / STRAGGLER ---
    ind_name <- g_inds[1]
    
    base_vec <- scale(raw_data %>% filter(Year %in% baseline_years) %>% pull(!!sym(ind_name))) %>% as.vector()
    fit_single <- MARSS(base_vec, fit = FALSE, silent = TRUE)
    fit_single$par <- fit_single$start
    
    marss_fits_locked[[paste0(g_name, "_smoothed")]] <- fit_single
    
    ext_vec <- scale_baseline_lock(raw_data[[ind_name]], full_years)
    ext_proj <- MARSS(ext_vec, fit = FALSE, silent = TRUE)
    ext_proj$par <- fit_single$par
    
    ext_trend <- as.numeric(t(MARSSkf(ext_proj)$xtT))
    extended_results[[paste0(g_name, "_smoothed")]] <- tibble(Year = full_years, val = ext_trend)
  }
}

# Save locked MARSS container
saveRDS(marss_fits_locked, file.path(outputDir, "marss_fits_locked_1998_2021.rds"))

# Pivot wide to create unshifted clusDataDFA matrix
clus_wide_dfaCols <- bind_rows(
  lapply(names(extended_results), function(nm) {
    extended_results[[nm]] %>% mutate(dfa_cols = nm)
  })
) %>%
  pivot_wider(id_cols = Year, names_from = dfa_cols, values_from = val) %>%
  arrange(Year)

# Export Unshifted Unmasked clusDataDFA (Alphanumerically ordered columns)
dfa_cols_sorted <- c("Year", sort(setdiff(names(clus_wide_dfaCols), "Year")))
clus_wide_dfaCols <- clus_wide_dfaCols %>% select(all_of(dfa_cols_sorted))

write.csv(clus_wide_dfaCols, file.path(outputDir, "clusDataDFA_wide_1995_2025_extended_dfaCols.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: RE-IMPOSE BOUNDARY NAs FROM COMPLETENESS.CSV BEFORE SHIFTING
# ------------------------------------------------------------------------------
cat("\nRe-imposing true missing/boundary NAs from completeness.csv...\n")

complete_df <- complete_df %>%
  mutate(CalendarYear = year(ymd(date)))

comp_cols <- setdiff(names(complete_df), c("date", "complete", "Year", "year", "X", "CalendarYear"))
masked_clus_wide <- clus_wide_dfaCols

for (c_col in comp_cols) {
  clean_c_col <- gsub("^X", "", c_col)
  matched_cols <- names(masked_clus_wide)[gsub("^X", "", names(masked_clus_wide)) == clean_c_col]
  
  if (length(matched_cols) > 0) {
    missing_cal_years <- complete_df$CalendarYear[is.na(complete_df[[c_col]]) | complete_df[[c_col]] == FALSE]
    
    if (length(missing_cal_years) > 0) {
      for (m_col in matched_cols) {
        idx <- which(masked_clus_wide$Year %in% missing_cal_years)
        masked_clus_wide[idx, m_col] <- NA
      }
    }
  }
}

# ------------------------------------------------------------------------------
# STEP 2: APPLY 2-YEAR COHORT SHIFTS ON MASKED MATRIX
# ------------------------------------------------------------------------------
cat("Applying 2-year cohort shifts (lead 2) for smoltyr and adultyr alignment...\n")

shifted_dfaCols <- masked_clus_wide %>%
  # 1. Adult Predators (X10, X15): Original column is _adultyr. Derive _smoltyr via lead(2)
  rename_with(~ paste0(.x, "_adultyr"), contains("10.PredMammal") | contains("15.PredFish")) %>%
  mutate(across((contains("10.PredMammal") | contains("15.PredFish")) & ends_with("_adultyr"), 
                ~ dplyr::lead(.x, 2), 
                .names = "{gsub('_adultyr', '_smoltyr', .col)}")) %>%
  
  # 2. Alternate Prey (X05, X13, X14): Original column is _smoltyr. Derive _adultyr via lead(2)
  rename_with(~ paste0(.x, "_smoltyr"), contains("05.Forage") | contains("13.FishPrey") | contains("14.CompAK")) %>%
  mutate(across((contains("05.Forage") | contains("13.FishPrey") | contains("14.CompAK")) & ends_with("_smoltyr"), 
                ~ dplyr::lead(.x, 2), 
                .names = "{gsub('_smoltyr', '_adultyr', .col)}"))

# Sort dfaCols alphanumerically
shifted_dfa_sorted <- c("Year", sort(setdiff(names(shifted_dfaCols), "Year")))
shifted_dfaCols <- shifted_dfaCols %>% select(all_of(shifted_dfa_sorted))

# Export Shifted dfaCols Matrix (1996-2025)
out_dfaCols <- file.path(outputDir, "clusDataDFA_wide_1996_2025_extended_shifted_dfaCols.csv")
write.csv(shifted_dfaCols %>% filter(Year >= 1996 & Year <= 2025), out_dfaCols, row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 3: MAP HEADERS TO LISANAMES & SORT ALPHANUMERICALLY
# ------------------------------------------------------------------------------
cat("Translating headers to LisaNames and sorting alphanumerically...\n")

crosswalk_df <- crosswalk_df %>%
  mutate(LisaName = gsub("_+$", "", trimws(LisaName)))

lookup_map <- setNames(crosswalk_df$LisaName, crosswalk_df$dfa_cols)

current_cols <- names(shifted_dfaCols)
new_cols     <- current_cols

for (i in seq_along(current_cols)) {
  col <- current_cols[i]
  if (col == "Year") next
  
  suffix <- ""
  base_col <- col
  if (grepl("_smoltyr$", col)) {
    suffix <- "_smoltyr"
    base_col <- gsub("_smoltyr$", "", col)
  } else if (grepl("_adultyr$", col)) {
    suffix <- "_adultyr"
    base_col <- gsub("_adultyr$", "", col)
  }
  
  base_col_x <- if (!startsWith(base_col, "X")) paste0("X", base_col) else base_col
  
  if (base_col_x %in% names(lookup_map)) {
    new_cols[i] <- paste0(lookup_map[[base_col_x]], suffix)
  } else if (base_col %in% names(lookup_map)) {
    new_cols[i] <- paste0(lookup_map[[base_col]], suffix)
  }
}

shifted_LisaNames <- shifted_dfaCols
names(shifted_LisaNames) <- new_cols

# Sort LisaNames alphanumerically: Year first, then X01_... through X16_SAR
lisa_cols_sorted <- c("Year", sort(setdiff(names(shifted_LisaNames), "Year")))
shifted_LisaNames <- shifted_LisaNames %>% select(all_of(lisa_cols_sorted))

# Export Final LisaName Matrix (1996-2025)
out_lisa <- file.path(outputDir, "clusDataDFA_wide_1996_2025_extended_shifted_LisaNames.csv")
write.csv(shifted_LisaNames %>% filter(Year >= 1996 & Year <= 2025), out_lisa, row.names = FALSE)

# ------------------------------------------------------------------------------
# VERIFICATION DIAGNOSTIC
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("PROCESS COMPLETE! EXPORTED PRODUCTS SUMMARY:\n")
cat("==============================================================================\n")
cat("1. Unshifted Smoothed DFA Matrix:\n   ", file.path(outputDir, "clusDataDFA_wide_1995_2025_extended_dfaCols.csv"), "\n")
cat("2. Shifted DFA Matrix (dfaCols):\n   ", out_dfaCols, "\n")
cat("3. Shifted DFA Matrix (LisaNames - Alphanumeric Order):\n   ", out_lisa, "\n")
cat(sprintf("\nTotal NAs in Final Shifted LisaName Matrix: %d\n", sum(is.na(shifted_LisaNames))))
cat("First 5 Columns in Output:\n")
print(names(shifted_LisaNames)[1:5])
cat("==============================================================================\n")