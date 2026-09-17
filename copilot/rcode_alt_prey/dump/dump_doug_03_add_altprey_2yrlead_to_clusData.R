

# ==============================================================================
# Script: doug_03_add_altprey_2yrlead_to_clusData.R
# Purpose: Extend 8 AltPrey series (including Sitka Herring) through 2023 
#          using locked 1998–2021 scaling and append their 2-year leads.
# Output: copilot/outputs_altprey/clusDataDFA_wide_1998_2021_altprey_extended.csv
# ==============================================================================

library(tidyverse)
library(MARSS)
library(lubridate)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & TARGET SPECIFICATIONS
# ------------------------------------------------------------------------------
rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"
outputDir <- file.path(rootdir, "copilot/outputs_altprey")

dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: LOAD RAW OBSERVATIONS (1998–2023) & LOCK 1998–2021 SCALING
# ------------------------------------------------------------------------------
raw_wide_path <- file.path(outputDir, "datWide_1996_2025_qualified.csv")
if (!file.exists(raw_wide_path)) {
  stop("CRITICAL ERROR: Cannot find datWide_1996_2025_qualified.csv in ", outputDir)
}

raw_data <- read.csv(raw_wide_path) %>%
  filter(Year >= 1998 & Year <= 2023)

years_full <- raw_data$Year  # 1998 to 2023

# Function to z-score extended series using ONLY 1998–2021 mean and SD
scale_with_baseline_lock <- function(vec, years) {
  baseline_idx <- which(years >= 1998 & years <= 2021)
  base_mean    <- mean(vec[baseline_idx], na.rm = TRUE)
  base_sd      <- sd(vec[baseline_idx], na.rm = TRUE)
  
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

# Load saved MARSS model fits from baseline 1998-2021 run
marss_fits_path <- file.path(outputDir, "marss_fits.rds")
if (!file.exists(marss_fits_path)) {
  stop("CRITICAL ERROR: Cannot find marss_fits.rds in ", outputDir)
}
marss_fits <- readRDS(marss_fits_path)

extended_results <- list()

# ------------------------------------------------------------------------------
# STEP 2A: BLOCK A - SINGLE-SERIES STRAGGLERS (Squid, Anchovy, Pollock, Sitka Herring)
# ------------------------------------------------------------------------------
straggler_targets <- tribble(
  ~LisaName,               ~raw_col,
  "X04_marketsquid_GAM",   "marketsquid_GAM_2025",
  "X05_anchovy_GAM",       "anchovy_GAM_2025",
  "X13_pollock_age1plus",  "pollockBiomassAIage1plus_predAK_2026",
  "X13_sitkaHerring_EGoA", "sitkaHerring_EGoA"
)

for (i in seq_len(nrow(straggler_targets))) {
  l_name  <- straggler_targets$LisaName[i]
  r_col   <- straggler_targets$raw_col[i]
  
  cat("Processing Straggler Extension for:", l_name, "-> Raw Col:", r_col, "\n")
  
  if (r_col %in% names(raw_data)) {
    raw_vec  <- raw_data[[r_col]]
    
    # Standardize using locked 1998-2021 baseline mean/SD
    z_scaled <- scale_with_baseline_lock(raw_vec, years_full)
    
    # Kalman Smooth via MARSS (fit = FALSE)
    fit_single <- MARSS(matrix(z_scaled, nrow = 1), fit = FALSE, silent = TRUE)
    fit_single$par <- fit_single$start
    kfList <- MARSSkf(fit_single)
    
    extended_results[[l_name]] <- tibble(
      Year = years_full,
      LisaName = l_name,
      value = as.numeric(kfList$xtT)
    )
  } else {
    stop("CRITICAL ERROR: Straggler raw column missing from raw_data: ", r_col)
  }
}

# ------------------------------------------------------------------------------
# STEP 2B: BLOCK B - MULTI-INDICATOR DFAs (Sardine DFA, ChinAbund DFA, Hake DFA, Euph DFA)
# ------------------------------------------------------------------------------
dfa_targets <- tribble(
  ~LisaName,                     ~DFAname,
  "X05_DFA_abundSardine",        "05.ForageFishNCC_DFA1",
  "X09_DFA_ChinAbundSnakeFall",  "09.PredFishNCC_DFA1",
  "X09_DFA_HakeAge5Plus",        "09.PredFishNCC_b_DFA1",
  "X12_DFA_biomassEuphShelfSum", "12.ZooPreyAK_DFA1"
)

for (i in seq_len(nrow(dfa_targets))) {
  l_name <- dfa_targets$LisaName[i]
  d_name <- dfa_targets$DFAname[i]
  
  cat("Processing DFA Projection for:", l_name, "(Fit key:", d_name, ")\n")
  
  fit_obj <- marss_fits[[d_name]]
  if (is.null(fit_obj)) {
    stop("CRITICAL ERROR: MARSS model fit not found in marss_fits.rds for key: ", d_name)
  }
  
  fitted_rows <- rownames(fit_obj$model$data)
  
  # Build matrix for 26 years (1998–2023)
  mat_data <- matrix(NA_real_, nrow = length(fitted_rows), ncol = length(years_full))
  rownames(mat_data) <- fitted_rows
  colnames(mat_data) <- years_full
  
  for (r_idx in seq_along(fitted_rows)) {
    r_name <- fitted_rows[r_idx]
    
    target_col <- case_when(
      r_name %in% names(raw_data) ~ r_name,
      r_name == "abundSardine" && "sardine_GAM_2025" %in% names(raw_data) ~ "sardine_GAM_2025",
      r_name == "abundSardine" && "sardine_NCC" %in% names(raw_data) ~ "sardine_NCC",
      r_name == "abundHerring" && "herring_GAM_2025" %in% names(raw_data) ~ "herring_GAM_2025",
      r_name == "abundHerring" && "herring_NCC" %in% names(raw_data) ~ "herring_NCC",
      TRUE ~ NA_character_
    )
    
    if (!is.na(target_col) && target_col %in% names(raw_data)) {
      mat_data[r_idx, ] <- scale_with_baseline_lock(raw_data[[target_col]], years_full)
    }
  }
  
  # Build model shell over 26-year matrix and assign locked parameters
  marss_model_26 <- MARSS(mat_data, model = list(m = 1), form = "dfa", fit = FALSE, silent = TRUE)
  marss_model_26$par <- fit_obj$par
  
  # Run Kalman Smoother
  kf_ext <- MARSSkf(marss_model_26)
  ext_trend <- as.numeric(kf_ext$xtT[1, ])
  
  # Preserve factor sign convention
  Z_matrix <- fit_obj$par$Z
  maxIndex <- which.max(abs(Z_matrix[, 1]))
  if (Z_matrix[maxIndex, 1] < 0) {
    ext_trend <- -ext_trend
  }
  
  extended_results[[l_name]] <- tibble(
    Year = years_full,
    LisaName = l_name,
    value = ext_trend
  )
}

# Combine all 8 extended series
all_extended_altprey <- bind_rows(extended_results)

# ------------------------------------------------------------------------------
# STEP 3: CONSTRUCT IN-YEAR AND 2-YEAR LEADS & MERGE ONTO BASELINE MATRIX
# ------------------------------------------------------------------------------
cat("\nStep 3: Constructing in-year and 2-year lead columns for target series...\n")

# 1. Prepare In-Year (t) dataframe for 1998–2021
in_year_df <- all_extended_altprey %>%
  filter(Year >= 1998 & Year <= 2021) %>%
  select(Year, LisaName, value)

# 2. Prepare 2-Year Lead (t-2) dataframe for 1998–2021 smolt entry
lead_df <- all_extended_altprey %>%
  mutate(
    Year_SmoltEntry = Year - 2,
    LisaName        = paste0(LisaName, "_2yrLead")
  ) %>%
  filter(Year_SmoltEntry >= 1998 & Year_SmoltEntry <= 2021) %>%
  select(Year = Year_SmoltEntry, LisaName, value)

# 3. Combine both in-year and lead series and pivot to wide format
targets_wide <- bind_rows(in_year_df, lead_df) %>%
  pivot_wider(names_from = LisaName, values_from = value)

# 4. Load baseline 1998-2021 DFA wide matrix
baseline_dfa_path <- file.path(outputDir, "clusDataDFA_wide_1998_2021.csv")
baseline_dfa_wide <- read.csv(baseline_dfa_path)

# Remove any existing duplicate target columns from baseline before joining
clean_baseline <- baseline_dfa_wide %>%
  select(-any_of(names(targets_wide)[names(targets_wide) != "Year"]))

# Left join both in-year and 2yrLead columns onto baseline
final_extended_matrix <- clean_baseline %>%
  left_join(targets_wide, by = "Year")

# Save final dataset
out_csv_path <- file.path(outputDir, "clusDataDFA_wide_1998_2021_altprey_extended.csv")
write.csv(final_extended_matrix, out_csv_path, row.names = FALSE)

cat("\nPipeline Step Complete!")
cat("\nSaved baseline matrix with both in-year and 2yrLead AltPrey columns to:\n", out_csv_path, "\n")



