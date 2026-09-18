
# ==============================================================================
# Script: doug_03b_add_altprey_2yrlead_akpred_smoltYear_to_clusData.R
# Purpose: Extend straggler/DFA series, perform 2-year lead/lag shifts, and 
#          apply interior-only interpolation to stragglers (preserving edge NAs).
# Output: copilot/outputs_altprey/clusDataDFA_wide_1998_2021_altprey_extended_smoltyear_shifted_LisaNames.csv
# ==============================================================================

library(tidyverse)
library(MARSS)
library(imputeTS)

rootdir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"
outputDir <- file.path(rootdir, "copilot/outputs_altprey")
metaDir   <- file.path(rootdir, "Rcode_for_paper/metadata")

# ------------------------------------------------------------------------------
# STEP 1: LOAD STANDARDIZED BASELINE MATRIX & FULL 1996+ RAW OBSERVATIONS
# ------------------------------------------------------------------------------
baseline_path <- file.path(outputDir, "clusDataDFA_wide_1998_2021_LisaNames.csv")
if (!file.exists(baseline_path)) stop("Missing standardized baseline file. Run doug_03a first!")
baseline_dfa_wide <- read.csv(baseline_path, stringsAsFactors = FALSE)

# Read raw dataset starting in 1996 to compute shifts across full window
raw_data <- read.csv(file.path(outputDir, "datWide_1996_2025_qualified.csv")) %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  arrange(Year)

years_full <- raw_data$Year

scale_with_baseline_lock <- function(vec, years) {
  baseline_idx <- which(years >= 1998 & years <= 2021)
  base_mean    <- mean(vec[baseline_idx], na.rm = TRUE)
  base_sd      <- sd(vec[baseline_idx], na.rm = TRUE)
  if (is.na(base_sd) || base_sd == 0) return(vec - base_mean)
  return((vec - base_mean) / base_sd)
}

# ------------------------------------------------------------------------------
# STEP 2: PROCESS EXTENDED DFA TRENDS VIA MARSS (1996–2025)
# ------------------------------------------------------------------------------
marss_fits <- readRDS(file.path(outputDir, "marss_fits.rds"))
extended_results <- list()

dfa_targets <- tribble(
  ~LisaName,                      ~DFAname,
  "X05_DFA_abundSardine",        "05.ForageFishNCC_DFA1",
  "X09_DFA_ChinAbundSnakeFall",  "09.PredFishNCC_DFA1",
  "X09_DFA_HakeAge5Plus",        "09.PredFishNCC_b_DFA1",
  "X12_DFA_biomassEuphShelfSum", "12.ZooPreyAK_DFA1"
)

for (i in seq_len(nrow(dfa_targets))) {
  l_name <- dfa_targets$LisaName[i]
  d_name <- dfa_targets$DFAname[i]
  fit_obj <- marss_fits[[d_name]]
  if (!is.null(fit_obj)) {
    fitted_rows <- rownames(fit_obj$model$data)
    mat_data <- matrix(NA_real_, nrow = length(fitted_rows), ncol = length(years_full))
    rownames(mat_data) <- fitted_rows
    colnames(mat_data) <- years_full
    for (r_idx in seq_along(fitted_rows)) {
      r_name <- fitted_rows[r_idx]
      target_col <- case_when(
        r_name %in% names(raw_data) ~ r_name,
        r_name == "abundSardine" && "sardine_GAM_2025" %in% names(raw_data) ~ "sardine_GAM_2025",
        r_name == "abundHerring" && "herring_GAM_2025" %in% names(raw_data) ~ "herring_GAM_2025",
        TRUE ~ NA_character_
      )
      if (!is.na(target_col) && target_col %in% names(raw_data)) {
        mat_data[r_idx, ] <- scale_with_baseline_lock(raw_data[[target_col]], years_full)
      }
    }
    marss_model_26 <- MARSS(mat_data, model = list(m = 1), form = "dfa", fit = FALSE, silent = TRUE)
    marss_model_26$par <- fit_obj$par
    ext_trend <- as.numeric(MARSSkf(marss_model_26)$xtT[1, ])
    if (fit_obj$par$Z[which.max(abs(fit_obj$par$Z[, 1])), 1] < 0) ext_trend <- -ext_trend
    extended_results[[l_name]] <- tibble(Year = years_full, LisaName = l_name, value = ext_trend)
  }
}

all_extended_altprey <- bind_rows(extended_results)

# ------------------------------------------------------------------------------
# STEP 3: COMPUTE SHIFTS & INTERIOR-ONLY INTERPOLATION FOR STRAGGLERS
# ------------------------------------------------------------------------------
extended_targets_wide <- all_extended_altprey %>%
  filter(Year >= 1996 & Year <= 2025) %>%
  pivot_wider(names_from = LisaName, values_from = value) %>%
  arrange(Year)

prey_lead_cols <- intersect(
  c("X09_DFA_ChinAbundSnakeFall", "X09_DFA_HakeAge5Plus", 
    "X13_pollock_age1plus", "X05_DFA_abundSardine", 
    "X05_anchovy_GAM", "X13_sitkaHerring_EGoA"),
  names(extended_targets_wide)
)

prey_leads_wide <- extended_targets_wide %>%
  mutate(across(all_of(prey_lead_cols), ~ dplyr::lead(.x, 2), .names = "{.col}_2yrLead"))

# Interpolate ONLY interior missing values (e.g. 2010), keeping boundary NAs intact
interpolate_interior_only <- function(vec) {
  if (all(is.na(vec)) || !any(is.na(vec))) return(vec)
  
  valid_indices <- which(!is.na(vec))
  first_valid   <- min(valid_indices)
  last_valid    <- max(valid_indices)
  
  if (first_valid < last_valid) {
    interior_subvec <- vec[first_valid:last_valid]
    if (any(is.na(interior_subvec))) {
      vec[first_valid:last_valid] <- imputeTS::na_interpolation(interior_subvec, option = "linear")
    }
  }
  return(vec)
}

get_raw_straggler_df <- function(data, pattern, col_name_base) {
  matched_col <- names(data)[grep(pattern, names(data), ignore.case = TRUE)][1]
  if (!is.na(matched_col)) {
    scaled_vec <- scale_with_baseline_lock(data[[matched_col]], data$Year)
    clean_vec  <- interpolate_interior_only(scaled_vec)
    return(tibble(Year = data$Year, !!col_name_base := clean_vec))
  }
  return(tibble(Year = data$Year, !!col_name_base := NA_real_))
}

murre_df  <- get_raw_straggler_df(raw_data, "commonmurre|murre_jsoes", "X08_commonMurre_JSOES")
canary_df <- get_raw_straggler_df(raw_data, "canaryrockfish|canary_rockfish", "X09_canaryRockfish")
hseal_df  <- get_raw_straggler_df(raw_data, "harbor_seal_cr", "X10_Harbor_seal_CR_2yrLead") %>%
  arrange(Year) %>%
  mutate(X10_Harbor_seal_CR_smoltYear = dplyr::lag(X10_Harbor_seal_CR_2yrLead, 2))

kw_col <- names(raw_data)[grep("killer.*whale", names(raw_data), ignore.case = TRUE)][1]
if (!is.na(kw_col)) {
  kw_df <- tibble(
    Year = raw_data$Year,
    kw_z = scale_with_baseline_lock(raw_data[[kw_col]], raw_data$Year)
  ) %>%
    arrange(Year) %>%
    mutate(
      killer.whales.nr.bc_2yrlead_2026   = dplyr::lead(kw_z, 2),
      killer.whales.nr.bc_smoltYear_2026 = dplyr::lag(kw_z, 2)
    ) %>%
    select(-kw_z)
} else {
  kw_df <- tibble(Year = raw_data$Year, killer.whales.nr.bc_2yrlead_2026 = NA_real_, killer.whales.nr.bc_smoltYear_2026 = NA_real_)
}

raw_shifts_wide <- murre_df %>%
  left_join(canary_df, by = "Year") %>%
  left_join(hseal_df,  by = "Year") %>%
  left_join(kw_df,     by = "Year")

# Merge baseline data starting from 1996 timeline to compute adult predator lags cleanly
baseline_1996_extended <- raw_data %>%
  select(Year) %>%
  left_join(baseline_dfa_wide, by = "Year") %>%
  arrange(Year) %>%
  mutate(
    X10_DFA_ssl.est.wholerange_smoltYear = dplyr::lag(X10_DFA_ssl.est.wholerange_2yrLead, 2),
    X10_Californian_s_l_smoltYear_WS      = dplyr::lag(X10_Californian_s_l_2yrLead_WS, 2),
    X10_Harbour_s_smoltYear_WS            = dplyr::lag(X10_Harbour_s_2yrLead_WS, 2),
    X10_Northern_f_s_smoltYear_WS         = dplyr::lag(X10_Northern_f_s_2yrLead_WS, 2),
    
    X15_ArrowtoothFlounderBiomass_predAK_smoltYear = dplyr::lag(X15_ArrowtoothFlounderBiomass_predAK, 2),
    X15_sablefishBiomass_predAK_smoltYear          = dplyr::lag(X15_sablefishBiomass_predAK, 2),
    X15_sablefishRecruitment_predAK_smoltYear      = dplyr::lag(X15_sablefishRecruitment_predAK, 2),
    X15_PacificCodBiomass_predAK_smoltYear         = dplyr::lag(X15_PacificCodBiomass_predAK, 2),
    X15_spinyDogfishGoA_predAK_smoltYear          = dplyr::lag(X15_spinyDogfishGoA_predAK, 2),
    X15_DFA_sleeperSharks_smoltYear               = dplyr::lag(X15_DFA_sleeperSharks, 2),
    X15_salmonSharkGoA_predAK_smoltYear           = dplyr::lag(X15_salmonSharkGoA_predAK, 2),
    X15_halibutBiomassAge8plus_smoltYear_predAK   = dplyr::lag(X15_halibutBiomassAge8plus_2yrLead_predAK, 2)
  )

clean_baseline <- baseline_1996_extended %>%
  select(-any_of(names(prey_leads_wide)[names(prey_leads_wide) != "Year"])) %>%
  select(-any_of(names(raw_shifts_wide)[names(raw_shifts_wide) != "Year"]))

final_extended_matrix <- clean_baseline %>%
  left_join(prey_leads_wide, by = "Year") %>%
  left_join(raw_shifts_wide, by = "Year") %>%
  filter(Year >= 1998 & Year <= 2021)

# ------------------------------------------------------------------------------
# STEP 4: EXPORT STANDARDIZED SHIFTED MATRIX
# ------------------------------------------------------------------------------
out_altprey_path <- file.path(outputDir, "clusDataDFA_wide_1998_2021_altprey_extended_smoltyear_shifted_LisaNames.csv")
out_meta_path    <- file.path(metaDir,   "clusDataDFA_wide_1998_2021_altprey_extended_smoltyear_shifted_LisaNames.csv")

write.csv(final_extended_matrix, out_altprey_path, row.names = FALSE)
write.csv(final_extended_matrix, out_meta_path,    row.names = FALSE)

cat("\nPipeline Complete!\nSaved clean standardized matrix to:\n - ", out_altprey_path, "\n - ", out_meta_path, "\n")


