# ==============================================================================
# Script: 11_filter_valid_ak_top_models.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Filter sweep CSV results to retain strictly valid Alaska (AK) 
#          predator models (SAR ~ CPUE + Peff) meeting Delta AIC < -4 and p < 0.05.
# Output: Rcode_for_paper/Routput_for_paper/data/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

# Select sweep file paths
sweep_24yr_path <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA_akferris.csv")
sweep_19yr_path <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_19yr_noNA_dougsmoothed.csv")

# ------------------------------------------------------------------------------
# STEP 1: DEFINE FILTERING RULES
# ------------------------------------------------------------------------------
blacklist <- c("canaryRockfish", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")
blacklist_pattern <- paste(blacklist, collapse = "|")

# Helper function to classify valid Alaska post-CPUE predators
is_valid_ak_predator <- function(pred_name) {
  # Drop blacklist
  if (grepl(blacklist_pattern, pred_name)) return(FALSE)
  
  # Valid AK predators start with X15 or are AK sea lion series (X10_DFA_ssl / X11_ssl)
  is_ak <- grepl("^X15_", pred_name) || grepl("^(X10_DFA_ssl|X11_ssl)", pred_name)
  
  # Ensure NCC smolt/avian/fish drivers are excluded
  is_ncc <- grepl("^(X08_|X09_|X10_|X11_DFA_Harbour_p_WS)", pred_name) && !grepl("ssl", pred_name)
  
  return(is_ak && !is_ncc)
}

# ------------------------------------------------------------------------------
# STEP 2: FILTER & RANK FUNCTION
# ------------------------------------------------------------------------------
process_ak_sweep <- function(file_path, dataset_label) {
  if (!file.exists(file_path)) {
    warning("Sweep file not found: ", file_path)
    return(NULL)
  }
  
  df <- read.csv(file_path, stringsAsFactors = FALSE) %>%
    # Deduplicate predator entries
    distinct(predator, .keep_all = TRUE) %>%
    # Filter for valid AK predators only
    filter(sapply(predator, is_valid_ak_predator)) %>%
    # Quality thresholds: Delta AIC <= -4 (or delta_aic < -4) AND p < 0.05
    filter(delta_aic <= -2.0, pvalue < 0.05) %>%
    # Sort by best Delta AIC
    arrange(desc(final_r2)) %>%
    mutate(dataset = dataset_label) %>%
    as_tibble()
  
  return(df)
}

# ------------------------------------------------------------------------------
# STEP 3: EXECUTE FILTERING
# ------------------------------------------------------------------------------
top_ak_24yr <- process_ak_sweep(sweep_24yr_path, "24-Year Series")
top_ak_19yr <- process_ak_sweep(sweep_19yr_path, "19-Year Series")

cat("==============================================================================\n")
cat("            TOP VALID ALASKA MODELS (19-YEAR DATASET)                         \n")
cat("==============================================================================\n\n")
if (!is.null(top_ak_19yr)) {
  print(top_ak_19yr %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, delta_aic, final_r2), n = 50)
}

cat("\n==============================================================================\n")
cat("            TOP VALID ALASKA MODELS (24-YEAR DATASET)                         \n")
cat("==============================================================================\n\n")
if (!is.null(top_ak_24yr)) {
  print(top_ak_24yr %>% select(predator, selected_prey, fitted_weights, beta_peff, pvalue, delta_aic, final_r2), n = 50)
}

# Combine both into a single summary
combined_ak_top <- bind_rows(top_ak_24yr,top_ak_19yr)

# Save to output CSV
out_csv <- file.path(data_out_dir, "Top_Valid_AK_Models_Filtered_Summary.csv")
write_csv(combined_ak_top, out_csv)

cat(sprintf("\nFiltered Alaska top models successfully written to: %s\n", out_csv))



#results----------
==============================================================================
  > cat("            TOP VALID ALASKA MODELS (19-YEAR DATASET)                         \n")
TOP VALID ALASKA MODELS (19-YEAR DATASET)                         
# A tibble: 14 × 7
predator                                          selected_prey                                                 fitted_weights beta_peff       pvalue delta_aic final_r2
<chr>                                             <chr>                                                         <chr>              <dbl>        <dbl>     <dbl>    <dbl> 1 X15_sablefishBiomass_predAK_smoltyr               X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr  0.9; 1            -7.07  0.0000000497    -17.3     0.718
2 X15_salmonSharkBSAI_predAK_smoltyr                X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr  0.7; 0.9          -2.47  0.000000457     -15.2     0.690
3 X15_halibutBiomassAge8plus_2yrLead_predAK_adultyr X21_sst_egoa_junjulaug_adultyr + X13_capelin_WGoA_adultyr     0.4; 1            -0.962 0.0000118        -9.74    0.640
4 X15_ArrowtoothFlounderBiomass_predAK_adultyr      X21_sst_egoa_junjulaug_adultyr + X13_capelin_WGoA_adultyr     0.2; 1            -0.731 0.0000213       -11.5     0.629
5 X15_DFA_sleeperSharks_adultyr                     X14_pinkSalmonNorthAmerica_adultyr + X13_capelin_WGoA_adultyr 0.4; 1            -0.735 0.0000345        -4.85    0.619
6 X15_salmonSharkGoA_predAK_adultyr                 X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       0.9; 1            -3.68  0.000371         -8.30    0.566
7 X15_spinyDogfishBSAI_predAK_smoltyr               X13_capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr             1; 1              -0.962 0.000787         -8.10    0.545
8 X15_halibutBiomassAge8plus_2yrLead_predAK_smoltyr X13_sitkaHerring_EGoA_smoltyr + X13_capelin_WGoA_smoltyr      1; 1              -0.569 0.00582          -4.78    0.483
9 X15_ArrowtoothFlounderBiomass_predAK_smoltyr      X13_capelin_WGoA_smoltyr                                      1                 -0.383 0.00851          -3.65    0.469
10 X15_salmonSharkBSAI_predAK_adultyr                X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       1; 1              -0.711 0.00928          -5.12    0.466
11 X15_salmonSharkGoA_predAK_smoltyr                 X13_sitkaHerring_EGoA_smoltyr + X13_capelin_WGoA_smoltyr      1; 1              -1.38  0.0108           -5.58    0.460
12 X15_sablefishBiomass_predAK_adultyr               X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       1; 1              -1.98  0.0214           -3.39    0.434
13 X15_PacificCodBiomass_predAK_adultyr              X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       0.7; 1            -0.623 0.0239           -4.52    0.429
14 X15_spinyDogfishBSAI_predAK_adultyr               X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr       0.7; 1            -0.705 0.0352           -3.86    0.413
> cat("\n==============================================================================\n")

==============================================================================
  > cat("            TOP VALID ALASKA MODELS (24-YEAR DATASET)                         \n")
TOP VALID ALASKA MODELS (24-YEAR DATASET)                         
# A tibble: 7 × 7
predator                                     selected_prey                                                fitted_weights beta_peff   pvalue delta_aic final_r2
<chr>                                        <chr>                                                        <chr>              <dbl>    <dbl>     <dbl>    <dbl>
  1 X15_ArrowtoothFlounderBiomass_predAK_adultyr X13_mid_il_capelin_adultyr                                   1                 -0.479 0.000343     -8.75    0.544
2 X15_sablefishBiomass_predAK_smoltyr          X13_pollock_age1plus_smoltyr + X13_sitkaHerring_EGoA_smoltyr 0.8; 1            -2.53  0.00392      -6.62    0.480
3 X15_ArrowtoothFlounderBiomass_predAK_smoltyr X13_mid_il_capelin_smoltyr                                   1                 -0.270 0.0246       -2.37    0.422
4 X15_sablefishBiomass_predAK_adultyr          X13_pollock_age1plus_adultyr + X13_mid_il_capelin_adultyr    1; 1              -1.70  0.0248       -3.45    0.422
5 X11_ssl_seak_pup_pred                        X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr   1; 1              -0.366 0.0250       -3.05    0.421
6 X10_DFA_ssl.est.wholerange_2yrLead_smoltyr   X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr         1; 0.8            -1.61  0.0393       -2.28    0.406
7 X15_salmonSharkGoA_predAK_smoltyr            X13_sitkaHerring_EGoA_smoltyr + X13_mid_il_capelin_smoltyr   1; 1              -1.09  0.0404       -3.81    0.404