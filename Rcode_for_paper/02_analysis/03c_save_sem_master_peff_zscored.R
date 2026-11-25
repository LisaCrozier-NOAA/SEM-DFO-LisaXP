# ==============================================================================
# Script: 02_analysis/03c_save_sem_master_peff_zscored.R
# Purpose: Add Z-score standardized columns (z_Peff_*) to peff_ts_df
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")

peff_ts_path <- file.path(data_out_dir, "obj3_data_peff_timeseries_24yr.csv")
if (!file.exists(peff_ts_path)) {
  stop("Missing Objective 3 time series file at: ", peff_ts_path)
}

peff_ts_df <- read_csv(peff_ts_path, show_col_types = FALSE)

# Identify all Peff metric columns
peff_cols <- names(peff_ts_df)[grepl("^Peff_", names(peff_ts_df))]

# Add z-scaled columns: z_Peff_*
peff_ts_df_z <- peff_ts_df %>%
  mutate(across(all_of(peff_cols), ~ as.numeric(scale(.x)), .names = "z_{.col}"))

# Save updated dataset overwriting or as extended version
 write_csv(peff_ts_df_z, peff_ts_path)
cat(sprintf("Successfully added %d z-standardized Peff columns to %s\n", 
            length(peff_cols), peff_ts_path))

# Quick inspection of original vs z-scaled values
peff_ts_df_z %>%
  select(any_of(c("year", "smoltyr", peff_cols[1], paste0("z_", peff_cols[1])))) %>%
  head()

# ==============================================================================
# Precursor Script 1: Audit df_raw NAs & Merge Clarified Panel A Peff Columns
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
})

proj_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_dir     <- file.path(proj_dir, "Outputs/data")
alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

# ------------------------------------------------------------------------------
# 1. LOAD & AUDIT df_raw (NO CROPPING)
# ------------------------------------------------------------------------------
df_raw <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)

cat("df_raw dimensions:", nrow(df_raw), "rows x", ncol(df_raw), "cols\n")
cat("Years covered in df_raw:", paste(range(df_raw$year, na.rm = TRUE), collapse = " - "), "\n\n")

# Audit NA counts across every column in df_raw
na_audit_raw <- tibble(
  column      = names(df_raw),
  n_valid     = map_int(df_raw, ~ sum(!is.na(.x))),
  n_missing   = map_int(df_raw, ~ sum(is.na(.x))),
  first_year  = map_dbl(df_raw, ~ if (all(is.na(.x))) NA_real_ else min(df_raw$year[!is.na(.x)])),
  last_year   = map_dbl(df_raw, ~ if (all(is.na(.x))) NA_real_ else max(df_raw$year[!is.na(.x)]))
) %>%
  arrange(n_missing, column)

print(as.data.frame(na_audit_raw), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. CLARIFY PANEL A (24-YR) Peff COLUMNS WITH EXPLICIT MECHANISM SUFFIXES
#    Suffix format: "__peff_<Mechanism>" (easy to strip via sub("__peff_.*$", "", x))
#     SELECT z_Peff_* COLUMNS
# ------------------------------------------------------------------------------
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

if ("Year" %in% names(peff_df) && !"year" %in% names(peff_df)) {
  peff_df <- peff_df %>% rename(year = Year)
}

# Exact Panel A mapping from Table 3 (no "Joint" or "Only")
panelA_zpeff_rename <- c(
  "z_Peff_ssl_seak_pup_pred_smoltyr"                  = "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z",
  "z_Peff_harbour_seal_DFO_WCVI_smoltyr"              = "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z",
  "z_Peff_sablefishAK_biomass3plus_smoltyr"           = "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z",
  "z_Peff_Killer_whales_NR_BC_adultyr"                = "Killer_whales_NR_BC_adultyr__peff_Temp_z",
  "z_Peff_DFA_Steller_s_l_WS_adultyr"                 = "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z",
  "z_Peff_DFA_sleeperSharkBSAI_predAK_smoltyr"        = "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z",
  "z_Peff_pollockGOA_SSB_smoltyr"                     = "pollockGOA_SSB_smoltyr__peff_TempPrey2_z",
  "z_Peff_Harbour_s_WS_adultyr"                       = "Harbour_s_WS_adultyr__peff_Base_z",
  "z_Peff_Northern_f_s_WS_adultyr"                    = "Northern_f_s_WS_adultyr__peff_Base_z",
  "z_Peff_harbour_seal_DFO_Coastwide_smoltyr"         = "harbour_seal_DFO_Coastwide_smoltyr__peff_TempPrey1_z",
  "z_Peff_DFA_AllSeaLionsEMB_adultyr"                 = "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z",
  "z_Peff_salmonSharkGoA_predAK_smoltyr"              = "salmonSharkGoA_predAK_smoltyr__peff_TempPrey1_z",
  "z_Peff_DFA_halibutBiomassAge8plus_predAK_adultyr"  = "DFA_halibutBiomassAge8plus_predAK_adultyr__peff_TempPrey1_z"
)

peff_clean <- peff_df %>%
  select(year, all_of(names(panelA_zpeff_rename))) %>%
  rename_with(~ panelA_zpeff_rename[.x], .cols = all_of(names(panelA_zpeff_rename)))

# Merge with full df_raw (retaining ALL rows and ALL time series) and ensure z-scoring on df_raw columns
sem_master_data <- df_raw %>%
  mutate(across(-year, ~ as.numeric(scale(.x)))) %>%
  full_join(peff_clean, by = "year") %>%
  arrange(year)

master_out_path <- file.path(data_dir, "sem_master_peff_zscored.csv")
write_csv(sem_master_data, master_out_path)
cat("\nSaved complete master dataset (", nrow(sem_master_data), "rows x", ncol(sem_master_data),
    "cols) to:\n", master_out_path, "\n")
